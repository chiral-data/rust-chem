use crate::buffers::BufferManager;
use crate::context::GpuContext;
use crate::error::GpuError;
use bytemuck::{Pod, Zeroable};
use chemcore::molecule::Molecule;
use wgpu;

// Every kernel in the batch pipeline (init_invariants_batch,
// morgan_iteration_batch, dedup_environments_batch, copy_invariants,
// copy_neighborhoods) is @workgroup_size(256), and each is dispatched with
// ceil(N / 256) workgroups, where N is num_molecules (dedup), total_atoms
// (init/iteration/copy), or total_neighborhood_words (copy_neighborhoods).
// WebGPU caps dispatch group count at 65,535 per dimension, so a single GPU
// dispatch cannot let any of those three counts exceed 256 * 65,535 (see #27
// for num_molecules, #33 for total_atoms/total_neighborhood_words).
// generate_fingerprints_batch transparently splits larger batches across
// multiple dispatches, keeping every chunk under all three caps at once,
// instead of erroring.
const WORKGROUP_SIZE: usize = 256;
const MAX_MOLECULES_PER_DISPATCH: usize = WORKGROUP_SIZE * 65_535;
const MAX_ATOMS_PER_DISPATCH: u32 = (WORKGROUP_SIZE * 65_535) as u32;
const MAX_NEIGHBORHOOD_WORDS_PER_DISPATCH: u32 = (WORKGROUP_SIZE * 65_535) as u32;

/// Per-dispatch caps a single GPU dispatch must stay under simultaneously
/// (see the constants above). Bundled into one struct so
/// `generate_fingerprints_batch_chunked` takes one parameter instead of
/// three, and so tests can shrink individual caps independently to exercise
/// each chunking trigger without needing tens of millions of atoms.
#[derive(Clone, Copy)]
struct DispatchLimits {
    max_molecules: usize,
    max_atoms: u32,
    max_neighborhood_words: u32,
}

const DISPATCH_LIMITS: DispatchLimits = DispatchLimits {
    max_molecules: MAX_MOLECULES_PER_DISPATCH,
    max_atoms: MAX_ATOMS_PER_DISPATCH,
    max_neighborhood_words: MAX_NEIGHBORHOOD_WORDS_PER_DISPATCH,
};

/// Atom count and neighborhood-word count `mol` contributes to a single
/// dispatch's totals (mirrors the per-molecule accumulation in
/// `generate_fingerprints_batch_single_dispatch`).
fn molecule_dispatch_cost(mol: &Molecule) -> (u32, u32) {
    let num_atoms = mol.num_atoms() as u32;
    let bond_words = (mol.num_bonds() as u32).div_ceil(32).max(1);
    (num_atoms, num_atoms * bond_words)
}

// ============================================================================
// BATCH GPU DATA STRUCTURES
// ============================================================================

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct BatchParams {
    num_molecules: u32,
    max_atoms: u32,
    radius: u32,
    fp_size: u32,
    fp_words: u32,
    use_chirality: u32,
    use_bond_types: u32,
    only_nonzero_invariants: u32,
    // Aggregate sizes needed by the shader to compute the combined `scratch`
    // buffer's region offsets (see morgan.wgsl) — not derivable from the
    // fields above alone, since they depend on each molecule's own bond
    // count.
    total_neighborhood_words: u32,
    total_seen_words: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct MoleculeOffset {
    atom_start: u32,
    atom_count: u32,
    bond_start: u32,
    bond_count: u32,
    adj_start: u32,
    fp_start: u32,
    /// Start word index in the (ping-ponged) neighborhood-bitset buffers,
    /// (bond_count + 31) / 32 words per atom.
    neighborhood_start: u32,
    /// Start word index in the persistent seen_neighborhoods buffer,
    /// capacity atom_count * radius * bond_words for this molecule.
    seen_start: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct GpuAtom {
    atomic_number: u32,
    degree: u32,
    is_aromatic: u32,
    chirality: u32,
    formal_charge: i32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct GpuBond {
    atom1: u32,
    atom2: u32,
    order: u32,
    stereo: u32,
    is_aromatic: u32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct GpuNeighbor {
    atom_idx: u32,
    bond_idx: u32,
}

// Create iteration params buffer
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct IterationParams {
    current_layer: u32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

// ============================================================================
// GPU MORGAN FINGERPRINT GENERATOR
// Processes multiple molecules in a single GPU dispatch for maximum performance
// ============================================================================

/// Cheaply `Clone`: every field is an `Arc`-backed wgpu handle (or a
/// `GpuContext`, itself cheaply `Clone` for the same reason), so cloning
/// shares the same pipelines/device rather than rebuilding them.
#[derive(Clone)]
pub struct GpuMorganFingerprint {
    ctx: GpuContext,
    init_pipeline: wgpu::ComputePipeline,
    iteration_pipeline: wgpu::ComputePipeline,
    dedup_pipeline: wgpu::ComputePipeline,
    copy_pipeline: wgpu::ComputePipeline,
    copy_neighborhoods_pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}

impl GpuMorganFingerprint {
    pub fn new() -> Result<Self, GpuError> {
        Self::from_context(GpuContext::new()?)
    }

    /// Async twin of [`Self::new`] for callers that can't block the current
    /// thread while waiting on adapter/device requests (namely wasm32).
    pub async fn new_async() -> Result<Self, GpuError> {
        Self::from_context(GpuContext::new_async().await?)
    }

    /// Build from an already-initialized [`GpuContext`] instead of creating a
    /// new one. Pipeline/shader setup is entirely synchronous (only
    /// obtaining the `GpuContext` itself needs async adapter/device
    /// requests), so both [`Self::new`] and [`Self::new_async`] share this
    /// once they have one.
    pub(crate) fn from_context(ctx: GpuContext) -> Result<Self, GpuError> {
        // The batch bind group layout below uses this many storage-buffer
        // bindings in its one compute stage. Some WebGPU backends (browsers)
        // only guarantee wgpu's own baseline (8) — real hardware has been
        // observed reporting exactly that (Firefox/Apple Silicon via Metal,
        // see #50) — so this has to be checked before attempting to build
        // the pipeline: wgpu's bind-group/pipeline creation calls aren't
        // fallible at the Rust API level on the WebGPU backend (errors
        // surface asynchronously, e.g. via the browser's console or a lost
        // device), so requesting more bindings than the adapter supports
        // doesn't error here — it silently creates invalid pipelines that
        // hang forever the first time something tries to actually read back
        // their output.
        const REQUIRED_STORAGE_BUFFERS_PER_STAGE: u32 = 6;
        let max_storage_buffers = ctx.limits().max_storage_buffers_per_shader_stage;
        if max_storage_buffers < REQUIRED_STORAGE_BUFFERS_PER_STAGE {
            return Err(GpuError::OperationFailed(format!(
                "adapter supports only {} storage buffers per compute stage, but the Morgan batch shader needs {}",
                max_storage_buffers, REQUIRED_STORAGE_BUFFERS_PER_STAGE
            )));
        }

        let shader_source = include_str!("shaders/morgan.wgsl");

        let shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("morgan_batch_shader"),
                source: wgpu::ShaderSource::Wgsl(shader_source.into()),
            });

        // Create bind group layout (8 bindings for batch Morgan: 2 uniform +
        // 6 storage — see morgan.wgsl's binding-count comment for why the
        // 9 originally-separate read-write buffers are now packed into one
        // combined `scratch` binding).
        let bind_group_layout =
            ctx.device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("morgan_batch_bind_group_layout"),
                    entries: &[
                        // @binding(0): params (uniform)
                        wgpu::BindGroupLayoutEntry {
                            binding: 0,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Uniform,
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(1): molecule_offsets
                        wgpu::BindGroupLayoutEntry {
                            binding: 1,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(2): atoms
                        wgpu::BindGroupLayoutEntry {
                            binding: 2,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(3): bonds
                        wgpu::BindGroupLayoutEntry {
                            binding: 3,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(4): adjacency
                        wgpu::BindGroupLayoutEntry {
                            binding: 4,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(5): adjacency_offsets
                        wgpu::BindGroupLayoutEntry {
                            binding: 5,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(6): scratch (combined read-write arena —
                        // current/next invariants, fingerprints, dead_atoms,
                        // current/next/seen neighborhoods, seen_count,
                        // neighbor_scratch; see morgan.wgsl)
                        wgpu::BindGroupLayoutEntry {
                            binding: 6,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(7): iter_params (uniform)
                        wgpu::BindGroupLayoutEntry {
                            binding: 7,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Uniform,
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                    ],
                });

        let pipeline_layout = ctx
            .device
            .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("morgan_batch_pipeline_layout"),
                bind_group_layouts: &[&bind_group_layout],
                push_constant_ranges: &[],
            });

        let init_pipeline = ctx
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("morgan_batch_init"),
                layout: Some(&pipeline_layout),
                module: &shader,
                entry_point: Some("init_invariants_batch"),
                cache: None,
                compilation_options: Default::default(),
            });

        let iteration_pipeline =
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some("morgan_batch_iter"),
                    layout: Some(&pipeline_layout),
                    module: &shader,
                    entry_point: Some("morgan_iteration_batch"),
                    cache: None,
                    compilation_options: Default::default(),
                });

        let dedup_pipeline = ctx
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("morgan_batch_dedup"),
                layout: Some(&pipeline_layout),
                module: &shader,
                entry_point: Some("dedup_environments_batch"),
                cache: None,
                compilation_options: Default::default(),
            });

        let copy_pipeline = ctx
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("morgan_batch_copy"),
                layout: Some(&pipeline_layout),
                module: &shader,
                entry_point: Some("copy_invariants"),
                cache: None,
                compilation_options: Default::default(),
            });

        let copy_neighborhoods_pipeline =
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some("morgan_batch_copy_neighborhoods"),
                    layout: Some(&pipeline_layout),
                    module: &shader,
                    entry_point: Some("copy_neighborhoods"),
                    cache: None,
                    compilation_options: Default::default(),
                });

        Ok(Self {
            ctx,
            init_pipeline,
            iteration_pipeline,
            dedup_pipeline,
            copy_pipeline,
            copy_neighborhoods_pipeline,
            bind_group_layout,
        })
    }

    /// Generate fingerprints for multiple molecules in a single batch
    ///
    /// This is MUCH faster than calling generate_fingerprint() in a loop because:
    /// 1. Single buffer allocation for all molecules
    /// 2. Single GPU dispatch processes all atoms in parallel
    /// 3. Single sync point at the end
    pub fn generate_fingerprints_batch(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
        use_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
    ) -> Result<Vec<Vec<u32>>, GpuError> {
        pollster::block_on(self.generate_fingerprints_batch_async(
            molecules,
            radius,
            fp_size,
            use_chirality,
            use_bond_types,
            only_nonzero_invariants,
        ))
    }

    /// Async twin of [`Self::generate_fingerprints_batch`] for callers that
    /// can't block the current thread while waiting on the GPU (namely
    /// wasm32, where blocking the browser's single JS thread would deadlock
    /// — nothing could drive the readback future forward).
    pub async fn generate_fingerprints_batch_async(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
        use_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
    ) -> Result<Vec<Vec<u32>>, GpuError> {
        self.generate_fingerprints_batch_chunked(
            molecules,
            radius,
            fp_size,
            use_chirality,
            use_bond_types,
            only_nonzero_invariants,
            DISPATCH_LIMITS,
        )
        .await
    }

    /// Splits `molecules` into chunks that each stay under every cap in
    /// `limits` simultaneously (molecule count, total atoms, and total
    /// neighborhood words — see #27 and #33) and runs each chunk through its
    /// own GPU dispatch. Each molecule's GPU state (buffers, dedup
    /// bookkeeping) is scoped to its own chunk, so chunking never changes a
    /// molecule's resulting fingerprint — only how many molecules one GPU
    /// dispatch processes at a time.
    ///
    /// `limits` is a parameter (rather than always using `DISPATCH_LIMITS`)
    /// purely so tests can shrink individual caps to exercise each chunking
    /// trigger without needing tens of millions of atoms/molecules.
    #[allow(clippy::too_many_arguments)]
    async fn generate_fingerprints_batch_chunked(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
        use_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
        limits: DispatchLimits,
    ) -> Result<Vec<Vec<u32>>, GpuError> {
        let mut result = Vec::with_capacity(molecules.len());
        let mut start = 0usize;
        while start < molecules.len() {
            // Greedily grow [start, end) while it stays under every cap.
            // `end > start` guards the first molecule in: a single molecule
            // whose own atom/neighborhood-word cost already exceeds a cap
            // still gets its own (unsplittable) chunk, which
            // single_dispatch's bounds check will then reject.
            let mut end = start;
            let mut atoms_acc = 0u32;
            let mut neighborhood_words_acc = 0u32;
            while end < molecules.len() {
                let (atoms, neighborhood_words) = molecule_dispatch_cost(&molecules[end]);
                let next_atoms = atoms_acc + atoms;
                let next_neighborhood_words = neighborhood_words_acc + neighborhood_words;
                let next_count = end - start + 1;
                if end > start
                    && (next_count > limits.max_molecules
                        || next_atoms > limits.max_atoms
                        || next_neighborhood_words > limits.max_neighborhood_words)
                {
                    break;
                }
                atoms_acc = next_atoms;
                neighborhood_words_acc = next_neighborhood_words;
                end += 1;
            }

            result.extend(
                self.generate_fingerprints_batch_single_dispatch(
                    &molecules[start..end],
                    radius,
                    fp_size,
                    use_chirality,
                    use_bond_types,
                    only_nonzero_invariants,
                )
                .await?,
            );
            start = end;
        }
        Ok(result)
    }

    async fn generate_fingerprints_batch_single_dispatch(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
        use_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
    ) -> Result<Vec<Vec<u32>>, GpuError> {
        if molecules.is_empty() {
            return Ok(Vec::new());
        }

        let fp_words = fp_size.div_ceil(32) as usize;
        let num_molecules = molecules.len();

        if num_molecules > MAX_MOLECULES_PER_DISPATCH {
            return Err(GpuError::OperationFailed(format!(
                "batch of {num_molecules} molecules exceeds the maximum supported by a \
                 single GPU dispatch ({MAX_MOLECULES_PER_DISPATCH}); this should have been \
                 chunked by generate_fingerprints_batch"
            )));
        }

        // Encode all molecules and compute offsets
        let mut all_atoms: Vec<GpuAtom> = Vec::new();
        let mut all_bonds: Vec<GpuBond> = Vec::new();
        let mut all_adjacency: Vec<GpuNeighbor> = Vec::new();
        let mut all_adjacency_offsets: Vec<u32> = Vec::new();
        let mut molecule_offsets: Vec<MoleculeOffset> = Vec::new();

        let mut total_atoms = 0u32;
        let mut total_bonds = 0u32;
        let mut total_adj = 0u32;
        let mut total_neighborhood_words = 0u32;
        let mut total_seen_words = 0u32;

        for (mol_idx, mol) in molecules.iter().enumerate() {
            let num_atoms = mol.num_atoms();
            let num_bonds = mol.num_bonds();
            let bond_words = (num_bonds as u32).div_ceil(32).max(1);
            let neighborhood_words_for_mol = num_atoms as u32 * bond_words;
            // Capacity for the cumulative "seen neighborhoods" set: at most
            // atom_count new environments can be accepted per layer, across
            // `radius` layers (round 0 has no redundancy check, see the
            // shader's init_invariants_batch).
            let seen_words_for_mol = num_atoms as u32 * radius * bond_words;

            let offset = MoleculeOffset {
                atom_start: total_atoms,
                atom_count: num_atoms as u32,
                bond_start: total_bonds,
                bond_count: num_bonds as u32,
                adj_start: total_adj,
                fp_start: (mol_idx * fp_words) as u32,
                neighborhood_start: total_neighborhood_words,
                seen_start: total_seen_words,
            };
            molecule_offsets.push(offset);
            total_neighborhood_words += neighborhood_words_for_mol;
            total_seen_words += seen_words_for_mol;

            // Encode atoms
            for (idx, atom) in mol.atoms().iter().enumerate() {
                all_atoms.push(GpuAtom {
                    atomic_number: atom.atomic_number() as u32,
                    degree: mol.degree(idx) as u32,
                    is_aromatic: atom.is_aromatic() as u32,
                    chirality: match atom.chirality() {
                        chemcore::atom::Chirality::None => 0,
                        chemcore::atom::Chirality::Unspecified => 1,
                        chemcore::atom::Chirality::CounterClockwise => 2,
                        chemcore::atom::Chirality::Clockwise => 3,
                    },
                    formal_charge: atom.formal_charge() as i32,
                    _padding1: 0,
                    _padding2: 0,
                    _padding3: 0,
                });
            }

            // Encode bonds
            for bond in mol.bonds() {
                all_bonds.push(GpuBond {
                    atom1: bond.atom1() as u32,
                    atom2: bond.atom2() as u32,
                    order: match bond.order() {
                        chemcore::bond::BondOrder::Single => 1,
                        chemcore::bond::BondOrder::Double => 2,
                        chemcore::bond::BondOrder::Triple => 3,
                        chemcore::bond::BondOrder::Quadruple => 4,
                        chemcore::bond::BondOrder::Aromatic => 12,
                    },
                    stereo: match bond.stereo() {
                        chemcore::bond::BondStereo::None => 0,
                        chemcore::bond::BondStereo::E => 1,
                        chemcore::bond::BondStereo::Z => 2,
                        chemcore::bond::BondStereo::Unspecified => 0,
                    },
                    is_aromatic: bond.is_aromatic() as u32,
                    _padding1: 0,
                    _padding2: 0,
                    _padding3: 0,
                });
            }

            // Build adjacency list for this molecule
            for atom_idx in 0..num_atoms {
                all_adjacency_offsets.push(all_adjacency.len() as u32);
                for neighbor in mol.neighbors(atom_idx) {
                    all_adjacency.push(GpuNeighbor {
                        atom_idx: neighbor.atom_idx as u32,
                        bond_idx: neighbor.bond_idx as u32,
                    });
                }
            }
            // Final offset for last atom's end
            all_adjacency_offsets.push(all_adjacency.len() as u32);

            total_atoms += num_atoms as u32;
            total_bonds += num_bonds as u32;
            total_adj = all_adjacency.len() as u32;
        }

        // Defensive fallback: generate_fingerprints_batch_chunked always
        // keeps every chunk it builds under these same caps, so this should
        // be unreachable via the public API — except for a single molecule
        // whose own atom/neighborhood-word count already exceeds a cap,
        // which cannot be split any further (see #33).
        if total_atoms > MAX_ATOMS_PER_DISPATCH {
            return Err(GpuError::OperationFailed(format!(
                "batch has {total_atoms} total atoms, exceeding the maximum supported by a \
                 single GPU dispatch ({MAX_ATOMS_PER_DISPATCH}); if this is a single molecule, \
                 it cannot be chunked any further"
            )));
        }
        if total_neighborhood_words > MAX_NEIGHBORHOOD_WORDS_PER_DISPATCH {
            return Err(GpuError::OperationFailed(format!(
                "batch has {total_neighborhood_words} total neighborhood words, exceeding the \
                 maximum supported by a single GPU dispatch \
                 ({MAX_NEIGHBORHOOD_WORDS_PER_DISPATCH}); if this is a single molecule, it \
                 cannot be chunked any further"
            )));
        }

        // Remove the extra adjacency offset we added per molecule
        // Actually, we need to fix the adjacency_offsets to be contiguous
        // Let me recalculate...
        all_adjacency_offsets.clear();
        let mut adj_idx = 0u32;
        for mol in molecules.iter() {
            for atom_idx in 0..mol.num_atoms() {
                all_adjacency_offsets.push(adj_idx);
                adj_idx += mol.neighbors(atom_idx).len() as u32;
            }
        }
        all_adjacency_offsets.push(adj_idx); // Final terminator

        let params = BatchParams {
            num_molecules: num_molecules as u32,
            max_atoms: total_atoms,
            radius,
            fp_size,
            fp_words: fp_words as u32,
            use_chirality: use_chirality as u32,
            use_bond_types: use_bond_types as u32,
            only_nonzero_invariants: only_nonzero_invariants as u32,
            total_neighborhood_words,
            total_seen_words,
        };

        let manager = BufferManager::new(&self.ctx.device, &self.ctx.queue);

        // Pad vectors if empty (wgpu requires non-zero buffer sizes)
        if all_atoms.is_empty() {
            all_atoms.push(GpuAtom {
                atomic_number: 0,
                degree: 0,
                is_aromatic: 0,
                chirality: 0,
                formal_charge: 0,
                _padding1: 0,
                _padding2: 0,
                _padding3: 0,
            });
        }
        if all_bonds.is_empty() {
            all_bonds.push(GpuBond {
                atom1: 0,
                atom2: 0,
                order: 0,
                stereo: 0,
                is_aromatic: 0,
                _padding1: 0,
                _padding2: 0,
                _padding3: 0,
            });
        }
        if all_adjacency.is_empty() {
            all_adjacency.push(GpuNeighbor {
                atom_idx: 0,
                bond_idx: 0,
            });
        }

        // Create all buffers ONCE
        let params_buffer = manager.create_uniform_buffer("batch_params", 40);
        let offsets_buffer = self.create_buffer(&manager, &molecule_offsets);
        let atoms_buffer = self.create_buffer(&manager, &all_atoms);
        let bonds_buffer = self.create_buffer(&manager, &all_bonds);
        let adjacency_buffer = self.create_buffer(&manager, &all_adjacency);
        let adjacency_offsets_buffer = self.create_buffer(&manager, &all_adjacency_offsets);

        let fp_total_words = num_molecules * fp_words;

        // All read-write working state (previously 9 separate buffers: the
        // current/next invariants, fingerprints, dead_atoms,
        // current/next/seen neighborhoods, seen_count, and neighbor_scratch)
        // is packed into one combined `scratch` buffer, at the word offsets
        // below — see morgan.wgsl's matching `*_base()` functions, which
        // must agree with this exactly. current_invariants starts at word 0
        // and next_invariants at `total_atoms`; only fingerprints_offset
        // onward is needed on the Rust side (for the final readback below).
        let fingerprints_offset = 2 * total_atoms;
        let dead_atoms_offset = fingerprints_offset + fp_total_words as u32;
        let current_neighborhoods_offset = dead_atoms_offset + total_atoms;
        let next_neighborhoods_offset = current_neighborhoods_offset + total_neighborhood_words;
        let seen_neighborhoods_offset = next_neighborhoods_offset + total_neighborhood_words;
        let seen_count_offset = seen_neighborhoods_offset + total_seen_words;
        let neighbor_scratch_offset = seen_count_offset + num_molecules as u32;
        // neighbor_scratch holds one (bond_inv, neighbor_inv) pair — 2 words
        // — per adjacency entry, indexed identically to `adjacency` itself
        // (see #12), hence the *2.
        let scratch_total_words = neighbor_scratch_offset + total_adj * 2;

        let scratch_buffer = manager.create_storage_buffer(
            "scratch",
            std::cmp::max((scratch_total_words as usize * 4) as u64, 4),
            true,
        );

        let iter_params_buffer = manager.create_uniform_buffer("iter_params", 16);

        // Zero the whole combined arena in one write — every region needs a
        // zeroed starting state except seen_neighborhoods (only ever read up
        // to seen_count, which itself starts at 0, so its initial content
        // never matters), which zeroing doesn't hurt either.
        manager.write_buffer(&scratch_buffer, &vec![0u32; scratch_total_words as usize]);
        manager.write_buffer(&params_buffer, &[params]);

        // Create bind group
        let bind_group = self
            .ctx
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("morgan_batch_bind_group"),
                layout: &self.bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: params_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: offsets_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: atoms_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: bonds_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 4,
                        resource: adjacency_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 5,
                        resource: adjacency_offsets_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 6,
                        resource: scratch_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 7,
                        resource: iter_params_buffer.as_entire_binding(),
                    },
                ],
            });

        // Create command encoder
        let mut encoder = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("morgan_batch_encoder"),
            });

        let workgroups = total_atoms.div_ceil(256);

        // Pass 1: Initialize invariants for ALL molecules
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("morgan_batch_init_pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.init_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(workgroups, 1, 1);
        }

        // Submit the init pass now: queue.write_buffer calls are ordered
        // relative to queue.submit calls, not to when compute passes are
        // *recorded* into an encoder. Below, the per-layer `iter_params`
        // buffer is rewritten once per layer via write_buffer — if every
        // layer's passes were recorded into one shared encoder and submitted
        // only once at the very end (as this used to do), every layer's
        // iteration pass would see whichever `current_layer` was written
        // *last*, not its own. Submitting once per layer keeps each write
        // correctly ordered before the pass that depends on it.
        self.ctx.queue.submit(Some(encoder.finish()));

        // Pass 2..N: Iterate through radius layers for ALL molecules
        for layer in 0..radius {
            let iter_params = IterationParams {
                current_layer: layer,
                _padding1: 0,
                _padding2: 0,
                _padding3: 0,
            };
            manager.write_buffer(&iter_params_buffer, &[iter_params]);

            let mut layer_encoder =
                self.ctx
                    .device
                    .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                        label: Some(&format!("morgan_batch_layer_{}_encoder", layer)),
                    });

            {
                let mut pass = layer_encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some(&format!("morgan_batch_iter_layer_{}", layer)),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&self.iteration_pipeline);
                pass.set_bind_group(0, &bind_group, &[]);
                pass.dispatch_workgroups(workgroups, 1, 1);
            }

            // Redundant-environment dedup: decides which of this round's
            // freshly computed environments are unique (contribute a
            // fingerprint bit, get recorded as "seen") vs. redundant
            // (freeze that atom via dead_atoms, no bit). Must run after the
            // iteration pass (needs this round's next_invariants/
            // next_neighborhoods) and before the copy passes (which promote
            // next -> current for the following round).
            {
                let mut pass = layer_encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some(&format!("morgan_batch_dedup_layer_{}", layer)),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&self.dedup_pipeline);
                pass.set_bind_group(0, &bind_group, &[]);
                pass.dispatch_workgroups(
                    (num_molecules as u32).div_ceil(WORKGROUP_SIZE as u32),
                    1,
                    1,
                );
            }

            // Copy next -> current using GPU (no host sync needed beyond the
            // per-layer submit above)
            {
                let mut pass = layer_encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some(&format!("morgan_batch_copy_layer_{}", layer)),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&self.copy_pipeline);
                pass.set_bind_group(0, &bind_group, &[]);
                pass.dispatch_workgroups(workgroups, 1, 1);
            }

            {
                let neighborhood_workgroups = total_neighborhood_words.div_ceil(256).max(1);
                let mut pass = layer_encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some(&format!("morgan_batch_copy_neighborhoods_layer_{}", layer)),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&self.copy_neighborhoods_pipeline);
                pass.set_bind_group(0, &bind_group, &[]);
                pass.dispatch_workgroups(neighborhood_workgroups, 1, 1);
            }

            self.ctx.queue.submit(Some(layer_encoder.finish()));
        }

        // Single readback of ALL fingerprints, from their region within the
        // combined scratch buffer.
        let all_fps: Vec<u32> = manager
            .read_buffer_at(
                &scratch_buffer,
                (fingerprints_offset as u64) * 4,
                fp_total_words,
            )
            .await?;

        // Split into individual fingerprints
        let result: Vec<Vec<u32>> = all_fps
            .chunks(fp_words)
            .map(|chunk| chunk.to_vec())
            .collect();

        Ok(result)
    }

    fn create_buffer<T: bytemuck::Pod>(&self, manager: &BufferManager, data: &[T]) -> wgpu::Buffer {
        let size = std::mem::size_of_val(data) as u64;
        let buffer = manager.create_storage_buffer("temp", size, false);
        manager.write_buffer(&buffer, data);
        buffer
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chemcore::atom::{Atom, Element};
    use chemcore::bond::{Bond, BondOrder};

    /// A single `GpuMorganFingerprint`, lazily built once from the crate's
    /// shared test `GpuContext` and reused by every test below.
    ///
    /// `GpuMorganFingerprint::new()` builds its own `GpuContext`, so three tests
    /// calling it raced to create three devices against one GPU — the hang #19
    /// diagnosed and fixed for `context.rs`, `buffers.rs`, `pipeline.rs` and
    /// `tanimoto.rs`, which left this file out. It is what still hung
    /// `cargo test --workspace` after #134 fixed the same mistake in `chem-app`.
    ///
    /// Cloned per test rather than lent: `wgpu`'s device, queue and pipelines are
    /// refcounted handles, so a clone shares the device instead of requesting
    /// another.
    fn shared_test_morgan() -> Option<GpuMorganFingerprint> {
        static GPU: std::sync::OnceLock<Option<GpuMorganFingerprint>> = std::sync::OnceLock::new();
        GPU.get_or_init(|| {
            let ctx = crate::context::shared_test_context()?.clone();
            match GpuMorganFingerprint::from_context(ctx) {
                Ok(gpu) => Some(gpu),
                Err(e) => {
                    println!("GPU morgan pipeline unavailable, skipping: {}", e);
                    None
                }
            }
        })
        .clone()
    }

    fn chain_molecule(num_atoms: u32) -> Molecule {
        let mut mol = Molecule::new();
        let first = mol.add_atom(Atom::new(Element::carbon()));
        let mut prev = first;
        for _ in 1..num_atoms {
            let next = mol.add_atom(Atom::new(Element::carbon()));
            mol.add_bond(Bond::new(prev, next, BondOrder::Single))
                .unwrap();
            prev = next;
        }
        mol.calculate_implicit_hydrogens();
        mol
    }

    fn no_bond_molecule(num_atoms: u32) -> Molecule {
        let mut mol = Molecule::new();
        for _ in 0..num_atoms {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        mol.calculate_implicit_hydrogens();
        mol
    }

    /// Regression test for #27's chunking follow-up: molecules must produce
    /// the same fingerprints whether they land in one GPU dispatch or are
    /// split across several, since chunking is purely a workaround for
    /// WebGPU's dispatch-group cap and must not change results. Uses a tiny
    /// artificial `max_molecules` so the multi-chunk path is exercised
    /// without needing millions of molecules.
    #[test]
    fn test_chunked_dispatch_matches_single_dispatch() {
        let Some(gpu) = shared_test_morgan() else {
            return;
        };

        // Varying atom/bond counts per molecule so chunk boundaries can't
        // accidentally line up with identical, trivially-matching molecules.
        let molecules: Vec<Molecule> = (1..=10u32).map(chain_molecule).collect();

        let radius = 2;
        let fp_size = 128;

        let single = pollster::block_on(gpu.generate_fingerprints_batch_chunked(
            &molecules,
            radius,
            fp_size,
            false,
            true,
            false,
            DispatchLimits {
                max_molecules: molecules.len(),
                max_atoms: u32::MAX,
                max_neighborhood_words: u32::MAX,
            },
        ))
        .unwrap();
        let chunked = pollster::block_on(gpu.generate_fingerprints_batch_chunked(
            &molecules,
            radius,
            fp_size,
            false,
            true,
            false,
            DispatchLimits {
                max_molecules: 3,
                max_atoms: u32::MAX,
                max_neighborhood_words: u32::MAX,
            },
        ))
        .unwrap();

        assert_eq!(chunked, single);
    }

    /// Regression test for #33: even when a batch's molecule count is well
    /// under `MAX_MOLECULES_PER_DISPATCH`, a small `max_atoms` cap (standing
    /// in for the real `256 * 65,535`, which would need tens of millions of
    /// atoms to hit) must still force chunking, since `init`/`iteration`/
    /// `copy` dispatch `ceil(total_atoms / 256)` workgroups.
    #[test]
    fn test_chunking_triggered_by_atom_count() {
        let Some(gpu) = shared_test_morgan() else {
            return;
        };

        // No bonds, so neighborhood_words == atom_count for each molecule —
        // isolates the atom-count cap from the neighborhood-word cap.
        let molecules: Vec<Molecule> = [5u32, 7, 9, 11].into_iter().map(no_bond_molecule).collect();

        let radius = 2;
        let fp_size = 128;
        let generous = DispatchLimits {
            max_molecules: usize::MAX,
            max_atoms: u32::MAX,
            max_neighborhood_words: u32::MAX,
        };
        let atom_capped = DispatchLimits {
            max_molecules: usize::MAX,
            max_atoms: 15,
            max_neighborhood_words: u32::MAX,
        };

        let single = pollster::block_on(gpu.generate_fingerprints_batch_chunked(
            &molecules, radius, fp_size, false, true, false, generous,
        ))
        .unwrap();
        let chunked = pollster::block_on(gpu.generate_fingerprints_batch_chunked(
            &molecules,
            radius,
            fp_size,
            false,
            true,
            false,
            atom_capped,
        ))
        .unwrap();

        assert_eq!(chunked, single);
    }

    /// Regression test for #33: a small `max_neighborhood_words` cap must
    /// force chunking even when molecule count and atom count are both well
    /// under their own caps, since `copy_neighborhoods` dispatches
    /// `ceil(total_neighborhood_words / 256)` workgroups. Uses chain
    /// molecules with >32 bonds each (bond_words = 2) so
    /// `neighborhood_words = atoms * 2` differs from the atom count, proving
    /// this cap is checked independently rather than piggybacking on the
    /// atom-count check.
    #[test]
    fn test_chunking_triggered_by_neighborhood_word_count() {
        let Some(gpu) = shared_test_morgan() else {
            return;
        };

        let molecules: Vec<Molecule> = [35u32, 40, 45, 50]
            .into_iter()
            .map(chain_molecule)
            .collect();

        let radius = 1;
        let fp_size = 128;
        let generous = DispatchLimits {
            max_molecules: molecules.len(),
            max_atoms: u32::MAX,
            max_neighborhood_words: u32::MAX,
        };
        let neighborhood_capped = DispatchLimits {
            max_molecules: usize::MAX,
            max_atoms: u32::MAX,
            max_neighborhood_words: 150,
        };

        let single = pollster::block_on(gpu.generate_fingerprints_batch_chunked(
            &molecules, radius, fp_size, false, true, false, generous,
        ))
        .unwrap();
        let chunked = pollster::block_on(gpu.generate_fingerprints_batch_chunked(
            &molecules,
            radius,
            fp_size,
            false,
            true,
            false,
            neighborhood_capped,
        ))
        .unwrap();

        assert_eq!(chunked, single);
    }
}
