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
        let ctx = GpuContext::new()?;

        let shader_source = include_str!("shaders/morgan.wgsl");

        let shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("morgan_batch_shader"),
                source: wgpu::ShaderSource::Wgsl(shader_source.into()),
            });

        // Create bind group layout (10 bindings for batch Morgan)
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
                        // @binding(6): current_invariants
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
                        // @binding(7): next_invariants
                        wgpu::BindGroupLayoutEntry {
                            binding: 7,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(8): fingerprints
                        wgpu::BindGroupLayoutEntry {
                            binding: 8,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(9): iter_params (uniform)
                        wgpu::BindGroupLayoutEntry {
                            binding: 9,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Uniform,
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(10): dead_atoms
                        wgpu::BindGroupLayoutEntry {
                            binding: 10,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(11): current_neighborhoods
                        wgpu::BindGroupLayoutEntry {
                            binding: 11,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(12): next_neighborhoods
                        wgpu::BindGroupLayoutEntry {
                            binding: 12,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(13): seen_neighborhoods
                        wgpu::BindGroupLayoutEntry {
                            binding: 13,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(14): seen_count
                        wgpu::BindGroupLayoutEntry {
                            binding: 14,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(15): neighbor_scratch
                        wgpu::BindGroupLayoutEntry {
                            binding: 15,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
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
        self.generate_fingerprints_batch_chunked(
            molecules,
            radius,
            fp_size,
            use_chirality,
            use_bond_types,
            only_nonzero_invariants,
            DISPATCH_LIMITS,
        )
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
    fn generate_fingerprints_batch_chunked(
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

            result.extend(self.generate_fingerprints_batch_single_dispatch(
                &molecules[start..end],
                radius,
                fp_size,
                use_chirality,
                use_bond_types,
                only_nonzero_invariants,
            )?);
            start = end;
        }
        Ok(result)
    }

    fn generate_fingerprints_batch_single_dispatch(
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
        let params_buffer = manager.create_uniform_buffer("batch_params", 32);
        let offsets_buffer = self.create_buffer(&manager, &molecule_offsets);
        let atoms_buffer = self.create_buffer(&manager, &all_atoms);
        let bonds_buffer = self.create_buffer(&manager, &all_bonds);
        let adjacency_buffer = self.create_buffer(&manager, &all_adjacency);
        let adjacency_offsets_buffer = self.create_buffer(&manager, &all_adjacency_offsets);

        let invariants_size = (total_atoms as usize * 4) as u64;
        let current_inv_buffer = manager.create_storage_buffer(
            "current_invariants",
            std::cmp::max(invariants_size, 4),
            true,
        );
        let next_inv_buffer = manager.create_storage_buffer(
            "next_invariants",
            std::cmp::max(invariants_size, 4),
            true,
        );

        let fp_total_words = num_molecules * fp_words;
        let fp_buffer =
            manager.create_storage_buffer("fingerprints", (fp_total_words * 4) as u64, true);

        let iter_params_buffer = manager.create_uniform_buffer("iter_params", 16);

        // Redundant-environment dedup buffers (see chemgpu/src/shaders/morgan.wgsl
        // dedup_environments_batch for the full explanation).
        let dead_atoms_buffer = manager.create_storage_buffer(
            "dead_atoms",
            std::cmp::max((total_atoms as usize * 4) as u64, 4),
            true,
        );
        let neighborhood_size = std::cmp::max((total_neighborhood_words as usize * 4) as u64, 4);
        let current_neigh_buffer =
            manager.create_storage_buffer("current_neighborhoods", neighborhood_size, true);
        let next_neigh_buffer =
            manager.create_storage_buffer("next_neighborhoods", neighborhood_size, true);
        let seen_neigh_buffer = manager.create_storage_buffer(
            "seen_neighborhoods",
            std::cmp::max((total_seen_words as usize * 4) as u64, 4),
            true,
        );
        let seen_count_buffer = manager.create_storage_buffer(
            "seen_count",
            std::cmp::max((num_molecules * 4) as u64, 4),
            true,
        );
        // One vec2<u32> (bond_inv, neighbor_inv) slot per adjacency entry —
        // indexed identically to `adjacency` itself, so each atom's own
        // neighbor candidate list is unbounded (no fixed-size cap; see #12).
        let neighbor_scratch_buffer = manager.create_storage_buffer(
            "neighbor_scratch",
            std::cmp::max((total_adj as usize * 8) as u64, 8),
            true,
        );

        // Initialize fingerprints, invariants, and dedup state to zero
        manager.write_buffer(&fp_buffer, &vec![0u32; fp_total_words]);
        manager.write_buffer(&params_buffer, &[params]);
        manager.write_buffer(&dead_atoms_buffer, &vec![0u32; total_atoms as usize]);
        manager.write_buffer(
            &current_neigh_buffer,
            &vec![0u32; total_neighborhood_words as usize],
        );
        manager.write_buffer(
            &next_neigh_buffer,
            &vec![0u32; total_neighborhood_words as usize],
        );
        manager.write_buffer(&seen_count_buffer, &vec![0u32; num_molecules]);
        manager.write_buffer(
            &neighbor_scratch_buffer,
            &vec![0u32; total_adj as usize * 2],
        );

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
                        resource: current_inv_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 7,
                        resource: next_inv_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 8,
                        resource: fp_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 9,
                        resource: iter_params_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 10,
                        resource: dead_atoms_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 11,
                        resource: current_neigh_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 12,
                        resource: next_neigh_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 13,
                        resource: seen_neigh_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 14,
                        resource: seen_count_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 15,
                        resource: neighbor_scratch_buffer.as_entire_binding(),
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

        // Single readback of ALL fingerprints
        let all_fps: Vec<u32> = manager.read_buffer_blocking(&fp_buffer, fp_total_words)?;

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
        let Ok(gpu) = GpuMorganFingerprint::new() else {
            println!("GPU not available, skipping");
            return;
        };

        // Varying atom/bond counts per molecule so chunk boundaries can't
        // accidentally line up with identical, trivially-matching molecules.
        let molecules: Vec<Molecule> = (1..=10u32).map(chain_molecule).collect();

        let radius = 2;
        let fp_size = 128;

        let single = gpu
            .generate_fingerprints_batch_chunked(
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
            )
            .unwrap();
        let chunked = gpu
            .generate_fingerprints_batch_chunked(
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
            )
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
        let Ok(gpu) = GpuMorganFingerprint::new() else {
            println!("GPU not available, skipping");
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

        let single = gpu
            .generate_fingerprints_batch_chunked(
                &molecules, radius, fp_size, false, true, false, generous,
            )
            .unwrap();
        let chunked = gpu
            .generate_fingerprints_batch_chunked(
                &molecules,
                radius,
                fp_size,
                false,
                true,
                false,
                atom_capped,
            )
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
        let Ok(gpu) = GpuMorganFingerprint::new() else {
            println!("GPU not available, skipping");
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

        let single = gpu
            .generate_fingerprints_batch_chunked(
                &molecules, radius, fp_size, false, true, false, generous,
            )
            .unwrap();
        let chunked = gpu
            .generate_fingerprints_batch_chunked(
                &molecules,
                radius,
                fp_size,
                false,
                true,
                false,
                neighborhood_capped,
            )
            .unwrap();

        assert_eq!(chunked, single);
    }
}
