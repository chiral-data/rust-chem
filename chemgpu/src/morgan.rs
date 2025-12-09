use crate::buffers::BufferManager;
use crate::context::GpuContext;
use crate::error::GpuError;
use crate::pipeline::ComputePipeline;
use bytemuck::{Pod, Zeroable};
use chemcore::molecule::Molecule;
use wgpu;

// ============================================================================
// GPU DATA STRUCTURES (must match WGSL layout)
// ============================================================================

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct GpuMoleculeData {
    num_atoms: u32,
    num_bonds: u32,
    radius: u32,
    fp_size: u32,
    use_chirality: u32,
    use_bond_types: u32,
    only_nonzero_invariants: u32,
    _padding: u32,
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
// ENCODER: Molecule → GPU format
// ============================================================================

pub struct MoleculeEncoder;

impl MoleculeEncoder {
    pub fn encode_molecule(
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
        use_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
    ) -> EncodedMolecule {
        let num_atoms = mol.num_atoms();
        let num_bonds = mol.num_bonds();

        // Encode metadata
        let metadata = GpuMoleculeData {
            num_atoms: num_atoms as u32,
            num_bonds: num_bonds as u32,
            radius,
            fp_size,
            use_chirality: use_chirality as u32,
            use_bond_types: use_bond_types as u32,
            only_nonzero_invariants: only_nonzero_invariants as u32,
            _padding: 0,
        };

        // Encode atoms
        let mut atoms = Vec::with_capacity(num_atoms);
        for (idx, atom) in mol.atoms().iter().enumerate() {
            atoms.push(GpuAtom {
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
        let mut bonds = Vec::with_capacity(num_bonds);
        for bond in mol.bonds() {
            bonds.push(GpuBond {
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

        // Build flat adjacency list
        let mut adjacency = Vec::new();
        let mut adjacency_offsets = vec![0u32];

        for atom_idx in 0..num_atoms {
            for neighbor in mol.neighbors(atom_idx) {
                adjacency.push(GpuNeighbor {
                    atom_idx: neighbor.atom_idx as u32,
                    bond_idx: neighbor.bond_idx as u32,
                });
            }
            adjacency_offsets.push(adjacency.len() as u32);
        }

        EncodedMolecule {
            metadata,
            atoms,
            bonds,
            adjacency,
            adjacency_offsets,
        }
    }
}

pub struct EncodedMolecule {
    pub metadata: GpuMoleculeData,
    pub atoms: Vec<GpuAtom>,
    pub bonds: Vec<GpuBond>,
    pub adjacency: Vec<GpuNeighbor>,
    pub adjacency_offsets: Vec<u32>,
}

// ============================================================================
// GPU MORGAN FINGERPRINT GENERATOR
// ============================================================================

pub struct GpuMorganFingerprint {
    ctx: GpuContext,
    init_pipeline: wgpu::ComputePipeline, // ComputePipeline,
    iteration_pipeline: wgpu::ComputePipeline, // ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}

impl GpuMorganFingerprint {
    pub fn new() -> Result<Self, GpuError> {
        let ctx = GpuContext::new()?;

        let shader_source = include_str!("shaders/morgan.wgsl");

        // Create shader module once
        let shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("morgan_shader"),
                source: wgpu::ShaderSource::Wgsl(shader_source.into()),
            });

        // Create bind group layout (8 bindings for Morgan)
        let bind_group_layout =
            ctx.device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("morgan_bind_group_layout"),
                    entries: &[
                        // @binding(0): molecule_data
                        wgpu::BindGroupLayoutEntry {
                            binding: 0,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(1): atoms
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
                        // @binding(2): bonds
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
                        // @binding(3): adjacency
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
                        // @binding(4): adjacency_offsets
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
                        // @binding(5): current_invariants (read_write)
                        wgpu::BindGroupLayoutEntry {
                            binding: 5,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(6): next_invariants (read_write)
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
                        // @binding(7): fingerprint (read_write)
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
                        wgpu::BindGroupLayoutEntry {
                            binding: 8,
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
                label: Some("morgan_pipeline_layout"),
                bind_group_layouts: &[&bind_group_layout],
                push_constant_ranges: &[],
            });

        let init_pipeline = ctx
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("morgan_init"),
                layout: Some(&pipeline_layout),
                module: &shader,
                entry_point: Some("init_invariants"),
                cache: None,
                compilation_options: Default::default(),
            });

        let iteration_pipeline =
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some("morgan_iter"),
                    layout: Some(&pipeline_layout),
                    module: &shader,
                    entry_point: Some("morgan_iteration"),
                    cache: None,
                    compilation_options: Default::default(),
                });

        //let init_pipeline =
        //    ComputePipeline::new(&ctx.device, shader_source, "init_invariants", "morgan_init")?;
        //
        //let iteration_pipeline = ComputePipeline::new(
        //    &ctx.device,
        //    shader_source,
        //    "morgan_iteration",
        //    "morgan_iter",
        //)?;

        Ok(Self {
            ctx,
            init_pipeline,
            iteration_pipeline,
            bind_group_layout,
        })
    }

    pub fn generate_fingerprint(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
        use_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
    ) -> Result<Vec<u32>, GpuError> {
        let encoded = MoleculeEncoder::encode_molecule(
            mol,
            radius,
            fp_size,
            use_chirality,
            use_bond_types,
            only_nonzero_invariants,
        );

        let num_atoms = encoded.atoms.len();
        let fp_words = ((fp_size + 31) / 32) as usize;

        let manager = BufferManager::new(&self.ctx.device, &self.ctx.queue);

        // Create buffers
        let metadata_buffer = self.create_buffer(&manager, &[encoded.metadata]);
        let atoms_buffer = self.create_buffer(&manager, &encoded.atoms);
        let bonds_buffer = self.create_buffer(&manager, &encoded.bonds);
        let adjacency_buffer = self.create_buffer(&manager, &encoded.adjacency);
        let adjacency_offsets_buffer = self.create_buffer(&manager, &encoded.adjacency_offsets);

        let current_inv_buffer =
            manager.create_storage_buffer("current_invariants", (num_atoms * 4) as u64, true);
        let next_inv_buffer =
            manager.create_storage_buffer("next_invariants", (num_atoms * 4) as u64, true);
        let fp_buffer = manager.create_storage_buffer("fingerprint", (fp_words * 4) as u64, true);
        let iter_params_buffer = manager.create_uniform_buffer("iter_params", 16);

        // Initialize fingerprint to zero
        manager.write_buffer(&fp_buffer, &vec![0u32; fp_words]);

        // Create bind group layout for Morgan pipeline
        let bind_group = self
            .ctx
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("morgan_bind_group"),
                layout: &self.bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: metadata_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: atoms_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: bonds_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: adjacency_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 4,
                        resource: adjacency_offsets_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 5,
                        resource: current_inv_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 6,
                        resource: next_inv_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 7,
                        resource: fp_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 8,
                        resource: iter_params_buffer.as_entire_binding(),
                    },
                ],
            });

        // Execute compute passes
        let mut encoder = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("morgan_encoder"),
            });

        // Pass 1: Initialize invariants (layer 0)
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("morgan_init_pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.init_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups((num_atoms as u32 + 63) / 64, 1, 1);
        }

        // Pass 2..N: Iterate through radius layers
        for layer in 0..radius {
            let params = IterationParams {
                current_layer: layer,
                _padding1: 0,
                _padding2: 0,
                _padding3: 0,
            };
            manager.write_buffer(&iter_params_buffer, &[params]);

            println!("GPU computing layer {}", layer);

            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some(&format!("morgan_iter_layer_{}", layer)),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.iteration_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups((num_atoms as u32 + 63) / 64, 1, 1);
            drop(pass);

            // Swap buffers (copy next → current)
            encoder.copy_buffer_to_buffer(
                &next_inv_buffer,
                0,
                &current_inv_buffer,
                0,
                (num_atoms * 4) as u64,
            );
        }

        self.ctx.queue.submit(Some(encoder.finish()));

        let gpu_invariants: Vec<u32> =
            manager.read_buffer_blocking(&current_inv_buffer, num_atoms)?;
        println!("\nGPU invariants after all layers:");
        for (i, &inv) in gpu_invariants.iter().enumerate() {
            let expected = if radius >= 2 {
                501318127
            } else if radius >= 1 {
                2837459962
            } else {
                37324
            };
            println!(
                "  Atom {}: GPU={}, CPU expected={}, match={}",
                i,
                inv,
                expected,
                inv == expected
            );
        }

        // Read back fingerprint
        let gpu_fp = manager.read_buffer_blocking(&fp_buffer, fp_words)?;

        // ADD THIS DEBUG:
        println!("\nGPU fingerprint bits set:");
        for (word_idx, &word) in gpu_fp.iter().enumerate() {
            if word != 0 {
                for bit in 0..32 {
                    if (word & (1 << bit)) != 0 {
                        let global_bit = word_idx * 32 + bit;
                        println!("  Bit {}: word {} bit {}", global_bit, word_idx, bit);
                    }
                }
            }
        }

        Ok(gpu_fp)
    }

    fn create_buffer<T: bytemuck::Pod>(&self, manager: &BufferManager, data: &[T]) -> wgpu::Buffer {
        let size = (data.len() * std::mem::size_of::<T>()) as u64;
        let buffer = manager.create_storage_buffer("temp", size, false);
        manager.write_buffer(&buffer, data);
        buffer
    }
}
