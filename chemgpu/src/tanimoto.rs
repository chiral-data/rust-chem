use crate::buffers::BufferManager;
use crate::context::GpuContext;
use crate::error::GpuError;
use bytemuck::{Pod, Zeroable};
use wgpu;

// ============================================================================
// GPU DATA STRUCTURES
// ============================================================================

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct TanimotoParams {
    num_queries: u32,
    num_targets: u32,
    fp_words: u32,
    _padding: u32,
}

// ============================================================================
// GPU TANIMOTO CALCULATOR
// ============================================================================

pub struct GpuTanimoto {
    ctx: GpuContext,
    all_pairs_pipeline: wgpu::ComputePipeline,
    single_query_pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}

impl GpuTanimoto {
    pub fn new() -> Result<Self, GpuError> {
        let ctx = GpuContext::new()?;
        let shader_source = include_str!("shaders/tanimoto.wgsl");

        let shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("tanimoto_shader"),
                source: wgpu::ShaderSource::Wgsl(shader_source.into()),
            });

        // Create bind group layout
        let bind_group_layout =
            ctx.device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("tanimoto_bind_group_layout"),
                    entries: &[
                        // @binding(0): queries
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
                        // @binding(1): targets
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
                        // @binding(2): results
                        wgpu::BindGroupLayoutEntry {
                            binding: 2,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // @binding(3): params
                        wgpu::BindGroupLayoutEntry {
                            binding: 3,
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
                label: Some("tanimoto_pipeline_layout"),
                bind_group_layouts: &[&bind_group_layout],
                push_constant_ranges: &[],
            });

        let all_pairs_pipeline =
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some("tanimoto_all_pairs"),
                    layout: Some(&pipeline_layout),
                    module: &shader,
                    entry_point: Some("compute_tanimoto"),
                    cache: None,
                    compilation_options: Default::default(),
                });

        let single_query_pipeline =
            ctx.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some("tanimoto_single_query"),
                    layout: Some(&pipeline_layout),
                    module: &shader,
                    entry_point: Some("compute_tanimoto_single_query"),
                    cache: None,
                    compilation_options: Default::default(),
                });

        Ok(Self {
            ctx,
            all_pairs_pipeline,
            single_query_pipeline,
            bind_group_layout,
        })
    }

    /// Compute Tanimoto similarity between one query and multiple targets
    ///
    /// # Arguments
    /// * `query` - Single fingerprint (as u32 words)
    /// * `targets` - Multiple fingerprints (flattened)
    /// * `fp_size` - Fingerprint size in bits
    ///
    /// # Returns
    /// Vector of similarities, one per target
    pub fn compute_single_query(
        &self,
        query: &[u32],
        targets: &[u32],
    ) -> Result<Vec<f32>, GpuError> {
        let fp_words = query.len();
        if fp_words == 0 {
            return Err(GpuError::BufferError(
                "Query fingerprint is empty".to_string(),
            ));
        }

        let num_targets = targets.len() / fp_words;

        if !targets.len().is_multiple_of(fp_words) {
            return Err(GpuError::BufferError(
                "Target fingerprints not aligned".to_string(),
            ));
        }

        // No targets means no results — return early rather than creating
        // 0-byte GPU buffers, which wgpu rejects as an invalid buffer size.
        if num_targets == 0 {
            return Ok(Vec::new());
        }

        let params = TanimotoParams {
            num_queries: 1,
            num_targets: num_targets as u32,
            fp_words: fp_words as u32,
            _padding: 0,
        };

        let manager = BufferManager::new(&self.ctx.device, &self.ctx.queue);

        // Create buffers
        let query_buffer = self.create_buffer(&manager, query);
        let target_buffer = self.create_buffer(&manager, targets);
        let params_buffer = manager.create_uniform_buffer("params", 16);
        let results_buffer =
            manager.create_storage_buffer("results", (num_targets * 4) as u64, true);

        manager.write_buffer(&params_buffer, &[params]);

        // Create bind group
        let bind_group = self
            .ctx
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("tanimoto_bind_group"),
                layout: &self.bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: query_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: target_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: results_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: params_buffer.as_entire_binding(),
                    },
                ],
            });

        // Execute
        let mut encoder = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("tanimoto_encoder"),
            });

        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("tanimoto_pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.single_query_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups((num_targets as u32).div_ceil(256), 1, 1);
        }

        self.ctx.queue.submit(Some(encoder.finish()));

        // Read results
        manager.read_buffer_blocking(&results_buffer, num_targets)
    }

    /// Compute all-pairs Tanimoto similarity matrix
    ///
    /// # Arguments
    /// * `queries` - Query fingerprints (flattened)
    /// * `targets` - Target fingerprints (flattened)
    /// * `fp_size` - Fingerprint size in bits
    ///
    /// # Returns
    /// Similarity matrix [num_queries × num_targets] in row-major order
    pub fn compute_all_pairs(
        &self,
        queries: &[u32],
        targets: &[u32],
        fp_size: u32,
    ) -> Result<Vec<f32>, GpuError> {
        let fp_words = fp_size.div_ceil(32) as usize;
        if fp_words == 0 {
            return Err(GpuError::BufferError("fp_size must be > 0".to_string()));
        }

        let num_queries = queries.len() / fp_words;
        let num_targets = targets.len() / fp_words;

        if !queries.len().is_multiple_of(fp_words) || !targets.len().is_multiple_of(fp_words) {
            return Err(GpuError::BufferError(
                "Fingerprints not aligned".to_string(),
            ));
        }

        // No queries or no targets means no results — return early rather
        // than creating 0-byte GPU buffers, which wgpu rejects as an invalid
        // buffer size.
        if num_queries == 0 || num_targets == 0 {
            return Ok(Vec::new());
        }

        let params = TanimotoParams {
            num_queries: num_queries as u32,
            num_targets: num_targets as u32,
            fp_words: fp_words as u32,
            _padding: 0,
        };

        let manager = BufferManager::new(&self.ctx.device, &self.ctx.queue);

        let query_buffer = self.create_buffer(&manager, queries);
        let target_buffer = self.create_buffer(&manager, targets);
        let params_buffer = manager.create_uniform_buffer("params", 16);
        let results_buffer =
            manager.create_storage_buffer("results", (num_queries * num_targets * 4) as u64, true);

        manager.write_buffer(&params_buffer, &[params]);

        let bind_group = self
            .ctx
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("tanimoto_bind_group"),
                layout: &self.bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: query_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: target_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: results_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: params_buffer.as_entire_binding(),
                    },
                ],
            });

        let mut encoder = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("tanimoto_encoder"),
            });

        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("tanimoto_all_pairs_pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.all_pairs_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);

            // 16x16 workgroup size
            let workgroups_x = (num_queries as u32).div_ceil(16);
            let workgroups_y = (num_targets as u32).div_ceil(16);
            pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
        }

        self.ctx.queue.submit(Some(encoder.finish()));

        manager.read_buffer_blocking(&results_buffer, num_queries * num_targets)
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

    #[test]
    fn test_compute_single_query_empty_targets_returns_empty_not_panic() {
        env_logger::try_init().ok();

        let Ok(gpu) = GpuTanimoto::new() else {
            println!("GPU not available, skipping");
            return;
        };

        let query = vec![0u32; 64];
        let result = gpu.compute_single_query(&query, &[]).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_compute_single_query_empty_query_is_error() {
        env_logger::try_init().ok();

        let Ok(gpu) = GpuTanimoto::new() else {
            println!("GPU not available, skipping");
            return;
        };

        let targets = vec![0u32; 64];
        assert!(gpu.compute_single_query(&[], &targets).is_err());
    }

    #[test]
    fn test_compute_all_pairs_empty_inputs_return_empty_not_panic() {
        env_logger::try_init().ok();

        let Ok(gpu) = GpuTanimoto::new() else {
            println!("GPU not available, skipping");
            return;
        };

        let fp_size = 2048u32;
        let one_fp = vec![0u32; (fp_size / 32) as usize];

        assert!(gpu.compute_all_pairs(&[], &one_fp, fp_size).unwrap().is_empty());
        assert!(gpu.compute_all_pairs(&one_fp, &[], fp_size).unwrap().is_empty());
        assert!(gpu.compute_all_pairs(&[], &[], fp_size).unwrap().is_empty());
    }

    #[test]
    fn test_compute_all_pairs_zero_fp_size_is_error() {
        env_logger::try_init().ok();

        let Ok(gpu) = GpuTanimoto::new() else {
            println!("GPU not available, skipping");
            return;
        };

        assert!(gpu.compute_all_pairs(&[1u32], &[1u32], 0).is_err());
    }
}
