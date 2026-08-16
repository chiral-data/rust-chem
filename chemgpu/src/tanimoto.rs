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

/// A target fingerprint set already uploaded to the GPU.
///
/// Search callers upload this once per dataset via [`GpuTanimoto::upload_targets`]
/// and reuse it across many queries with [`GpuTanimoto::compute_single_query_against`],
/// instead of re-uploading the whole dataset on every query.
pub struct GpuTargetSet {
    buffer: wgpu::Buffer,
    count: usize,
    fp_words: usize,
}

impl GpuTargetSet {
    /// Number of target fingerprints held in this set.
    pub fn count(&self) -> usize {
        self.count
    }

    /// Words per fingerprint in this set.
    pub fn fp_words(&self) -> usize {
        self.fp_words
    }
}

impl GpuTanimoto {
    pub fn new() -> Result<Self, GpuError> {
        Self::from_context(GpuContext::new()?)
    }

    /// Build from an already-initialized [`GpuContext`] instead of creating a
    /// new one. Lets callers (notably this crate's own test suite) share a
    /// single context across many `GpuTanimoto` instances instead of each
    /// racing to create its own `wgpu::Instance`/adapter/device concurrently.
    pub(crate) fn from_context(ctx: GpuContext) -> Result<Self, GpuError> {
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
    ///
    /// Uploads `targets` to a fresh GPU buffer on every call. When the same
    /// target set is queried repeatedly (e.g. an interactive search box),
    /// prefer [`Self::upload_targets`] once plus
    /// [`Self::compute_single_query_against`] per query instead.
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

        if !targets.len().is_multiple_of(fp_words) {
            return Err(GpuError::BufferError(
                "Target fingerprints not aligned".to_string(),
            ));
        }

        // No targets means no results — return early rather than creating
        // 0-byte GPU buffers, which wgpu rejects as an invalid buffer size.
        if targets.is_empty() {
            return Ok(Vec::new());
        }

        let manager = BufferManager::new(&self.ctx.device, &self.ctx.queue);
        let target_buffer = self.create_buffer(&manager, targets);
        let target_set = GpuTargetSet {
            buffer: target_buffer,
            count: targets.len() / fp_words,
            fp_words,
        };

        self.compute_single_query_against(query, &target_set)
    }

    /// Upload a target fingerprint set to the GPU once, for reuse across
    /// many queries via [`Self::compute_single_query_against`].
    ///
    /// # Arguments
    /// * `targets` - Multiple fingerprints (flattened, as u32 words)
    /// * `fp_words` - Words per fingerprint
    pub fn upload_targets(
        &self,
        targets: &[u32],
        fp_words: usize,
    ) -> Result<GpuTargetSet, GpuError> {
        if fp_words == 0 {
            return Err(GpuError::BufferError("fp_words must be > 0".to_string()));
        }

        if !targets.len().is_multiple_of(fp_words) {
            return Err(GpuError::BufferError(
                "Target fingerprints not aligned".to_string(),
            ));
        }

        let manager = BufferManager::new(&self.ctx.device, &self.ctx.queue);
        let buffer = self.create_buffer(&manager, targets);

        Ok(GpuTargetSet {
            buffer,
            count: targets.len() / fp_words,
            fp_words,
        })
    }

    /// Compute Tanimoto similarity between one query and a previously
    /// uploaded target set, without re-uploading the targets.
    pub fn compute_single_query_against(
        &self,
        query: &[u32],
        targets: &GpuTargetSet,
    ) -> Result<Vec<f32>, GpuError> {
        let fp_words = query.len();
        if fp_words == 0 {
            return Err(GpuError::BufferError(
                "Query fingerprint is empty".to_string(),
            ));
        }
        if fp_words != targets.fp_words {
            return Err(GpuError::BufferError(format!(
                "Query fingerprint size ({} words) does not match uploaded target set ({} words)",
                fp_words, targets.fp_words
            )));
        }

        // No targets means no results — return early rather than creating
        // 0-byte GPU buffers, which wgpu rejects as an invalid buffer size.
        if targets.count == 0 {
            return Ok(Vec::new());
        }

        let params = TanimotoParams {
            num_queries: 1,
            num_targets: targets.count as u32,
            fp_words: fp_words as u32,
            _padding: 0,
        };

        let manager = BufferManager::new(&self.ctx.device, &self.ctx.queue);

        let query_buffer = self.create_buffer(&manager, query);
        let params_buffer = manager.create_uniform_buffer("params", 16);
        let results_buffer =
            manager.create_storage_buffer("results", (targets.count * 4) as u64, true);

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
                        resource: targets.buffer.as_entire_binding(),
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
            pass.dispatch_workgroups((targets.count as u32).div_ceil(256), 1, 1);
        }

        self.ctx.queue.submit(Some(encoder.finish()));

        // Read results
        manager.read_buffer_blocking(&results_buffer, targets.count)
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
    use crate::context::shared_test_context;

    /// A single `GpuTanimoto`, lazily built once from the crate's shared test
    /// `GpuContext` and reused by every test below — building each test's own
    /// `GpuTanimoto::new()` would each create an independent `GpuContext`,
    /// reintroducing the concurrent-context-creation hang this crate's test
    /// suite hit under the default parallel runner (see #19).
    fn shared_test_tanimoto() -> Option<&'static GpuTanimoto> {
        static GPU: std::sync::OnceLock<Option<GpuTanimoto>> = std::sync::OnceLock::new();
        GPU.get_or_init(|| {
            let ctx = shared_test_context()?.clone();
            match GpuTanimoto::from_context(ctx) {
                Ok(gpu) => Some(gpu),
                Err(e) => {
                    println!("GPU tanimoto pipeline unavailable, skipping: {}", e);
                    None
                }
            }
        })
        .as_ref()
    }

    #[test]
    fn test_compute_single_query_empty_targets_returns_empty_not_panic() {
        let Some(gpu) = shared_test_tanimoto() else {
            println!("GPU not available, skipping");
            return;
        };

        let query = vec![0u32; 64];
        let result = gpu.compute_single_query(&query, &[]).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_compute_single_query_empty_query_is_error() {
        let Some(gpu) = shared_test_tanimoto() else {
            println!("GPU not available, skipping");
            return;
        };

        let targets = vec![0u32; 64];
        assert!(gpu.compute_single_query(&[], &targets).is_err());
    }

    #[test]
    fn test_compute_all_pairs_empty_inputs_return_empty_not_panic() {
        let Some(gpu) = shared_test_tanimoto() else {
            println!("GPU not available, skipping");
            return;
        };

        let fp_size = 2048u32;
        let one_fp = vec![0u32; (fp_size / 32) as usize];

        assert!(
            gpu.compute_all_pairs(&[], &one_fp, fp_size)
                .unwrap()
                .is_empty()
        );
        assert!(
            gpu.compute_all_pairs(&one_fp, &[], fp_size)
                .unwrap()
                .is_empty()
        );
        assert!(gpu.compute_all_pairs(&[], &[], fp_size).unwrap().is_empty());
    }

    #[test]
    fn test_compute_all_pairs_zero_fp_size_is_error() {
        let Some(gpu) = shared_test_tanimoto() else {
            println!("GPU not available, skipping");
            return;
        };

        assert!(gpu.compute_all_pairs(&[1u32], &[1u32], 0).is_err());
    }
}
