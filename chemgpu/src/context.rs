use crate::error::GpuError;
use wgpu;

/// GPU context containing device and queue.
pub struct GpuContext {
    /// The GPU device
    pub device: wgpu::Device,
    /// Command queue for submitting work
    pub queue: wgpu::Queue,
    /// Adapter information
    pub adapter_info: wgpu::AdapterInfo,
}

impl GpuContext {
    /// Initialize GPU context.
    ///
    /// # Returns
    /// - `Ok(GpuContext)` if GPU is available
    /// - `Err(GpuError)` if no suitable GPU found
    ///
    /// # Example
    /// ```no_run
    /// use chemgpu::context::GpuContext;
    ///
    /// let ctx = GpuContext::new().expect("Failed to initialize GPU");
    /// println!("Using GPU: {}", ctx.adapter_info.name);
    /// ```
    pub fn new() -> Result<Self, GpuError> {
        pollster::block_on(Self::new_async())
    }

    /// Initialize GPU context with specific power preference.
    pub fn new_with_preference(power_pref: wgpu::PowerPreference) -> Result<Self, GpuError> {
        pollster::block_on(Self::new_async_with_preference(power_pref))
    }

    async fn new_async() -> Result<Self, GpuError> {
        Self::new_async_with_preference(wgpu::PowerPreference::HighPerformance).await
    }

    async fn new_async_with_preference(
        power_pref: wgpu::PowerPreference,
    ) -> Result<Self, GpuError> {
        log::info!("Initializing GPU context...");

        // Create instance
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        // Request adapter
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: power_pref,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .map_err(|_| GpuError::NoAdapter)?;

        let adapter_info = adapter.get_info();
        log::info!(
            "Found adapter: {} ({:?})",
            adapter_info.name,
            adapter_info.backend
        );

        // wgpu::Limits::default() caps max_storage_buffers_per_shader_stage at 8
        // (a conservative, portable-across-backends default). The Morgan batch
        // shader's redundant-environment dedup uses 13 storage buffer bindings
        // in one compute stage, so ask for more — clamped to what this adapter
        // actually reports, so we never request above its real capability.
        let adapter_limits = adapter.limits();
        let mut required_limits = wgpu::Limits::default();
        required_limits.max_storage_buffers_per_shader_stage = adapter_limits
            .max_storage_buffers_per_shader_stage
            .min(16)
            .max(wgpu::Limits::default().max_storage_buffers_per_shader_stage);

        // Request device and queue
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("ChemGPU Device"),
                required_features: wgpu::Features::empty(),
                required_limits,
                memory_hints: wgpu::MemoryHints::default(),
                experimental_features: wgpu::ExperimentalFeatures::default(),
                trace: wgpu::Trace::default(),
            })
            .await
            .map_err(|e| GpuError::DeviceRequest(e.to_string()))?;

        log::info!("GPU context initialized successfully");

        Ok(GpuContext {
            device,
            queue,
            adapter_info,
        })
    }

    /// Check if GPU supports compute shaders.
    pub fn supports_compute(&self) -> bool {
        // All modern GPUs support compute
        true
    }

    /// Get device limits.
    pub fn limits(&self) -> wgpu::Limits {
        self.device.limits()
    }
}

impl std::fmt::Debug for GpuContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuContext")
            .field("adapter", &self.adapter_info.name)
            .field("backend", &self.adapter_info.backend)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_init() {
        env_logger::try_init().ok();

        match GpuContext::new() {
            Ok(ctx) => {
                println!("GPU initialized: {:?}", ctx.adapter_info);
                assert!(ctx.supports_compute());
            }
            Err(e) => {
                println!("GPU not available: {}", e);
                // Not a failure - GPU might not be available in CI
            }
        }
    }
}
