use crate::error::GpuError;
use wgpu;

/// GPU context containing device and queue.
///
/// Cheaply `Clone`: `wgpu::Device`/`Queue` are thin `Arc`-backed handles, so
/// cloning shares the same underlying GPU device rather than creating a new
/// one.
#[derive(Clone)]
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

    /// Async twin of [`Self::new`] for callers that can't block the current
    /// thread while waiting on adapter/device requests (namely wasm32,
    /// where blocking the browser's single JS thread would deadlock).
    pub async fn new_async() -> Result<Self, GpuError> {
        match Self::new_async_with_preference(wgpu::PowerPreference::HighPerformance).await {
            Ok(ctx) => Ok(ctx),
            // Observed on Firefox/macOS (Apple Silicon): requestAdapter()
            // returns null when asked for "high-performance" specifically,
            // even though the same call with no power preference at all
            // succeeds and returns a real adapter — there's only one
            // unified GPU on Apple Silicon, so "high-performance" doesn't
            // clearly map to anything for that implementation to return.
            // Retry with no preference rather than give up.
            Err(GpuError::NoAdapter) => {
                log::warn!(
                    "No adapter found with HighPerformance preference; retrying with no power preference"
                );
                Self::new_async_with_preference(wgpu::PowerPreference::None).await
            }
            Err(e) => Err(e),
        }
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
        let required_limits = wgpu::Limits {
            max_storage_buffers_per_shader_stage: adapter_limits
                .max_storage_buffers_per_shader_stage
                .min(16)
                .max(wgpu::Limits::default().max_storage_buffers_per_shader_stage),
            ..Default::default()
        };

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

/// A single `GpuContext`, lazily created once and shared by every GPU test
/// across the crate (`context`, `buffers`, `pipeline`).
///
/// Each test previously called `GpuContext::new()` independently, which under
/// the default parallel test runner meant multiple threads racing to create
/// a `wgpu::Instance`/adapter/device concurrently against the same physical
/// GPU — that race hangs indefinitely on at least one driver (see #19).
/// Sharing one context removes the concurrent-creation race entirely; the
/// shared `Device`/`Queue` are `Send + Sync` and safe to use from multiple
/// test threads afterward.
#[cfg(test)]
pub(crate) fn shared_test_context() -> Option<&'static GpuContext> {
    static CTX: std::sync::OnceLock<Option<GpuContext>> = std::sync::OnceLock::new();
    CTX.get_or_init(|| {
        env_logger::try_init().ok();
        match GpuContext::new() {
            Ok(ctx) => {
                println!("GPU initialized: {:?}", ctx.adapter_info);
                Some(ctx)
            }
            Err(e) => {
                println!("GPU not available: {}", e);
                // Not a failure - GPU might not be available in CI
                None
            }
        }
    })
    .as_ref()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_init() {
        if let Some(ctx) = shared_test_context() {
            assert!(ctx.supports_compute());
        }
    }
}
