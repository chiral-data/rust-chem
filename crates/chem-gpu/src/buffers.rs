use crate::error::GpuError;
use bytemuck;
use wgpu;

/// Manager for GPU buffers.
pub struct BufferManager<'a> {
    device: &'a wgpu::Device,
    queue: &'a wgpu::Queue,
}

impl<'a> BufferManager<'a> {
    /// Create a new buffer manager.
    pub fn new(device: &'a wgpu::Device, queue: &'a wgpu::Queue) -> Self {
        BufferManager { device, queue }
    }

    /// Create a storage buffer (read/write on GPU).
    ///
    /// # Arguments
    /// - `label`: Debug label
    /// - `size`: Size in bytes
    /// - `readable`: Whether CPU can read from this buffer
    pub fn create_storage_buffer(&self, label: &str, size: u64, readable: bool) -> wgpu::Buffer {
        let mut usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST;
        if readable {
            usage |= wgpu::BufferUsages::COPY_SRC;
        }

        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size,
            usage,
            mapped_at_creation: false,
        })
    }

    /// Create a uniform buffer (read-only constants on GPU).
    pub fn create_uniform_buffer(&self, label: &str, size: u64) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    }

    /// Create a staging buffer (for reading results back to CPU).
    pub fn create_staging_buffer(&self, label: &str, size: u64) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    }

    /// Write data to a buffer.
    pub fn write_buffer<T: bytemuck::Pod>(&self, buffer: &wgpu::Buffer, data: &[T]) {
        self.queue
            .write_buffer(buffer, 0, bytemuck::cast_slice(data));
    }

    /// Write data to a buffer at a given byte offset — for populating one
    /// region of a larger buffer shared by several logically distinct
    /// arrays (see chem-gpu/src/shaders/morgan.wgsl's combined `scratch`
    /// binding).
    pub fn write_buffer_at<T: bytemuck::Pod>(
        &self,
        buffer: &wgpu::Buffer,
        offset: u64,
        data: &[T],
    ) {
        self.queue
            .write_buffer(buffer, offset, bytemuck::cast_slice(data));
    }

    /// Read data from a buffer (async operation).
    pub async fn read_buffer<T: bytemuck::Pod>(
        &self,
        buffer: &wgpu::Buffer,
        size: usize,
    ) -> Result<Vec<T>, GpuError> {
        self.read_buffer_at(buffer, 0, size).await
    }

    /// Read data from a given byte offset within a buffer (async operation)
    /// — for reading back one region of a larger combined buffer.
    pub async fn read_buffer_at<T: bytemuck::Pod>(
        &self,
        buffer: &wgpu::Buffer,
        offset: u64,
        size: usize,
    ) -> Result<Vec<T>, GpuError> {
        let staging =
            self.create_staging_buffer("read_staging", (size * std::mem::size_of::<T>()) as u64);

        // Copy GPU buffer to staging
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("read_encoder"),
            });
        encoder.copy_buffer_to_buffer(buffer, offset, &staging, 0, staging.size());
        self.queue.submit(Some(encoder.finish()));

        // Map staging buffer and read
        let buffer_slice = staging.slice(..);
        let (sender, receiver) = futures::channel::oneshot::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            sender.send(result).ok();
        });

        let _ = self.device.poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: None,
        });

        receiver
            .await
            .map_err(|_| GpuError::BufferError("Channel closed".to_string()))?
            .map_err(|e| GpuError::BufferError(e.to_string()))?;

        let data = buffer_slice.get_mapped_range();
        let result: Vec<T> = bytemuck::cast_slice(&data).to_vec();

        drop(data);
        staging.unmap();

        Ok(result)
    }

    /// Read buffer synchronously (blocking).
    pub fn read_buffer_blocking<T: bytemuck::Pod>(
        &self,
        buffer: &wgpu::Buffer,
        size: usize,
    ) -> Result<Vec<T>, GpuError> {
        pollster::block_on(self.read_buffer(buffer, size))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::shared_test_context;

    #[test]
    fn test_buffer_creation() {
        let Some(ctx) = shared_test_context() else {
            return; // Skip if no GPU
        };

        let manager = BufferManager::new(&ctx.device, &ctx.queue);

        let storage = manager.create_storage_buffer("test_storage", 1024, true);
        assert_eq!(storage.size(), 1024);

        let uniform = manager.create_uniform_buffer("test_uniform", 256);
        assert_eq!(uniform.size(), 256);
    }

    #[test]
    fn test_buffer_read_write() {
        let Some(ctx) = shared_test_context() else {
            return;
        };

        let manager = BufferManager::new(&ctx.device, &ctx.queue);

        // Write data
        let data: Vec<u32> = vec![1, 2, 3, 4, 5];
        let buffer = manager.create_storage_buffer("test_rw", (data.len() * 4) as u64, true);
        manager.write_buffer(&buffer, &data);

        // Read back
        let result: Vec<u32> = manager.read_buffer_blocking(&buffer, data.len()).unwrap();
        assert_eq!(result, data);
    }
}
