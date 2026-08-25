use crate::error::GpuError;
use wgpu;

/// A compute pipeline with shader and bind group layout.
pub struct ComputePipeline {
    // The compiled shader pipeline (contains GPU shader code)
    pub pipeline: wgpu::ComputePipeline,
    // The layout describing how buffers connect to the shader
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl ComputePipeline {
    /// Create a compute pipeline from WGSL shader source.
    pub fn new(
        device: &wgpu::Device,
        shader_source: &str,
        entry_point: &str,
        label: &str,
    ) -> Result<ComputePipeline, GpuError> {
        // Create shader module (WGSL → SPIR-V → GPU machine code)
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some(&format!("{}_shader", label)),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });

        // Create bind group layout (what the GPU expects from the CPU)
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some(&format!("{}_bind_group_layout", label)),
            entries: &[
                // Input buffer
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
                // Output buffer
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
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

        // Create pipeline layout (container describing all bind groups the shader uses)
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some(&format!("{}_pipeline_layout", label)),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        // Create compute pipeline (WGSL → GPU executable kernel happens)
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(label),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some(entry_point),
            cache: None,
            compilation_options: Default::default(),
        });

        Ok(ComputePipeline {
            pipeline,
            bind_group_layout,
        })
    }

    /// Create a bind group for this pipeline.
    pub fn create_bind_group(
        &self,
        device: &wgpu::Device,
        input_buffer: &wgpu::Buffer,
        output_buffer: &wgpu::Buffer,
        label: &str,
    ) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some(label),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: input_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: output_buffer.as_entire_binding(),
                },
            ],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::shared_test_context;

    #[test]
    fn test_pipeline_creation() {
        let Some(ctx) = shared_test_context() else {
            return;
        };

        let shader = include_str!("shaders/hello.wgsl");
        let pipeline = ComputePipeline::new(&ctx.device, shader, "main", "test_pipeline");

        assert!(pipeline.is_ok());
    }
}
