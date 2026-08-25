use chem::gpu::{BufferManager, ComputePipeline, GpuContext, init_logging};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    init_logging();

    println!("=== GPU Hello World ===\n");

    // Initialize GPU
    println!("1. Initializing GPU...");
    let ctx = GpuContext::new()?;
    println!(
        "   Using: {} ({:?})\n",
        ctx.adapter_info.name, ctx.adapter_info.backend
    );

    // Create buffer manager
    let buffer_mgr = BufferManager::new(&ctx.device, &ctx.queue);

    // Input data
    let input_data: Vec<u32> = vec![1, 2, 3, 4, 5, 6, 7, 8];
    println!("2. Input data: {:?}\n", input_data);

    // Create buffers
    println!("3. Creating GPU buffers...");
    let input_size = (input_data.len() * std::mem::size_of::<u32>()) as u64;
    let input_buffer = buffer_mgr.create_storage_buffer("input", input_size, false);
    let output_buffer = buffer_mgr.create_storage_buffer("output", input_size, true);

    // Write input data
    buffer_mgr.write_buffer(&input_buffer, &input_data);
    println!("   Buffers created and data uploaded\n");

    // Load shader
    println!("4. Loading compute shader...");
    let shader_source = include_str!("../src/gpu/shaders/hello.wgsl");
    let pipeline = ComputePipeline::new(&ctx.device, shader_source, "main", "hello_pipeline")?;
    println!("   Shader compiled\n");

    // Create bind group
    let bind_group = pipeline.create_bind_group(
        &ctx.device,
        &input_buffer,
        &output_buffer,
        "hello_bind_group",
    );

    // Execute compute shader
    println!("5. Executing compute shader...");
    let mut encoder = ctx
        .device
        .create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("compute_encoder"),
        });

    {
        let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("compute_pass"),
            timestamp_writes: None,
        });
        compute_pass.set_pipeline(&pipeline.pipeline);
        compute_pass.set_bind_group(0, &bind_group, &[]);

        // Dispatch: (num_elements + 64 - 1) / 64 workgroups
        let num_workgroups = (input_data.len() as u32 + 63).div_ceil(64);
        compute_pass.dispatch_workgroups(num_workgroups, 1, 1);
    }

    ctx.queue.submit(Some(encoder.finish()));
    println!("   Compute shader dispatched\n");

    // Read results
    println!("6. Reading results from GPU...");
    let output_data: Vec<u32> =
        buffer_mgr.read_buffer_blocking(&output_buffer, input_data.len())?;
    println!("   Output data: {:?}\n", output_data);

    // Verify
    println!("7. Verification:");
    let mut all_correct = true;
    for (i, (&input, &output)) in input_data.iter().zip(output_data.iter()).enumerate() {
        let expected = input * 2;
        let correct = output == expected;
        all_correct &= correct;
        println!(
            "   [{}] {} * 2 = {} {}",
            i,
            input,
            output,
            if correct { "✓" } else { "✗" }
        );
    }

    println!(
        "\n{}",
        if all_correct {
            "SUCCESS: All values correct!"
        } else {
            "FAILURE: Some values incorrect"
        }
    );

    Ok(())
}
