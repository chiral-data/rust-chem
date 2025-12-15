use bitvec::prelude::*;
use chemcore::atom::{Atom, Element};
use chemcore::bond::{Bond, BondOrder};
use chemcore::molecule::Molecule;
use chemfp::morgan::MorganFingerprint;
use chemfp::tanimoto;
use chemgpu::{GpuMorganFingerprint, GpuTanimoto};
use std::time::Instant; // ADD THIS

fn create_molecules(count: usize) -> Vec<Molecule> {
    let mut molecules = Vec::new();

    for i in 0..count {
        let mut mol = Molecule::new();

        // Create different molecules by varying ring size
        let ring_size = 5 + (i % 3); // 5, 6, or 7-membered rings

        // Add carbons
        for _ in 0..ring_size {
            mol.add_atom(Atom::new(Element::carbon()).with_aromatic(true));
        }

        // Add ring bonds
        for j in 0..ring_size {
            let next = (j + 1) % ring_size;
            mol.add_bond(Bond::new(j, next, BondOrder::Aromatic).with_aromatic(true))
                .unwrap();
        }

        mol.calculate_implicit_hydrogens();
        molecules.push(mol);
    }

    molecules
}

fn bitvec_to_u32_vec(bv: &BitVec, fp_size: u32) -> Vec<u32> {
    // CHANGED: removed bitvec::prelude::
    let num_words = ((fp_size + 31) / 32) as usize;
    let mut result = vec![0u32; num_words];

    for word_idx in 0..num_words {
        let mut word = 0u32;
        for bit_idx in 0..32 {
            let global_bit = word_idx * 32 + bit_idx;
            if global_bit < fp_size as usize && bv[global_bit] {
                word |= 1u32 << bit_idx;
            }
        }
        result[word_idx] = word;
    }

    result
}

fn main() {
    println!("=== GPU Morgan + Tanimoto Full Benchmark ===\n");

    let num_molecules = 100;
    let radius = 2;
    let fp_size = 2048;

    println!("Configuration:");
    println!("  Molecules: {}", num_molecules);
    println!("  Radius: {}", radius);
    println!("  Fingerprint size: {} bits\n", fp_size);

    // Create test molecules
    println!("Creating {} test molecules...", num_molecules);
    let molecules = create_molecules(num_molecules);

    // ========================================================================
    // PART 1: Morgan Fingerprint Generation
    // ========================================================================

    println!("\n--- PART 1: Morgan Fingerprint Generation ---");

    // CPU Morgan
    println!("\nCPU Morgan fingerprints...");
    let cpu_morgan_start = Instant::now();
    let cpu_fps: Vec<_> = molecules
        .iter()
        .map(|mol| {
            MorganFingerprint::get_fingerprint_as_bitvec(
                mol, radius, fp_size, None, None, false, true, false,
            )
            .unwrap()
            // REMOVED: .into()
        })
        .collect();
    let cpu_morgan_time = cpu_morgan_start.elapsed();
    println!("  Time: {:?}", cpu_morgan_time);

    // Convert to u32 vectors
    let cpu_fps_u32: Vec<Vec<u32>> = cpu_fps
        .iter()
        .map(|fp| bitvec_to_u32_vec(fp, fp_size))
        .collect();

    // GPU Morgan (batched)
    println!("\nGPU Morgan fingerprints (batched)...");
    let gpu_morgan_gen = GpuMorganFingerprint::new().unwrap();
    let gpu_morgan_start = Instant::now();
    let gpu_fps = gpu_morgan_gen
        .generate_fingerprints_batch(&molecules, radius, fp_size, false, true, false)
        .unwrap();
    let gpu_morgan_time = gpu_morgan_start.elapsed();
    println!("  Time: {:?}", gpu_morgan_time);

    // Validate Morgan
    println!("\nValidating Morgan fingerprints...");
    let mut fp_matches = 0;
    let mut fp_mismatches = 0;

    for (i, (cpu_fp, gpu_fp)) in cpu_fps_u32.iter().zip(gpu_fps.iter()).enumerate() {
        if cpu_fp == gpu_fp {
            fp_matches += 1;
        } else {
            fp_mismatches += 1;
            println!("  Molecule {}: MISMATCH", i);
        }
    }

    println!("  Matches: {}/{}", fp_matches, num_molecules);
    if fp_mismatches > 0 {
        println!("  ✗ FAILED: {} fingerprints differ", fp_mismatches);
        std::process::exit(1);
    }
    println!("  ✓ All Morgan fingerprints match!");
    println!(
        "  Morgan speedup: {:.2}x",
        cpu_morgan_time.as_secs_f64() / gpu_morgan_time.as_secs_f64()
    );

    // ========================================================================
    // PART 2: Tanimoto Similarity Computation
    // ========================================================================

    println!("\n--- PART 2: Tanimoto Similarity Computation ---");

    // Use first molecule as query, rest as targets
    let query_fp = &cpu_fps[0];
    let target_fps = &cpu_fps[1..];

    println!("\nComputing 1 query vs {} targets", target_fps.len());

    // CPU Tanimoto
    println!("\nCPU Tanimoto...");
    let cpu_tanimoto_start = Instant::now();
    let cpu_similarities = tanimoto::batch_tanimoto(query_fp, target_fps).unwrap();
    let cpu_tanimoto_time = cpu_tanimoto_start.elapsed();
    println!("  Time: {:?}", cpu_tanimoto_time);

    // GPU Tanimoto
    println!("\nGPU Tanimoto...");
    let gpu_tanimoto = GpuTanimoto::new().unwrap();
    let query_u32 = &cpu_fps_u32[0];
    let targets_u32: Vec<u32> = cpu_fps_u32[1..].iter().flatten().copied().collect();

    let gpu_tanimoto_start = Instant::now();
    let gpu_similarities = gpu_tanimoto
        .compute_single_query(query_u32, &targets_u32, fp_size)
        .unwrap();
    let gpu_tanimoto_time = gpu_tanimoto_start.elapsed();
    println!("  Time: {:?}", gpu_tanimoto_time);

    // Validate Tanimoto
    println!("\nValidating Tanimoto similarities...");
    let mut sim_matches = 0;
    let mut sim_mismatches = 0;
    let tolerance = 1e-6;

    for (i, (cpu_sim, gpu_sim)) in cpu_similarities
        .iter()
        .zip(gpu_similarities.iter())
        .enumerate()
    {
        let diff = (cpu_sim - *gpu_sim as f64).abs();
        if diff < tolerance {
            sim_matches += 1;
        } else {
            sim_mismatches += 1;
            println!(
                "  Target {}: CPU={:.6}, GPU={:.6}, diff={:.6}",
                i, cpu_sim, gpu_sim, diff
            );
        }
    }

    println!("  Matches: {}/{}", sim_matches, target_fps.len());
    if sim_mismatches > 0 {
        println!("  ✗ FAILED: {} similarities differ", sim_mismatches);
        std::process::exit(1);
    }
    println!("  ✓ All similarities match!");
    println!(
        "  Tanimoto speedup: {:.2}x",
        cpu_tanimoto_time.as_secs_f64() / gpu_tanimoto_time.as_secs_f64()
    );

    // ========================================================================
    // PART 3: All-pairs similarity (if small dataset)
    // ========================================================================

    if num_molecules <= 20 {
        println!(
            "\n--- PART 3: All-pairs Similarity ({}x{}) ---",
            num_molecules, num_molecules
        );

        // CPU all-pairs
        println!("\nCPU all-pairs...");
        let cpu_allpairs_start = Instant::now();
        let mut cpu_matrix = Vec::new();
        for query in &cpu_fps {
            for target in &cpu_fps {
                let sim = tanimoto::tanimoto_similarity(query, target).unwrap();
                cpu_matrix.push(sim);
            }
        }
        let cpu_allpairs_time = cpu_allpairs_start.elapsed();
        println!("  Time: {:?}", cpu_allpairs_time);

        // GPU all-pairs
        println!("\nGPU all-pairs...");
        let all_fps_u32: Vec<u32> = cpu_fps_u32.iter().flatten().copied().collect();

        let gpu_allpairs_start = Instant::now();
        let gpu_matrix = gpu_tanimoto
            .compute_all_pairs(&all_fps_u32, &all_fps_u32, fp_size)
            .unwrap();
        let gpu_allpairs_time = gpu_allpairs_start.elapsed();
        println!("  Time: {:?}", gpu_allpairs_time);

        // Validate
        println!("\nValidating all-pairs matrix...");
        let mut matrix_matches = 0;
        let mut matrix_mismatches = 0;

        for (i, (cpu_val, gpu_val)) in cpu_matrix.iter().zip(gpu_matrix.iter()).enumerate() {
            let diff = (cpu_val - *gpu_val as f64).abs();
            if diff < tolerance {
                matrix_matches += 1;
            } else {
                matrix_mismatches += 1;
                if matrix_mismatches <= 5 {
                    println!("  Position {}: CPU={:.6}, GPU={:.6}", i, cpu_val, gpu_val);
                }
            }
        }

        println!("  Matches: {}/{}", matrix_matches, cpu_matrix.len());
        if matrix_mismatches > 0 {
            println!("  ✗ FAILED: {} values differ", matrix_mismatches);
            std::process::exit(1);
        }
        println!("  ✓ All-pairs matrix matches!");
        println!(
            "  All-pairs speedup: {:.2}x",
            cpu_allpairs_time.as_secs_f64() / gpu_allpairs_time.as_secs_f64()
        );
    }

    // ========================================================================
    // SUMMARY
    // ========================================================================

    println!("\n=== SUMMARY ===");
    println!("✓ Morgan fingerprints: GPU matches CPU");
    println!("✓ Tanimoto similarity: GPU matches CPU");
    println!("\nOverall times:");
    println!(
        "  Morgan:   CPU={:?}, GPU={:?}",
        cpu_morgan_time, gpu_morgan_time
    );
    println!(
        "  Tanimoto: CPU={:?}, GPU={:?}",
        cpu_tanimoto_time, gpu_tanimoto_time
    );
    println!("\nSpeedups:");
    println!(
        "  Morgan GPU vs CPU: {:.2}x",
        cpu_morgan_time.as_secs_f64() / gpu_morgan_time.as_secs_f64()
    );
    println!("\nALL TESTS PASSED!");
}
