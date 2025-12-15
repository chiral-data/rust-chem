use chemcore::atom::{Atom, Element};
use chemcore::bond::{Bond, BondOrder};
use chemcore::molecule::Molecule;
use chemfp::morgan::MorganFingerprint;
use chemgpu::GpuMorganFingerprint;
use std::time::Instant;

fn create_test_molecule() -> Molecule {
    // Benzene ring
    let mut mol = Molecule::new();

    // 6 carbons
    for _ in 0..6 {
        mol.add_atom(Atom::new(Element::carbon()).with_aromatic(true));
    }

    // Ring bonds
    mol.add_bond(Bond::new(0, 1, BondOrder::Aromatic).with_aromatic(true))
        .unwrap();
    mol.add_bond(Bond::new(1, 2, BondOrder::Aromatic).with_aromatic(true))
        .unwrap();
    mol.add_bond(Bond::new(2, 3, BondOrder::Aromatic).with_aromatic(true))
        .unwrap();
    mol.add_bond(Bond::new(3, 4, BondOrder::Aromatic).with_aromatic(true))
        .unwrap();
    mol.add_bond(Bond::new(4, 5, BondOrder::Aromatic).with_aromatic(true))
        .unwrap();
    mol.add_bond(Bond::new(5, 0, BondOrder::Aromatic).with_aromatic(true))
        .unwrap();

    mol.calculate_implicit_hydrogens();
    mol
}

fn bitvec_to_u32_vec(bv: &bitvec::prelude::BitVec, fp_size: u32) -> Vec<u32> {
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

fn compare_fingerprints(cpu: &[u32], gpu: &[u32]) -> bool {
    if cpu.len() != gpu.len() {
        println!("Length mismatch: CPU={}, GPU={}", cpu.len(), gpu.len());
        return false;
    }

    let mut matches = 0;
    let mut mismatches = 0;

    for (i, (c, g)) in cpu.iter().zip(gpu.iter()).enumerate() {
        if c != g {
            println!("Word {}: CPU={:032b}, GPU={:032b}", i, c, g);
            mismatches += 1;
        } else {
            matches += 1;
        }
    }

    println!("Matches: {}, Mismatches: {}", matches, mismatches);
    mismatches == 0
}

fn main() {
    println!("=== GPU Morgan Fingerprint Validation ===\n");

    let mol = create_test_molecule();
    let radius = 2;
    let fp_size = 2048;

    println!(
        "Molecule: {} atoms, {} bonds",
        mol.num_atoms(),
        mol.num_bonds()
    );
    println!("Parameters: radius={}, fp_size={}\n", radius, fp_size);

    // CPU fingerprint
    println!("Computing CPU fingerprint...");
    let cpu_start = Instant::now();
    let cpu_fp = MorganFingerprint::get_fingerprint_as_bitvec(
        &mol, radius, fp_size, None, None, false, true, false,
    )
    .unwrap();
    let cpu_time = cpu_start.elapsed();
    let cpu_vec = bitvec_to_u32_vec(&cpu_fp, fp_size);

    println!("CPU time: {:?}", cpu_time);
    println!("CPU on-bits: {}\n", cpu_fp.count_ones());
    println!("CPU hash values:");
    println!("  Layer 0: 37324 (all atoms should be this)");
    println!("  Layer 1: 2837459962");
    println!("  Layer 2: 501318127");
    println!(
        "  Bit positions: {}, {}, {}",
        37324 % fp_size,
        2837459962 % fp_size,
        501318127 % fp_size
    );
    println!("CPU on-bits: {}", cpu_fp.count_ones());
    println!("CPU on-bit positions: {:?}", {
        let mut positions = Vec::new();
        for i in 0..fp_size as usize {
            if cpu_fp[i] {
                positions.push(i);
            }
        }
        positions
    });

    // GPU fingerprint (using batched API with single molecule)
    println!("Computing GPU fingerprint...");
    let gpu_generator = GpuMorganFingerprint::new().unwrap();

    let gpu_start = Instant::now();
    let gpu_vecs = gpu_generator
        .generate_fingerprints_batch(&[mol.clone()], radius, fp_size, false, true, false)
        .unwrap();
    let gpu_vec = gpu_vecs[0].clone();
    let gpu_time = gpu_start.elapsed();

    let gpu_on_bits: usize = gpu_vec.iter().map(|w| w.count_ones() as usize).sum();

    println!("GPU time: {:?}", gpu_time);
    println!("GPU on-bits: {}\n", gpu_on_bits);

    // Validation
    println!("=== VALIDATION ===");
    let valid = compare_fingerprints(&cpu_vec, &gpu_vec);

    if valid {
        println!("\n✓ SUCCESS: GPU output matches CPU exactly!");
        println!(
            "Speedup: {:.2}x",
            cpu_time.as_secs_f64() / gpu_time.as_secs_f64()
        );
    } else {
        println!("\n✗ FAILED: GPU output differs from CPU");
        std::process::exit(1);
    }
}
