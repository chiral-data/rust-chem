use chem::core::atom::{Atom, Element};
use chem::core::bond::{Bond, BondOrder};
use chem::core::molecule::Molecule;
use chem::fp::morgan::MorganFingerprint;
use chem::gpu::GpuMorganFingerprint;

/// Regression test for issue #12: the GPU Morgan shader used to collect each
/// atom's neighbors into a fixed-size `array<vec2<u32>, 16>`, silently
/// truncating (and diverging from the CPU reference, which has no such cap)
/// any atom bonded to more than 16 neighbors.
#[test]
fn test_gpu_matches_cpu_for_high_degree_atom() {
    let Ok(gpu) = GpuMorganFingerprint::new() else {
        println!("GPU not available, skipping");
        return;
    };

    // A central atom bonded to 20 leaves — exceeds the old 16-neighbor cap.
    // Chemically nonsensical, but a valid graph for stress-testing degree.
    let degree = 20;
    let mut mol = Molecule::new();
    let center = mol.add_atom(Atom::new(Element::carbon()));
    for i in 0..degree {
        let leaf = mol.add_atom(Atom::new(Element::carbon()));
        let order = if i % 3 == 0 {
            BondOrder::Double
        } else {
            BondOrder::Single
        };
        mol.add_bond(Bond::new(center, leaf, order)).unwrap();
    }
    mol.calculate_implicit_hydrogens();

    let radius = 1;
    let fp_size = 2048;

    let cpu_fp = MorganFingerprint::get_fingerprint_as_bitvec(
        &mol, radius, fp_size, None, None, false, true, false,
    )
    .unwrap();

    let gpu_fps = gpu
        .generate_fingerprints_batch(&[mol], radius, fp_size, false, true, false)
        .unwrap();

    let cpu_bits: Vec<usize> = cpu_fp
        .iter()
        .enumerate()
        .filter(|(_, b)| **b)
        .map(|(i, _)| i)
        .collect();
    let gpu_bits: Vec<usize> = (0..fp_size as usize)
        .filter(|&i| (gpu_fps[0][i / 32] >> (i % 32)) & 1 != 0)
        .collect();

    assert_eq!(
        cpu_bits, gpu_bits,
        "GPU fingerprint diverged from CPU for a {}-neighbor atom",
        degree
    );
}
