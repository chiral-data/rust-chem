use chem::core::atom::{Atom, Element};
use chem::core::molecule::Molecule;
use chem::gpu::GpuMorganFingerprint;

/// Regression test for issue #27: `dedup_environments_batch` used to dispatch
/// one WebGPU workgroup per molecule with no chunking
/// (`pass.dispatch_workgroups(num_molecules as u32, 1, 1)` against a
/// `@workgroup_size(1)` shader), so any batch over WebGPU's 65,535
/// dispatch-group-count cap panicked instead of completing. This drives a
/// batch just past that old cap to confirm it now completes without
/// panicking.
#[test]
fn test_batch_above_old_dispatch_cap_does_not_panic() {
    let Ok(gpu) = GpuMorganFingerprint::new() else {
        println!("GPU not available, skipping");
        return;
    };

    // One atom each: cheap per-molecule GPU work, but still exercises the
    // dedup dispatch at a molecule count above the old 65,535 cap.
    let num_molecules = 65_536 + 1000;
    let molecules: Vec<Molecule> = (0..num_molecules)
        .map(|_| {
            let mut mol = Molecule::new();
            mol.add_atom(Atom::new(Element::carbon()));
            mol.calculate_implicit_hydrogens();
            mol
        })
        .collect();

    let radius = 1;
    let fp_size = 128;

    let fps = gpu
        .generate_fingerprints_batch(&molecules, radius, fp_size, false, true, false)
        .unwrap();

    assert_eq!(fps.len(), num_molecules);
}
