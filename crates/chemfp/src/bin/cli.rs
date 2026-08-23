use chemfp::fingerprint::ROMol;
use chemfp::morgan::MorganFingerprint;
use chemio::smiles::parse_smiles;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <SMILES> [radius] [nbits]", args[0]);
        eprintln!("Examples:");
        eprintln!("  {} 'CCO' 2 2048", args[0]);
        eprintln!("  {} 'c1ccccc1' 3 1024", args[0]);
        std::process::exit(1);
    }

    let smiles = &args[1];
    let radius = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2);
    let nbits = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2048);

    eprintln!("=== Morgan Fingerprint Generator ===");
    eprintln!("SMILES: {}", smiles);
    eprintln!("Radius: {}", radius);
    eprintln!("Bits:   {}", nbits);
    eprintln!();

    // Parse SMILES
    let mol: ROMol = match parse_smiles(smiles) {
        Ok(mut mol) => {
            mol.calculate_implicit_hydrogens();

            eprintln!("Molecule Info:");
            eprintln!("  Atoms: {}", mol.num_atoms());
            eprintln!("  Bonds: {}", mol.num_bonds());
            eprintln!("  Formula: {}", mol.formula());

            // Print atom details
            eprintln!("\n  Atom Details:");
            for (idx, atom) in mol.atoms().iter().enumerate() {
                eprintln!(
                    "    Atom {}: {} (aromatic={}, degree={})",
                    idx,
                    atom.element().symbol(),
                    atom.is_aromatic(),
                    mol.degree(idx)
                );
            }
            eprintln!();

            mol
        }
        Err(e) => {
            eprintln!("Error parsing SMILES: {}", e);
            std::process::exit(1);
        }
    };

    // Generate fingerprint with additional output
    eprintln!("Generating fingerprint...");

    let fp = match MorganFingerprint::get_fingerprint_as_bitvec(
        &mol, radius, nbits, None, None, false, true, false,
    ) {
        Ok(fp) => fp,
        Err(e) => {
            eprintln!("ERROR: {:?}", e);
            std::process::exit(1);
        }
    };

    // Analyze results
    let on_bits: Vec<usize> = fp.iter_ones().collect();
    let density = on_bits.len() as f64 / nbits as f64;

    eprintln!("Results:");
    eprintln!("  Set bits: {}", on_bits.len());
    eprintln!("  Density:  {:.4}%", density * 100.0);
    eprintln!();

    if on_bits.len() <= 20 {
        eprintln!("All set bits:");
        for &bit in &on_bits {
            eprintln!("  {}", bit);
        }
    } else {
        eprintln!("First 20 set bits:");
        for &bit in on_bits.iter().take(20) {
            eprintln!("  {}", bit);
        }
        eprintln!("  ... and {} more", on_bits.len() - 20);
    }
    eprintln!();

    // Output for comparison (CSV format)
    println!(
        "{}",
        on_bits
            .iter()
            .map(|b| b.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );
}
