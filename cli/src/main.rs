use chemfp::morgan::{morgan_fingerprint, MorganOptions};
use chemio::smiles::parse_smiles;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 4 {
        eprintln!("Usage: {} <smiles> <radius> <nbits>", args[0]);
        eprintln!("Example: {} CCO 2 2048", args[0]);
        std::process::exit(1);
    }

    let smiles = &args[1];
    let radius: u32 = args[2].parse().expect("Invalid radius");
    let nbits: usize = args[3].parse().expect("Invalid nbits");

    // Parse SMILES
    let mol = match parse_smiles(smiles) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("ERROR: Failed to parse SMILES: {}", e);
            std::process::exit(1);
        }
    };

    // Generate fingerprint
    let options = MorganOptions {
        radius,
        nbits,
        bits_per_feature: 1,
    };

    let fp = match morgan_fingerprint(&mol, &options) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("ERROR: Failed to generate fingerprint: {}", e);
            std::process::exit(1);
        }
    };

    // Output bit positions as comma-separated list
    let mut positions = Vec::new();
    for i in 0..fp.size() {
        if fp.get(i) {
            positions.push(i.to_string());
        }
    }

    println!("{}", positions.join(","));
}
