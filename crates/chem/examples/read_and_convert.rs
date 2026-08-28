//! Reading a file, perceiving aromaticity, and writing it back out.
//!
//! ```sh
//! cargo run --example read_and_convert -- test.smi
//! ```
//!
//! [`chem::io::reader`] is the file-level layer: it splits a `.smi` by line or
//! an `.sdf` on `$$$$`, names records the file left unnamed, and — importantly
//! — hands back the records it could *not* read rather than logging them, so a
//! caller can count them, report them, or fail on them.

use chem::core::layout::ensure_coords;
use chem::io::aromaticity::detect_aromaticity;
use chem::io::reader::{self, Format};
use chem::io::sdf::write_sdf_all;
use chem::io::smiles_writer::write_smiles_for_molecule;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).unwrap_or_else(|| "test.smi".into());
    let text = std::fs::read_to_string(&path)?;

    // Only `.sdf` is SDF; `.smi`, `.txt` and no extension all read as SMILES.
    let outcome = reader::read(&text, Format::from_filename(&path));
    println!("{}: {} molecules read", path, outcome.len());

    // Failures come back rather than vanishing into a log, which is what lets
    // a batch tool exit non-zero on them.
    for skipped in outcome.skipped.iter().take(3) {
        println!("  line {} skipped: {}", skipped.position, skipped.error);
    }
    if outcome.skipped.len() > 3 {
        println!("  ...and {} more", outcome.skipped.len() - 3);
    }

    let mut molecules = Vec::new();
    let mut newly_aromatic = 0;
    for record in &outcome.records {
        let mut molecule = record.molecule.clone();

        // Kekulé input — explicit alternating double bonds — carries no
        // aromatic flags until something perceives them. Input already written
        // with lowercase atoms is flagged by the parser and unaffected.
        let before = molecule.atoms().iter().filter(|a| a.is_aromatic()).count();
        detect_aromaticity(&mut molecule);
        if molecule.atoms().iter().filter(|a| a.is_aromatic()).count() != before {
            newly_aromatic += 1;
        }

        // SMILES cannot carry coordinates, so they have to be computed before
        // anything can draw or store geometry.
        ensure_coords(&mut molecule);
        molecule.set_name(record.name.clone());
        molecules.push(molecule);
    }
    println!("aromaticity newly perceived in {newly_aromatic} molecules");

    if let Some(first) = molecules.first() {
        println!(
            "first molecule as SMILES: {}",
            write_smiles_for_molecule(first)
        );
    }

    // SDF is the only format here that can hold the coordinates just computed.
    let sdf = write_sdf_all(&molecules);
    println!(
        "written as SDF: {} bytes, {} records",
        sdf.len(),
        sdf.matches("$$$$").count()
    );
    Ok(())
}
