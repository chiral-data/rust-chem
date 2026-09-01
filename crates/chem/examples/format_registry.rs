//! The format registry: what this build supports, and every way to ask.
//!
//! ```sh
//! cargo run -p chem --example format_registry
//! ```
//!
//! [`chem::io::format`] is the table a `-L formats` listing will iterate.
//! Nothing in the CLI reaches most of it yet — `--format` still goes through
//! clap's own `smiles|sdf` mirror — so this is where the registry can actually
//! be seen working.

use chem::core::molecule::Molecule;
use chem::io::format::{self, Format};
use chem::io::reader;

fn main() {
    // What a format listing prints. Note `codes` is a slice: SDF answers to
    // four of them, which is why lookup is not a single-string compare.
    println!(
        "{:<16} {:<18} {:<22} {:<4} CATEGORY",
        "CODES", "NAME", "EXTENSIONS", "R/W"
    );
    for f in format::all() {
        let exts: Vec<_> = f.extensions().iter().map(|e| format!(".{e}")).collect();
        let rw = format!(
            "{}{}",
            if f.can_read() { "r" } else { "-" },
            if f.can_write() { "w" } else { "-" },
        );
        println!(
            "{:<16} {:<18} {:<22} {:<4} {}",
            f.codes().join(","),
            f.name(),
            exts.join(","),
            rw,
            f.category().label(),
        );
    }

    // Lookup by code, case-insensitively, across every alias.
    println!("\n-- from_code --");
    for code in ["smi", "smiles", "sdf", "sd", "mol", "mdl", "SDF", "pdb"] {
        match Format::from_code(code) {
            Some(f) => println!("  {code:>6} -> {f:?}  ({f})"),
            None => println!("  {code:>6} -> none, not registered in this build"),
        }
    }

    // The SMILES fallback, which is the rule that used to exist three times.
    println!("\n-- from_filename --");
    for name in [
        "a.smi", "a.sdf", "A.SDF", "a.txt", "a.mol", "a.pdb", "a.sdf.gz", "a", "-",
    ] {
        println!("  {name:>10} -> {}", Format::from_filename(name).label());
    }

    // Round-tripping through the descriptor's own writer, then reading back.
    println!("\n-- write, then read back --");
    let outcome = reader::read("CCO ethanol\nc1ccccc1 benzene\n", Format::SMILES);
    let records: Vec<(String, Molecule)> = outcome
        .records
        .iter()
        .map(|r| (r.name.clone(), r.molecule.clone()))
        .collect();
    for f in format::all() {
        let text = f
            .write(&records)
            .expect("every registered format writes today");
        let back = reader::read(&text, f);
        println!(
            "  {:<6} {:>6} bytes -> {} molecules back",
            f.label(),
            text.len(),
            back.len()
        );
    }

    // Pinned, because a drifted constant would silently swap two formats.
    assert_eq!(Format::from_code("mol"), Some(Format::SDF));
    assert_eq!(Format::SMILES.label(), "SMILES");
    assert_eq!(Format::SDF.name(), "MDL MOL format");
    // `mol` is a *code*, not an extension — so a `.mol` file still reads as SMILES.
    assert_eq!(Format::from_filename("a.mol"), Format::SMILES);
    println!("\nall assertions held");
}
