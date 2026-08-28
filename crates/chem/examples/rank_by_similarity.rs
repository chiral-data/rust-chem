//! Fingerprint a small library and rank it against a query.
//!
//! ```sh
//! cargo run --example rank_by_similarity
//! ```
//!
//! This is the whole search workflow in one file: parse, fingerprint, compare.
//! For a real library use [`chem::search::FingerprintSearch`], which picks a
//! backend and batches the work; this shows what it does underneath.

use chem::fp::morgan::MorganFingerprint;
use chem::fp::tanimoto::tanimoto_similarity;
use chem::io::smiles::parse_smiles;

const RADIUS: u32 = 2;
const BITS: u32 = 2048;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let library = [
        ("phenol", "c1ccccc1O"),
        // The same molecule, written from the other end. Morgan fingerprints
        // describe atom environments rather than the string, so these must
        // score exactly 1.000 — which they did not until #167, because the
        // parser gave one spelling an aromatic C-O bond.
        ("phenol, written backwards", "Oc1ccccc1"),
        ("p-cresol", "Cc1ccc(O)cc1"),
        ("toluene", "Cc1ccccc1"),
        ("benzene", "c1ccccc1"),
        ("cyclohexane", "C1CCCCC1"),
        ("ethanol", "CCO"),
    ];

    let fingerprint = |smiles: &str| -> Result<_, Box<dyn std::error::Error>> {
        let molecule = parse_smiles(smiles)?;
        Ok(MorganFingerprint::get_fingerprint_as_bitvec(
            &molecule, RADIUS, BITS, None, None, false, true, false,
        )?)
    };

    let query = fingerprint("c1ccccc1O")?;

    let mut ranked = Vec::new();
    for (name, smiles) in library {
        let similarity = tanimoto_similarity(&query, &fingerprint(smiles)?)?;
        ranked.push((similarity, name));
    }

    // Descending, and stable — so molecules that tie keep the order the library
    // listed them in rather than an arbitrary one.
    ranked.sort_by(|a, b| b.0.partial_cmp(&a.0).expect("similarities are never NaN"));

    println!("query: phenol (c1ccccc1O), Morgan radius {RADIUS}, {BITS} bits\n");
    for (similarity, name) in ranked {
        println!("  {similarity:.3}  {name}");
    }
    Ok(())
}
