//! SMILES in, SVG out.
//!
//! ```sh
//! cargo run --example smiles_to_svg -- 'CC(=O)Oc1ccccc1C(=O)O' aspirin.svg
//! ```
//!
//! The whole depiction pipeline is three steps: parse, lay out, describe. Only
//! the last one knows about drawing, which is why an SVG writer and a GUI
//! painter can share it.

use chem::core::layout::ensure_coords;
use chem::draw::Vec2;
use chem::draw::structure::{StructureOptions, StructureTheme};
use chem::draw::svg::structure_to_svg;
use chem::io::smiles::parse_smiles;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let smiles = args.next().unwrap_or_else(|| "c1ccccc1O".to_string());
    let path = args.next().unwrap_or_else(|| "structure.svg".to_string());

    let mut molecule = parse_smiles(&smiles)?;

    // SMILES carries no coordinates, so a molecule parsed from one has no
    // geometry until something computes it. Drawing without this gives a
    // placeholder rather than a structure.
    if !ensure_coords(&mut molecule) {
        return Err(format!("could not lay out {smiles}").into());
    }

    // The theme is explicit rather than ambient: an SVG bound for a document
    // should not carry a dark background's colours because of where it was
    // generated.
    let svg = structure_to_svg(
        &molecule,
        Vec2::new(360.0, 300.0),
        &StructureOptions::default(),
        &StructureTheme::light(),
    );

    std::fs::write(&path, &svg)?;
    println!("{smiles} -> {path} ({} bytes)", svg.len());
    Ok(())
}
