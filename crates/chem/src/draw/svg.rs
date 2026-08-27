//! Writing a structure out as SVG.
//!
//! A second consumer of [`describe_structure`], which is what made this possible
//! — before it, the renderer computed geometry and emitted egui shapes in one
//! pass, with nothing in between to serialise.
//!
//! # Text
//!
//! Labels are emitted as `<text>` naming a font stack, and the bond insets are
//! computed from an *estimate* for that stack rather than from egui's measured
//! metrics. That matters: bonds stop at the edge of their labels, so insets
//! measured with one font and drawn with another leave bonds striking through
//! text or stopping short of nothing. Estimating for the font the file declares
//! keeps the file internally consistent — which is the property worth having,
//! since it cannot be pixel-identical to the screen anyway.
//!
//! The residue is that a viewer substituting a metrically different face shows
//! slight drift. Converting glyphs to outlines would remove even that, at the
//! cost of a font-parsing dependency and labels that are no longer text.

use crate::core::molecule::Molecule;
use crate::draw::structure::{
    MeasureText, StructureOptions, StructureShape, StructureTheme, describe_structure,
};
use ecolor::Color32;
use emath::{Align2, Rect, Vec2, pos2};

/// The font the file asks for, most-preferred first.
///
/// A stack rather than one name, so a viewer without the first has somewhere to
/// fall back to before reaching its own default.
const FONT_STACK: &str = "Helvetica, Arial, 'DejaVu Sans', sans-serif";

/// Width of a glyph as a fraction of the font size, averaged over the
/// characters that appear in atom labels.
///
/// Element symbols are one or two capitals and digits, which in a humanist
/// sans sit near 0.6 em; this is not a substitute for real metrics, and does not
/// need to be, because the same number decides both where a label is drawn and
/// how far the bond stops short of it. Being wrong about the font makes the
/// drawing slightly loose, not inconsistent.
const GLYPH_ASPECT: f32 = 0.6;

/// Cap height as a fraction of the font size, for centring a label vertically.
const CAP_HEIGHT: f32 = 0.72;

/// Estimates a label's box for [`FONT_STACK`] at a size.
fn estimate(text: &str, size: f32) -> Vec2 {
    Vec2::new(text.chars().count() as f32 * size * GLYPH_ASPECT, size)
}

/// Renders a molecule as a standalone SVG document.
///
/// `size` is the viewport in points, the same units the on-screen renderer uses,
/// so a structure exported at the size it was displayed comes out at the same
/// proportions.
pub fn structure_to_svg(
    molecule: &Molecule,
    size: Vec2,
    options: &StructureOptions,
    theme: &StructureTheme,
) -> String {
    let rect = Rect::from_min_size(pos2(0.0, 0.0), size);
    let measure: MeasureText<'_> = &estimate;
    let shapes = describe_structure(molecule, rect, options, theme, theme.foreground, measure);

    let mut out = String::with_capacity(1024 + shapes.len() * 96);
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" \
         viewBox=\"0 0 {:.0} {:.0}\">\n",
        size.x, size.y, size.x, size.y
    ));
    // No background rect: a structure dropped into a document should take the
    // document's paper, not carry its own.
    out.push_str("  <g fill=\"none\" stroke-linecap=\"round\">\n");

    for shape in &shapes {
        match shape {
            StructureShape::Line {
                from,
                to,
                width,
                color,
            } => {
                out.push_str(&format!(
                    "    <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                     stroke=\"{}\" stroke-width=\"{:.2}\"{}/>\n",
                    from.x,
                    from.y,
                    to.x,
                    to.y,
                    hex(*color),
                    width,
                    opacity(*color)
                ));
            }
            StructureShape::DashedLine {
                from,
                to,
                width,
                color,
                dash,
            } => {
                out.push_str(&format!(
                    "    <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                     stroke=\"{}\" stroke-width=\"{:.2}\" stroke-dasharray=\"{:.2} {:.2}\"{}/>\n",
                    from.x,
                    from.y,
                    to.x,
                    to.y,
                    hex(*color),
                    width,
                    dash,
                    dash,
                    opacity(*color)
                ));
            }
            StructureShape::Disc {
                center,
                radius,
                color,
            } => {
                out.push_str(&format!(
                    "    <circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"{}\"{}/>\n",
                    center.x,
                    center.y,
                    radius,
                    hex(*color),
                    opacity(*color)
                ));
            }
            StructureShape::Text {
                pos,
                align,
                text,
                size,
                color,
            } => {
                // SVG anchors text on a baseline, so a centred label has to be
                // shifted down by half its cap height. `dominant-baseline`
                // would say this more directly but is unevenly supported, and
                // an export is no place to rely on that.
                let (anchor, dy) = match *align {
                    Align2::CENTER_CENTER => ("middle", size * CAP_HEIGHT / 2.0),
                    _ => ("start", size * CAP_HEIGHT / 2.0),
                };
                out.push_str(&format!(
                    "    <text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"{}\" \
                     font-family=\"{}\" font-size=\"{:.2}\" fill=\"{}\"{}>{}</text>\n",
                    pos.x,
                    pos.y + dy,
                    anchor,
                    FONT_STACK,
                    size,
                    hex(*color),
                    opacity(*color),
                    escape(text)
                ));
            }
        }
    }

    out.push_str("  </g>\n</svg>\n");
    out
}

fn hex(color: Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", color.r(), color.g(), color.b())
}

/// Alpha as a separate attribute, since `#rrggbbaa` is CSS Color 4 and older
/// viewers ignore it — silently drawing a transparent shape opaque.
fn opacity(color: Color32) -> String {
    if color.a() == 255 {
        String::new()
    } else {
        format!(" opacity=\"{:.3}\"", color.a() as f32 / 255.0)
    }
}

/// Escapes the five characters that cannot appear literally in XML content.
///
/// Atom labels are element symbols and digits today, so nothing here needs it —
/// but an export that produces invalid XML on some future label is a bug waiting
/// rather than a risk worth taking.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

/// A filename for a molecule, safe to write or to offer to a save dialog.
///
/// Molecule names come from a file's own records and can hold anything; a slash
/// would silently redirect where the file lands, and a leading dot would hide
/// it. Sanitising here rather than at each call site means the app's dialog and
/// a command-line `-o` get the same guarantee.
pub fn suggested_filename(molecule_name: &str) -> String {
    let cleaned: String = molecule_name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    let trimmed = cleaned.trim_matches('_');
    if trimmed.is_empty() {
        "structure.svg".to_owned()
    } else {
        format!("{trimmed}.svg")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::layout::layout;
    use crate::core::prelude::*;
    use emath::vec2;

    fn ring(order: BondOrder) -> Molecule {
        let mut mol = Molecule::new();
        for _ in 0..6 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)] {
            mol.add_bond(Bond::new(a, b, order)).unwrap();
        }
        layout(&mut mol);
        mol
    }

    fn render(mol: &Molecule) -> String {
        structure_to_svg(
            mol,
            vec2(200.0, 160.0),
            &StructureOptions::default(),
            &StructureTheme::light(),
        )
    }

    #[test]
    fn test_document_declares_its_size_and_namespace() {
        let svg = render(&ring(BondOrder::Single));

        assert!(svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\""));
        assert!(svg.contains("width=\"200\" height=\"160\""));
        assert!(svg.contains("viewBox=\"0 0 200 160\""));
        assert!(svg.trim_end().ends_with("</svg>"));
    }

    #[test]
    fn test_a_hexagon_is_six_lines() {
        let svg = render(&ring(BondOrder::Single));
        assert_eq!(svg.matches("<line ").count(), 6);
        // No carbon is labelled, so a plain ring carries no text at all.
        assert_eq!(svg.matches("<text ").count(), 0);
    }

    #[test]
    fn test_aromatic_bonds_become_dashed_strokes() {
        let svg = render(&ring(BondOrder::Aromatic));

        // Six solid perimeter lines and six dashed inner ones. The dash length
        // survives as an attribute rather than being cut into segments, which is
        // the point of carrying it through the description.
        assert_eq!(svg.matches("<line ").count(), 12);
        assert_eq!(svg.matches("stroke-dasharray=").count(), 6);
    }

    #[test]
    fn test_a_heteroatom_becomes_text_in_its_element_colour() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));
        mol.add_bond(Bond::new(0, 1, BondOrder::Single)).unwrap();
        layout(&mut mol);

        let svg = render(&mol);
        assert!(svg.contains(">O</text>"));
        assert!(svg.contains("text-anchor=\"middle\""));
        // Oxygen is drawn red by the light theme, and the export must not fall
        // back to the live theme or to black.
        let oxygen = hex(StructureTheme::light().oxygen);
        assert!(svg.contains(&format!("fill=\"{oxygen}\"")), "{svg}");
    }

    #[test]
    fn test_no_background_is_drawn() {
        // A structure dropped into a document should take the document's paper.
        // A white rect would be invisible on white and wrong everywhere else.
        let svg = render(&ring(BondOrder::Single));
        assert!(!svg.contains("<rect"));
    }

    #[test]
    fn test_a_molecule_without_coordinates_exports_a_message_not_an_empty_file() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));

        let svg = render(&mol);
        assert!(svg.contains("No 2D coordinates"));
    }

    #[test]
    fn test_xml_special_characters_are_escaped() {
        // Element symbols and digits need none of this today. An export that
        // emits invalid XML on some future label is a bug waiting rather than a
        // risk worth running.
        assert_eq!(escape("A&B"), "A&amp;B");
        assert_eq!(escape("<x>"), "&lt;x&gt;");
        assert_eq!(escape("q\"'"), "q&quot;&apos;");
        assert_eq!(escape("O"), "O");
    }

    #[test]
    fn test_transparency_is_a_separate_attribute() {
        // `#rrggbbaa` is CSS Color 4; an older viewer ignores the alpha and
        // draws the shape opaque, which is worse than a slightly longer file.
        assert_eq!(opacity(Color32::from_rgb(1, 2, 3)), "");
        assert!(opacity(Color32::from_rgba_unmultiplied(1, 2, 3, 128)).contains("opacity="));
        assert_eq!(hex(Color32::from_rgb(0x12, 0x34, 0x56)), "#123456");
    }

    #[test]
    fn test_a_filename_cannot_escape_the_dialog() {
        // Dataset names come from a file's own records and can hold anything. A
        // slash would silently redirect the save; a leading dot would hide the
        // file.
        assert_eq!(suggested_filename("Phenol"), "Phenol.svg");
        assert_eq!(suggested_filename("../../etc/passwd"), "etc_passwd.svg");
        assert_eq!(suggested_filename(".hidden"), "hidden.svg");
        assert_eq!(suggested_filename("c1ccccc1O"), "c1ccccc1O.svg");
        assert_eq!(suggested_filename("mol 1 (a)"), "mol_1__a.svg");
    }

    #[test]
    fn test_an_unnameable_molecule_still_gets_a_filename() {
        // Everything stripped leaves nothing to build a name from, and offering
        // ".svg" or an empty box to a save dialog is worse than a generic name.
        assert_eq!(suggested_filename(""), "structure.svg");
        assert_eq!(suggested_filename("///"), "structure.svg");
        assert_eq!(suggested_filename("___"), "structure.svg");
    }

    #[test]
    fn test_the_estimate_scales_with_text_and_size() {
        // The insets and the label positions come from this same number, so
        // being wrong about the font loosens the drawing rather than making it
        // inconsistent.
        assert!(estimate("Cl", 10.0).x > estimate("C", 10.0).x);
        assert!(estimate("C", 20.0).x > estimate("C", 10.0).x);
        assert_eq!(estimate("C", 12.0).y, 12.0);
    }
}
