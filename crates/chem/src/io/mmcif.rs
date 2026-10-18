//! mmCIF — tag-value text (`_category.item value`), with `loop_` blocks for
//! repeated records. The modern PDB successor (#224).
//!
//! Data-model-wise this is #223 (PDB) again: [`crate::core::site::AtomSite`],
//! [`crate::core::residue::Chain`]/[`crate::core::residue::Residue`]
//! (already carrying mmCIF's own dual `auth_*`/`label_*` numbering — built
//! for exactly this), [`crate::core::cell::UnitCell`]/
//! [`crate::core::cell::SpaceGroup`] all already exist and need no changes.
//! What's new is the *grammar*, not the model. The actual atom-site/chain/
//! cell interpretation now lives in [`crate::io::cif_model`], shared with
//! BinaryCIF (#319) — this module is the text-specific tokenizer/formatter
//! around it.
//!
//! **No bonds by default.** With
//! [`MmcifReadOptions::residue_templates`](crate::io::options::MmcifReadOptions::residue_templates)
//! on (#403), standard residues are bonded from
//! [`crate::core::residue_templates`], then the covalent rows of
//! `_struct_conn` (`disulf`, `covale*`), each partner matched by its
//! `label_*` IDs or, failing that, its `auth_*` IDs. Hydrogen bonds, metal
//! coordination, links to a symmetry mate and unresolvable partners are
//! skipped. BinaryCIF shares all of this through [`crate::io::cif_model`].
//!
//! With
//! [`MmcifReadOptions::chem_comp_bonds`](crate::io::options::MmcifReadOptions::chem_comp_bonds)
//! on, each residue is bonded by the block's own `_chem_comp_bond` table,
//! matched on atom name, with the stated orders (#406); atoms the table
//! names but the residue lacks are skipped. BinaryCIF reads the same way.
//! Neither writer emits `_chem_comp_bond`, so bonds don't survive a write.
//!
//! **Grammar subset**: bare (whitespace-delimited) tokens, single- and
//! double-quoted tokens. Multi-line semicolon (`;...;`) text fields are
//! **not** parsed — real files never need one inside `_atom_site`/`_cell`/
//! `_symmetry` data, only for free-text descriptions this module doesn't
//! read, so the limitation costs nothing in practice.
//!
//! **Multiple `data_` blocks** in one file each become a separate record —
//! splitting is the reader's/supplier's job
//! ([`crate::io::reader::read_mmcif_with_options`],
//! [`crate::io::supplier::MmcifSupplier`]), the same division every prior
//! format has. **A varying `pdbx_PDB_model_num` within one block's
//! `_atom_site` loop** is mmCIF's NMR-ensemble convention: models share one
//! loop rather than living in separate blocks. [`parse_mmcif`] returns the
//! first model. The reader returns the whole ensemble as one
//! [`crate::io::reader::Payload::Frames`] record under
//! [`crate::io::options::MultiFrameMode::Frames`], and otherwise reports
//! each later model as skipped (#404). A later model whose atoms differ
//! from the first model's is always skipped, never merged.

use crate::core::molecule::Molecule;
use crate::io::cif_model::{
    AtomSiteRows, CifLoops, CifModels, RowPrecision, build_models, build_rows, category_of,
};
use crate::io::errors::MmcifError;
use crate::io::options::MmcifReadOptions;

/// Splits one line into its whitespace-delimited or quoted tokens. A quote
/// is a delimiter only where it starts a token (preceded by whitespace or
/// nothing) and ends one (followed by whitespace or nothing) — the
/// standard CIF quoting rule, so `O5'` (a common atom name containing an
/// apostrophe) is not mistaken for the start of a quoted string.
///
/// `pub(crate)`: the CIF grammar this tokenizes is shared with
/// [`crate::io::cif_core`] (#320), which has no dictionary-specific
/// knowledge baked into it.
pub(crate) fn tokenize_line(line: &str) -> Vec<String> {
    let chars: Vec<char> = line.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        while i < chars.len() && chars[i].is_whitespace() {
            i += 1;
        }
        if i >= chars.len() {
            break;
        }
        if chars[i] == '\'' || chars[i] == '"' {
            let quote = chars[i];
            i += 1;
            let start = i;
            while i < chars.len()
                && !(chars[i] == quote && chars.get(i + 1).is_none_or(|c| c.is_whitespace()))
            {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
            i += 1; // the closing quote
        } else {
            let start = i;
            while i < chars.len() && !chars[i].is_whitespace() {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
        }
    }
    tokens
}

/// Parses one mmCIF `data_` block (or a whole file with only one).
pub fn parse_mmcif(text: &str) -> Result<Molecule, MmcifError> {
    parse_mmcif_with_options(text, &MmcifReadOptions::default())
}

/// [`parse_mmcif`] with mmCIF's own read options. Returns the first model
/// of an NMR ensemble.
pub fn parse_mmcif_with_options(
    text: &str,
    options: &MmcifReadOptions,
) -> Result<Molecule, MmcifError> {
    parse_mmcif_models(text, options).map(|models| models.molecule)
}

/// One block's first model plus its later models (#404).
pub(crate) fn parse_mmcif_models(
    text: &str,
    options: &MmcifReadOptions,
) -> Result<CifModels, MmcifError> {
    let mut singles: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut atom_site: Option<AtomSiteRows> = None;
    let mut loops = CifLoops::new();

    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        let trimmed = lines[i].trim();
        let lower = trimmed.to_ascii_lowercase();
        if trimmed.is_empty() || trimmed.starts_with('#') || lower.starts_with("data_") {
            i += 1;
            continue;
        }

        if lower == "loop_" {
            i += 1;
            let mut tags = Vec::new();
            while i < lines.len() && lines[i].trim_start().starts_with('_') {
                tags.push(lines[i].trim().to_string());
                i += 1;
            }
            let mut rows: Vec<Vec<String>> = Vec::new();
            while i < lines.len() {
                let row_trimmed = lines[i].trim();
                let row_lower = row_trimmed.to_ascii_lowercase();
                if row_trimmed.is_empty()
                    || row_trimmed.starts_with('#')
                    || row_trimmed.starts_with('_')
                    || row_lower == "loop_"
                    || row_lower.starts_with("data_")
                {
                    break;
                }
                rows.push(tokenize_line(lines[i]));
                i += 1;
            }

            if tags.first().is_some_and(|t| t.starts_with("_atom_site.")) {
                atom_site = Some((tags, rows));
            } else if let Some(first) = tags.first() {
                loops.insert(category_of(first), (tags, rows));
            }
            continue;
        }

        if trimmed.starts_with('_') {
            let tokens = tokenize_line(lines[i]);
            if tokens.len() >= 2 {
                singles.insert(tokens[0].clone(), tokens[1].clone());
            }
        }
        i += 1;
    }

    build_models(
        &singles,
        atom_site
            .as_ref()
            .map(|(t, r)| (t.as_slice(), r.as_slice())),
        &loops,
        options,
    )
}

/// Writes one mmCIF `data_` block.
pub fn write_mmcif(mol: &Molecule) -> String {
    let (singles, (tags, rows)) = build_rows(mol, RowPrecision::Text);

    let mut out = String::from("data_chem\n");
    // One block, written together or not at all -- `build_rows` only ever
    // inserts all six `_cell.*` keys as a unit, so checking one is checking
    // all of them.
    if singles.contains_key("_cell.length_a") {
        for key in [
            "_cell.length_a",
            "_cell.length_b",
            "_cell.length_c",
            "_cell.angle_alpha",
            "_cell.angle_beta",
            "_cell.angle_gamma",
        ] {
            out.push_str(&format!("{key} {}\n", singles[key]));
        }
    }
    if let Some(sg) = singles.get("_symmetry.space_group_name_H-M") {
        out.push_str(&format!("_symmetry.space_group_name_H-M '{sg}'\n"));
    }

    out.push_str("loop_\n");
    for tag in &tags {
        out.push_str(tag);
        out.push('\n');
    }
    for row in &rows {
        out.push_str(&row.join(" "));
        out.push('\n');
    }

    out
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::core::geometry::Point3;

    const WATER_MMCIF: &str = "\
data_water
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
HETATM 1 O O . HOH A 1 A HOH . 0.000 0.000 0.000 1.00 20.00 1
HETATM 2 H H1 . HOH A 1 A HOH . 0.759 0.000 0.504 1.00 20.00 1
HETATM 3 H H2 . HOH A 1 A HOH . 0.759 0.000 -0.504 1.00 20.00 1
";

    #[test]
    fn test_a_single_block_round_trips() {
        let mol = parse_mmcif(WATER_MMCIF).expect("valid mmCIF");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.residues().len(), 1);
        assert!(mol.residues()[0].is_hetero);
        assert_eq!(mol.residues()[0].name, "HOH");
        assert_eq!(mol.chains()[0].id, "A");
        assert_eq!(mol.site(0).unwrap().occupancy, Some(1.0));

        let written = write_mmcif(&mol);
        let back = parse_mmcif(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.residues()[0].name, "HOH");
        assert_eq!(back.chains()[0].id, "A");
    }

    #[test]
    fn test_cell_and_space_group_round_trip() {
        let text = format!(
            "data_x\n_cell.length_a 10.000\n_cell.length_b 20.000\n_cell.length_c 30.000\n\
             _cell.angle_alpha 90.00\n_cell.angle_beta 90.00\n_cell.angle_gamma 90.00\n\
             _symmetry.space_group_name_H-M 'P 21 21 21'\n{}",
            &WATER_MMCIF["data_water\n".len()..]
        );
        let mol = parse_mmcif(&text).expect("valid mmCIF");
        let cell = mol.cell().expect("has a cell");
        assert_eq!((cell.a, cell.b, cell.c), (10.0, 20.0, 30.0));
        assert_eq!(
            mol.space_group().and_then(|g| g.symbol.as_deref()),
            Some("P 21 21 21")
        );

        let written = write_mmcif(&mol);
        let back = parse_mmcif(&written).expect("round trips");
        assert_eq!(back.cell(), mol.cell());
    }

    #[test]
    fn test_a_quoted_value_containing_a_space_parses_correctly() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_comp_id\n\
                    _atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
                    C 'LIG AND' 0.0 0.0 0.0\n";
        let mol = parse_mmcif(text).expect("valid mmCIF");
        assert_eq!(mol.residues()[0].name, "LIG AND");
    }

    #[test]
    fn test_question_mark_and_dot_both_mean_absent() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n\
                    _atom_site.B_iso_or_equiv\nC 0.0 0.0 0.0 ? .\n";
        let mol = parse_mmcif(text).expect("valid mmCIF");
        assert_eq!(mol.site(0).unwrap().occupancy, None);
        assert_eq!(mol.site(0).unwrap().b_factor, None);
    }

    #[test]
    fn test_multiple_data_blocks_are_this_modules_callers_job() {
        // parse_mmcif itself only ever sees one already-isolated block --
        // confirm it does not choke on a bare `data_` line at the top,
        // which is all that distinguishes "isolated" from "not".
        let mol = parse_mmcif(WATER_MMCIF).expect("valid mmCIF");
        assert_eq!(mol.num_atoms(), 3);
    }

    #[test]
    fn test_parse_mmcif_returns_the_first_model_of_an_ensemble() {
        // The later model reaches only the reader (#404); `parse_mmcif`'s
        // signature has room for one molecule.
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n\
                    C 0.0 0.0 0.0 1\nC 1.0 1.0 1.0 2\n";
        let mol = parse_mmcif(text).expect("valid mmCIF");
        assert_eq!(mol.num_atoms(), 1);
        assert_eq!(mol.coord3(0), Some(Point3::new(0.0, 0.0, 0.0)));

        let models = parse_mmcif_models(text, &MmcifReadOptions::default()).unwrap();
        assert_eq!(models.later.len(), 1);
        assert_eq!(
            models.later[0].positions,
            Ok(vec![Point3::new(1.0, 1.0, 1.0)])
        );
    }

    const ENSEMBLE_HEADER: &str = "data_x\nloop_\n_atom_site.type_symbol\n\
        _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n\
        _atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
        _atom_site.pdbx_PDB_model_num\n";

    fn ensemble(rows: &[&str]) -> CifModels {
        let text = format!("{ENSEMBLE_HEADER}{}\n", rows.join("\n"));
        parse_mmcif_models(&text, &MmcifReadOptions::default()).expect("valid mmCIF")
    }

    #[test]
    fn test_a_three_model_block_gives_one_molecule_and_two_later_models() {
        let models = ensemble(&[
            "N N GLY A 1 0.0 0.0 0.0 1",
            "C CA GLY A 1 1.5 0.0 0.0 1",
            "N N GLY A 1 0.1 0.0 0.0 2",
            "C CA GLY A 1 1.6 0.0 0.0 2",
            "N N GLY A 1 0.2 0.0 0.0 3",
            "C CA GLY A 1 1.7 0.0 0.0 3",
        ]);
        assert_eq!(models.molecule.num_atoms(), 2);
        assert_eq!(models.molecule.residues().len(), 1);
        let numbers: Vec<&str> = models.later.iter().map(|m| m.number.as_str()).collect();
        assert_eq!(numbers, ["2", "3"]);
        assert_eq!(
            models.later[1].positions,
            Ok(vec![Point3::new(0.2, 0.0, 0.0), Point3::new(1.7, 0.0, 0.0)])
        );
    }

    #[test]
    fn test_a_later_model_with_one_atom_fewer_is_rejected_with_a_reason() {
        let models = ensemble(&[
            "N N GLY A 1 0.0 0.0 0.0 1",
            "C CA GLY A 1 1.5 0.0 0.0 1",
            "N N GLY A 1 0.1 0.0 0.0 2",
        ]);
        assert_eq!(models.molecule.num_atoms(), 2);
        let reason = models.later[0].positions.as_ref().unwrap_err();
        assert_eq!(reason, "1 atoms where model 1 has 2");
    }

    #[test]
    fn test_a_later_model_naming_a_different_atom_is_rejected() {
        // Same count, different atom: what a count-only check would merge.
        let models = ensemble(&[
            "N N GLY A 1 0.0 0.0 0.0 1",
            "C CA GLY A 1 1.5 0.0 0.0 1",
            "N N GLY A 1 0.1 0.0 0.0 2",
            "C CB GLY A 1 1.6 0.0 0.0 2",
        ]);
        let reason = models.later[0].positions.as_ref().unwrap_err();
        assert!(reason.starts_with("atom 2 is C CB GLY A 1"), "{reason}");
    }

    #[test]
    fn test_a_single_model_block_has_no_later_models() {
        let models = ensemble(&["N N GLY A 1 0.0 0.0 0.0 1", "C CA GLY A 1 1.5 0.0 0.0 1"]);
        assert_eq!(models.molecule.num_atoms(), 2);
        assert!(models.later.is_empty());
    }

    #[test]
    fn test_a_malformed_row_is_a_clear_error_not_a_panic() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\nC not-a-number 0.0 0.0\n";
        let err = parse_mmcif(text).unwrap_err();
        assert!(matches!(err, MmcifError::InvalidAtomRow(_)), "{err}");
    }

    #[test]
    fn test_a_second_write_keeps_every_column_aligned() {
        // #260. CIF has no empty token, so a chain whose id is the empty
        // string used to reach the file as *nothing* -- 17 values under a
        // header declaring 19, shifting B_iso_or_equiv into the model-number
        // column. Only the second write showed it, because the empty id came
        // from reading back our own `.`, which is why one round trip in the
        // oracle harness never caught it.
        let mut mol = crate::io::smiles::parse_smiles("CC").expect("valid SMILES");
        let mut site = crate::core::site::AtomSite::empty();
        site.b_factor = Some(42.0);
        mol.set_sites(vec![site, crate::core::site::AtomSite::empty()])
            .expect("one per atom");

        let tags = |text: &str| {
            text.lines()
                .filter(|l| l.trim_start().starts_with("_atom_site."))
                .count()
        };

        let mut current = mol;
        for pass in 1..=3 {
            let text = write_mmcif(&current);
            let columns = tags(&text);
            for row in text.lines().filter(|l| l.starts_with("ATOM")) {
                assert_eq!(
                    row.split_whitespace().count(),
                    columns,
                    "pass {pass}: {row:?} does not fill its {columns} declared tags"
                );
            }
            current = parse_mmcif(&text).expect("our own output reads back");
            assert_eq!(
                current.site(0).and_then(|s| s.b_factor),
                Some(42.0),
                "pass {pass}: the b-factor moved out of its column"
            );
        }
    }

    #[test]
    fn test_an_unrecognised_element_is_a_clear_error() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\nXx 0.0 0.0 0.0\n";
        let err = parse_mmcif(text).unwrap_err();
        assert!(matches!(err, MmcifError::InvalidElement(_)), "{err}");
    }

    #[test]
    fn test_all_zero_coordinates_are_not_read_as_a_conformer() {
        // #270.
        let text = "data_x\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C UNK A 1 A 0.000 0.000 0.000\n\
ATOM 2 O O UNK A 1 A 0.000 0.000 0.000\n";
        let mol = parse_mmcif(text).expect("valid mmCIF");
        assert_eq!(mol.num_atoms(), 2);
        assert!(
            !mol.has_coords3(),
            "all-zero is a placeholder, not a conformer"
        );
    }

    #[test]
    fn test_garbage_text_reports_no_atoms_instead_of_an_empty_molecule() {
        for input in [
            "",
            "data_x\n_cell.length_a 10.000\n",
            "not an mmcif file at all\n",
        ] {
            assert!(
                matches!(parse_mmcif(input), Err(MmcifError::NoAtoms)),
                "{input:?} should report NoAtoms"
            );
        }
    }

    #[test]
    fn test_an_atom_site_loop_with_zero_rows_is_a_legitimately_empty_structure() {
        // Distinct from the garbage case above: this loop header genuinely
        // says "here are atom sites" (#268) -- it just lists none, which is
        // legal mmCIF for a deposited empty structure and must not error.
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\n";
        let mol = parse_mmcif(text).expect("a present-but-empty atom_site loop is legal mmCIF");
        assert_eq!(mol.num_atoms(), 0);
    }

    /// Partners named by `label_*` IDs, one ligand only by `auth_*` (its
    /// `label_seq_id` is `.`), and three rows that must not bond: metal
    /// coordination, a hydrogen bond, and a link to a symmetry mate.
    pub(crate) const LINKED_MMCIF: &str = "\
data_links
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 S SG . CYS A 1 3 A 0.0 0.0 0.0
ATOM 2 S SG . CYS A 2 40 A 2.0 0.0 0.0
ATOM 3 N ND2 . ASN A 3 60 A 9.0 0.0 0.0
ATOM 4 N NE2 . HIS A 4 50 A 5.0 0.0 0.0
HETATM 5 C C1 . NAG B . 102 A 10.4 0.0 0.0
HETATM 6 Zn ZN . ZN C . 101 A 7.0 0.0 0.0
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_symmetry
disulf1 disulf A 1 SG A 3 A 2 SG A 40 1_555 1_555
covale1 covale A 3 ND2 A 60 B . C1 A 102 1_555 1_555
metalc1 metalc A 4 NE2 A 50 C . ZN A 101 1_555 1_555
hydrog1 hydrog A 1 SG A 3 A 4 NE2 A 50 1_555 1_555
covale2 covale A 1 SG A 3 A 3 ND2 A 60 1_555 2_555
";

    #[test]
    fn test_struct_conn_bonds_only_covalent_rows_inside_the_model() {
        let options = MmcifReadOptions {
            residue_templates: true,
            ..Default::default()
        };
        let mol = parse_mmcif_with_options(LINKED_MMCIF, &options).expect("reads");
        let bonded = |a: usize, b: usize| mol.graph().has_edge(a, b);
        assert!(bonded(0, 1), "disulf, by label IDs");
        assert!(bonded(2, 4), "covale, falling back to auth IDs");
        assert_eq!(
            mol.num_bonds(),
            2,
            "metalc, hydrog and a symmetry mate stay out"
        );
        assert_eq!(mol.atom(0).hydrogens(), Some(1));
        assert_eq!(mol.atom(3).hydrogens(), None, "bondless stays unknown");
    }

    #[test]
    fn test_a_single_struct_conn_row_written_without_a_loop_still_bonds() {
        let (atoms, _) = LINKED_MMCIF.split_once("loop_\n_struct_conn").unwrap();
        let text = format!(
            "{atoms}_struct_conn.id disulf1\n_struct_conn.conn_type_id disulf\n\
             _struct_conn.ptnr1_label_asym_id A\n_struct_conn.ptnr1_label_seq_id 1\n\
             _struct_conn.ptnr1_label_atom_id SG\n_struct_conn.ptnr2_label_asym_id A\n\
             _struct_conn.ptnr2_label_seq_id 2\n_struct_conn.ptnr2_label_atom_id SG\n"
        );
        let options = MmcifReadOptions {
            residue_templates: true,
            ..Default::default()
        };
        let mol = parse_mmcif_with_options(&text, &options).expect("reads");
        assert!(mol.graph().has_edge(0, 1));
        assert_eq!(mol.num_bonds(), 1);
    }

    #[test]
    fn test_struct_conn_is_ignored_by_default() {
        let mol = parse_mmcif(LINKED_MMCIF).expect("reads");
        assert_eq!(mol.num_bonds(), 0);
        assert!(mol.atoms().iter().all(|a| a.hydrogens().is_none()));
    }
}
