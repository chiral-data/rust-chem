//! The committed corpus, checked without a toolkit in sight.
//!
//! The oracle harness in `tools/oracle/` needs Docker, RDKit and OpenBabel.
//! This does not: it reads the same files and asserts what they promise, so the
//! corpus stays meaningful in an ordinary `cargo test` run.
//!
//! Without this the fixtures would only matter inside a container, and would
//! rot the first time someone changed the parser without building the image.

use std::fs;
use std::path::{Path, PathBuf};

use chem::io::format::Format;
use chem::io::reader;

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus")
}

/// The `SMILES<space>name` lines of a corpus file, comments and blanks dropped.
fn records(path: &Path) -> Vec<(String, String)> {
    let text = fs::read_to_string(path).expect("corpus file is readable");
    text.lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let mut parts = line.splitn(2, char::is_whitespace);
            let smiles = parts.next().unwrap_or_default().to_string();
            let name = parts.next().unwrap_or_default().trim().to_string();
            (smiles, name)
        })
        .collect()
}

/// Whether this entry is a pinned gap.
///
/// Marked by the entry's *name* rather than a nearby comment. The first
/// version of this scanned for a `KNOWN GAP` comment and armed itself on the
/// phrase appearing in the file's own header, then swept up every line below
/// it. A name carries the marker into every listing and cannot mis-arm.
///
/// Gaps run in both directions: input that should be rejected and is not, and
/// valid input the parser cannot yet read. Both are pinned so the corpus
/// describes the parser as it behaves, and closing either fails the test that
/// pins it — which is the prompt to delete the marker.
fn is_known_gap(name: &str) -> bool {
    name.starts_with("known-gap-")
}

fn themed_files() -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(corpus_dir())
        .expect("corpus directory exists")
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|e| e == "smi"))
        .filter(|path| path.file_stem().is_some_and(|s| s != "invalid"))
        .collect();
    files.sort();
    files
}

#[test]
fn test_every_themed_fixture_parses() {
    // The corpus is only worth having if it is readable. A fixture that stopped
    // parsing would silently drop out of the oracle run too, taking whatever it
    // was probing with it.
    //
    // Except the pinned gaps, which are valid SMILES this parser cannot read
    // yet (#191) and must keep failing until it can.
    let files = themed_files();
    assert!(!files.is_empty(), "no themed corpus files found");

    for path in files {
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        for (smiles, entry) in records(&path) {
            let outcome = reader::read(&format!("{smiles} {entry}\n"), Format::SMILES);
            let parsed = !outcome.records.is_empty();

            if is_known_gap(&entry) {
                assert!(
                    !parsed,
                    "{name}: {entry} is pinned as a gap but now parses — \
                     the gap is closed, so rename it and drop its comment"
                );
            } else {
                assert!(
                    parsed,
                    "{name}: {entry} ({smiles}) must parse; {:?}",
                    outcome
                        .skipped
                        .iter()
                        .map(|s| s.error.clone())
                        .collect::<Vec<_>>()
                );
            }
        }
    }
}

#[test]
fn test_the_rejection_fixture_is_rejected() {
    // The file the parse check exists for. A suite of valid molecules proves
    // what a parser accepts and never what it refuses.
    let path = corpus_dir().join("invalid.smi");
    let text = fs::read_to_string(&path).expect("readable");
    let outcome = reader::read(&text, Format::SMILES);

    let accepted: Vec<&str> = outcome
        .records
        .iter()
        .map(|record| record.name.as_str())
        .collect();

    // Anything that parsed must be a pinned gap, and nothing else.
    for name in &accepted {
        assert!(
            is_known_gap(name),
            "invalid.smi: {name} was accepted but is not pinned as a gap — \
             either it is valid and belongs in a themed file, or the parser is wrong"
        );
    }

    // And every pinned gap must still be accepted. When one starts being
    // rejected the parser was fixed, and the marker has to go with it.
    for (_, name) in records(&path) {
        if is_known_gap(&name) {
            assert!(
                accepted.contains(&name.as_str()),
                "invalid.smi: {name} is pinned as a gap but is now rejected — \
                 the gap is closed, so rename it and drop its comment"
            );
        }
    }

    assert!(
        !outcome.skipped.is_empty(),
        "invalid.smi rejected nothing at all, which cannot be right"
    );
}

#[test]
fn test_promoted_regressions_still_parse() {
    // Promoted fixtures are the harness's whole output. If the directory exists
    // its contents have to stay loadable, or a finding silently stops being
    // checked.
    let dir = corpus_dir().join("regressions");
    if !dir.exists() {
        return; // nothing promoted yet
    }
    for entry in fs::read_dir(&dir).expect("readable").filter_map(Result::ok) {
        let path = entry.path();
        if path.extension().is_some_and(|e| e == "smi") {
            let text = fs::read_to_string(&path).expect("readable");
            let outcome = reader::read(&text, Format::SMILES);
            assert!(
                !outcome.records.is_empty(),
                "{:?} promoted nothing readable",
                path.file_name().unwrap()
            );
        }
    }
}

/// One bond end: chain, residue number, atom name.
type NamedAtom = (String, i32, String);

/// Each bond as its two (chain, residue number, atom name) ends, sorted, so
/// two formats' reads compare by name rather than by atom index.
fn named_bonds(mol: &chem::core::molecule::Molecule) -> Vec<(NamedAtom, NamedAtom)> {
    let end = |i: usize| {
        let residue = mol.residue_of(i).expect("every atom has a residue");
        let chain = mol.chain_of(i).expect("every atom has a chain");
        let name = mol.site(i).and_then(|s| s.name.clone()).unwrap_or_default();
        (chain.id.clone(), residue.sequence, name)
    };
    let mut bonds: Vec<_> = mol
        .bonds()
        .iter()
        .map(|b| {
            let (a, c) = (end(b.atom1()), end(b.atom2()));
            if a <= c { (a, c) } else { (c, a) }
        })
        .collect();
    bonds.sort();
    bonds
}

#[test]
fn test_crambin_gets_the_same_template_bonds_from_pdb_and_mmcif() {
    let pdb_text = fs::read_to_string(corpus_dir().join("pdb/crambin-1crn.pdb")).unwrap();
    let cif_text = fs::read_to_string(corpus_dir().join("mmcif/crambin-1crn.cif")).unwrap();
    let pdb = chem::io::pdb::parse_pdb_with_options(
        &pdb_text,
        &chem::io::options::PdbReadOptions {
            residue_templates: true,
            ..Default::default()
        },
    )
    .unwrap();
    let cif = chem::io::mmcif::parse_mmcif_with_options(
        &cif_text,
        &chem::io::options::MmcifReadOptions {
            residue_templates: true,
            ..Default::default()
        },
    )
    .unwrap();

    // 327 heavy atoms in one chain, plus one ring closure each for 5 PRO,
    // 1 PHE, 2 TYR and 3 disulfides: fully connected, no break.
    assert_eq!(pdb.num_atoms(), 327);
    assert_eq!(pdb.num_bonds(), 326 + 11);
    assert_eq!(named_bonds(&pdb), named_bonds(&cif));

    let disulfides = named_bonds(&cif)
        .into_iter()
        .filter(|(a, b)| a.2 == "SG" && b.2 == "SG")
        .count();
    assert_eq!(disulfides, 3);

    // Hydrogens fill from the template bonds: crambin's own formula.
    assert_eq!(pdb.formula(), "C202H315N55O64S6");
    assert_eq!(cif.formula(), pdb.formula());

    // The default stays the file's own statement: CONECT for PDB, nothing for mmCIF.
    assert_eq!(chem::io::pdb::parse_pdb(&pdb_text).unwrap().num_bonds(), 3);
    assert_eq!(
        chem::io::mmcif::parse_mmcif(&cif_text).unwrap().num_bonds(),
        0
    );
}

#[test]
fn test_templates_and_chem_comp_bonds_together_read_like_templates_alone() {
    // 1CRN's own `_chem_comp_bond` table repeats the templates' intra-residue
    // bonds, so turning it on as well must change nothing -- in particular
    // the hydrogen fill must run once, after every bond step.
    let text = fs::read_to_string(corpus_dir().join("mmcif/crambin-1crn.cif")).unwrap();
    let read = |chem_comp_bonds| {
        chem::io::mmcif::parse_mmcif_with_options(
            &text,
            &chem::io::options::MmcifReadOptions {
                residue_templates: true,
                chem_comp_bonds,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let (alone, both) = (read(false), read(true));
    assert_eq!(named_bonds(&both), named_bonds(&alone));
    let orders = |m: &chem::prelude::Molecule| {
        let mut v: Vec<_> = m
            .bonds()
            .iter()
            .map(|b| {
                (
                    b.atom1().min(b.atom2()),
                    b.atom1().max(b.atom2()),
                    b.order(),
                )
            })
            .collect();
        v.sort_by_key(|&(a, b, _)| (a, b));
        v
    };
    assert_eq!(orders(&both), orders(&alone));
    assert_eq!(both.formula(), "C202H315N55O64S6");
}
