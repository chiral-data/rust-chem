//! Rendering structures to SVG files.
//!
//! The depiction itself is `chemdraw`; what lives here is the part a directory
//! of files needs and a single export does not — turning molecule names into
//! filenames that cannot collide, escape, or overwrite each other.

use anyhow::{Context, Result};
use chem::draw::svg::suggested_filename;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Chooses a unique filename for every molecule.
///
/// [`suggested_filename`] sanitises one name in isolation, which is all the
/// app's save dialog needs. A directory export needs more: molecule names come
/// from a file's own records and are not unique, so two records called
/// `Phenol` — or `Phenol` and `Phenol!`, which sanitise to the same thing —
/// would silently write to one file and leave the export one structure short.
///
/// Duplicates get `-2`, `-3` and so on, in the order the file listed them, so a
/// rerun over the same input produces the same names.
pub fn unique_filenames(names: &[String]) -> Vec<String> {
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut out = Vec::with_capacity(names.len());

    for name in names {
        let base = suggested_filename(name);
        let count = seen.entry(base.clone()).or_insert(0);
        *count += 1;
        if *count == 1 {
            out.push(base);
        } else {
            // Split on the extension rather than appending, so the result is
            // still `something.svg` and not `something.svg-2`.
            let stem = base.strip_suffix(".svg").unwrap_or(&base);
            out.push(format!("{stem}-{count}.svg"));
        }
    }
    out
}

/// Writes one SVG per molecule into `dir`, creating it if needed.
pub fn write_directory(dir: &Path, files: &[(String, String)]) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    for (name, svg) in files {
        let path: PathBuf = dir.join(name);
        std::fs::write(&path, svg).with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distinct_names_are_left_alone() {
        let names = vec!["Phenol".to_string(), "Benzene".to_string()];
        assert_eq!(unique_filenames(&names), ["Phenol.svg", "Benzene.svg"]);
    }

    #[test]
    fn test_duplicate_names_do_not_overwrite_each_other() {
        // The failure this prevents is silent: without it the export finishes
        // reporting three files and the directory holds one.
        let names = vec!["Phenol".to_string(); 3];
        assert_eq!(
            unique_filenames(&names),
            ["Phenol.svg", "Phenol-2.svg", "Phenol-3.svg"]
        );
    }

    #[test]
    fn test_names_that_sanitise_to_the_same_thing_also_collide() {
        // `Phenol!` and `Phenol/` both become `Phenol`, so uniqueness has to be
        // decided after sanitising rather than before.
        let names = vec![
            "Phenol".to_string(),
            "Phenol!".to_string(),
            "Phenol/".to_string(),
        ];
        assert_eq!(
            unique_filenames(&names),
            ["Phenol.svg", "Phenol-2.svg", "Phenol-3.svg"]
        );
    }

    #[test]
    fn test_the_suffix_keeps_the_extension() {
        let names = vec!["x".to_string(), "x".to_string()];
        assert!(unique_filenames(&names).iter().all(|n| n.ends_with(".svg")));
    }

    #[test]
    fn test_unnameable_molecules_still_get_distinct_files() {
        // `suggested_filename` maps anything unusable to `structure.svg`, so a
        // file of unnamed records would otherwise be one file.
        let names = vec!["".to_string(), "!!!".to_string(), "///".to_string()];
        assert_eq!(
            unique_filenames(&names),
            ["structure.svg", "structure-2.svg", "structure-3.svg"]
        );
    }

    #[test]
    fn test_a_name_cannot_escape_the_output_directory() {
        // Molecule names come from a file's own records, so this is untrusted
        // input reaching a filesystem path.
        for hostile in ["../../etc/passwd", "/etc/passwd", "..", "a/b"] {
            let picked = &unique_filenames(&[hostile.to_string()])[0];
            assert!(!picked.contains('/'), "{hostile:?} -> {picked:?}");
            assert!(!picked.contains(".."), "{hostile:?} -> {picked:?}");
        }
    }

    #[test]
    fn test_ordering_is_stable_across_runs() {
        let names: Vec<String> = ["a", "b", "a", "c", "a"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(unique_filenames(&names), unique_filenames(&names));
        assert_eq!(
            unique_filenames(&names),
            ["a.svg", "b.svg", "a-2.svg", "c.svg", "a-3.svg"]
        );
    }
}
