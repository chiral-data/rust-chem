use std::fs;
use std::path::Path;

const README: &str = include_str!("../README.md");

/// Regression test for #343: the README described a version three releases
/// stale, twice, because nothing checked it against the crate's own manifest.
#[test]
fn test_the_readme_version_matches_the_crate_version() {
    let version = env!("CARGO_PKG_VERSION");
    // `chem = "0.8"` style snippets state only major.minor, not the patch —
    // matching that prefix is the README's actual promise.
    let major_minor = version.rsplit_once('.').map_or(version, |(mm, _)| mm);

    let occurrences = README.matches(major_minor).count();
    assert!(
        occurrences >= 2,
        "expected at least 2 version snippets in README.md matching \
         crate version {version} (prefix {major_minor}), found {occurrences} — \
         either the README fell behind Cargo.toml, or a snippet was removed \
         and this count should shrink with it"
    );
    assert!(
        !README.contains("= \"0.6\"") && !README.contains("= \"0.7\""),
        "README.md still names a version older than the crate's own {version}"
    );
}

/// Regression test for #343: two of the crate's eight examples were never
/// added to the README's table, so `cargo run --example format_registry`
/// worked but nobody reading docs.rs would know it existed.
#[test]
fn test_every_example_is_mentioned_in_the_readme() {
    let examples_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let mut missing = Vec::new();

    for entry in fs::read_dir(&examples_dir).expect("examples dir exists") {
        let path = entry.expect("readable dir entry").path();
        if path.extension().and_then(|e| e.to_str()) != Some("rs") {
            continue;
        }
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .expect("utf-8 example filename")
            .to_string();
        if !README.contains(&name) {
            missing.push(name);
        }
    }

    assert!(
        missing.is_empty(),
        "these examples exist under crates/chem/examples/ but are not named \
         anywhere in README.md: {missing:?}"
    );
}
