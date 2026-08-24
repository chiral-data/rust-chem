//! Behavioural tests for the `chem` binary.
//!
//! Run as a process rather than by calling `run()` directly, because what these
//! assert *is* the process boundary: which stream a line went to, and what the
//! exit code was. A unit test cannot observe either.
//!
//! `CARGO_BIN_EXE_chem` is set by cargo for integration tests of a crate with a
//! `[[bin]]`, so this needs no test-harness dependency.

use std::io::Write;
use std::process::{Command, Stdio};

const EXE: &str = env!("CARGO_BIN_EXE_chem");

struct Run {
    code: i32,
    stdout: String,
    stderr: String,
}

fn run(args: &[&str], stdin: Option<&str>) -> Run {
    let mut child = Command::new(EXE)
        .args(args)
        .stdin(if stdin.is_some() {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn chem");

    if let Some(text) = stdin {
        child
            .stdin
            .as_mut()
            .expect("piped stdin")
            .write_all(text.as_bytes())
            .expect("write stdin");
    }

    let out = child.wait_with_output().expect("wait");
    Run {
        code: out.status.code().expect("exited normally"),
        stdout: String::from_utf8_lossy(&out.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
    }
}

/// Writes a fixture beside the test binary, so tests need no fixture directory
/// and cannot collide with each other's files.
fn fixture(name: &str, contents: &str) -> std::path::PathBuf {
    let path = std::env::temp_dir().join(format!("chem-cli-test-{name}"));
    std::fs::write(&path, contents).expect("write fixture");
    path
}

const GOOD: &str = "CCO ethanol\nc1ccccc1 benzene\n# comment\n\nCC(=O)O acetic\n";
const MIXED: &str = "CCO ethanol\nnot-a-smiles bad\nc1ccccc1 benzene\n";
const ALL_BAD: &str = "xxx\nyyy\n";

#[test]
fn test_data_goes_to_stdout_and_everything_else_to_stderr() {
    let path = fixture("streams.smi", GOOD);
    let r = run(&["info", path.to_str().unwrap()], None);

    assert_eq!(r.code, 0);
    // stdout is the table and nothing but: a header plus one row per molecule.
    let lines: Vec<&str> = r.stdout.lines().collect();
    assert_eq!(lines.len(), 4, "header + 3 molecules, got {:?}", lines);
    assert!(lines[0].starts_with("name\t"));
    // This is the property that makes `chem a | chem b` possible at all.
    assert!(
        !r.stdout.contains("read 3 molecules"),
        "progress must not reach stdout"
    );
    assert!(r.stderr.contains("read 3 molecules"));
}

#[test]
fn test_stdin_is_read_when_no_file_is_named() {
    let r = run(&["info"], Some(GOOD));
    assert_eq!(r.code, 0);
    assert_eq!(r.stdout.lines().count(), 4);
    assert!(r.stderr.contains("from -"));
}

#[test]
fn test_a_lone_dash_also_means_stdin() {
    let r = run(&["info", "-"], Some(GOOD));
    assert_eq!(r.code, 0);
    assert_eq!(r.stdout.lines().count(), 4);
}

#[test]
fn test_output_flag_writes_a_file_and_leaves_stdout_empty() {
    let path = fixture("out-src.smi", GOOD);
    let dest = std::env::temp_dir().join("chem-cli-test-out.tsv");
    let _ = std::fs::remove_file(&dest);

    let r = run(
        &["info", path.to_str().unwrap(), "-o", dest.to_str().unwrap()],
        None,
    );

    assert_eq!(r.code, 0);
    assert!(r.stdout.is_empty(), "stdout was {:?}", r.stdout);
    let written = std::fs::read_to_string(&dest).expect("output file");
    assert_eq!(written.lines().count(), 4);
}

#[test]
fn test_a_partly_bad_file_still_succeeds_and_names_what_it_skipped() {
    let path = fixture("mixed.smi", MIXED);
    let r = run(&["info", path.to_str().unwrap()], None);

    // Two good molecules out of three lines: throwing them away would help
    // nobody, so this is a success that reports a loss.
    assert_eq!(r.code, 0);
    assert_eq!(r.stdout.lines().count(), 3);
    assert!(r.stderr.contains("skipped line 2"));
    assert!(r.stderr.contains("not-a-smiles"));
}

#[test]
fn test_strict_turns_a_skipped_record_into_a_failure() {
    let path = fixture("strict.smi", MIXED);
    let r = run(&["--strict", "info", path.to_str().unwrap()], None);
    assert_eq!(r.code, 4);
}

#[test]
fn test_strict_is_silent_when_nothing_was_skipped() {
    let path = fixture("strict-clean.smi", GOOD);
    let r = run(&["--strict", "info", path.to_str().unwrap()], None);
    assert_eq!(r.code, 0);
}

#[test]
fn test_an_unreadable_file_is_distinct_from_an_unusable_one() {
    let unusable = fixture("allbad.smi", ALL_BAD);
    let r = run(&["info", unusable.to_str().unwrap()], None);
    assert_eq!(r.code, 2, "nothing parsed");

    let r = run(&["info", "/nonexistent/nope.smi"], None);
    assert_eq!(r.code, 1, "the file itself could not be read");
    assert!(r.stderr.contains("error:"));
    // The cause is kept, not flattened into a one-line summary.
    assert!(r.stderr.contains("caused by:"));
}

#[test]
fn test_format_is_detected_from_the_extension() {
    let path = fixture("detect.sdf", SDF);
    let r = run(&["info", path.to_str().unwrap()], None);
    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("(SDF)"));
    // SDF carries coordinates; SMILES does not. The column proves the format
    // was actually honoured rather than merely reported.
    assert!(r.stdout.contains("\tyes"));
}

#[test]
fn test_format_can_be_forced_for_stdin_which_has_no_name() {
    let r = run(&["info", "--format", "sdf"], Some(SDF));
    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("(SDF)"));
    assert!(r.stdout.contains("ethanol"));
}

#[test]
fn test_an_unknown_extension_reads_as_smiles() {
    let path = fixture("nameless.dat", GOOD);
    let r = run(&["info", path.to_str().unwrap()], None);
    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("(SMILES)"));
}

#[test]
fn test_the_backend_choice_is_echoed_so_a_typo_is_visible() {
    let path = fixture("backend.smi", GOOD);
    let r = run(&["--backend", "cpu", "info", path.to_str().unwrap()], None);
    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("backend: cpu"));

    // An invalid value is clap's to reject, and it must not look like a
    // successful run with a default.
    let r = run(
        &["--backend", "quantum", "info", path.to_str().unwrap()],
        None,
    );
    assert_ne!(r.code, 0);
}

const SDF: &str = "\
ethanol
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 C   0  0
    2.2500    1.3000    0.0000 O   0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
";
