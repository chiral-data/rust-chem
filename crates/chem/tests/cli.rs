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

// ---- chem fp / chem search (#140) ------------------------------------------

/// Every test below pins `--backend cpu`. A GPU is not present on every machine
/// that runs these, and a test that quietly changes which code path it exercises
/// depending on the host is a test that proves less than it appears to.
/// CPU/GPU agreement is asserted once, separately, and skipped where there is
/// no device.
const CPU: &str = "--backend";

#[test]
fn test_fp_writes_a_self_describing_file() {
    let path = fixture("fp-basic.smi", GOOD);
    let r = run(&[CPU, "cpu", "fp", path.to_str().unwrap()], None);

    assert_eq!(r.code, 0);
    let lines: Vec<&str> = r.stdout.lines().collect();
    // The header is what lets `chem search` refuse a mismatched query rather
    // than silently comparing incomparable fingerprints.
    assert_eq!(lines[0], "# chem-fingerprints 1");
    assert!(lines.contains(&"# radius 2"));
    assert!(lines.contains(&"# size 2048"));
    assert_eq!(lines[3], "name\tfingerprint");
    assert_eq!(lines.len(), 7, "3 header + 1 column row + 3 molecules");
    assert!(r.stderr.contains("fingerprinted 3 molecules"));
}

#[test]
fn test_fp_honours_the_requested_parameters() {
    let path = fixture("fp-params.smi", GOOD);
    let r = run(
        &[
            CPU,
            "cpu",
            "fp",
            path.to_str().unwrap(),
            "--radius",
            "3",
            "--size",
            "512",
        ],
        None,
    );
    assert_eq!(r.code, 0);
    assert!(r.stdout.contains("# radius 3"));
    assert!(r.stdout.contains("# size 512"));
    // 512 bits is 128 hex characters.
    let row = r.stdout.lines().nth(4).expect("a molecule row");
    let hex = row.split('\t').nth(1).expect("a fingerprint column");
    assert_eq!(hex.len(), 128);
}

#[test]
fn test_fp_and_search_compose_through_a_pipe() {
    // The reason these are two commands rather than one. If either put its
    // progress on stdout, this would feed garbage into the second process.
    let path = fixture("pipe.smi", GOOD);
    let fp = run(&[CPU, "cpu", "fp", path.to_str().unwrap()], None);
    assert_eq!(fp.code, 0);

    let search = run(
        &[CPU, "cpu", "search", "--query", "CCO", "--top", "1"],
        Some(&fp.stdout),
    );
    assert_eq!(search.code, 0);
    let rows: Vec<&str> = search.stdout.lines().collect();
    assert_eq!(rows[0], "rank\tname\tsimilarity");
    // Ethanol against itself.
    assert!(rows[1].starts_with("1\tethanol\t1.000000"), "{:?}", rows[1]);
}

#[test]
fn test_search_fingerprints_the_query_with_the_files_parameters() {
    // Not with its own defaults. A query at radius 2 compared against targets
    // at radius 3 produces plausible-looking similarities rather than an error,
    // which is the failure the header exists to prevent.
    let path = fixture("search-params.smi", GOOD);
    let fp = run(
        &[
            CPU,
            "cpu",
            "fp",
            path.to_str().unwrap(),
            "--radius",
            "3",
            "--size",
            "512",
        ],
        None,
    );
    let search = run(
        &[CPU, "cpu", "search", "--query", "CCO", "--top", "1"],
        Some(&fp.stdout),
    );
    assert_eq!(search.code, 0);
    assert!(search.stderr.contains("radius 3"));
    assert!(search.stderr.contains("512 bits"));
    // An identical molecule must still score 1.0 — which only holds if the
    // query used the same parameters.
    assert!(search.stdout.contains("1.000000"), "{:?}", search.stdout);
}

#[test]
fn test_search_ranks_by_descending_similarity() {
    let path = fixture("rank.smi", "CCO ethanol\nCCCO propanol\nc1ccccc1 benzene\n");
    let fp = run(&[CPU, "cpu", "fp", path.to_str().unwrap()], None);
    let search = run(
        &[CPU, "cpu", "search", "--query", "CCO", "--top", "0"],
        Some(&fp.stdout),
    );

    assert_eq!(search.code, 0);
    let scores: Vec<f64> = search
        .stdout
        .lines()
        .skip(1)
        .map(|l| l.split('\t').nth(2).unwrap().parse().unwrap())
        .collect();
    assert_eq!(scores.len(), 3, "--top 0 returns everything");
    assert!(
        scores.windows(2).all(|w| w[0] >= w[1]),
        "not descending: {scores:?}"
    );
    assert_eq!(scores[0], 1.0, "the query itself ranks first");
}

#[test]
fn test_ties_come_out_in_dataset_order() {
    // Two identical molecules under different names score identically. The sort
    // is stable, so they must come out in the order the file listed them —
    // otherwise a diff of two runs is noisy for no reason. This is currently
    // guaranteed only by `sort_by` being stable, which is exactly why it is
    // worth pinning: switching to `sort_unstable_by` would break it silently.
    let path = fixture("ties.smi", "CCO first\nCCO second\nCCO third\n");
    let fp = run(&[CPU, "cpu", "fp", path.to_str().unwrap()], None);
    let search = run(
        &[CPU, "cpu", "search", "--query", "CCO", "--top", "0"],
        Some(&fp.stdout),
    );

    assert_eq!(search.code, 0);
    let names: Vec<&str> = search
        .stdout
        .lines()
        .skip(1)
        .map(|l| l.split('\t').nth(1).unwrap())
        .collect();
    assert_eq!(names, ["first", "second", "third"]);
}

#[test]
fn test_search_refuses_a_file_that_is_not_a_fingerprint_file() {
    // `chem search mols.smi` instead of `chem search fps.tsv` is the ordinary
    // mistake, and the error should name it rather than describe a bad row.
    let path = fixture("not-fps.smi", GOOD);
    let r = run(
        &[
            CPU,
            "cpu",
            "search",
            path.to_str().unwrap(),
            "--query",
            "CCO",
        ],
        None,
    );
    assert_eq!(r.code, 1);
    assert!(
        r.stderr.contains("not a chem fingerprint file"),
        "{:?}",
        r.stderr
    );
}

#[test]
fn test_fp_on_an_unusable_file_exits_before_touching_a_backend() {
    let path = fixture("fp-allbad.smi", ALL_BAD);
    let r = run(&[CPU, "cpu", "fp", path.to_str().unwrap()], None);
    assert_eq!(r.code, 2);
    assert!(r.stdout.is_empty());
}

#[test]
fn test_fp_reports_timing_on_stderr_not_stdout() {
    let path = fixture("fp-timing.smi", GOOD);
    let r = run(&[CPU, "cpu", "fp", path.to_str().unwrap()], None);
    assert!(r.stderr.contains(" ms ("));
    assert!(!r.stdout.contains(" ms "));
}

#[test]
fn test_a_gpu_request_fails_loudly_when_there_is_no_gpu() {
    // Either outcome is correct depending on the host; what must never happen
    // is exit 0 having quietly used the CPU. A batch job pinned to the GPU that
    // silently ran on the CPU reports success and takes the slowdown.
    let path = fixture("gpu-req.smi", GOOD);
    let r = run(&[CPU, "gpu", "fp", path.to_str().unwrap()], None);
    if r.code == 0 {
        assert!(r.stderr.contains("backend: gpu"), "{:?}", r.stderr);
    } else {
        assert!(r.stderr.contains("no GPU is usable") || r.stderr.contains("did not engage"));
    }
}

// ---- chem aromatic / chem coords (#141) ------------------------------------

/// Kekulé SMILES — explicit alternating bonds, uppercase atoms. Lowercase
/// aromatic input is already flagged by the parser, so a test written that way
/// would pass whether or not `chem aromatic` did anything.
const KEKULE: &str = "C1=CC=NC=C1 pyridine\nC1=CC=CC=C1 benzene\nCCO ethanol\n";

#[test]
fn test_aromatic_perceives_rings_and_writes_them_lowercase() {
    let path = fixture("arom.smi", KEKULE);
    let r = run(&["aromatic", path.to_str().unwrap()], None);

    assert_eq!(r.code, 0);
    let lines: Vec<&str> = r.stdout.lines().collect();
    assert_eq!(lines[0], "c1ccncc1 pyridine", "pyridine must be perceived");
    assert_eq!(lines[1], "c1ccccc1 benzene");
    // Not everything is aromatic, and the command must not claim otherwise.
    assert_eq!(lines[2], "CCO ethanol");
}

#[test]
fn test_aromatic_reports_how_many_molecules_it_changed() {
    // "it did nothing" and "nothing needed doing" produce identical output
    // files, and only one of them is a problem.
    let path = fixture("arom-count.smi", KEKULE);
    let r = run(&["aromatic", path.to_str().unwrap()], None);
    assert!(
        r.stderr.contains("2 of 3 molecules changed"),
        "{:?}",
        r.stderr
    );
}

#[test]
fn test_aromatic_output_round_trips_back_to_aromatic() {
    // The point of defaulting to SMILES here: the perception survives, because
    // lowercase atoms carry it. If the writer emitted uppercase, this command
    // would compute a result and then throw it away.
    let path = fixture("arom-rt.smi", KEKULE);
    let first = run(&["aromatic", path.to_str().unwrap()], None);
    let second = run(&["aromatic"], Some(&first.stdout));

    assert_eq!(second.code, 0);
    assert_eq!(
        first.stdout, second.stdout,
        "already-aromatic input is stable"
    );
    assert!(second.stderr.contains("0 of 3 molecules changed"));
}

#[test]
fn test_coords_defaults_to_sdf_because_smiles_cannot_hold_coordinates() {
    let path = fixture("coords.smi", KEKULE);
    let r = run(&["coords", path.to_str().unwrap()], None);

    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("writing SDF"));
    assert!(r.stdout.contains("V2000"));
    assert!(r.stdout.trim_end().ends_with("$$$$"));
}

#[test]
fn test_coords_output_actually_carries_the_coordinates() {
    // The whole reason the default is SDF. Reading it back and finding `yes` in
    // the coords column is what proves the result was not discarded — these
    // molecules came from SMILES, which has nowhere to put geometry.
    let path = fixture("coords-carry.smi", KEKULE);
    let coords = run(&["coords", path.to_str().unwrap()], None);
    let info = run(&["info", "--format", "sdf"], Some(&coords.stdout));

    assert_eq!(info.code, 0);
    let rows: Vec<&str> = info.stdout.lines().skip(1).collect();
    assert_eq!(rows.len(), 3);
    for row in rows {
        // The 2D column specifically: `chem coords` computes a depiction, and
        // a layout written to SDF comes back a layout rather than a conformer.
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[3], "yes", "no layout in {row:?}");
        assert_eq!(fields[4], "no", "a layout must not become a conformer");
    }
}

#[test]
fn test_coords_on_a_3d_file_writes_the_depiction_it_computed() {
    // A conformer outranks a layout in the SDF writer, so without an explicit
    // drop `chem coords` would report "laid out 1" and then emit the input
    // untouched — doing the work and discarding it. The depiction has to win
    // here, because producing one is the entire command.
    let sdf = "\
Acetone
APtclcactv06051922463D 0   0.00000     0.00000

  4  3  0  0  0  0  0  0  0  0999 V2000
    1.3051    0.6772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.0763    0.8900 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3051    0.6772   -0.8900 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.2839    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
M  END
$$$$
";
    let path = fixture("coords-3d.sdf", sdf);
    let r = run(&["coords", path.to_str().unwrap()], None);

    assert_eq!(r.code, 0);
    assert!(
        r.stderr
            .contains("discarded the 3D conformer of 1 molecules"),
        "dropping geometry has to be announced: {:?}",
        r.stderr
    );

    // Every z is zero now, and the header says 2D rather than 3D.
    let z: Vec<&str> = r
        .stdout
        .lines()
        .skip(4)
        .take(4)
        .map(|line| line.split_whitespace().nth(2).unwrap())
        .collect();
    assert_eq!(
        z,
        vec!["0.0000", "0.0000", "0.0000", "0.0000"],
        "{:?}",
        r.stdout
    );
    assert!(r.stdout.lines().nth(1).unwrap().contains("2D"));

    // And reading it back finds a layout, not a conformer.
    let info = run(&["info", "--format", "sdf"], Some(&r.stdout));
    let row = info.stdout.lines().nth(1).unwrap();
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(fields[3], "yes", "no layout in {row:?}");
    assert_eq!(
        fields[4], "no",
        "conformer should have been dropped: {row:?}"
    );
}

#[test]
fn test_coords_warns_rather_than_silently_discarding_the_result() {
    let path = fixture("coords-warn.smi", KEKULE);

    // Explicitly asking for SMILES is allowed — connectivity only is a real
    // want — but it must say what it is costing. The wording comes from the
    // format registry now rather than a hand-written string, so this asserts
    // on the attribute that was lost rather than on a sentence.
    let r = run(
        &["coords", path.to_str().unwrap(), "--out-format", "smiles"],
        None,
    );
    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("cannot carry"), "{:?}", r.stderr);
    assert!(r.stderr.contains("coords_2d"), "{:?}", r.stderr);

    // And an output name that is not .sdf is the same mistake made implicitly.
    // Named `-dest` rather than reusing the fixture stem: the first draft of
    // this test picked the same path as its own input, and the clobber guard
    // refused it — correctly, but that is a different assertion.
    let dest = std::env::temp_dir().join("chem-cli-test-coords-warn-dest.smi");
    let r = run(
        &[
            "coords",
            path.to_str().unwrap(),
            "-o",
            dest.to_str().unwrap(),
        ],
        None,
    );
    assert!(r.stderr.contains("cannot carry"), "{:?}", r.stderr);
    assert!(r.stderr.contains("coords_2d"), "{:?}", r.stderr);
}

/// Two nitro compounds, whose charges an SDF write drops.
///
/// The registry says `Format::SDF` does not carry `formal_charge`, because
/// `write_sdf` documents that it does not write the field. Real files in
/// `test.smi` hit this, so it is not a contrived case.
/// Two molecules that lose the *same* attribute when written.
///
/// **Coordinates on the way to SMILES**, which is a structural limit rather
/// than a gap: a SMILES string has nowhere to put a position and never will.
///
/// The earlier choices both stopped being losses as the writers improved —
/// formal charge until #197 taught V2000 `M  CHG`, then double-bond geometry
/// until #198 taught the layout to draw it — and each time these tests failed
/// by reporting *nothing*, which was the drop report being honest rather than
/// broken. Picking a loss that cannot be fixed stops that recurring.
const LAID_OUT: &str = "F/C=C/F trans-difluoroethene\nF/C=C\\F cis-difluoroethene\n";

#[test]
fn test_a_write_that_loses_nothing_says_nothing() {
    // The regression a general mechanism most easily introduces: noise on the
    // path that was previously quiet. SMILES in, SMILES out drops nothing.
    let path = fixture("drops-none.smi", GOOD);
    let r = run(&["aromatic", path.to_str().unwrap()], None);

    assert_eq!(r.code, 0);
    assert!(
        !r.stderr.contains("cannot carry"),
        "a lossless write reported a loss: {:?}",
        r.stderr
    );
}

#[test]
fn test_the_summary_names_each_lost_attribute_once() {
    // Two molecules lose the same attribute, and the default report is one
    // line with a count rather than one line per molecule — the whole reason
    // the per-molecule form is behind a flag.
    let path = fixture("drops-charge.smi", LAID_OUT);
    let dest = std::env::temp_dir().join("chem-cli-test-drops-charge-out.smi");
    let _ = std::fs::remove_file(&dest);

    // `coords` computes a layout, and SMILES has nowhere to keep it.
    let r = run(
        &[
            "coords",
            path.to_str().unwrap(),
            "-o",
            dest.to_str().unwrap(),
        ],
        None,
    );

    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("coords_2d (2)"), "{:?}", r.stderr);
    // The molecules are not named without the flag.
    assert!(!r.stderr.contains("difluoroethene"), "{:?}", r.stderr);
    assert!(r.stderr.contains("--explain-drops"), "{:?}", r.stderr);
}

#[test]
fn test_explain_drops_names_the_molecules() {
    let path = fixture("drops-explain.smi", LAID_OUT);
    let dest = std::env::temp_dir().join("chem-cli-test-drops-explain-out.smi");
    let _ = std::fs::remove_file(&dest);

    let r = run(
        &[
            "coords",
            path.to_str().unwrap(),
            "-o",
            dest.to_str().unwrap(),
            "--explain-drops",
        ],
        None,
    );

    assert_eq!(r.code, 0);
    assert!(
        r.stderr.contains("trans-difluoroethene: coords_2d"),
        "{:?}",
        r.stderr
    );
    assert!(
        r.stderr.contains("cis-difluoroethene: coords_2d"),
        "{:?}",
        r.stderr
    );
    // The pointer to the flag is pointless once the flag is on.
    assert!(!r.stderr.contains("--explain-drops"), "{:?}", r.stderr);
}

#[test]
fn test_an_sdf_name_is_honoured_without_the_flag() {
    let path = fixture("coords-ext.smi", KEKULE);
    let dest = std::env::temp_dir().join("chem-cli-test-coords-ext.sdf");
    let _ = std::fs::remove_file(&dest);

    let r = run(
        &[
            "coords",
            path.to_str().unwrap(),
            "-o",
            dest.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(r.code, 0);
    assert!(!r.stderr.contains("warning:"), "{:?}", r.stderr);
    assert!(
        std::fs::read_to_string(&dest)
            .expect("written")
            .contains("V2000")
    );
}

#[test]
fn test_writing_over_the_input_is_refused_unless_forced() {
    // It would happen to work today, because the input is read fully before the
    // write begins — which is worse than failing, since it would keep working
    // until the day reading became streamed.
    let path = fixture("clobber.smi", KEKULE);
    let arg = path.to_str().unwrap();

    let r = run(&["coords", arg, "-o", arg], None);
    assert_eq!(r.code, 1);
    assert!(r.stderr.contains("is the input file"), "{:?}", r.stderr);
    // And the file is untouched.
    assert_eq!(std::fs::read_to_string(&path).expect("intact"), KEKULE);

    let r = run(&["aromatic", arg, "-o", arg, "--force"], None);
    assert_eq!(r.code, 0);
    assert!(
        std::fs::read_to_string(&path)
            .expect("rewritten")
            .contains("c1ccncc1")
    );
}

#[test]
fn test_relayout_distinguishes_kept_from_recomputed() {
    let smi = fixture("relayout.smi", KEKULE);
    let sdf = run(&["coords", smi.to_str().unwrap()], None).stdout;

    let kept = run(&["coords", "--format", "sdf"], Some(&sdf));
    assert!(
        kept.stderr.contains("laid out 0, kept 3"),
        "{:?}",
        kept.stderr
    );

    let redone = run(&["coords", "--format", "sdf", "--relayout"], Some(&sdf));
    assert!(
        redone.stderr.contains("laid out 3, kept 0"),
        "{:?}",
        redone.stderr
    );
}

#[test]
fn test_a_cpu_only_command_notes_the_backend_flag_instead_of_failing() {
    // A script pinning --backend globally should not have to special-case which
    // subcommands can use it.
    let path = fixture("cpu-note.smi", KEKULE);
    let r = run(
        &["--backend", "gpu", "aromatic", path.to_str().unwrap()],
        None,
    );
    assert_eq!(r.code, 0, "must not fail: {:?}", r.stderr);
    assert!(r.stderr.contains("no GPU path"), "{:?}", r.stderr);
}

#[test]
fn test_aromaticity_changes_the_fingerprint() {
    // Not a formatting detail: perceiving aromaticity changes the bonds, so it
    // changes the Morgan environments. This is why `chem aromatic` is worth
    // having in a pipeline rather than being cosmetic.
    let path = fixture("arom-fp.smi", KEKULE);
    let plain = run(&[CPU, "cpu", "fp", path.to_str().unwrap()], None);
    let perceived = run(&["aromatic", path.to_str().unwrap()], None);
    let after = run(&[CPU, "cpu", "fp"], Some(&perceived.stdout));

    assert_eq!(plain.code, 0);
    assert_eq!(after.code, 0);
    assert_ne!(
        plain.stdout, after.stdout,
        "aromatic perception must reach the fingerprint"
    );
}

// ---- chem draw (#142) ------------------------------------------------------

/// Parses far enough to tell a real SVG from a string that merely starts with
/// `<svg`. Counts the drawable elements, since a document with none is a blank
/// page that still passes a "does it contain `<svg`" check.
fn svg_marks(text: &str) -> usize {
    text.matches("<line").count() + text.matches("<text").count() + text.matches("<circle").count()
}

#[test]
fn test_a_single_molecule_goes_to_stdout_as_one_document() {
    let r = run(&["draw"], Some("c1ccccc1O phenol\n"));
    assert_eq!(r.code, 0);
    assert!(r.stdout.starts_with("<svg"));
    assert!(r.stdout.trim_end().ends_with("</svg>"));
    assert_eq!(r.stdout.matches("<svg").count(), 1);
    assert!(
        svg_marks(&r.stdout) > 6,
        "phenol should draw more than 6 marks"
    );
}

#[test]
fn test_many_molecules_without_outdir_is_refused_with_the_fix_named() {
    // Concatenated SVG documents are not a valid SVG, so this cannot be done
    // silently. The error has to name the flag, or the user is left guessing.
    let path = fixture("draw-many.smi", GOOD);
    let r = run(&["draw", path.to_str().unwrap()], None);
    assert_eq!(r.code, 1);
    assert!(r.stderr.contains("--outdir"), "{:?}", r.stderr);
    assert!(r.stdout.is_empty(), "no partial document may escape");
}

#[test]
fn test_outdir_writes_one_file_per_molecule() {
    let path = fixture("draw-dir.smi", GOOD);
    let dir = std::env::temp_dir().join("chem-cli-test-drawdir");
    let _ = std::fs::remove_dir_all(&dir);

    let r = run(
        &[
            "draw",
            path.to_str().unwrap(),
            "--outdir",
            dir.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(r.code, 0);

    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .expect("directory created")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert_eq!(names, ["acetic.svg", "benzene.svg", "ethanol.svg"]);

    for name in &names {
        let text = std::fs::read_to_string(dir.join(name)).expect(name);
        assert!(text.starts_with("<svg"), "{name}");
        assert!(svg_marks(&text) > 0, "{name} drew nothing");
    }
}

#[test]
fn test_duplicate_molecule_names_do_not_overwrite_each_other() {
    // Molecule names come from a file's own records and are not unique. Without
    // suffixing, this reports three files and leaves one — silently short.
    let path = fixture("draw-dup.smi", "CCO Phenol\nc1ccccc1O Phenol\nCC Phenol!\n");
    let dir = std::env::temp_dir().join("chem-cli-test-drawdup");
    let _ = std::fs::remove_dir_all(&dir);

    let r = run(
        &[
            "draw",
            path.to_str().unwrap(),
            "--outdir",
            dir.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(r.code, 0);
    assert_eq!(std::fs::read_dir(&dir).expect("dir").count(), 3);
    // And it says so, since otherwise the caller believes the filenames match
    // the molecule names.
    assert!(r.stderr.contains("duplicate names"), "{:?}", r.stderr);
}

#[test]
fn test_coordinates_are_generated_and_the_command_says_so() {
    // SMILES carries no coordinates, so drawing one means doing work that was
    // not asked for. Doing it silently would make `chem coords` look redundant.
    let r = run(&["draw"], Some("CCO ethanol\n"));
    assert_eq!(r.code, 0);
    assert!(r.stderr.contains("generated coordinates"), "{:?}", r.stderr);
    assert!(
        r.stderr.contains("chem coords"),
        "should point at the explicit form"
    );
}

#[test]
fn test_an_sdf_with_coordinates_is_drawn_without_relaying_out() {
    let smi = fixture("draw-sdf-src.smi", "CCO ethanol\n");
    let sdf = run(&["coords", smi.to_str().unwrap()], None).stdout;

    let r = run(&["draw", "--format", "sdf"], Some(&sdf));
    assert_eq!(r.code, 0);
    assert!(
        !r.stderr.contains("generated coordinates"),
        "existing coordinates must be used as-is: {:?}",
        r.stderr
    );
}

#[test]
fn test_the_theme_flag_changes_the_palette_and_defaults_to_light() {
    // An SVG is bound for a document or a slide, so light is the honest
    // default — a CLI has no ambient theme to follow.
    let light = run(&["draw", "--theme", "light"], Some("c1ccccc1O x\n"));
    let dark = run(&["draw", "--theme", "dark"], Some("c1ccccc1O x\n"));
    let default = run(&["draw"], Some("c1ccccc1O x\n"));

    assert!(light.stdout.contains("#222222"), "light foreground");
    assert!(dark.stdout.contains("#ffffff"), "dark foreground");
    assert_eq!(default.stdout, light.stdout, "default must be light");
    // The heteroatom hues are shared between the palettes, so a test that only
    // compared those would pass whatever the flag did.
    assert!(light.stdout.contains("#e74c3c") && dark.stdout.contains("#e74c3c"));
}

#[test]
fn test_the_size_flags_reach_the_document() {
    let r = run(
        &["draw", "--width", "800", "--height", "600"],
        Some("CCO x\n"),
    );
    assert_eq!(r.code, 0);
    assert!(r.stdout.contains("width=\"800\""), "{:?}", &r.stdout[..80]);
    assert!(r.stdout.contains("height=\"600\""));
    assert!(r.stdout.contains("viewBox=\"0 0 800 600\""));
}

#[test]
fn test_output_flag_writes_the_single_structure_to_a_file() {
    let dest = std::env::temp_dir().join("chem-cli-test-draw-out.svg");
    let _ = std::fs::remove_file(&dest);

    let r = run(&["draw", "-o", dest.to_str().unwrap()], Some("CCO x\n"));
    assert_eq!(r.code, 0);
    assert!(r.stdout.is_empty());
    assert!(
        std::fs::read_to_string(&dest)
            .expect("written")
            .starts_with("<svg")
    );
}

#[test]
fn test_drawing_an_unusable_file_exits_before_writing_anything() {
    let path = fixture("draw-allbad.smi", ALL_BAD);
    let dir = std::env::temp_dir().join("chem-cli-test-drawbad");
    let _ = std::fs::remove_dir_all(&dir);

    let r = run(
        &[
            "draw",
            path.to_str().unwrap(),
            "--outdir",
            dir.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(r.code, 2);
    assert!(!dir.exists(), "no directory should be created for nothing");
}

#[test]
fn test_convert_round_trips_smiles_to_sdf_and_back() {
    let path = fixture("convert-in.smi", GOOD);
    let sdf = std::env::temp_dir().join("chem-cli-test-convert-out.sdf");
    let _ = std::fs::remove_file(&sdf);

    let r = run(
        &[
            "convert",
            path.to_str().unwrap(),
            "-o",
            sdf.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(
        r.stderr.contains("converted 3, skipped 0"),
        "{:?}",
        r.stderr
    );
    let written = std::fs::read_to_string(&sdf).expect("output file");
    assert!(written.contains("ethanol"));
    assert!(written.contains("$$$$"));

    // And back, to stdout as SMILES via --to.
    let r = run(&["convert", sdf.to_str().unwrap(), "--to", "smi"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(r.stdout.contains("ethanol"));
    assert!(r.stdout.contains("benzene"));
}

#[test]
fn test_convert_from_and_to_override_the_extension() {
    // Named .txt, so the extension alone would resolve to neither format.
    let path = fixture("convert-wrong-ext.txt", GOOD);
    let dest = std::env::temp_dir().join("chem-cli-test-convert-wrong-ext-out.txt");
    let _ = std::fs::remove_file(&dest);

    let r = run(
        &[
            "convert",
            path.to_str().unwrap(),
            "--from",
            "smi",
            "--to",
            "sdf",
            "-o",
            dest.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    let written = std::fs::read_to_string(&dest).expect("output file");
    assert!(written.contains("$$$$"));
}

#[test]
fn test_convert_literal_reads_the_string_directly() {
    let r = run(&["convert", "--literal", "c1ccccc1O", "--to", "sdf"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(r.stdout.contains("$$$$"));
    assert!(
        r.stderr.contains("converted 1, skipped 0"),
        "{:?}",
        r.stderr
    );
}

#[test]
fn test_convert_gen3d_is_rejected_with_a_clear_message() {
    let r = run(
        &["convert", "--literal", "CCO", "--to", "sdf", "--gen3d"],
        None,
    );
    assert_ne!(r.code, 0);
    assert!(
        r.stderr.contains("does not generate 3D coordinates"),
        "{:?}",
        r.stderr
    );
}

#[test]
fn test_convert_output_format_is_ambiguous_without_to_or_a_named_extension() {
    let r = run(&["convert", "--literal", "CCO"], None);
    assert_ne!(r.code, 0);
    assert!(r.stderr.contains("ambiguous"), "{:?}", r.stderr);
}

#[test]
fn test_convert_unrecognized_format_code_is_an_error() {
    let r = run(&["convert", "--literal", "CCO", "--to", "pdbqt"], None);
    assert_ne!(r.code, 0);
    assert!(
        r.stderr.contains("unrecognized format code"),
        "{:?}",
        r.stderr
    );
}

#[test]
fn test_convert_strict_exits_partial_on_a_skipped_record() {
    let path = fixture("convert-mixed.smi", MIXED);
    let r = run(
        &["--strict", "convert", path.to_str().unwrap(), "--to", "sdf"],
        None,
    );
    assert_eq!(r.code, 4);
}

#[test]
fn test_convert_refuses_to_clobber_its_own_input() {
    let path = fixture("convert-clobber.sdf", SDF);
    let r = run(
        &[
            "convert",
            path.to_str().unwrap(),
            "--to",
            "smi",
            "-o",
            path.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(r.code, 1);
    assert!(r.stderr.contains("pass --force"), "{:?}", r.stderr);
}

#[test]
fn test_format_bogus_is_an_error_not_a_clap_panic() {
    // #215: --format used to be a clap ValueEnum, so an unknown value was
    // rejected by clap itself. Now it's a registry lookup -- same outcome,
    // different mechanism, and the message should say so plainly.
    let r = run(&["info", "--format", "bogus"], Some(GOOD));
    assert_ne!(r.code, 0);
    assert!(
        r.stderr.contains("unrecognized format code"),
        "{:?}",
        r.stderr
    );
}

#[test]
fn test_convert_list_formats_lists_every_registered_code() {
    let r = run(&["convert", "-L", "formats"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(r.stdout.contains("smi"), "{:?}", r.stdout);
    assert!(r.stdout.contains("sdf"), "{:?}", r.stdout);
}

#[test]
fn test_convert_list_one_format_shows_its_detail() {
    let r = run(&["convert", "-L", "sdf"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(r.stdout.contains("sdf"), "{:?}", r.stdout);
    assert!(r.stdout.contains("codes:"), "{:?}", r.stdout);
    assert!(r.stdout.contains("extensions:"), "{:?}", r.stdout);
    assert!(r.stdout.contains("category:"), "{:?}", r.stdout);
    assert!(r.stdout.contains("carries:"), "{:?}", r.stdout);
}

#[test]
fn test_convert_list_unrecognized_code_is_an_error() {
    let r = run(&["convert", "-L", "bogus"], None);
    assert_ne!(r.code, 0);
    assert!(
        r.stderr.contains("unrecognized format code"),
        "{:?}",
        r.stderr
    );
}

#[test]
fn test_convert_describe_reports_no_options_exposed_yet() {
    let r = run(&["convert", "-H", "sdf"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(
        r.stdout.contains("no format-specific options are exposed"),
        "{:?}",
        r.stdout
    );

    let r = run(&["convert", "-H", "smi"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(
        r.stdout.contains("no format-specific options are exposed"),
        "{:?}",
        r.stdout
    );
}

#[test]
fn test_convert_list_and_describe_are_pure_queries_that_run_no_conversion() {
    // Neither flag should require --literal/input/--to -- they're pure
    // queries against the registry, resolved before any of that is checked.
    let r = run(&["convert", "-L", "formats"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    let r = run(&["convert", "-H", "sdf"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
}

#[test]
fn test_convert_writes_an_enhanced_stereo_group_as_cxsmiles() {
    // #221: --to cxsmiles needed no main.rs changes at all, since #215
    // already made every format-code flag a registry lookup.
    let r = run(
        &[
            "convert",
            "--literal",
            "N[C@@H](C)C(=O)O",
            "--from",
            "smi",
            "--to",
            "cxsmiles",
        ],
        None,
    );
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    // Plain SMILES has no stereo groups, so nothing to write -- confirms
    // the block is genuinely omitted rather than always present.
    assert!(!r.stdout.contains('|'), "{:?}", r.stdout);
}

#[test]
fn test_convert_round_trips_an_enhanced_stereo_group_through_cxsmiles() {
    let r = run(
        &[
            "convert",
            "--literal",
            "N[C@@H](C)C(=O)O |&1:1|",
            "--from",
            "cxsmiles",
            "--to",
            "cxsmiles",
        ],
        None,
    );
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(r.stdout.contains("|&1:1|"), "{:?}", r.stdout);
}

#[test]
fn test_convert_reports_dropping_a_stereo_group_when_writing_plain_smiles() {
    // Plain SMILES's carries mask does not include STEREO_GROUP, so the
    // drop report (#214's DropTracker) should fire rather than silently
    // discard the group.
    let r = run(
        &[
            "convert",
            "--literal",
            "N[C@@H](C)C(=O)O |&1:1|",
            "--from",
            "cxsmiles",
            "--to",
            "smi",
        ],
        None,
    );
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(r.stderr.contains("stereo_group"), "{:?}", r.stderr);
}

#[test]
fn test_convert_list_shows_cxsmiles_is_registered() {
    let r = run(&["convert", "-L", "cxsmiles"], None);
    assert_eq!(r.code, 0, "{:?}", r.stderr);
    assert!(r.stdout.contains("stereo_group"), "{:?}", r.stdout);
}
