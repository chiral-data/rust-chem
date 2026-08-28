//! Where the GPU starts paying for itself on similarity search.
//!
//! ```text
//! cargo run --release --example search_cpu_vs_gpu -p chem --features gpu -- [path.smi]
//! ```
//!
//! Fingerprint *generation* is setup here, not the measurement: it runs once on
//! the CPU per batch and is not timed. What is timed is ranking one query
//! against N targets, which is the operation a workbench repeats and the one a
//! larger dataset makes more expensive.
//!
//! Two GPU numbers, because they answer different questions:
//!
//! - **warm** — targets already uploaded. This is the steady state: the app
//!   uploads a dataset once and then searches it as the user types.
//! - **cold** — upload included. This is a single `chem search` invocation,
//!   which gets no second query to amortise the upload over.
//!
//! A single figure hides the gap between those, and the gap is the whole
//! decision about when the GPU is worth reaching for.
//!
//! Run in `--release`. A debug CPU build is several times slower, which moves
//! the crossover to a batch size that does not exist.
//!
//! Molecules come from a SMILES file — `test.smi` at the repo root by default —
//! repeated to reach each batch size. Real molecules rather than generated
//! rings: fingerprint density depends on branching and ring structure, and
//! Tanimoto cost depends on density.

use chem::search::{FingerprintSearch, SearchResult};
use std::time::{Duration, Instant};

const RADIUS: u32 = 2;
const SIZE: u32 = 2048;
const TOP_K: usize = 10;

/// Stops at 200,000. Fingerprint *generation* panics inside wgpu above roughly
/// 250,000 molecules of this size (#158), and while search itself uploads a
/// smaller buffer, staying under the known-bad range keeps this example about
/// search rather than about that bug.
const BATCHES: [usize; 8] = [100, 500, 1_000, 5_000, 20_000, 50_000, 100_000, 200_000];

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).unwrap_or_else(|| "test.smi".into());
    let content =
        std::fs::read_to_string(&path).map_err(|e| anyhow::anyhow!("reading {path}: {e}"))?;
    let outcome = chem::io::reader::read_smiles(&content);
    if outcome.is_empty() {
        anyhow::bail!("no molecules in {path}");
    }
    let corpus: Vec<_> = outcome.records.iter().map(|r| r.molecule.clone()).collect();
    println!("corpus: {} molecules from {path}", corpus.len());
    println!("params: radius {RADIUS}, {SIZE} bits, top {TOP_K}");
    println!("fingerprint generation is setup here, and not timed\n");

    let cpu = FingerprintSearch::new_cpu_only();

    // Timed separately because it is paid once per process, and smearing it
    // across whichever batch ran first would misreport it as a per-query cost.
    let started = Instant::now();
    let mut gpu = FingerprintSearch::new_cpu_only();
    let gpu_available = gpu.retry_gpu_init().is_ok() && gpu.is_using_gpu();
    let init = started.elapsed();

    if gpu_available {
        println!("GPU init: {:.0} ms (once per process)\n", ms(init));
    } else {
        println!("no usable GPU on this machine — reporting CPU timings only");
        if let Some(e) = gpu.gpu_init_error() {
            println!("  ({e})");
        }
        println!();
    }

    println!(
        "{:>9}  {:>10}  {:>10}  {:>10}  {:>8}  agree",
        "targets", "cpu (ms)", "gpu warm", "gpu cold", "warm x"
    );

    for n in BATCHES {
        let batch: Vec<_> = corpus.iter().cycle().take(n).cloned().collect();

        // Setup, deliberately untimed.
        let targets =
            pollster::block_on(cpu.generate_fingerprints_batch_async(&batch, RADIUS, SIZE))?;
        let query = pollster::block_on(cpu.generate_fingerprint_async(&corpus[0], RADIUS, SIZE))?;

        let mut cpu_engine = FingerprintSearch::new_cpu_only();
        let (cpu_results, cpu_time) =
            time(|| pollster::block_on(cpu_engine.search_async(&query, &targets, TOP_K)))?;

        if !gpu_available {
            println!(
                "{n:>9}  {:>10.2}  {:>10}  {:>10}  {:>8}  -",
                ms(cpu_time),
                "-",
                "-",
                "-"
            );
            continue;
        }

        // Cold: what a single `chem search` pays, upload included. Dropping the
        // cache first, since a previous iteration may have left one.
        gpu.invalidate_target_dataset();
        let (_, cold) = time(|| -> anyhow::Result<Vec<SearchResult>> {
            gpu.set_target_dataset(&targets)?;
            pollster::block_on(gpu.search_async(&query, &targets, TOP_K))
        })?;

        // Warm: targets already resident, which is the app's steady state. One
        // untimed pass first, so pipeline setup is not charged to the first
        // measured query.
        let _ = pollster::block_on(gpu.search_async(&query, &targets, TOP_K))?;
        let (gpu_results, warm) =
            time(|| pollster::block_on(gpu.search_async(&query, &targets, TOP_K)))?;

        // A benchmark that skips this can report a fine speedup for wrong
        // answers. Compared with a tolerance because the GPU accumulates in f32
        // and the CPU in f64.
        let agree = same_ranking(&cpu_results, &gpu_results);

        println!(
            "{n:>9}  {:>10.2}  {:>10.2}  {:>10.2}  {:>7.1}x  {}",
            ms(cpu_time),
            ms(warm),
            ms(cold),
            ms(cpu_time) / ms(warm),
            if agree { "yes" } else { "NO — MISMATCH" }
        );

        if !agree {
            anyhow::bail!("CPU and GPU disagreed at {n} targets");
        }
    }

    if gpu_available {
        println!(
            "\n`warm x` is the speedup once the dataset is resident — the workbench's\n\
             case, where one upload serves every keystroke. `gpu cold` is a single\n\
             `chem search` run, which has no second query to amortise the upload\n\
             over, and additionally pays the {:.0} ms of init above.",
            ms(init)
        );
    }
    Ok(())
}

/// Same molecules in the same order, with similarities within f32 precision.
fn same_ranking(a: &[SearchResult], b: &[SearchResult]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| x.index == y.index && (x.similarity - y.similarity).abs() < 1e-6)
}

fn time<T, E>(mut f: impl FnMut() -> Result<T, E>) -> Result<(T, Duration), E> {
    let started = Instant::now();
    let value = f()?;
    Ok((value, started.elapsed()))
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}
