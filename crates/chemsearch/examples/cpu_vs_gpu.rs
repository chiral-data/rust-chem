//! Where the GPU starts paying for itself.
//!
//! ```text
//! cargo run --release --example cpu_vs_gpu -p chemsearch -- [path.smi]
//! ```
//!
//! A single speedup figure is the least useful thing to report here, because
//! the answer is a function of batch size. The GPU has fixed costs per call —
//! uploading molecules, dispatching, reading results back — that the CPU does
//! not, so below some batch size the CPU wins and above it the GPU does. The
//! crossover is the number worth knowing: it is what decides whether
//! `--backend gpu` is worth passing.
//!
//! Run in `--release`. A debug CPU build is several times slower, which moves
//! the crossover to a batch size that does not exist.
//!
//! Molecules come from a SMILES file — `test.smi` at the repo root by default —
//! and are repeated to reach each batch size. Real molecules rather than
//! generated rings: fingerprinting cost depends on branching and ring
//! structure, so a corpus of identical synthetic rings measures a workload
//! nobody has.

use chemsearch::FingerprintSearch;
use std::time::{Duration, Instant};

const RADIUS: u32 = 2;
const SIZE: u32 = 2048;
/// Stops at 50,000 deliberately. Above roughly 200,000 molecules of this size
/// the GPU path panics inside wgpu — the scratch buffer outgrows
/// `max_storage_buffer_binding_size`, which the existing dispatch caps do not
/// bound (#158). The ceiling is a function of total atoms rather than molecule
/// count, so a corpus of larger molecules reaches it sooner and there is no
/// molecule count that is safe in general.
const BATCHES: [usize; 8] = [1, 10, 100, 500, 1_000, 5_000, 20_000, 50_000];

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).unwrap_or_else(|| "test.smi".into());
    let content =
        std::fs::read_to_string(&path).map_err(|e| anyhow::anyhow!("reading {path}: {e}"))?;
    let outcome = chemio::reader::read_smiles(&content);
    if outcome.is_empty() {
        anyhow::bail!("no molecules in {path}");
    }
    let corpus: Vec<_> = outcome.records.iter().map(|r| r.molecule.clone()).collect();
    println!("corpus: {} molecules from {path}", corpus.len());
    println!("params: radius {RADIUS}, {SIZE} bits\n");

    let cpu = FingerprintSearch::new_cpu_only();

    // Timed separately, because it is paid once per process and would otherwise
    // be smeared across whichever batch happened to run first. It is also the
    // reason a one-molecule GPU run looks terrible.
    let started = Instant::now();
    let mut gpu = FingerprintSearch::new_cpu_only();
    let gpu_available = gpu.retry_gpu_init().is_ok() && gpu.is_using_gpu();
    let init = started.elapsed();

    if !gpu_available {
        println!("no usable GPU on this machine — reporting CPU timings only");
        if let Some(e) = gpu.gpu_init_error() {
            println!("  ({e})\n");
        }
    } else {
        println!("GPU init: {:.0} ms (once per process)\n", ms(init));
    }

    println!(
        "{:>9}  {:>11}  {:>11}  {:>9}  agree",
        "molecules", "cpu (ms)", "gpu (ms)", "speedup"
    );

    for n in BATCHES {
        let batch: Vec<_> = corpus.iter().cycle().take(n).cloned().collect();

        let cpu_time = time(|| {
            pollster::block_on(cpu.generate_fingerprints_batch_async(&batch, RADIUS, SIZE))
        })?;

        if !gpu_available {
            println!(
                "{n:>9}  {:>11.1}  {:>11}  {:>9}  -",
                ms(cpu_time.1),
                "-",
                "-"
            );
            continue;
        }

        // One warm-up pass: the first dispatch pays for pipeline setup that
        // every later one reuses, and charging that to the smallest batch would
        // misreport it as a per-molecule cost.
        let _ = pollster::block_on(gpu.generate_fingerprints_batch_async(&batch, RADIUS, SIZE))?;
        let gpu_time = time(|| {
            pollster::block_on(gpu.generate_fingerprints_batch_async(&batch, RADIUS, SIZE))
        })?;

        // A benchmark that does not check agreement can report a fine speedup
        // for wrong answers. This is the same property the parity tests assert,
        // re-checked here because it is cheap and the alternative is a
        // meaningless number.
        let agree = cpu_time.0 == gpu_time.0;

        println!(
            "{n:>9}  {:>11.1}  {:>11.1}  {:>8.1}x  {}",
            ms(cpu_time.1),
            ms(gpu_time.1),
            ms(cpu_time.1) / ms(gpu_time.1),
            if agree { "yes" } else { "NO — MISMATCH" }
        );

        if !agree {
            anyhow::bail!("CPU and GPU disagreed at {n} molecules");
        }
    }

    if gpu_available {
        // The column above measures throughput once a device exists. A CLI
        // invocation does not get that for free: it pays init every run, and
        // that is what decides whether `--backend gpu` helps in a script.
        println!(
            "\nThe speedup column excludes the {:.0} ms of init, which a single\n\
             `chem fp` run pays every time. Including it, the GPU only wins once\n\
             the CPU time it saves exceeds {:.0} ms — far beyond the sizes above.\n\
             The GPU is for a long-lived process doing many operations, which is\n\
             what the workbench is and what a one-shot CLI call is not.",
            ms(init),
            ms(init)
        );
    }
    Ok(())
}

fn time<T, E>(mut f: impl FnMut() -> Result<T, E>) -> Result<(T, Duration), E> {
    let started = Instant::now();
    let value = f()?;
    Ok((value, started.elapsed()))
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}
