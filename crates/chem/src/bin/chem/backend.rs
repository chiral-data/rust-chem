//! Choosing between GPU and CPU.
//!
//! A global flag rather than a per-subcommand one: a batch script pins the
//! backend once for a whole pipeline, and having `chem fp --gpu | chem search
//! --gpu` be the spelling would invite one half of it to be forgotten.
//!
//! Probing costs a device creation, so only the commands that actually compute
//! pay it — [`Backend::open`] is where that happens, and the commands that just
//! read or describe never call it.

use anyhow::{Result, bail};
use chem::search::FingerprintSearch;
use clap::ValueEnum;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, Default)]
pub enum Backend {
    /// GPU if one is usable, CPU otherwise. Reports which on stderr.
    #[default]
    Auto,
    /// GPU, and fail rather than fall back — so a batch job cannot silently
    /// take a hundredfold slowdown.
    Gpu,
    /// CPU, skipping device probing entirely.
    Cpu,
}

impl Backend {
    pub fn label(&self) -> &'static str {
        match self {
            Backend::Auto => "auto",
            Backend::Gpu => "gpu",
            Backend::Cpu => "cpu",
        }
    }

    /// Builds a search engine on the requested backend, reporting on stderr
    /// which one it got.
    ///
    /// `Gpu` fails rather than falling back. That is the whole point of asking
    /// for it explicitly: a batch job pinned to the GPU that quietly ran on the
    /// CPU would take a large slowdown and report success, and the operator
    /// would find out from the wall clock rather than from the tool.
    pub fn open(&self) -> Result<FingerprintSearch> {
        match self {
            Backend::Cpu => {
                eprintln!("backend: cpu (by request)");
                Ok(FingerprintSearch::new_cpu_only())
            }
            Backend::Gpu => {
                let mut search = FingerprintSearch::new_cpu_only();
                if let Err(e) = search.retry_gpu_init() {
                    bail!("--backend gpu was requested but no GPU is usable: {e}");
                }
                if !search.is_using_gpu() {
                    bail!("--backend gpu was requested but the GPU did not engage");
                }
                eprintln!("backend: gpu");
                Ok(search)
            }
            Backend::Auto => {
                let search = FingerprintSearch::new();
                if search.is_using_gpu() {
                    eprintln!("backend: gpu");
                } else {
                    // Say *why*, so "it was slow" and "there is no GPU here"
                    // are distinguishable without guessing.
                    match search.gpu_init_error() {
                        Some(e) => eprintln!("backend: cpu (no usable GPU: {e})"),
                        None => eprintln!("backend: cpu"),
                    }
                }
                Ok(search)
            }
        }
    }
}
