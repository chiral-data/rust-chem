//! Choosing between GPU and CPU.
//!
//! A global flag rather than a per-subcommand one: a batch script pins the
//! backend once for a whole pipeline, and having `chem fp --gpu | chem search
//! --gpu` be the spelling would invite one half of it to be forgotten.
//!
//! Resolution is deliberately not done here. Probing for a device costs a
//! device creation, and the commands that do not compute — reading, describing,
//! writing — should not pay it. Each command that needs a backend resolves
//! [`Backend`] itself, which is also where a fallback message belongs, since
//! only that command knows what it fell back *to*.

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
}
