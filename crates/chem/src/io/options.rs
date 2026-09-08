//! Per-format read and write options.
//!
//! `ReadFn`/`WriteFn` in [`crate::io::format`] thread these through to every
//! registered reader and writer. Most formats have nothing to configure yet
//! — the types exist so a format's signature only has to widen once, per
//! the plan already recorded on those type aliases before any format needed
//! a real option (#212).

/// No options exist yet. A real field lands here the first time a reader
/// needs one to configure, rather than being invented ahead of that need.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReadOptions;

/// One field per format that has a write option today.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct WriteOptions {
    pub sdf: SdfWriteOptions,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfWriteOptions {
    pub version: MolfileVersion,
}

/// Which molfile dialect to write. Only `V2000` is implemented — this crate
/// has no V3000 writer yet, a different block structure entirely, not just
/// a flag. The variant exists because the choice itself is the reason this
/// type exists (#197 deferred V3000 until this option could be expressed);
/// implementing it is a separate story.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum MolfileVersion {
    #[default]
    V2000,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_write_options_pick_v2000() {
        assert_eq!(WriteOptions::default().sdf.version, MolfileVersion::V2000);
    }
}
