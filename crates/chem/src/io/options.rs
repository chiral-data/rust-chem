//! Per-format read and write options.
//!
//! `ReadFn`/`WriteFn` in [`crate::io::format`] thread these through to every
//! registered reader and writer. Most formats have nothing to configure yet
//! — the types exist so a format's signature only has to widen once, per
//! the plan already recorded on those type aliases before any format needed
//! a real option (#212).

/// One field per format that has a read option today (#324's `lammps` is
/// the first). A real field lands here the first time a reader needs one
/// to configure, rather than being invented ahead of that need.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReadOptions {
    pub lammps: LammpsReadOptions,
}

/// LAMMPS data's own read option: which `Atoms`-section column layout to
/// assume when the section header carries no recognized `# style` comment
/// and the column count alone is ambiguous (#324) -- `charge` and
/// `molecular` (which also covers `bond`/`angle`'s identical shape) share
/// a column count both with and without the optional image-flag triplet.
/// `None` (the default) means: resolve from the comment or an unambiguous
/// column count, refusing rather than guessing when neither settles it.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LammpsReadOptions {
    pub atom_style: Option<AtomStyle>,
}

/// A LAMMPS `atom_style`, restricted to the four column layouts this crate
/// models (#324). `molecular`, `bond` and `angle` share one shape --
/// `id, molecule-id, type, x, y, z` -- so `Molecular` covers all three;
/// this crate never needs to tell them apart. Any other real style
/// (`sphere`, `ellipsoid`, `electron`, ...) is a clear, refusing error
/// rather than a guessed-at column layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomStyle {
    /// `id, type, x, y, z`.
    Atomic,
    /// `id, type, charge, x, y, z`.
    Charge,
    /// `id, molecule-id, type, x, y, z` -- also `bond` and `angle`.
    Molecular,
    /// `id, molecule-id, type, charge, x, y, z`.
    Full,
}

/// One field per format that has a write option today.
///
/// No longer `Eq` since #325's `xtc: XtcWriteOptions` carries an `f32` --
/// nothing in this crate ever needed `WriteOptions` as a hash key or in an
/// exhaustive-equality context, only `assert_eq!`, which only needs
/// `PartialEq`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WriteOptions {
    pub sdf: SdfWriteOptions,
    pub xtc: XtcWriteOptions,
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

/// XTC's own write option: the coordinate precision, in inverse nanometres
/// (#325) -- a position rounds to the nearest `1.0 / precision` nm. `1000.0`
/// (0.001 nm resolution) is the real format's own documented default and
/// what every GROMACS-written file uses unless told otherwise.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XtcWriteOptions {
    pub precision: f32,
}

impl Default for XtcWriteOptions {
    fn default() -> Self {
        Self { precision: 1000.0 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_write_options_pick_v2000() {
        assert_eq!(WriteOptions::default().sdf.version, MolfileVersion::V2000);
    }

    #[test]
    fn test_default_xtc_precision_is_the_formats_own_documented_default() {
        assert_eq!(WriteOptions::default().xtc.precision, 1000.0);
    }
}
