//! Physical unit conversions shared across formats.

/// Nanometres to Ångström. GROMACS's whole ecosystem (GRO, TOP, TRR, XTC)
/// uses nm; this crate's data model, and every other format registered so
/// far, is Å.
pub(crate) const NM_TO_ANGSTROM: f64 = 10.0;

/// Bohr radii to Ångström (CODATA). Gaussian's CUBE format (#331) states
/// coordinates and grid vectors in Bohr by default.
pub(crate) const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;
