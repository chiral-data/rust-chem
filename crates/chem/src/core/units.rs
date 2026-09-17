//! Physical unit conversions shared across formats.

/// Nanometres to Ångström. GROMACS's whole ecosystem (GRO, TOP, TRR, XTC)
/// uses nm; this crate's data model, and every other format registered so
/// far, is Å.
pub(crate) const NM_TO_ANGSTROM: f64 = 10.0;
