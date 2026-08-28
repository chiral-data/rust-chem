//! Exit codes.
//!
//! Distinguishable rather than merely non-zero, because the caller of a batch
//! job needs to tell "the file was bad" from "the file was fine and nothing
//! matched" — those want different responses, and a bare `1` collapses them.
//!
//! Codes start at 2: the shell uses 126 and 127 for its own failures, and 1 is
//! what a panic or an unhandled error produces, so leaving it unclaimed keeps
//! "the tool broke" distinct from "the tool ran and says no".

/// Everything worked.
pub const OK: i32 = 0;

/// Nothing in the input parsed — the file is unusable rather than partly bad.
pub const NO_INPUT: i32 = 2;

/// The input parsed but the operation produced nothing: a search with no hits,
/// a filter that matched nothing. Not an error, and callers who treat it as one
/// can, but distinct from a failure.
///
/// Unused so far — `chem search` in #140 is its first caller. It is declared
/// here anyway, because the point of this module is that the codes are decided
/// once, in one table, rather than invented by whichever command needs one
/// first and then quietly reused with a different meaning.
#[allow(dead_code)]
pub const EMPTY_RESULT: i32 = 3;

/// Some records were skipped, and the caller asked for that to be fatal.
pub const PARTIAL: i32 = 4;
