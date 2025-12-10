//pub mod bitvec;
//pub mod errors;
//pub mod fingerprint;
//pub mod morgan;
//pub mod tanimoto;
//
//pub use bitvec::BitVec;
//pub use errors::FingerprintError;
//pub use morgan::{
//    morgan_fingerprint, morgan_fingerprint_with_args, morgan_fingerprint_with_info, BitInfo,
//    MorganArguments, MorganOptions,
//};
//pub use tanimoto::tanimoto_similarity;

pub mod bitvec;
pub mod errors;
pub mod fingerprint;
pub mod hash;
pub mod morgan;
pub mod tanimoto;

pub use morgan::{MorganArguments, MorganFingerprint};
