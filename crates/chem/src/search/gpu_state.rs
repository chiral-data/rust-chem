//! The GPU-resident half of [`super::FingerprintSearch`], in two shapes.
//!
//! With the `gpu` feature, three refcounted wgpu handles. Without it, a
//! zero-sized struct whose every query answers "no".
//!
//! One field on `FingerprintSearch` either way, which is the point: the struct,
//! its `Clone`, `new_cpu_only`, and every CPU code path are written once and
//! compile identically in both configurations. The alternative — cfg-ing three
//! fields individually — forks the constructor and every method that touches
//! them.

#[cfg(feature = "gpu")]
mod imp {
    use crate::gpu::{GpuMorganFingerprint, GpuTanimoto, GpuTargetSet};

    #[derive(Clone, Default)]
    pub(crate) struct GpuState {
        pub(crate) morgan: Option<GpuMorganFingerprint>,
        pub(crate) tanimoto: Option<GpuTanimoto>,
        /// Target dataset already uploaded, reused across searches so only the
        /// small query fingerprint round-trips per search.
        pub(crate) targets: Option<GpuTargetSet>,
    }

    impl GpuState {
        pub(crate) fn new() -> Self {
            Self::default()
        }
        pub(crate) fn has_morgan(&self) -> bool {
            self.morgan.is_some()
        }
        pub(crate) fn has_tanimoto(&self) -> bool {
            self.tanimoto.is_some()
        }
        pub(crate) fn is_available(&self) -> bool {
            self.has_morgan() && self.has_tanimoto()
        }
        pub(crate) fn clear_targets(&mut self) {
            self.targets = None;
        }
    }
}

#[cfg(not(feature = "gpu"))]
mod imp {
    #[derive(Clone, Default)]
    pub(crate) struct GpuState;

    // Every caller of these sits behind `#[cfg(feature = "gpu")]`, so without
    // the feature they are unreachable — but they must still exist for the
    // shared code paths to name them.
    #[allow(dead_code)]
    impl GpuState {
        pub(crate) fn new() -> Self {
            Self
        }
        pub(crate) fn has_morgan(&self) -> bool {
            false
        }
        pub(crate) fn has_tanimoto(&self) -> bool {
            false
        }
        pub(crate) fn is_available(&self) -> bool {
            false
        }
        pub(crate) fn clear_targets(&mut self) {}
    }
}

pub(crate) use imp::GpuState;
