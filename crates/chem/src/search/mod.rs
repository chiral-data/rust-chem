//! Fingerprint generation and similarity search, across CPU and GPU backends.
//!
//! [`FingerprintSearch`] picks a backend and runs the work on it: Morgan
//! fingerprints from `chem::fp` or `chem::gpu`, Tanimoto similarity from either,
//! with the GPU used where one is usable and the CPU otherwise.
//!
//! # Why this is a crate rather than part of `chem::fp`
//!
//! `chem::fp`'s library is pure CPU and depends on nothing heavy. Putting the
//! backend selection there would make `wgpu` a permanent dependency of anyone
//! who only wants to fingerprint on a CPU — which is every wasm build and every
//! machine without a usable device. Keeping the orchestration separate lets
//! `chem::fp` stay small and lets a consumer opt into the GPU by depending on
//! this instead.
//!
//! # Why the API is async
//!
//! Not for the GUI's benefit. `wgpu` device creation and buffer readback are
//! inherently async, and on the web they cannot block at all — the browser has
//! no thread to park. A native caller that wants to block does so at its own
//! edge with `pollster::block_on`, which is one line and keeps a single code
//! path for both.
//!
//! # Where the GPU is *not* used
//!
//! Similarity is computed on the GPU only when a device initialised and the
//! caller has not forced the CPU. Both paths must agree, which is what
//! `chem::fp`'s parity tests check; the CPU path is not a fallback of last
//! resort but the reference the GPU one is measured against.

use crate::core::molecule::Molecule;
use crate::fp::morgan::MorganFingerprint;
use crate::fp::tanimoto::tanimoto_similarity;
#[cfg(feature = "gpu")]
use crate::gpu::{GpuMorganFingerprint, GpuTanimoto};

mod gpu_state;
use bitvec::prelude::BitVec;
use gpu_state::GpuState;

#[derive(Clone, Debug)]
pub struct SearchResult {
    pub index: usize,
    pub similarity: f64,
}

/// Cheaply `Clone` (see [`GpuMorganFingerprint`]/[`GpuTanimoto`]) — the
/// wasm32 async call sites in `chem-app`'s app module clone a snapshot to move into
/// `spawn_local` rather than sharing `&mut` access across the await
/// boundary. Cloned snapshots don't propagate their `gpu_targets` cache
/// back, so those call sites re-upload the target dataset per search rather
/// than reusing the cache #16 added — an accepted, wasm32-only tradeoff.
#[derive(Clone)]
pub struct FingerprintSearch {
    /// Every GPU handle, in one field so the struct has the same shape with
    /// and without the `gpu` feature. See [`gpu_state`].
    gpu: GpuState,
    use_gpu: bool,
    /// Set when a GPU init attempt (initial or retried) fails, so callers
    /// can distinguish "on CPU because GPU was never available/failed" from
    /// "on CPU by choice" (see [`Self::force_cpu`]) — e.g. to show a
    /// distinct error indicator instead of the plain CPU-mode one. Cleared
    /// on any successful init (see [`Self::install_gpu`]).
    gpu_init_error: Option<String>,
}

impl FingerprintSearch {
    // On native, GPU init can block the calling thread — it's only ever
    // called from a plain (non-UI) thread at startup. On wasm32, `new()`
    // runs inside eframe's synchronous app-construction callback, so it
    // can't block OR run the async adapter/device request itself; start
    // CPU-only and upgrade to GPU shortly after via `try_init_gpu_async` +
    // `install_gpu`, which `chem-app`'s app module drives with `spawn_local` and
    // polls once per frame (mirroring the file-load pattern from #37).
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    pub fn new() -> Self {
        let mut search = Self::new_cpu_only();
        if let Err(e) = search.retry_gpu_init() {
            log::warn!("GPU initialization failed, using CPU fallback: {e}");
        }
        search
    }

    /// wasm32 cannot block the browser's only thread, so it starts CPU-only and
    /// upgrades through [`Self::try_init_gpu_async`]. A build without the `gpu`
    /// feature has nothing to upgrade to, and records that as the reason —
    /// otherwise the downgrade is invisible and someone benchmarks this crate
    /// concluding the GPU work does not exist.
    /// wasm32 cannot block the browser's only thread, so it starts CPU-only and
    /// upgrades through [`Self::try_init_gpu_async`].
    #[cfg(all(feature = "gpu", target_arch = "wasm32"))]
    pub fn new() -> Self {
        log::info!("Web build: starting CPU-only, upgrading to GPU asynchronously if available");
        Self::new_cpu_only()
    }

    /// Built without the `gpu` feature, so there is no device to find. The
    /// reason is recorded rather than left implicit: otherwise the downgrade is
    /// invisible, and someone benchmarks this crate concluding the GPU work
    /// does not exist. [`Self::gpu_init_error`] reports it, and the CLI prints
    /// it.
    #[cfg(not(feature = "gpu"))]
    pub fn new() -> Self {
        Self {
            gpu_init_error: Some("built without the `gpu` feature".to_string()),
            ..Self::new_cpu_only()
        }
    }

    /// Bypasses GPU init entirely and always uses the CPU fingerprinting/
    /// search path. Used on wasm32 (see [`Self::new`]) and by tests that need
    /// deterministic CPU-path coverage regardless of whether the machine
    /// running them has a GPU; also public for callers that want to force
    /// CPU mode on native (e.g. benchmarking against the GPU path).
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    pub fn new_cpu_only() -> Self {
        Self {
            gpu: GpuState::new(),
            use_gpu: false,
            gpu_init_error: None,
        }
    }

    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    fn init_gpu() -> anyhow::Result<(GpuMorganFingerprint, GpuTanimoto)> {
        let morgan = GpuMorganFingerprint::new()?;
        let tanimoto = GpuTanimoto::new()?;
        Ok((morgan, tanimoto))
    }

    /// (Re)attempts GPU init synchronously, e.g. from a "retry GPU" UI
    /// action. Installs and switches to it on success; on failure, records
    /// the reason (see [`Self::gpu_init_error`]) and leaves the instance on
    /// whatever it was using before (CPU, unless a GPU was already active).
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    pub fn retry_gpu_init(&mut self) -> Result<(), String> {
        match Self::init_gpu() {
            Ok((m, t)) => {
                log::info!("GPU acceleration enabled");
                self.install_gpu(m, t);
                Ok(())
            }
            Err(e) => {
                let msg = e.to_string();
                self.gpu_init_error = Some(msg.clone());
                Err(msg)
            }
        }
    }

    /// Attempts GPU init without blocking the calling thread — for wasm32,
    /// where blocking (e.g. from `new()`, or a "retry GPU" UI action) would
    /// deadlock the browser's single JS thread. An associated function
    /// rather than a method since it doesn't touch any existing instance —
    /// callers install the result themselves via [`Self::install_gpu`] on
    /// success, or [`Self::record_gpu_init_failure`] on failure.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    #[cfg(feature = "gpu")]
    pub async fn try_init_gpu_async() -> Result<(GpuMorganFingerprint, GpuTanimoto), String> {
        let morgan = GpuMorganFingerprint::new_async().await.map_err(|e| {
            let msg = format!("Morgan: {}", e);
            log::warn!("Async GPU initialization failed ({})", msg);
            msg
        })?;
        let tanimoto = GpuTanimoto::new_async().await.map_err(|e| {
            let msg = format!("Tanimoto: {}", e);
            log::warn!("Async GPU initialization failed ({})", msg);
            msg
        })?;
        log::info!("GPU acceleration enabled (async)");
        Ok((morgan, tanimoto))
    }

    /// Upgrades an instance in place once GPU init succeeds (see
    /// [`Self::try_init_gpu_async`]/[`Self::retry_gpu_init`]). Drops any
    /// cached GPU target upload — there wasn't one while CPU-only, but this
    /// keeps the invariant explicit rather than assuming it — and clears
    /// [`Self::gpu_init_error`].
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    #[cfg(feature = "gpu")]
    pub fn install_gpu(&mut self, gpu_morgan: GpuMorganFingerprint, gpu_tanimoto: GpuTanimoto) {
        self.gpu.morgan = Some(gpu_morgan);
        self.gpu.tanimoto = Some(gpu_tanimoto);
        self.use_gpu = true;
        self.gpu.clear_targets();
        self.gpu_init_error = None;
    }

    /// Records that a GPU init attempt (e.g. the wasm32 async one kicked off
    /// by `chem-app`'s app module) failed, for [`Self::gpu_init_error`] to report.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    pub fn record_gpu_init_failure(&mut self, error: String) {
        self.gpu_init_error = Some(error);
    }

    /// Reason the most recent GPU init attempt failed, if any — `None` if
    /// GPU is currently active, or if no attempt has completed yet (e.g.
    /// wasm32's async init is still in flight).
    pub fn gpu_init_error(&self) -> Option<&str> {
        self.gpu_init_error.as_deref()
    }

    /// True if a GPU context was successfully initialized at some point,
    /// whether or not it's the one currently in use (see [`Self::force_cpu`]/
    /// [`Self::force_gpu`]).
    pub fn has_gpu_available(&self) -> bool {
        self.gpu.is_available()
    }

    /// Switches to the CPU fingerprinting/search path without discarding an
    /// already-initialized GPU context, so a later [`Self::force_gpu`] call
    /// can switch back instantly instead of re-initializing.
    pub fn force_cpu(&mut self) {
        self.use_gpu = false;
    }

    /// Switches back to the GPU path if one was already initialized.
    /// Returns `false` (no-op) if no GPU context is available yet — the
    /// caller should trigger a (re)init attempt instead (see
    /// [`Self::retry_gpu_init`]/[`Self::try_init_gpu_async`]).
    pub fn force_gpu(&mut self) -> bool {
        if self.has_gpu_available() {
            self.use_gpu = true;
            true
        } else {
            false
        }
    }

    // Unused on wasm32: chem_app::app only calls the _async twin there
    // (blocking would deadlock the browser's single JS thread — see
    // Self::new). Still part of the public API and exercised by tests.
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    /// Generates a fingerprint, on the GPU where one is available.
    ///
    /// Async on every platform: wasm32 can't block, and native callers go
    /// through `task::Task`, which does the blocking in one place rather than
    /// each operation carrying a sync twin of itself.
    pub async fn generate_fingerprint_async(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        #[cfg(feature = "gpu")]
        if self.use_gpu && self.gpu.has_morgan() {
            return self
                .generate_fingerprint_gpu_async(mol, radius, fp_size)
                .await;
        }
        self.generate_fingerprint_cpu(mol, radius, fp_size)
    }

    /// Batch form of [`Self::generate_fingerprint_async`].
    pub async fn generate_fingerprints_batch_async(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        #[cfg(feature = "gpu")]
        if self.use_gpu && self.gpu.has_morgan() {
            return self
                .generate_fingerprints_gpu_batch_async(molecules, radius, fp_size)
                .await;
        }
        self.generate_fingerprints_cpu_batch(molecules, radius, fp_size)
    }

    /// Upload `target_fps` to the GPU once, so subsequent
    /// [`Self::search_async`] calls reuse the same buffer instead of
    /// re-uploading the whole dataset on every query. Call this whenever the
    /// dataset changes; `search_async` will still lazily upload on its own if
    /// this was skipped.
    pub fn set_target_dataset(&mut self, target_fps: &[BitVec]) -> anyhow::Result<()> {
        #[cfg(feature = "gpu")]
        if self.use_gpu && self.gpu.has_tanimoto() && !target_fps.is_empty() {
            let gpu = self.gpu.tanimoto.as_ref().unwrap();
            let fp_words = Self::bitvec_to_words(&target_fps[0]).len();
            let flattened: Vec<u32> = target_fps.iter().flat_map(Self::bitvec_to_words).collect();
            self.gpu.targets = Some(gpu.upload_targets(&flattened, fp_words)?);
            return Ok(());
        }
        // Named so the argument is used in every configuration.
        let _ = target_fps;
        self.gpu.clear_targets();
        Ok(())
    }

    /// Drop the cached GPU target upload, e.g. when the dataset is cleared.
    pub fn invalidate_target_dataset(&mut self) {
        self.gpu.clear_targets();
    }

    /// Ranks `target_fps` against `query_fp`, on the GPU where one is
    /// available. Async for the same reason as
    /// [`Self::generate_fingerprint_async`].
    pub async fn search_async(
        &mut self,
        query_fp: &BitVec,
        target_fps: &[BitVec],
        top_k: usize,
    ) -> anyhow::Result<Vec<SearchResult>> {
        #[cfg(feature = "gpu")]
        let similarities = if self.use_gpu && self.gpu.has_tanimoto() {
            self.compute_similarities_gpu_async(query_fp, target_fps)
                .await?
        } else {
            self.compute_similarities_cpu(query_fp, target_fps)?
        };
        #[cfg(not(feature = "gpu"))]
        let similarities = self.compute_similarities_cpu(query_fp, target_fps)?;

        let mut results: Vec<SearchResult> = similarities
            .into_iter()
            .enumerate()
            .map(|(index, similarity)| SearchResult { index, similarity })
            .collect();

        results.sort_by(|a, b| b.similarity.partial_cmp(&a.similarity).unwrap());
        results.truncate(top_k);

        Ok(results)
    }

    #[cfg(feature = "gpu")]
    async fn generate_fingerprint_gpu_async(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        let gpu = self.gpu.morgan.as_ref().unwrap();
        let fp_words = gpu
            .generate_fingerprints_batch_async(
                std::slice::from_ref(mol),
                radius,
                fp_size,
                false,
                true,
                false,
            )
            .await?;

        Ok(Self::words_to_bitvec(&fp_words[0], fp_size as usize))
    }

    #[cfg(feature = "gpu")]
    async fn generate_fingerprints_gpu_batch_async(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        let gpu = self.gpu.morgan.as_ref().unwrap();
        let fp_words_batch = gpu
            .generate_fingerprints_batch_async(molecules, radius, fp_size, false, true, false)
            .await?;

        Ok(fp_words_batch
            .iter()
            .map(|words| Self::words_to_bitvec(words, fp_size as usize))
            .collect())
    }

    #[cfg(feature = "gpu")]
    async fn compute_similarities_gpu_async(
        &mut self,
        query_fp: &BitVec,
        target_fps: &[BitVec],
    ) -> anyhow::Result<Vec<f64>> {
        // Word-count divisibility alone can't catch a wrong-but-divisible target
        // size (e.g. 1024-bit targets against a 2048-bit query), since the GPU
        // path only sees flattened u32 words, not each fingerprint's bit length.
        // Check exact bit lengths here, same as the CPU path's tanimoto_similarity.
        if let Some(mismatch_idx) = target_fps.iter().position(|fp| fp.len() != query_fp.len()) {
            anyhow::bail!(
                "Target fingerprint size mismatch at index {}: expected {} bits, got {} bits",
                mismatch_idx,
                query_fp.len(),
                target_fps[mismatch_idx].len()
            );
        }

        // No targets means no results — mirrors the empty-buffer guard the
        // GPU path itself uses, and sidesteps caching an empty upload.
        if target_fps.is_empty() {
            return Ok(Vec::new());
        }

        let query_words = Self::bitvec_to_words(query_fp);

        // Reuse the cached upload when it already matches this target set's
        // shape; only re-upload (the whole point of #16) when the caller
        // never called `set_target_dataset` or the dataset actually changed.
        let cache_valid =
            self.gpu.targets.as_ref().is_some_and(|t| {
                t.count() == target_fps.len() && t.fp_words() == query_words.len()
            });
        if !cache_valid {
            self.set_target_dataset(target_fps)?;
        }

        let gpu = self.gpu.tanimoto.as_ref().unwrap();
        let targets = self
            .gpu
            .targets
            .as_ref()
            .expect("set_target_dataset populates gpu_targets for a non-empty dataset");
        let similarities_f32 = gpu
            .compute_single_query_against_async(&query_words, targets)
            .await?;
        Ok(similarities_f32.into_iter().map(|s| s as f64).collect())
    }

    fn generate_fingerprint_cpu(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        MorganFingerprint::get_fingerprint_as_bitvec(
            mol, radius, fp_size, None, None, false, true, false,
        )
        .map_err(|e| anyhow::anyhow!("CPU fingerprint generation failed: {}", e))
    }

    fn generate_fingerprints_cpu_batch(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        molecules
            .iter()
            .map(|mol| self.generate_fingerprint_cpu(mol, radius, fp_size))
            .collect()
    }

    fn compute_similarities_cpu(
        &self,
        query_fp: &BitVec,
        target_fps: &[BitVec],
    ) -> anyhow::Result<Vec<f64>> {
        target_fps
            .iter()
            .map(|target| {
                tanimoto_similarity(query_fp, target)
                    .map_err(|e| anyhow::anyhow!("Tanimoto calculation failed: {}", e))
            })
            .collect()
    }

    #[cfg(feature = "gpu")]
    fn words_to_bitvec(words: &[u32], fp_size: usize) -> BitVec {
        let mut bv = BitVec::repeat(false, fp_size);
        for (word_idx, &word) in words.iter().enumerate() {
            for bit_idx in 0..32 {
                if word & (1 << bit_idx) != 0 {
                    let pos = word_idx * 32 + bit_idx;
                    if pos < fp_size {
                        bv.set(pos, true);
                    }
                }
            }
        }
        bv
    }

    #[cfg(feature = "gpu")]
    fn bitvec_to_words(bv: &BitVec) -> Vec<u32> {
        let num_words = bv.len().div_ceil(32);
        let mut words = vec![0u32; num_words];

        for (i, bit) in bv.iter().enumerate() {
            if *bit {
                let word_idx = i / 32;
                let bit_idx = i % 32;
                words[word_idx] |= 1 << bit_idx;
            }
        }

        words
    }

    pub fn is_using_gpu(&self) -> bool {
        self.use_gpu
    }
}

impl Default for FingerprintSearch {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The tests are sync, so they do the blocking the app delegates to
    /// `task::Task`.
    fn block<T>(future: impl std::future::Future<Output = T>) -> T {
        pollster::block_on(future)
    }

    /// One GPU-initialised search, built once and cloned per test.
    ///
    /// Every test that called `FingerprintSearch::new()` requested its own
    /// adapter and device, and two of those racing against one physical GPU
    /// hangs indefinitely on at least one driver — the deadlock #19 diagnosed
    /// for `chem::gpu` and #134 for here. `chem::gpu` solved it by sharing a context
    /// through a `OnceLock`; this is the same move one layer up.
    ///
    /// Clones rather than lending a reference because these tests mutate the
    /// search — `set_target_dataset` fills its upload cache. Cloning is what
    /// makes that safe *and* cheap: `wgpu`'s device, queue and pipelines are
    /// refcounted handles, so a clone shares the device instead of asking for
    /// another one, while getting a cache of its own.
    ///
    /// `None` means no usable GPU, which is not a failure — it is what every
    /// CI runner reports.
    ///
    /// Native only, and there is nothing to share on wasm32 anyway:
    /// `FingerprintSearch::new()` starts CPU-only there, since it cannot block
    /// the browser's single thread to request a device, so a GPU-path test has
    /// nothing to run against. A `static` would not compile either — wgpu's
    /// WebGPU types hold a `RefCell` and raw pointers, the web being
    /// single-threaded, so they are not `Sync`.
    #[cfg(not(target_arch = "wasm32"))]
    fn shared_gpu_search() -> Option<FingerprintSearch> {
        static SEARCH: std::sync::OnceLock<Option<FingerprintSearch>> = std::sync::OnceLock::new();
        SEARCH
            .get_or_init(|| {
                let search = FingerprintSearch::new();
                if search.is_using_gpu() {
                    Some(search)
                } else {
                    println!("GPU not available; GPU-path tests will skip");
                    None
                }
            })
            .clone()
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn test_gpu_search_rejects_mismatched_fingerprint_size() {
        let Some(mut search) = shared_gpu_search() else {
            return;
        };

        // Two 1024-bit targets flatten to exactly the same word count as one
        // 2048-bit query (64 words), so a word-count-divisibility-only check
        // would silently accept this and treat the two targets as one —
        // exactly the "2*N-1 index" corruption the size check must prevent.
        let query_fp = BitVec::repeat(false, 2048);
        let mismatched_targets = vec![BitVec::repeat(false, 1024), BitVec::repeat(false, 1024)];

        let result = block(search.search_async(&query_fp, &mismatched_targets, 2));
        assert!(
            result.is_err(),
            "expected a size-mismatch error, got {:?}",
            result
        );
    }

    // The two tests above only assert once GPU acceleration is actually in
    // use, so on any machine/CI runner without a usable GPU adapter they
    // silently skip and verify nothing. The CPU path is what every wasm32
    // (web) build and any GPU-less native machine actually runs, so give it
    // its own deterministic coverage that never skips.

    #[test]
    fn test_cpu_search_rejects_mismatched_fingerprint_size() {
        let mut search = FingerprintSearch::new_cpu_only();
        assert!(!search.is_using_gpu());

        let query_fp = BitVec::repeat(false, 2048);
        let mismatched_targets = vec![BitVec::repeat(false, 1024), BitVec::repeat(false, 1024)];

        let result = block(search.search_async(&query_fp, &mismatched_targets, 2));
        assert!(
            result.is_err(),
            "expected a size-mismatch error, got {:?}",
            result
        );
    }

    #[test]
    fn test_cpu_search_similarity_is_correct() {
        let mut search = FingerprintSearch::new_cpu_only();

        let query_fp = BitVec::repeat(true, 2048);

        let dataset_a = vec![BitVec::repeat(true, 2048)];
        let results_a = block(search.search_async(&query_fp, &dataset_a, 1)).unwrap();
        assert!(
            (results_a[0].similarity - 1.0).abs() < 1e-6,
            "expected perfect self-similarity, got {:?}",
            results_a
        );

        let dataset_b = vec![BitVec::repeat(false, 2048)];
        let results_b = block(search.search_async(&query_fp, &dataset_b, 1)).unwrap();
        assert!(
            results_b[0].similarity.abs() < 1e-6,
            "expected zero similarity against an all-zero target, got {:?}",
            results_b
        );
    }

    #[test]
    fn test_cpu_fingerprint_generation_produces_requested_size() {
        let search = FingerprintSearch::new_cpu_only();
        let mol = crate::io::smiles::parse_smiles("c1ccccc1").expect("valid SMILES");

        let fp = block(search.generate_fingerprint_async(&mol, 2, 2048))
            .expect("CPU fingerprint generation should succeed");
        assert_eq!(fp.len(), 2048);

        let batch =
            block(search.generate_fingerprints_batch_async(std::slice::from_ref(&mol), 2, 2048))
                .expect("CPU batch fingerprint generation should succeed");
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].len(), 2048);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn test_gpu_search_reflects_new_dataset_after_reupload() {
        let Some(mut search) = shared_gpu_search() else {
            return;
        };

        // Same shape (1 target, 64 words) for both datasets, so a cache keyed
        // only on count/fp_words would wrongly keep serving dataset A's
        // upload after dataset B replaces it via `set_target_dataset`.
        let query_fp = BitVec::repeat(true, 2048);

        let dataset_a = vec![BitVec::repeat(true, 2048)];
        search.set_target_dataset(&dataset_a).unwrap();
        let results_a = block(search.search_async(&query_fp, &dataset_a, 1)).unwrap();
        assert!(
            (results_a[0].similarity - 1.0).abs() < 1e-6,
            "expected perfect self-similarity against dataset A, got {:?}",
            results_a
        );

        let dataset_b = vec![BitVec::repeat(false, 2048)];
        search.set_target_dataset(&dataset_b).unwrap();
        let results_b = block(search.search_async(&query_fp, &dataset_b, 1)).unwrap();
        assert!(
            results_b[0].similarity.abs() < 1e-6,
            "expected zero similarity against dataset B, got a stale result from dataset A: {:?}",
            results_b
        );
    }
}
