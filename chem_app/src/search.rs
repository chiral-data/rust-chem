use bitvec::prelude::BitVec;
use chemcore::molecule::Molecule;
use chemfp::morgan::MorganFingerprint;
use chemfp::tanimoto::tanimoto_similarity;
use chemgpu::{GpuMorganFingerprint, GpuTanimoto, GpuTargetSet};

#[derive(Clone, Debug)]
pub struct SearchResult {
    pub index: usize,
    pub similarity: f64,
}

/// Cheaply `Clone` (see [`GpuMorganFingerprint`]/[`GpuTanimoto`]) — the
/// wasm32 async call sites in `chem_app::app` clone a snapshot to move into
/// `spawn_local` rather than sharing `&mut` access across the await
/// boundary. Cloned snapshots don't propagate their `gpu_targets` cache
/// back, so those call sites re-upload the target dataset per search rather
/// than reusing the cache #16 added — an accepted, wasm32-only tradeoff.
#[derive(Clone)]
pub struct FingerprintSearch {
    gpu_morgan: Option<GpuMorganFingerprint>,
    gpu_tanimoto: Option<GpuTanimoto>,
    use_gpu: bool,
    /// Target dataset already uploaded to the GPU, reused across searches
    /// so only the (small) query fingerprint round-trips per search instead
    /// of re-uploading the whole dataset every time.
    gpu_targets: Option<GpuTargetSet>,
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
    // `install_gpu`, which `chem_app::app` drives with `spawn_local` and
    // polls once per frame (mirroring the file-load pattern from #37).
    #[cfg(target_arch = "wasm32")]
    pub fn new() -> Self {
        log::info!("Web build: starting CPU-only, upgrading to GPU asynchronously if available");
        Self::new_cpu_only()
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn new() -> Self {
        let mut search = Self::new_cpu_only();
        if let Err(e) = search.retry_gpu_init() {
            log::warn!("GPU initialization failed, using CPU fallback: {}", e);
        }
        search
    }

    /// Bypasses GPU init entirely and always uses the CPU fingerprinting/
    /// search path. Used on wasm32 (see [`Self::new`]) and by tests that need
    /// deterministic CPU-path coverage regardless of whether the machine
    /// running them has a GPU; also public for callers that want to force
    /// CPU mode on native (e.g. benchmarking against the GPU path).
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    pub fn new_cpu_only() -> Self {
        Self {
            gpu_morgan: None,
            gpu_tanimoto: None,
            use_gpu: false,
            gpu_targets: None,
            gpu_init_error: None,
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn init_gpu() -> anyhow::Result<(GpuMorganFingerprint, GpuTanimoto)> {
        let morgan = GpuMorganFingerprint::new()?;
        let tanimoto = GpuTanimoto::new()?;
        Ok((morgan, tanimoto))
    }

    /// (Re)attempts GPU init synchronously, e.g. from a "retry GPU" UI
    /// action. Installs and switches to it on success; on failure, records
    /// the reason (see [`Self::gpu_init_error`]) and leaves the instance on
    /// whatever it was using before (CPU, unless a GPU was already active).
    #[cfg(not(target_arch = "wasm32"))]
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
    pub fn install_gpu(&mut self, gpu_morgan: GpuMorganFingerprint, gpu_tanimoto: GpuTanimoto) {
        self.gpu_morgan = Some(gpu_morgan);
        self.gpu_tanimoto = Some(gpu_tanimoto);
        self.use_gpu = true;
        self.gpu_targets = None;
        self.gpu_init_error = None;
    }

    /// Records that a GPU init attempt (e.g. the wasm32 async one kicked off
    /// by `chem_app::app`) failed, for [`Self::gpu_init_error`] to report.
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
        self.gpu_morgan.is_some() && self.gpu_tanimoto.is_some()
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
    pub fn generate_fingerprint(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        pollster::block_on(self.generate_fingerprint_async(mol, radius, fp_size))
    }

    /// Async twin of [`Self::generate_fingerprint`] for callers that can't
    /// block the current thread (namely wasm32).
    pub async fn generate_fingerprint_async(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        if self.use_gpu && self.gpu_morgan.is_some() {
            self.generate_fingerprint_gpu_async(mol, radius, fp_size)
                .await
        } else {
            self.generate_fingerprint_cpu(mol, radius, fp_size)
        }
    }

    // Unused on wasm32 — see the note on generate_fingerprint above.
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    pub fn generate_fingerprints_batch(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        pollster::block_on(self.generate_fingerprints_batch_async(molecules, radius, fp_size))
    }

    /// Async twin of [`Self::generate_fingerprints_batch`] for callers that
    /// can't block the current thread (namely wasm32).
    pub async fn generate_fingerprints_batch_async(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        if self.use_gpu && self.gpu_morgan.is_some() {
            self.generate_fingerprints_gpu_batch_async(molecules, radius, fp_size)
                .await
        } else {
            self.generate_fingerprints_cpu_batch(molecules, radius, fp_size)
        }
    }

    /// Upload `target_fps` to the GPU once, so subsequent [`Self::search`]
    /// calls reuse the same buffer instead of re-uploading the whole
    /// dataset on every query. Call this whenever the dataset changes;
    /// `search` will still lazily upload on its own if this was skipped.
    pub fn set_target_dataset(&mut self, target_fps: &[BitVec]) -> anyhow::Result<()> {
        if !(self.use_gpu && self.gpu_tanimoto.is_some()) || target_fps.is_empty() {
            self.gpu_targets = None;
            return Ok(());
        }

        let gpu = self.gpu_tanimoto.as_ref().unwrap();
        let fp_words = Self::bitvec_to_words(&target_fps[0]).len();
        let flattened: Vec<u32> = target_fps.iter().flat_map(Self::bitvec_to_words).collect();

        self.gpu_targets = Some(gpu.upload_targets(&flattened, fp_words)?);
        Ok(())
    }

    /// Drop the cached GPU target upload, e.g. when the dataset is cleared.
    pub fn invalidate_target_dataset(&mut self) {
        self.gpu_targets = None;
    }

    // Unused on wasm32 — see the note on generate_fingerprint above.
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    pub fn search(
        &mut self,
        query_fp: &BitVec,
        target_fps: &[BitVec],
        top_k: usize,
    ) -> anyhow::Result<Vec<SearchResult>> {
        pollster::block_on(self.search_async(query_fp, target_fps, top_k))
    }

    /// Async twin of [`Self::search`] for callers that can't block the
    /// current thread (namely wasm32).
    pub async fn search_async(
        &mut self,
        query_fp: &BitVec,
        target_fps: &[BitVec],
        top_k: usize,
    ) -> anyhow::Result<Vec<SearchResult>> {
        let similarities = if self.use_gpu && self.gpu_tanimoto.is_some() {
            self.compute_similarities_gpu_async(query_fp, target_fps)
                .await?
        } else {
            self.compute_similarities_cpu(query_fp, target_fps)?
        };

        let mut results: Vec<SearchResult> = similarities
            .into_iter()
            .enumerate()
            .map(|(index, similarity)| SearchResult { index, similarity })
            .collect();

        results.sort_by(|a, b| b.similarity.partial_cmp(&a.similarity).unwrap());
        results.truncate(top_k);

        Ok(results)
    }

    async fn generate_fingerprint_gpu_async(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        let gpu = self.gpu_morgan.as_ref().unwrap();
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

    async fn generate_fingerprints_gpu_batch_async(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        let gpu = self.gpu_morgan.as_ref().unwrap();
        let fp_words_batch = gpu
            .generate_fingerprints_batch_async(molecules, radius, fp_size, false, true, false)
            .await?;

        Ok(fp_words_batch
            .iter()
            .map(|words| Self::words_to_bitvec(words, fp_size as usize))
            .collect())
    }

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
        let cache_valid = self
            .gpu_targets
            .as_ref()
            .is_some_and(|t| t.count() == target_fps.len() && t.fp_words() == query_words.len());
        if !cache_valid {
            self.set_target_dataset(target_fps)?;
        }

        let gpu = self.gpu_tanimoto.as_ref().unwrap();
        let targets = self
            .gpu_targets
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

    #[test]
    fn test_gpu_search_rejects_mismatched_fingerprint_size() {
        let mut search = FingerprintSearch::new();
        if !search.is_using_gpu() {
            println!("GPU not available, skipping GPU-path size-mismatch test");
            return;
        }

        // Two 1024-bit targets flatten to exactly the same word count as one
        // 2048-bit query (64 words), so a word-count-divisibility-only check
        // would silently accept this and treat the two targets as one —
        // exactly the "2*N-1 index" corruption the size check must prevent.
        let query_fp = BitVec::repeat(false, 2048);
        let mismatched_targets = vec![BitVec::repeat(false, 1024), BitVec::repeat(false, 1024)];

        let result = search.search(&query_fp, &mismatched_targets, 2);
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

        let result = search.search(&query_fp, &mismatched_targets, 2);
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
        let results_a = search.search(&query_fp, &dataset_a, 1).unwrap();
        assert!(
            (results_a[0].similarity - 1.0).abs() < 1e-6,
            "expected perfect self-similarity, got {:?}",
            results_a
        );

        let dataset_b = vec![BitVec::repeat(false, 2048)];
        let results_b = search.search(&query_fp, &dataset_b, 1).unwrap();
        assert!(
            results_b[0].similarity.abs() < 1e-6,
            "expected zero similarity against an all-zero target, got {:?}",
            results_b
        );
    }

    #[test]
    fn test_cpu_fingerprint_generation_produces_requested_size() {
        let search = FingerprintSearch::new_cpu_only();
        let mol = chemio::smiles::parse_smiles("c1ccccc1").expect("valid SMILES");

        let fp = search
            .generate_fingerprint(&mol, 2, 2048)
            .expect("CPU fingerprint generation should succeed");
        assert_eq!(fp.len(), 2048);

        let batch = search
            .generate_fingerprints_batch(std::slice::from_ref(&mol), 2, 2048)
            .expect("CPU batch fingerprint generation should succeed");
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].len(), 2048);
    }

    #[test]
    fn test_gpu_search_reflects_new_dataset_after_reupload() {
        let mut search = FingerprintSearch::new();
        if !search.is_using_gpu() {
            println!("GPU not available, skipping GPU target-cache test");
            return;
        }

        // Same shape (1 target, 64 words) for both datasets, so a cache keyed
        // only on count/fp_words would wrongly keep serving dataset A's
        // upload after dataset B replaces it via `set_target_dataset`.
        let query_fp = BitVec::repeat(true, 2048);

        let dataset_a = vec![BitVec::repeat(true, 2048)];
        search.set_target_dataset(&dataset_a).unwrap();
        let results_a = search.search(&query_fp, &dataset_a, 1).unwrap();
        assert!(
            (results_a[0].similarity - 1.0).abs() < 1e-6,
            "expected perfect self-similarity against dataset A, got {:?}",
            results_a
        );

        let dataset_b = vec![BitVec::repeat(false, 2048)];
        search.set_target_dataset(&dataset_b).unwrap();
        let results_b = search.search(&query_fp, &dataset_b, 1).unwrap();
        assert!(
            results_b[0].similarity.abs() < 1e-6,
            "expected zero similarity against dataset B, got a stale result from dataset A: {:?}",
            results_b
        );
    }
}
