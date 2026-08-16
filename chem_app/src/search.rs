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

pub struct FingerprintSearch {
    gpu_morgan: Option<GpuMorganFingerprint>,
    gpu_tanimoto: Option<GpuTanimoto>,
    use_gpu: bool,
    /// Target dataset already uploaded to the GPU, reused across searches
    /// so only the (small) query fingerprint round-trips per search instead
    /// of re-uploading the whole dataset every time.
    gpu_targets: Option<GpuTargetSet>,
}

impl FingerprintSearch {
    pub fn new() -> Self {
        let (gpu_morgan, gpu_tanimoto, use_gpu) = match Self::init_gpu() {
            Ok((m, t)) => {
                log::info!("GPU acceleration enabled");
                (Some(m), Some(t), true)
            }
            Err(e) => {
                log::warn!("GPU initialization failed, using CPU fallback: {}", e);
                (None, None, false)
            }
        };

        Self {
            gpu_morgan,
            gpu_tanimoto,
            use_gpu,
            gpu_targets: None,
        }
    }

    fn init_gpu() -> anyhow::Result<(GpuMorganFingerprint, GpuTanimoto)> {
        let morgan = GpuMorganFingerprint::new()?;
        let tanimoto = GpuTanimoto::new()?;
        Ok((morgan, tanimoto))
    }

    pub fn generate_fingerprint(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        if self.use_gpu && self.gpu_morgan.is_some() {
            self.generate_fingerprint_gpu(mol, radius, fp_size)
        } else {
            self.generate_fingerprint_cpu(mol, radius, fp_size)
        }
    }

    pub fn generate_fingerprints_batch(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        if self.use_gpu && self.gpu_morgan.is_some() {
            self.generate_fingerprints_gpu_batch(molecules, radius, fp_size)
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

    pub fn search(
        &mut self,
        query_fp: &BitVec,
        target_fps: &[BitVec],
        top_k: usize,
    ) -> anyhow::Result<Vec<SearchResult>> {
        let similarities = if self.use_gpu && self.gpu_tanimoto.is_some() {
            self.compute_similarities_gpu(query_fp, target_fps)?
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

    fn generate_fingerprint_gpu(
        &self,
        mol: &Molecule,
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<BitVec> {
        let gpu = self.gpu_morgan.as_ref().unwrap();
        let fp_words = gpu.generate_fingerprints_batch(
            std::slice::from_ref(mol),
            radius,
            fp_size,
            false,
            true,
            false,
        )?;

        Ok(Self::words_to_bitvec(&fp_words[0], fp_size as usize))
    }

    fn generate_fingerprints_gpu_batch(
        &self,
        molecules: &[Molecule],
        radius: u32,
        fp_size: u32,
    ) -> anyhow::Result<Vec<BitVec>> {
        let gpu = self.gpu_morgan.as_ref().unwrap();
        let fp_words_batch =
            gpu.generate_fingerprints_batch(molecules, radius, fp_size, false, true, false)?;

        Ok(fp_words_batch
            .iter()
            .map(|words| Self::words_to_bitvec(words, fp_size as usize))
            .collect())
    }

    fn compute_similarities_gpu(
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
        let similarities_f32 = gpu.compute_single_query_against(&query_words, targets)?;
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
