use bitvec::prelude::BitVec;
use chemcore::molecule::Molecule;
use chemfp::morgan::MorganFingerprint;
use chemfp::tanimoto::tanimoto_similarity;
use chemgpu::{GpuMorganFingerprint, GpuTanimoto};

#[derive(Clone, Debug)]
pub struct SearchResult {
    pub index: usize,
    pub similarity: f64,
}

pub struct FingerprintSearch {
    gpu_morgan: Option<GpuMorganFingerprint>,
    gpu_tanimoto: Option<GpuTanimoto>,
    use_gpu: bool,
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

    pub fn search(
        &self,
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
        let fp_words =
            gpu.generate_fingerprints_batch(&[mol.clone()], radius, fp_size, false, true, false)?;

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
        &self,
        query_fp: &BitVec,
        target_fps: &[BitVec],
    ) -> anyhow::Result<Vec<f64>> {
        let gpu = self.gpu_tanimoto.as_ref().unwrap();

        let query_words = Self::bitvec_to_words(query_fp);
        let target_words: Vec<u32> = target_fps
            .iter()
            .flat_map(|fp| Self::bitvec_to_words(fp))
            .collect();

        let similarities_f32 = gpu.compute_single_query(&query_words, &target_words)?;
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
        let num_words = (bv.len() + 31) / 32;
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
