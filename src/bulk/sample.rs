//! Profile-driven statistical samplers for bulk generation.
//!
//! Every sampler draws i.i.d. from the [`crate::bulk::profile::Fitted`]
//! distributions in a [`crate::bulk::profile::Profile`] — no LD, no
//! haplotype copying, no coalescent simulation. This module draws each
//! record's *target* allele count `ac` from the fitted SFS
//! ([`Samplers::allele_count`]); [`crate::bulk::gen::gen_record`] then places
//! exactly `ac` alt alleles among the record's genotypes (not an i.i.d.
//! Bernoulli draw at the implied frequency — see that function's doc for
//! why), which keeps genotypes HWE conditional on `ac`; see
//! `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`.
//!
//! [`Samplers::new`] precomputes cumulative weights once so that per-record
//! sampling ([`Samplers::gap`], [`Samplers::allele_count`], etc.) is a
//! binary search rather than a re-scan — this runs roughly 265k times per
//! bulk-generation run, so the precomputation matters.

use rand::Rng;

use crate::bulk::profile::{ClassMix, Fitted, Histogram};
use crate::bulk::BulkError;

/// The structural class of a generated variant record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass {
    Snp,
    Insertion,
    Deletion,
    Mnp,
    Complex,
    Symbolic,
}

/// A histogram sampler with a precomputed CDF over its bins.
#[derive(Debug, Clone)]
struct HistSampler {
    edges: Vec<f64>,
    cdf: Vec<f64>,
}

impl HistSampler {
    fn new(h: &Histogram) -> Result<HistSampler, BulkError> {
        h.validate()?;
        let total: f64 = h.weights.iter().sum();
        let mut cdf = Vec::with_capacity(h.weights.len());
        let mut acc = 0.0;
        for w in &h.weights {
            acc += w / total;
            cdf.push(acc);
        }
        Ok(HistSampler {
            edges: h.edges.clone(),
            cdf,
        })
    }

    /// Draw a bin by CDF binary search, then a value uniformly within it.
    ///
    /// Callers that need an integer quantize with `.floor()`, not
    /// `.round()`: flooring a value uniform on `[lo, hi)` stays within the
    /// bin's own integer range and so preserves the bin's probability mass
    /// exactly, whereas rounding straddles the bin edge and bleeds ~half of
    /// a bin's mass into its neighbor (this is most visible on the
    /// `sfs` histogram's `[1, 2)` singleton bin, which must map to exactly
    /// `1`, never `2`).
    fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        let u: f64 = rng.gen();
        let bin = self.cdf.partition_point(|c| *c < u).min(self.cdf.len() - 1);
        let (lo, hi) = (self.edges[bin], self.edges[bin + 1]);
        rng.gen_range(lo..hi)
    }
}

/// Precomputed samplers for one [`Fitted`] profile.
///
/// Construction validates and precomputes cumulative weights once; sampling
/// afterwards touches no allocation and is a binary search over a small
/// slice.
#[derive(Debug, Clone)]
pub struct Samplers {
    gap: HistSampler,
    sfs: HistSampler,
    indel: HistSampler,
    class_cdf: [f64; 6],
    ti_frac: f64,
}

impl Samplers {
    /// Build samplers from a fitted profile, precomputing CDFs once.
    pub fn new(fitted: &Fitted) -> Result<Samplers, BulkError> {
        fitted.variant_classes.validate()?;
        if fitted.titv <= 0.0 {
            return Err(BulkError::Invalid("titv must be > 0".into()));
        }
        let m: &ClassMix = &fitted.variant_classes;
        let mut acc = 0.0;
        let mut class_cdf = [0.0f64; 6];
        for (i, w) in [m.snp, m.insertion, m.deletion, m.mnp, m.complex, m.symbolic]
            .iter()
            .enumerate()
        {
            acc += w;
            class_cdf[i] = acc;
        }
        Ok(Samplers {
            gap: HistSampler::new(&fitted.gap_dist)?,
            sfs: HistSampler::new(&fitted.sfs)?,
            indel: HistSampler::new(&fitted.indel_length)?,
            class_cdf,
            ti_frac: fitted.titv / (fitted.titv + 1.0),
        })
    }

    /// Draw a gap (in bp) to the next variant. Always `>= 1`.
    pub fn gap<R: Rng>(&self, rng: &mut R) -> u64 {
        (self.gap.sample(rng).floor() as u64).max(1)
    }

    /// Draw an allele count from the fitted site-frequency spectrum,
    /// clamped to `1..=n_alleles`.
    pub fn allele_count<R: Rng>(&self, rng: &mut R, n_alleles: u64) -> u64 {
        (self.sfs.sample(rng).floor() as u64).clamp(1, n_alleles.max(1))
    }

    /// Draw an indel length in bases. Always `>= 1`.
    pub fn indel_len<R: Rng>(&self, rng: &mut R) -> usize {
        (self.indel.sample(rng).floor() as usize).max(1)
    }

    /// Draw a structural variant class from the fitted class mix.
    pub fn class<R: Rng>(&self, rng: &mut R) -> VariantClass {
        let u: f64 = rng.gen();
        let i = self.class_cdf.partition_point(|c| *c < u).min(5);
        [
            VariantClass::Snp,
            VariantClass::Insertion,
            VariantClass::Deletion,
            VariantClass::Mnp,
            VariantClass::Complex,
            VariantClass::Symbolic,
        ][i]
    }

    /// Draw a uniformly random base, one of `b"ACGT"`.
    pub fn base<R: Rng>(&self, rng: &mut R) -> u8 {
        b"ACGT"[rng.gen_range(0..4)]
    }

    /// Draw a SNP ALT base `!= ref_base`, with transitions drawn at
    /// `titv / (titv + 1)` of SNPs.
    pub fn snp_alt<R: Rng>(&self, rng: &mut R, ref_base: u8) -> u8 {
        let transition = match ref_base {
            b'A' => b'G',
            b'G' => b'A',
            b'C' => b'T',
            b'T' => b'C',
            _ => b'A',
        };
        if rng.gen::<f64>() < self.ti_frac {
            transition
        } else {
            let transversions: [u8; 2] = match ref_base {
                b'A' | b'G' => [b'C', b'T'],
                _ => [b'A', b'G'],
            };
            transversions[rng.gen_range(0..2)]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bulk::profile::Profile;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn samplers() -> Samplers {
        let p = Profile::builtin("germline-1kgp").unwrap();
        Samplers::new(&p.fitted).unwrap()
    }

    #[test]
    fn gaps_are_at_least_one() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        for _ in 0..1000 {
            assert!(s.gap(&mut rng) >= 1);
        }
    }

    #[test]
    fn allele_count_is_in_range() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        for _ in 0..1000 {
            let ac = s.allele_count(&mut rng, 6404);
            assert!((1..=6404).contains(&ac), "ac out of range: {ac}");
        }
    }

    #[test]
    fn sfs_reproduces_singleton_fraction() {
        // Targets `germline-1kgp-unphased`, not `germline-1kgp`: phasing
        // drops unphaseable singletons, so the phased panel's SFS is ~0%
        // singletons by construction (see
        // `germline_variants_differ_in_phasing` in `profile.rs`). The raw
        // unphased callset puts ~35.8% of weight in the [1, 2) bin -- this
        // is the whole point of fitting an empirical SFS (a neutral 1/x SFS
        // would give ~12%), so guard it here.
        let p = Profile::builtin("germline-1kgp-unphased").unwrap();
        let s = Samplers::new(&p.fitted).unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let n = 20_000;
        let singletons = (0..n)
            .filter(|_| s.allele_count(&mut rng, 6404) == 1)
            .count();
        let frac = singletons as f64 / n as f64;
        assert!(
            (frac - 0.358).abs() < 0.02,
            "singleton fraction {frac} != ~0.358"
        );
    }

    #[test]
    fn class_mix_is_reproduced() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(4);
        let n = 20_000;
        let snps = (0..n)
            .filter(|_| matches!(s.class(&mut rng), VariantClass::Snp))
            .count();
        let frac = snps as f64 / n as f64;
        assert!((frac - 0.87).abs() < 0.02, "snp fraction {frac} != ~0.87");
    }

    #[test]
    fn snp_alt_never_equals_ref_and_respects_titv() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let n = 20_000;
        let mut ti = 0usize;
        for _ in 0..n {
            let alt = s.snp_alt(&mut rng, b'A');
            assert_ne!(alt, b'A');
            if alt == b'G' {
                ti += 1;
            } // A<->G is a transition
        }
        // titv = 2.05 => transitions are 2.05 / 3.05 of SNPs
        let frac = ti as f64 / n as f64;
        assert!((frac - 2.05 / 3.05).abs() < 0.02, "ti fraction {frac}");
    }

    #[test]
    fn sampling_is_deterministic_for_a_seed() {
        let s = samplers();
        let mut a = ChaCha8Rng::seed_from_u64(7);
        let mut b = ChaCha8Rng::seed_from_u64(7);
        let xs: Vec<u64> = (0..100).map(|_| s.gap(&mut a)).collect();
        let ys: Vec<u64> = (0..100).map(|_| s.gap(&mut b)).collect();
        assert_eq!(xs, ys);
    }
}
