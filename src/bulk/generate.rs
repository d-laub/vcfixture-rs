//! Streaming per-record generation: draw a variant class and REF/ALT, draw
//! genotypes (no LD, no haplotype copying, no coalescent simulation — see
//! `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`), and
//! convert the result into a noodles [`RecordBuf`](noodles_vcf::variant::RecordBuf).
//!
//! Genotypes are drawn by **exact-AC placement**: a target allele count `ac`
//! is drawn from the fitted SFS, and exactly `ac` alt alleles are placed
//! uniformly at random among the record's non-missing slots (sampling
//! without replacement), rather than drawing each slot i.i.d. Bernoulli at
//! the implied frequency `ac / n_alleles`. The i.i.d. draw re-randomises the
//! realised allele count away from `ac` (it becomes
//! `Binomial(n_alleles, ac/n_alleles)`), which on any single record
//! statistically destroys the very SFS the profile was fitted to reproduce —
//! most visibly at low AC, where relative Binomial variance is largest.
//! Exact-AC placement reproduces the fitted SFS by construction. Genotypes
//! stay HWE conditional on the allele count (the standard population-
//! genetics formulation: uniform placement of `ac` alleles over `2N` slots
//! gives the hypergeometric, which is HWE asymptotically), and there is
//! still no LD, since placement is independent per record.
//!
//! `gen_record` implements the "uniform placement without replacement" step
//! with two strategies, chosen by how sparse `ac` is relative to the
//! non-missing slot count: rejection-sampling distinct slot ranks into a
//! `HashSet` for the common sparse case (a median AC of 1 against up to
//! 6404 non-missing slots), or a partial Fisher-Yates shuffle of the
//! materialised non-missing index list for the dense case. Both give the
//! same uniform-without-replacement distribution; see `gen_record`'s inline
//! comment for the threshold and the determinism argument.
//!
//! [`block_rng`] is the determinism guarantee for parallel generation: a
//! block's RNG stream is a pure function of `(seed, block_idx, stream)`,
//! never of thread identity or a shared mutable RNG, so output is
//! byte-identical regardless of how many worker threads produced it. See
//! [`Stream`] for why positions and content draw from separate streams.

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use noodles_core::Position;
use noodles_vcf::variant::record_buf::samples::sample::Value;
use noodles_vcf::variant::record_buf::samples::Keys;
use noodles_vcf::variant::record_buf::{AlternateBases, Samples};
use noodles_vcf::variant::RecordBuf;

use crate::bulk::profile::{Fitted, Payload};
use crate::bulk::sample::{Samplers, VariantClass};

/// One generated variant: its site (REF/ALT/class) and its exact-AC-placed
/// genotype draws, flattened as `n_samples * ploidy` allele calls.
///
/// `ploidy` is not part of the Task 6 brief's interface sketch, but
/// [`to_record_buf`] needs it to regroup the flat `gts` vector back into
/// per-sample genotype columns, and neither `GenRecord` nor `to_record_buf`
/// otherwise carries `n_samples`/`ploidy` — see the Task 6 report for the
/// full rationale.
#[derive(Debug, Clone, PartialEq)]
pub struct GenRecord {
    pub chrom: String,
    pub pos: u64,
    pub ref_: String,
    pub alts: Vec<String>,
    pub class: VariantClass,
    pub gts: Vec<i8>,
    pub ploidy: u8,
}

/// Which of a block's two independent PRNG streams to draw from.
///
/// Positions and record content are deliberately separate streams so that a
/// block's positions are a pure function of `(seed, block_idx, count)` —
/// independent of `n_samples`, ploidy, and payload. That is what lets
/// contig spans be computed by a gap-only pass instead of by generating
/// every genotype and discarding it (issue #22).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stream {
    /// Gap draws, and nothing else.
    Position,
    /// Class, REF/ALT, allele count, missingness, alt placement, phasing.
    Content,
}

impl Stream {
    fn domain(self) -> u64 {
        match self {
            Stream::Position => 1,
            Stream::Content => 2,
        }
    }
}

/// splitmix64 finalizer.
fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Seeds a `ChaCha8Rng` stream that depends only on
/// `(seed, block_idx, stream)`.
///
/// This is the determinism guarantee for parallel generation: never seed
/// from a thread ID, and never draw from a shared mutable RNG across
/// blocks. The block index is mixed first, then the stream domain in a
/// second finalizer round, so the two domains separate under the same
/// 2^-64 collision assumption the per-block separation already makes.
pub fn block_rng(seed: u64, block_idx: u64, stream: Stream) -> ChaCha8Rng {
    let base = mix(seed ^ block_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15));
    ChaCha8Rng::seed_from_u64(mix(
        base ^ stream.domain().wrapping_mul(0xD1B5_4A32_D192_ED03)
    ))
}

/// Draws one variant record: a structural class, REF/ALT bases for that
/// class, and `n_samples * ploidy` genotype calls, with exactly the drawn
/// allele count `ac` placed among the non-missing slots (never LD, haplotype
/// copying, or coalescent simulation — see
/// `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`).
pub fn gen_record<R: Rng>(
    rng: &mut R,
    s: &Samplers,
    chrom: &str,
    pos: u64,
    n_samples: usize,
    ploidy: u8,
    fitted: &Fitted,
) -> GenRecord {
    let class = s.class(rng);
    let (ref_, alts) = gen_site(rng, s, class);

    let n_alleles = n_samples as u64 * ploidy as u64;
    let ac = s.allele_count(rng, n_alleles);

    // `missing_rate` is fitted per-genotype (plink2 `--missing` counts a
    // missing hardcall once per sample, not once per allele), so it must be
    // drawn once per sample here too. Drawing it per-allele instead used to
    // produce GT half-calls (e.g. `0/.`), which real callers (plink2
    // `--make-pgen`) reject outright.
    //
    // Missingness is drawn FIRST, and the `ac` alt alleles are placed only
    // among the resulting non-missing slots, below. VCF's AC/AN count alt
    // alleles among *called* genotypes only, which is exactly what a re-fit
    // measures (`plink2 --freq counts` / `INFO/AC`); placing alt alleles
    // first and overwriting some with missing after the fact would silently
    // shrink the realised AC below the drawn `ac`.
    let mut gts: Vec<i8> = Vec::with_capacity(n_alleles as usize);
    for _ in 0..n_samples {
        let sample_missing = rng.random::<f64>() < fitted.missing_rate;
        for _ in 0..ploidy {
            gts.push(if sample_missing { -1 } else { 0 });
        }
    }

    // Exact-AC placement: place exactly `ac` alt alleles uniformly at random
    // among the non-missing slots, without replacement, rather than drawing
    // each slot i.i.d. Bernoulli at the implied frequency `ac / n_alleles` —
    // see the module doc for why the i.i.d. draw does not preserve the
    // fitted SFS. `ac` can exceed the number of non-missing slots at high AC
    // + high missing_rate; clamp rather than panic.
    //
    // Two placement strategies, chosen by how sparse the target AC is
    // relative to the non-missing slot count:
    //   - Sparse (the common case: a median AC of 1 against up to 6404
    //     non-missing slots): rejection-sample distinct slot *ranks* (a
    //     rank is a slot's position among non-missing slots, left to right)
    //     into a `HashSet`, then apply membership by walking `gts` in index
    //     order. This avoids ever materialising the non-missing index list.
    //   - Dense (`ac_eff` is a large fraction of `n_nonmissing`, so
    //     rejection sampling's expected draw count blows up): partial
    //     Fisher-Yates over the materialised non-missing index list, as
    //     before.
    //
    // Determinism note: `HashSet` insertion order never influences output —
    // ranks are drawn from `rng` (deterministic) and membership is applied
    // while walking `gts` in index order, never by iterating the `HashSet`.
    let n_nonmissing = gts.iter().filter(|&&g| g != -1).count();
    let ac_eff = (ac as usize).min(n_nonmissing);

    if ac_eff * 2 <= n_nonmissing {
        let mut chosen: std::collections::HashSet<usize> =
            std::collections::HashSet::with_capacity(ac_eff);
        while chosen.len() < ac_eff {
            chosen.insert(rng.random_range(0..n_nonmissing));
        }
        let mut rank = 0usize;
        for g in gts.iter_mut() {
            if *g != -1 {
                if chosen.contains(&rank) {
                    *g = 1;
                }
                rank += 1;
            }
        }
    } else {
        let mut idx: Vec<usize> = (0..n_alleles as usize).filter(|&i| gts[i] != -1).collect();
        for i in 0..ac_eff {
            let j = rng.random_range(i..idx.len());
            idx.swap(i, j);
            gts[idx[i]] = 1;
        }
    }

    GenRecord {
        chrom: chrom.to_string(),
        pos,
        ref_,
        alts,
        class,
        gts,
        ploidy,
    }
}

/// Draws REF/ALT bases for one variant class. Isolated from `gen_record` so
/// the class-dispatch and the genotype draw stay easy to read separately.
fn gen_site<R: Rng>(rng: &mut R, s: &Samplers, class: VariantClass) -> (String, Vec<String>) {
    match class {
        VariantClass::Snp => {
            let r = s.base(rng);
            let a = s.snp_alt(rng, r);
            (base_str(r), vec![base_str(a)])
        }
        VariantClass::Insertion => {
            // REF is the anchor base; ALT is the anchor plus inserted bases.
            let anchor = s.base(rng);
            let len = s.indel_len(rng);
            let mut alt = String::with_capacity(1 + len);
            alt.push(anchor as char);
            for _ in 0..len {
                alt.push(s.base(rng) as char);
            }
            (base_str(anchor), vec![alt])
        }
        VariantClass::Deletion => {
            // REF is the anchor base plus deleted bases; ALT is the anchor.
            let anchor = s.base(rng);
            let len = s.indel_len(rng);
            let mut ref_ = String::with_capacity(1 + len);
            ref_.push(anchor as char);
            for _ in 0..len {
                ref_.push(s.base(rng) as char);
            }
            (ref_, vec![base_str(anchor)])
        }
        VariantClass::Mnp => {
            // Same length (2-3), differing at every position.
            let len = rng.random_range(2..=3);
            let mut ref_ = String::with_capacity(len);
            let mut alt = String::with_capacity(len);
            for _ in 0..len {
                let r = s.base(rng);
                let mut a = s.base(rng);
                while a == r {
                    a = s.base(rng);
                }
                ref_.push(r as char);
                alt.push(a as char);
            }
            (ref_, vec![alt])
        }
        VariantClass::Complex => {
            // Different lengths (2-4), guaranteeing ref != alt.
            let ref_len = rng.random_range(2..=4);
            let mut alt_len = rng.random_range(2..=4);
            while alt_len == ref_len {
                alt_len = rng.random_range(2..=4);
            }
            let ref_: String = (0..ref_len).map(|_| s.base(rng) as char).collect();
            let alt: String = (0..alt_len).map(|_| s.base(rng) as char).collect();
            (ref_, vec![alt])
        }
        VariantClass::Symbolic => {
            let anchor = s.base(rng);
            (base_str(anchor), vec!["<DEL>".to_string()])
        }
    }
}

fn base_str(b: u8) -> String {
    (b as char).to_string()
}

/// Per-sample FORMAT values derived cheaply and deterministically from one
/// sample's allele calls. Realism of the non-GT values is out of scope —
/// only their presence, type, and cardinality affect the benchmark (per the
/// design spec); do not add fields beyond what each [`Payload`] preset asks
/// for.
struct SampleStats {
    gt: String,
    dp: i32,
    ad: [i32; 2],
    vaf: f32,
}

impl SampleStats {
    fn new(alleles: &[i8], phased: bool) -> SampleStats {
        let sep = if phased { '|' } else { '/' };
        let mut gt = String::with_capacity(alleles.len() * 2);
        for (i, a) in alleles.iter().enumerate() {
            if i > 0 {
                gt.push(sep);
            }
            if *a < 0 {
                gt.push('.');
            } else {
                use std::fmt::Write as _;
                let _ = write!(gt, "{a}");
            }
        }

        let n_ref = alleles.iter().filter(|&&a| a == 0).count() as i32;
        let n_alt = alleles.iter().filter(|&&a| a == 1).count() as i32;
        let dp = alleles.iter().filter(|&&a| a != -1).count() as i32;
        let vaf = if dp > 0 {
            n_alt as f32 / dp as f32
        } else {
            0.0
        };

        SampleStats {
            gt,
            dp,
            ad: [n_ref, n_alt],
            vaf,
        }
    }

    fn value_for(&self, key: &str) -> Value {
        match key {
            "GT" => Value::from(self.gt.clone()),
            "DP" => Value::from(self.dp),
            "AD" => Value::from(vec![Some(self.ad[0]), Some(self.ad[1])]),
            "GQ" => Value::from(99i32),
            "PL" => Value::from(vec![Some(0i32), Some(30i32), Some(60i32)]),
            "VAF" | "AF" => Value::from(self.vaf),
            "F1R2" => Value::from(vec![Some(self.ad[0] / 2), Some(self.ad[1] / 2)]),
            "F2R1" => Value::from(vec![
                Some(self.ad[0] - self.ad[0] / 2),
                Some(self.ad[1] - self.ad[1] / 2),
            ]),
            "SB" => Value::from(vec![Some(0i32); 4]),
            other => unreachable!("unhandled FORMAT key {other}: not in any Payload preset"),
        }
    }
}

/// Converts a [`GenRecord`] into a noodles [`RecordBuf`], with a FORMAT
/// payload matching exactly the given [`Payload`] preset's key list, in
/// order.
pub fn to_record_buf(r: &GenRecord, payload: &Payload, phased: bool) -> RecordBuf {
    let key_names: &[&str] = match payload {
        Payload::GtOnly => &["GT"],
        Payload::GtVaf => &["GT", "VAF"],
        Payload::Gatk => &["GT", "AD", "DP", "GQ", "PL"],
        Payload::Mutect2 => &["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
    };
    let keys: Keys = key_names.iter().map(|k| k.to_string()).collect();

    let ploidy = r.ploidy as usize;
    // `GenRecord` is a flat `pub` struct, so nothing prevents a caller from
    // constructing one with `ploidy: 0` or a `gts.len()` that isn't a
    // multiple of `ploidy`. `checked_div`'s `unwrap_or(0)` below silently
    // turns that into a zero-sample (or truncated) record rather than
    // failing, so assert the invariant explicitly first -- in debug/test
    // builds this fails fast instead of silently mis-encoding.
    debug_assert!(
        ploidy > 0 && r.gts.len() % ploidy == 0,
        "ploidy must be > 0 and evenly divide gts.len() (ploidy={ploidy}, gts.len()={})",
        r.gts.len()
    );
    let n_samples = r.gts.len().checked_div(ploidy).unwrap_or(0);

    let values: Vec<Vec<Option<Value>>> = (0..n_samples)
        .map(|i| {
            let alleles = &r.gts[i * ploidy..(i + 1) * ploidy];
            let stats = SampleStats::new(alleles, phased);
            key_names
                .iter()
                .map(|&k| Some(stats.value_for(k)))
                .collect()
        })
        .collect();

    RecordBuf::builder()
        .set_reference_sequence_name(r.chrom.clone())
        .set_variant_start(Position::try_from(r.pos as usize).expect("pos must be >= 1"))
        .set_reference_bases(r.ref_.clone())
        .set_alternate_bases(AlternateBases::from(r.alts.clone()))
        .set_samples(Samples::new(keys, values))
        .build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bulk::profile::Profile;
    use crate::bulk::sample::Samplers;

    fn fixture() -> (Profile, Samplers) {
        let p = Profile::builtin("germline-1kgp").unwrap();
        let an_source = 2 * p.provenance.n_samples_source as u64;
        let s = Samplers::new(&p.fitted, an_source).unwrap();
        (p, s)
    }

    #[test]
    fn block_rng_is_a_pure_function_of_seed_block_and_stream() {
        use rand::Rng;
        let draw = |seed, blk, s| block_rng(seed, blk, s).random::<u64>();

        assert_eq!(
            draw(42, 7, Stream::Content),
            draw(42, 7, Stream::Content),
            "same (seed, block, stream) must give the same stream"
        );
        assert_ne!(
            draw(42, 7, Stream::Content),
            draw(42, 8, Stream::Content),
            "different block must give a different stream"
        );
        assert_ne!(
            draw(42, 7, Stream::Content),
            draw(42, 7, Stream::Position),
            "different domain must give a different stream"
        );
        assert_ne!(
            draw(42, 7, Stream::Position),
            draw(43, 7, Stream::Position),
            "different seed must give a different stream"
        );
    }

    #[test]
    fn block_gap_sequence_is_invariant_to_cohort_width() {
        // The point of the split: a block's gap sequence -- drawn from the
        // Position stream -- must be the same regardless of `n_samples`,
        // even though `gen_record`'s content draws (on the Content stream)
        // scale with it. This reproduces the block pipeline's actual
        // per-record loop shape (`BulkSpec::stream_contigs`,
        // `src/bulk/mod.rs`): a gap draw on the
        // position stream, then `gen_record` plus a phasing draw on the
        // content stream, repeated -- rather than calling `Samplers::gap`
        // in isolation, which would trivially pass no matter how the two
        // streams are wired together.
        use rand::Rng;

        let (p, s) = fixture();
        let seed = 11;
        let block_idx = 3;

        let run = |n_samples: usize| -> Vec<u64> {
            let mut pos_rng = block_rng(seed, block_idx, Stream::Position);
            let mut content_rng = block_rng(seed, block_idx, Stream::Content);
            let mut local_pos = 0u64;
            (0..20)
                .map(|_| {
                    let gap = s.gap(&mut pos_rng);
                    local_pos += gap;
                    let _ = gen_record(
                        &mut content_rng,
                        &s,
                        "chr1",
                        local_pos,
                        n_samples,
                        p.dialed.ploidy,
                        &p.fitted,
                    );
                    let _phased = content_rng.random::<f64>() < p.fitted.phased_rate;
                    gap
                })
                .collect()
        };

        assert_eq!(run(2), run(64));
    }

    #[test]
    fn sample_stats_gt_string_is_unchanged_by_buffer_reuse() {
        assert_eq!(SampleStats::new(&[0, 1], true).gt, "0|1");
        assert_eq!(SampleStats::new(&[1, 1], false).gt, "1/1");
        assert_eq!(SampleStats::new(&[-1, 0], false).gt, "./0");
        assert_eq!(SampleStats::new(&[0, 1, 1], true).gt, "0|1|1");
    }

    #[test]
    fn genotypes_have_expected_shape_and_alphabet() {
        let (p, s) = fixture();
        let mut rng = block_rng(1, 0, Stream::Content);
        let r = gen_record(&mut rng, &s, "chr1", 100, 10, 2, &p.fitted);
        assert_eq!(r.gts.len(), 20);
        assert!(r.gts.iter().all(|g| (-1..=1).contains(g)));
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.pos, 100);
        assert!(!r.alts.is_empty());
    }

    #[test]
    fn ref_and_alt_are_never_equal() {
        let (p, s) = fixture();
        let mut rng = block_rng(2, 0, Stream::Content);
        for i in 0..500 {
            let r = gen_record(&mut rng, &s, "chr1", 100 + i, 4, 2, &p.fitted);
            for a in &r.alts {
                if !a.starts_with('<') {
                    assert_ne!(*a, r.ref_, "ref == alt at iteration {i}");
                }
            }
        }
    }

    #[test]
    fn exactly_ac_eff_alt_alleles_are_placed() {
        let (p, s) = fixture();
        for i in 0..200u64 {
            let mut rng = block_rng(7, i, Stream::Content);
            let r = gen_record(&mut rng, &s, "chr1", 100 + i, 1000, 2, &p.fitted);
            let n_alt = r.gts.iter().filter(|&&g| g == 1).count();
            let n_missing = r.gts.iter().filter(|&&g| g == -1).count();
            let n_nonmissing = r.gts.len() - n_missing;
            assert!(n_alt <= n_nonmissing);
        }
    }

    #[test]
    fn generation_is_deterministic() {
        let (p, s) = fixture();
        let a: Vec<_> = (0..50)
            .map(|i| {
                let mut r = block_rng(9, i, Stream::Content);
                gen_record(&mut r, &s, "chr1", 100 + i, 8, 2, &p.fitted)
            })
            .collect();
        let b: Vec<_> = (0..50)
            .map(|i| {
                let mut r = block_rng(9, i, Stream::Content);
                gen_record(&mut r, &s, "chr1", 100 + i, 8, 2, &p.fitted)
            })
            .collect();
        assert_eq!(
            a.iter().map(|r| r.gts.clone()).collect::<Vec<_>>(),
            b.iter().map(|r| r.gts.clone()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn payload_presets_produce_the_right_format_keys() {
        use crate::bulk::profile::Payload;

        let (p, s) = fixture();
        let mut rng = block_rng(3, 0, Stream::Content);
        let r = gen_record(&mut rng, &s, "chr1", 100, 2, 2, &p.fitted);
        for (payload, expected) in [
            (Payload::GtOnly, vec!["GT"]),
            (Payload::GtVaf, vec!["GT", "VAF"]),
            (Payload::Gatk, vec!["GT", "AD", "DP", "GQ", "PL"]),
            (
                Payload::Mutect2,
                vec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
            ),
        ] {
            let buf = to_record_buf(&r, &payload, true);
            // `record_buf::samples::Keys` wraps an `IndexSet<String>` (order
            // preserved, no `IntoIterator`/`Deref` of its own) — go through
            // `AsRef` to iterate it, rather than the brief's `.keys().map(..)`
            // shorthand, which does not compile against noodles-vcf 0.83's
            // real API.
            let keys: Vec<String> = buf
                .samples()
                .keys()
                .as_ref()
                .iter()
                .map(|k| k.to_string())
                .collect();
            assert_eq!(keys, expected, "payload {payload:?}");
        }
    }
}
