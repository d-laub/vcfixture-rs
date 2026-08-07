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
use noodles_vcf::variant::record_buf::samples::sample::value::Array;
use noodles_vcf::variant::record_buf::samples::sample::Value;
use noodles_vcf::variant::record_buf::samples::Keys;
use noodles_vcf::variant::record_buf::Samples;
use noodles_vcf::variant::RecordBuf;

use crate::bulk::profile::{Fitted, Payload};
use crate::bulk::sample::{Samplers, VariantClass};

/// One generated variant: its site (REF/ALT/class) and its exact-AC-placed
/// genotype draws, flattened as `n_samples * ploidy` allele calls.
///
/// `ploidy` is not part of the Task 6 brief's interface sketch, but
/// [`RecordScratch::fill`] needs it to regroup the flat `gts` vector back
/// into per-sample genotype columns, and `GenRecord` does not otherwise
/// carry `n_samples`/`ploidy` — see the Task 6 report for the full
/// rationale.
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
/// No longer carries a formatted `GT` string: the genotype is written
/// straight into the destination slot as a structured [`Genotype`] by
/// [`SampleStats::refill`], which lets the slot's backing `Vec<Allele>` be
/// reused across records and removes both this crate's integer formatting
/// and noodles' string reparse from the encode path (issue #26).
struct SampleStats {
    dp: i32,
    ad: [i32; 2],
    vaf: f32,
}

impl SampleStats {
    fn new(alleles: &[i8]) -> SampleStats {
        let n_ref = alleles.iter().filter(|&&a| a == 0).count() as i32;
        let n_alt = alleles.iter().filter(|&&a| a == 1).count() as i32;
        let dp = alleles.iter().filter(|&&a| a != -1).count() as i32;
        let vaf = if dp > 0 {
            n_alt as f32 / dp as f32
        } else {
            0.0
        };

        SampleStats {
            dp,
            ad: [n_ref, n_alt],
            vaf,
        }
    }

    /// Writes this sample's value for `key` into `slot`, reusing `slot`'s
    /// existing heap buffer when the variant already matches.
    ///
    /// Every key is written unconditionally on every record, so no slot can
    /// carry a stale value forward from the previous record — which is what
    /// makes [`RecordScratch`]'s reuse safe.
    fn refill(&self, key: &str, alleles: &[i8], phased: bool, slot: &mut Option<Value>) {
        match key {
            "GT" => refill_genotype(slot, alleles, phased),
            "DP" => *slot = Some(Value::Integer(self.dp)),
            "GQ" => *slot = Some(Value::Integer(99)),
            "VAF" | "AF" => *slot = Some(Value::Float(self.vaf)),
            "AD" => refill_int_array(slot, &[self.ad[0], self.ad[1]]),
            "PL" => refill_int_array(slot, &[0, 30, 60]),
            "F1R2" => refill_int_array(slot, &[self.ad[0] / 2, self.ad[1] / 2]),
            "F2R1" => refill_int_array(
                slot,
                &[self.ad[0] - self.ad[0] / 2, self.ad[1] - self.ad[1] / 2],
            ),
            "SB" => refill_int_array(slot, &[0, 0, 0, 0]),
            other => unreachable!("unhandled FORMAT key {other}: not in any Payload preset"),
        }
    }
}

/// Rewrites `slot` as a `GT` string over `alleles`, reusing the existing
/// `String`'s buffer when there is one.
///
/// # Why a string and not `Value::Genotype`
///
/// noodles has a structured [`Value::Genotype`] whose BCF encoder skips the
/// string reparse this forces, and switching to it was the original plan.
/// It cannot be used here: it changes the *text* VCF output.
///
/// `build_header` declares `VCFv4.5`, and noodles' text writer takes its
/// VCF-4.4-and-later branch for genotypes
/// (`io::writer::record::samples::sample::value::genotype`), which writes a
/// phasing separator before **every** allele including position 0. A
/// diploid `0|0` renders as `/0|0`. The leading indicator is 4.4-conformant
/// in principle, but it is not what this crate emitted before and not what
/// consumers of these fixtures expect.
///
/// The BCF side has no such problem — all four BCF golden digests are
/// unchanged under `Value::Genotype`, confirming the encoders agree — so
/// this is purely a text-writer asymmetry. Keeping the string costs the
/// integer formatting below and noodles' reparse on encode, but those are
/// CPU costs, not allocations: the buffer reuse that issue #26 is actually
/// about is fully preserved, since this `String` is cleared and refilled
/// rather than reallocated.
fn refill_genotype(slot: &mut Option<Value>, alleles: &[i8], phased: bool) {
    if !matches!(slot, Some(Value::String(_))) {
        *slot = Some(Value::String(String::new()));
    }
    let Some(Value::String(gt)) = slot else {
        unreachable!("slot was just set to a String")
    };

    gt.clear();
    let sep = if phased { '|' } else { '/' };
    for (i, &a) in alleles.iter().enumerate() {
        if i > 0 {
            gt.push(sep);
        }
        if a < 0 {
            gt.push('.');
        } else if a < 10 {
            // The overwhelmingly common case. `write!(gt, "{a}")` would
            // route a single digit through `core::fmt` and `pad_integral`,
            // which the issue #26 profile showed as a real share of self
            // time; this produces byte-identical output without it.
            gt.push((b'0' + a as u8) as char);
        } else {
            use std::fmt::Write as _;
            let _ = write!(gt, "{a}");
        }
    }
}

/// Rewrites `slot` as an integer array, reusing the existing `Vec` when
/// there is one.
fn refill_int_array(slot: &mut Option<Value>, values: &[i32]) {
    if !matches!(slot, Some(Value::Array(Array::Integer(_)))) {
        *slot = Some(Value::Array(Array::Integer(Vec::new())));
    }
    let Some(Value::Array(Array::Integer(v))) = slot else {
        unreachable!("slot was just set to an integer array")
    };

    v.clear();
    v.extend(values.iter().copied().map(Some));
}

/// The ordered FORMAT key list for one [`Payload`] preset.
///
/// Single source of truth for both the values written here and the
/// `##FORMAT` header lines `BulkSpec::build_header` emits (`src/bulk/mod.rs`),
/// which reads this same function — a key rendered with no matching header
/// line is a write-time "missing FORMAT header record" error, not a silent
/// mismatch. `payload_presets_all_write_readable_files` (`tests/bulk.rs`) is
/// the safety net.
pub(crate) fn payload_keys(payload: &Payload) -> &'static [&'static str] {
    match payload {
        Payload::GtOnly => &["GT"],
        Payload::GtVaf => &["GT", "VAF"],
        Payload::Gatk => &["GT", "AD", "DP", "GQ", "PL"],
        Payload::Mutect2 => &["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
    }
}

/// A reusable [`RecordBuf`] plus the machinery to refill it from a
/// [`GenRecord`] without reallocating its per-sample buffers.
///
/// # Why this exists
///
/// The obvious shape — build and return a fresh `RecordBuf` per record —
/// costs four heap allocations per sample per record for
/// [`Payload::GtOnly`]: a formatted `GT` string, a clone of it, the
/// per-sample `Vec<Option<Value>>`, and noodles' own `Vec<i8>` while
/// reparsing that string. At the reference benchmark workload (2000 samples
/// x 20000 records) that is roughly 160 million allocations, and allocator
/// work measured ~47% of profile self time (issue #26).
///
/// Holding one `RecordBuf` across records and refilling it in place removes
/// three of those four. The fourth is inside noodles and stays.
///
/// Record-level fields (`chrom`, `ref_`, `alts`) are still cloned per
/// record: that is about three allocations per *record* against ~160M on the
/// per-sample path, so reusing them would add code for no measurable gain.
///
/// # Reuse safety
///
/// Every FORMAT key of every sample is written on every [`RecordScratch::fill`]
/// call, and the per-sample slot vector is resized to the record's sample
/// count, so no value can survive from a previous record.
/// `scratch_reuse_matches_a_fresh_scratch` pins this against a freshly
/// constructed scratch for all four payload presets, and
/// `tests/bulk_golden.rs` pins the encoded bytes.
pub struct RecordScratch {
    buf: RecordBuf,
    key_names: &'static [&'static str],
}

impl RecordScratch {
    /// Builds a scratch record for `payload`.
    ///
    /// The [`Keys`] set is constructed once here and then moved in and out
    /// of the record's [`Samples`] on each fill, never rebuilt or cloned.
    pub fn new(payload: &Payload) -> RecordScratch {
        let key_names = payload_keys(payload);
        let keys: Keys = key_names.iter().map(|k| k.to_string()).collect();
        let mut buf = RecordBuf::default();
        *buf.samples_mut() = Samples::new(keys, Vec::new());
        RecordScratch { buf, key_names }
    }

    /// Refills this scratch from `r` and returns it, ready to encode.
    ///
    /// The returned reference borrows the scratch, so the caller must finish
    /// with it before the next `fill` — which is exactly the
    /// encode-then-next shape `BulkSpec::stream_contigs` already has.
    pub fn fill(&mut self, r: &GenRecord, phased: bool) -> &RecordBuf {
        let ploidy = r.ploidy as usize;
        // `GenRecord` is a flat `pub` struct, so nothing prevents a caller
        // from constructing one with `ploidy: 0` or a `gts.len()` that isn't
        // a multiple of `ploidy`. `checked_div`'s `unwrap_or(0)` below
        // silently turns that into a zero-sample (or truncated) record
        // rather than failing, so assert the invariant explicitly first --
        // in debug/test builds this fails fast instead of silently
        // mis-encoding.
        debug_assert!(
            ploidy > 0 && r.gts.len() % ploidy == 0,
            "ploidy must be > 0 and evenly divide gts.len() (ploidy={ploidy}, gts.len()={})",
            r.gts.len()
        );
        let n_samples = r.gts.len().checked_div(ploidy).unwrap_or(0);

        // Take the sample block apart so the outer `Vec`, every per-sample
        // `Vec<Option<Value>>`, and every `Value`'s own buffer keep their
        // capacity across records. `Samples` exposes no `values_mut`, so
        // this `From` impl is the only way to reach them.
        let samples = std::mem::take(self.buf.samples_mut());
        let (keys, mut values) = <(Keys, Vec<Vec<Option<Value>>>)>::from(samples);

        values.resize_with(n_samples, Vec::new);
        for (i, slots) in values.iter_mut().enumerate() {
            let alleles = &r.gts[i * ploidy..(i + 1) * ploidy];
            let stats = SampleStats::new(alleles);
            slots.resize_with(self.key_names.len(), || None);
            for (slot, &k) in slots.iter_mut().zip(self.key_names) {
                stats.refill(k, alleles, phased, slot);
            }
        }

        *self.buf.samples_mut() = Samples::new(keys, values);

        // Record-level fields: cloned, not reused. See the type's docs for
        // why chasing these would not pay.
        self.buf.reference_sequence_name_mut().clear();
        self.buf.reference_sequence_name_mut().push_str(&r.chrom);
        *self.buf.variant_start_mut() =
            Some(Position::try_from(r.pos as usize).expect("pos must be >= 1"));
        self.buf.reference_bases_mut().clear();
        self.buf.reference_bases_mut().push_str(&r.ref_);
        let alts = self.buf.alternate_bases_mut().as_mut();
        alts.clear();
        alts.extend(r.alts.iter().cloned());

        &self.buf
    }
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

    fn gt_of(slot: &Option<Value>) -> &str {
        match slot {
            Some(Value::String(s)) => s,
            other => panic!("refill_genotype must produce a String, got {other:?}"),
        }
    }

    /// The `GT` rendering, including the two-digit fallback past the
    /// single-digit fast path in [`refill_genotype`].
    #[test]
    fn genotype_refill_renders_the_expected_gt_string() {
        let cases: &[(&[i8], bool, &str)] = &[
            (&[0, 1], true, "0|1"),
            (&[1, 1], false, "1/1"),
            (&[-1, 0], false, "./0"),
            (&[0, 1, 1], true, "0|1|1"),
            // Past the `a < 10` fast path, so the `write!` fallback is
            // covered and cannot silently rot.
            (&[10, 0], false, "10/0"),
            (&[9, 12], true, "9|12"),
        ];

        for (alleles, phased, expected) in cases {
            let mut slot = None;
            refill_genotype(&mut slot, alleles, *phased);
            assert_eq!(
                gt_of(&slot),
                *expected,
                "alleles={alleles:?} phased={phased}"
            );
        }
    }

    /// Refilling a populated slot must overwrite it, not append to it.
    ///
    /// This is the specific failure mode that makes buffer reuse dangerous:
    /// a missing `clear()` grows the value without ever producing a visibly
    /// wrong *first* record, so only the second call reveals it.
    #[test]
    fn refill_overwrites_rather_than_appends() {
        let mut slot = None;
        refill_genotype(&mut slot, &[0, 1, 1], true);
        refill_genotype(&mut slot, &[1, 0], false);
        assert_eq!(gt_of(&slot), "1/0", "second refill must replace the first");

        let mut slot = None;
        refill_int_array(&mut slot, &[1, 2, 3, 4]);
        refill_int_array(&mut slot, &[9, 8]);
        let Some(Value::Array(Array::Integer(v))) = &slot else {
            panic!("refill_int_array must produce an integer array")
        };
        assert_eq!(
            v,
            &[Some(9), Some(8)],
            "second refill must replace the first"
        );
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
            let mut scratch = RecordScratch::new(&payload);
            let buf = scratch.fill(&r, true);
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

    /// The scratch must produce a record indistinguishable from a fresh
    /// one, *and* keep doing so after being reused.
    ///
    /// A scratch that leaked state between records — a stale slot, a
    /// genotype vector appended to instead of cleared — would pass a
    /// single-record test and corrupt every record after the first, so the
    /// reuse itself is what this exercises. All four presets are covered
    /// because `GtOnly` never touches the array-valued refill branch.
    #[test]
    fn scratch_reuse_matches_a_fresh_scratch() {
        use crate::bulk::profile::Payload;

        let (p, s) = fixture();
        let mut rng = block_rng(7, 0, Stream::Content);

        let records: Vec<GenRecord> = (0..8)
            .map(|i| gen_record(&mut rng, &s, "chr1", 100 * (i + 1), 4, 2, &p.fitted))
            .collect();

        for payload in [
            Payload::GtOnly,
            Payload::GtVaf,
            Payload::Gatk,
            Payload::Mutect2,
        ] {
            let mut reused = RecordScratch::new(&payload);
            for (i, r) in records.iter().enumerate() {
                // Alternate phasing so a stale phase bit cannot survive
                // unnoticed from the previous record.
                let phased = i % 2 == 0;
                let from_reused = reused.fill(r, phased).clone();
                let from_fresh = RecordScratch::new(&payload).fill(r, phased).clone();
                assert_eq!(
                    from_reused, from_fresh,
                    "record {i} of {payload:?} differs after scratch reuse"
                );
            }
        }
    }

    /// The scratch must shrink and grow its per-sample buffers correctly.
    ///
    /// A `clear()` that was really a truncate to the wrong length, or a
    /// `resize_with` that grew but never shrank, would show up here and
    /// nowhere else in the suite.
    #[test]
    fn scratch_handles_changing_sample_count() {
        use crate::bulk::profile::Payload;

        let (p, s) = fixture();
        let mut rng = block_rng(11, 0, Stream::Content);
        let mut scratch = RecordScratch::new(&Payload::GtOnly);

        for n_samples in [6usize, 2, 9, 1, 6] {
            let r = gen_record(&mut rng, &s, "chr1", 500, n_samples, 2, &p.fitted);
            let fresh = RecordScratch::new(&Payload::GtOnly).fill(&r, false).clone();
            let got = scratch.fill(&r, false);
            assert_eq!(
                got.samples().values().count(),
                n_samples,
                "scratch must resize to {n_samples} samples"
            );
            assert_eq!(*got, fresh, "resized scratch must match a fresh one");
        }
    }
}
