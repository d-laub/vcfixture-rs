//! Streaming per-record generation: draw a variant class and REF/ALT, draw
//! genotypes i.i.d. from HWE (no LD, no haplotype copying, no coalescent
//! simulation — see `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`),
//! and convert the result into a noodles [`RecordBuf`](noodles_vcf::variant::RecordBuf).
//!
//! [`block_rng`] is the determinism guarantee for parallel generation: a
//! block's RNG stream is a pure function of `(seed, block_idx)`, never of
//! thread identity or a shared mutable RNG, so output is byte-identical
//! regardless of how many worker threads produced it.

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use noodles_core::Position;
use noodles_vcf::variant::record_buf::samples::sample::Value;
use noodles_vcf::variant::record_buf::samples::Keys;
use noodles_vcf::variant::record_buf::{AlternateBases, Samples};
use noodles_vcf::variant::RecordBuf;

use crate::bulk::profile::{Fitted, Payload};
use crate::bulk::sample::{Samplers, VariantClass};

/// One generated variant: its site (REF/ALT/class) and its i.i.d.-HWE
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

/// Seeds a `ChaCha8Rng` stream that depends only on `(seed, block_idx)`.
///
/// This is the determinism guarantee for parallel generation: never seed
/// from a thread ID, and never draw from a shared mutable RNG across
/// blocks. A splitmix64-style finalizer keeps adjacent block indices'
/// streams well-separated despite differing from `seed` by a single
/// multiplication.
pub fn block_rng(seed: u64, block_idx: u64) -> ChaCha8Rng {
    let mut z = seed ^ block_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    ChaCha8Rng::seed_from_u64(z)
}

/// Draws one variant record: a structural class, REF/ALT bases for that
/// class, and `n_samples * ploidy` genotype calls drawn i.i.d. from HWE
/// (never LD, haplotype copying, or coalescent simulation — see
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
    let p = ac as f64 / (n_alleles.max(1) as f64);

    let gts: Vec<i8> = (0..n_alleles)
        .map(|_| {
            if rng.gen::<f64>() < fitted.missing_rate {
                -1
            } else if rng.gen::<f64>() < p {
                1
            } else {
                0
            }
        })
        .collect();

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
            let len = rng.gen_range(2..=3);
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
            let ref_len = rng.gen_range(2..=4);
            let mut alt_len = rng.gen_range(2..=4);
            while alt_len == ref_len {
                alt_len = rng.gen_range(2..=4);
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
        let gt = alleles
            .iter()
            .map(|a| {
                if *a < 0 {
                    ".".to_string()
                } else {
                    a.to_string()
                }
            })
            .collect::<Vec<_>>()
            .join(&sep.to_string());

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
pub fn to_record_buf(r: &GenRecord, payload: Payload, phased: bool) -> RecordBuf {
    let key_names: &[&str] = match payload {
        Payload::GtOnly => &["GT"],
        Payload::GtVaf => &["GT", "VAF"],
        Payload::Gatk => &["GT", "AD", "DP", "GQ", "PL"],
        Payload::Mutect2 => &["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
    };
    let keys: Keys = key_names.iter().map(|k| k.to_string()).collect();

    let ploidy = r.ploidy as usize;
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
        let s = Samplers::new(&p.fitted).unwrap();
        (p, s)
    }

    #[test]
    fn block_rng_is_a_pure_function_of_seed_and_block() {
        use rand::Rng;
        let mut a = block_rng(42, 7);
        let mut b = block_rng(42, 7);
        let mut c = block_rng(42, 8);
        let xa: u64 = a.gen();
        let xb: u64 = b.gen();
        let xc: u64 = c.gen();
        assert_eq!(xa, xb, "same (seed, block) must give the same stream");
        assert_ne!(xa, xc, "different block must give a different stream");
    }

    #[test]
    fn genotypes_have_expected_shape_and_alphabet() {
        let (p, s) = fixture();
        let mut rng = block_rng(1, 0);
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
        let mut rng = block_rng(2, 0);
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
    fn generation_is_deterministic() {
        let (p, s) = fixture();
        let a: Vec<_> = (0..50)
            .map(|i| {
                let mut r = block_rng(9, i);
                gen_record(&mut r, &s, "chr1", 100 + i, 8, 2, &p.fitted)
            })
            .collect();
        let b: Vec<_> = (0..50)
            .map(|i| {
                let mut r = block_rng(9, i);
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
        let mut rng = block_rng(3, 0);
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
            let buf = to_record_buf(&r, payload.clone(), true);
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
