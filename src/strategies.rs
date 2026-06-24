//! Proptest strategies for generating valid VCF documents.
//!
//! All strategies are correct-by-construction: every generated `Document`
//! passes `VcfBuilder` validation. The public `documents()` entry-point is the
//! main one consumed by property tests in downstream crates.

use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

use crate::allele::{Allele, SvType};
use crate::build::{RecordSpec, VcfBuilder};
use crate::model::Document;
use crate::reference::{DrawOpts, ReferenceBuilder, ReferenceSpec, VariantKlass};
use crate::spec::field::{FieldDef, FieldKind};
use crate::spec::number::Number;
use crate::spec::types::Type;
use crate::spec::version::{VcfVersion, LATEST};
use crate::truth::GroundTruth;
use crate::value::{FieldValue, Scalar};

/// All sequence-level variant classes that the reference-free `documents`
/// strategy can generate.
pub const ALL_VARIANT_CLASSES: [VariantKlass; 6] = [
    VariantKlass::Snp,
    VariantKlass::Mnp,
    VariantKlass::Ins,
    VariantKlass::Del,
    VariantKlass::Delins,
    VariantKlass::SpanningDel,
];

const BASES: [&str; 4] = ["A", "C", "G", "T"];

/// All classic Number × Type combinations as `(Number, Type, FieldKind)`.
/// `Flag` appears only for `FieldKind::Info` (VCF spec constraint).
pub fn all_number_type_combos() -> Vec<(Number, Type, FieldKind)> {
    let numbers = [
        Number::ONE,
        Number::fixed(2).unwrap(),
        Number::A,
        Number::R,
        Number::G,
        Number::DOT,
    ];
    let mut combos = Vec::new();
    for kind in [FieldKind::Info, FieldKind::Format] {
        let allowed: Vec<Type> = match kind {
            FieldKind::Info => Type::info_allowed()
                .into_iter()
                .filter(|t| *t != Type::Flag)
                .collect(),
            FieldKind::Format => Type::format_allowed().to_vec(),
        };
        for n in numbers {
            for t in &allowed {
                combos.push((n, *t, kind));
            }
        }
        if kind == FieldKind::Info {
            combos.push((Number::FLAG, Type::Flag, FieldKind::Info));
        }
    }
    combos
}

fn next_base(b: &str, off: usize) -> String {
    let i = BASES.iter().position(|&x| x == b).unwrap_or(0);
    BASES[(i + off) % 4].to_string()
}

/// Draw a `(ref_seq, alt_seq)` pair that is consistent with the given klass.
///
/// This ensures REF and ALT lengths agree with the biological variant class
/// (e.g. Del has a longer REF than ALT). Using this guarantees the builder's
/// per-class validation never fires.
pub fn ref_alt_strategy(klass: VariantKlass) -> impl Strategy<Value = (String, String)> {
    let base = prop::sample::select(BASES.to_vec());
    let base2 = prop::sample::select(BASES.to_vec());
    let ins_regex = prop::string::string_regex("[ACGT]{1,3}").unwrap();
    (base, base2, ins_regex, 0usize..3).prop_map(move |(b, b2, tail, snp_off)| match klass {
        VariantKlass::Snp => (b.to_string(), next_base(b, 1 + snp_off)),
        VariantKlass::Mnp => {
            let r = format!("{b}{b2}");
            let a = format!("{}{}", next_base(b, 1), next_base(b2, 1));
            (r, a)
        }
        VariantKlass::Ins => (b.to_string(), format!("{b}{tail}")),
        VariantKlass::Del => (format!("{b}{tail}"), b.to_string()),
        VariantKlass::Delins => {
            // REF is 2-base; ALT is the random tail (drawn as [ACGT]{1,3}).
            (format!("{b}{b2}"), tail)
        }
        VariantKlass::SpanningDel => (b.to_string(), "*".to_string()),
    })
}

/// Generate a single-sample GT string with `ploidy` slots, allele indices in
/// `0..=n_alt`, and each slot independently set to missing (`.`) with
/// probability `missing_rate`.
pub fn genotype_strategy(
    ploidy: usize,
    n_alt: usize,
    missing_rate: f64,
) -> impl Strategy<Value = String> {
    let slots = prop::collection::vec((0.0f64..1.0, 0u32..=(n_alt as u32)), ploidy);
    (slots, any::<bool>()).prop_map(move |(slots, phased)| {
        let sep = if phased { "|" } else { "/" };
        let tokens: Vec<String> = slots
            .into_iter()
            .map(|(r, idx)| {
                if r < missing_rate {
                    ".".to_string()
                } else {
                    idx.to_string()
                }
            })
            .collect();
        tokens.join(sep)
    })
}

/// Generate a scalar value for the given VCF field type.
pub fn scalar_strategy(t: Type) -> BoxedStrategy<Scalar> {
    match t {
        Type::Integer => (-1000i64..=1000).prop_map(Scalar::Int).boxed(),
        Type::Float => (-1.0e6f64..1.0e6)
            .prop_map(|f| Scalar::Float(f as f32 as f64))
            .boxed(),
        Type::Character => prop::sample::select(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
                .chars()
                .collect::<Vec<_>>(),
        )
        .prop_map(Scalar::Char)
        .boxed(),
        Type::String => prop::string::string_regex("[A-Za-z0-9]{1,6}")
            .unwrap()
            .prop_map(Scalar::Str)
            .boxed(),
        // Flag handled separately in field_value_strategy; this branch is unused.
        Type::Flag => Just(Scalar::Int(1)).boxed(),
    }
}

/// Generate a spec-valid `FieldValue` for the given field definition and record
/// dimensions. Cardinality is resolved exactly (`Number::cardinality`), and Dot
/// fields draw 1–3 values.
pub fn field_value_strategy(
    fd: &FieldDef,
    n_alt: usize,
    ploidy: usize,
) -> BoxedStrategy<FieldValue> {
    if fd.type_ == Type::Flag {
        return Just(FieldValue::Flag).boxed();
    }
    let card = fd.number.cardinality(n_alt, ploidy);
    let t = fd.type_;
    match card {
        Some(c) => prop::collection::vec(scalar_strategy(t), c)
            .prop_map(|xs| FieldValue::List(xs.into_iter().map(Some).collect()))
            .boxed(),
        None => (1usize..=3)
            .prop_flat_map(move |c| prop::collection::vec(scalar_strategy(t), c))
            .prop_map(|xs| FieldValue::List(xs.into_iter().map(Some).collect()))
            .boxed(),
    }
}

/// Options controlling the `documents` strategy.
#[derive(Debug, Clone)]
pub struct DocumentOpts {
    /// Maximum number of samples per document.
    pub max_samples: usize,
    /// Maximum number of records per document.
    pub max_records: usize,
    /// Maximum number of ALT alleles per record.
    pub max_alt: usize,
    /// VCF spec version for generated documents.
    pub version: VcfVersion,
}

impl Default for DocumentOpts {
    fn default() -> Self {
        DocumentOpts {
            max_samples: 3,
            max_records: 4,
            max_alt: 1,
            version: LATEST,
        }
    }
}

/// Resolve a single record REF and its list of ALT alleles from per-alt drawn
/// `(ref, alt)` pairs.
///
/// For `n_alt == 1` the drawn pair is used directly (REF and ALT are always
/// class-consistent). For multi-allelic records a single common REF is required
/// (VCF semantics), so REF is taken from the first allele's pair and the other
/// ALTs are rebuilt deterministically off the common REF anchor base. The
/// builder does not enforce REF/ALT length agreement for sequence alleles, so
/// the result always validates.
fn build_ref_alts(
    n_alt: usize,
    klasses: &[VariantKlass],
    ref_alts: &[(String, String)],
) -> (String, Vec<Allele>) {
    let final_ref = ref_alts[0].0.clone();
    if n_alt == 1 {
        return (
            final_ref,
            vec![Allele::parse(&ref_alts[0].1).expect("strategy generates a valid allele")],
        );
    }
    let refbase = &final_ref[..1]; // single-char anchor
    let mut alts: Vec<Allele> = Vec::with_capacity(n_alt);
    for (i, k) in klasses.iter().enumerate() {
        let alt_str = match k {
            VariantKlass::Snp => next_base(refbase, 1),
            VariantKlass::Mnp => {
                format!("{}{}", next_base(refbase, 1), next_base(refbase, 2))
            }
            VariantKlass::Ins => format!("{refbase}T"),
            VariantKlass::Del | VariantKlass::Delins => ref_alts[i].1.clone(),
            VariantKlass::SpanningDel => "*".to_string(),
        };
        alts.push(Allele::parse(&alt_str).expect("strategy generates a valid allele"));
    }
    (final_ref, alts)
}

/// A reference-free document strategy over a single synthetic `chr1`.
///
/// Every drawn `Document` is valid-by-construction: the builder's eager
/// validation never fires because:
/// - REF/ALT pairs are drawn via `ref_alt_strategy`, which guarantees
///   biological consistency for each class.
/// - SPANNING_DEL (`*`) is only placed in the last ALT slot; earlier slots are
///   downgraded to SNP.
/// - Genotype allele indices are bounded by `n_alt` via `genotype_strategy`.
pub fn documents(opts: DocumentOpts) -> impl Strategy<Value = Document> {
    (
        1usize..=opts.max_samples,
        1usize..=2usize, // ploidy
        1usize..=opts.max_records,
    )
        .prop_flat_map(move |(n_samples, ploidy, n_rec)| {
            let max_alt = opts.max_alt;
            let version = opts.version;

            // Per-record strategy: n_alt, one (ref, alt) pair per alt allele,
            // genotypes for all samples, and a position gap.
            let rec_strat = (1usize..=max_alt).prop_flat_map(move |n_alt| {
                let klasses = prop::collection::vec(
                    prop::sample::select(ALL_VARIANT_CLASSES.to_vec()),
                    n_alt,
                );
                let ref_alts = prop::collection::vec(
                    prop::sample::select(ALL_VARIANT_CLASSES.to_vec())
                        .prop_flat_map(ref_alt_strategy),
                    n_alt,
                );
                let gts = prop::collection::vec(genotype_strategy(ploidy, n_alt, 0.1), n_samples);
                let gap = 1u64..=50u64;
                (Just(n_alt), klasses, ref_alts, gts, gap)
            });

            prop::collection::vec(rec_strat, n_rec).prop_map(move |recs| {
                let samples: Vec<String> = (0..n_samples).map(|i| format!("s{i}")).collect();
                let mut b =
                    VcfBuilder::new(samples, [("chr1", Some(100_000u64))], version).format("GT");

                let mut pos = 1000u64;
                for (n_alt, mut klasses, ref_alts, gts, gap) in recs {
                    // SPANNING_DEL is only valid as the last ALT; downgrade any
                    // earlier ones to SNP.
                    let last = n_alt - 1;
                    for (j, k) in klasses.iter_mut().enumerate() {
                        if *k == VariantKlass::SpanningDel && j != last {
                            *k = VariantKlass::Snp;
                        }
                    }
                    let (final_ref, final_alts) = build_ref_alts(n_alt, &klasses, &ref_alts);

                    let spec = RecordSpec::at("chr1", pos)
                        .ref_(final_ref)
                        .alt(final_alts)
                        .gt(gts);
                    b = b.record(spec);
                    pos += gap;
                }
                b.build().expect("valid document")
            })
        })
}

/// A curated, valid subset of extra INFO/FORMAT fields exercised by
/// `documents_with_fields`. Each entry is `(id, Number, Type, FieldKind)` with
/// an explicit Number/Type so we control population. Flag is Info-only.
fn extra_field_defs() -> Vec<FieldDef> {
    vec![
        // INFO Number=A Float — one value per alt.
        FieldDef::new(
            "AF",
            Number::A,
            Type::Float,
            "Allele frequency",
            FieldKind::Info,
        )
        .unwrap(),
        // INFO Number=1 Integer.
        FieldDef::new(
            "NS",
            Number::ONE,
            Type::Integer,
            "Samples with data",
            FieldKind::Info,
        )
        .unwrap(),
        // INFO Flag.
        FieldDef::new(
            "DB",
            Number::FLAG,
            Type::Flag,
            "dbSNP membership",
            FieldKind::Info,
        )
        .unwrap(),
        // FORMAT Number=1 Integer.
        FieldDef::new(
            "GQ",
            Number::ONE,
            Type::Integer,
            "Genotype quality",
            FieldKind::Format,
        )
        .unwrap(),
        // FORMAT Number=R Integer — one value per allele (ref + alts).
        FieldDef::new(
            "AD",
            Number::R,
            Type::Integer,
            "Allelic depths",
            FieldKind::Format,
        )
        .unwrap(),
    ]
}

/// Like `documents`, but also declares a curated set of extra INFO and FORMAT
/// fields and populates each record's values via `field_value_strategy` so
/// cardinality always matches the record's `n_alt`/`ploidy`. Valid-by-
/// construction: `field_value_strategy` resolves the exact cardinality for the
/// declared Number/Type and returns `FieldValue::Flag` for Flag fields.
pub fn documents_with_fields(opts: DocumentOpts) -> impl Strategy<Value = Document> {
    (
        1usize..=opts.max_samples,
        1usize..=2usize, // ploidy
        1usize..=opts.max_records,
    )
        .prop_flat_map(move |(n_samples, ploidy, n_rec)| {
            let max_alt = opts.max_alt;
            let version = opts.version;
            let field_defs = extra_field_defs();

            // Per-record strategy: n_alt, per-alt (ref, alt) pairs, genotypes,
            // a position gap, and the populated INFO/FORMAT field values.
            let rec_strat = (1usize..=max_alt).prop_flat_map(move |n_alt| {
                let klasses = prop::collection::vec(
                    prop::sample::select(ALL_VARIANT_CLASSES.to_vec()),
                    n_alt,
                );
                let ref_alts = prop::collection::vec(
                    prop::sample::select(ALL_VARIANT_CLASSES.to_vec())
                        .prop_flat_map(ref_alt_strategy),
                    n_alt,
                );
                let gts = prop::collection::vec(genotype_strategy(ploidy, n_alt, 0.1), n_samples);
                let gap = 1u64..=50u64;

                // INFO values: one FieldValue per INFO field.
                let info_defs: Vec<FieldDef> = extra_field_defs()
                    .into_iter()
                    .filter(|fd| fd.kind == FieldKind::Info)
                    .collect();
                let info_strats: Vec<BoxedStrategy<FieldValue>> = info_defs
                    .iter()
                    .map(|fd| field_value_strategy(fd, n_alt, ploidy))
                    .collect();

                // FORMAT values: per field, one FieldValue per sample.
                let fmt_defs: Vec<FieldDef> = extra_field_defs()
                    .into_iter()
                    .filter(|fd| fd.kind == FieldKind::Format)
                    .collect();
                let fmt_strats: Vec<BoxedStrategy<Vec<FieldValue>>> = fmt_defs
                    .iter()
                    .map(|fd| {
                        prop::collection::vec(field_value_strategy(fd, n_alt, ploidy), n_samples)
                            .boxed()
                    })
                    .collect();

                (
                    Just(n_alt),
                    klasses,
                    ref_alts,
                    gts,
                    gap,
                    info_strats,
                    fmt_strats,
                )
            });

            prop::collection::vec(rec_strat, n_rec).prop_map(move |recs| {
                let samples: Vec<String> = (0..n_samples).map(|i| format!("s{i}")).collect();
                let mut b =
                    VcfBuilder::new(samples, [("chr1", Some(100_000u64))], version).format("GT");
                // Declare the curated extra fields.
                for fd in &field_defs {
                    b = match fd.kind {
                        FieldKind::Info => b.info(
                            crate::build::Field::typed(&fd.id, fd.number, fd.type_)
                                .description(fd.description.clone()),
                        ),
                        FieldKind::Format => b.format(
                            crate::build::Field::typed(&fd.id, fd.number, fd.type_)
                                .description(fd.description.clone()),
                        ),
                    };
                }

                let info_ids: Vec<String> = field_defs
                    .iter()
                    .filter(|fd| fd.kind == FieldKind::Info)
                    .map(|fd| fd.id.clone())
                    .collect();
                let fmt_ids: Vec<String> = field_defs
                    .iter()
                    .filter(|fd| fd.kind == FieldKind::Format)
                    .map(|fd| fd.id.clone())
                    .collect();

                let mut pos = 1000u64;
                for (n_alt, mut klasses, ref_alts, gts, gap, info_vals, fmt_vals) in recs {
                    let last = n_alt - 1;
                    for (j, k) in klasses.iter_mut().enumerate() {
                        if *k == VariantKlass::SpanningDel && j != last {
                            *k = VariantKlass::Snp;
                        }
                    }
                    let (final_ref, final_alts) = build_ref_alts(n_alt, &klasses, &ref_alts);

                    let mut spec = RecordSpec::at("chr1", pos)
                        .ref_(final_ref)
                        .alt(final_alts)
                        .gt(gts);
                    for (id, val) in info_ids.iter().zip(info_vals) {
                        spec = spec.info(id.clone(), val);
                    }
                    for (id, per_sample) in fmt_ids.iter().zip(fmt_vals) {
                        spec = spec.format(id.clone(), per_sample);
                    }
                    b = b.record(spec);
                    pos += gap;
                }
                b.build().expect("valid document")
            })
        })
}

// ── Symbolic SV documents ────────────────────────────────────────────────────

/// SV types that require SVCLAIM at VCF >= 4.4.
fn svclaim_required(t: SvType) -> bool {
    matches!(t, SvType::Del | SvType::Dup)
}

/// First allowed SVCLAIM for the type (used to make records valid).
fn svclaim_for(t: SvType) -> &'static str {
    match t {
        SvType::Del | SvType::Dup => "D",
        SvType::Ins | SvType::Inv => "J",
        SvType::Cnv => "D",
    }
}

fn sv_allele(t: SvType) -> Allele {
    match t {
        SvType::Del => Allele::deletion(Vec::<&str>::new()),
        SvType::Ins => Allele::insertion(Vec::<&str>::new()),
        SvType::Dup => Allele::duplication(Vec::<&str>::new()),
        SvType::Inv => Allele::inversion(Vec::<&str>::new()),
        SvType::Cnv => Allele::cnv(Vec::<&str>::new()),
    }
}

const ALL_SV_TYPES: [SvType; 5] = [
    SvType::Del,
    SvType::Ins,
    SvType::Dup,
    SvType::Inv,
    SvType::Cnv,
];

/// Generates documents whose records each contain exactly one symbolic SV
/// allele (DEL/INS/DUP/INV/CNV). All generated documents are valid:
/// - Single-base REF padding is enforced.
/// - SVLEN (Number=A) is declared and populated.
/// - SVCLAIM is declared and populated for types that require it at >= 4.4.
pub fn symbolic_documents(opts: DocumentOpts) -> impl Strategy<Value = Document> {
    let version = opts.version;
    (1usize..=opts.max_samples, 1usize..=opts.max_records).prop_flat_map(
        move |(n_samples, n_rec)| {
            let base_strings: Vec<String> = BASES.iter().map(|s| s.to_string()).collect();
            let rec_strat = (
                prop::sample::select(base_strings),
                prop::sample::select(ALL_SV_TYPES.to_vec()),
                100i64..=10_000i64, // SVLEN
                genotype_strategy(2, 1, 0.0),
            );
            prop::collection::vec(rec_strat, n_rec).prop_map(
                move |recs: Vec<(String, SvType, i64, String)>| {
                    let samples: Vec<String> = (0..n_samples).map(|i| format!("s{i}")).collect();

                    let mut b =
                        VcfBuilder::new(samples.clone(), [("chr1", Some(100_000u64))], version)
                            .format("GT")
                            .info("SVLEN")
                            .info("SVCLAIM");

                    let mut pos = 1000u64;
                    for (refbase, sv_type, svlen, gt) in &recs {
                        let allele = sv_allele(*sv_type);
                        let gt_vals: Vec<String> = samples.iter().map(|_| gt.clone()).collect();
                        let need_svclaim = version >= crate::spec::version::VcfVersion::V4_4
                            && svclaim_required(*sv_type);
                        let svclaim_val = svclaim_for(*sv_type);

                        let mut spec = RecordSpec::at("chr1", pos)
                            .ref_(refbase.clone())
                            .alt([allele])
                            .gt(gt_vals)
                            .info("SVLEN", FieldValue::ints([*svlen]));
                        if need_svclaim {
                            spec = spec.info("SVCLAIM", FieldValue::strings([svclaim_val]));
                        }
                        b = b.record(spec);
                        pos += 100;
                    }
                    b.build().expect("valid symbolic document")
                },
            )
        },
    )
}

// ── Reference strategies ─────────────────────────────────────────────────────

/// Generates a small `ReferenceSpec` with 1–3 contigs and optional tandem
/// repeats.
pub fn references(n_contigs: usize, contig_len: usize) -> impl Strategy<Value = ReferenceSpec> {
    let seeds = prop::collection::vec(any::<u64>(), 1..=n_contigs);
    let repeat_count = 0usize..=3;
    (seeds, repeat_count).prop_map(move |(seeds, n_repeats)| {
        let mut rb = ReferenceBuilder::new(seeds[0]);
        for (i, _seed) in seeds.iter().enumerate() {
            let name = format!("chr{}", i + 1);
            rb.add_contig(&name, contig_len).expect("add_contig");
        }
        // Plant a few tandem repeats on chr1.
        let motifs = ["CAG", "AC", "T"];
        for j in 0..n_repeats.min(3) {
            let pos0 = (j * 20).min(contig_len.saturating_sub(20));
            let motif = motifs[j % motifs.len()];
            let _ = rb.tandem_repeat("chr1", pos0, motif, 3);
        }
        rb.build()
    })
}

const REF_CONTIG_LEN: usize = 500;
/// Margin left at the end of each contig so `draw_ref_alt` (which may read up to
/// `pos0 + 2` for Del/Mnp under `DrawOpts::default`) never runs out of bounds.
const REF_DRAW_MARGIN: usize = 8;

/// Generates a `(ReferenceSpec, Document, GroundTruth)` triple whose records are
/// reference-consistent: each record draws a 0-based position within a contig
/// (with margin) and derives its `(REF, ALT)` pair from the reference sequence
/// via `ReferenceSpec::draw_ref_alt`. Positions are emitted as 1-based VCF POS.
///
/// Valid-by-construction: positions are bounded so `draw_ref_alt` stays in
/// range, and the derived REF/ALT pairs are class-consistent.
pub fn reference_and_documents(
    opts: DocumentOpts,
) -> impl Strategy<Value = (ReferenceSpec, Document, GroundTruth)> {
    let version = opts.version;
    let ref_strat = references(2, REF_CONTIG_LEN);
    ref_strat
        .prop_flat_map(move |ref_spec| {
            let n_contigs = ref_spec.contigs.len();
            // Per-record draw: contig index, 0-based position, class.
            let max_pos = REF_CONTIG_LEN.saturating_sub(REF_DRAW_MARGIN).max(1);
            let rec_strat = (
                0usize..n_contigs,
                0usize..max_pos,
                prop::sample::select(ALL_VARIANT_CLASSES.to_vec()),
            );
            (
                Just(ref_spec),
                1usize..=opts.max_samples,
                1usize..=2usize, // ploidy
                prop::collection::vec(rec_strat, 1..=opts.max_records),
            )
        })
        .prop_map(move |(ref_spec, n_samples, ploidy, recs)| {
            let samples: Vec<String> = (0..n_samples).map(|i| format!("s{i}")).collect();
            // Declare contigs matching the reference spec so records on any
            // contig validate against a declared contig.
            let contig_pairs: Vec<(String, Option<u64>)> = ref_spec
                .contigs
                .iter()
                .map(|(id, seq)| (id.clone(), Some(seq.len() as u64)))
                .collect();
            let mut b = VcfBuilder::new(samples, contig_pairs, version).format("GT");

            for (ci, pos0, klass) in recs {
                let (contig_id, _seq) = &ref_spec.contigs[ci];
                // draw_ref_alt is bounded by REF_DRAW_MARGIN; it can only fail if
                // pos0 is too close to the end, which we have excluded.
                let (ref_seq, alt_seqs) = ref_spec
                    .draw_ref_alt(contig_id, pos0, klass, &DrawOpts::default())
                    .expect("reference draw in bounds");
                let alts: Vec<Allele> = alt_seqs
                    .iter()
                    .map(|a| Allele::parse(a).expect("strategy generates a valid allele"))
                    .collect();
                // Simple valid phased GT: (ploidy-1) ref slots then one alt slot.
                // Allele index 1 is always <= n_alt (>= 1 here).
                let gt = {
                    let mut slots: Vec<&str> = vec!["0"; ploidy.saturating_sub(1)];
                    slots.push("1");
                    slots.join("|")
                };
                let gts: Vec<String> = (0..n_samples).map(|_| gt.clone()).collect();
                // VCF POS is 1-based.
                let pos = pos0 as u64 + 1;
                let spec = RecordSpec::at(contig_id.clone(), pos)
                    .ref_(ref_seq)
                    .alt(alts)
                    .gt(gts);
                b = b.record(spec);
            }
            let doc = b.build().expect("valid reference-consistent document");
            let truth = doc.truth();
            (ref_spec, doc, truth)
        })
}
