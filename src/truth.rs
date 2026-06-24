//! [`GroundTruth`] — the decoded oracle derived from a [`Document`]: position,
//! genotype, phasing, and per-allele arrays a parser test asserts against.

use std::collections::{BTreeSet, HashMap};

use ndarray::{Array1, Array2, Array3};

use crate::allele::{Allele, SvType};
use crate::model::Document;
use crate::value::{FieldValue, Scalar};
use crate::variants::{classify_seq, record_class, VariantClass};

/// Coarse classification of an ALT allele in the ground-truth oracle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleKind {
    Snp,
    Mnp,
    Ins,
    Del,
    Delins,
    SpanningDel,
    Symbolic,
    Unspecified,
    Bnd,
}

/// Per-allele ground-truth metadata derived from a [`crate::model::Document`].
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleTruth {
    /// Coarse variant class for this allele.
    pub kind: AlleleKind,
    /// `true` for sequence (`Seq`) alleles; `false` for symbolic/breakend/special.
    pub is_sequence: bool,
    /// Symbolic SV type string (e.g. `"DEL"`, `"DUP:TANDEM"`), if applicable.
    pub sv_type: Option<String>,
    /// Absolute SVLEN value, if present in INFO.
    pub svlen: Option<i64>,
    /// Computed end position (`pos + svlen`) for spanning SV types (DEL, DUP, INV, CNV).
    pub sv_end: Option<i64>,
}

/// Decoded oracle derived from a [`crate::model::Document`].
///
/// All arrays are indexed `[record, sample, ploidy]` unless noted.
#[derive(Debug, Clone)]
pub struct GroundTruth {
    /// Sample names in declaration order.
    pub samples: Vec<String>,
    /// Contig IDs in declaration order.
    pub contigs: Vec<String>,
    /// 1-based POS for each record; shape `[n_records]`.
    pub pos: Array1<i64>,
    /// REF allele string for each record.
    pub ref_: Vec<String>,
    /// Rendered ALT strings per record (`alts[record][alt_index]`).
    pub alts: Vec<Vec<String>>,
    /// Record-level variant class.
    pub variant_class: Vec<VariantClass>,
    /// Allele indices; shape `[n_records, n_samples, ploidy]`; `-1` for missing.
    pub genotypes: Array3<i32>,
    /// `true` if the genotype for that record+sample is phased; shape `[n_records, n_samples]`.
    pub phasing: Array2<bool>,
    /// INFO field values per record.
    pub info: Vec<HashMap<String, FieldValue>>,
    /// FORMAT field values per record per sample (excludes GT, which is in `genotypes`).
    pub format: Vec<Vec<HashMap<String, FieldValue>>>,
    /// Arbitrary labels attached via `RecordSpec::labels`.
    pub labels: Vec<BTreeSet<String>>,
    /// Per-allele truth for each record (`alts_truth[record][alt_index]`).
    pub alts_truth: Vec<Vec<AlleleTruth>>,
    /// Boolean mask: `true` for sequence alleles; shape `[n_records][n_alts]`.
    pub is_sequence_mask: Vec<Array1<bool>>,
}

/// Symbolic SV types that have a reference span (=> computed end). Excludes INS.
fn sv_spanning(t: SvType) -> bool {
    matches!(t, SvType::Del | SvType::Dup | SvType::Inv | SvType::Cnv)
}

fn seq_kind_to_allele_kind(c: VariantClass) -> AlleleKind {
    match c {
        VariantClass::Snp => AlleleKind::Snp,
        VariantClass::Mnp => AlleleKind::Mnp,
        VariantClass::Ins => AlleleKind::Ins,
        VariantClass::Del => AlleleKind::Del,
        VariantClass::Delins => AlleleKind::Delins,
        VariantClass::SpanningDel => AlleleKind::SpanningDel,
        _ => AlleleKind::Delins, // unreachable for sequence pairs
    }
}

/// Compute AlleleTruth for non-Seq alleles only.
/// Seq alleles are handled inline in `derive()` with the real record REF.
fn allele_truth(pos: u64, allele: &Allele, svlen_val: Option<i64>) -> AlleleTruth {
    match allele {
        Allele::Star => AlleleTruth {
            kind: AlleleKind::SpanningDel,
            is_sequence: false,
            sv_type: None,
            svlen: None,
            sv_end: None,
        },
        Allele::Unspecified => AlleleTruth {
            kind: AlleleKind::Unspecified,
            is_sequence: false,
            sv_type: None,
            svlen: None,
            sv_end: None,
        },
        Allele::Breakend { .. } => AlleleTruth {
            kind: AlleleKind::Bnd,
            is_sequence: false,
            sv_type: None,
            svlen: None,
            sv_end: None,
        },
        Allele::Symbolic { first_type, .. } => {
            let svlen = svlen_val.map(|v| v.abs());
            let end = match (svlen, sv_spanning(*first_type)) {
                (Some(l), true) => Some(pos as i64 + l),
                _ => None,
            };
            AlleleTruth {
                kind: AlleleKind::Symbolic,
                is_sequence: false,
                sv_type: allele.symbolic_type_str(),
                svlen,
                sv_end: end,
            }
        }
        // Seq is never passed here; handled inline in derive().
        Allele::Seq(_) => unreachable!("Seq allele must be handled inline in derive()"),
    }
}

/// Derive the [`GroundTruth`] oracle from a validated [`Document`].
pub fn derive(doc: &Document) -> GroundTruth {
    let n_rec = doc.records.len();
    let n_smp = doc.samples.len();
    let ploidy = doc.max_ploidy();

    // Resolution B: n_rec is usize; no .max(0) needed.
    let mut genos = Array3::<i32>::from_elem((n_rec, n_smp, ploidy), -1);
    let mut phasing = Array2::<bool>::from_elem((n_rec, n_smp), false);
    let mut pos = Array1::<i64>::zeros(n_rec);
    let mut ref_ = Vec::with_capacity(n_rec);
    let mut alts = Vec::with_capacity(n_rec);
    let mut vclass = Vec::with_capacity(n_rec);
    let mut info = Vec::with_capacity(n_rec);
    let mut fmt = Vec::with_capacity(n_rec);
    let mut labels = Vec::with_capacity(n_rec);
    let mut alts_truth = Vec::with_capacity(n_rec);
    let mut seq_mask = Vec::with_capacity(n_rec);

    for (ri, rec) in doc.records.iter().enumerate() {
        pos[ri] = rec.pos as i64;
        ref_.push(rec.ref_.clone());
        alts.push(rec.alts.iter().map(|a| a.render()).collect());
        vclass.push(record_class(&rec.ref_, &rec.alts));
        info.push(
            rec.info
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect(),
        );

        let mut per_sample = Vec::with_capacity(n_smp);
        for (si, sample) in rec.samples.iter().enumerate() {
            if let Some(gt) = &sample.gt {
                for (ai, allele) in gt.alleles.iter().enumerate() {
                    genos[[ri, si, ai]] = match allele {
                        Some(v) => *v as i32,
                        None => -1,
                    };
                }
                phasing[[ri, si]] = gt.is_phased();
            }
            // GT lives in genotypes/phasing; exclude it from the format map.
            let m: HashMap<String, FieldValue> = sample
                .values
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect();
            per_sample.push(m);
        }
        fmt.push(per_sample);
        labels.push(rec.labels.clone());

        // Per-ALT truth. SVLEN is Number=A => per-allele list of ints.
        let svlen = rec.info.get("SVLEN");
        let mut per_alt = Vec::with_capacity(rec.alts.len());
        for (ai, allele) in rec.alts.iter().enumerate() {
            let sv = svlen_at(svlen, ai);
            // Resolution A: handle Seq inline with the real record REF.
            let at = match allele {
                Allele::Seq(bases) => AlleleTruth {
                    kind: seq_kind_to_allele_kind(classify_seq(&rec.ref_, bases)),
                    is_sequence: true,
                    sv_type: None,
                    svlen: None,
                    sv_end: None,
                },
                _ => allele_truth(rec.pos, allele, sv),
            };
            per_alt.push(at);
        }
        let mask = Array1::from(per_alt.iter().map(|a| a.is_sequence).collect::<Vec<_>>());
        seq_mask.push(mask);
        alts_truth.push(per_alt);
    }

    GroundTruth {
        samples: doc.samples.clone(),
        contigs: doc.contigs.iter().map(|c| c.id.clone()).collect(),
        pos,
        ref_,
        alts,
        variant_class: vclass,
        genotypes: genos,
        phasing,
        info,
        format: fmt,
        labels,
        alts_truth,
        is_sequence_mask: seq_mask,
    }
}

fn svlen_at(value: Option<&FieldValue>, i: usize) -> Option<i64> {
    match value {
        Some(FieldValue::List(v)) => v.get(i).and_then(|x| match x {
            Some(Scalar::Int(n)) => Some(*n),
            Some(Scalar::Float(f)) => Some(*f as i64),
            _ => None,
        }),
        Some(FieldValue::Scalar(Scalar::Int(n))) if i == 0 => Some(*n),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;
    use crate::spec::version::LATEST;
    use crate::value::FieldValue;
    use crate::variants::VariantClass;
    use crate::RecordSpec;
    use crate::VcfBuilder;

    #[test]
    fn genotypes_phasing_and_missing() {
        let t = VcfBuilder::new(["s1", "s2"], [("chr1", Some(1000u64))], LATEST)
            .format("GT")
            .record(
                RecordSpec::at("chr1", 10)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "./."]),
            )
            .build()
            .unwrap()
            .truth();
        assert_eq!(t.genotypes[[0, 0, 0]], 0);
        assert_eq!(t.genotypes[[0, 0, 1]], 1);
        assert_eq!(t.genotypes[[0, 1, 0]], -1);
        assert!(t.phasing[[0, 0]]);
        assert!(!t.phasing[[0, 1]]);
        assert_eq!(t.variant_class[0], VariantClass::Snp);
    }

    #[test]
    fn symbolic_svlen_and_end() {
        let t = VcfBuilder::new(["s1"], [("chr1", Some(100_000u64))], LATEST)
            .format("GT")
            .info("SVLEN")
            .info("SVCLAIM")
            .record(
                RecordSpec::at("chr1", 100)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .info("SVLEN", FieldValue::ints([250]))
                    .info("SVCLAIM", FieldValue::strings(["D"]))
                    .gt(["0|1"]),
            )
            .build()
            .unwrap()
            .truth();
        let at = &t.alts_truth[0][0];
        assert_eq!(at.kind, AlleleKind::Symbolic);
        assert_eq!(at.svlen, Some(250));
        assert_eq!(at.sv_end, Some(350)); // pos + svlen
        assert!(!t.is_sequence_mask[0][0]);
    }

    #[test]
    fn info_excludes_gt_from_format() {
        let t = VcfBuilder::new(["s1"], [("chr1", Some(1000u64))], LATEST)
            .format("GT")
            .format("GQ")
            .record(
                RecordSpec::at("chr1", 10)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1"])
                    .format("GQ", [FieldValue::ints([42])]),
            )
            .build()
            .unwrap()
            .truth();
        assert!(t.format[0][0].contains_key("GQ"));
        assert!(!t.format[0][0].contains_key("GT"));
    }
}
