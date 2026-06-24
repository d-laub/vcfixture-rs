use std::collections::BTreeSet;

use indexmap::IndexMap;

use crate::allele::{Allele, SvType};
use crate::error::BuildError;
use crate::genotype::Genotype;
use crate::model::{AltDef, ContigDef, Document, Record, SampleValues};
use crate::spec::field::{FieldDef, FieldKind};
use crate::spec::number::{Number, NumberKind};
use crate::spec::reserved::reserved;
use crate::spec::types::Type;
use crate::spec::version::VcfVersion;
use crate::value::FieldValue;

/// Allowed SVCLAIM tokens per first-level SV type.
fn svclaim_allowed(t: SvType) -> &'static [&'static str] {
    match t {
        SvType::Del | SvType::Dup => &["D", "J", "DJ"],
        SvType::Cnv => &["D"],
        SvType::Ins | SvType::Inv => &["J"],
    }
}

fn svclaim_required(t: SvType) -> bool {
    matches!(t, SvType::Del | SvType::Dup)
}

fn cn_svlen_type(t: SvType) -> bool {
    matches!(t, SvType::Cnv | SvType::Del | SvType::Dup)
}

/// A record's spec, before validation/appending.
#[derive(Debug, Clone, Default)]
pub struct RecordSpec {
    chrom: String,
    pos: u64,
    ref_: String,
    alts: Vec<Allele>,
    ids: Option<Vec<String>>,
    qual: Option<f64>,
    filters: Option<Vec<String>>,
    gt: Option<Vec<String>>,
    info: IndexMap<String, FieldValue>,
    fmt: IndexMap<String, Vec<FieldValue>>,
    labels: BTreeSet<String>,
}

impl RecordSpec {
    pub fn at(chrom: impl Into<String>, pos: u64) -> RecordSpec {
        RecordSpec {
            chrom: chrom.into(),
            pos,
            ..Default::default()
        }
    }
    pub fn ref_(mut self, r: impl Into<String>) -> Self {
        self.ref_ = r.into();
        self
    }
    pub fn alt(mut self, alts: impl IntoIterator<Item = Allele>) -> Self {
        self.alts = alts.into_iter().collect();
        self
    }
    pub fn ids(mut self, ids: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.ids = Some(ids.into_iter().map(Into::into).collect());
        self
    }
    pub fn qual(mut self, q: f64) -> Self {
        self.qual = Some(q);
        self
    }
    pub fn filter(mut self, f: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.filters = Some(f.into_iter().map(Into::into).collect());
        self
    }
    pub fn gt(mut self, gts: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.gt = Some(gts.into_iter().map(Into::into).collect());
        self
    }
    pub fn info(mut self, id: impl Into<String>, value: FieldValue) -> Self {
        self.info.insert(id.into(), value);
        self
    }
    pub fn format(
        mut self,
        id: impl Into<String>,
        per_sample: impl IntoIterator<Item = FieldValue>,
    ) -> Self {
        self.fmt.insert(id.into(), per_sample.into_iter().collect());
        self
    }
    pub fn labels(mut self, labels: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.labels = labels.into_iter().map(Into::into).collect();
        self
    }
}

pub struct VcfBuilder {
    samples: Vec<String>,
    contigs: Vec<ContigDef>,
    version: VcfVersion,
    info_defs: IndexMap<String, FieldDef>,
    format_defs: IndexMap<String, FieldDef>,
    filter_defs: Vec<(String, String)>,
    alt_defs: IndexMap<String, String>,
    records: Vec<Record>,
}

impl VcfBuilder {
    pub fn new(
        samples: impl IntoIterator<Item = impl Into<String>>,
        contigs: impl IntoIterator<Item = (impl Into<String>, Option<u64>)>,
        version: VcfVersion,
    ) -> VcfBuilder {
        VcfBuilder {
            samples: samples.into_iter().map(Into::into).collect(),
            contigs: contigs
                .into_iter()
                .map(|(id, length)| ContigDef {
                    id: id.into(),
                    length,
                })
                .collect(),
            version,
            info_defs: IndexMap::new(),
            format_defs: IndexMap::new(),
            filter_defs: Vec::new(),
            alt_defs: IndexMap::new(),
            records: Vec::new(),
        }
    }

    fn make_def(
        &self,
        id: &str,
        number: Option<Number>,
        type_: Option<Type>,
        description: Option<String>,
        kind: FieldKind,
    ) -> Result<FieldDef, BuildError> {
        match (number, type_) {
            (Some(n), Some(t)) => FieldDef::new(
                id,
                n,
                t,
                description.unwrap_or_else(|| id.to_string()),
                kind,
            ),
            _ => reserved(id, kind, self.version),
        }
    }

    pub fn info(
        mut self,
        id: impl AsRef<str>,
        number: Option<Number>,
        type_: Option<Type>,
        description: Option<String>,
    ) -> Result<VcfBuilder, BuildError> {
        let id = id.as_ref();
        let def = self.make_def(id, number, type_, description, FieldKind::Info)?;
        self.info_defs.insert(id.to_string(), def);
        Ok(self)
    }

    pub fn format(
        mut self,
        id: impl AsRef<str>,
        number: Option<Number>,
        type_: Option<Type>,
        description: Option<String>,
    ) -> Result<VcfBuilder, BuildError> {
        let id = id.as_ref();
        let def = self.make_def(id, number, type_, description, FieldKind::Format)?;
        self.format_defs.insert(id.to_string(), def);
        Ok(self)
    }

    pub fn filter(mut self, id: impl Into<String>, description: impl Into<String>) -> VcfBuilder {
        self.filter_defs.push((id.into(), description.into()));
        self
    }

    pub fn alt(mut self, id: impl Into<String>, description: impl Into<String>) -> VcfBuilder {
        self.alt_defs.insert(id.into(), description.into());
        self
    }

    pub fn record(mut self, spec: RecordSpec) -> Result<VcfBuilder, BuildError> {
        let n_alt = spec.alts.len();
        self.validate_alleles(&spec)?;

        let mut fmt_keys: Vec<String> = Vec::new();
        let mut samples: Vec<SampleValues> = vec![SampleValues::default(); self.samples.len()];

        // GT
        if let Some(gts) = &spec.gt {
            if !self.format_defs.contains_key("GT") {
                return Err(BuildError::GtNotDeclared);
            }
            fmt_keys.push("GT".to_string());
            for (si, s) in gts.iter().enumerate() {
                let geno = Genotype::parse(s);
                for a in geno.alleles.iter().flatten() {
                    if *a as usize > n_alt {
                        return Err(BuildError::AlleleIndexOutOfRange { index: *a, n_alt });
                    }
                }
                samples[si].gt = Some(geno);
            }
        }

        let ploidy = samples
            .iter()
            .filter_map(|s| s.gt.as_ref().map(|g| g.ploidy()))
            .max()
            .unwrap_or(2);

        // FORMAT (non-GT)
        for (key, per_sample) in &spec.fmt {
            let fdef = self
                .format_defs
                .get(key)
                .ok_or_else(|| BuildError::UndeclaredField {
                    kind: "FORMAT".into(),
                    id: key.clone(),
                })?;
            fmt_keys.push(key.clone());
            let card = fdef.number.cardinality(n_alt, ploidy);
            for (si, val) in per_sample.iter().enumerate() {
                check_cardinality(key, fdef.number.kind, card, val)?;
                samples[si].values.insert(key.clone(), val.clone());
            }
        }

        // INFO
        let mut info: IndexMap<String, FieldValue> = IndexMap::new();
        for (key, val) in &spec.info {
            let fdef = self
                .info_defs
                .get(key)
                .ok_or_else(|| BuildError::UndeclaredField {
                    kind: "INFO".into(),
                    id: key.clone(),
                })?;
            let card = fdef.number.cardinality(n_alt, ploidy);
            if fdef.number.kind != NumberKind::Flag {
                check_cardinality(key, fdef.number.kind, card, val)?;
            }
            info.insert(key.clone(), val.clone());
        }

        // FORMAT CN requires equal SVLEN across CNV/DEL/DUP alleles.
        if fmt_keys.iter().any(|k| k == "CN") {
            let svlen = spec.info.get("SVLEN");
            let mut seen: Vec<Option<i64>> = Vec::new();
            for (i, a) in spec.alts.iter().enumerate() {
                if let Allele::Symbolic { first_type, .. } = a {
                    if cn_svlen_type(*first_type) {
                        seen.push(per_allele_int(svlen, i));
                    }
                }
            }
            seen.dedup();
            if seen.len() > 1 {
                return Err(BuildError::CnSvlenMismatch);
            }
        }

        self.records.push(Record {
            chrom: spec.chrom,
            pos: spec.pos,
            ids: spec.ids,
            ref_: spec.ref_,
            alts: spec.alts,
            qual: spec.qual,
            filters: spec.filters,
            info,
            fmt_keys,
            samples,
            labels: spec.labels,
        });
        Ok(self)
    }

    fn validate_alleles(&self, spec: &RecordSpec) -> Result<(), BuildError> {
        let svlen = spec.info.get("SVLEN");
        let svclaim = spec.info.get("SVCLAIM");
        let needs_padding = spec
            .alts
            .iter()
            .any(|a| matches!(a, Allele::Symbolic { .. } | Allele::Breakend { .. }));
        if needs_padding && spec.ref_.len() != 1 {
            return Err(BuildError::MissingRefPadding(spec.ref_.clone()));
        }
        for (i, a) in spec.alts.iter().enumerate() {
            let sv = per_allele_int(svlen, i);
            let cl = per_allele_str(svclaim, i);
            match a {
                Allele::Symbolic { first_type, .. } => {
                    if sv.is_none() {
                        return Err(BuildError::MissingSvlen(a.render()));
                    }
                    let allowed = svclaim_allowed(*first_type);
                    if let Some(c) = &cl {
                        if !allowed.contains(&c.as_str()) {
                            return Err(BuildError::BadSvclaim {
                                claim: c.clone(),
                                allele: a.render(),
                                allowed: allowed.iter().map(|s| s.to_string()).collect(),
                            });
                        }
                    }
                    if self.version >= VcfVersion::V4_4
                        && svclaim_required(*first_type)
                        && cl.is_none()
                    {
                        return Err(BuildError::SvclaimRequired(a.render()));
                    }
                }
                Allele::Breakend { .. } | Allele::Unspecified | Allele::Star => {
                    if sv.is_some() {
                        return Err(BuildError::SvlenMustBeMissing(a.render()));
                    }
                }
                Allele::Seq(_) => {}
            }
        }
        Ok(())
    }

    pub fn build(self) -> Result<Document, BuildError> {
        // Auto-describe symbolic ALT types; explicit .alt() descriptions win.
        let mut alt_ids: IndexMap<String, String> = IndexMap::new();
        for rec in &self.records {
            for a in &rec.alts {
                if let Some(ts) = a.symbolic_type_str() {
                    alt_ids
                        .entry(ts.clone())
                        .or_insert_with(|| format!("{ts} structural variant"));
                }
            }
        }
        for (id, desc) in &self.alt_defs {
            alt_ids.insert(id.clone(), desc.clone());
        }
        let alt_defs = alt_ids
            .into_iter()
            .map(|(id, description)| AltDef { id, description })
            .collect();

        Ok(Document {
            version: self.version,
            info_defs: self.info_defs.into_values().collect(),
            format_defs: self.format_defs.into_values().collect(),
            filter_defs: self.filter_defs,
            contigs: self.contigs,
            samples: self.samples,
            records: self.records,
            alt_defs,
        })
    }
}

/// Resolve the i-th per-allele integer of a Number=A info value.
fn per_allele_int(value: Option<&FieldValue>, i: usize) -> Option<i64> {
    match value {
        Some(FieldValue::List(v)) => v.get(i).and_then(|x| match x {
            Some(crate::value::Scalar::Int(n)) => Some(*n),
            Some(crate::value::Scalar::Float(f)) => Some(*f as i64),
            _ => None,
        }),
        Some(FieldValue::Scalar(crate::value::Scalar::Int(n))) if i == 0 => Some(*n),
        _ => None,
    }
}

fn per_allele_str(value: Option<&FieldValue>, i: usize) -> Option<String> {
    match value {
        Some(FieldValue::List(v)) => v.get(i).and_then(|x| match x {
            Some(crate::value::Scalar::Str(s)) => Some(s.clone()),
            _ => None,
        }),
        Some(FieldValue::Scalar(crate::value::Scalar::Str(s))) if i == 0 => Some(s.clone()),
        _ => None,
    }
}

/// Enforce that a list value's length equals the resolved cardinality.
/// Flag and unbounded (`None`) cardinalities are not checked. Lone scalars and
/// Flag values bypass the length check (they have no list length).
fn check_cardinality(
    id: &str,
    kind: NumberKind,
    card: Option<usize>,
    val: &FieldValue,
) -> Result<(), BuildError> {
    if kind == NumberKind::Flag {
        return Ok(());
    }
    if let (Some(expected), Some(got)) = (card, val.list_len()) {
        if expected != got {
            return Err(BuildError::Cardinality {
                id: id.to_string(),
                expected,
                got,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;
    use crate::spec::number::Number;
    use crate::spec::types::Type;
    use crate::spec::version::LATEST;
    use crate::value::FieldValue;

    fn base() -> VcfBuilder {
        VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
    }

    #[test]
    fn happy_path_builds() {
        let doc = base()
            .info("AF", None, None, None)
            .unwrap()
            .format("GT", None, None, None)
            .unwrap()
            .format("DS", Some(Number::A), Some(Type::Float), None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1000)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1"])
                    .info("AF", FieldValue::floats([0.25]))
                    .format("DS", [FieldValue::floats([0.4]), FieldValue::floats([1.9])]),
            )
            .unwrap()
            .build()
            .unwrap();
        assert_eq!(doc.records.len(), 1);
        assert_eq!(
            doc.records[0].samples[0].gt.as_ref().unwrap().render(),
            "0|1"
        );
    }

    #[test]
    fn undeclared_field_errs() {
        let r = base().format("GT", None, None, None).unwrap().record(
            RecordSpec::at("chr1", 1)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .info("AF", FieldValue::floats([0.1])),
        );
        assert!(matches!(
            r,
            Err(crate::error::BuildError::UndeclaredField { .. })
        ));
    }

    #[test]
    fn cardinality_checked() {
        let r = base()
            .format("GT", None, None, None)
            .unwrap()
            .info("AF", None, None, None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()]) // n_alt = 1, AF is Number::A
                    .info("AF", FieldValue::floats([0.1, 0.2])),
            ); // 2 values -> mismatch
        assert!(matches!(
            r,
            Err(crate::error::BuildError::Cardinality { .. })
        ));
    }

    #[test]
    fn symbolic_requires_svlen_and_padding() {
        // missing SVLEN
        let r = base()
            .format("GT", None, None, None)
            .unwrap()
            .info("SVLEN", None, None, None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())]),
            );
        assert!(matches!(r, Err(crate::error::BuildError::MissingSvlen(_))));

        // multi-base REF padding violation
        let r = base()
            .format("GT", None, None, None)
            .unwrap()
            .info("SVLEN", None, None, None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("AC")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .info("SVLEN", FieldValue::ints([100])),
            );
        assert!(matches!(
            r,
            Err(crate::error::BuildError::MissingRefPadding(_))
        ));
    }

    #[test]
    fn gt_index_out_of_range() {
        let r = base().format("GT", None, None, None).unwrap().record(
            RecordSpec::at("chr1", 1)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|2", "0|0"]),
        ); // index 2 > n_alt 1
        assert!(matches!(
            r,
            Err(crate::error::BuildError::AlleleIndexOutOfRange { .. })
        ));
    }

    #[test]
    fn gt_not_declared_errs() {
        // .gt(...) used but the builder never declared FORMAT "GT".
        let r = base().record(
            RecordSpec::at("chr1", 1)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|1", "1|1"]),
        );
        assert!(matches!(r, Err(crate::error::BuildError::GtNotDeclared)));
    }

    #[test]
    fn svlen_must_be_missing_for_unspecified() {
        // <*> (Unspecified) ALT with single-base REF padding but SVLEN set.
        let r = base()
            .format("GT", None, None, None)
            .unwrap()
            .info("SVLEN", None, None, None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::unspecified()])
                    .info("SVLEN", FieldValue::ints([100])),
            );
        assert!(matches!(
            r,
            Err(crate::error::BuildError::SvlenMustBeMissing(_))
        ));
    }

    #[test]
    fn bad_svclaim_errs() {
        // DEL allows SVCLAIM D/J/DJ; "Z" is invalid.
        let r = base()
            .format("GT", None, None, None)
            .unwrap()
            .info("SVLEN", None, None, None)
            .unwrap()
            .info("SVCLAIM", None, None, None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .info("SVLEN", FieldValue::ints([100]))
                    .info("SVCLAIM", FieldValue::strings(["Z"])),
            );
        assert!(matches!(
            r,
            Err(crate::error::BuildError::BadSvclaim { .. })
        ));
    }

    #[test]
    fn svclaim_required_for_del_at_4_5() {
        // At LATEST (V4_5 >= 4.4), DEL requires SVCLAIM; SVLEN present but no SVCLAIM.
        let r = base()
            .format("GT", None, None, None)
            .unwrap()
            .info("SVLEN", None, None, None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .info("SVLEN", FieldValue::ints([100])),
            );
        assert!(matches!(
            r,
            Err(crate::error::BuildError::SvclaimRequired(_))
        ));
    }

    #[test]
    fn cn_svlen_mismatch_errs() {
        // Two CNV alleles (CN-relevant, no SVCLAIM-required rule) with differing
        // SVLEN, plus a per-sample CN FORMAT value -> CnSvlenMismatch.
        let r = base()
            .format("GT", None, None, None)
            .unwrap()
            .format("CN", None, None, None)
            .unwrap()
            .info("SVLEN", None, None, None)
            .unwrap()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([
                        Allele::cnv(Vec::<&str>::new()),
                        Allele::cnv(Vec::<&str>::new()),
                    ])
                    .info("SVLEN", FieldValue::ints([100, 200]))
                    .format("CN", [FieldValue::floats([2.0]), FieldValue::floats([3.0])]),
            );
        assert!(matches!(r, Err(crate::error::BuildError::CnSvlenMismatch)));
    }
}
