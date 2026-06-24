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

/// A field declaration, resolved to a `FieldDef` at build time.
#[derive(Debug, Clone)]
pub struct Field {
    id: String,
    decl: Decl,
    description: Option<String>,
}

#[derive(Debug, Clone)]
enum Decl {
    /// Look the field up in the reserved registry for the document version.
    Reserved,
    /// Explicit `Number` and `Type`.
    Typed(Number, Type),
    /// Flag field: `Number=0`, `Type=Flag` (INFO only; enforced at build).
    Flag,
}

impl Field {
    /// Resolve `id` via the reserved registry at build time.
    pub fn reserved(id: impl Into<String>) -> Field {
        Field {
            id: id.into(),
            decl: Decl::Reserved,
            description: None,
        }
    }

    /// Declare `id` with an explicit `number` and `type_`.
    pub fn typed(id: impl Into<String>, number: Number, type_: Type) -> Field {
        Field {
            id: id.into(),
            decl: Decl::Typed(number, type_),
            description: None,
        }
    }

    /// Declare a Flag field (`Number=0`, `Type=Flag`). Valid for INFO only.
    pub fn flag(id: impl Into<String>) -> Field {
        Field {
            id: id.into(),
            decl: Decl::Flag,
            description: None,
        }
    }

    /// Set the `Description=` header text. Defaults to the field id.
    pub fn description(mut self, d: impl Into<String>) -> Field {
        self.description = Some(d.into());
        self
    }

    /// Resolve to a concrete `FieldDef` for the given kind and version.
    fn resolve(&self, kind: FieldKind, version: VcfVersion) -> Result<FieldDef, BuildError> {
        let desc = || self.description.clone().unwrap_or_else(|| self.id.clone());
        match &self.decl {
            Decl::Reserved => reserved(&self.id, kind, version),
            Decl::Typed(number, type_) => {
                FieldDef::new(self.id.as_str(), *number, *type_, desc(), kind)
            }
            Decl::Flag => FieldDef::new(self.id.as_str(), Number::FLAG, Type::Flag, desc(), kind),
        }
    }
}

impl From<&str> for Field {
    fn from(id: &str) -> Field {
        Field::reserved(id)
    }
}

impl From<String> for Field {
    fn from(id: String) -> Field {
        Field::reserved(id)
    }
}

pub struct VcfBuilder {
    samples: Vec<String>,
    contigs: Vec<ContigDef>,
    version: VcfVersion,
    info_fields: Vec<Field>,
    format_fields: Vec<Field>,
    filter_defs: Vec<(String, String)>,
    alt_defs: IndexMap<String, String>,
    records: Vec<RecordSpec>,
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
            info_fields: Vec::new(),
            format_fields: Vec::new(),
            filter_defs: Vec::new(),
            alt_defs: IndexMap::new(),
            records: Vec::new(),
        }
    }

    pub fn info(mut self, field: impl Into<Field>) -> VcfBuilder {
        self.info_fields.push(field.into());
        self
    }

    pub fn format(mut self, field: impl Into<Field>) -> VcfBuilder {
        self.format_fields.push(field.into());
        self
    }

    pub fn filter(mut self, id: impl Into<String>, description: impl Into<String>) -> VcfBuilder {
        self.filter_defs.push((id.into(), description.into()));
        self
    }

    pub fn alt(mut self, id: impl Into<String>, description: impl Into<String>) -> VcfBuilder {
        self.alt_defs.insert(id.into(), description.into());
        self
    }

    pub fn record(mut self, spec: RecordSpec) -> VcfBuilder {
        self.records.push(spec);
        self
    }

    pub fn render(self) -> Result<String, BuildError> {
        Ok(self.build()?.render())
    }

    pub fn write(
        self,
        path: impl AsRef<std::path::Path>,
        opts: crate::write::WriteOpts,
    ) -> Result<std::path::PathBuf, BuildError> {
        self.build()?.write(path, opts)
    }

    pub fn truth(self) -> Result<crate::truth::GroundTruth, BuildError> {
        Ok(self.build()?.truth())
    }

    pub fn build(self) -> Result<Document, BuildError> {
        // 1. Resolve field declarations to concrete defs (reserved lookup,
        //    explicit FieldDef::new). Last declaration of an id wins.
        let mut info_defs: IndexMap<String, FieldDef> = IndexMap::new();
        for field in &self.info_fields {
            let def = field.resolve(FieldKind::Info, self.version)?;
            info_defs.insert(def.id.clone(), def);
        }
        let mut format_defs: IndexMap<String, FieldDef> = IndexMap::new();
        for field in &self.format_fields {
            let def = field.resolve(FieldKind::Format, self.version)?;
            format_defs.insert(def.id.clone(), def);
        }

        // 2. Validate and convert each record, tagging errors with their index.
        let mut records: Vec<Record> = Vec::with_capacity(self.records.len());
        for (index, spec) in self.records.into_iter().enumerate() {
            let rec = build_record(spec, &self.samples, self.version, &info_defs, &format_defs)
                .map_err(|source| BuildError::InRecord {
                    index,
                    source: Box::new(source),
                })?;
            records.push(rec);
        }

        // 3. Auto-describe symbolic ALT types; explicit .alt() descriptions win.
        let mut alt_ids: IndexMap<String, String> = IndexMap::new();
        for rec in &records {
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
            info_defs: info_defs.into_values().collect(),
            format_defs: format_defs.into_values().collect(),
            filter_defs: self.filter_defs,
            contigs: self.contigs,
            samples: self.samples,
            records,
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

fn validate_alleles(spec: &RecordSpec, version: VcfVersion) -> Result<(), BuildError> {
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
                if version >= VcfVersion::V4_4 && svclaim_required(*first_type) && cl.is_none() {
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

/// Validate a `RecordSpec` against the resolved field defs and convert it to a
/// `Record`. This is the per-record pipeline run by `VcfBuilder::build`.
fn build_record(
    spec: RecordSpec,
    samples: &[String],
    version: VcfVersion,
    info_defs: &IndexMap<String, FieldDef>,
    format_defs: &IndexMap<String, FieldDef>,
) -> Result<Record, BuildError> {
    let n_alt = spec.alts.len();
    validate_alleles(&spec, version)?;

    let mut fmt_keys: Vec<String> = Vec::new();
    let mut sample_vals: Vec<SampleValues> = vec![SampleValues::default(); samples.len()];

    // GT
    if let Some(gts) = &spec.gt {
        if !format_defs.contains_key("GT") {
            return Err(BuildError::GtNotDeclared);
        }
        if gts.len() != samples.len() {
            return Err(BuildError::SampleCountMismatch {
                kind: "GT".into(),
                expected: samples.len(),
                got: gts.len(),
            });
        }
        fmt_keys.push("GT".to_string());
        for (si, s) in gts.iter().enumerate() {
            let geno = Genotype::parse(s)?;
            for a in geno.alleles.iter().flatten() {
                if *a as usize > n_alt {
                    return Err(BuildError::AlleleIndexOutOfRange { index: *a, n_alt });
                }
            }
            sample_vals[si].gt = Some(geno);
        }
    }

    let ploidy = sample_vals
        .iter()
        .filter_map(|s| s.gt.as_ref().map(|g| g.ploidy()))
        .max()
        .unwrap_or(2);

    // FORMAT (non-GT)
    for (key, per_sample) in &spec.fmt {
        let fdef = format_defs
            .get(key)
            .ok_or_else(|| BuildError::UndeclaredField {
                kind: "FORMAT".into(),
                id: key.clone(),
            })?;
        if per_sample.len() != samples.len() {
            return Err(BuildError::SampleCountMismatch {
                kind: key.clone(),
                expected: samples.len(),
                got: per_sample.len(),
            });
        }
        fmt_keys.push(key.clone());
        let card = fdef.number.cardinality(n_alt, ploidy);
        for (si, val) in per_sample.iter().enumerate() {
            check_cardinality(key, fdef.number.kind, card, val)?;
            sample_vals[si].values.insert(key.clone(), val.clone());
        }
    }

    // INFO
    let mut info: IndexMap<String, FieldValue> = IndexMap::new();
    for (key, val) in &spec.info {
        let fdef = info_defs
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

    Ok(Record {
        chrom: spec.chrom,
        pos: spec.pos,
        ids: spec.ids,
        ref_: spec.ref_,
        alts: spec.alts,
        qual: spec.qual,
        filters: spec.filters,
        info,
        fmt_keys,
        samples: sample_vals,
        labels: spec.labels,
    })
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

    /// Match a per-record `BuildError` produced during `build()`.
    macro_rules! assert_in_record {
        ($result:expr, $pat:pat) => {
            match $result {
                Err(BuildError::InRecord { source, .. }) => {
                    assert!(
                        matches!(*source, $pat),
                        "unexpected inner error: {source:?}"
                    );
                }
                other => panic!("expected InRecord error, got {other:?}"),
            }
        };
    }

    // --- Field resolution (Task 1) ---

    #[test]
    fn field_reserved_resolves_from_registry() {
        let def = Field::reserved("AF")
            .resolve(FieldKind::Info, LATEST)
            .unwrap();
        assert_eq!(def.id, "AF");
        assert_eq!(def.number, Number::A);
        assert_eq!(def.type_, Type::Float);
    }

    #[test]
    fn field_typed_resolves_explicitly_with_description() {
        let def = Field::typed("DP", Number::ONE, Type::Integer)
            .description("read depth")
            .resolve(FieldKind::Info, LATEST)
            .unwrap();
        assert_eq!(def.id, "DP");
        assert_eq!(def.number, Number::ONE);
        assert_eq!(def.type_, Type::Integer);
        assert_eq!(def.description, "read depth");
    }

    #[test]
    fn field_typed_defaults_description_to_id() {
        let def = Field::typed("DP", Number::ONE, Type::Integer)
            .resolve(FieldKind::Info, LATEST)
            .unwrap();
        assert_eq!(def.description, "DP");
    }

    #[test]
    fn field_flag_resolves_as_info_flag() {
        let def = Field::flag("SOMATIC")
            .resolve(FieldKind::Info, LATEST)
            .unwrap();
        assert_eq!(def.type_, Type::Flag);
        assert_eq!(def.number, Number::FLAG);
    }

    #[test]
    fn field_from_str_is_reserved() {
        let from_str: Field = "AF".into();
        let explicit = Field::reserved("AF");
        assert_eq!(
            from_str.resolve(FieldKind::Info, LATEST).unwrap(),
            explicit.resolve(FieldKind::Info, LATEST).unwrap()
        );
    }

    // --- Builder happy path + deferred validation ---

    #[test]
    fn happy_path_builds() {
        let doc = base()
            .info("AF")
            .format("GT")
            .format(Field::typed("DS", Number::A, Type::Float))
            .record(
                RecordSpec::at("chr1", 1000)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1"])
                    .info("AF", FieldValue::floats([0.25]))
                    .format("DS", [FieldValue::floats([0.4]), FieldValue::floats([1.9])]),
            )
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
        let r = base()
            .format("GT")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .info("AF", FieldValue::floats([0.1])),
            )
            .build();
        assert_in_record!(r, BuildError::UndeclaredField { .. });
    }

    #[test]
    fn cardinality_checked() {
        let r = base()
            .format("GT")
            .info("AF")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()]) // n_alt = 1, AF is Number::A
                    .info("AF", FieldValue::floats([0.1, 0.2])),
            )
            .build();
        assert_in_record!(r, BuildError::Cardinality { .. });
    }

    #[test]
    fn symbolic_requires_svlen_and_padding() {
        // missing SVLEN
        let r = base()
            .format("GT")
            .info("SVLEN")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())]),
            )
            .build();
        assert_in_record!(r, BuildError::MissingSvlen(_));

        // multi-base REF padding violation
        let r = base()
            .format("GT")
            .info("SVLEN")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("AC")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .info("SVLEN", FieldValue::ints([100])),
            )
            .build();
        assert_in_record!(r, BuildError::MissingRefPadding(_));
    }

    #[test]
    fn gt_index_out_of_range() {
        let r = base()
            .format("GT")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|2", "0|0"]),
            )
            .build(); // index 2 > n_alt 1
        assert_in_record!(r, BuildError::AlleleIndexOutOfRange { .. });
    }

    #[test]
    fn gt_not_declared_errs() {
        // .gt(...) used but the builder never declared FORMAT "GT".
        let r = base()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1"]),
            )
            .build();
        assert_in_record!(r, BuildError::GtNotDeclared);
    }

    #[test]
    fn svlen_must_be_missing_for_unspecified() {
        let r = base()
            .format("GT")
            .info("SVLEN")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::unspecified()])
                    .info("SVLEN", FieldValue::ints([100])),
            )
            .build();
        assert_in_record!(r, BuildError::SvlenMustBeMissing(_));
    }

    #[test]
    fn bad_svclaim_errs() {
        let r = base()
            .format("GT")
            .info("SVLEN")
            .info("SVCLAIM")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .info("SVLEN", FieldValue::ints([100]))
                    .info("SVCLAIM", FieldValue::strings(["Z"])),
            )
            .build();
        assert_in_record!(r, BuildError::BadSvclaim { .. });
    }

    #[test]
    fn svclaim_required_for_del_at_4_5() {
        let r = base()
            .format("GT")
            .info("SVLEN")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .info("SVLEN", FieldValue::ints([100])),
            )
            .build();
        assert_in_record!(r, BuildError::SvclaimRequired(_));
    }

    #[test]
    fn too_many_genotypes_errs() {
        let r = base()
            .format("GT")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1", "0|0"]),
            )
            .build();
        assert_in_record!(r, BuildError::SampleCountMismatch { .. });
    }

    #[test]
    fn too_many_format_values_errs() {
        let r = base()
            .format("GT")
            .format(Field::typed("DS", Number::A, Type::Float))
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1"])
                    .format(
                        "DS",
                        [
                            FieldValue::floats([0.1]),
                            FieldValue::floats([0.2]),
                            FieldValue::floats([0.3]),
                        ],
                    ),
            )
            .build();
        assert_in_record!(r, BuildError::SampleCountMismatch { .. });
    }

    #[test]
    fn malformed_gt_errors() {
        let r = base()
            .format("GT")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|x", "0|0"]),
            )
            .build();
        assert_in_record!(r, BuildError::BadGenotype(_));
    }

    #[test]
    fn cn_svlen_mismatch_errs() {
        let r = base()
            .format("GT")
            .format("CN")
            .info("SVLEN")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([
                        Allele::cnv(Vec::<&str>::new()),
                        Allele::cnv(Vec::<&str>::new()),
                    ])
                    .info("SVLEN", FieldValue::ints([100, 200]))
                    .format("CN", [FieldValue::floats([2.0]), FieldValue::floats([3.0])]),
            )
            .build();
        assert_in_record!(r, BuildError::CnSvlenMismatch);
    }

    // --- New guarantees (Task 3) ---

    #[test]
    fn flag_on_format_errs() {
        // Field::flag is INFO-only; using it on FORMAT must fail at build().
        let r = base().format(Field::flag("SOMATIC")).build();
        assert!(matches!(r, Err(BuildError::FlagNotInfo)));
    }

    #[test]
    fn declaration_order_independent() {
        // record() appears before the .format("GT") that it depends on.
        let doc = base()
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1"]),
            )
            .format("GT")
            .build()
            .unwrap();
        assert_eq!(doc.records.len(), 1);
    }

    #[test]
    fn info_str_shorthand_matches_reserved() {
        // .info("AF") and .info(Field::reserved("AF")) produce the same header.
        let a = base().info("AF").build().unwrap();
        let b = base().info(Field::reserved("AF")).build().unwrap();
        assert_eq!(a.info_defs, b.info_defs);
    }

    // --- Record index tagging ---

    #[test]
    fn record_index_in_error() {
        // Second record (index 1) is the bad one.
        let r = base()
            .format("GT")
            .record(
                RecordSpec::at("chr1", 1)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1"]),
            )
            .record(
                RecordSpec::at("chr1", 2)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|2", "0|0"]), // out of range
            )
            .build();
        match r {
            Err(BuildError::InRecord { index, .. }) => assert_eq!(index, 1),
            other => panic!("expected InRecord, got {other:?}"),
        }
    }
}
