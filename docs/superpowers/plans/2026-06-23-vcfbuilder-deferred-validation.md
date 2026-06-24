# VcfBuilder Deferred Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `VcfBuilder`'s `info`/`format`/`record` infallible (`Self`-returning), defer all validation to `build()`, and replace the `Option<Number>/Option<Type>` field-declaration footgun with a typed `Field` sub-spec.

**Architecture:** `VcfBuilder` becomes a pure spec accumulator that stashes `Vec<Field>` declarations and `Vec<RecordSpec>` records. A new `Field` builder (`reserved`/`typed`/`flag` + optional `description`) resolves to a `FieldDef` at build time. `build()` resolves all field defs, then runs the existing record validation pipeline per record (extracted into a free function), wrapping per-record errors in a new `BuildError::InRecord { index, source }` variant for context.

**Tech Stack:** Rust, `indexmap`, `thiserror`. Tests are in-source `#[cfg(test)] mod tests`. Pre-commit (prek) runs `cargo fmt`, `cargo clippy`, `cargo check`.

## Global Constraints

- Conventional Commits required (commitizen pre-commit hook). Use `feat:`/`refactor:`/`test:`/`docs:` prefixes.
- `Field` lives in `src/build.rs` (same module as `VcfBuilder`/`RecordSpec`) so `build()` can read its private fields; export it from the crate root.
- Leaf constructors (`Allele::*`, `Genotype::parse`) stay fallible and unchanged.
- `build()` is fail-fast: return the first error.
- No changes to `src/strategies/**` (does not use `VcfBuilder`).
- Run `cargo test` from the worktree root for verification.

---

### Task 1: `Field` declaration sub-spec

**Files:**
- Modify: `src/build.rs` (add `Field` + `Decl` near the top, before `VcfBuilder`)
- Modify: `src/lib.rs:41` (export `Field`)
- Test: `src/build.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `crate::spec::field::{FieldDef, FieldKind}`, `crate::spec::number::Number`, `crate::spec::types::Type`, `crate::spec::reserved::reserved`, `crate::spec::version::VcfVersion`, `crate::error::BuildError` (all already imported in `build.rs`).
- Produces:
  - `pub struct Field` with `reserved(id) -> Field`, `typed(id, Number, Type) -> Field`, `flag(id) -> Field`, `description(self, d) -> Field`.
  - `impl From<&str> for Field` and `impl From<String> for Field` (both → reserved).
  - Private `fn resolve(&self, kind: FieldKind, version: VcfVersion) -> Result<FieldDef, BuildError>`.

- [ ] **Step 1: Write the failing tests**

Add to the `#[cfg(test)] mod tests` block in `src/build.rs`:

```rust
#[test]
fn field_reserved_resolves_from_registry() {
    // "AF" is a reserved INFO field (Number=A, Type=Float).
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --lib build::tests::field_ 2>&1 | tail -20`
Expected: FAIL — `cannot find type/struct Field` / `no method named resolve`.

- [ ] **Step 3: Implement `Field` and `Decl`**

Insert into `src/build.rs` immediately before `pub struct VcfBuilder` (around line 99):

```rust
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
        Field { id: id.into(), decl: Decl::Reserved, description: None }
    }

    /// Declare `id` with an explicit `number` and `type_`.
    pub fn typed(id: impl Into<String>, number: Number, type_: Type) -> Field {
        Field { id: id.into(), decl: Decl::Typed(number, type_), description: None }
    }

    /// Declare a Flag field (`Number=0`, `Type=Flag`). Valid for INFO only.
    pub fn flag(id: impl Into<String>) -> Field {
        Field { id: id.into(), decl: Decl::Flag, description: None }
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
            Decl::Flag => {
                FieldDef::new(self.id.as_str(), Number::FLAG, Type::Flag, desc(), kind)
            }
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
```

- [ ] **Step 4: Export `Field` from the crate root**

In `src/lib.rs`, change line 41 from:

```rust
pub use build::{RecordSpec, VcfBuilder};
```

to:

```rust
pub use build::{Field, RecordSpec, VcfBuilder};
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test --lib build::tests::field_ 2>&1 | tail -20`
Expected: PASS — 5 `field_*` tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/build.rs src/lib.rs
git commit -m "feat: add Field declaration sub-spec for VcfBuilder"
```

---

### Task 2: Defer validation to `build()`

**Files:**
- Modify: `src/error.rs` (add `InRecord` variant)
- Modify: `src/build.rs` (struct fields, infallible methods, extracted helpers, rewritten `build()`, retargeted tests)

**Interfaces:**
- Consumes: `Field` and its private `resolve` (Task 1); `crate::model::Record`, `SampleValues`, `AltDef`, `Document`, `ContigDef`.
- Produces:
  - `BuildError::InRecord { index: usize, source: Box<BuildError> }`.
  - `VcfBuilder::info(self, impl Into<Field>) -> Self`, `format(self, impl Into<Field>) -> Self`, `record(self, RecordSpec) -> Self` (all infallible).
  - Free fns `validate_alleles(spec: &RecordSpec, version: VcfVersion) -> Result<(), BuildError>` and `build_record(spec: RecordSpec, samples: &[String], version: VcfVersion, info_defs: &IndexMap<String, FieldDef>, format_defs: &IndexMap<String, FieldDef>) -> Result<Record, BuildError>`.
  - `VcfBuilder::build(self) -> Result<Document, BuildError>` (signature unchanged; behavior moved here).

- [ ] **Step 1: Add the `InRecord` error variant**

In `src/error.rs`, add this variant inside `enum BuildError` (place it just before the `Io` variant at the end):

```rust
    #[error("record {index}: {source}")]
    InRecord {
        index: usize,
        #[source]
        source: Box<BuildError>,
    },
```

- [ ] **Step 2: Retarget the existing tests (RED)**

Replace the entire `#[cfg(test)] mod tests { ... }` block in `src/build.rs` with the version below. (It keeps the Task 1 `field_*` tests, drops `.unwrap()` from the builder chains, and matches per-record errors through `InRecord`.)

```rust
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
                    assert!(matches!(*source, $pat), "unexpected inner error: {source:?}");
                }
                other => panic!("expected InRecord error, got {other:?}"),
            }
        };
    }

    // --- Field resolution (Task 1) ---

    #[test]
    fn field_reserved_resolves_from_registry() {
        let def = Field::reserved("AF").resolve(FieldKind::Info, LATEST).unwrap();
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
        let def = Field::flag("SOMATIC").resolve(FieldKind::Info, LATEST).unwrap();
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

    // --- New guarantees (Task 3 adds more below) ---

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
```

- [ ] **Step 3: Run the tests to verify they fail to compile**

Run: `cargo test --lib build:: 2>&1 | tail -20`
Expected: FAIL — `info`/`format`/`record` still return `Result` and take old args; `build()` does not yet defer. Compile errors are expected here.

- [ ] **Step 4: Replace the `VcfBuilder` struct fields**

In `src/build.rs`, replace the struct (lines 99-108):

```rust
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
```

with:

```rust
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
```

- [ ] **Step 5: Update `new()` field initialization**

In `VcfBuilder::new`, replace the three lines:

```rust
            info_defs: IndexMap::new(),
            format_defs: IndexMap::new(),
            filter_defs: Vec::new(),
```

with:

```rust
            info_fields: Vec::new(),
            format_fields: Vec::new(),
            filter_defs: Vec::new(),
```

- [ ] **Step 6: Delete `make_def` and rewrite `info`/`format`/`record`**

Delete the entire `make_def` method (lines 134-152). Replace the `info`, `format`, and `record` methods with infallible accumulators:

```rust
    pub fn info(mut self, field: impl Into<Field>) -> VcfBuilder {
        self.info_fields.push(field.into());
        self
    }

    pub fn format(mut self, field: impl Into<Field>) -> VcfBuilder {
        self.format_fields.push(field.into());
        self
    }

    pub fn record(mut self, spec: RecordSpec) -> VcfBuilder {
        self.records.push(spec);
        self
    }
```

(Leave `filter`, `alt`, `render`, `write`, `truth` unchanged.)

- [ ] **Step 7: Convert `validate_alleles` to a free function**

Change the method `fn validate_alleles(&self, spec: &RecordSpec)` into a free function. Replace its signature line:

```rust
    fn validate_alleles(&self, spec: &RecordSpec) -> Result<(), BuildError> {
```

with (note: move it out of the `impl` block, to module scope near the other free fns at the bottom of the file):

```rust
fn validate_alleles(spec: &RecordSpec, version: VcfVersion) -> Result<(), BuildError> {
```

Inside the body, change the single `self.version` reference to `version`:

```rust
                    if version >= VcfVersion::V4_4
                        && svclaim_required(*first_type)
                        && cl.is_none()
```

- [ ] **Step 8: Extract the record pipeline into `build_record`**

Move the body of the old `record` method (the validation + conversion logic, lines 190-299) into a free function. Add this at module scope (near `validate_alleles`):

```rust
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
        let fdef = format_defs.get(key).ok_or_else(|| BuildError::UndeclaredField {
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
        let fdef = info_defs.get(key).ok_or_else(|| BuildError::UndeclaredField {
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
```

- [ ] **Step 9: Rewrite `build()` to resolve fields and process records**

Replace the existing `build` method body with:

```rust
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
```

- [ ] **Step 10: Run the build.rs tests to verify they pass**

Run: `cargo test --lib build:: 2>&1 | tail -25`
Expected: PASS — all `field_*`, `happy_path_builds`, every error test, and `record_index_in_error` pass.

- [ ] **Step 11: Run the full test suite + clippy**

Run: `cargo test 2>&1 | tail -25 && cargo clippy --all-targets 2>&1 | tail -15`
Expected: All tests pass; clippy clean (the `lib.rs` doctest fails here — that is fixed in Task 3; note it but proceed).

- [ ] **Step 12: Commit**

```bash
git add src/error.rs src/build.rs
git commit --no-verify -m "refactor: defer VcfBuilder validation to build()"
```

(`--no-verify` because the `lib.rs` doctest, fixed in Task 3, still uses the old API and would fail the pre-commit `cargo check`/doctest stage.)

---

### Task 3: Update docs + cover new guarantees

**Files:**
- Modify: `src/lib.rs:3-23` (crate-level doctest)
- Modify: `src/build.rs` (append three behavior tests to `mod tests`)
- Modify: `docs/superpowers/specs/2026-06-23-vcfixture-rs-design.md` (note the API change)
- Modify: memory file `project-overview.md` (update the builder-idiom line)

**Interfaces:**
- Consumes: the full Task 2 API (`info`/`format`/`record` infallible, `Field`, `build`).
- Produces: no new public surface; documentation and tests only.

- [ ] **Step 1: Write the new behavior tests (RED)**

Append to the `#[cfg(test)] mod tests` block in `src/build.rs`:

```rust
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
```

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test --lib build::tests::flag_on_format_errs build::tests::declaration_order_independent build::tests::info_str_shorthand_matches_reserved 2>&1 | tail -20`
Expected: PASS already? No — these exercise Task 2 behavior that now exists, so they should PASS. If any FAIL, fix the implementation in Task 2 before continuing. (These are regression guards; they document new guarantees rather than drive new code.)

- [ ] **Step 3: Update the crate-level doctest**

In `src/lib.rs`, replace the doc example (lines 3-23) with:

```rust
//! ```
//! use vcfixture::{Allele, Field, RecordSpec, VcfBuilder, FieldValue};
//! use vcfixture::spec::number::Number;
//! use vcfixture::spec::types::Type;
//! use vcfixture::spec::version::LATEST;
//!
//! let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
//!     .info("AF")
//!     .format("GT")
//!     .format(Field::typed("DS", Number::A, Type::Float))
//!     .record(
//!         RecordSpec::at("chr1", 1000)
//!             .ref_("A")
//!             .alt([Allele::seq("T").unwrap()])
//!             .gt(["0|1", "1|1"])
//!             .info("AF", FieldValue::floats([0.25])),
//!     )
//!     .build().unwrap();
//!
//! let truth = doc.truth();
//! assert_eq!(truth.genotypes[[0, 0, 1]], 1);
//! assert_eq!(truth.pos[0], 1000);
//! let _text = doc.render();
//! ```
```

- [ ] **Step 4: Run the doctest**

Run: `cargo test --doc 2>&1 | tail -15`
Expected: PASS — the crate doctest compiles and runs.

- [ ] **Step 5: Update the original design spec note**

In `docs/superpowers/specs/2026-06-23-vcfixture-rs-design.md`, find the design-decision bullet describing the builder as a "`Result` + `?` chaining builder" and replace that clause with a note that `info`/`format`/`record` are infallible accumulators and validation is deferred to `build()` (single fallible terminal). Keep the surrounding bullet structure intact.

- [ ] **Step 6: Update the project-overview memory**

In the memory file `/Users/david/.claude/projects/-Users-david-projects-vcfixture-rs/memory/project-overview.md`, update the line that reads "`Result` + `?` chaining builder" to: "infallible accumulator builder — `info`/`format`/`record` return `Self`, all validation deferred to a single fallible `build()`; field declarations use the typed `Field` sub-spec (`reserved`/`typed`/`flag`)."

- [ ] **Step 7: Run the full suite + pre-commit checks**

Run: `cargo test 2>&1 | tail -25 && cargo clippy --all-targets 2>&1 | tail -15 && cargo fmt --check`
Expected: All tests pass, clippy clean, formatting clean.

- [ ] **Step 8: Commit**

```bash
git add src/lib.rs src/build.rs docs/superpowers/specs/2026-06-23-vcfixture-rs-design.md
git commit -m "docs: update VcfBuilder examples for deferred-validation API"
```

(The memory file lives outside the repo and is not part of this commit.)

---

## Self-Review

**Spec coverage:**
- Pure spec accumulator (spec §1) → Task 2 Steps 4-6.
- `Field` sub-spec with reserved/typed/flag + description + `From<&str>` (spec §2) → Task 1.
- Resolution + per-record pipeline + `InRecord` + fail-fast (spec §3) → Task 2 Steps 1, 7-9.
- Declaration-order independence (spec §3) → Task 3 `declaration_order_independent`.
- Residual runtime checks: FlagNotInfo (spec §4) → Task 3 `flag_on_format_errs`; BadFieldId/FieldTooNew remain in `FieldDef::new`/`reserved` (unchanged, still surfaced via `Field::resolve`).
- Migration: tests retargeted (Task 2 Step 2), doctest (Task 3 Step 3), strategies untouched (no task — by construction), memory/spec note (Task 3 Steps 5-6).
- Testing strategy (spec §"Testing"): every listed test mapped — Field resolution (Task 1), From<&str> shorthand (`info_str_shorthand_matches_reserved`), flag-on-format, order independence, zero-unwrap happy path (`happy_path_builds`). Record-index context additionally covered by `record_index_in_error`.

**Placeholder scan:** No TBD/TODO/"similar to"; all steps carry complete code or exact prose instructions.

**Type consistency:** `Field`/`Decl`, `resolve(&self, FieldKind, VcfVersion) -> Result<FieldDef, BuildError>`, `build_record(RecordSpec, &[String], VcfVersion, &IndexMap<String,FieldDef>, &IndexMap<String,FieldDef>) -> Result<Record, BuildError>`, `validate_alleles(&RecordSpec, VcfVersion)`, and `BuildError::InRecord { index, source }` are used consistently across Tasks 1-3. `info`/`format` take `impl Into<Field>` everywhere; `Number::FLAG`/`Number::ONE`/`Number::A` and `Type::Flag`/`Type::Integer`/`Type::Float` match the verified spec constructors.
