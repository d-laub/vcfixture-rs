# vcfixture-rs design

A Rust port of [vcfixture](https://github.com/d-laub/vcfixture): generate small,
spec-conformant VCF (v4.1–4.5) test data **and return the decoded ground truth
alongside**, so parser tests assert against a known oracle instead of
hand-maintained literals. The immediate consumer is property-based testing of
**svar2** (the planned Rust rewrite of genoray's SparseVar `.svar` reader/writer).

## Goals

- Full feature parity with the Python library's *capabilities*, expressed with an
  idiomatic Rust API (enums, `Result`, ownership, proptest).
- One immutable model is the hub: a builder (or proptest strategies) produces it;
  `render()`/`write()` serialize it; `truth()` derives the oracle. The bytes and
  the oracle come from the same object, so they can never diverge.
- Drive svar2's property tests: draw a document, write it, read it back with an
  independent reader, assert the decode equals the document's `truth()`.

## Non-goals

- **No cross-language byte or oracle parity with the Python library.** The two
  ports are independent. The Rust port is *internally* consistent (its bytes and
  its truth agree), which is all svar2's tests need. Exact byte output and
  float formatting will differ from Python (different serializer, different RNG).
  If Python↔Rust oracle snapshots are wanted later, that is a serde add-on we can
  revisit; it is out of scope here.
- No VCF *parsing* into our model (we generate, we don't ingest). Reading back is
  the consumer's job, done with an independent reader (noodles-vcf) in tests.

## Architecture

The Python design's central virtue is that **one immutable model is the hub**.
We preserve that exactly. Our own idiomatic model is the single source of truth;
noodles is an output adapter only.

```
                      VcfBuilder ──┐
proptest strategies ──────────────┤
                                   ▼
                            Document (model.rs)   ← the hub, owned by us
                              │            │
                   truth()    │            │   render()/write()
                              ▼            ▼
                       GroundTruth     noodles adapter (write.rs)
                       (ndarray)       → text / .vcf.gz + .csi via noodles
```

**Decision: own model as hub, noodles as output adapter (chosen).** Keep
idiomatic Rust structs/enums as the canonical representation. `truth()` derives
from them directly; the writer converts `Document → noodles RecordBuf + Header`
and uses noodles for text rendering, bgzip compression, and CSI/tabix indexing.
Truth derivation stays decoupled from noodles and independently unit-testable.

Rejected alternative: use noodles `RecordBuf` *as* the model. It couples truth
derivation to noodles' API and drops the typed-allele distinction the builder
relies on for validation.

If noodles' record writer cannot express a deliberately edge-case construct, a
small custom text-rendering fallback is the escape hatch (documented, expected to
be unused).

## Crate layout

```
vcfixture/  (crate "vcfixture")
  Cargo.toml
  src/
    lib.rs            re-exports + crate docs
    error.rs          BuildError (thiserror)
    spec/
      mod.rs
      number.rs       Number, NumberKind, cardinality()
      types.rs        Type
      version.rs      VcfVersion (Ord-derived, chronological); LATEST
      field.rs        FieldDef (+ header_line, validation)
      reserved.rs     reserved-field registry + version gating
      genotype_order.rs   Number=G ordering
    allele.rs         Allele enum + ctors + parse + classify
    genotype.rs       Genotype (parse/render/ploidy/is_phased)
    variants.rs       classify() / record_class() / class ctors
    model.rs          Document, Record, ContigDef, AltDef (the hub)
    build.rs          VcfBuilder + eager validation
    truth.rs          GroundTruth (ndarray) + AlleleTruth + derive()
    write.rs          render() -> String; write(path, bgzip, index) via noodles
    reference.rs      ReferenceSpec, ReferenceBuilder, RepeatFeature, FASTA write
    strategies.rs     proptest strategies (feature = "proptest")
  tests/
    roundtrip.rs      proptest: draw → write → read (noodles) → assert vs truth
    snapshots.rs      insta snapshots of rendered VCF text
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `ndarray` | GroundTruth arrays (`Array1/Array2/Array3`) |
| `noodles-vcf`, `noodles-bgzf`, `noodles-csi`, `noodles-tabix`, `noodles-fasta`, `noodles-core` | render / bgzip / index / FASTA write |
| `indexmap` | order-preserving INFO/FORMAT maps in the model |
| `thiserror` | `BuildError` |
| `rand`, `rand_chacha` | seeded reference fill (core: `ReferenceBuilder`) |
| `proptest` *(optional)* | strategies; behind feature `proptest` |

**Features.** `proptest` (optional) pulls in the `proptest` crate. It is **off by
default**, so a consumer that only wants the builder + oracle stays light.
genoray's svar2 tests enable it (`vcfixture = { version = "...", features =
["proptest"] }`). noodles I/O and the reference subsystem (with its seeded
`rand_chacha` fill) are always on (core to `render`/`write`/`ReferenceBuilder`).

## Core types (idiomatic redesign of the Python model)

### Allele (`allele.rs`)

One enum replaces the Python class hierarchy and its `Seq/Sym/Star/...` aliases:

```rust
pub enum SvType { Del, Ins, Dup, Inv, Cnv }

pub enum Allele {
    Seq(String),                                   // validated [ACGTNacgtn]+
    Star,                                          // spanning deletion "*"
    Symbolic { first_type: SvType, subtypes: Vec<String> },
    Unspecified,                                   // "<*>"
    Breakend { raw: String, single: bool },
}
```

- Constructors: `Allele::seq("A") -> Result<Allele, BuildError>` (validates
  `[ACGTNacgtn]+`); `Allele::deletion(subtypes)`, `::insertion`, `::duplication`,
  `::inversion`, `::cnv`; `Allele::breakend_parse("T[chr2:5[") -> Result<...>`
  (the paired/single regex dispatch); `Allele::star()`, `Allele::unspecified()`.
- `render(&self) -> String` and `Allele::parse(&str) -> Allele` (the
  `classify_allele` syntactic dispatcher: `*`, `<*>`, `<...>`, breakend brackets,
  single-breakend `.t`/`t.`, else sequence).
- `Symbolic::type_str()` → inner token, e.g. `DEL` or `DUP:TANDEM`.

### Genotype (`genotype.rs`)

```rust
pub struct Genotype { pub alleles: Vec<Option<u32>>, pub phased: Vec<bool> }
```

`Genotype::parse("0|1")`, `render()`, `ploidy()`, `is_phased()` (true iff there
is at least one separator and all separators are `|`). `None` = missing allele.

### Spec (`spec/`)

- `NumberKind { Fixed, A, R, G, Dot, Flag }`; `Number { kind, count: Option<u32> }`
  with `Number::ONE/A/R/G/DOT/FLAG` consts, `Number::fixed(n)`, `header_str()`,
  and `cardinality(n_alt, ploidy) -> Option<u32>` (G via binomial
  `C(n_alleles + ploidy - 1, ploidy)`).
- `Type { Integer, Float, Flag, Character, String }` with `info_allowed()` /
  `format_allowed()` (Flag excluded from FORMAT).
- `VcfVersion { V4_1..=V4_5 }`, `#[derive(PartialOrd, Ord)]` in chronological
  order; `LATEST = V4_5`.
- `FieldDef { id, number, type, description, kind: FieldKind }` with the
  ID-regex check and the Flag↔Number=0 cross-validation; `header_line()`.
  `FieldKind { Info, Format }`.
- `reserved(id, kind, version) -> Result<FieldDef, BuildError>`: the registry
  (AA/AC/AF/AN/DP/DB/H2/END/SVTYPE/SVLEN/SVCLAIM/CIPOS/CIEND/CILEN/MATEID/PARID/
  IMPRECISE; GT/GQ/DP/AD/PL/GL/PS/CN/LEN) with version gating (SVCLAIM & LEN since
  4.4) and the SVLEN 4.3→4.4 form switch (Number=. signed → Number=A unsigned).
- `genotype_ordering(ploidy, n_alleles) -> Vec<Vec<u32>>` (Number=G order).

### Decoded values (`truth.rs`)

```rust
pub enum Scalar { Int(i64), Float(f64), Char(char), Str(String) }
pub enum FieldValue { Flag, Scalar(Scalar), List(Vec<Option<Scalar>>) }
```

`FieldValue` is what INFO/FORMAT decode to in the oracle (a Flag, a lone scalar,
or a list with possible missing entries).

## Model (`model.rs`) — the hub

```rust
pub struct ContigDef { pub id: String, pub length: Option<u64> }
pub struct AltDef    { pub id: String, pub description: String }

pub struct Record {
    pub chrom: String,
    pub pos: u64,                          // 1-based
    pub ids: Option<Vec<String>>,          // None -> "."
    pub ref_: String,
    pub alts: Vec<Allele>,
    pub qual: Option<f64>,
    pub filters: Option<Vec<String>>,      // None -> "."; empty -> "PASS"
    pub info: IndexMap<String, FieldValue>,
    pub fmt_keys: Vec<String>,             // FORMAT column order
    pub samples: Vec<SampleValues>,        // per-sample key -> value / Genotype
    pub labels: BTreeSet<String>,          // test tags, not serialized
}

pub struct Document {
    pub version: VcfVersion,
    pub info_defs: Vec<FieldDef>,
    pub format_defs: Vec<FieldDef>,
    pub filter_defs: Vec<(String, String)>,
    pub contigs: Vec<ContigDef>,
    pub samples: Vec<String>,
    pub records: Vec<Record>,
    pub alt_defs: Vec<AltDef>,
}
```

`Document::max_ploidy()`, `Document::render()`, `Document::truth()`,
`Document::write(path, opts)`. Per-sample values keep insertion order (GT plus
declared FORMAT keys); `SampleValues` holds an optional `Genotype` and an ordered
map of other FORMAT values. Use `indexmap::IndexMap` to preserve field order in
INFO/FORMAT for deterministic rendering (add `indexmap` to deps).

## Builder (`build.rs`)

Infallible-accumulator design: `info`/`format`/`record` return `Self` (no
`Result`), accumulating declarations without any validation. All validation is
deferred to the single fallible terminal `build() -> Result<Document, BuildError>`.
Field declarations use the typed `Field` sub-spec (`Field::reserved(id)`,
`Field::typed(id, number, type_)`, `Field::flag(id)`); a bare `&str` or `String`
converts via `From<&str> for Field` (resolves as reserved).

```rust
let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000))], VcfVersion::LATEST)
    .info("AF")                                      // reserved → Number/Type resolved at build()
    .format("GT")
    .format(Field::typed("DS", Number::A, Type::Float))
    .record(
        RecordSpec::at("chr1", 1000)
            .ref_("A")
            .alt([Allele::seq("T").unwrap()])
            .gt(["0|1", "1|1"])
            .info("AF", FieldValue::floats([0.25])),
    )
    .build().unwrap();

let text  = doc.render();
let truth = doc.truth();
doc.write("x.vcf.gz", WriteOpts { bgzip: true, index: true }).unwrap();
```

- `info(impl Into<Field>)` / `format(impl Into<Field>)`: accept a `Field` or
  any `Into<Field>` (bare `&str`/`String` → `Field::reserved`). `filter(id, desc)`,
  `alt(id, desc)`.
- Declaration order is independent of record order: `record()` may appear before
  the `format("GT")` it depends on; `build()` resolves everything in one pass.
- A small **`Record` sub-builder** stands in for Python's keyword arguments
  (Rust has no kwargs). `Record::at(chrom, pos)` then chained setters
  `.ref_()/.alt()/.ids()/.qual()/.filter()/.gt()/.info()/.format()/.labels()`.
  `VcfBuilder::record(rec)` accumulates a record spec; validation runs later in `build()`.
- **Deferred validation** (run in `build()`) ports every Python check, surfaced as `BuildError`:
  - reserved resolution failure when number/type omitted for an unknown id;
  - symbolic/breakend ALT requires a single-base REF padding base;
  - symbolic allele requires `SVLEN`; SVCLAIM per-type allow-list
    (`DEL/DUP: D/J/DJ`, `CNV: D`, `INS/INV: J`); SVCLAIM required for `DEL`/`DUP`
    at version ≥ 4.4; SVLEN must be absent for breakend/unspecified/spanning;
  - INFO/FORMAT field must be declared before use;
  - value count must equal the resolved `Number` cardinality (A per-ALT, R
    per-allele incl. REF, G per-genotype, Fixed exact);
  - GT allele index must not exceed `n_alt`;
  - `FORMAT CN` requires equal SVLEN across `<CNV>/<DEL>/<DUP>` alleles.
- `build()` auto-registers an `AltDef` for each symbolic ALT type encountered
  (default description), with explicit `.alt()` descriptions overriding.

Inside proptest strategies, generation is correct by construction, so the
builder's `Result` is `.expect()`-ed there (mirrors Hypothesis: no rejection).

## Ground truth (`truth.rs`) — ndarray oracle

Record-first indexing, `-1` = missing, mirroring the Python field set:

```rust
pub struct AlleleTruth {
    pub kind: AlleleKind,          // Snp|Mnp|Ins|Del|Delins|SpanningDel|Symbolic|Unspecified|Bnd
    pub is_sequence: bool,         // literal DNA a tool may splice
    pub sv_type: Option<String>,   // e.g. "DEL"/"DUP:TANDEM"
    pub svlen: Option<i64>,        // resolved per-allele, absolute
    pub sv_end: Option<i64>,       // 1-based inclusive end = pos + svlen for DEL/DUP/INV/CNV
}

pub struct GroundTruth {
    pub samples: Vec<String>,
    pub contigs: Vec<String>,                       // contig id per record
    pub pos: Array1<i64>,                           // (records,), 1-based
    pub ref_: Vec<String>,
    pub alts: Vec<Vec<String>>,                     // rendered per ALT
    pub variant_class: Vec<VariantClass>,
    pub genotypes: Array3<i32>,                     // (records, samples, ploidy), -1 missing
    pub phasing: Array2<bool>,                      // (records, samples)
    pub info: Vec<HashMap<String, FieldValue>>,     // per record
    pub format: Vec<Vec<HashMap<String, FieldValue>>>, // per record, per sample (GT excluded)
    pub labels: Vec<BTreeSet<String>>,
    pub alts_truth: Vec<Vec<AlleleTruth>>,
    pub is_sequence_mask: Vec<Array1<bool>>,        // per record, over ALTs
}
```

`VariantClass` is an enum covering `Snp Mnp Ins Del Delins SpanningDel Unspecified
Bnd Multiallelic SvDel SvIns SvDup SvInv Cnv` (port of `classify` + `record_class`:
multiple ALTs → `Multiallelic`; single symbolic → `Cnv` for CNV else `Sv*`).
`derive(&Document) -> GroundTruth` ports `derive_truth` + `_allele_truth`
(SVLEN→END only for the spanning SV types DEL/DUP/INV/CNV; INS has length but no
span).

## Reference subsystem (`reference.rs`)

- `RepeatFeature { contig, pos0, motif, count }` with `length()`.
- `ReferenceBuilder::new(seed)` — seeded `rand_chacha::ChaCha8Rng` fill of
  contigs from `ACGT`; `add_contig(id, length)`, `set_base`, `set_seq`,
  `tandem_repeat(contig, pos0, motif, n)`; `build() -> ReferenceSpec`.
- `ReferenceSpec { contigs: Vec<(String, String)>, repeats: Vec<RepeatFeature> }`
  with `base`/`seq`/`length`/`draw_ref_alt(contig, pos0, klass, opts)` and
  `write(path, bgzip, index)` (60-column FASTA; `.fai` via noodles-fasta, `.gzi`
  via noodles-bgzf when bgzipped).
- Note: seeded fill uses a Rust RNG, so generated sequences differ from Python's
  numpy RNG. Acceptable — no cross-language parity goal.

## I/O (`write.rs`)

- `Document::render() -> String`: drive the noodles VCF writer into an in-memory
  `Vec<u8>` buffer, return as `String`.
- `Document::write(path, WriteOpts { bgzip, index }) -> Result<PathBuf>`:
  - plain text `.vcf` when `!bgzip`;
  - bgzipped `.vcf.gz` (noodles-bgzf) when `bgzip` (append `.gz` if absent);
  - `.csi` index alongside (noodles-csi/tabix) when `index` (ignored if `!bgzip`).
- Conversion `Document → noodles Header + RecordBuf` lives here. The custom
  text-fallback escape hatch lives here too if any construct proves
  inexpressible via noodles' writer.

## proptest strategies (`strategies.rs`, feature `proptest`)

Each Hypothesis composite becomes a function returning `impl Strategy<Value =
...>`, built from proptest combinators (`prop_flat_map` for dependent dimensions
such as n_alt→cardinality; `prop::collection::vec` for variable-length
collections) so shrinking works end-to-end through the same builder.

- `genotypes(ploidy, n_alt, missing_rate) -> impl Strategy<Value = String>`
- `field_value(&FieldDef, n_alt, ploidy) -> impl Strategy<Value = FieldValue>`
- `documents(DocumentOpts) -> impl Strategy<Value = Document>` (reference-free
  body and reference-consistent body, with the `violations` /
  `label_overrides` options: `multiallelic`, `non_atomic`, `non_left_aligned`)
- `documents_with_fields(...)` — full Number×Type matrix as INFO + FORMAT (+ Flag
  INFO)
- `symbolic_documents(...)` — symbolic SV / `<*>` records with consistent
  SVLEN/SVCLAIM and the version-correct SVLEN sign
- `references(...)` — synthetic `ReferenceSpec` with non-overlapping planted
  repeats
- `reference_and_documents(...)` — consistent `(ReferenceSpec, Document,
  GroundTruth)` triple
- Coverage tables as `const`/`static`: `ALL_VARIANT_CLASSES`,
  `ALL_NUMBER_TYPE_COMBOS`, for `prop::sample::select` parametrization; reserve
  proptest for values within a chosen combo.

Options structs (`DocumentOpts`, etc.) replace Python's keyword defaults, with a
`Default` impl matching the Python defaults (e.g. `max_samples=3`,
`max_records=4`, `max_alt=1`).

## Error handling (`error.rs`)

A single `BuildError` enum via `thiserror`, with variants for each validation
failure (bad allele bases, unknown reserved field, missing/invalid SVLEN/SVCLAIM,
undeclared field, cardinality mismatch, GT index out of range, CN/SVLEN
mismatch, missing REF padding, version-unavailable field, I/O error). Public
fallible API returns `Result<_, BuildError>`.

## Testing strategy

- **Unit tests** per module: allele parse/render round-trips; `Number`
  cardinality (incl. G binomial); genotype parse/render/`is_phased`; reserved
  version-gating and SVLEN form switch; `derive` on hand-built documents
  (genotypes/phasing/info/format/variant-class/allele-truth).
- **Round-trip property tests** (`tests/roundtrip.rs`): draw `documents()` (and
  the other strategies), `write()`, read back with **noodles-vcf as an
  independent reader**, assert the decode equals `truth()`. This is exactly the
  loop svar2 will run against its own reader.
- **Snapshots** (`tests/snapshots.rs`): `insta` snapshots of a few rendered VCFs
  to lock byte-level output (add `insta` as a dev-dependency).

## Dev environment

### pixi (`pixi.toml`)

conda-forge provides the Rust toolchain and prek; mirror the genoray conventions.

```toml
[workspace]
name = "vcfixture-rs"
channels = ["conda-forge"]
platforms = ["osx-arm64", "linux-64"]

[dependencies]
rust = "*"            # cargo, rustc, rustfmt, clippy
prek = "*"
commitizen = "*"

[tasks]
build       = "cargo build"
test        = "cargo test --all-features"
check       = "cargo check --all-features"
fmt         = "cargo fmt"
fmt-check   = "cargo fmt --check"
clippy      = "cargo clippy --all-features -- -D warnings"
prek-install = "prek install -t pre-commit -t commit-msg"
bump-dry    = "cz bump --dry-run"
```

(If the conda-forge `rust` package turns out not to bundle `clippy`/`rustfmt` on
a target platform, add `rust-std`/the component package during implementation —
verify with `pixi run cargo clippy --version` at setup time.)

### prek (`.pre-commit-config.yaml`)

Run cargo check, fmt, and clippy on every commit. Cargo hooks are `language:
system` (pixi supplies cargo); commitizen and the generic hooks use their public
repos, matching the Python projects.

```yaml
default_install_hook_types: [pre-commit, commit-msg]

repos:
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v6.0.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-toml
      - id: check-merge-conflict
      - id: check-case-conflict
      - id: check-added-large-files
      - id: mixed-line-ending
        args: [--fix=lf]

  - repo: https://github.com/commitizen-tools/commitizen
    rev: v4.16.3
    hooks:
      - id: commitizen
        stages: [commit-msg]

  - repo: local
    hooks:
      - id: cargo-fmt
        name: cargo fmt
        entry: cargo fmt --check
        language: system
        types: [rust]
        pass_filenames: false
      - id: cargo-clippy
        name: cargo clippy
        entry: cargo clippy --all-features -- -D warnings
        language: system
        types: [rust]
        pass_filenames: false
      - id: cargo-check
        name: cargo check
        entry: cargo check --all-features
        language: system
        types: [rust]
        pass_filenames: false
```

## Implementation order

1. Scaffold: `Cargo.toml`, `pixi.toml`, `.pre-commit-config.yaml`, `lib.rs`,
   `error.rs`; `pixi install`; `prek install`.
2. `spec/` (version, type, number, field, reserved, genotype_order) + unit tests.
3. `allele.rs`, `genotype.rs`, `variants.rs` + unit tests.
4. `model.rs` (the hub).
5. `build.rs` (`VcfBuilder` + `Record` sub-builder + validation) + unit tests.
6. `truth.rs` (`derive` + `GroundTruth` + `AlleleTruth`) + unit tests.
7. `write.rs` (noodles render/write/index) + snapshot tests.
8. `reference.rs` + unit tests.
9. `strategies.rs` (feature `proptest`) + `tests/roundtrip.rs`.
10. Docs pass on the public API.
