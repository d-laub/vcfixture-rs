# Documentation & Examples Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the stale README, add a runnable `examples/` directory, write an mdBook guide deployed to GitHub Pages, and enrich rustdoc so docs.rs reads well.

**Architecture:** Five compile-checked `examples/*.rs` files are the single source of truth for code. An mdBook in `docs/book/` pulls snippets from those examples via `{{#rustdoc_include}}` (anchored regions), so book code never drifts. docs.rs carries the API reference (with the `proptest` feature enabled via Cargo metadata). A new CI workflow builds + tests + deploys the book to GitHub Pages on push to `main`.

**Tech Stack:** Rust 1.86+, mdBook, GitHub Actions (Pages), pixi (conda-forge `mdbook`).

## Global Constraints

- **MSRV:** Rust `1.86` (`rust-version` in `Cargo.toml`). Do not use newer language features.
- **Edition:** 2021.
- **Crate name:** `vcfixture`. Public re-exports live in `src/lib.rs`.
- **Current public API (verified against source):**
  - `VcfBuilder::new(samples: IntoIterator<Item: Into<String>>, contigs: IntoIterator<Item = (Into<String>, Option<u64>)>, version: VcfVersion)`
  - Builder methods (all consuming `self`, returning `VcfBuilder`): `.info(impl Into<Field>)`, `.format(impl Into<Field>)`, `.filter(id, description)`, `.alt(id, description)`, `.record(RecordSpec)`. Terminal: `.build() -> Result<Document, BuildError>` (also `.render()/.write()/.truth()` shortcuts).
  - `Field::reserved(id)`, `Field::typed(id, Number, Type)`, `Field::flag(id)`, `.description(d)`. `&str`/`String` convert to `Field::reserved` via `Into<Field>`.
  - `RecordSpec::at(chrom, pos)` then `.ref_(s)`, `.alt(IntoIterator<Item = Allele>)`, `.ids(..)`, `.qual(f64)`, `.filter(..)`, `.gt(IntoIterator<Item: Into<String>>)`, `.info(id, FieldValue)`, `.format(id, IntoIterator<Item = FieldValue>)`, `.labels(..)`.
  - `Allele::seq(bases) -> Result<Allele, BuildError>`, `Allele::star()`, `Allele::unspecified()`, `Allele::deletion(subtypes)`, `Allele::insertion(subtypes)`, `Allele::duplication(subtypes)`, `Allele::inversion(subtypes)`, `Allele::cnv(subtypes)`, `Allele::breakend_parse(&str) -> Result<Allele, BuildError>`. Symbolic/breakend alleles require a single-base REF pad. Symbolic alleles require `SVLEN`; DEL/DUP require `SVCLAIM` at VCF ≥ 4.4. Breakend/Unspecified/Star must NOT carry `SVLEN`.
  - `FieldValue::ints([i64])`, `FieldValue::floats([f64])`, `FieldValue::strings([Into<String>])`, `FieldValue::Flag`. (`FieldValue`/`Scalar` derive `PartialEq`.)
  - `Number::ONE`, `Number::A`, `Number::fixed(n) -> Result<..>`; `Type::Integer/Float/String/Character/Flag`; `vcfixture::spec::version::LATEST` = `VcfVersion::V4_5` (renders `##fileformat=VCFv4.5`).
  - `Document::render() -> String`, `Document::write(path, WriteOpts) -> Result<PathBuf, BuildError>`, `Document::truth() -> GroundTruth`.
  - `WriteOpts::text()`, `WriteOpts::bgzipped()`, `WriteOpts::bgzipped_indexed()`.
  - `GroundTruth` public fields used here: `genotypes: Array3<i32>` `[records, samples, ploidy]` (allele index, `-1` = missing/pad), `phasing: Array2<bool>` `[records, samples]`, `pos: Array1<i64>` (1-based), `info: Vec<HashMap<String, FieldValue>>`, `format: Vec<Vec<HashMap<String, FieldValue>>>`, `alts_truth: Vec<Vec<AlleleTruth>>`. `AlleleTruth { kind, is_sequence, sv_type: Option<String>, svlen: Option<i64>, sv_end: Option<i64> }`. For a symbolic DEL with SVLEN=100: `sv_type == Some("DEL")`, `svlen == Some(100)`, `is_sequence == false`.
  - Strategies (feature `proptest`): `vcfixture::strategies::{documents, documents_with_fields, symbolic_documents, reference_and_documents, DocumentOpts}`. `documents(DocumentOpts) -> impl Strategy<Value = Document>`. `DocumentOpts::default()` = `{ max_samples: 3, max_records: 4, max_alt: 1, version: LATEST }`.
- **Exact rendered header/data format (from snapshot):**
  ```
  ##fileformat=VCFv4.5
  ##contig=<ID=chr1,length=100000>
  ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
  chr1	1000	.	A	T	.	.	AF=0.25	GT	0|1	1|1
  ```
  Floats render verbatim (`0.25` → `0.25`, `0.5` → `0.5`). Symbolic ALT types are auto-described as `##ALT=<ID=DEL,...>`.
- **Repo conventions:** Conventional Commits (`cz check` in CI). pixi for tooling. prek runs `cargo fmt`/`clippy`/`check` on commit. `actionlint` lints workflows.

---

### Task 1: `examples/core.rs` — build → ground truth → render

**Files:**
- Create: `examples/core.rs`

**Interfaces:**
- Consumes: public API from Global Constraints.
- Produces: anchors `build`, `truth`, `render` (used by Task 7's `building-documents.md` and `ground-truth.md`).

- [ ] **Step 1: Write the example**

Create `examples/core.rs`:

```rust
//! Core workflow: build a tiny VCF document, read its ground-truth oracle, and
//! render it to text.
//!
//! Run with: `cargo run --example core`

// ANCHOR: build
use vcfixture::spec::version::LATEST;
use vcfixture::{Allele, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    // Two samples, one contig, the latest VCF version.
    let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
        .info("AF") // reserved INFO field: Number=A, Type=Float
        .format("GT") // reserved FORMAT field
        .record(
            RecordSpec::at("chr1", 1000)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|1", "1|1"])
                .info("AF", FieldValue::floats([0.25])),
        )
        .build()
        .expect("document is valid");
    // ANCHOR_END: build

    // ANCHOR: truth
    let truth = doc.truth();

    // `genotypes` is an [records, samples, ploidy] array of allele indices
    // (0 = REF, 1 = first ALT, -1 = missing/padding).
    assert_eq!(truth.genotypes[[0, 0, 0]], 0); // record 0, sample s1, allele 0 = REF
    assert_eq!(truth.genotypes[[0, 0, 1]], 1); // s1 allele 1 = first ALT
    assert_eq!(truth.genotypes[[0, 1, 1]], 1); // s2 allele 1 = first ALT

    // `phasing` is [records, samples]: both genotypes used '|'.
    assert!(truth.phasing[[0, 0]]);
    assert!(truth.phasing[[0, 1]]);

    // `pos` is 1-based, per the VCF spec.
    assert_eq!(truth.pos[0], 1000);
    // ANCHOR_END: truth

    // ANCHOR: render
    let text = doc.render();
    assert!(text.starts_with("##fileformat=VCFv"));
    assert!(text.contains("AF=0.25"));
    assert!(text.contains("0|1"));
    print!("{text}");
    // ANCHOR_END: render
}
```

- [ ] **Step 2: Run the example and verify it succeeds**

Run: `cargo run --example core`
Expected: prints the VCF text ending with `chr1\t1000\t.\tA\tT\t.\t.\tAF=0.25\tGT\t0|1\t1|1`, exit code 0 (all asserts pass).

- [ ] **Step 3: Verify it compiles under the test harness**

Run: `cargo test --example core`
Expected: builds; the example's `main` is not a test, so output is `running 0 tests ... ok` after a successful compile.

- [ ] **Step 4: Commit**

```bash
git add examples/core.rs
git commit -m "docs: add core build-to-truth example"
```

---

### Task 2: `examples/fields.rs` — field declarations and typed values

**Files:**
- Create: `examples/fields.rs`

**Interfaces:**
- Consumes: public API; `Field`, `Number`, `Type`, `FieldValue::Flag`.
- Produces: anchor `fields` (used by Task 7's `fields-and-values.md`).

- [ ] **Step 1: Write the example**

Create `examples/fields.rs`:

```rust
//! Declaring INFO/FORMAT fields and constructing typed values.
//!
//! Run with: `cargo run --example fields`

// ANCHOR: fields
use vcfixture::spec::number::Number;
use vcfixture::spec::types::Type;
use vcfixture::spec::version::LATEST;
use vcfixture::{Allele, Field, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
        // Reserved: looked up in the spec registry (AF => Number=A, Type=Float).
        .info("AF")
        // Reserved via the explicit constructor (identical to the &str form).
        .info(Field::reserved("DP"))
        // Typed: you choose Number and Type, plus an optional description.
        .info(Field::typed("AC", Number::A, Type::Integer).description("allele count"))
        // Flag: INFO-only, Number=0, Type=Flag.
        .info(Field::flag("SOMATIC"))
        .format("GT")
        // Per-allele FORMAT field (Number=A).
        .format(Field::typed("DS", Number::A, Type::Float))
        .record(
            RecordSpec::at("chr1", 2000)
                .ref_("G")
                .alt([Allele::seq("C").unwrap()])
                .gt(["0|0", "0|1"])
                .info("AF", FieldValue::floats([0.5])) // list of floats
                .info("DP", FieldValue::ints([42])) // single int (as a 1-elem list)
                .info("AC", FieldValue::ints([1]))
                .info("SOMATIC", FieldValue::Flag) // flag present
                // FORMAT DS: one value per sample.
                .format("DS", [FieldValue::floats([0.4]), FieldValue::floats([1.9])]),
        )
        .build()
        .expect("document is valid");

    // The decoded oracle exposes INFO and FORMAT per record/sample.
    let truth = doc.truth();
    assert_eq!(truth.info[0]["DP"], FieldValue::ints([42]));
    assert_eq!(truth.format[0][1]["DS"], FieldValue::floats([1.9]));

    let text = doc.render();
    // Typed and flag fields render deterministic headers.
    assert!(text.contains("##INFO=<ID=AC,Number=A,Type=Integer,"));
    assert!(text.contains("##INFO=<ID=SOMATIC,Number=0,Type=Flag,"));
    assert!(text.contains("##FORMAT=<ID=DS,Number=A,Type=Float,"));
    // INFO column joins fields with ';' in declaration order.
    assert!(text.contains("AF=0.5;DP=42;AC=1;SOMATIC"));
    print!("{text}");
}
// ANCHOR_END: fields
```

- [ ] **Step 2: Run and verify**

Run: `cargo run --example fields`
Expected: prints VCF text with the four INFO fields and `GT:DS` sample columns; exit 0.

- [ ] **Step 3: Commit**

```bash
git add examples/fields.rs
git commit -m "docs: add fields-and-values example"
```

---

### Task 3: `examples/symbolic.rs` — symbolic SVs, breakends, allele truth

**Files:**
- Create: `examples/symbolic.rs`

**Interfaces:**
- Consumes: `Allele::deletion/insertion/breakend_parse`, `AlleleTruth`.
- Produces: anchor `symbolic` (used by Task 7's `alleles-and-svs.md`).

- [ ] **Step 1: Write the example**

Create `examples/symbolic.rs`:

```rust
//! Symbolic structural variants, breakends, and per-allele truth.
//!
//! Symbolic/breakend ALTs require a single-base REF pad. Symbolic alleles
//! require SVLEN; DEL/DUP additionally require SVCLAIM at VCF >= 4.4. Breakends
//! must NOT carry SVLEN.
//!
//! Run with: `cargo run --example symbolic`

// ANCHOR: symbolic
use vcfixture::spec::version::LATEST;
use vcfixture::{Allele, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    let doc = VcfBuilder::new(["s1"], [("chr1", Some(1_000_000u64))], LATEST)
        .info("SVLEN")
        .info("SVCLAIM")
        .format("GT")
        // Symbolic deletion: single-base REF pad, SVLEN required, and at
        // VCF >= 4.4 a DEL also requires SVCLAIM.
        .record(
            RecordSpec::at("chr1", 5000)
                .ref_("A")
                .alt([Allele::deletion(Vec::<&str>::new())])
                .gt(["0|1"])
                .info("SVLEN", FieldValue::ints([100]))
                .info("SVCLAIM", FieldValue::strings(["D"])),
        )
        // Symbolic insertion: SVLEN required, no SVCLAIM requirement.
        .record(
            RecordSpec::at("chr1", 8000)
                .ref_("C")
                .alt([Allele::insertion(Vec::<&str>::new())])
                .gt(["1|1"])
                .info("SVLEN", FieldValue::ints([50])),
        )
        // Paired breakend: the raw replacement string carries the mate locus.
        // Breakends must NOT carry SVLEN.
        .record(
            RecordSpec::at("chr1", 9000)
                .ref_("G")
                .alt([Allele::breakend_parse("G]chr2:321]").unwrap()])
                .gt(["0|1"]),
        )
        .build()
        .expect("document is valid");

    let truth = doc.truth();
    // The deletion's per-allele truth is decoded for you.
    let del = &truth.alts_truth[0][0];
    assert_eq!(del.sv_type.as_deref(), Some("DEL"));
    assert_eq!(del.svlen, Some(100));
    assert!(!del.is_sequence);

    let text = doc.render();
    assert!(text.contains("<DEL>"));
    assert!(text.contains("<INS>"));
    assert!(text.contains("G]chr2:321]"));
    // Symbolic ALT types are auto-described in the header.
    assert!(text.contains("##ALT=<ID=DEL,"));
    print!("{text}");
}
// ANCHOR_END: symbolic
```

- [ ] **Step 2: Run and verify**

Run: `cargo run --example symbolic`
Expected: prints three data lines (`<DEL>`, `<INS>`, `G]chr2:321]`) plus an `##ALT=<ID=DEL,...>` header; exit 0.

- [ ] **Step 3: Commit**

```bash
git add examples/symbolic.rs
git commit -m "docs: add symbolic SV and breakend example"
```

---

### Task 4: `examples/writing.rs` — write a bgzipped, indexed file

**Files:**
- Create: `examples/writing.rs`

**Interfaces:**
- Consumes: `Document::write`, `WriteOpts::bgzipped_indexed`.
- Produces: anchor `writing` (used by Task 7's `rendering-and-writing.md`).

- [ ] **Step 1: Write the example**

Create `examples/writing.rs`:

```rust
//! Writing a document to a bgzipped, CSI-indexed `.vcf.gz` file.
//!
//! Run with: `cargo run --example writing`

// ANCHOR: writing
use std::env;
use std::fs;

use vcfixture::spec::version::LATEST;
use vcfixture::write::WriteOpts;
use vcfixture::{Allele, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
        .info("AF")
        .format("GT")
        .record(
            RecordSpec::at("chr1", 1000)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|1", "1|1"])
                .info("AF", FieldValue::floats([0.25])),
        )
        .build()
        .expect("document is valid");

    // Render in-memory text whenever you just need a string:
    let _text = doc.render();

    // Or write a bgzipped + CSI-indexed file to disk. The `.gz` extension is
    // ensured for you; `write` returns the final path.
    let dir = env::temp_dir().join("vcfixture_example");
    fs::create_dir_all(&dir).expect("create temp dir");
    let out = doc
        .write(dir.join("out.vcf"), WriteOpts::bgzipped_indexed())
        .expect("write succeeds");

    assert!(out.exists());
    assert_eq!(out.extension().and_then(|e| e.to_str()), Some("gz"));
    // The CSI index sits next to the data file.
    let csi = out.with_extension("gz.csi");
    assert!(csi.exists());
    println!("wrote {} and {}", out.display(), csi.display());
}
// ANCHOR_END: writing
```

- [ ] **Step 2: Run and verify**

Run: `cargo run --example writing`
Expected: prints `wrote /tmp/.../vcfixture_example/out.vcf.gz and /tmp/.../out.vcf.gz.csi`; exit 0. (The writer appends `.gz` to the path, then writes the index as `<path>.gz.csi`; `out.with_extension("gz.csi")` on `out.vcf.gz` yields `out.vcf.gz.csi`.)

- [ ] **Step 3: Commit**

```bash
git add examples/writing.rs
git commit -m "docs: add file-writing example"
```

---

### Task 5: `examples/proptest_fuzz.rs` — fuzzing a parser against the oracle

**Files:**
- Create: `examples/proptest_fuzz.rs`
- Modify: `Cargo.toml` (add `[[example]]` with `required-features`)

**Interfaces:**
- Consumes: `vcfixture::strategies::{documents, DocumentOpts}`, `proptest::test_runner`.
- Produces: anchor `proptest` (used by Task 7's `property-testing.md`).

- [ ] **Step 1: Register the feature-gated example in `Cargo.toml`**

Append to `Cargo.toml` (after the `[features]` block):

```toml
[[example]]
name = "proptest_fuzz"
required-features = ["proptest"]
```

- [ ] **Step 2: Write the example**

Create `examples/proptest_fuzz.rs`:

```rust
//! Fuzzing a VCF parser against the ground-truth oracle.
//!
//! Requires the `proptest` feature:
//!   cargo run --example proptest_fuzz --features proptest
//!
//! In your own test suite you would write this as a `proptest!` block (see the
//! note at the bottom). Here we drive the strategy directly so the example has
//! a runnable `main`.

// ANCHOR: proptest
use proptest::strategy::{Strategy, ValueTree};
use proptest::test_runner::TestRunner;

use vcfixture::strategies::{documents, DocumentOpts};

fn main() {
    let mut runner = TestRunner::default();
    let strategy = documents(DocumentOpts::default());

    for _ in 0..32 {
        // Draw one valid-by-construction document.
        let doc = strategy
            .new_tree(&mut runner)
            .expect("strategy produces a value")
            .current();

        // Derive the oracle and render the document.
        let truth = doc.truth();
        let text = doc.render();

        // A real parser test would parse `text` and compare against `truth`.
        // Here we assert the structural invariant the oracle guarantees: one
        // data line per record in the genotype matrix.
        let data_lines = text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .count();
        assert_eq!(data_lines, truth.pos.len());
        assert_eq!(truth.genotypes.shape()[0], truth.pos.len());
    }

    println!("checked 32 generated documents against their oracle");
}

// In your crate's tests, the idiomatic form is:
//
//   use proptest::prelude::*;
//   use vcfixture::strategies::{documents, DocumentOpts};
//
//   proptest! {
//       #[test]
//       fn parser_matches_oracle(doc in documents(DocumentOpts::default())) {
//           let truth = doc.truth();
//           let text = doc.render();
//           // parse `text` with your parser and prop_assert_eq! against `truth`.
//       }
//   }
// ANCHOR_END: proptest
```

- [ ] **Step 3: Run and verify (feature enabled)**

Run: `cargo run --example proptest_fuzz --features proptest`
Expected: prints `checked 32 generated documents against their oracle`; exit 0.

- [ ] **Step 4: Verify it is correctly gated (feature disabled)**

Run: `cargo build --example proptest_fuzz`
Expected: Cargo reports the example requires the `proptest` feature and does NOT attempt to compile it (no error about missing `proptest` crate).

- [ ] **Step 5: Verify the full test build still passes**

Run: `cargo test --all-features --locked`
Expected: all existing tests pass; every example (including `proptest_fuzz`) compiles.

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml examples/proptest_fuzz.rs
git commit -m "docs: add proptest fuzzing example"
```

---

### Task 6: Enrich rustdoc and enable the `proptest` feature on docs.rs

**Files:**
- Modify: `Cargo.toml` (add `[package.metadata.docs.rs]`)
- Modify: `src/lib.rs` (crate-level overview)
- Modify: `src/build.rs`, `src/value.rs`, `src/allele.rs`, `src/truth.rs`, `src/write.rs` (module-level `//!` docs + doc comments on key public items missing them)

**Interfaces:**
- Consumes: nothing new.
- Produces: nothing other tasks depend on (pure documentation).

- [ ] **Step 1: Enable all-features builds on docs.rs**

Append to `Cargo.toml`:

```toml
[package.metadata.docs.rs]
all-features = true
```

- [ ] **Step 2: Expand the crate-level overview in `src/lib.rs`**

Replace the existing top-of-file `//!` block (lines 1–26) with an overview that explains the oracle concept and links the main types, keeping a compiling doctest. New content:

```rust
//! Generate small VCF test data with a decoded ground-truth oracle.
//!
//! `vcfixture` builds a VCF [`Document`] in code, renders it to text (or a
//! bgzipped, indexed file), and derives a [`GroundTruth`] — arrays of
//! positions, genotypes, and per-allele metadata — so parser tests assert
//! against a known oracle instead of hand-coded literals.
//!
//! # Workflow
//!
//! 1. [`VcfBuilder`] accumulates samples, contigs, field declarations, and
//!    records. It is infallible until [`VcfBuilder::build`], which validates
//!    everything at once and returns a [`Document`] or a [`BuildError`].
//! 2. [`Document::render`] produces VCF text; [`Document::write`] writes a file;
//!    [`Document::truth`] derives the [`GroundTruth`] oracle.
//!
//! Property-test strategies for fuzzing a parser live in [`strategies`], behind
//! the `proptest` feature (off by default).
//!
//! # Example
//!
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

- [ ] **Step 3: Add module-level docs where missing**

For each of `src/build.rs`, `src/value.rs`, `src/allele.rs`, `src/truth.rs`, `src/write.rs`: if the file does not already begin with a `//!` block, add a one-paragraph module overview as the first lines. Examples:

`src/build.rs` (insert at top, before `use`):
```rust
//! The [`VcfBuilder`] accumulator and its `build()` validation pipeline.
//!
//! Field declarations and records are collected infallibly; all validation is
//! deferred to [`VcfBuilder::build`], which resolves reserved fields, checks
//! cardinality and structural-variant rules, and returns a [`Document`].
```

`src/value.rs`:
```rust
//! [`FieldValue`] and [`Scalar`]: decoded INFO/FORMAT values.
```

`src/allele.rs`:
```rust
//! [`Allele`] — sequence, symbolic structural-variant, and breakend ALTs, with
//! constructors and VCF-string rendering/parsing.
```

`src/truth.rs`:
```rust
//! [`GroundTruth`] — the decoded oracle derived from a [`Document`]: position,
//! genotype, phasing, and per-allele arrays a parser test asserts against.
```

`src/write.rs`:
```rust
//! VCF text serialization ([`render`]) and file output ([`write`], optionally
//! bgzipped and CSI-indexed) plus [`WriteOpts`].
```

Use intra-doc links (`[`Type`]`) only for items in scope; bare backticks otherwise.

- [ ] **Step 4: Add doc comments to undocumented public items**

For public items lacking a doc comment in the files above (e.g. `RecordSpec` builder methods, `WriteOpts` constructors, `GroundTruth` fields), add a one-line `///` describing each. Do not document private helpers.

- [ ] **Step 5: Build docs with the same config docs.rs uses**

Run: `cargo doc --all-features --no-deps`
Expected: builds with no warnings; `strategies` appears in the generated docs (`target/doc/vcfixture/strategies/index.html` exists).

- [ ] **Step 6: Verify doctests still pass**

Run: `cargo test --doc --all-features`
Expected: the crate-level example and any other doctests pass.

- [ ] **Step 7: Commit**

```bash
git add Cargo.toml src/lib.rs src/build.rs src/value.rs src/allele.rs src/truth.rs src/write.rs
git commit -m "docs: enrich rustdoc and enable proptest on docs.rs"
```

---

### Task 7: mdBook guide pulling snippets from examples

**Files:**
- Create: `docs/book/book.toml`
- Create: `docs/book/src/SUMMARY.md`
- Create: `docs/book/src/introduction.md`, `building-documents.md`, `fields-and-values.md`, `alleles-and-svs.md`, `ground-truth.md`, `rendering-and-writing.md`, `property-testing.md`
- Modify: `.gitignore` (ignore `docs/book/book/`)
- Modify: `pixi.toml` (add `mdbook` dependency + `docs` tasks)

**Interfaces:**
- Consumes: example anchors `build`, `truth`, `render`, `fields`, `symbolic`, `writing`, `proptest` from Tasks 1–5.
- Produces: built site at `docs/book/book/` (deployed by Task 9).

- [ ] **Step 1: Create `docs/book/book.toml`**

```toml
[book]
title = "vcfixture"
authors = ["David Laub"]
description = "Generate small VCF test data with a decoded ground-truth oracle."
src = "src"
language = "en"

[output.html]
git-repository-url = "https://github.com/d-laub/vcfixture-rs"
edit-url-template = "https://github.com/d-laub/vcfixture-rs/edit/main/docs/book/{path}"
```

- [ ] **Step 2: Create `docs/book/src/SUMMARY.md`**

```markdown
# Summary

- [Introduction](./introduction.md)
- [Building documents](./building-documents.md)
- [Fields and values](./fields-and-values.md)
- [Alleles and structural variants](./alleles-and-svs.md)
- [Ground truth](./ground-truth.md)
- [Rendering and writing](./rendering-and-writing.md)
- [Property testing](./property-testing.md)
```

- [ ] **Step 3: Create the chapter files**

`docs/book/src/introduction.md`:
```markdown
# Introduction

`vcfixture` generates small, spec-conformant VCF (v4.x) test data and returns a
decoded **ground-truth oracle** alongside it. Parser tests assert against the
oracle instead of hand-coded expected arrays.

The primary consumer is property-based testing of VCF/SparseVar parsers. You
build a document in code, render it to text (or a bgzipped, indexed file), and
derive a `GroundTruth` of positions, genotypes, and per-allele metadata.

See the [API reference on docs.rs](https://docs.rs/vcfixture) for full type
signatures. Every code block in this guide is taken verbatim from a compiled
example in the crate's `examples/` directory.
```

`docs/book/src/building-documents.md`:
```markdown
# Building documents

`VcfBuilder` accumulates samples, contigs, field declarations, and records. It
is infallible until `build()`, which validates everything at once.

```rust
{{#rustdoc_include ../../../examples/core.rs:build}}
```

Validation is deferred: declaration order does not matter, and a single
`BuildError` (tagged with the offending record index) is returned from `build()`
if anything is inconsistent.
```

`docs/book/src/fields-and-values.md`:
```markdown
# Fields and values

Declare INFO/FORMAT fields as reserved (looked up in the spec registry), typed
(you choose `Number` and `Type`), or flag (INFO-only). Values are built with
`FieldValue`.

```rust
{{#rustdoc_include ../../../examples/fields.rs:fields}}
```
```

`docs/book/src/alleles-and-svs.md`:
```markdown
# Alleles and structural variants

ALTs are typed `Allele` values: sequence, symbolic SVs (`<DEL>`, `<INS>`, ...),
breakends, `*`, and `<*>`. Symbolic and breakend alleles require a single-base
REF pad; symbolic alleles require `SVLEN`, and DEL/DUP require `SVCLAIM` at
VCF ≥ 4.4.

```rust
{{#rustdoc_include ../../../examples/symbolic.rs:symbolic}}
```
```

`docs/book/src/ground-truth.md`:
```markdown
# Ground truth

`Document::truth()` derives the oracle. `genotypes` is an
`[records, samples, ploidy]` array of allele indices (`-1` = missing/padding),
`phasing` is `[records, samples]`, and `pos` is 1-based.

```rust
{{#rustdoc_include ../../../examples/core.rs:truth}}
```

INFO and FORMAT are decoded per record (and per sample for FORMAT) into maps of
field id to `FieldValue`; per-allele structural metadata lives in `alts_truth`.
```

`docs/book/src/rendering-and-writing.md`:
```markdown
# Rendering and writing

`Document::render()` returns VCF text. `Document::write()` writes a file,
optionally bgzipped and CSI-indexed via `WriteOpts`.

```rust
{{#rustdoc_include ../../../examples/writing.rs:writing}}
```
```

`docs/book/src/property-testing.md`:
```markdown
# Property testing

With the `proptest` feature enabled, `vcfixture::strategies` provides
valid-by-construction `Document` strategies for fuzzing a parser against the
oracle.

```rust
{{#rustdoc_include ../../../examples/proptest_fuzz.rs:proptest}}
```
```

- [ ] **Step 4: Ignore the build output**

Add to `.gitignore`:
```
docs/book/book/
```

- [ ] **Step 5: Add mdBook to pixi**

In `pixi.toml`, add to `[dependencies]`:
```toml
mdbook = "*"
```
And add to `[tasks]`:
```toml
docs-build   = "mdbook build docs/book"
docs-test    = "mdbook test docs/book"
docs-serve   = "mdbook serve docs/book"
```

- [ ] **Step 6: Install the new dependency**

Run: `pixi install`
Expected: resolves and installs `mdbook` from conda-forge.

- [ ] **Step 7: Build the book**

Run: `pixi run docs-build`
Expected: builds with no `unresolved` include errors; produces `docs/book/book/index.html`. Each `{{#rustdoc_include}}` resolves (anchors `build`, `truth`, `fields`, `symbolic`, `writing`, `proptest` exist in the example files).

- [ ] **Step 8: Verify included snippets compile against the crate**

Run: `pixi run docs-test`
Expected: `mdbook test` compiles the rendered Rust blocks. Note: `{{#rustdoc_include}}` regions are shown but `mdbook test` only compiles them as standalone snippets; the authoritative compile check is `cargo test --all-features` (Task 5, Step 5). If `mdbook test` reports failures for the included `main`-bearing snippets, this is expected for `fn main` examples — rely on `cargo run --example` / `cargo test` for correctness and confirm `mdbook build` succeeds.

- [ ] **Step 9: Commit**

```bash
git add docs/book .gitignore pixi.toml pixi.lock
git commit -m "docs: add mdBook guide sourced from examples"
```

---

### Task 8: Rewrite the README

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: current API (the example matches `src/lib.rs`'s doctest).
- Produces: nothing.

- [ ] **Step 1: Replace `README.md` with the current API**

```markdown
# vcfixture

[![CI](https://github.com/d-laub/vcfixture-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/d-laub/vcfixture-rs/actions/workflows/ci.yml)
[![docs.rs](https://img.shields.io/docsrs/vcfixture)](https://docs.rs/vcfixture)
[![Guide](https://img.shields.io/badge/guide-d--laub.github.io-blue)](https://d-laub.github.io/vcfixture-rs/)

Generate small VCF (v4.x) test data with a decoded ground-truth oracle, for
property-testing VCF parsers. Build a `Document` in code, render it to VCF text
(or a bgzipped, indexed file), and get back a `GroundTruth` with arrays of
positions, genotypes, and per-allele metadata — no hand-coded expected arrays.

```rust
use vcfixture::{Allele, Field, RecordSpec, VcfBuilder, FieldValue};
use vcfixture::spec::number::Number;
use vcfixture::spec::types::Type;
use vcfixture::spec::version::LATEST;

let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
    .info("AF")
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

let truth = doc.truth();
assert_eq!(truth.genotypes[[0, 0, 1]], 1);
assert_eq!(truth.pos[0], 1000);
let _text = doc.render();
```

## Examples

Runnable examples live in [`examples/`](examples/):

```bash
cargo run --example core       # build -> truth -> render
cargo run --example fields     # field declarations and typed values
cargo run --example symbolic   # symbolic SVs and breakends
cargo run --example writing    # write a bgzipped, indexed file
cargo run --example proptest_fuzz --features proptest
```

## Proptest strategies

Hypothesis-style strategies for fuzzing a VCF parser are available behind the
`proptest` feature:

```toml
[dev-dependencies]
vcfixture = { version = "0.1", features = ["proptest"] }
```

## Documentation

- [User guide](https://d-laub.github.io/vcfixture-rs/) (mdBook)
- [API reference](https://docs.rs/vcfixture) (docs.rs)
- [Design spec](docs/superpowers/specs/2026-06-23-vcfixture-rs-design.md)
```

- [ ] **Step 2: Verify the README example matches the doctest**

Run: `cargo test --doc --all-features`
Expected: passes (the README snippet mirrors the `src/lib.rs` doctest verified in Task 6).

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: rewrite README for current API"
```

---

### Task 9: CI workflow to build, test, and deploy the book to GitHub Pages

**Files:**
- Create: `.github/workflows/docs.yml`

**Interfaces:**
- Consumes: `docs/book/` (Task 7).
- Produces: a deployed GitHub Pages site.

- [ ] **Step 1: Create `.github/workflows/docs.yml`**

```yaml
name: Docs

on:
  push:
    branches: [main]
  workflow_dispatch:

permissions:
  contents: read
  pages: write
  id-token: write

concurrency:
  group: pages
  cancel-in-progress: false

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Install mdBook
        uses: taiki-e/install-action@v2
        with:
          tool: mdbook
      - name: Test book snippets
        run: mdbook test docs/book
      - name: Build book
        run: mdbook build docs/book
      - uses: actions/upload-pages-artifact@v3
        with:
          path: docs/book/book

  deploy:
    needs: build
    runs-on: ubuntu-latest
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    steps:
      - id: deployment
        uses: actions/deploy-pages@v4
```

- [ ] **Step 2: Lint the workflow**

Run: `pixi run lint-actions`
Expected: `actionlint` reports no errors for `.github/workflows/docs.yml`.

- [ ] **Step 3: Verify the book builds exactly as CI will**

Run: `mdbook build docs/book` (or `pixi run docs-build`)
Expected: produces `docs/book/book/index.html` with no errors.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/docs.yml
git commit -m "ci: build and deploy mdBook to GitHub Pages"
```

- [ ] **Step 5: One-time manual repository setting (document, do not automate)**

In the GitHub repo: **Settings → Pages → Build and deployment → Source → GitHub Actions**. The first push to `main` after this setting is enabled publishes the site at `https://d-laub.github.io/vcfixture-rs/`. Note this in the PR description.

---

## Final verification (run after all tasks)

- [ ] `cargo test --all-features --locked` — all tests + doctests pass; all examples compile.
- [ ] `cargo run --example core && cargo run --example fields && cargo run --example symbolic && cargo run --example writing && cargo run --example proptest_fuzz --features proptest` — each exits 0.
- [ ] `cargo doc --all-features --no-deps` — builds cleanly; `strategies` present.
- [ ] `pixi run docs-build` — book builds; all `{{#rustdoc_include}}` anchors resolve.
- [ ] `cargo fmt --all -- --check` and `cargo clippy --all-features -- -D warnings` — clean (prek runs these on commit).
- [ ] `pixi run lint-actions` — workflows lint clean.
```
