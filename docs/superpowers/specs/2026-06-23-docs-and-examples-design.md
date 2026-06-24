# Documentation & Examples — Design

**Date:** 2026-06-23
**Status:** Approved

## Goal

Bring `vcfixture-rs` documentation up to date and complete:

1. Fix the stale `README.md` (its example uses the old fallible builder API).
2. Add a Rust `examples/` directory of runnable, compile-checked usage examples.
3. Write an mdBook prose guide and publish it to GitHub Pages.
4. Make docs.rs (rustdoc) carry the API reference, with the `proptest` feature
   visible.

The split is fully Rust-native: **mdBook** for narrative guides, **docs.rs** for
the API reference. No second documentation toolchain (the Python sibling uses
zensical/mkdocstrings; we deliberately diverge to idiomatic Rust tooling).

## Background

- The crate (`vcfixture` v0.1.0, unpublished) exposes an infallible accumulator
  builder: `VcfBuilder::new(...).info(...).format(...).record(...).build()`
  returning `Result<Document, BuildError>`, with all validation deferred to
  `build()`. `Document` then offers `render()`, `write()`, and `truth()`.
- Feature-gated `proptest` strategies live in `src/strategies.rs` (entry point
  `documents()`), off by default.
- The current `README.md` example is **stale**: it calls
  `.info("AF", None, None, None).unwrap()` (the old fallible signature). The
  current API is `.info("AF")` (infallible, `impl Into<Field>`).
- There is **no `examples/` directory**.
- The Python `vcfixture` ships a zensical docs site (Home / Guide /
  API Reference) on GitHub Pages with a badged README. We mirror the *shape* of
  that experience, not the tooling.
- Existing CI (`.github/workflows/ci.yml`) runs lint + a test matrix
  (`cargo test --all-features --locked`). No docs deployment exists.

## Components

### 1. `examples/` directory — the source of truth for snippets

Five runnable files. Each has a header doc-comment explaining the workflow, a
`fn main()` that prints output, and closing `assert!`s so it doubles as a smoke
test. The existing `cargo test --all-features` already compiles every example,
so they are CI-verified without extra wiring.

| File | Covers |
|------|--------|
| `examples/core.rs` | `VcfBuilder` with samples/contigs → declare `AF` (INFO) and `GT` (FORMAT) → one record → `build()` → inspect `GroundTruth` (`genotypes`, `pos`, `phasing`) → `render()` to text. The "hello world". |
| `examples/fields.rs` | `Field::reserved` vs `Field::typed` vs `Field::flag`; `FieldValue` construction (ints/floats/strings/lists); multi-sample FORMAT (`DS`, `Number=A`); cardinality. |
| `examples/symbolic.rs` | `Allele::seq`, `deletion`/`dup`/`cnv`/`ins`/`inv` symbolic SVs, breakends; `SVLEN`/`SVCLAIM`; ref-padding rules; version differences (4.4+ requires SVCLAIM for DEL/DUP). |
| `examples/writing.rs` | `Document::write()` to a bgzipped + CSI-indexed `.vcf.gz` via `WriteOpts`, into a temp directory. |
| `examples/proptest_fuzz.rs` | `strategies::documents()` driving a `proptest!` block that round-trips a generated document against its oracle. Gated with `[[example]] name = "proptest_fuzz" required-features = ["proptest"]` in `Cargo.toml`. |

### 2. mdBook guide (`docs/book/`)

Layout:

```
docs/book/
  book.toml
  src/
    SUMMARY.md
    introduction.md
    building-documents.md
    fields-and-values.md
    alleles-and-svs.md
    ground-truth.md
    rendering-and-writing.md
    property-testing.md
  book/        # build output, gitignored
```

Chapters:

- **Introduction** — what the crate is, the ground-truth-oracle idea, the
  primary consumer (svar2 property tests).
- **Building documents** — `VcfBuilder`, `RecordSpec`, deferred validation,
  `build()` and its `BuildError`.
- **Fields & values** — reserved/typed/flag `Field`s, `FieldValue`.
- **Alleles & structural variants** — sequence + symbolic alleles, SV rules.
- **Ground truth** — the `GroundTruth` arrays as the parser oracle.
- **Rendering & writing** — `render()` vs `write()` (+ bgzip/index).
- **Property testing** — the `proptest` feature and `strategies`.

**Snippet strategy:** chapter code is pulled from `examples/*.rs` using mdBook's
`{{#rustdoc_include ../../../examples/<file>.rs:<anchor>}}`. Because the example
files are compiled by `cargo test`, every snippet shown in the book is real,
compiled code — no drift. Examples use `// ANCHOR` / `// ANCHOR_END` comments to
mark the regions a chapter includes. `mdbook test docs/book` additionally checks
any inline ```` ```rust ```` blocks against the crate.

Each chapter links out to docs.rs for the relevant API symbols rather than
duplicating signatures.

`book.toml` essentials: `title`, `authors`, `git-repository-url`, and
`[output.html]` with `git-repository-url` / `edit-url-template` pointing at the
GitHub repo so the "edit" and "view source" actions work.

### 3. README rewrite

- Replace the stale example with the current infallible-API example (kept in
  sync with `lib.rs`'s module doc-test).
- Add badges: CI status, docs.rs, and the Book (GitHub Pages).
- Quick example + links to the Book and docs.rs, mirroring the Python sibling's
  README shape.
- Keep the `proptest` feature note and the design-spec link.

### 4. docs.rs (rustdoc)

- Add to `Cargo.toml`:

  ```toml
  [package.metadata.docs.rs]
  all-features = true
  ```

  so the `proptest`-gated `strategies` module is built and shown on docs.rs.
- Fill gaps in module-level (`//!`) and public-item doc comments across the main
  surface: `build`, `truth`, `value`, `allele`, `write`, `strategies`,
  `reference`. Enough that docs.rs reads well end to end. Private helpers are not
  exhaustively documented.

### 5. CI: build, test, and deploy the book to GitHub Pages

New `.github/workflows/docs.yml`:

- Trigger: `push` to `main` (plus `workflow_dispatch` for manual runs).
- Permissions: `contents: read`, `pages: write`, `id-token: write`.
- Concurrency group `pages` with `cancel-in-progress: false` (standard Pages
  pattern).
- Steps: checkout → install mdBook (pinned version, e.g. via
  `taiki-e/install-action` or a cached cargo-binstall) → `mdbook test docs/book`
  → `mdbook build docs/book` → `actions/upload-pages-artifact` (path
  `docs/book/book`) → `actions/deploy-pages`.
- Examples remain compiled by the existing `test` job (`cargo test
  --all-features` builds all examples); no change needed there.

**One-time manual step (documented, not automated):** set the repository's
Pages **Source → GitHub Actions** in repo settings.

### 6. `.gitignore`

Add `docs/book/book/` (the mdBook build output).

## Out of scope (YAGNI)

- crates.io publishing.
- Versioned / multi-version docs.
- Custom mdBook themes or plugins beyond defaults.
- Documentation coverage gating in CI.

## Verification

- `cargo test --all-features` passes (compiles all examples + doc-tests).
- Each example runs: `cargo run --example core` (etc.), proptest one with
  `--features proptest`.
- `mdbook build docs/book` and `mdbook test docs/book` succeed locally.
- README example compiles (it matches the `lib.rs` doc-test, which `cargo test`
  checks).
- After merge to `main`, the docs workflow deploys and the GitHub Pages site is
  reachable.
