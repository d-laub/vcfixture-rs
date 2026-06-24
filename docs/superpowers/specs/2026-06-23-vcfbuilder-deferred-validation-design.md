# VcfBuilder: deferred validation + typed `Field` declarations

**Date:** 2026-06-23
**Status:** Design approved, pending implementation plan
**Branch/worktree:** `builder-fallibility`

## Problem

`VcfBuilder`'s per-step methods are fallible: `info()`, `format()`, and `record()`
each return `Result<Self, BuildError>`. In test-fixture code — the crate's entire
reason to exist — this forces a `.unwrap()` (or `?`) after **every** chained call:

```rust
let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
    .info("AF", None, None, None).unwrap()
    .format("GT", None, None, None).unwrap()
    .format("DS", Some(Number::A), Some(Type::Float), None).unwrap()
    .record(/* ... */).unwrap()
    .build().unwrap();
```

Two problems:

1. **Ergonomics.** The `.unwrap()` noise dominates fixture code. The natural shape
   for a builder is infallible chaining with a single fallible terminal.
2. **A silent footgun.** `info(id, number, type_, description)` resolves via
   `make_def`: only `(Some(n), Some(t))` builds an explicit `FieldDef`; **every
   other combination** — including `(Some(Number::A), None)` — silently falls back
   to a reserved-registry lookup, ignoring the partially-supplied arguments.

Most VCF validation is irreducibly value-dependent and cross-field (SVLEN ↔ SVCLAIM
↔ SvType; cardinality depends on `n_alt` × `ploidy`; GT allele indices depend on
`alts.len()`), so it cannot move to compile time. But the *failure surface* can be
collapsed to one point, and the field-declaration footgun can be designed away.

## Goals

- **Primary:** infallible per-step chaining; `build()` is the only fallible point.
- **Secondary (cheap correctness only):** make impossible field-declaration combos
  unrepresentable, without contorting the fluent API.

Non-goals: typestate/ordering enforcement at compile time; aggregated multi-error
reporting; making every static check a compile error.

## Design

### 1. `VcfBuilder` becomes a pure spec accumulator

The builder stops resolving and validating eagerly. It stashes raw specs and does
all work in `build()`.

```rust
pub struct VcfBuilder {
    samples: Vec<String>,
    contigs: Vec<ContigDef>,
    version: VcfVersion,
    info_fields:   Vec<Field>,     // was IndexMap<String, FieldDef> (resolved eagerly)
    format_fields: Vec<Field>,     // was IndexMap<String, FieldDef>
    filter_defs: Vec<(String, String)>,
    alt_defs: IndexMap<String, String>,
    records: Vec<RecordSpec>,      // was Vec<Record> (converted + validated eagerly)
}
```

Method fallibility:

| method | before | after |
|---|---|---|
| `info` / `format` | `Result<Self, BuildError>` | **`Self`** |
| `record` | `Result<Self, BuildError>` | **`Self`** |
| `new` / `filter` / `alt` | `Self` | `Self` (unchanged) |
| `build` | `Result<Document, BuildError>` | `Result<Document, BuildError>` (unchanged) |
| `render` / `write` / `truth` | `Result<_, BuildError>` | `Result<_, BuildError>` (delegate to `build`) |

`RecordSpec` and all its setters are unchanged (already infallible).

### 2. The `Field` declaration sub-spec

Replaces the `info(id, Option<Number>, Option<Type>, Option<String>)` signature
with a small builder that mirrors the existing `RecordSpec` idiom.

```rust
pub struct Field {
    id: String,
    decl: Decl,
    description: Option<String>,
}

enum Decl {
    Reserved,                // resolve via the reserved registry at build()
    Typed(Number, Type),     // explicit FieldDef::new at build()
    Flag,                    // Number=Flag, Type=Flag
}

impl Field {
    pub fn reserved(id: impl Into<String>) -> Field;
    pub fn typed(id: impl Into<String>, number: Number, type_: Type) -> Field;
    pub fn flag(id: impl Into<String>) -> Field;
    pub fn description(self, d: impl Into<String>) -> Field;  // optional chain
}

impl From<&str>   for Field { /* => Field::reserved */ }
impl From<String> for Field { /* => Field::reserved */ }
```

`info` and `format` take `impl Into<Field>`, so the common reserved case stays terse:

```rust
.info("AF")                                   // reserved (via From<&str>)
.info(Field::reserved("AF"))                  // explicit reserved
.format(Field::typed("DS", Number::A, Type::Float))
.info(Field::typed("DP", Number::ONE, Type::Integer).description("read depth"))
.info(Field::flag("SOMATIC"))                 // Number=Flag, Type=Flag
```

Method signatures:

```rust
pub fn info(mut self,   field: impl Into<Field>) -> Self;
pub fn format(mut self, field: impl Into<Field>) -> Self;
```

**Footgun eliminated:** there is no longer any partial `(Some, None)` combination to
silently ignore. A declaration is exactly one of reserved / typed / flag.

### 3. Resolution and validation move into `build()`

`build()` performs, in order:

1. **Resolve field defs.** For each `Field` in `info_fields` / `format_fields`:
   - `Decl::Reserved` → `reserved(id, kind, self.version)`
   - `Decl::Typed(n, t)` → `FieldDef::new(id, n, t, description.unwrap_or(id), kind)`
   - `Decl::Flag` → `FieldDef::new(id, Number{Flag}, Type::Flag, description.unwrap_or(id), kind)`

   into `IndexMap<String, FieldDef>` for info and format (preserving today's
   structure for the rest of the pipeline). Any error here is returned bare.

2. **Process records.** For each `RecordSpec` at index `i`, run the existing
   record validation+conversion pipeline (the current body of `record()`):
   allele validation, GT parsing + allele-index range, sample-count match,
   cardinality, INFO/FORMAT cross-field SV rules, CN/SVLEN. Convert to `Record`.

3. **Auto-describe symbolic ALTs** and assemble `Document` (unchanged from today's
   `build()`).

**Fail-fast:** the first error wins. Record-processing errors are wrapped with
their index for context via one new `BuildError` variant:

```rust
#[error("record {index}: {source}")]
InRecord {
    index: usize,
    #[source]
    source: Box<BuildError>,
},
```

Field-resolution errors (step 1) are returned bare (no record context applies).

**Side effect — declaration-order independence.** Because all declarations are
collected before any record is validated, a record may reference a field declared
*after* its `record()` call. Today `GtNotDeclared` fires if GT is not declared *at
the time `record()` runs*; under this design it fires only if GT is never declared.
This is strictly more permissive and removes an ordering gotcha.

### 4. Residual runtime checks (scope boundary)

These remain `BuildError`s surfaced at `build()`, **not** compile errors. Making
them static would require splitting `Field` into distinct `InfoField` / `FormatField`
types, contorting the fluent API — rejected per the "ergonomics first" goal:

- `FlagNotInfo` — `Field::flag(...)` (or a `Type::Flag` typed decl) passed to
  `format()`; already enforced inside `FieldDef::new`.
- `BadFieldId` — VCF key regex.
- `FieldTooNew` — reserved field gated by `version`.

## Error handling summary

- One new variant: `BuildError::InRecord { index, source }`.
- All existing variants retained; the per-record ones now arrive wrapped in
  `InRecord` when produced during `build()`.
- `build()` is fail-fast (first error). Aggregated reporting is an explicit
  non-goal; the wrapping variant leaves the door open to revisit later.

## Migration impact

- **In-crate tests** (`src/build.rs`): drop `.unwrap()` from `info`/`format`/`record`
  chains; keep one at `build()`. Error-assertion tests change from matching on the
  `record(...)` result to matching on the `build()` result — per-record errors now
  match `BuildError::InRecord { source, .. }` wrapping the original variant.
- **Doctest** (`src/lib.rs`): update the example to the new infallible chaining and
  the `Field` / `"AF"` declaration form.
- **`strategies` / proptest:** no impact — the strategy layer does not use
  `VcfBuilder`.
- **Public API:** `Field` is newly exported from the crate root alongside
  `RecordSpec`. `info`/`format` signatures change (breaking, pre-1.0, acceptable).
- **Memory / prior spec:** `project-overview` memory and the original design spec
  record "Result + `?` chaining builder" as the chosen idiom; update both to reflect
  the infallible-chaining + single-terminal decision.

## Testing strategy

Reuse the existing `build.rs` error-path tests (they already exercise every
per-record `BuildError`), retargeted to `build()` + `InRecord`. Add:

- `Field` resolution: reserved, typed (with and without `.description()`), flag.
- `From<&str>` shorthand: `.info("AF")` ≡ `.info(Field::reserved("AF"))`.
- flag-on-format → `FlagNotInfo` at `build()`.
- Declaration-order independence: declare a FORMAT field *after* the `record()`
  call that uses it; `build()` succeeds.
- Happy-path chain has zero `.unwrap()` before `build()`.

## Open questions

None blocking. Decisions locked during design:

- **(a)** Fail-fast (first error), not collect-all. Easy to revisit via the
  `InRecord` wrapper.
- **(b)** flag-on-format stays a runtime `BuildError`, not a compile error.
- **(c)** `Field::typed` uses existing `Number` constructors (`Number::A`,
  `Number::ONE`, …); a `Number::fixed(n)` helper may be added if call sites want it.
