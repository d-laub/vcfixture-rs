# Bulk-generation follow-ups (#6–#10) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close GitHub issues #6–#10 (IDE/edition hygiene, profile-schema
correctness, `Size::Target` performance, and `fit_profile.py` genome-scale
memory) as four independent workstreams, ending with a re-fitted genome-wide
somatic profile.

**Architecture:** Four workstreams over mostly-disjoint files. **A** (rand 0.9 +
`gen`→`generate` rename) and **B** (move `ploidy` to `Dialed`, add
`provenance.supplied`, parse-time payload/ploidy check, fresh-fit CI gate) are
independent and land first. **C** (two-point byte-target calibration that
promotes the measured file instead of regenerating; per-record allocation
fixes; `n_variants` split weight) and **D** (`scan_csv` `skip_lines`; gap
shift+mask; Ti/Tv direct comparisons; sequential `collect`; single-bin
histogram fix; somatic re-fit) build on B's settled schema and land second.

**Tech Stack:** Rust (edition 2021, `bulk`/`cli` cargo features), noodles-vcf
0.83, rayon, rand 0.9, tempfile; Python 3.11 + polars 1.42 + numpy + plink2 +
bcftools under `pixi run -e fit`.

## Global Constraints

- **Rust edition stays 2021.** Do NOT migrate to edition 2024 (out of scope);
  the `gen`→`generate` rename is what unblocks rust-analyzer.
- **`bulk` is a cargo feature, not default.** All Rust build/test/lint commands
  in this plan use `--all-features` so the `bulk`/`cli` code is compiled.
- **Every Rust gate is:** `cargo test --all-features`,
  `cargo clippy --all-features -- -D warnings`, `cargo fmt --check` — all green.
- **Every fit-script gate is:** `pixi run -e fit test-fit` (unit) and, where a
  task says so, `pixi run -e fit test-fidelity`.
- **prek hooks must be installed** (`pixi run prek-install`) before the first
  commit; they run `cargo fmt`/`clippy`/`check` and commitizen on commit.
- **Commit messages** use Conventional Commits (commitizen `commit-msg` hook
  enforces it): `type(scope): summary`. End every commit body with
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.
- **`fitted` must contain only values measured from data.** `ploidy` moves to
  `Dialed`; every non-measured field must be named in `provenance.supplied`.
- **Determinism:** generation output is a pure function of `(seed, block_idx)`.
  Tasks that change the RNG draw pattern (C2) must re-baseline any golden
  output explicitly and keep the self-consistency determinism tests passing.
- **Real source data** for D's re-fit and benchmarks (present on this node):
  somatic `/carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt`
  (`.pgen/.pvar/.psam`, 348,259,675 rows); 32 GB SLURM allocation.

**Sequencing:** Do **A ∥ B first, then C ∥ D**. Within D, task D6 (re-fit) runs
last and depends on B (schema) + D1–D5 (memory).

---

## File map

| File | Responsibility | Workstreams |
|---|---|---|
| `Cargo.toml` | rand 0.9 / rand_chacha 0.9 version bump; `validate-profile` bin | A, B |
| `src/reference.rs` | rand call-site migration | A |
| `src/bulk/sample.rs` | rand call-site migration | A |
| `src/bulk/gen.rs` → `src/bulk/generate.rs` | rand migration; module rename; GT-alloc fix (C1); alt rejection-sample (C2) | A, C |
| `src/bulk/mod.rs` | `pub mod generate`; rand migration; ploidy reader moves (B1); split weight (C3); calibration+promote (C4) | A, B, C |
| `src/bulk/profile.rs` | `ploidy` Fitted→Dialed; `provenance.supplied`; parse-time payload check | B |
| `src/bin/validate_profile.rs` | CI validator entry point | B |
| `profiles/germline-1kgp.json`, `germline-1kgp-unphased.json`, `somatic-gdc.json` | schema move; somatic re-fit | B, D |
| `tests/bulk.rs` | inline profile JSON schema move; payload-guard test moves to validate; overshoot bound unchanged | B, C |
| `scripts/fit/fit_profile.py` | `read_pvar` scan; `_gap_bins_lazy`; `_titv_lazy`; `compute_pvar_stats`; `_bucket_index_expr`; `build_profile` schema; `supplied`; enum-drift guard | B, D |
| `scripts/fit/test_fit_profile.py` | tests for the above | B, D |
| `.github/workflows/*.yml` | fresh-fit validation gate | B |

---

# Workstream A — hygiene (#6)

## Task A1: Upgrade rand 0.8 → 0.9 and migrate all call sites

**Files:**
- Modify: `Cargo.toml` (deps `rand`, `rand_chacha`)
- Modify: `src/reference.rs:204`
- Modify: `src/bulk/sample.rs:67,70,169,183,196,203`
- Modify: `src/bulk/gen.rs:106,122,171,187,188,190,335,336,337`
- Modify: `src/bulk/mod.rs:465`
- Test: existing `cargo test --all-features` suite (re-baseline if needed)

**Interfaces:**
- Consumes: nothing new.
- Produces: no signature changes. rand 0.9 method names (`random`,
  `random_range`) replace 0.8 (`gen`, `gen_range`) crate-wide. `block_rng`
  (`src/bulk/gen.rs`) keeps its signature `fn block_rng(seed: u64, block_idx: u64) -> ChaCha8Rng`.

- [ ] **Step 1: Capture the current generation baseline BEFORE any change**

Run and save output — this is the oracle for whether output shifts:

```bash
cargo test --all-features 2>&1 | tee /tmp/pre_rand_tests.txt | tail -20
```

Expected: all pass. Keep `/tmp/pre_rand_tests.txt`.

- [ ] **Step 2: Bump versions in `Cargo.toml`**

```toml
rand = "0.9"
rand_chacha = "0.9"
```

- [ ] **Step 3: Migrate `src/reference.rs`**

Line 204, change `gen_range` → `random_range`:

```rust
            .map(|_| BASES[self.rng.random_range(0..4)])
```

- [ ] **Step 4: Migrate `src/bulk/sample.rs`**

```rust
// line 67:
        let u: f64 = rng.random();
// line 70:
        rng.random_range(lo..hi)
// line 169:
        let u: f64 = rng.random();
// line 183:
        b"ACGT"[rng.random_range(0..4)]
// line 196:
        if rng.random::<f64>() < self.ti_frac {
// line 203:
            transversions[rng.random_range(0..2)]
```

- [ ] **Step 5: Migrate `src/bulk/gen.rs`**

```rust
// line 106:
        let sample_missing = rng.random::<f64>() < fitted.missing_rate;
// line 122:
        let j = rng.random_range(i..idx.len());
// line 171:
            let len = rng.random_range(2..=3);
// line 187:
            let ref_len = rng.random_range(2..=4);
// line 188:
            let mut alt_len = rng.random_range(2..=4);
// line 190:
                alt_len = rng.random_range(2..=4);
// lines 335-337 (in tests):
        let xa: u64 = a.random();
        let xb: u64 = b.random();
        let xc: u64 = c.random();
```

- [ ] **Step 6: Migrate `src/bulk/mod.rs`**

```rust
// line 465:
                        let phased = rng.random::<f64>() < fitted.phased_rate;
```

- [ ] **Step 7: Compile and fix any residual 0.9 breakage**

```bash
cargo build --all-features 2>&1 | tail -30
```

Expected: clean. If the compiler flags removed 0.8 imports (e.g. a `use
rand::distributions::Standard` — none exist in this tree per grep, but verify),
follow its suggestions. `rng.random::<f64>()` requires `use rand::Rng;` which is
already imported at every call site.

- [ ] **Step 8: Run the full suite and diff against the baseline**

```bash
cargo test --all-features 2>&1 | tee /tmp/post_rand_tests.txt | tail -20
diff /tmp/pre_rand_tests.txt /tmp/post_rand_tests.txt || true
```

Expected outcomes, and what each means:
- **All green, no failures:** `ChaCha8Rng::seed_from_u64` and the uniform
  samplers produced identical streams; nothing to re-baseline. Proceed.
- **A snapshot/golden test fails** (e.g. an `insta` snapshot in `reference.rs`
  or a golden hash in `tests/bulk.rs`): rand 0.9's integer-uniform sampling
  shifted the stream. This is expected and acceptable — re-baseline in Step 9.

- [ ] **Step 9: Re-baseline shifted expectations (only if Step 8 showed failures)**

For `insta` snapshots: `cargo insta review` (accept the new values) or
`INSTA_UPDATE=always cargo test --all-features`. For hand-coded golden
constants: update them to the new observed values. In each changed test, add a
one-line comment: `// re-baselined for rand 0.9 (uniform-int stream change)`.
Do NOT re-baseline a test that asserts a *structural invariant* (shape,
alphabet, determinism self-consistency) — those must still pass unchanged; if
one of those fails, the migration has a real bug, stop and debug.

- [ ] **Step 10: Lint and commit**

```bash
cargo clippy --all-features -- -D warnings && cargo fmt --check
git add Cargo.toml Cargo.lock src/ tests/
git commit -m "build(deps): upgrade rand 0.8 -> 0.9, migrate call sites

gen()/gen_range() -> random()/random_range(). Removes the rand-0.8 half of
the edition-2024 gen-keyword collision (#6). Any re-baselined golden output
is noted inline."
```

---

## Task A2: Rename module `gen` → `generate`

**Files:**
- Rename: `src/bulk/gen.rs` → `src/bulk/generate.rs`
- Modify: `src/bulk/mod.rs:12` (`pub mod gen;`), `:35` (`use gen::{...}`), and
  doc comments referencing `[`gen`]` / `gen::` (`:8`, `:783`, `:808`, and any
  others surfaced by grep)

**Interfaces:**
- Consumes: nothing new.
- Produces: the module is `crate::bulk::generate`; public items
  (`GenRecord`, `block_rng`, `gen_record`, `to_record_buf`) keep their names —
  only the module path changes. (`gen_record` is a function name, not the
  reserved word `gen`, so it stays; the reserved word only bites `mod gen` and
  `rng.gen()`, the latter already gone in A1.)

- [ ] **Step 1: Rename the file with git**

```bash
git mv src/bulk/gen.rs src/bulk/generate.rs
```

- [ ] **Step 2: Update the module declaration and import in `src/bulk/mod.rs`**

```rust
// line 12:
pub mod generate;
// line 35:
use generate::{block_rng, gen_record, to_record_buf, GenRecord};
```

- [ ] **Step 3: Update every remaining `gen` reference**

Find them all:

```bash
rg -n '\bgen\b|gen::|`gen`|crate::bulk::gen' src/ tests/ examples/ docs/book
```

Update doc-comment references (`src/bulk/mod.rs:8` "the record generator
([`gen`])" → "([`generate`])"; `:783` and `:808` "`gen::to_record_buf`" /
"`gen::SampleStats`" → "`generate::to_record_buf`" / "`generate::SampleStats`").
In `src/bulk/generate.rs`, update any self-referential doc mention of the module
name. Leave `gen_record`, `GenRecord`, and `gen_site` identifiers untouched —
they are not the reserved word.

- [ ] **Step 4: Verify no bare `gen` token survives in `src/`**

```bash
rg -n '\bgen\b' src/            # module/keyword only; expect NO hits
rg -n 'mod gen|rng\.gen' src/   # expect NO hits
```

Expected: no matches for `mod gen` or `rng.gen`. (`gen_record`/`GenRecord`/
`gen_site` are fine and will not match `\bgen\b`.)

- [ ] **Step 5: Build, test, lint**

```bash
cargo build --all-features && cargo test --all-features \
  && cargo clippy --all-features -- -D warnings && cargo fmt --check
```

Expected: all green.

- [ ] **Step 6: Manual rust-analyzer check**

Open `src/bulk/generate.rs` and `src/bulk/mod.rs` in the editor; confirm no
"Syntax Error" diagnostics and that `profile.rs` is no longer reported as "not
included in any crates." (This is the acceptance criterion #6 exists for.)

- [ ] **Step 7: Commit**

```bash
git add -A src/ tests/ examples/ docs/
git commit -m "refactor(bulk): rename module gen -> generate (#6)

'gen' is reserved in edition 2024; rust-analyzer resolved src/bulk/ under
2024 rules and broke IDE services across the module. Combined with the
rand 0.9 upgrade this removes every 'gen' token. Closes #6."
```

---

# Workstream B — schema (#9, #10 schema items)

## Task B1: Move `ploidy` from `Fitted` to `Dialed`

**Files:**
- Modify: `src/bulk/profile.rs` (struct `Fitted` `:49`, struct `Dialed` `:80`,
  `validate` `:119`, test `:221`)
- Modify: `src/bulk/mod.rs:308` (payload guard reader), `:439` (generate reader)
- Modify: `profiles/germline-1kgp.json`, `profiles/germline-1kgp-unphased.json`,
  `profiles/somatic-gdc.json`
- Modify: `tests/bulk.rs` inline JSON (`:409`, `:520`)
- Modify: `scripts/fit/fit_profile.py:452` (`build_profile`)
- Test: `src/bulk/profile.rs` tests, `cargo test --all-features`

**Interfaces:**
- Consumes: nothing new.
- Produces: `Dialed` gains `pub ploidy: u8`; `Fitted` loses `ploidy`. Readers
  use `self.profile.dialed.ploidy` (in `BulkSpec`) / `dialed.ploidy`.

- [ ] **Step 1: Write a failing test for the new location**

In `src/bulk/profile.rs` tests module, add:

```rust
    #[test]
    fn ploidy_lives_in_dialed_not_fitted() {
        let p = Profile::builtin("germline-1kgp").unwrap();
        assert_eq!(p.dialed.ploidy, 2);
    }
```

- [ ] **Step 2: Run it — expect a compile error**

```bash
cargo test --all-features profile:: 2>&1 | tail -15
```

Expected: FAIL to compile — `no field ploidy on Dialed`.

- [ ] **Step 3: Move the struct field**

`src/bulk/profile.rs`: remove `pub ploidy: u8,` from `Fitted` (line 49) and add
it to `Dialed`:

```rust
pub struct Dialed {
    pub payload: Payload,
    pub ploidy: u8,
}
```

- [ ] **Step 4: Update `validate` to read `dialed.ploidy`**

`src/bulk/profile.rs:119`:

```rust
        if self.dialed.ploidy == 0 {
            return Err(BulkError::Invalid("ploidy must be >= 1".into()));
        }
```

- [ ] **Step 5: Update the two readers in `src/bulk/mod.rs`**

Line 308 (payload guard) — read `self.profile.dialed.ploidy`. Because this
block already borrows `fitted`, bind a local first for clarity:

```rust
        let keys = payload_keys(&self.payload);
        let ploidy = self.profile.dialed.ploidy;
        if (keys.contains(&"PL") || keys.contains(&"AD")) && ploidy != 2 {
            return Err(BulkError::Invalid(format!(
                "payload {:?} declares PL and/or AD, which are hard-coded for \
                 diploid (ploidy 2) genotype calls, but the profile's ploidy is {}",
                self.payload, ploidy
            )));
        }
```

Line 439 (`generate_contig`):

```rust
        let ploidy = self.profile.dialed.ploidy;
```

- [ ] **Step 6: Update the existing `profile.rs` test at line 221**

```rust
        assert_eq!(p.dialed.ploidy, 2);
```

- [ ] **Step 7: Move the key in all three committed profiles**

For each of `profiles/germline-1kgp.json`, `germline-1kgp-unphased.json`,
`somatic-gdc.json`: delete `"ploidy": N` from the `"fitted"` object and add it
to `"dialed"`. Use this jq (adjust value per file — all three are currently `2`):

```bash
for f in profiles/germline-1kgp.json profiles/germline-1kgp-unphased.json profiles/somatic-gdc.json; do
  jq '.dialed.ploidy = .fitted.ploidy | del(.fitted.ploidy)' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done
```

- [ ] **Step 8: Move the key in the two inline test profiles**

`tests/bulk.rs`: in the JSON string constants, move `"ploidy": 2` (line ~409)
and `"ploidy": 3` (line ~520) out of `"fitted"` and into `"dialed"`. Example
for the first:

```rust
  "dialed": { "payload": "gt-only", "ploidy": 2 }
```

and for the ploidy-3 fixture:

```rust
  "dialed": { "payload": "gt-only", "ploidy": 3 }
```

(delete the `"ploidy": N` line from each `"fitted"` block).

- [ ] **Step 9: Update `build_profile` in the fit script**

`scripts/fit/fit_profile.py`: remove `"ploidy": ploidy,` from the `fitted` dict
(line 452) and add it to the `dialed` dict (line 464):

```python
        "dialed": {"payload": payload, "ploidy": ploidy},
```

- [ ] **Step 10: Run all Rust tests and the fit unit tests**

```bash
cargo test --all-features 2>&1 | tail -15
pixi run -e fit test-fit 2>&1 | tail -15
```

Expected: all green (the new `ploidy_lives_in_dialed_not_fitted` passes; any
fit-script test asserting `fitted.ploidy` now asserts `dialed.ploidy` — update
it if present).

- [ ] **Step 11: Commit**

```bash
git add -A src/ tests/ profiles/ scripts/fit/fit_profile.py
git commit -m "refactor(bulk): move ploidy from fitted to dialed (#9)

ploidy is a generation choice, never derived from source data, so it
belongs in Dialed alongside payload -- not in Fitted, whose invariant is
'measured from data only'. Closes #9's core ask."
```

---

## Task B2: Add `provenance.supplied` and scope validation to it

**Files:**
- Modify: `src/bulk/profile.rs` (struct `Provenance` `:31`)
- Modify: `scripts/fit/fit_profile.py` (`build_profile` `:408`, both callers
  `:1006`, `:1042`)
- Modify: `profiles/*.json` (add the array)
- Test: `src/bulk/profile.rs`, `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: nothing new.
- Produces: `Provenance` gains `pub supplied: Vec<String>` (serde default
  `Vec::new()` so old JSON still parses). `build_profile` gains a keyword arg
  `supplied: list[str]`.

- [ ] **Step 1: Write a failing Rust test**

In `src/bulk/profile.rs` tests:

```rust
    #[test]
    fn sites_only_profile_marks_phased_rate_supplied() {
        // germline-1kgp-unphased is fitted from a sites-only VCF, so
        // phased_rate and n_samples are supplied, not measured.
        let p = Profile::builtin("germline-1kgp-unphased").unwrap();
        assert!(p.provenance.supplied.contains(&"ploidy".to_string()));
        assert!(p.provenance.supplied.contains(&"phased_rate".to_string()));
    }
```

- [ ] **Step 2: Run it — expect compile failure**

```bash
cargo test --all-features profile:: 2>&1 | tail -12
```

Expected: FAIL — `no field supplied on Provenance`.

- [ ] **Step 3: Add the field with a serde default**

`src/bulk/profile.rs`, struct `Provenance`:

```rust
pub struct Provenance {
    pub source: String,
    pub n_samples_source: usize,
    pub n_variants_source: u64,
    pub fitted_on: String,
    pub fit_tool_version: String,
    #[serde(default)]
    pub supplied: Vec<String>,
}
```

- [ ] **Step 4: Add `supplied` to `build_profile` and both callers**

`scripts/fit/fit_profile.py`, `build_profile` signature — add `supplied: list[str]`
(no default; force callers to be explicit), and emit it in the provenance dict:

```python
        "provenance": {
            "source": source,
            "n_samples_source": n_samples,
            "n_variants_source": n_variants_source,
            "fitted_on": _dt.date.today().isoformat(),
            "fit_tool_version": __version__,
            "supplied": sorted(supplied),
        },
```

`_fit_from_pgen` caller (phased_rate IS measured here):

```python
        ploidy=args.ploidy,
        payload=args.payload,
        supplied=["ploidy"],
    )
```

`_fit_from_sites_vcf` caller (phased_rate and n_samples supplied):

```python
        ploidy=args.ploidy,
        payload=args.payload,
        supplied=["ploidy", "phased_rate", "n_samples"],
    )
```

- [ ] **Step 5: Add `supplied` to the three committed profiles**

```bash
# germline-1kgp is fitted from --pgen:
jq '.provenance.supplied = ["ploidy"]' profiles/germline-1kgp.json > t && mv t profiles/germline-1kgp.json
# germline-1kgp-unphased is fitted from --sites-vcf:
jq '.provenance.supplied = ["n_samples","phased_rate","ploidy"]' profiles/germline-1kgp-unphased.json > t && mv t profiles/germline-1kgp-unphased.json
# somatic-gdc is fitted from --pgen (will be re-emitted in D6, but set it now so it loads):
jq '.provenance.supplied = ["ploidy"]' profiles/somatic-gdc.json > t && mv t profiles/somatic-gdc.json
```

- [ ] **Step 6: Write a fit-script test asserting the two paths' `supplied`**

In `scripts/fit/test_fit_profile.py`:

```python
def test_build_profile_records_supplied_fields():
    prof = build_profile(
        name="t", source="x", n_samples=10,
        contigs=[{"id": "chr1", "n_variants": 100, "density_per_kb": 40.0}],
        gap_dist={"edges": [1.0, 2.0], "weights": [1.0]},
        sfs={"edges": [1.0, 2.0], "weights": [1.0]},
        indel_length={"edges": [1.0, 2.0], "weights": [1.0]},
        class_counts={n: 1 for n in CLASS_NAMES},
        titv=2.0, multiallelic_rate=0.1, missing_rate=0.0,
        phased_rate=1.0, ploidy=2, supplied=["ploidy", "phased_rate"],
    )
    assert prof["provenance"]["supplied"] == ["phased_rate", "ploidy"]
    assert "ploidy" not in prof["fitted"]
    assert prof["dialed"]["ploidy"] == 2
```

- [ ] **Step 7: Run both test suites**

```bash
cargo test --all-features profile:: 2>&1 | tail -12
pixi run -e fit test-fit 2>&1 | tail -12
```

Expected: green.

- [ ] **Step 8: Commit**

```bash
git add -A src/ profiles/ scripts/fit/
git commit -m "feat(bulk): add provenance.supplied naming non-measured fields (#9)

phased_rate on the sites-vcf path is hand-supplied but cannot move to
Dialed (it IS measured on --pgen). provenance.supplied makes every
non-measured field (ploidy always; phased_rate/n_samples on sites-vcf)
auditable from the JSON alone."
```

---

## Task B3: Reject `ploidy != 2` with AD/PL payload at parse time

**Files:**
- Modify: `src/bulk/profile.rs` (`validate`, add `Payload::needs_diploid`)
- Modify: `src/bulk/mod.rs:307-314` (remove the runtime guard)
- Modify: `tests/bulk.rs` (`:461` test moves from `write()` to `validate()`)

**Interfaces:**
- Consumes: `Dialed.payload`, `Dialed.ploidy` (both now in `Dialed` after B1).
- Produces: `Profile::validate` returns `Err` when `dialed.payload` uses AD/PL
  and `dialed.ploidy != 2`. The runtime guard in `write()` is deleted. New
  `pub(crate) fn Payload::needs_diploid(&self) -> bool`.

- [ ] **Step 1: Add a `needs_diploid` accessor on `Payload`**

`src/bulk/profile.rs`:

```rust
impl Payload {
    /// True if this preset emits AD or PL, both hard-coded diploid in
    /// `generate::SampleStats`.
    pub(crate) fn needs_diploid(&self) -> bool {
        matches!(self, Payload::Gatk | Payload::Mutect2)
    }
}
```

(Gatk keys include AD+PL; Mutect2 includes AD; GtOnly/GtVaf include neither —
matches `to_record_buf`'s `key_names`.)

- [ ] **Step 2: Write the failing test — validate rejects ploidy 3 + Gatk**

In `src/bulk/profile.rs` tests:

```rust
    #[test]
    fn validate_rejects_diploid_payload_with_nondiploid_ploidy() {
        let mut p = Profile::builtin("germline-1kgp").unwrap();
        p.dialed.payload = Payload::Gatk;
        p.dialed.ploidy = 3;
        let err = p.validate().unwrap_err();
        assert!(
            format!("{err}").contains("diploid"),
            "expected a diploid-ploidy rejection, got: {err:?}"
        );
    }
```

- [ ] **Step 3: Run it — expect failure (validate accepts it today)**

```bash
cargo test --all-features validate_rejects_diploid 2>&1 | tail -12
```

Expected: FAIL — `validate` currently returns Ok, so `unwrap_err` panics.

- [ ] **Step 4: Add the check to `validate`**

`src/bulk/profile.rs`, in `validate`, after the `ploidy == 0` check:

```rust
        if self.dialed.payload.needs_diploid() && self.dialed.ploidy != 2 {
            return Err(BulkError::Invalid(format!(
                "payload {:?} emits AD and/or PL, which are hard-coded for \
                 diploid (ploidy 2) calls, but ploidy is {}",
                self.dialed.payload, self.dialed.ploidy
            )));
        }
```

- [ ] **Step 5: Remove the now-redundant runtime guard in `write()`**

Delete `src/bulk/mod.rs` lines 307-314 (the `let keys = payload_keys(...)` +
`if (keys.contains(&"PL") ...)` block). `validate()` is already called at the
top of `write()`, so the constraint is enforced earlier and this is dead. If
`payload_keys` becomes unused as a result, leave it (it is still used by
`build_header` at `:394`).

- [ ] **Step 6: Move the `tests/bulk.rs` guard test to assert via `validate`**

The existing test around `tests/bulk.rs:461` builds a ploidy-3 profile + Gatk
payload and expects `write()` (or `BulkSpec` construction) to error. Change it
to assert `Profile::validate()` errors directly (parse-time), keeping the
ploidy-3 fixture. Keep any assertion that GtOnly/GtVaf with ploidy 3 is still
*accepted* by validate.

- [ ] **Step 7: Test, lint**

```bash
cargo test --all-features 2>&1 | tail -15
cargo clippy --all-features -- -D warnings
```

Expected: green.

- [ ] **Step 8: Commit**

```bash
git add -A src/ tests/
git commit -m "fix(bulk): reject non-diploid AD/PL payloads at validate time (#10)

AD/PL are hard-coded diploid in SampleStats; validate now rejects
ploidy != 2 with a Gatk/Mutect2 payload instead of write() catching it at
runtime, making the invalid state unrepresentable at construction."
```

---

## Task B4: CI gate — validate every freshly-written profile through Rust

**Files:**
- Create: `src/bin/validate_profile.rs`
- Modify: `Cargo.toml` (new bin target)
- Modify: `scripts/fit/fit_profile.py` (post-write self-check in `main`)
- Modify: `.github/workflows/ci.yml` (or the existing fit workflow)
- Test: `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: `Profile::from_json` + `Profile::validate` (already `pub`,
  re-exported from `vcfixture::bulk`).
- Produces: a `validate-profile` binary (behind `bulk`) exiting non-zero on any
  invalid profile at an arbitrary path. Writing a bad profile fails the fit run.

- [ ] **Step 1: Add the bin target to `Cargo.toml`**

```toml
[[bin]]
name = "validate-profile"
path = "src/bin/validate_profile.rs"
required-features = ["bulk"]
```

- [ ] **Step 2: Create `src/bin/validate_profile.rs`**

```rust
//! Validate a profile JSON file through the same `Profile::validate` the
//! embedded profiles pass, so a freshly-fitted profile fails here in CI
//! rather than later at `include_str!` time.
use std::process::ExitCode;
use vcfixture::bulk::Profile;

fn main() -> ExitCode {
    let path = match std::env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("usage: validate-profile <path.json>");
            return ExitCode::FAILURE;
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("{path}: read error: {e}");
            return ExitCode::FAILURE;
        }
    };
    match Profile::from_json(&text).and_then(|p| p.validate().map(|_| p)) {
        Ok(p) => {
            println!("{path}: OK ({} contigs)", p.fitted.contigs.len());
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("{path}: INVALID: {e}");
            ExitCode::FAILURE
        }
    }
}
```

- [ ] **Step 3: Build it and validate the committed profiles**

```bash
cargo build --features bulk --bin validate-profile
for f in profiles/*.json; do cargo run -q --features bulk --bin validate-profile -- "$f"; done
```

Expected: each prints `OK`. (Confirms `Profile`, `from_json`, `validate`, and
`fitted.contigs` are reachable from `vcfixture::bulk`.)

- [ ] **Step 4: Write a test that a broken profile is rejected**

In `scripts/fit/test_fit_profile.py`:

```python
import json, subprocess

def test_validate_profile_binary_rejects_nan(tmp_path):
    prof = build_profile(
        name="bad", source="x", n_samples=10,
        contigs=[{"id": "chr1", "n_variants": 100, "density_per_kb": 40.0}],
        gap_dist={"edges": [1.0, 2.0], "weights": [1.0]},
        sfs={"edges": [1.0, 2.0], "weights": [1.0]},
        indel_length={"edges": [1.0, 2.0], "weights": [1.0]},
        class_counts={n: 1 for n in CLASS_NAMES},
        titv=2.0, multiallelic_rate=0.1, missing_rate=0.0,
        phased_rate=1.0, ploidy=2, supplied=["ploidy"],
    )
    prof["fitted"]["variant_classes"]["snp"] = float("nan")  # poison
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(prof))
    r = subprocess.run(
        ["cargo", "run", "--quiet", "--features", "bulk",
         "--bin", "validate-profile", "--", str(p)],
        capture_output=True, text=True,
    )
    assert r.returncode != 0, r.stdout + r.stderr
```

Mark it with `@pytest.mark.skipif(shutil.which("cargo") is None, ...)` so a
Python-only sandbox skips it; it is primarily a CI check.

- [ ] **Step 5: Add the post-write self-check to `fit_profile.py`**

After `main()` writes the JSON (`fit_profile.py:1152`,
`args.out.write_text(...)`), run the validator so a bad fit fails immediately.
Prefer the binary (single source of truth), warn-and-continue if it cannot be
built:

```python
    args.out.write_text(json.dumps(profile, indent=2) + "\n")
    print(f"wrote {args.out}")
    _validate_with_rust(args.out)
```

`_validate_with_rust(path)` runs `cargo run -q --features bulk --bin
validate-profile -- <path>`; on non-zero exit it raises `SystemExit(1)` with the
validator's stderr; if `cargo` is absent it prints a warning and returns. Put
the exact behavior in the function docstring.

- [ ] **Step 6: Wire the gate into CI**

In the workflow that runs `test-fit` (or a new job), after building the crate:

```yaml
      - name: Validate committed profiles through Rust
        run: |
          cargo build --features bulk --bin validate-profile
          for f in profiles/*.json; do
            cargo run --quiet --features bulk --bin validate-profile -- "$f"
          done
```

- [ ] **Step 7: Run locally**

```bash
for f in profiles/*.json; do cargo run -q --features bulk --bin validate-profile -- "$f"; done
pixi run -e fit test-fit 2>&1 | tail -12
```

Expected: all three profiles `OK`; fit tests pass.

- [ ] **Step 8: Commit**

```bash
git add -A Cargo.toml src/bin/validate_profile.rs scripts/fit/ .github/
git commit -m "ci(bulk): gate fresh fits through the Rust profile validator (#10)

A new profile emitting titv<=0, ploidy 0, or a NaN bin previously wrote
happily and only failed at include_str! time. validate-profile runs
Profile::validate over any path; CI runs it over profiles/*.json and the
fit script self-checks after writing."
```

---

## Task B5: Enforce `CLASS_NAMES` / `--payload` choices against the Rust enums

**Files:**
- Modify: `scripts/fit/fit_profile.py` (`CLASS_NAMES` `:69`, `--payload` choices)
- Test: `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: the Rust `ClassMix` field names and `Payload` variants as truth.
- Produces: `_PAYLOAD_CHOICES`/`_payload_choices()` module constant; a test that
  fails loudly if `CLASS_NAMES` or the `--payload` choices drift from the Rust
  enums, replacing the current "must match exactly" comment.

- [ ] **Step 1: Expose the `--payload` choices as a named list**

In `fit_profile.py`, replace the inline `choices=[...]` on `--payload` with a
module-level constant the test can import:

```python
_PAYLOAD_CHOICES = ("gt-only", "gt-vaf", "gatk", "mutect2")

def _payload_choices() -> tuple[str, ...]:
    return _PAYLOAD_CHOICES
```

and `... "--payload", choices=_PAYLOAD_CHOICES, ...`.

- [ ] **Step 2: Write the drift tests**

`ClassMix` fields are `pub snp: f64`, etc.; `Payload` variants
`GtOnly`/`GtVaf`/`Gatk`/`Mutect2` with `#[serde(rename_all = "kebab-case")]`.
Add to `scripts/fit/test_fit_profile.py`:

```python
import re
from pathlib import Path

_PROFILE_RS = Path(__file__).resolve().parents[2] / "src" / "bulk" / "profile.rs"

def _rust_classmix_fields() -> list[str]:
    src = _PROFILE_RS.read_text()
    block = re.search(r"pub struct ClassMix \{(.*?)\}", src, re.S).group(1)
    return re.findall(r"pub (\w+): f64", block)

def _rust_payload_variants() -> list[str]:
    src = _PROFILE_RS.read_text()
    block = re.search(r"pub enum Payload \{(.*?)\}", src, re.S).group(1)
    camel = re.findall(r"\b([A-Z]\w+),", block)
    return [re.sub(r"(?<!^)(?=[A-Z])", "-", c).lower() for c in camel]

def test_class_names_match_rust_classmix():
    assert list(CLASS_NAMES) == _rust_classmix_fields()

def test_payload_choices_match_rust_enum():
    assert sorted(_payload_choices()) == sorted(_rust_payload_variants())
```

- [ ] **Step 3: Run the tests — expect green (no drift today)**

```bash
pixi run -e fit test-fit -k "class_names_match or payload_choices_match" 2>&1 | tail -12
```

Expected: PASS. (If it fails, the enums already drifted — reconcile the Python
constants to the Rust truth.)

- [ ] **Step 4: Replace the "must match exactly" comment**

Delete the hand-mirror comment near `CLASS_NAMES` and the `--payload` choices;
replace with `# Enforced against the Rust enums by test_fit_profile.py`.

- [ ] **Step 5: Commit**

```bash
git add -A scripts/fit/
git commit -m "test(fit): enforce CLASS_NAMES/payload against Rust enums (#10)

Replaces the unenforced 'must match exactly' comment with tests that parse
ClassMix/Payload out of profile.rs and fail on drift."
```

---

# Workstream C — generator perf (#8, #10 perf items)

> **Note:** C tasks edit `src/bulk/generate.rs` (renamed from `gen.rs` in A2).
> If A has not landed yet, the file is still `src/bulk/gen.rs`; adjust paths.

## Task C1: Reuse a `String` buffer for GT instead of `Vec<String>` + `join`

**Files:**
- Modify: `src/bulk/generate.rs` (`SampleStats::new` `:220`)
- Test: `src/bulk/generate.rs` tests

**Interfaces:**
- Consumes: nothing new.
- Produces: identical GT strings (`"0|1"`, `"./."`, etc.); no output change,
  purely fewer allocations. `SampleStats.gt` stays `String`.

- [ ] **Step 1: Write a test pinning the exact GT strings**

In `src/bulk/generate.rs` tests (in-module, so it can call the private ctor):

```rust
    #[test]
    fn sample_stats_gt_string_is_unchanged_by_buffer_reuse() {
        assert_eq!(SampleStats::new(&[0, 1], true).gt, "0|1");
        assert_eq!(SampleStats::new(&[1, 1], false).gt, "1/1");
        assert_eq!(SampleStats::new(&[-1, 0], false).gt, "./0");
        assert_eq!(SampleStats::new(&[0, 1, 1], true).gt, "0|1|1");
    }
```

- [ ] **Step 2: Run it — expect PASS on current code (baseline)**

```bash
cargo test --all-features generate::tests::sample_stats_gt 2>&1 | tail -8
```

Expected: PASS (pins current behavior before the refactor).

- [ ] **Step 3: Rewrite the GT construction to a single reused `String`**

`src/bulk/generate.rs`, in `SampleStats::new`, replace the
`.map().collect::<Vec<_>>().join()`:

```rust
        let sep = if phased { '|' } else { '/' };
        let mut gt = String::with_capacity(alleles.len() * 2);
        for (i, a) in alleles.iter().enumerate() {
            if i > 0 {
                gt.push(sep);
            }
            if *a < 0 {
                gt.push('.');
            } else {
                use std::fmt::Write as _;
                let _ = write!(gt, "{a}");
            }
        }
```

(`write!` into the existing buffer avoids the per-allele `a.to_string()` and
stays correct if allele indices ever exceed 9.)

- [ ] **Step 4: Run the test again — still PASS**

```bash
cargo test --all-features generate:: 2>&1 | tail -10
```

Expected: PASS, including `payload_presets_produce_the_right_format_keys`.

- [ ] **Step 5: Lint and commit**

```bash
cargo clippy --all-features -- -D warnings && cargo fmt --check
git add src/bulk/generate.rs
git commit -m "perf(bulk): build GT into a reused String, not Vec<String>+join (#10)

SampleStats::new allocated ~3 strings per sample per record; at 100 MB
(265k records x 3202 samples) that is ~2.5B allocations. Write digits and
the separator straight into one preallocated buffer. Output unchanged."
```

---

## Task C2: Rejection-sample alt-allele placement when `ac << n_alleles`

**Files:**
- Modify: `src/bulk/generate.rs` (`gen_record` `:112-125`)
- Test: `src/bulk/generate.rs` tests
- Re-baseline: any golden output in `tests/bulk.rs` (RNG draw pattern changes)

**Interfaces:**
- Consumes: nothing new.
- Produces: same *invariant* — exactly `ac_eff = min(ac, n_nonmissing)` alt
  alleles placed among non-missing slots — but a different RNG consumption
  pattern, so byte-level output for a given seed changes. `gen_record`'s
  signature is unchanged.

- [ ] **Step 1: Write the invariant test (survives the algorithm change)**

In `src/bulk/generate.rs` tests:

```rust
    #[test]
    fn exactly_ac_eff_alt_alleles_are_placed() {
        let (p, s) = fixture();
        for i in 0..200u64 {
            let mut rng = block_rng(7, i);
            let r = gen_record(&mut rng, &s, "chr1", 100 + i, 1000, 2, &p.fitted);
            let n_alt = r.gts.iter().filter(|&&g| g == 1).count();
            let n_missing = r.gts.iter().filter(|&&g| g == -1).count();
            let n_nonmissing = r.gts.len() - n_missing;
            assert!(n_alt <= n_nonmissing);
        }
    }
```

- [ ] **Step 2: Capture the current golden output baseline**

```bash
cargo test --all-features 2>&1 | tee /tmp/pre_c2.txt | tail -5
```

- [ ] **Step 3: Replace the partial Fisher-Yates with a two-branch placement**

`src/bulk/generate.rs`, `gen_record`, replace lines 119-125 (the `let mut idx
... for i in 0..ac_eff {...}` block):

```rust
    // Number of non-missing slots, and the exact alt count to place.
    let n_nonmissing = gts.iter().filter(|&&g| g != -1).count();
    let ac_eff = (ac as usize).min(n_nonmissing);

    // When ac_eff is small relative to n_nonmissing (the common case: 36% of
    // records are singletons), rejection-sample distinct slot ranks instead
    // of materialising and shuffling all n_nonmissing indices.
    if ac_eff * 2 <= n_nonmissing {
        let mut chosen: std::collections::HashSet<usize> =
            std::collections::HashSet::with_capacity(ac_eff);
        while chosen.len() < ac_eff {
            chosen.insert(rng.random_range(0..n_nonmissing));
        }
        // rank = position among non-missing slots, left to right
        let mut rank = 0usize;
        for g in gts.iter_mut() {
            if *g != -1 {
                if chosen.contains(&rank) {
                    *g = 1;
                }
                rank += 1;
            }
        }
    } else {
        // dense: partial Fisher-Yates over the non-missing index list
        let mut idx: Vec<usize> = (0..n_alleles as usize).filter(|&i| gts[i] != -1).collect();
        for i in 0..ac_eff {
            let j = rng.random_range(i..idx.len());
            idx.swap(i, j);
            gts[idx[i]] = 1;
        }
    }
```

Update the module doc / inline comment above to describe the two-branch
placement (sparse rejection vs dense shuffle), preserving the existing
explanation of *why* exact-AC placement is used at all.

> **Determinism note:** `HashSet` insertion order never influences output — we
> draw ranks from `rng` (deterministic) and apply membership while walking
> `gts` in index order. No hash-map *iteration* touches output.

- [ ] **Step 4: Run the invariant and determinism tests**

```bash
cargo test --all-features generate:: 2>&1 | tail -15
```

Expected: `exactly_ac_eff_alt_alleles_are_placed` passes;
`generation_is_deterministic` still passes (same seed twice → same output).

- [ ] **Step 5: Re-baseline shifted golden output**

```bash
cargo test --all-features 2>&1 | tee /tmp/post_c2.txt | tail -15
diff /tmp/pre_c2.txt /tmp/post_c2.txt || true
```

Any `tests/bulk.rs` golden hash / snapshot that changed is expected (draw
pattern changed). Re-baseline and annotate:
`// re-baselined: alt placement switched to rejection sampling (#10)`.
Structural/statistical assertions (SFS fidelity, determinism) must stay green
unchanged — if one breaks, stop and debug.

- [ ] **Step 6: Lint and commit**

```bash
cargo clippy --all-features -- -D warnings && cargo fmt --check
git add src/bulk/generate.rs tests/bulk.rs
git commit -m "perf(bulk): rejection-sample sparse alt placement (#10)

gen_record allocated a Vec<usize> of n_alleles (up to 6404) every record
to place a median of 1 alt allele. Rejection-sample distinct ranks when
ac<<n_alleles; keep Fisher-Yates for the dense case. Exact-AC invariant
preserved; golden output re-baselined for the new draw pattern."
```

---

## Task C3: Split `Size::Records` by `n_variants`, not `density_per_kb`

**Files:**
- Modify: `src/bulk/mod.rs` (`distribute_by_density` `:733-738`, its doc, and
  the `Size` doc mentioning "density"); call sites `:329`, `:535`, `:548`
- Test: `src/bulk/mod.rs` tests (private-fn access) or `tests/bulk.rs`

**Interfaces:**
- Consumes: `ContigStat.n_variants` (already fitted, previously unread).
- Produces: `distribute_by_density` renamed to `distribute_by_n_variants`,
  weighting on `n_variants`; per-contig split reproduces the source's variant
  distribution. Same signature: `fn distribute_by_n_variants(fitted: &Fitted,
  contig_ids: &[String], total: u64) -> Vec<u64>`.

- [ ] **Step 1: Write a discriminating test (n_variants vs density disagree)**

Add to `src/bulk/mod.rs`'s `#[cfg(test)] mod tests` (so it can see the private
fn):

```rust
    // A profile where n_variants order != density order, so the two split
    // strategies give different answers.
    const DISCRIMINATING_PROFILE: &str = r#"{
      "name": "disc", "provenance": {"source":"x","n_samples_source":10,
        "n_variants_source":11000,"fitted_on":"2026-01-01","fit_tool_version":"t",
        "supplied":["ploidy"]},
      "fitted": { "contigs": [
          { "id": "big", "n_variants": 10000, "density_per_kb": 10.0 },
          { "id": "small", "n_variants": 1000, "density_per_kb": 90.0 }
        ],
        "gap_dist": {"edges":[1.0,2.0],"weights":[1.0]},
        "sfs": {"edges":[1.0,2.0],"weights":[1.0]},
        "variant_classes": {"snp":1.0,"insertion":0.0,"deletion":0.0,"mnp":0.0,"complex":0.0,"symbolic":0.0},
        "indel_length": {"edges":[1.0,2.0],"weights":[1.0]},
        "titv": 2.0, "multiallelic_rate": 0.0, "missing_rate": 0.0, "phased_rate": 1.0
      },
      "dialed": { "payload": "gt-only", "ploidy": 2 }
    }"#;

    #[test]
    fn records_split_follows_n_variants_not_density() {
        let p = Profile::from_json(DISCRIMINATING_PROFILE).unwrap();
        let ids = vec!["big".to_string(), "small".to_string()];
        let counts = distribute_by_n_variants(&p.fitted, &ids, 11_000);
        assert_eq!(counts, vec![10_000, 1_000]);
    }
```

- [ ] **Step 2: Run it — expect failure (fn not yet renamed/reweighted)**

```bash
cargo test --all-features records_split_follows_n_variants 2>&1 | tail -12
```

Expected: FAIL to compile (`distribute_by_n_variants` undefined) — that is the
red state; it becomes a value mismatch only if you rename without reweighting.

- [ ] **Step 3: Rename and reweight**

`src/bulk/mod.rs:733`, rename `distribute_by_density` → `distribute_by_n_variants`
and weight on `n_variants`:

```rust
fn distribute_by_n_variants(fitted: &Fitted, contig_ids: &[String], total: u64) -> Vec<u64> {
    let weights: Vec<f64> = contig_ids
        .iter()
        .enumerate()
        .map(|(i, id)| resolve_contig_stat(fitted, i, id).n_variants as f64)
        .collect();
    // ... rest unchanged (largest-remainder, even-split fallback) ...
```

Update the call sites: `:329` (`Size::Records`), and `:535`/`:548` inside
`resolve_target_counts` (C4 rewrites that function; if C4 lands first, its
calibration split already calls `distribute_by_n_variants`). Update
`distribute_by_density`'s doc and the `Size::Records`/`Size::Target` doc
comments: "proportional to fitted density" → "proportional to fitted per-contig
variant count (`n_variants`)".

- [ ] **Step 4: Run tests**

```bash
cargo test --all-features 2>&1 | tail -15
```

Expected: the new test passes; existing split tests pass (or are updated where
they asserted density-specific numbers — note the change inline).

- [ ] **Step 5: Lint and commit**

```bash
cargo clippy --all-features -- -D warnings && cargo fmt --check
git add src/bulk/mod.rs tests/bulk.rs
git commit -m "fix(bulk): split records by n_variants, not density_per_kb (#10)

Output density is 1/mean(gap) globally (gap_dist isn't per-contig), so
fitted per-contig density was never reproduced. Weighting the split by the
fitted-but-unread n_variants reproduces the source's per-contig variant
distribution and stops the MT outlier (350/kb, 12x) skewing the split."
```

---

## Task C4: Replace the byte-target search with two-point calibration + promote

**Files:**
- Modify: `src/bulk/mod.rs` (`resolve_target_counts` `:523-557`; add
  `write_to_temp`; `write` `:327-381` Target branch; add `promote_temp` +
  `move_file`)
- Test: `tests/bulk.rs` (existing overshoot bound `:159-186`, plus a new
  determinism/byte-exactness test)

**Interfaces:**
- Consumes: `measure_compressed_bytes` (`:581`, unchanged, bytes-only, for the
  two cheap calibration points), `distribute_by_n_variants` (C3), `Summary`,
  `BulkWriter`, `tempfile::NamedTempFile`.
- Produces:
  - `fn resolve_target_counts(...) -> Result<(Vec<u64>, tempfile::NamedTempFile, u64, Summary), BulkError>`
    — returns the winning counts, the already-written temp file (byte-len
    `>= target`), and its summary.
  - `fn write_to_temp(&self, pool, samplers, fitted, counts) -> Result<(NamedTempFile, u64, Summary), BulkError>`.
  - `fn promote_temp(tmp: NamedTempFile, dest: &Path, format: Format) -> Result<(), BulkError>`
    and free `fn move_file(src: &Path, dst: &Path) -> Result<(), BulkError>`.

- [ ] **Step 1: Write the acceptance test — target reached, overshoot bounded, deterministic**

Keep the existing `tests/bulk.rs:159-186` overshoot test (`got >= target`,
`got <= target * 1.25`). Add a determinism test proving the calibrate+promote
path is reproducible (weaker than resolved-count exactness but sufficient, and
needs no new `Summary` accessor):

```rust
#[test]
fn target_size_is_byte_identical_across_runs() {
    let dir = tempfile::tempdir().unwrap();
    let profile = Profile::from_json(NONUNIFORM_DENSITY_PROFILE).unwrap();
    let run = |name: &str| {
        let out = dir.path().join(name);
        BulkSpec::new(profile.clone(), 4, vec!["1".into(), "2".into(), "3".into()])
            .size(Size::Target(512 * 1024))
            .seed(99)
            .write(&out)
            .unwrap();
        std::fs::read(&out).unwrap()
    };
    let a = run("a.vcf.gz");
    let b = run("b.vcf.gz");
    assert!(a.len() as u64 >= 512 * 1024, "must reach target");
    assert_eq!(a, b, "calibrate+promote must be byte-identical run to run");
}
```

(Match `BulkSpec`'s real builder API — adjust `::new(...).size(...).seed(...)`
to whatever the crate exposes; the existing overshoot test shows the exact form.)

- [ ] **Step 2: Run it — expect the current code to pass overshoot but be slow**

```bash
cargo test --all-features target_size 2>&1 | tail -12
```

The determinism test may already pass under the old code (it is deterministic
too); its purpose is to guard the rewrite. Proceed regardless.

- [ ] **Step 3: Add `write_to_temp` (write pass building the Summary, keeps the file)**

In `src/bulk/mod.rs`, alongside `measure_compressed_bytes`:

```rust
    /// Like `measure_compressed_bytes`, but builds the `Summary` during the
    /// write pass and returns the live temp file instead of deleting it, so
    /// the caller can promote it to the real destination. Byte-exact:
    /// identical header, records, and (absent) flush cadence as `write()`.
    fn write_to_temp(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        per_contig_count: &[u64],
    ) -> Result<(tempfile::NamedTempFile, u64, Summary), BulkError> {
        let spans: Vec<u64> = self
            .contig_ids
            .iter()
            .zip(per_contig_count)
            .enumerate()
            .map(|(i, (id, &n))| {
                let recs = self.generate_contig(pool, samplers, fitted, id, i as u64, n);
                contig_span(&recs)
            })
            .collect();

        let header = self.build_header(&spans);
        let tmp = tempfile::NamedTempFile::new()?;
        let tmp_path = tmp.path().to_path_buf();

        let mut w = BulkWriter::create(
            &tmp_path, self.format, &header, self.compression_level, self.workers,
        )?;
        let mut summary = Summary::new(self.n_samples);
        for (i, (id, &n)) in self.contig_ids.iter().zip(per_contig_count).enumerate() {
            let recs = self.generate_contig(pool, samplers, fitted, id, i as u64, n);
            for r in &recs {
                let buf = to_record_buf(&r.g, self.payload.clone(), r.phased);
                w.write(&header, &buf)?;
                summary.observe(id, r.g.pos, r.g.class, &r.g.gts);
            }
        }
        w.finish_and_index(&tmp_path)?;
        let bytes = std::fs::metadata(&tmp_path)?.len();
        Ok((tmp, bytes, summary))
    }
```

- [ ] **Step 4: Add `promote_temp` and `move_file`**

```rust
    /// Move a written temp file (and, for BCF, its `.csi` companion) onto the
    /// real destination. Rename when possible; fall back to copy across
    /// filesystems (TMPDIR may differ from the output dir).
    fn promote_temp(
        tmp: tempfile::NamedTempFile,
        dest: &Path,
        format: Format,
    ) -> Result<(), BulkError> {
        let tmp_path = tmp.path().to_path_buf();
        if matches!(format, Format::Bcf) {
            let mut src_csi = tmp_path.as_os_str().to_os_string();
            src_csi.push(".csi");
            let mut dst_csi = dest.as_os_str().to_os_string();
            dst_csi.push(".csi");
            move_file(Path::new(&src_csi), Path::new(&dst_csi))?;
        }
        match tmp.persist(dest) {
            Ok(_) => Ok(()),
            Err(e) => {
                std::fs::copy(e.file.path(), dest)?;
                Ok(())
            }
        }
    }
```

Free function (place near the other module-level helpers):

```rust
fn move_file(src: &Path, dst: &Path) -> Result<(), BulkError> {
    if std::fs::rename(src, dst).is_ok() {
        return Ok(());
    }
    std::fs::copy(src, dst)?;
    let _ = std::fs::remove_file(src);
    Ok(())
}
```

- [ ] **Step 5: Rewrite `resolve_target_counts` as two-point calibration**

Replace the whole 25-round loop:

```rust
    fn resolve_target_counts(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        target_bytes: u64,
    ) -> Result<(Vec<u64>, tempfile::NamedTempFile, u64, Summary), BulkError> {
        let n_contigs = self.contig_ids.len() as u64;

        // Two calibration points; c2 = 2*c1 so the slope is well-conditioned.
        let split1 = distribute_by_n_variants(fitted, &self.contig_ids, 1_000 * n_contigs);
        let split2 = distribute_by_n_variants(fitted, &self.contig_ids, 2_000 * n_contigs);
        let bytes1 = self.measure_compressed_bytes(pool, samplers, fitted, &split1)?;
        let bytes2 = self.measure_compressed_bytes(pool, samplers, fitted, &split2)?;

        let r1 = split1.iter().sum::<u64>() as f64;
        let r2 = split2.iter().sum::<u64>() as f64;
        // bytes ~= b0 + k*records ; k bytes/record, b0 the fixed header cost.
        let k = ((bytes2 as f64 - bytes1 as f64) / (r2 - r1)).max(1e-9);
        let b0 = bytes1 as f64 - k * r1;

        // Direct count; never below the larger calibration (a known-good
        // measurement) and never below 1 record/contig.
        let want = (((target_bytes as f64 - b0) / k).ceil() as i64).max(r2 as i64) as u64;
        let mut counts = distribute_by_n_variants(fitted, &self.contig_ids, want);

        // Slope-based correction; converges in 1-2 rounds.
        const MAX_CORRECTIONS: usize = 4;
        for _ in 0..MAX_CORRECTIONS {
            let (tmp, bytes, summary) = self.write_to_temp(pool, samplers, fitted, &counts)?;
            if bytes >= target_bytes {
                return Ok((counts, tmp, bytes, summary));
            }
            let shortfall = (target_bytes - bytes) as f64;
            let extra = ((shortfall / k) * 1.02).ceil() as u64 + 1;
            let extra_split = distribute_by_n_variants(fitted, &self.contig_ids, extra);
            for (c, e) in counts.iter_mut().zip(&extra_split) {
                *c += e;
            }
            drop(tmp); // discard the under-target temp before regenerating
        }

        Err(BulkError::Invalid(format!(
            "could not reach target size {target_bytes} bytes within {MAX_CORRECTIONS} corrective rounds"
        )))
    }
```

Update the function's doc comment (`:485-522`) to describe calibration +
promote, dropping the "up to 25 rounds" narrative.

- [ ] **Step 6: Rewrite the `Size::Target` branch of `write` to promote**

`src/bulk/mod.rs`, in `write()`, at the `match self.size` (`:327`):

```rust
        let counts: Vec<u64> = match self.size {
            Size::RecordsPerContig(n) => vec![n; self.contig_ids.len()],
            Size::Records(total) => distribute_by_n_variants(fitted, &self.contig_ids, total),
            Size::Target(target_bytes) => {
                let (_counts, tmp, _bytes, summary) =
                    self.resolve_target_counts(&pool, &samplers, fitted, target_bytes)?;
                Self::promote_temp(tmp, path, self.format)?;
                let json = summary.to_json()?;
                let mut summary_path = path.as_os_str().to_os_string();
                summary_path.push(".summary.json");
                std::fs::write(&summary_path, json)?;
                return Ok(summary);
            }
        };
        // ... existing span pass + write pass for Records/RecordsPerContig ...
```

Leave the non-Target path below unchanged (span pass, `build_header`,
`BulkWriter::create`, write pass with `summary.observe`, `finish_and_index`,
summary json).

- [ ] **Step 7: Build, run the target tests and full suite**

```bash
cargo test --all-features target 2>&1 | tail -20
cargo test --all-features 2>&1 | tail -10
```

Expected: overshoot bound and byte-identical determinism tests pass; nothing
else regresses.

- [ ] **Step 8: Benchmark against the 226 s baseline**

```bash
cargo build --release --features cli
# #8 baseline: --target-size 8MB was 226 s, 50,590 records.
# Flags per src/bin/vcfixture.rs: --output (short -o), --samples,
# --target-size, --contigs (comma-delimited), --format (default bcf).
/usr/bin/time -v ./target/release/vcfixture bulk \
  --profile germline-1kgp --samples 100 --target-size 8MB \
  --contigs chr1,chr2,chr3 --seed 1 --format vcf-gz --output /tmp/t8.vcf.gz 2>&1 \
  | grep -E "Elapsed|Maximum resident"
ls -l /tmp/t8.vcf.gz
```

Expected: wall time low-tens-of-seconds (down from 226 s); output `>= 8 MB`.
Record the number for the commit.

- [ ] **Step 9: Lint and commit**

```bash
cargo clippy --all-features -- -D warnings && cargo fmt --check
git add src/bulk/mod.rs tests/bulk.rs
git commit -m "perf(bulk): calibrate byte target in 2 points, promote the file (#8)

The old search did up to 25 rounds, each generating every contig TWICE
(~50 generations), and its bytes/record divided by header-inclusive bytes
so every round undershot. Fit b0+k*records from two small measurements,
compute the count directly, correct at most once, and promote the winning
byte-exact temp file instead of regenerating. 8MB: 226 s -> <NN> s."
```

---

# Workstream D — fit-script memory (#7, #10 histogram bug)

## Task D1: `read_pvar` — skip `##` lines with `skip_lines`, drop `comment_prefix`

**Files:**
- Modify: `scripts/fit/fit_profile.py` (`read_pvar` `:317-323`; add
  `_count_leading_meta_lines`)
- Test: `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: a `.pvar`/`.pvar.zst` path.
- Produces: `read_pvar` returns the same 5-column lazy frame, but its `scan_csv`
  uses `skip_lines=n` (counted `##` lines) instead of `comment_prefix="##"` —
  removing the whole-file materialization that `comment_prefix` forces (measured
  6.2 GB just to count rows, 3× slower). New helper
  `_count_leading_meta_lines(path) -> int`.

- [ ] **Step 1: Write a test on a small synthetic pvar**

In `scripts/fit/test_fit_profile.py`:

```python
def test_read_pvar_skips_meta_and_reads_all_rows(tmp_path):
    p = tmp_path / "t.pvar"
    p.write_text(
        "##fileformat=PVARv1.0\n"
        "##contig=<ID=1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\t.\tA\tG\n"
        "1\t200\t.\tC\tT\n"
        "2\t50\t.\tG\tA\n"
    )
    lf = read_pvar(p)
    df = lf.collect()
    assert df.height == 3
    assert df["CHROM"].to_list() == ["1", "1", "2"]
    assert df["POS"].to_list() == [100, 200, 50]
```

- [ ] **Step 2: Run it — expect PASS on current code (baseline behavior)**

```bash
pixi run -e fit test-fit -k read_pvar_skips_meta 2>&1 | tail -10
```

Expected: PASS (current `comment_prefix` code reads the same rows; the test
pins behavior across the change).

- [ ] **Step 3: Add `_count_leading_meta_lines`, handling plain/gz/zst**

`scripts/fit/fit_profile.py`:

```python
def _count_leading_meta_lines(path: str | Path) -> int:
    """Count leading ``##`` metadata lines so scan_csv can `skip_lines` past
    them instead of `comment_prefix='##'` (which forces polars to materialise
    the whole file -- measured 6.2 GB / 3x slower on a 5.9 GB somatic pvar).

    Handles the pvar's optional .zst compression: the meta block is a small
    text prefix, so only the first chunk needs decompressing.
    """
    p = str(path)
    if p.endswith(".zst"):
        import zstandard  # in the `fit` env
        n = 0
        with open(p, "rb") as raw:
            dctx = zstandard.ZstdDecompressor()
            with dctx.stream_reader(raw) as r:
                buf = r.read(1 << 16).split(b"\n")
                for line in buf:
                    if line.startswith(b"##"):
                        n += 1
                    else:
                        break
        return n
    n = 0
    with open(p, "rb") as fh:
        for line in fh:
            if line.startswith(b"##"):
                n += 1
            else:
                break
    return n
```

> Verify `zstandard` is importable: `pixi run -e fit python -c "import zstandard"`.
> If absent, add it to `[feature.fit.dependencies]` in `pixi.toml`, or read the
> zst prefix by shelling `zstd -dc <file> | head` — but prefer the module.

- [ ] **Step 4: Rewrite `read_pvar` to use `skip_lines`**

```python
    n_meta = _count_leading_meta_lines(path)
    lf = pl.scan_csv(
        str(path),
        separator="\t",
        skip_lines=n_meta,
        schema_overrides={"#CHROM": pl.Utf8, "ID": pl.Utf8, "REF": pl.Utf8, "ALT": pl.Utf8},
    ).rename({"#CHROM": "CHROM"})
    # ID is carried only to keep read_pvar/read_sites_vcf schemas identical
    # for compute_pvar_stats; nothing reads it, and in this file it is "."
    # (REF/ALT carry the wide strings), so dropping it saves ~0 memory --
    # keep it for schema parity with read_sites_vcf's fabricated ID.
    return lf.select(["CHROM", "POS", "ID", "REF", "ALT"])
```

- [ ] **Step 5: Run the test and confirm the real header count**

```bash
pixi run -e fit test-fit -k read_pvar 2>&1 | tail -10
pixi run -e fit python -c "
import sys; sys.path.insert(0,'scripts/fit'); import fit_profile as fp
print(fp._count_leading_meta_lines('/carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt.pvar'))
"
```

Expected: test PASS; header count = `43`.

- [ ] **Step 6: Measure the scan floor on real data (commit evidence)**

```bash
/usr/bin/time -v pixi run -e fit python -c "
import sys; sys.path.insert(0,'scripts/fit'); import fit_profile as fp
import polars as pl
print(fp.read_pvar('/carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt.pvar').select(pl.len()).collect(engine='streaming').to_dicts())
" 2>&1 | grep -E "len|Maximum resident|Elapsed"
```

Expected: 348,259,675 rows, peak RSS ~6.2 GB, ~3× faster than the
`comment_prefix` scan.

- [ ] **Step 7: Commit**

```bash
git add scripts/fit/fit_profile.py scripts/fit/test_fit_profile.py
git commit -m "perf(fit): scan pvar with skip_lines, not comment_prefix (#7)

comment_prefix='##' forces polars to materialise the whole 5.9 GB CSV
(6.2 GB just to count rows, 3x slower); bcftools output on the sites-vcf
path has no ## lines, which is the real reason that path was leaner.
Count the ## header once and skip_lines past it."
```

---

## Task D2: `_gap_bins_lazy` — drop the sort, replace the window with shift+mask

**Files:**
- Modify: `scripts/fit/fit_profile.py` (`_gap_bins_lazy` `:577-594`; add
  `assert_pvar_sorted`, called in `compute_pvar_stats`)
- Test: `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: the per-record `lf` (coordinate-sorted, as plink2 emits).
- Produces: the same `histogram_lazy` bin-count frame, computed without
  `.sort()` or `.over()`. New `assert_pvar_sorted(lf) -> None`. Measured
  genome-wide: current sort+window ~44-57 B/row (~20 GB at 348M rows);
  shift+mask ~0.6 B/row (220 MB), bit-identical (chr21 gaps n=4,107,307,
  total=41,668,608; genome-wide 318,372,446 gaps).

- [ ] **Step 1: Write a test asserting shift+mask == sort+window on sorted input**

In `scripts/fit/test_fit_profile.py`:

```python
import polars as pl

def test_gap_bins_matches_sorted_window_reference():
    lf = pl.LazyFrame({
        "CHROM": ["1","1","1","2","2"],
        "POS":   [100, 250, 251, 5, 30],
        "ID": ["."]*5, "REF": ["A"]*5, "ALT": ["T"]*5,
    })
    got = _gap_bins_lazy(lf).collect()
    ref_gaps = (lf.sort(["CHROM","POS"])
        .select(pl.col("POS").diff().over("CHROM").alias("gap"))
        .filter(pl.col("gap").is_not_null() & (pl.col("gap") > 0)))
    ref = histogram_lazy(ref_gaps, pl.col("gap"), _gap_edges()).collect()
    assert got.equals(ref)
```

- [ ] **Step 2: Run it — expect PASS on current code (baseline)**

```bash
pixi run -e fit test-fit -k gap_bins_matches 2>&1 | tail -10
```

Expected: PASS (compares the function to its own current formula).

- [ ] **Step 3: Rewrite `_gap_bins_lazy`**

```python
def _gap_bins_lazy(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Inter-variant gap (bp) histogram bin counts, within each contig.

    `lf` is per-RECORD and is assumed coordinate-sorted within each contig
    (plink2 emits pvar sorted by CHROM, POS). Gaps are a straight
    `POS.diff()` masked to same-contig adjacent rows
    (`CHROM == CHROM.shift(1)`), NOT `sort().diff().over("CHROM")`: the sort
    is a full-frame pipeline breaker and `.over()` is a non-streaming window,
    which together cost ~20 GB at genome scale (348M rows). shift+mask
    streams in ~220 MB and is bit-identical on sorted input (guarded by
    `assert_pvar_sorted`).
    """
    gaps = (
        lf.select(
            pl.col("POS").diff().alias("gap"),
            (pl.col("CHROM") == pl.col("CHROM").shift(1)).alias("same_contig"),
        )
        .filter(
            pl.col("same_contig")
            & pl.col("gap").is_not_null()
            & (pl.col("gap") > 0)
        )
        .select("gap")
    )
    return histogram_lazy(gaps, pl.col("gap"), _gap_edges())
```

- [ ] **Step 4: Add `assert_pvar_sorted` and call it in `compute_pvar_stats`**

```python
def assert_pvar_sorted(lf: pl.LazyFrame) -> None:
    """Fail if any within-contig POS is out of order (the precondition
    `_gap_bins_lazy` relies on after dropping its sort). Streams: one boolean
    reduction, bounded memory."""
    bad = (
        lf.select(
            (pl.col("CHROM") == pl.col("CHROM").shift(1)).alias("same"),
            (pl.col("POS") < pl.col("POS").shift(1)).alias("descending"),
        )
        .select((pl.col("same") & pl.col("descending")).sum().alias("n_desc"))
        .collect(engine="streaming")
        .item()
    )
    if bad:
        raise ValueError(
            f"pvar is not sorted within contigs ({bad} descending steps); "
            "gap fitting requires coordinate-sorted input"
        )
```

Call `assert_pvar_sorted(lf)` once near the top of `compute_pvar_stats`, before
building the six plans.

- [ ] **Step 5: Add a test that the precondition fires on unsorted input**

```python
def test_assert_pvar_sorted_rejects_unsorted():
    import pytest
    lf = pl.LazyFrame({"CHROM": ["1","1"], "POS": [200, 100],
                       "ID": [".","."], "REF": ["A","A"], "ALT": ["T","T"]})
    with pytest.raises(ValueError, match="not sorted"):
        assert_pvar_sorted(lf)
```

- [ ] **Step 6: Run the tests**

```bash
pixi run -e fit test-fit -k "gap_bins or pvar_sorted" 2>&1 | tail -12
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add scripts/fit/fit_profile.py scripts/fit/test_fit_profile.py
git commit -m "perf(fit): compute gaps by shift+mask, not sort+window (#7)

sort(['CHROM','POS']) is a full-frame pipeline breaker feeding a
non-streaming .over() window; together ~20 GB at 348M rows. The pvar is
already coordinate-sorted, so POS.diff() masked to same-contig adjacency
is bit-identical and streams in ~220 MB. Guarded by assert_pvar_sorted."
```

---

## Task D3: `_titv_lazy` — replace `is_in(TRANSITION_PAIRS)` with direct comparisons

**Files:**
- Modify: `scripts/fit/fit_profile.py` (`_titv_lazy` `:630-636`; maybe remove
  `TRANSITION_PAIRS`)
- Test: `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: the per-allele `alleles` frame (`class`, `REF`, `ALT`).
- Produces: the same `(n_snps, n_ts)` one-row frame. Measured genome-wide:
  `concat_str().is_in(...)` cost 16.7 GB; direct comparisons cost 6.4 GB (the
  scan floor), identical `n_ts = 163,208,320`.

- [ ] **Step 1: Write an equivalence test**

```python
def test_titv_direct_matches_is_in_reference():
    alleles = pl.LazyFrame({
        "class": ["snp","snp","snp","snp","insertion"],
        "REF":   ["A","G","C","A","A"],
        "ALT":   ["G","A","T","C","AT"],  # A>G, G>A, C>T ts; A>C tv
    })
    got = _titv_lazy(alleles).collect().row(0)  # (n_snps, n_ts)
    assert got == (4, 3)
```

- [ ] **Step 2: Run it — expect PASS on current code (baseline)**

```bash
pixi run -e fit test-fit -k titv_direct 2>&1 | tail -8
```

Expected: PASS.

- [ ] **Step 3: Rewrite `_titv_lazy`**

```python
def _titv_lazy(alleles: pl.LazyFrame) -> pl.LazyFrame:
    """SNP count and transition count (both scalars) over a per-ALLELE frame.

    Transitions are the four purine<->purine / pyrimidine<->pyrimidine pairs,
    expressed as direct (REF,ALT) comparisons rather than
    `concat_str([REF,ALT]).is_in(TRANSITION_PAIRS)`: measured, the `is_in`
    against a literal list costs ~10 GB extra at 348M rows (16.7 GB vs the
    6.4 GB scan floor) for a bit-identical result.
    """
    snps = alleles.filter(pl.col("class") == "snp")
    r, a = pl.col("REF"), pl.col("ALT")
    is_ts = (
        ((r == "A") & (a == "G"))
        | ((r == "G") & (a == "A"))
        | ((r == "C") & (a == "T"))
        | ((r == "T") & (a == "C"))
    )
    return snps.select(n_snps=pl.len(), n_ts=is_ts.sum())
```

If `TRANSITION_PAIRS` is now unused (`rg TRANSITION_PAIRS scripts/fit/`), remove it.

- [ ] **Step 4: Run tests**

```bash
pixi run -e fit test-fit -k titv 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/fit/fit_profile.py scripts/fit/test_fit_profile.py
git commit -m "perf(fit): compute Ti/Tv by direct base compares, not is_in (#7)

concat_str([REF,ALT]).is_in(TRANSITION_PAIRS) cost ~10 GB extra at 348M
rows (16.7 GB vs the 6.4 GB scan floor). Four direct (REF,ALT) comparisons
give a bit-identical n_ts at the floor."
```

---

## Task D4: `compute_pvar_stats` — collect plans sequentially, not `pl.collect_all`

**Files:**
- Modify: `scripts/fit/fit_profile.py` (`compute_pvar_stats` `:684-687` +
  docstring `:671-673`)
- Test: `scripts/fit/test_fit_profile.py` (existing `compute_pvar_stats` test
  as regression guard)

**Interfaces:**
- Consumes: the six lazy plans already built in `compute_pvar_stats`.
- Produces: the same six result frames, collected one at a time. Measured:
  `pl.collect_all` CSE caches the shared 348M-row scan/explode output, +20 GB
  (26.4 GB) for no gain; sequential collect is 6.4 GB and faster.

- [ ] **Step 1: Note the current expected values of the stats test**

```bash
pixi run -e fit test-fit -k compute_pvar_stats 2>&1 | tail -12
```

Record the expected dict values; they must be identical after Step 2. (If no
`compute_pvar_stats` test exists, add one on a small synthetic frame asserting
`titv`, `multiallelic_rate`, contig count, and a gap-bin total.)

- [ ] **Step 2: Replace `pl.collect_all` with sequential `.collect`**

`scripts/fit/fit_profile.py:684-687`:

```python
    # Collect each plan on its own rather than pl.collect_all: collect_all's
    # common-subplan elimination *caches* the shared 348M-row scan/explode
    # output, costing +20 GB (26.4 GB vs 6.4 GB) and buying nothing -- a
    # re-scan from warm page cache is both leaner and faster than holding a
    # materialised copy. Measured on the 348M-row somatic pvar.
    contig_df = contig_lf.collect(engine="streaming")
    gap_bins_df = gap_bins_lf.collect(engine="streaming")
    multi_df = multi_lf.collect(engine="streaming")
    class_df = class_lf.collect(engine="streaming")
    indel_bins_df = indel_bins_lf.collect(engine="streaming")
    titv_df = titv_lf.collect(engine="streaming")
```

Replace the stale docstring paragraph (`:671-673`) claiming `collect_all` runs
the shared scans once "via common-subplan elimination" — it is measurably false
(that caching is the +20 GB). State the sequential rationale.

- [ ] **Step 3: Run the regression test — values identical**

```bash
pixi run -e fit test-fit -k compute_pvar_stats 2>&1 | tail -12
```

Expected: identical values to Step 1.

- [ ] **Step 4: End-to-end genome-wide smoke test (the payoff)**

```bash
/usr/bin/time -v pixi run -e fit python -c "
import sys; sys.path.insert(0,'scripts/fit'); import fit_profile as fp
s = fp.compute_pvar_stats(fp.read_pvar('/carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt.pvar'))
print('contigs', len(s['contigs']), 'titv', round(s['titv'],4))
" 2>&1 | grep -E "contigs|Maximum resident|Elapsed"
```

Expected: completes under 32 GB — peak RSS ~6.4 GB, ~2.5 min (vs OOM-killed
before D2+D3+D4). Record the numbers.

- [ ] **Step 5: Commit**

```bash
git add scripts/fit/fit_profile.py scripts/fit/test_fit_profile.py
git commit -m "perf(fit): collect pvar-stat plans sequentially, not collect_all (#7)

collect_all's CSE caches the shared 348M-row scan/explode output: +20 GB
(26.4 GB vs 6.4 GB) for no gain, since re-scanning from warm page cache is
leaner and faster. With D2+D3 this takes the genome-wide somatic fit from
OOM to 6.4 GB / ~2.5 min. Docstring corrected."
```

---

## Task D5: Fix `_bucket_index_expr` single-bin off-by-one

**Files:**
- Modify: `scripts/fit/fit_profile.py` (`_bucket_index_expr` `:152-180`)
- Test: `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: a value expression and `edges`.
- Produces: for `n_bins == 1`, a value equal to `edges[-1]` now lands in bin 0
  (closed on both ends), matching `numpy.histogram`; previously dropped.

- [ ] **Step 1: Write the failing test (numpy parity at n_bins == 1)**

```python
import numpy as np

def test_bucket_index_single_bin_includes_last_edge():
    edges = [1.0, 10.0]  # one bin, closed [1, 10]
    df = pl.DataFrame({"v": [1.0, 5.0, 10.0, 10.0, 0.5, 11.0]})
    got = df.select(_bucket_index_expr(pl.col("v"), edges).alias("b"))["b"].to_list()
    counts, _ = np.histogram(df["v"].to_numpy(), bins=edges)
    assert counts.tolist() == [4]
    assert got == [0, 0, 0, 0, None, None]
```

- [ ] **Step 2: Run it — expect failure (10.0 dropped)**

```bash
pixi run -e fit test-fit -k bucket_index_single_bin 2>&1 | tail -10
```

Expected: FAIL — current code gives `[0, 0, None, None, None, None]` (bin 0 is
half-open `[1,10)` and the closed-last clause is skipped when `n_bins == 1`).

- [ ] **Step 3: Fix the single-bin case**

`scripts/fit/fit_profile.py`, `_bucket_index_expr`:

```python
    edges = _validate_edges(edges)
    n_bins = len(edges) - 1
    if n_bins == 1:
        # single bin is closed on both ends, matching numpy.histogram
        return (
            pl.when((value >= edges[0]) & (value <= edges[1]))
            .then(pl.lit(0, dtype=pl.Int64))
            .otherwise(None)
        )
    expr = pl.when((value >= edges[0]) & (value < edges[1])).then(pl.lit(0, dtype=pl.Int64))
    for i in range(1, n_bins - 1):
        expr = expr.when((value >= edges[i]) & (value < edges[i + 1])).then(
            pl.lit(i, dtype=pl.Int64)
        )
    expr = expr.when((value >= edges[-2]) & (value <= edges[-1])).then(
        pl.lit(n_bins - 1, dtype=pl.Int64)
    )
    return expr.otherwise(None)
```

(With the early return, the trailing closed-bin clause runs only for
`n_bins >= 2`, where `edges[-2]` is a distinct edge — so its `if n_bins > 1`
guard is gone.)

- [ ] **Step 4: Run the test and the existing histogram tests**

```bash
pixi run -e fit test-fit -k "bucket_index or histogram" 2>&1 | tail -12
```

Expected: PASS, including existing multi-bin parity tests (unchanged path).

- [ ] **Step 5: Commit**

```bash
git add scripts/fit/fit_profile.py scripts/fit/test_fit_profile.py
git commit -m "fix(fit): count edges[-1] in single-bin histograms (#10)

_bucket_index_expr skipped the closed-last-bin clause when n_bins==1,
dropping any value equal to edges[-1] (reachable via _sfs_edges(1)),
contradicting its numpy-parity docstring. Close the single bin on both ends."
```

---

## Task D6: Re-fit the genome-wide somatic profile and commit it

**Files:**
- Modify: `profiles/somatic-gdc.json` (regenerated genome-wide)
- Modify: `scripts/fit/README.md`, docs referencing the somatic scope
- Depends on: B1, B2 (schema), B4 (validator), D1–D5 (memory)

**Interfaces:**
- Consumes: the fixed `fit_profile.py` and the real somatic `.pgen` fileset.
- Produces: a genome-wide `somatic-gdc.json` — `n_variants_source` ~7.9M →
  ~348M, ~25 contigs — carrying the new schema (`ploidy` under `dialed`,
  `provenance.supplied`), validating via the B4 gate.

- [ ] **Step 1: Run the genome-wide fit under budget**

```bash
/usr/bin/time -v pixi run -e fit python scripts/fit/fit_profile.py \
  --name somatic-gdc \
  --pgen /carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt \
  --payload gt-vaf --ploidy 2 \
  --out /tmp/somatic-gdc.genomewide.json 2>&1 \
  | grep -E "Maximum resident|Elapsed" | tail -4
```

(No `--contigs` restriction — the whole point is the full genome. Match the real
CLI flags in `fit_profile.py`'s argparser.) Expected: peak RSS well under 32 GB
(~7 GB), minutes.

- [ ] **Step 2: Validate the fresh profile through the Rust gate**

```bash
cargo run -q --features bulk --bin validate-profile -- /tmp/somatic-gdc.genomewide.json
```

Expected: `OK (~25 contigs)`.

- [ ] **Step 3: Sanity-check schema and scale**

```bash
python3 -c "
import json; d=json.load(open('/tmp/somatic-gdc.genomewide.json'))
print('n_variants_source', d['provenance']['n_variants_source'])
print('contigs', len(d['fitted']['contigs']))
print('supplied', d['provenance']['supplied'])
print('dialed', d['dialed'])
assert 'ploidy' not in d['fitted'], 'ploidy must be under dialed'
assert d['provenance']['n_variants_source'] > 300_000_000, 'expected genome-wide'
print('OK')
"
```

Expected: ~348M variants, ~25 contigs, `supplied` contains `"ploidy"`,
`dialed.ploidy == 2`, no `fitted.ploidy`.

- [ ] **Step 4: Replace the committed profile and re-run the Rust suite**

```bash
cp /tmp/somatic-gdc.genomewide.json profiles/somatic-gdc.json
cargo test --all-features 2>&1 | tail -12
```

Expected: `builtin_somatic_loads_and_validates` and all others pass. The
existing test asserts `n_samples_source == 16007` (same fileset → unchanged). If
any test asserts the old 2-contig scale, update it to the genome-wide reality
with a comment.

- [ ] **Step 5: Update docs referencing the somatic profile's scope**

```bash
rg -n "chr21|21 22|chr21\+22|two-contig|two contig" docs scripts profiles
```

Change "chr21+22" descriptions of `somatic-gdc` to "genome-wide" in
`scripts/fit/README.md` and any design/book doc.

- [ ] **Step 6: Commit**

```bash
git add profiles/somatic-gdc.json tests/bulk.rs scripts/fit/README.md docs/
git commit -m "data(fit): re-fit somatic-gdc genome-wide (#7)

The shipped profile described only chr21+22 because a genome-wide --pgen
fit OOMed at 25.4 GB. With the D1-D5 memory fixes the full 348M-row fit
peaks at ~6.4 GB, so somatic-gdc now describes the whole genome
(n_variants_source ~7.9M -> ~348M, ~25 contigs) with the new schema
(ploidy under dialed, provenance.supplied). Closes #7."
```

---

## Final integration check (after all tasks)

- [ ] **Rust gates:**

```bash
cargo test --all-features && cargo clippy --all-features -- -D warnings && cargo fmt --check
```

- [ ] **Fit gates:**

```bash
pixi run -e fit test-fit
pixi run -e fit test-fidelity   # if the somatic re-fit is in its fixture set
```

- [ ] **Profiles validate:**

```bash
for f in profiles/*.json; do cargo run -q --features bulk --bin validate-profile -- "$f"; done
```

- [ ] **rust-analyzer:** no syntax errors across `src/bulk/`; `profile.rs`
  linked into the crate graph (issue #6 acceptance).

- [ ] **Issue closure map:** #6 (A1+A2), #7 (D1–D6), #8 (C4), #9 (B1+B2),
  #10 (B3, B4, B5, C1, C2, C3, D5) — every checkbox above maps to one.
