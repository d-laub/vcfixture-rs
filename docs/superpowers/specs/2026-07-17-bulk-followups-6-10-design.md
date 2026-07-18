# Bulk-generation follow-ups (issues #6–#10)

**Date:** 2026-07-17
**Status:** approved, ready for planning
**Scope:** five open GitHub issues filed during the bulk-generation branch's
final review, tackled as one spec across four independent workstreams.

## Context

PR #5 (`feat(bulk): bulk VCF/BCF generation from data-fitted profiles`) merged a
system that fits a statistical profile from a real cohort and generates bulk
VCF/BCF fixtures from it. Its final review deferred five findings as issues
#6–#10. This spec addresses all five.

The findings split into four workstreams that touch mostly disjoint files:

- **A — hygiene** (#6): the `gen` module name and rand 0.8's `gen()` collide
  with the edition-2024 reserved keyword, breaking rust-analyzer.
- **B — schema** (#9, part of #10): `ploidy` is a generation choice living in
  the "measured-only" `Fitted` struct; provenance of hand-supplied fields is
  not auditable from the JSON; a diploid-only constraint is enforced at runtime
  rather than at parse time; nothing gates a fresh fit against the validator.
- **C — generator perf** (#8, part of #10): `Size::Target` does dozens of full
  regenerations to hit a byte target; two per-record allocation hot spots; the
  `Size::Records` split uses a weight that is never reproduced in the output.
- **D — fit-script memory** (#7, part of #10): `fit_profile.py --pgen` peak RSS
  scales with row count and OOMs at 348M rows; a single-bin histogram bug.

**Sequencing:** A and B are independent and cheap; C and D depend on B (the
somatic re-fit in D must emit B's new schema, and both should build on a settled
profile contract). Execute **A ∥ B first, then C ∥ D**.

### Empirical grounding

The memory and performance claims below were **measured on the real source
data** present on this node, not estimated:

- Somatic source: `/carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt`
  (`.pvar` 5.9 GB, 348,259,675 rows), 32 GB SLURM allocation (the "32 GB
  cgroup" in #7 is `SLURM_MEM_PER_NODE=32768`).
- polars 1.42.1, `pixi run -e fit`.

Two of the issues' stated leads turned out to be **wrong**, and two **additional
offenders the issues never named** were found. The measurements are recorded
inline in each workstream so the implementer builds on facts, not the issues'
hypotheses.

---

## Workstream A — hygiene (#6)

### Problem

`gen` is a reserved keyword in edition 2024. `src/bulk/mod.rs` declares
`pub mod gen;` and rand 0.8's API is `rng.gen()` / `rng.gen_range(a..b)`. The
crate is edition 2021 so **cargo is green** (`cargo test --all-features`,
`clippy -D warnings`, `fmt --check` all pass); only rust-analyzer breaks,
resolving these files under edition-2024 rules. Effect: no completion, no
go-to-definition, and a wall of false syntax errors across `src/bulk/`, plus
`profile.rs` reported as "not included in any crates."

### Design

Two changes that together remove **every** `gen` token:

1. **Rename module `gen` → `generate`.** Update `pub mod gen;` in
   `src/bulk/mod.rs` and every `bulk::gen::` / `gen::` path reference. Also
   forward-compatible with an edition-2024 migration (where `pub mod gen;`
   becomes a hard error needing `r#gen`).
2. **Upgrade rand 0.8 → 0.9** (and `rand_chacha` 0.3 → 0.9). API renames:
   `gen()` → `random()`, `gen_range()` → `random_range()`. Audit for other 0.9
   breaking changes (`Rng`/`SeedableRng` trait surface, `thread_rng` →
   `rng`, distribution imports).

**Not migrating to edition 2024 now.** It is the root cause but drags in an
unrelated breaking-change surface. The rename makes us forward-compatible;
defer the edition bump to its own change.

### Risk

The rand 0.9 upgrade **changes the RNG algorithm**, so generated output shifts
for a given seed. Profiles are statistical, so fits remain valid, but any
seed-pinned test expectation (insta snapshots, hard-coded genotype assertions)
must be re-baselined. This is the main risk in A and the reason to land it
before C touches the same tests.

### Acceptance

- `rg '\bgen\b' src/` returns nothing (module and method).
- `cargo test --all-features`, `cargo clippy --all-features -- -D warnings`,
  `cargo fmt --check` all green.
- rust-analyzer reports no syntax errors across `src/bulk/` (manual check).
- Any re-baselined snapshot/test documents the seed change in its commit.

---

## Workstream B — schema (#9, #10 schema items)

### Problem

1. **`fitted.ploidy` is a generation choice, not a measurement** (#9).
   `fit_profile.py`'s `--ploidy` (default 2) is written straight into
   `fitted.ploidy` on both paths and is never derived from source data. It is
   exactly like `payload`, which correctly lives in `Dialed`. The branch's
   headline invariant is "`fitted` contains only values measured from data";
   `ploidy` violates it.
2. **Hand-supplied fields aren't auditable from the JSON** (#9). `phased_rate`
   on `--sites-vcf` is also hand-supplied (a sites-only file has no genotypes),
   but on `--pgen` it *is* measured — so it cannot simply move to `Dialed`. A
   reader of the emitted JSON can't tell a supplied value from a measured one,
   and `Profile::validate` asserts the hand-typed number as an invariant.
3. **A diploid-only constraint is enforced at runtime** (#10). `validate`
   accepts `ploidy >= 1` but `write()` re-rejects `ploidy != 2` for AD/PL
   payloads, because AD/PL are hard-coded diploid in `SampleStats`.
4. **Nothing gates a fresh fit against the Rust validator** (#10). A new fit
   emitting `titv <= 0`, `ploidy 0`, or a NaN bin writes happily and only fails
   later at `include_str!` time.
5. **`CLASS_NAMES` / `--payload` choices hand-mirror the Rust enums** (#10) with
   a "must match exactly" comment and no enforcement.

### Design

**Move `ploidy` `Fitted` → `Dialed`:**
- `src/bulk/profile.rs`: remove `ploidy` from `Fitted`, add to `Dialed`; update
  validation and the `fitted_pairs`/`dialed` accessors.
- Re-emit the three committed `profiles/*.json` to move the key (germline-1kgp,
  germline-1kgp-unphased, somatic-gdc). The somatic one is re-emitted anyway by
  workstream D; the two germline ones are hand-moved or re-emitted.
- `scripts/fit/fit_profile.py`: write `ploidy` under `dialed`.
- Update any test reading `p.fitted.ploidy`.

**Add `provenance.supplied: Vec<String>`:** a list naming every non-measured
field, populated by the fit script per path:
- `--pgen`: `["ploidy"]` (phased_rate *is* measured here).
- `--sites-vcf`: `["ploidy", "phased_rate", "n_samples"]`.

This makes both `ploidy` and `phased_rate` auditable from the JSON alone, which
is the invariant #9 is actually defending. `Profile::validate` uses the list to
scope which fields it asserts as measured invariants vs. accepts as supplied.

**Reject `ploidy != 2` with AD/PL payload at parse time:** move the constraint
from `write()` into `Profile::validate`, making the invalid state
unrepresentable at construction rather than at generation. (Chosen over
deriving PL cardinality from ploidy: that invents a capability nothing has
asked for. Small, honest, and matches the current diploid-only reality.)

**CI gate for fresh fits:** a CI step (and/or a `fit_profile.py` post-write
self-check) that loads each freshly-written profile through `vcfixture`'s
validator before it can be committed. Closes the `include_str!`-time failure gap.

**Enforce `CLASS_NAMES` / `--payload` against the Rust enums:** replace the
"must match exactly" comment with an actual check. Options for the implementer:
a generated constant, a test that parses the Rust enum variants, or a shared
data file. Pick the lightest that fails loudly on drift.

### Acceptance

- `profiles/*.json` carry `ploidy` under `dialed`, and a `provenance.supplied`
  list; `Profile::validate` passes on all three.
- A profile with `ploidy: 3` and an AD/PL payload is rejected by `validate`
  (new test), not only by `write()`.
- A deliberately-broken fit (e.g. `titv: 0`) is caught by the CI gate, not at
  `include_str!` time (new test).
- A drift between `CLASS_NAMES` and the Rust `ClassMix` variants fails a test.

---

## Workstream C — generator perf (#8, #10 perf items)

### Problem

**`Size::Target` is ~7× slower than `--records-per-contig` for similar output**
(#8): `--target-size 8MB` took 226 s vs 34.5 s for `--records-per-contig 3000`.

Root cause, confirmed by reading `resolve_target_counts` and
`measure_compressed_bytes` (`src/bulk/mod.rs`):
- The search runs up to `MAX_ROUNDS = 25`, and **each round's
  `measure_compressed_bytes` generates every contig twice** (a span pass then a
  write pass), so 25 rounds ≈ **50 full generations**, and `write()` adds ~2
  more.
- The convergence rule is one-sided (`if bytes >= target_bytes { return }`)
  with additive growth `extra = ceil(shortfall / bytes_per_record * 1.15) + 1`.
- **`bytes_per_record = total_file_bytes / total_records` includes the header**,
  which at small counts (round 1 = 500/contig) with many samples inflates the
  per-record estimate, so every round *undershoots* and the loop iterates many
  times. The `*1.15` margin is a band-aid over this bias.
- Measurement is **byte-exact**: same `self.seed`, same `build_header`, same
  format/compression/workers, no mid-stream flush — the docs already assert the
  measured count *equals* the final file size. Yet the winning temp file is
  deleted and the final output regenerated from scratch.

**Two per-record allocation hot spots** (#10):
- `to_record_buf` (`src/bulk/generate.rs`): `SampleStats::new` builds GT via
  `map` → `Vec<String>` → `join`, ~3 string allocations per sample per record.
  At 100 MB (265k records × 3202 samples) that is ~2.5B allocations.
- `gen_record` allocates a `Vec<usize>` of `n_alleles` (6,404) every record to
  place a median of 1 alt allele (36% of records are singletons).

**A fitted-but-unreproduced split weight** (#10): `Size::Records` splits by
`density_per_kb`, but the output's actual density is `1/mean(gap)` globally
(`gap_dist` is not per-contig), so fitted per-contig density is never
reproduced. `ContigStat.n_variants` is fitted but never read.

### Design

**Replace the search with two-point linear calibration + at most one
correction, then promote the measured file:**

1. Measure compressed bytes at **two small calibration counts** (e.g. `c1` and
   `c2 = 2·c1` records total, split by the chosen weight). Solve
   `bytes ≈ b0 + k·records` for both intercept `b0` (the header/fixed cost) and
   slope `k` (bytes per record). Fitting the intercept **structurally removes
   the header-contamination bias** that drove the round count.
2. Compute the target count directly: `records = ceil((target − b0) / k)`.
3. Do **at most one corrective round**: measure at the computed count; if under
   target, top up once using the same slope. (Compression ratio is stable
   across records from one profile, so this lands within the existing overshoot
   tolerance.)
4. **Promote the winning measurement instead of regenerating.** Because the
   write is byte-exact, the final corrective measurement's temp file *is* the
   final output. Restructure so the summary (`summary.observe`) is computed
   during that write pass and the temp file is moved to the destination
   (`fs::rename`, falling back to copy across filesystems), rather than
   `write()` regenerating from the returned counts.

This turns ~50 generations into **~3–4** (two calibration + one corrective +
promote), targeting #8's "seconds to low-tens-of-seconds" expectation.

**Fix the two allocation hot spots:**
- GT: write into a reused `String::with_capacity(2 * ploidy)` (or a small
  stack buffer) instead of `Vec<String>` + `join`.
- Alt placement: **rejection-sample** the alt positions when `ac << n_alleles`
  instead of shuffling all `n_alleles` indices. Keep the shuffle path for the
  dense case.

**Switch the `Size::Records` split weight `density_per_kb` → `n_variants`:**
reproduces the source's per-contig variant distribution (the thing a user
splitting by contig expects) and stops the MT outlier (350/kb, 12×) skewing a
split it was never meant to drive. This makes the fitted-but-unread
`n_variants` field load-bearing. Note the interaction with workstream D's
`n_variants` re-fit — the values must be correct genome-wide.

### Risk / constraints

- Promotion via `fs::rename` must handle the temp file being on a different
  filesystem than the destination (`TMPDIR` vs output dir); fall back to copy.
- The BCF CSI index is written alongside; promotion must carry or regenerate
  the `.csi` for the final path.
- Calibration counts must be large enough that `k` is stable but small enough
  to stay fast; pick from measurement, don't hard-code blindly.

### Acceptance

- `tests/bulk.rs`'s existing overshoot bound (`MAX_OVERSHOOT_FRACTION = 0.25`,
  one-sided `got >= target`) still passes.
- `--target-size 8MB` wall time is benchmarked before/after against the 226 s
  baseline and lands in low-tens-of-seconds.
- The promoted output is byte-identical to what a from-scratch `write()` of the
  same counts would produce (regression test: compare against a
  known-good hash or a re-generated file).
- GT and alt-placement changes leave `truth()` and rendered output unchanged
  (existing generation tests pass; the RNG stream for alt placement must be
  preserved or the change explicitly re-baselined).

---

## Workstream D — fit-script memory (#7, #10 histogram bug)

### Problem

`fit_profile.py --pgen` peak RSS scales with source rows: 75.5M rows → 6.7 GB
(OK), 348M rows → 25.4 GB (**OOM-killed** at the 32 GB allocation). The shipped
`somatic-gdc.json` describes only chr21+22 because a genome-wide fit OOMs.

**The issue's stated lead is wrong.** There are no joins in `fit_profile.py`;
`.acount`/`.vmiss` are each scanned independently and reduced to bounded frames.
Every `.collect()` already returns a bounded result. Measurement located the
**three actual O(rows) offenders** (all confirmed on the 348M-row source):

| stage (genome-wide) | current peak RSS | wall |
|---|---|---|
| **whole `compute_pvar_stats`, current code** | **OOM / 32.2 GB killed** | — |
| `_gap_bins_lazy` alone (`sort` + `.over()` window) | ~44–57 B/row over a 6.2 GB floor → ~20 GB+ | slow |
| `_titv_lazy` alone (`concat_str().is_in(...)`) | **16.7 GB** | 50 s |
| `pl.collect_all` CSE cache (6 plans together) | +20 GB over sequential | — |

A **fixed ~6.2 GB scan floor** exists because `read_pvar`'s `scan_csv` uses
`comment_prefix="##"`, which forces polars to materialize the whole 5.9 GB CSV
(counting rows alone costs 6.2 GB; projecting only `POS` still costs 6.2 GB —
projection pushdown does not help). The `--sites-vcf` path is lean **not**
because of bcftools/narrow-TSV but because its `scan_csv` has **no
`comment_prefix`** (bcftools output has no `##` lines). This is the true reason
the two paths differ, and it corrects the issue's framing.

**Single-bin histogram bug** (#10): `_bucket_index_expr` disagrees with
`numpy.histogram` for single-bin histograms — the closed-last-bin clause is
skipped at `n_bins == 1`, dropping a value equal to `edges[-1]`; reachable via
`_sfs_edges(1)`. The docstring claims exact numpy parity.

### Design — four fixes, each measured

1. **`read_pvar` scan: drop `comment_prefix`, use `skip_lines`.** Count the
   leading `##` lines with a cheap byte loop, then
   `scan_csv(..., skip_lines=n)` with no `comment_prefix`. Removes the
   materialization forced by comment-scanning; also ~3× faster on the scan.
   (Also drop the unused `ID` column from the select — but note honestly: in
   this file `ID` is `.`, so this is tidiness, ~0 memory. The wide columns are
   REF/ALT, which are needed.)
2. **`_gap_bins_lazy`: drop the sort, replace the window.** The pvar is already
   coordinate-sorted (plink2 emits it so). Replace
   `sort(["CHROM","POS"]).select(POS.diff().over("CHROM"))` with
   `POS.diff()` masked by `CHROM == CHROM.shift(1)`. **Verified bit-identical**
   on chr21 (n=4,107,307, total=41,668,608) and genome-wide (318,372,446 gaps).
   Add a cheap streaming **monotonicity precondition check** so the dropped-sort
   assumption is asserted, not silent.
3. **`_titv_lazy`: replace `is_in` with direct base comparisons.** Replace
   `concat_str([REF,ALT]).is_in(TRANSITION_PAIRS)` with the four transition
   pairs as explicit `(REF==a)&(ALT==b)` disjunctions. **Verified identical**
   `n_ts = 163,208,320`; 16.7 GB → 6.4 GB.
4. **`compute_pvar_stats`: collect plans sequentially, not `pl.collect_all`.**
   CSE caches the shared 348M-row scan/explode output, costing +20 GB and
   buying nothing — re-scanning from warm page cache is *leaner and faster*.
   Correct the docstring (lines ~671–673) which currently asserts the opposite.

**Measured end-to-end result of fixes 2+3+4, genome-wide:**

| variant | peak RSS | wall |
|---|---|---|
| current | OOM (killed) | — |
| gap + titv, `collect_all` CSE on | 26.4 GB | 2:57 |
| gap + titv, CSE off | 18.2 GB | 3:00 |
| **gap + titv + sequential collect** | **6.4 GB** | **2:30** |

The genome-wide somatic fit goes from OOM to **6.4 GB, 2.5 min** — a 4× memory
reduction that also runs faster. All three histograms and scalars are
bit-identical to the current code's output on the contigs where the current
code can run.

5. **Fix `_bucket_index_expr` single-bin bug**: include the closed-last-bin
   clause at `n_bins == 1` so a value equal to `edges[-1]` is counted, matching
   `numpy.histogram`. Regression test against numpy for `n_bins == 1`.

### Re-fit deliverable

After B (schema) and D (memory) land, **re-run the genome-wide somatic fit and
commit the corrected `somatic-gdc.json`**, replacing the chr21+22 stopgap. This
is the payoff #7 exists to deliver.

- `n_variants_source` changes ~7.9M → ~348M; `contigs` grows from 2 to ~25.
- The re-emitted profile carries B's new schema (`ploidy` under `dialed`,
  `provenance.supplied`).
- `Profile::builtin("somatic-gdc")` must validate and round-trip.
- Because C switches the `Size::Records` split to `n_variants`, the genome-wide
  per-contig `n_variants` must be correct — this re-fit provides them.

### Acceptance

- Genome-wide `--pgen` fit of the somatic source completes **under the 32 GB
  allocation** (measured peak ≤ ~7 GB) and is benchmarked.
- Fixed `_gap_bins_lazy` / `_titv_lazy` outputs are asserted bit-identical to
  the current code on a small contig subset (unit test on real or synthetic
  sorted input).
- `_bucket_index_expr` matches `numpy.histogram` at `n_bins == 1` (regression
  test).
- `profiles/somatic-gdc.json` describes the whole genome and validates via the
  B CI gate.

---

## Cross-cutting acceptance

- All Rust: `cargo test --all-features`, `cargo clippy --all-features -- -D
  warnings`, `cargo fmt --check` green.
- Fit script: `pixi run -e fit test-fit` (and `test-fidelity` where relevant)
  green.
- prek hooks installed and passing before any commit (per project convention).
- The four workstreams land as reviewable commits grouped by workstream;
  A ∥ B merge before C ∥ D begin.

## Out of scope

- Edition-2024 migration (A makes us forward-compatible; the bump is deferred).
- Making AD/PL genuinely ploidy-generic (B rejects non-diploid instead).
- Re-fitting the two germline profiles from scratch (only their `ploidy` key
  moves; their fits are unchanged).
- Any change to the generation *statistics* beyond what the allocation and
  split-weight fixes necessarily touch.
