# Task 8 report: extraction script (`scripts/fit/fit_profile.py`)

## Status: DONE

## Commit
`7a140c2932b9d62f8b1c39468a5b3d8217048083` on branch `bulk-task-8`
("feat(fit): add profile extraction script for pgen sources")

## What was built

- `pixi.toml`: added a `fit` feature/environment (`python`, `polars`,
  `numpy`, `plink2`, `pytest`), restricted to `platforms = ["linux-64"]`
  because `plink2`'s only bioconda build is linux-64 (no osx-arm64 build
  exists — confirmed via `pixi search -c bioconda plink2`). Added the
  `bioconda` channel (needed for `plink2`). Tasks: `fit` and `test-fit`.
  `pixi.lock` regenerated via `pixi install -e fit` (517 additive lines,
  no changes to other environments).
- `scripts/fit/fit_profile.py`: implements every function the brief
  specifies —
  - `histogram(values, edges)` — `numpy.histogram`-based binning,
    normalized weights, validates edges are strictly increasing.
  - `class_mix_from_counts(counts)` — normalizes the 6 `ClassMix` fields.
  - `classify(ref, alt)` — scalar reference classifier (symbolic / snp /
    mnp / insertion / deletion / complex), plus `_classify_expr`, a
    vectorized polars expression doing the same classification over an
    entire pvar column (a Python loop over 10-100M+ pvar rows would be a
    non-starter — matches the "vectorized over Python loops" coding
    principle).
  - `read_pvar(path)` — `polars.scan_csv` (lazy) + streaming `.collect()`;
    handles `.pvar`/`.pvar.zst` (zstd decoded natively by polars from the
    path extension, verified empirically); `comment_prefix="##"` skips
    `##` metadata lines without swallowing the `#CHROM` header row.
  - `fit_sfs(pgen_prefix)` — shells out to `plink2 --freq counts`, reads
    `.acount`, returns `ALT_CTS`.
  - `fit_missing_rate` / `fit_phased_rate` — additional helpers (beyond
    the brief's explicit list) needed to populate `missing_rate` and
    `phased_rate`, which the brief's `build_profile` signature takes as
    plain scalars but doesn't say how to compute. `fit_missing_rate` uses
    `plink2 --missing`; `fit_phased_rate` estimates the phased-call
    fraction from a bounded (`--phase-sample-mb`, default 1 Mb) VCF export
    window, since pgen has no direct phase-fraction report.
  - `build_profile(...)` — assembles the schema dict; `provenance` is
    always populated from real inputs (source path, sample/variant
    counts, `datetime.date.today()`, `__version__`), never a placeholder.
  - `main()` — `argparse` CLI: `--pgen`, `--name`, `--out`, `--payload`
    (default `gt-only`), `--contigs` (default: all present), plus
    `--ploidy` and `--phase-sample-mb`.
  - SFS edges (`_sfs_edges`): doubles from 1.0 (`1, 2, 4, 8, ...`),
    capped at exactly `2 * n_samples`, so the first bin is **exactly**
    `[1, 2)` for any `n_samples` — the load-bearing detail from the task
    context. Verified: `_sfs_edges(3202)` → edges[0]==1.0,
    edges[1]==2.0, edges[-1]==6404.0.
  - Gap edges: `numpy.geomspace(1, 1e5, 25)` (24 bins). Indel edges: the
    literal `[1,2,3,4,5,6,10,20,50,100,1000]` from the brief.
- `scripts/fit/test_fit_profile.py`: the exact test file from the brief
  (verbatim), unmodified.
- `scripts/fit/README.md`: documents what the script does, the
  fitted-vs-dialed split, and the exact commands to re-fit both committed
  profiles (paths taken from
  `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`):
  - germline: `/carter/users/dlaub/data/1kGP/plink2/hg38.norm` (3,202
    samples, `gt-only`)
  - somatic: `/carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt`
    (16,007 samples, `gt-vaf`; note the source `.pvar` is uncompressed)
- `.gitignore`: added `__pycache__/`, `*.pyc`, `.pytest_cache/` (the repo
  had no Python ignore rules at all).

## TDD: failure then pass

Step 3 (test written, `fit_profile.py` not yet created):

```
$ pixi run -e fit test-fit
E   ModuleNotFoundError: No module named 'fit_profile'
1 error in 0.62s
```

Step 5, after implementation:

```
$ pixi run -e fit test-fit
...                                                                      [100%]
3 passed in 2.18s
```

Final re-run after all fixes below:

```
$ pixi run -e fit test-fit
...                                                                      [100%]
3 passed in 0.44s
```

## Beyond the brief's 3 tests: manual verification (not committed as tests)

The brief's given tests exercise `histogram`/`class_mix_from_counts`/
`build_profile` only with pre-computed arrays — none of the pvar-reading,
polars-vectorized, or plink2-shelling code paths are covered. Since the
brief explicitly forbids running against real `/carter` data, I verified
the rest against synthetic data (plink2 `--dummy` and a small hand-written
VCF converted to pgen via `plink2 --vcf ... --make-pgen`, all in `/tmp`,
never committed):

- Confirmed exact `plink2 --freq counts` (`.acount`: `#CHROM ID REF ALT
  ALT_CTS OBS_CT`) and `plink2 --missing` (`.vmiss`: `#CHROM ID
  MISSING_CT OBS_CT F_MISS`) column names against a real plink2 binary
  (had assumed these; confirmed correct).
- Ran `main()` end-to-end against a hand-built 6-variant/4-sample VCF
  covering every variant class (snp/insertion/deletion/mnp/symbolic),
  producing a full profile JSON.
- **Verified the emitted JSON round-trips through the actual Rust
  `Profile::from_json(...).validate()`** (via a scratch integration test,
  compiled with `cargo test --features bulk`, then deleted before
  committing — not part of the diff) — confirms the schema match is
  exact, not just visually inspected.
- Cross-checked `classify()` (scalar) against `_classify_expr` (vectorized)
  on 8 REF/ALT pairs covering every class — outputs matched.

### Bugs this caught and fixed before committing

1. **`read_pvar` let polars infer numeric-looking `CHROM` values (e.g.
   `"1"`) as an integer dtype.** This silently turned `ContigStat.id`
   into an int and broke every downstream `plink2 --chr` argument built
   from it (`TypeError: expected str, bytes or os.PathLike object, not
   int` from `subprocess.run`). Fixed by adding explicit
   `schema_overrides={"#CHROM": pl.Utf8, ...}` to the `scan_csv` call.
2. **`_indel_lengths` computed `ALT.str.len_chars() - REF.str.len_chars()`
   without casting off polars' `UInt32` dtype for `str.len_chars()`.**
   For deletions (alt shorter than ref) this **silently underflowed**
   (e.g. `1 - 3` wrapped to `4294967294` instead of `-2`), which
   `.abs()` then left as a wildly wrong huge value — this would have
   corrupted the `indel_length` histogram on any real deletion-containing
   data without raising any error. Fixed by casting both sides to
   `Int64` before subtracting.

Both were real correctness bugs a scale-test would have caught only via
downstream Rust histogram-validation failures or, worse, silently wrong
generated data — glad to have caught them here via manual synthetic-data
verification rather than relying solely on the brief's 3 given tests
(which don't exercise this code path at all).

## Deviations from the brief, with justification

1. **Added `numpy` to `feature.fit.dependencies`** (not listed in the
   brief's Step 1 snippet). Required for `numpy.histogram`/`geomspace`,
   which is the natural vectorized way to implement `histogram()` and the
   log-spaced edge helpers; without it the import fails immediately.
2. **Restricted the `fit` feature to `platforms = ["linux-64"]`** (the
   brief's Step 1 snippet has no `[feature.fit]` platforms line).
   `plink2` on bioconda has no `osx-arm64` build (confirmed via `pixi
   search`), so `pixi install -e fit` fails to solve without this on the
   workspace's `osx-arm64` platform.
3. **Added `.gitignore` entries** for `__pycache__/`, `*.pyc`,
   `.pytest_cache/` — the repo had no Python ignore rules at all, and
   running `pixi run -e fit test-fit` generates these. Small,
   uncontroversial hygiene addition (Boy Scout Rule).
4. **Added `fit_missing_rate` and `fit_phased_rate` helpers plus
   `--ploidy`/`--phase-sample-mb` CLI flags** — not named in the brief's
   function list, but required: `build_profile`'s signature (given
   verbatim in the brief's test) takes `missing_rate`/`phased_rate` as
   plain scalar arguments, and something has to compute them from
   `--pgen` for `main()` to be a working CLI (the brief's own Step 4 says
   "assemble ... filling `provenance` with the real source"). Documented
   the phased-rate approach (bounded VCF-export sampling) in code
   comments and the README since pgen has no native phase-fraction
   report and this is the least-bad honest way to measure it rather than
   hand-picking it (which would violate "never write a hand-chosen value
   into fitted").
5. **Working directory / branch note**: this agent's tool sandbox
   (Edit/Write) was initially pinned to a different worktree
   (`.claude/worktrees/bulk-generation`, branch `worktree-bulk-generation`)
   despite the task brief specifying `bulk-task-8`. I reverted the one
   accidental edit made there before discovering the mismatch (clean,
   confirmed via `git status`), then used `EnterWorktree` to switch into
   `.claude/worktrees/bulk-task-8` (branch `bulk-task-8`, confirmed via
   `git branch --show-current` throughout). The Edit/Write tools remained
   incorrectly pinned to `bulk-generation` even after the switch (a tool
   quirk, not a repo issue), so all file writes in `bulk-task-8` were done
   via `Bash` heredocs/`python3 -c` file writes instead, which did operate
   in the correct worktree (verified via `pwd`/`git branch
   --show-current`/`git log` before every write and again at commit time).
   The final commit `7a140c2` is confirmed on branch `bulk-task-8` in
   `.claude/worktrees/bulk-task-8`, with no stray changes left in
   `bulk-generation`.

## Files touched

- `scripts/fit/fit_profile.py` (new)
- `scripts/fit/test_fit_profile.py` (new, verbatim from brief)
- `scripts/fit/README.md` (new)
- `pixi.toml` (modified: `fit` feature/environment)
- `pixi.lock` (regenerated)
- `.gitignore` (modified: Python ignores)

No changes to `src/**` or `profiles/germline-1kgp.json`.

## Fixes: review round 1

Two review findings addressed, both reproduced end-to-end before and after
the fix. Only `scripts/fit/fit_profile.py`, `scripts/fit/test_fit_profile.py`,
and `scripts/fit/README.md` were touched -- no Rust source or
`profiles/germline-1kgp.json` changed.

### Critical: multiallelic pvar records crashed the script

**Root cause.** plink2 pvar retains multiallelic records natively (it does
NOT auto-split them the way `bcftools norm -m-` would). `fit_sfs()` shelled
out to `plink2 --freq counts`, which for a multiallelic site emits
comma-joined `ALT_CTS` aligned with the comma-joined `ALT` column (e.g.
`ALT="G,T"`, `ALT_CTS="1,1"`). That string went straight into
`histogram()`'s `np.asarray(values, dtype=np.float64)`, which raises
`ValueError: could not convert string to float: '1,1'`. Separately,
`classify()`/`_classify_expr` are only defined for a single REF/ALT pair;
run on a raw joined ALT they silently misclassify (`classify("A", "A,G")`
reads as `"insertion"`, not two independent SNP alleles), corrupting
`variant_classes` and `indel_length` on every multiallelic row that didn't
even crash.

**Fix.** Two units of observation now apply, and the code keeps them
separate:

- **Per-RECORD** (`contigs` density, `gap_dist` inputs, `multiallelic_rate`):
  computed on the raw, un-exploded pvar frame, exactly as before -- a
  triallelic site is still one record.
- **Per-ALLELE** (`variant_classes`, `indel_length`, `titv`): a new
  `_explode_alleles(df)` splits `ALT` on `,` and explodes into one row per
  REF/ALT pair (`df.with_columns(pl.col("ALT").str.split(",")).explode(...)`
  -- vectorized polars, no Python row loop, so it still streams over a
  500+ MB pvar). `classify`/`_classify_expr` then run on that per-allele
  frame, so a 2-ALT site now contributes exactly 2 class/indel-length/Ti-Tv
  observations, matching the brief's stated invariant. A new
  `compute_pvar_stats(df)` function bundles this ordering (record-level
  stats first, then explode+classify for allele-level stats) so `main()`
  can't accidentally reorder it and reintroduce the bug.
- On the plink2 side, `fit_sfs()` now reads `ALT_CTS` as a forced `Utf8`
  column and parses it with a new `_split_alt_cts()` helper
  (`Series.str.split(",").explode(...).cast(Float64)`), splitting the
  comma-joined counts into one float per allele, positionally aligned with
  the exploded ALT alleles -- so each ALT allele contributes its own count
  to the `sfs` histogram instead of crashing.

`multiallelic_rate`'s formula (`_multiallelic_rate`) is untouched -- it
already ran on the un-exploded frame and is now explicitly called before
any exploding happens, so a triallelic site still counts as 1 record in
both its numerator and denominator (verified by
`test_multiallelic_rate_counts_records_not_alleles`).

### Important: `histogram()` silently dropped partially out-of-range values

**Root cause.** `histogram()` only raised `ValueError` if *all* values fell
outside `edges`; a partial drop (e.g. `gap_dist` values beyond the 1e5 bp
cap, from real inter-variant gaps over centromeres/telomeres) vanished from
the normalization with zero diagnostic output.

**Fix.** `histogram()` now compares `arr.size` (total input) against
`counts.sum()` (in-range count) and, whenever any values are dropped, emits
a `UserWarning` reporting the exact count and fraction dropped, e.g.:

```
histogram: dropped 1/6 values (16.7%) outside the edge range [1.0, 100.0]
```

The dropping behavior itself is unchanged (values outside `[edges[0],
edges[-1]]` still don't contribute to the weights) -- I judged clipping to
the outermost bins to be a bigger, riskier behavioral change than the
brief asked for, and a loud diagnostic is enough for a human to decide
whether the edge range needs widening (e.g. bumping `GAP_HIGH` past 1e5)
on the next fit. The all-out-of-range case still raises `ValueError`, as
before.

### Reasoning on per-allele semantics

The task asked me to think through what "correct per-allele semantics"
means for each fitted stat, and keep it consistent:

| stat | unit of observation | why |
|---|---|---|
| `contigs` (density) | record | density is about how many *sites* occupy a span of genome, not how many alleles they carry |
| `gap_dist` | record | inter-variant gap is a property of two adjacent *positions*; exploding first would fabricate spurious zero-length gaps between alleles of the same site |
| `multiallelic_rate` | record | it's asking "what fraction of sites are multiallelic", which is definitionally a per-site question -- exploding first would make every multiallelic site's rate contribution >1, which no longer means "fraction of sites" |
| `variant_classes` | allele | a triallelic SNP/indel site really does draw its two ALT alleles from up to two different classes (e.g. one SNP ALT and one indel ALT is possible) -- classification is fundamentally a REF/ALT-pair property |
| `indel_length` | allele | same reasoning -- each ALT allele has its own length delta from REF |
| `titv` | allele | each SNP allele is an independent transition/transversion draw; a triallelic SNP with one transition and one transversion ALT should contribute one observation to each bucket, not be excluded or double/half-counted |
| `sfs` (via `ALT_CTS`) | allele | the site-frequency spectrum is inherently per-allele already (each ALT allele has its own frequency in the cohort) -- this was the crashing half of the bug |

### Covering tests

Added to `scripts/fit/test_fit_profile.py` (20 tests total, up from 3; all
pass with no `/carter` access and no network):

- `test_split_alt_cts_parses_plink2s_comma_joined_column` -- direct
  regression test for the exact string (`"1,1"`) that crashed
  `np.asarray(..., dtype=float)` before the fix.
- `test_fit_sfs_regression_comma_joined_alt_cts_no_longer_crashes` --
  exercises `fit_sfs()`'s `.acount`-parsing path with a mocked plink2
  subprocess (no plink2 binary needed), feeding it exactly the
  `ALT="G,T"`/`ALT_CTS="2,1"` shape plink2 emits for a real multiallelic
  site.
- `test_multiallelic_alt_corrupts_scalar_classify_but_not_the_split_pipeline`
  -- documents the reviewer's exact `classify("A", "A,G") == "insertion"`
  finding, then proves the real pipeline (`_explode_alleles` +
  `_classify_df`) classifies the same alleles correctly.
- `test_explode_alleles_splits_multiallelic_and_is_a_no_op_for_biallelic`
  -- unit test for the new splitting primitive.
- `test_multiallelic_rate_counts_records_not_alleles` -- proves a
  triallelic site (2 ALTs) contributes exactly 1, not 2, to
  `multiallelic_rate`.
- `test_compute_pvar_stats_gives_two_class_and_two_indel_observations_per_multiallelic_site`
  -- proves a 2-ALT insertion site contributes 2 class observations and 2
  indel-length observations (lengths 1 and 2) while still counting as 1
  record for `multiallelic_rate` (0.5, not 1.0, in a 2-record fixture) --
  this is the brief's stated invariant, tested directly.
- `test_histogram_warns_when_some_values_are_out_of_range` /
  `test_histogram_raises_when_all_values_are_out_of_range` /
  `test_histogram_does_not_warn_when_all_values_in_range` -- the
  `histogram()` fix, including the reviewer's exact repro values
  (`[1,1,2,5,50,99999]`, edges `[1,2,10,100]`), asserting the warning text
  contains the count and percentage dropped.
- `test_sfs_edges_first_bin_is_exactly_one_two` (parametrized over
  `n in {1, 2, 3, 10, 3202, 16007}`) -- explicit regression coverage for
  the singleton-bin invariant this task must not break.
- `test_end_to_end_multiallelic_pgen_does_not_crash_and_has_per_allele_semantics`
  and `test_fit_sfs_preserves_positional_alignment_between_alt_and_alt_cts`
  -- the reviewer's gold-standard repro: a synthetic VCF with a genuine
  `REF=A ALT=G,T` record, converted with `plink2 --vcf ... --make-pgen`
  (confirmed empirically that plink2 keeps it multiallelic, no special
  flags needed), run through the real `main()` CLI end-to-end. Both are
  `@pytest.mark.skipif(not PLINK2_AVAILABLE, ...)` (checked via
  `shutil.which("plink2")`) so CI (which lacks plink2) skips them cleanly
  instead of failing. The second test specifically checks *positional*
  alignment (rs1's two ALT alleles have different counts, G=2 vs T=1) --
  a stronger check than a symmetric `"1,1"` case, which wouldn't catch a
  bug that mixed up which count belongs to which allele.

`_replica_validate_profile()` in the test file mirrors every invariant
`Profile::validate`/`Histogram::validate`/`ClassMix::validate` enforce in
`src/bulk/profile.rs` (strictly increasing finite edges, non-negative
finite weights summing > 0, ClassMix summing to 1.0 within 1e-6, rates in
`[0, 1]`, >=1 contig) without depending on the Rust crate, since modifying
Rust source (even temporarily, to add a scratch `#[test]`) is out of scope
for this task -- the permission system in this environment actively
blocked an attempt to append a throwaway verification test to
`src/bulk/profile.rs`, which is the correct behavior here.

### Exact commands and output

```
$ pixi run -e fit test-fit
✨ Pixi task (test-fit in fit): pytest scripts/fit -q
....................                                                     [100%]
=============================== warnings summary ===============================
scripts/fit/test_fit_profile.py::test_histogram_raises_when_all_values_are_out_of_range
  .../test_fit_profile.py:86: UserWarning: histogram: dropped 2/2 values (100.0%) outside the edge range [1.0, 100.0]
    histogram([1000, 2000], edges=[1, 2, 10, 100])
-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
20 passed, 1 warning in 5.46s
```

Reviewer's exact repro, re-run directly against the fixed `histogram()`:

```
$ pixi run -e fit python -c "
from fit_profile import histogram
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    h = histogram([1,1,2,5,50,99999], edges=[1,2,10,100])
    print(h['weights'])
    print(w[0].message)
"
[0.4, 0.4, 0.2]
histogram: dropped 1/6 values (16.7%) outside the edge range [1.0, 100.0]
```

Singleton-bin invariant, re-verified explicitly for every n the reviewer
listed:

```
$ pixi run -e fit python -c "
from fit_profile import _sfs_edges
for n in [1, 2, 3, 10, 3202, 16007]:
    e = _sfs_edges(n)
    print(n, e[0], e[1], e[-1])
"
1 1.0 2.0 2.0
2 1.0 2.0 4.0
3 1.0 2.0 6.0
10 1.0 2.0 20.0
3202 1.0 2.0 6404.0
16007 1.0 2.0 32014.0
```

Every row: `e[0] == 1.0` and `e[1] == 2.0` -- the singleton bin `[1, 2)` is
untouched by this change.

End-to-end sanity run (synthetic 3-variant/3-sample pgen with the
`REF=A ALT=G,T` multiallelic site, via the real CLI, not just the test
suite):

```
$ pixi run -e fit python scripts/fit/fit_profile.py \
    --pgen /tmp/.../prefix --name synthtest --out /tmp/fitout/synth.json --ploidy 2
wrote /tmp/fitout/synth.json
```

emits a schema-shaped profile with `variant_classes.snp == 0.75`,
`variant_classes.deletion == 0.25` (3 SNP + 1 deletion allele observations
across 3 records, one of which is the multiallelic site contributing 2 SNP
alleles), `multiallelic_rate == 0.333...` (1 multiallelic record out of 3),
and `sfs.edges == [1.0, 2.0, 4.0, 6.0]` -- previously this run raised
`ValueError: could not convert string to float: '2,1'` inside `fit_sfs()`.
