# Investigating `vcfixture bulk`'s sub-linear scaling

**Date:** 2026-08-06
**Branch:** `worktree-fix-22-bulk-perf` (PR #28)
**Follows:** `2026-08-06-bulk-parallel-encode-design.md`, `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md`

## Problem

The parallel-encode work (PR #28) recorded an unresolved finding: wall clock
fell 26.103s → 7.502s from 1 to 8 workers, a 3.48x speedup for 8x the workers.
The measurements doc reported this as **43.5% parallel efficiency** with an
Amdahl-implied serial fraction of ~18.6%, and named two untested explanations
(glibc allocator arena contention; rayon + bgzf thread oversubscription)
without distinguishing them.

## The premise was wrong

`nproc` reports 8 and `SLURM_CPUS_PER_TASK=8`, so the efficiency denominator
was taken to be 8 cores. The allocation is not 8 cores:

```
Cpus_allowed_list:      0-3,48-51
cpu0 thread_siblings:   0,48
cpu1 thread_siblings:   1,49
cpu2 thread_siblings:   2,50
cpu3 thread_siblings:   3,51
```

`carter-cn-03` is a 4-socket Xeon E5-4650 v3: 48 physical cores, 2 threads per
core, 96 logical CPUs. Our cgroup holds **4 physical cores** (0–3) and both SMT
siblings of each. `std::thread::available_parallelism()` reports 8 because it
counts logical CPUs under the affinity mask, which is why `workers` defaulted
to 8 and why 8 became the denominator.

Against 4 physical cores, 3.48x is **~87% efficiency**, before crediting SMT
with anything. The reported 43.5% figure is an artifact of the wrong
denominator, not a measurement of the code.

Two consequences follow:

- The **speedup-vs-baseline numbers in PR #28 are unaffected** (2.51x–3.08x,
  avg 2.81x). Baseline and post-change runs used the same allocation, so the
  comparison is sound. Only the efficiency interpretation was wrong.
- The **1-worker baseline was never single-threaded.** `workers=1` sizes the
  bgzf pool to 1 as well, so that run had 1 rayon + 1 bgzf + main on 4 cores,
  overlapping compression with generation. The 3.48x is measured against an
  already-partly-parallel baseline, which understates the true parallel gain.

## Goal

Replace the "sub-linear scaling is unexplained" caveat with a measured
account, and correct the two places that publish the wrong figure: the
measurements doc and the PR #28 body.

Non-goal: optimizing the allocator hot spots. That is issue #26 and stays
filed.

## Design

Four stages, ordered so the cheapest and most-invalidating run first. **Each
stage is a decision gate** — Stage 1 may resolve the question outright, in
which case Stages 3 and 4 are dropped rather than run for completeness.

Reference workload throughout: `VCFIXTURE_BENCH_SAMPLES=2000
VCFIXTURE_BENCH_RECORDS=20000` (80M cells), the same one the original scaling
check used, so new numbers are directly comparable to the recorded ones.

### Stage 0 — Topology (complete)

Recorded above. Its output is a new "Hardware and allocation" section in the
measurements doc, stating the true core count and how `available_parallelism`
relates to it. Every table in that doc gains the context that 8 workers ran on
4 physical cores.

### Stage 1 — Scaling curve

The existing evidence is two points (1 and 8). A curve discriminates between
the candidate explanations by its *shape*.

Sweep `workers ∈ {1, 2, 3, 4, 6, 8}`, 3 repetitions each, reporting min and
median. Then repeat the sweep under `taskset -c 0-3` (physical cores only, SMT
siblings excluded) to price SMT directly rather than assuming a textbook
figure.

Predictions, stated before the run so they can fail:

| Shape observed | Reading |
|---|---|
| Near-linear to 4 (~3.5–3.9x), shallow 4→8 tail | Core count explains it; residual serial fraction is small |
| Clearly sub-linear by 4 (≤3.0x) | A real serial or contention effect remains — proceed to Stage 3 |
| Non-monotonic; 6 or 8 slower than 4 | Oversubscription confirmed — proceed to Stage 2 |

The `taskset -c 0-3` variant at `workers=4` versus the unpinned `workers=4`
isolates whether SMT siblings help or interfere. The unpinned-8 versus
pinned-4 comparison gives SMT's actual contribution on this workload.

**Gate:** if the curve is near-linear to 4 with a shallow SMT tail, the
question is answered — record it, correct the docs, and stop. Stages 3 and 4
do not run.

### Stage 2 — Writer-pool oversubscription

Under 4 physical cores this hypothesis is far stronger than it looked under
the assumed 8: `workers=8` spawns up to 8 rayon threads + 8 bgzf compression
threads + main ≈ 17 threads on 4 cores.

It is testable with **no library API change**, via two comparisons the
existing knobs already allow:

1. **Format ablation.** Re-run the Stage 1 sweep with uncompressed VCF output,
   which removes the bgzf pool entirely. If BCF turns over at 4 workers and
   VCF keeps climbing, the writer pool is implicated.
2. **Turnover in the BCF curve itself.** Oversubscription shows up directly as
   `workers=6`/`8` being no faster, or slower, than `workers=4`.

Only if these confirm oversubscription do we consider decoupling the two pool
sizes. That would be a `BulkSpec` API addition (a writer-worker count
defaulting to `workers`), and it is a separate decision made with data in
hand — not part of this investigation's committed scope.

### Stage 3 — Measure the serial fraction rather than infer it

This stage exists because of a third hypothesis the original write-up never
recorded. `stream_contigs` (`src/bulk/mod.rs:735-791`) is a **barrier per
chunk**:

```rust
let encoded: Vec<...> = pool.install(|| chunk.par_iter().map_init(...).collect());
for item in encoded {
    writer.write_encoded(&bytes)?;   // every rayon thread idle here
    summary.merge_block(id, &bs);
}
```

`collect()` drains the whole chunk before the serial loop starts, so the
`write_encoded` memcpy into the bgzf staging buffer runs with the entire rayon
pool parked. That memcpy is O(bytes), hence O(cells) — precisely the shape
Amdahl's law was pointing at, and a better-founded candidate than either
recorded hypothesis.

Measure it directly: a throwaway patch accumulating `Duration` across the
`pool.install` region and the drain loop separately, printed at the end of the
run. Compare the measured serial share against the Amdahl-implied figure
recomputed for 4 cores.

The patch is **not committed**. It is saved to scratch, its diff is quoted in
the measurements doc so the number is reproducible, and it is reverted before
any commit.

### Stage 4 — Allocator, only if a gap remains

Runs only if Stages 1–3 leave the curve materially unexplained.

The originally-proposed test — comparing glibc self-time share between a
1-worker and an 8-worker profile — is correlational and stays dropped. Use
causal interventions instead:

- `MALLOC_ARENA_MAX=1` versus default. If forcing a single arena sharply
  degrades the 8-worker run, arena parallelism is load-bearing and contention
  is real.
- A `#[global_allocator]` swap to mimalloc in the bench binary only. If
  efficiency moves, the allocator is confirmed *and* we know the remedy.

Either result feeds issue #26 rather than this branch.

## Harness changes

Two small additions to `examples/bulk_bench.rs`, both genuinely reusable:

- `VCFIXTURE_BENCH_REPS` (default 1) — repeat each sweep cell, report min and
  median. The recorded ~10% run-to-run spread is currently unquantified per
  point; a 3.48x figure derived from single shots deserves error bars.
- `VCFIXTURE_BENCH_FORMAT` (default `bcf`) — select output format, so the
  Stage 2 format ablation needs no separate binary.

No changes to `src/` are part of this investigation's committed scope.

## Operational constraints

- Foreground only, no detaching. NFS-blocked processes on this cluster are
  unkillable and drained a node on 2026-07-29.
- `TMPDIR=$CLAUDE_JOB_DIR/tmp` for every bench run. `bulk_bench` writes via
  `env::temp_dir()`, and parallel background jobs share `/tmp`.
- Total runtime is bounded: one sweep variant is ~4 minutes, and at most three
  variants run (unpinned BCF and pinned-to-physical BCF in Stage 1, plus the
  VCF ablation only if Stage 2 is reached), so this stays well under the
  threshold that would require `sbatch`.
- Confirm no bench process survives before reporting.

## Deliverables

1. `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md` — new
   "Hardware and allocation" section; the Hypothesis 2 verdict rewritten with
   the corrected denominator; the scaling-curve table; a resolution of the
   open question, or an explicit statement of what remains open and why.
2. PR #28 body — the "sub-linear scaling is not fully explained" bullet
   replaced with the measured account.
3. `examples/bulk_bench.rs` — the two knobs above.
4. Issue #26 — a comment recording whatever Stage 4 establishes, if it runs.

## Testing

The investigation is measurement work; its correctness gate is that the
committed tree still passes. `cargo test --all-features`, `cargo fmt --check`,
and `cargo clippy --all-features --all-targets -- -D warnings` must be clean
at the final commit, and the Stage 3 instrumentation patch must be verifiably
absent from it.
