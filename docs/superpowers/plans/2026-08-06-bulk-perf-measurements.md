# Bulk generation measurements (issue #22)

Harness: `cargo run --release --example bulk_bench --features bulk`
Machine: 8 cores, 4.18.0-553.36.1.el8_10.x86_64

## Baseline — v0.4.0 (`9e64a0c`), before any change

```
workers=8
 samples   records        cells       secs       s/cell peakRSS_MB
     500      5000      5000000      1.443     2.886e-7       27.2
     500     20000     20000000      5.606     2.803e-7       45.4
    2000      5000     20000000      5.511     2.755e-7       47.9
    2000     20000     80000000     21.746     2.718e-7       78.5
    8000      5000     80000000     21.316     2.664e-7       83.9
    8000     20000    320000000     87.470     2.733e-7      200.9
```

## After the change — parallel encode + gap-only span pass (`a222aae`, `cf18a72`)

Harness: `cargo run --release --example bulk_bench --features bulk` (worker-count control added in `7bb2f6d`)
Machine: same box (8 cores, 4.18.0-553.36.1.el8_10.x86_64)

```
workers=8
 samples   records        cells       secs       s/cell peakRSS_MB
     500      5000      5000000      0.576     1.152e-7       27.8
     500     20000     20000000      1.974     9.868e-8       40.1
    2000      5000     20000000      2.083     1.042e-7       50.0
    2000     20000     80000000      7.049     8.812e-8       77.0
    8000      5000     80000000      7.615     9.519e-8      117.3
    8000     20000    320000000     29.438     9.199e-8      140.3
```

(`real 0m58.317s`, `user 4m44.372s`, `sys 0m1.844s` for the full sweep, including a ~9.5s incremental compile.)

### Comparison to the v0.4.0 baseline

| samples | records | baseline s/cell | new s/cell | speedup | RSS before (MB) | RSS after (MB) |
|---|---|---|---|---|---|---|
| 500 | 5000 | 2.886e-7 | 1.152e-7 | 2.51x | 27.2 | 27.8 |
| 500 | 20000 | 2.803e-7 | 9.868e-8 | 2.84x | 45.4 | 40.1 |
| 2000 | 5000 | 2.755e-7 | 1.042e-7 | 2.65x | 47.9 | 50.0 |
| 2000 | 20000 | 2.718e-7 | 8.812e-8 | 3.08x | 78.5 | 77.0 |
| 8000 | 5000 | 2.664e-7 | 9.519e-8 | 2.80x | 83.9 | 117.3 |
| 8000 | 20000 | 2.733e-7 | 9.199e-8 | 2.97x | 200.9 | 140.3 |

Speedup range 2.51x–3.08x, average ~2.81x. This clears the "wins" bar plainly — not a 1.1x rounding-error result. RSS falls at the largest workload (200.9 → 140.3 MB, the redundant generation pass no longer holds a second full copy of record state) but *rises* at 8000×5000 (83.9 → 117.3 MB): the parallel fan-out now holds several blocks' `RecordBuf`s and encoded byte buffers in flight concurrently across workers, instead of one block's data alive at a time in the old serial encode stage. Net RSS stays well under 1 order of magnitude change either way and the largest workload is the best RSS case, so this is not a concern, just a real trade-off worth naming.

All (samples, records) points, workers, and workload sizes above are **single-shot measurements**, not averaged over repeated runs. The same (2000, 20000) workload was independently re-measured three times over the course of this task (7.049s in the full sweep, 7.502s in the Step 4 scaling run, 7.762s in the Step 5 profiling run) — a ~10% spread. The overall win is far larger than this noise, so no conclusion above changes, but individual row speedups should be read as accurate to roughly ±10%, not to the three significant figures the table prints.

## Hardware and allocation

```
$ hostname
carter-cn-03

$ nproc
8

$ nproc --all
96

$ grep Cpus_allowed_list /proc/self/status
Cpus_allowed_list:	0-3,48-51

$ cat /sys/devices/system/cpu/cpu0/topology/thread_siblings_list
0,48

$ lscpu | rg -i 'thread|core|socket|model name'
Thread(s) per core:  2
Core(s) per socket:  12
Socket(s):           4
Model name:          Intel(R) Xeon(R) CPU E5-4650 v3 @ 2.10GHz
```

> Every measurement in this document ran inside a Slurm allocation of
> `Cpus_allowed_list: 0-3,48-51` on `carter-cn-03` — **4 physical cores** plus
> their SMT siblings, not 8 cores. `nproc` and
> `std::thread::available_parallelism()` both report 8 because they count
> logical CPUs, which is why `workers` defaulted to 8. Parallel efficiency must
> be computed against 4, not 8.

### Hypothesis 1 — "removing the redundant generation pass roughly halves generation CPU"

**Measured directly with an isolated A/B comparison**, using the `VCFIXTURE_BENCH_WORKERS` knob added in this task to factor parallelism out entirely. A throwaway `git worktree` was checked out at the pre-change baseline commit `9e64a0c` (which predates the `bulk_bench` harness), the current `examples/bulk_bench.rs` was copied into it unmodified (it uses only public `BulkSpec`/`Payload`/`Profile`/`Size` API that is identical at `9e64a0c`), built `--release`, and run once at `WORKERS=1` for the same workload used in the Hypothesis 2 scaling check:

```
OLD code (9e64a0c), WORKERS=1: 34.179s  (4.272e-7 s/cell)
NEW code (HEAD),    WORKERS=1: 26.103s  (3.263e-7 s/cell)
```

At 1 worker there is no parallelism benefit to factor out — this isolates whatever the pipeline-shape change (gap-only span pass + moving encode/summary into a single interleaved pass) is worth on its own, end to end. Result: **1.31x speedup, a 23.6% reduction in end-to-end single-threaded wall clock** (26.103s / 34.179s = 0.7637; 34.179 − 26.103 = 8.076s saved).

This A/B is one run per side, not averaged, on top of the ~10% run-to-run spread already established above for this exact workload — so the ~1.31x/23.6% figures should be read as accurate in direction but not pinned to three significant figures. Even at a worst-case reading (old run 5% low, new run 5% high: 32.470s vs. 27.408s) the ratio is still ~1.18x, safely above 1x — the *direction* of the result is not noise, only its precise magnitude is uncertain by a few points. Also worth flagging: the two binaries were run at the same `seed=42`, but `block_rng` changed from `(seed, block_idx)` to `(seed, block_idx, stream)` in an earlier task in this same body of work, so the old and new binaries generate different variant *content* at that seed — same distributions, same record and cell counts (verified: both `n_records_total()` values and both `cells` denominators are 80,000,000), but not byte-identical records. At 20k records this is very unlikely to move timing, but it is an unstated difference between the two sides worth naming for completeness.

**Read against what H1 actually claims — generation CPU, not end-to-end wall clock:** H1 says removing the redundant pass roughly halves *generation CPU*, not that it halves total single-threaded time. At `9e64a0c`, the spans loop calls `generate_contig` to lay out the file, and the write loop calls it again to produce records — two identical full-genotype generation passes. Deleting one of two identical passes halves generation CPU **by construction**; that mechanism is not in question.

The A/B's 8.076s delta is the cost of removing *one* of those two passes (the pipeline-shape change also restructures encode/summary into the same interleaved pass, per the methodology paragraph above — the A/B cannot separate the two effects, so 8.076s is an **upper bound** on the removed generation pass's cost, not a clean measurement of it; the true generation-pass cost may be somewhat lower if any of the 8.076s is attributable to the encode/summary restructuring instead). Treating it as the removed pass's cost for the sake of a sanity check: the *old* binary paid for **two** such passes, so total generation work in the old run was ≈ 2 × 8.076s = 16.152s of the 34.179s total — **≈47%** of old single-threaded cost, an upper-bound estimate. In the *new* binary, one pass remains: ≈8.076s of 26.103s — **≈31%**, also an upper-bound estimate. Non-generation work (encoding, header/CSI, I/O) is the complementary ≈53% of the old total, not the ≈76% an earlier draft of this paragraph stated — that draft mistakenly treated the 8.076s delta as *all* of generation rather than as *one of two* generation passes, which halved the implied generation share and inflated the non-generation share to match.

With the corrected arithmetic, the observed 1.31x end-to-end result is not a shortfall against H1 — it is exactly what Amdahl's law predicts from halving a ≈47%-of-total component: removing one full generation pass out of two (34.179s − 8.076s = 26.103s) is definitionally what the measured end-to-end delta is, so the two numbers agreeing is not independent confirmation, just consistent bookkeeping. The substantive point is the ratio it implies: you cannot get a 2x *end-to-end* win by halving something that was under half the total to begin with, even when the halving itself is complete and exact.

Verdict: **held, as a claim about generation CPU specifically — not as a claim about end-to-end wall clock.** The mechanism H1 describes (deleting one of two identical generation passes) did occur and does halve generation CPU by construction; the corrected arithmetic above (generation ≈47% of old single-threaded cost, upper bound) makes the observed 1.31x end-to-end result the expected Amdahl's-law consequence of that halving, not evidence against it. What this A/B does rule out is reading H1 as implying a 2x *end-to-end* single-threaded speedup: generation was under half of total single-threaded cost, so halving it nets a real but smaller ~1.31x end-to-end improvement. Do not read the ~1.31x figure as H1 failing, and do not read it as licensing a claim of a 2x end-to-end win either — both would misstate what was measured.

### Hypothesis 2 — "moving encode/summary into the fan-out removes the serial O(cells) stage, so wall clock scales with `--threads`"

Measured directly at fixed workload (2000 samples × 20000 records, 80M cells):

```
VCFIXTURE_BENCH_WORKERS=1: 26.103s  (3.263e-7 s/cell)
VCFIXTURE_BENCH_WORKERS=8:  7.502s  (9.377e-8 s/cell)
```

Speedup = 3.48x for 8x the workers (43.5% parallel efficiency). Wall clock is **not flat** — it falls substantially, confirming the serial `O(cells)` stage is no longer the dominant cost. But scaling is clearly sub-linear, not the full 8x a purely-parallel fan-out with a cheap serial tail would predict.

Verdict: **held, partially.** The serial `O(cells)` encode stage that dominated before is gone — s/cell nearly triples in efficiency from 1→8 workers, it does not stay flat. But claiming full linear scaling would be dishonest: at ~44% parallel efficiency, something is still capping speedup below 8x (implied serial fraction ≈ 18.6% by Amdahl's law on the 3.48x/8-worker numbers).

The profile (below) shows glibc allocator-family symbols at ~47% of self time — the largest measured bucket, and large enough to plausibly explain a meaningful chunk of the sub-linear scaling. But the profile supports allocator **overhead**, not allocator **contention**: the `perf` event used is `cycles:u` (user-space cycles only), which cannot see kernel-side lock/futex wait time, so there is no direct evidence of threads blocking on the allocator's arena lock specifically — that would require a 1-worker profile to compare malloc self-time share against (not collected here) or a lock/futex-aware event, neither of which was gathered.

A second, untested, and at least equally plausible explanation is plain **thread oversubscription**: `workers` sizes both the rayon pool (`rayon::ThreadPoolBuilder::new().num_threads(self.workers.get())`, `src/bulk/mod.rs`) and the bgzf multithreaded writer's compression pool (`bgzf::io::multithreaded_writer::Builder::default().set_worker_count(workers)`, `src/bulk/writer.rs:149`) *independently* from the same `workers` value. At `workers=8` on this 8-core box, that is up to 8 rayon worker threads plus up to 8 bgzf compression threads plus the main thread — as many as 17 threads contending for 8 cores, which alone would depress parallel efficiency well below 100% regardless of any allocator effect.

Honest statement: **the profile establishes that glibc allocator overhead is the largest measured self-time bucket (~47%), not that it is the cause of the sub-linear scaling.** Both "allocator arena-lock contention" and "rayon+bgzf thread oversubscription" are plausible, untested explanations for the ~44% parallel efficiency ceiling, and this measurement cannot distinguish between them.

### Allocator interventions

Task 3 instrumented `stream_contigs`'s chunk barrier directly and exonerated
it: measured `serial_frac` (2.9%-3.1%) sits well under the Amdahl-implied
serial fraction (6.3% BCF, 10.0% VCF) needed to explain the measured `S(4)`.
The gap between measured and ideal-linear speedup is therefore still
unexplained, and the allocator — the other untested candidate from Hypothesis
2 above — is the next one to test. Comparing glibc self-time share between a
1-worker and an 8-worker profile would be **correlational**: a higher share
at 8 workers is equally consistent with "more allocation work" and with
"lock contention," and no such comparison was even collected (Hypothesis 2
notes only a single 8-worker profile was taken). The two probes below are
**interventions** instead — each changes the allocator's behavior and
observes whether `S(w)` changes, which a self-time comparison cannot do.

**Step 1 — arena-count intervention.** Reference workload (BCF, unpinned),
`workers=8`, `reps=3`, `MALLOC_ARENA_MAX=1` versus glibc's default, back to
back, same freshly rebuilt `bulk_bench` binary:

```
MALLOC_ARENA_MAX=1 TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

| arena config | min_s | med_s | max_s |
|---|---|---|---|
| `MALLOC_ARENA_MAX=1` | 77.558 | 77.678 | 78.908 |
| default | 8.467 | 9.021 | 9.460 |

Forcing a single arena makes `workers=8` **9.16x slower** (`77.558 / 8.467`)
— two orders of magnitude past the ~10% run-to-run spread this document
otherwise treats as noise, so this is unambiguously a real effect. (The
default arm here, 8.467s, reads ~18% above Task 2's originally recorded
BCF-unpinned `workers=8` figure of 7.153s from a separate session — wider
than the ~10% spread quoted elsewhere, and worth flagging as the true
session-to-session noise band on this box, but it changes nothing about the
Step 1 reading: the arena effect is ~9x, dwarfing either baseline by close to
an order of magnitude.)

Per the brief's reading rule, sharply degrading `min_s` means **per-thread
arenas are load-bearing and allocator contention is real** in the sense that
removing glibc's per-thread-arena mitigation causes catastrophic contention.
It does not by itself say whether that mitigation is *complete* under the
default configuration every other measurement in this document was taken
under — Step 2 tests that directly.

**Step 2 — allocator-swap intervention.** `mimalloc = { version = "0.1",
default-features = false }` added to `[dev-dependencies]`, and
`#[global_allocator] static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;`
added to `examples/bulk_bench.rs` after the `use` block — both temporary,
reverted immediately after the runs below (`git checkout -- Cargo.toml
examples/bulk_bench.rs Cargo.lock`, confirmed with a clean `git status
--short`; see the commit this section belongs to). Built cleanly against
crates.io (`mimalloc 0.1.52`, `libmimalloc-sys 0.1.49`, `cc 1.4.0`). Reference
workload, BCF, unpinned, `reps=3`, `workers ∈ {1,4,8}`:

```
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=1 ./target/release/examples/bulk_bench
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=4 ./target/release/examples/bulk_bench
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

| workers | mimalloc min_s | glibc min_s (Task 2, BCF-unpinned) | speedup |
|---|---|---|---|
| 1 | 15.416 | 27.577 | 1.79x |
| 4 | 4.736 | 8.205 | 1.73x |
| 8 | 4.375 | 7.153 | 1.63x |

mimalloc's own scaling, `S(w) = mimalloc_min_s(1) / mimalloc_min_s(w)`,
against the same BCF-unpinned series used throughout this document:

| workers | mimalloc `S(w)` | glibc `S(w)` (Task 2) | gap |
|---|---|---|---|
| 4 | 3.26 (`15.416/4.736`) | 3.36 | 3.2% lower |
| 8 | 3.52 (`15.416/4.375`) | 3.86 | 8.6% lower |

Two separate findings here, and they point in different directions. mimalloc
makes every worker count **uniformly faster in absolute terms** — 1.63x to
1.79x faster wall clock at `w=1`, `4`, and `8` alike, with no trend toward a
*larger* speedup at higher worker counts. That is consistent with the
profile's ~47% allocator self-time finding: allocation cost is a large, real,
per-call overhead, and a faster allocator removes a roughly constant fraction
of it independent of thread count — exactly the "swap that makes everything
uniformly faster without changing the shape of the curve" case, not one that
changes scaling. But mimalloc's *own* `S(w)` is not higher than glibc's: 3.2%
lower at `w=4` (inside the noise band) and 8.6% lower at `w=8` — smaller than
the ~18% session-to-session spread Step 1's own default-arm reading showed
for the identical configuration, so this cannot be read as a confirmed
degradation either, only as "not better." If arena-lock contention were
capping scaling under the default allocator, swapping to mimalloc's
fundamentally different, contention-avoiding design (thread-local heaps, no
shared arena lock) should have *raised* `S(w)`. It did not.

**Verdict.** Combining both interventions: allocator cost is large and real
— Step 2 confirms this causally, not just as a profile bucket, since a
different allocator saves 1.6x-1.8x of wall time at every worker count — but
it is **not the cause of the sub-linear scaling curves** in the Scaling
curves section below. Step 1 shows glibc's per-thread arenas are load-bearing
scaffolding that, if removed, produces catastrophic (9.16x) contention; but
that scaffolding is exactly what the default configuration under which every
scaling curve in this document was measured already has, and Step 2 shows a
qualitatively different, more contention-resistant allocator does not unlock
additional scaling headroom on top of it. At this thread count, allocator
behavior reads as **pure per-thread overhead, not a scaling limiter** — the
concern belongs to issue #26 (allocation count/size), not to this
investigation's search for what caps `S(w)` below linear. The gap between
measured `S(4)`/`S(8)` and the ideal 4x/8x — Amdahl-implied serial fractions
of 6.3% (BCF) to 10.0% (VCF) against a measured chunk-barrier serial fraction
of only ~3% (Task 3) — **remains open** after this task, with both
originally-proposed candidates (chunk-barrier serial drain, allocator
contention) now tested and exonerated.

### Scaling curves

Fixed workload throughout (2000 samples × 20000 records, 80M cells, `seed=42`), swept over `VCFIXTURE_BENCH_WORKERS ∈ {1,2,3,4,6,8}`, 3 reps per cell (`VCFIXTURE_BENCH_REPS=3`), reporting `min_s`/`med_s`/`max_s` per cell. `S(w) = min_s(1) / min_s(w)`; efficiency is `S(w) / min(w, 4)` — against **4 physical cores**, per the hardware section above (so `w=6` and `w=8` both divide by 4, not by `w`).

**BCF, unpinned** (`Cpus_allowed_list: 0-3,48-51`, all 8 logical CPUs available):

| workers | min_s | med_s | max_s | S(w) | S(w)/min(w,4) |
|---|---|---|---|---|---|
| 1 | 27.577 | 31.642 | 36.312 | 1.00 | 1.00 |
| 2 | 13.309 | 14.042 | 18.509 | 2.07 | 1.04 |
| 3 | 10.614 | 12.355 | 12.526 | 2.60 | 0.87 |
| 4 | 8.205 | 8.305 | 8.374 | 3.36 | 0.84 |
| 6 | 7.580 | 7.650 | 8.024 | 3.64 | 0.91 |
| 8 | 7.153 | 7.363 | 7.802 | 3.86 | 0.96 |

**BCF, pinned to physical cores only** (`taskset -c 0-3`, excludes the 48-51 SMT siblings):

| workers | min_s | med_s | max_s | S(w) | S(w)/min(w,4) |
|---|---|---|---|---|---|
| 1 | 23.996 | 24.468 | 25.856 | 1.00 | 1.00 |
| 2 | 12.916 | 12.981 | 13.085 | 1.86 | 0.93 |
| 3 | 9.305 | 11.030 | 11.223 | 2.58 | 0.86 |
| 4 | 7.565 | 7.709 | 8.083 | 3.17 | 0.79 |
| 6 | 7.741 | 7.876 | 7.930 | 3.10 | 0.77 |
| 8 | 7.303 | 7.312 | 7.457 | 3.29 | 0.82 |

**VCF, unpinned, no bgzf pool** (uncompressed `--format vcf`, so `workers` sizes only the rayon pool, not a second compression pool):

| workers | min_s | med_s | max_s | S(w) | S(w)/min(w,4) |
|---|---|---|---|---|---|
| 1 | 18.529 | 18.718 | 19.018 | 1.00 | 1.00 |
| 2 | 10.615 | 12.253 | 17.061 | 1.75 | 0.87 |
| 3 | 7.287 | 7.460 | 7.962 | 2.54 | 0.85 |
| 4 | 6.026 | 6.169 | 6.218 | 3.07 | 0.77 |
| 6 | 5.906 | 6.082 | 6.095 | 3.14 | 0.78 |
| 8 | 5.521 | 5.794 | 5.947 | 3.36 | 0.84 |

Two of the three curves rise monotonically from `w=1` through `w=8`; the BCF-pinned curve does not. BCF-unpinned reaches 3.36x by `w=4` (84% efficient against 4 physical cores) and keeps climbing more slowly through the SMT range to 3.86x at `w=8` (96% efficient) — consistent with the SMT siblings contributing real, if diminished, throughput on this workload rather than pure oversubscription noise. VCF-unpinned also climbs monotonically throughout, from 3.07x at `w=4` to 3.36x at `w=8`. BCF-pinned tracks BCF-unpinned reasonably closely through `w=4` (within the ~10% single-run spread already established elsewhere in this document for this exact workload, see Hypothesis 1) but shows the oversubscription signature the other two curves don't: `S(6)=3.10` dips slightly below `S(4)=3.17` before recovering to `S(8)=3.29` — running 6 or 8 rayon+bgzf threads on 4 pinned physical cores produces mild, non-monotonic thrashing rather than a clean plateau, expected once the process is confined to fewer cores than it requests threads for.

VCF's absolute speedup is lower than BCF-**unpinned**'s at every worker count `w≥4` (3.07 vs 3.36 at `w=4`, 3.14 vs 3.64 at `w=6`, 3.36 vs 3.86 at `w=8`), and its peak (3.36x at `w=8`) is not materially higher than BCF-unpinned's peak (3.86x) — this is the comparison Gate B's "VCF reaches a materially higher peak speedup than BCF" clause evaluates, using the brief's default BCF-unpinned series, so Gate B does not fire on this evidence. VCF is *not*, however, uniformly below BCF-**pinned**: it trails at `w=4` (3.07 vs 3.17) but overtakes it at `w=6` (3.14 vs 3.10) and `w=8` (3.36 vs 3.29) — VCF's peak is in fact ~2% above BCF-pinned's peak (3.36 vs 3.29). That crossover is the same order of magnitude as the ~10% single-run spread noted above, so it reads as directionally suggestive rather than conclusive proof that uncompressed rayon-only writes outrun bgzf-pinned writes past `w=4` — but it is exactly the comparison Gate B's peak-speedup clause is asking about, so it is worth flagging even though it does not change the verdict (the gate's default series is BCF-unpinned, against which VCF's peak is clearly lower, not higher). VCF's efficiency plateaus at 77%–84% from `w=4` onward (77% at `w=4`, 78% at `w=6`, 84% at `w=8`) rather than climbing toward 100%, indicating a real, if modest, sub-linear component in rayon-only fan-out scaling that is not explained by the bgzf pool at all.

### Serial-fraction measurement

Gate C fired: the VCF-unpinned series gave `S(4) = 3.07 < 3.2`, so a real serial
stage or contention effect is worth chasing directly rather than inferring from
Amdahl's law alone. `stream_contigs` barriers per chunk —
`chunk.par_iter().map_init(...).collect()` fully drains before the serial
`for item in encoded { writer.write_encoded(&bytes)?; ... }` loop starts, so the
`O(bytes)` memcpy into the bgzf staging buffer runs with the entire rayon pool
parked. That is an `O(cells)` serial stage, the shape Amdahl's law points at, and
a candidate neither originally-recorded hypothesis considered. This section
measures that stage's cost directly, rather than inferring it.

**Method:** `stream_contigs` was instrumented (temporarily — reverted
immediately after the two runs below, never committed) to accumulate wall time
spent inside `pool.install(...)` (`t_par`) separately from wall time spent in
the serial drain loop (`t_ser`), across all contigs and chunks, and to print
`serial_frac = t_ser / (t_par + t_ser)` once per run — a fraction of
*instrumented* time (`t_par + t_ser`), not of the run's total wall clock
(`min_s`); see the reconciliation below the results table. The instrumentation
patch, applied against `src/bulk/mod.rs` at the commit this document was
written against, was generated with `git diff -U0` (zero context lines): the
repo's `trailing-whitespace` pre-commit hook runs on all files with no
exclusions and strips the single leading space that marks a blank *context*
line in a unified diff, which corrupts an ordinary `-U3` patch quoted in a
committed file (confirmed: `git apply --check` on the previously quoted `-U3`
block failed with `error: corrupt patch`). `-U0` has no context lines, so
there is nothing for the hook to strip. Apply with `git apply --unidiff-zero`:

```diff
diff --git a/src/bulk/mod.rs b/src/bulk/mod.rs
index 40a90f4..bf5c0de 100644
--- a/src/bulk/mod.rs
+++ b/src/bulk/mod.rs
@@ -729,0 +730,2 @@ impl BulkSpec {
+        let mut t_par = std::time::Duration::ZERO;
+        let mut t_ser = std::time::Duration::ZERO;
@@ -737,0 +740 @@ impl BulkSpec {
+                let t_par0 = std::time::Instant::now();
@@ -784,0 +788 @@ impl BulkSpec {
+                t_par += t_par0.elapsed();
@@ -785,0 +790 @@ impl BulkSpec {
+                let t_ser0 = std::time::Instant::now();
@@ -790,0 +796 @@ impl BulkSpec {
+                t_ser += t_ser0.elapsed();
@@ -793,0 +800,7 @@ impl BulkSpec {
+        let p = t_par.as_secs_f64();
+        let s = t_ser.as_secs_f64();
+        eprintln!(
+            "[instr] parallel={p:.3}s serial={s:.3}s serial_frac={:.4}",
+            s / (p + s)
+        );
+
```

Run commands (default BCF format — the brief's Step 3 does not set
`VCFIXTURE_BENCH_FORMAT`; `VCFIXTURE_BENCH_REPS` is also unset, so `reps=1`,
**single-shot** measurements sitting on the ~10% run-to-run spread already
established elsewhere in this document, not averaged over repetitions):

```
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_WORKERS=1 ./target/release/examples/bulk_bench
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

Results:

| workers | parallel (t_par) | serial (t_ser) | `serial_frac` |
|---|---|---|---|
| 1 | 24.670s | 0.742s | 0.0292 |
| 8 |  7.135s | 0.229s | 0.0310 |

The `workers=1` value is the control: with only one worker there is nothing to
park at the barrier, so `serial_frac=0.0292` at `workers=1` is measuring the
drain loop's intrinsic cost (memcpy + bookkeeping), not contention. Going from
`workers=1` to `workers=8`, `serial_frac` barely moves — 0.0292 → 0.0310, a
0.18-point rise — which is consistent with no material barrier-contention
effect at 8 workers beyond the loop's fixed cost.

**Denominator reconciliation.** `serial_frac`'s denominator is
`t_par + t_ser` (instrumented time only, summed across `stream_contigs`'s
chunk loop), not the run's reported wall clock `min_s` — `stream_contigs` is
one step inside `BulkSpec::write`, which also does profile validation,
`compute_layouts`'s own separate `pool.install` gap-sum pass
(`src/bulk/mod.rs:647-655`), header construction, writer creation,
`finish_and_index`, and the summary JSON write, none of which the
instrumentation times. At `workers=1`: `24.670 + 0.742 = 25.412s` instrumented
vs `min_s = 25.485s` reported, a residual of `0.073s` (0.29% of `min_s`)
outside both timers. At `workers=8`: `7.135 + 0.229 = 7.364s` instrumented vs
`min_s = 7.448s` reported, a residual of `0.084s` (1.13% of `min_s`) outside
both timers. Worst-casing that entire residual as additional serial time —
i.e. assuming, pessimistically, that all of it is barrier-adjacent rather than
setup/teardown work outside `stream_contigs` entirely — gives
`(0.742 + 0.073) / 25.485 ≈ 0.032` at `workers=1` and
`(0.229 + 0.084) / 7.448 ≈ 0.042` at `workers=8`. Both worst-case figures are
still comfortably under the 5% exoneration threshold, so this reconciliation
strengthens rather than weakens the exoneration reading below.

**Amdahl comparison.** Solving `S = 1 / (s + (1-s)/p)` for `s` at `p = 4` gives
`s = (p - S) / (S(p - 1))`. Inverting for `s` is sensitive to the third
decimal of `S`, so the arithmetic below uses `S(4)` computed directly from the
underlying `min_s` ratio to 3 decimal places, not the Scaling curves table's
2-decimal display (`27.577 / 8.205 = 3.361`, `18.529 / 6.026 = 3.075`; the
table rounds these to `3.36` and `3.07` respectively — same numbers, tighter
precision here). The instrumented runs above use the default BCF format, so
the like-for-like comparison is the **BCF-unpinned** curve from Task 2, whose
`S(4) = 3.361` (min_s 8.205s vs 27.577s at `workers=1`, from the Scaling
curves table above):

```
s = (4 - 3.361) / (3.361 * (4 - 1)) = 0.639 / 10.083 ≈ 0.0634
```

i.e. Amdahl's law says a 6.34% serial fraction would fully explain
BCF-unpinned's `S(4) = 3.361`. For reference, since Gate C fired on the VCF
series, the same arithmetic against VCF-unpinned's `S(4) = 3.075` (min_s 6.026s
vs 18.529s):

```
s = (4 - 3.075) / (3.075 * (4 - 1)) = 0.925 / 9.225 ≈ 0.1003
```

i.e. VCF-unpinned's scaling would need a 10.03% serial fraction to be fully
explained by Amdahl's law alone — a larger implied fraction than BCF's, which
is expected since VCF has no bgzf compression pool to compete with rayon for
cores, so any residual sub-linearity has to come from somewhere other than that
oversubscription.

**Reading.** Measured `serial_frac` at `workers=8` (**0.031**, 3.1%) is under a
third of the BCF-unpinned Amdahl-implied fraction (0.063, 6.3%) — and under a
third of the VCF-unpinned one too (0.100, 10.0%). Both instrumented values sit
comfortably under the brief's 5% exoneration threshold, and the `workers=1`
control shows the loop's baseline cost (2.9%) accounts for nearly all of the
`workers=8` figure (3.1%), leaving at most ~0.2 points attributable to
contention at 8 workers. **The chunk barrier's serial drain is exonerated: it
is not the source of the sub-linear scaling gate C flagged.** The gap between
measured speedup and linear scaling lies elsewhere — consistent with
Hypothesis 2's two untested candidates (glibc allocator arena-lock contention,
rayon+bgzf thread oversubscription), which Task 4 investigates.

### Profiling

`samply` is installed but could not run: `perf_event_paranoid` is `2` on this box and `samply` requires it at `1` or lower; `sudo` is not available non-interactively (`sudo: a password is required`). Fell back to `perf record -g` per the brief's instructions, built with `CARGO_PROFILE_RELEASE_DEBUG=true` for symbolized frames (release binary otherwise ships no debug info; this env var was only used for the profiling run, not committed to `Cargo.toml`).

```
cd $CLAUDE_JOB_DIR/tmp
CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release --example bulk_bench --features bulk
VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 \
  perf record -g -F 999 -o perf.data -- <path-to-binary>
perf report -i perf.data --stdio -g none   # flat self-time table
```

Workload reduced to a single (2000, 20000) point (80M cells, ~7.8s) rather than the full sweep, to keep profiling time bounded.

Top self-time symbols (43,859 samples, `cycles:u`):

| self% | symbol |
|---|---|
| 18.43% | `_int_malloc` (libc) |
| 10.54% | `malloc` (libc) |
| 9.40% | `_int_free` (libc) |
| 4.38% | `malloc_consolidate` (libc) |
| 4.26% | `SampleStats::new` |
| 3.69% | `noodles_bcf::…encode_genotype_str` |
| 3.51% | `core::fmt::write` |
| 3.07% | `cfree` (libc) |
| 2.95% | `RandomState::hash_one::<&usize>` |
| 2.87% | `core::fmt::Formatter::pad_integral` |
| 2.49% | `gen_record::<ChaCha8Rng>` |
| 2.38% | `Vec<Option<Value>>::from_iter` |
| 2.35% | `noodles_bcf::…write_genotype_values` |
| 2.24% | `Hasher::write` (siphash) |
| 2.03% | `noodles_bcf::…write_samples` |
| 1.73% | `<i8 as Display>::fmt` |
| 1.59% | `unlink_chunk` (libc) |
| 1.48% | `String::write_str` |
| ... | (long tail, each <1.2%) |

glibc allocator-family symbols (`_int_malloc` + `malloc` + `_int_free` + `malloc_consolidate` + `cfree` + `unlink_chunk`) sum to **~47% of self time**. This is now the single largest cost bucket by a wide margin, and it is not attributable to any one Rust function — it is the aggregate cost of many small per-sample heap allocations.

**DSO- and compression-level sanity check** (re-derived from the same `perf.data`, `perf report --stdio -g none --sort=dso`, no re-run of the benchmark): self time splits ~49% `bulk_bench` (the binary, including statically-linked noodles/zlib-rs code), ~49% `libc-2.28.so` (almost entirely the allocator family above), ~1.5% `[unknown]`. There is no large unsymbolized DSO hiding time — the `[unknown]` bucket is small and, on inspection, is dominated by zero-self-time call-graph artifacts (garbage addresses from broken frame-pointer unwinding), not real cost. Compression (`zlib_rs::*`, the BCF writer's bgzf DEFLATE) sums to **~2% of self time across ~18 distinct call sites** (AVX2 match-finding, CRC32, portable fallbacks) — genuinely negligible, not a hidden cost, despite not appearing in the top-18 table above.

**Deferred question — `to_record_buf`'s per-sample `Vec<Option<Value>>`/`String` GT clone vs. noodles' BCF genotype encoder:**

Summing self time by side of the boundary, symbol-by-symbol at one decimal place so the arithmetic is auditable (the final totals are then rounded to whole percents — rounding each symbol first and summing the rounded values would itself introduce a ~1-point error, so precision is kept until the last step):

- vcfixture-side value construction: `SampleStats::new` 4.26 + `core::fmt::write` 3.51 + `<i8 as Display>::fmt` 1.73 + `pad_integral` 2.87 + `pad_integral::write_prefix` 0.90 + `SampleStats::value_for` 0.88 + `Vec<Option<Value>>::from_iter` 2.38 + `String::write_str` 1.48 + the `to_record_buf` closure fold 1.15 + `String::clone` 0.98 = 20.14 → **~20%** of self time.
- noodles BCF genotype encoder, excluding the shared hash lookup below: `encode_genotype_str` 3.69 + `encode_genotype_str::encode` 1.19 + `write_genotype_values` 2.35 + `write_samples` 2.03 + `Sample::get_index` 1.16 + the `Value` `From` conversion 0.75 = 11.17 → **~11%** of self time.
- Shared string-keyed lookup (`RandomState::hash_one::<&usize>` 2.95 + siphash `Hasher::write` 2.24 = 5.19 → ~5%): `src/bulk/` uses only `BTreeMap` internally, no `HashMap`, so this almost certainly belongs to noodles' own FORMAT-key index lookup (it appears immediately alongside `Sample::get_index` in the symbol list) rather than to anything in `vcfixture`. Counting it on the encoder side gives noodles' encoder 11.17 + 5.19 = 16.36 → **~16%** of self time.

So: **~20% (`to_record_buf`) vs. ~16% (noodles' encoder, hash lookup included) — `to_record_buf` is larger, but not decisively so.** An earlier draft of this analysis omitted the hash-lookup symbols from the encoder side entirely (reporting ~12%), which overstated the gap as ~1.7x; with the hash lookup correctly included, the ratio is ~1.25x. Read this as "roughly comparable, `to_record_buf` somewhat ahead," not as a clear-cut answer.

**Limitation of this attribution:** ~10 of the ~20 vcfixture-side points (`core::fmt::write`, `pad_integral`(+`write_prefix`), `Vec<Option<Value>>::from_iter`, `String::write_str`) are generic Rust standard-library symbols with no `vcfixture`-specific name — they were assigned to the vcfixture side by *plausibility* (the `write!(gt, "{a}")` call is at `src/bulk/generate.rs:306`, the `Vec<Option<Value>>` construction is at `src/bulk/generate.rs:372`), not by measurement. A call graph would settle this definitively, but `perf report --stdio --sort=self,symbol` (the flag combination needed to safely sort a call-graph view) reproducibly segfaults on this box's `perf 6.3.10` with Rust demangling, and the flat `-g none` list used throughout this profile discards caller information entirely. The vcfixture-vs-noodles split above should be read as a reasonable but unverified assignment, not a call-graph-confirmed measurement — this is on top of, not instead of, the general sampling-precision caveat above.

If this is pursued further, the lever is allocation count, not formatting or encoding logic: `SampleStats::new` allocates a `String` per sample per record for GT (with `write!` formatting each allele digit), and `to_record_buf` allocates a `Vec<Option<Value>>` per sample per record. For `Payload::GtOnly` that is 2 allocations × 40M samples-per-record-pairs in the (2000, 20000) profiled workload. **Not optimized in this task per the brief** — recorded as a finding only. Filed as [issue #26](https://github.com/d-laub/vcfixture-rs/issues/26), citing the shares above.

A separate, unrelated defect surfaced while reading `src/bulk/mod.rs` during this task: a mid-stream error from block encoding or the writer leaves a truncated, un-indexed file at the destination path with no cleanup (`stream_contigs`'s `?` skips `finish_and_index`). This predates the block-pipeline change and is not a regression from it. Filed as [issue #27](https://github.com/d-laub/vcfixture-rs/issues/27).

### Stop condition

The project's rule is "no optimization is kept unless the benchmark shows it wins and the oracle still passes": the benchmark shows a substantial, reproducible win (2.51x–3.08x wall-clock speedup, average ~2.81x, with wall clock falling substantially rather than staying flat as worker count rises), `cargo test --all-features` passes in full (153 tests, 0 failures, verified at this task's final commit), and the next hot spot found by profiling (glibc allocator overhead, ~47% of self time) is new, separate optimization work — not unfinished work from this change — so the change is kept and this measurement task stops here.
