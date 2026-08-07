# Bulk generation measurements (issue #22)

Harness: `cargo run --release --example bulk_bench --features bulk`
Machine: `carter-cn-03`, Slurm allocation of **4 physical cores** plus their SMT
siblings (see "Hardware and allocation" below — *not* 8 independent cores, even
though `nproc` reports 8), 4.18.0-553.36.1.el8_10.x86_64

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
Machine: same box and same allocation (`carter-cn-03`, 4 physical cores plus SMT
siblings; see "Hardware and allocation" below), 4.18.0-553.36.1.el8_10.x86_64

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

All (samples, records) points, workers, and workload sizes above are **single-shot measurements**, not averaged over repeated runs. The same (2000, 20000) workload at `workers=8` was independently re-measured three times over the course of this task (7.049s in the full sweep, 7.502s in the Step 4 scaling run, 7.762s in the Step 5 profiling run) — a 10% spread. That figure describes `workers=8` runs only and **must not be generalized**: the per-cell spreads measured later in this document (see "Measurement noise and anchor stability", below) run from 1% to 61% depending on worker count, and the low-worker-count cells are far noisier than the high ones. Every row in the sweep above is a `workers=8` row, so ±10% is the right band to read them at; the overall 2.5x–3.1x win is far larger than that, so no conclusion above changes, but individual row speedups should be read as accurate to roughly ±10%, not to the three significant figures the table prints.

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

This A/B is one run per side, not averaged, and both sides are `workers=1` runs — the noisiest regime in this study (see "Measurement noise and anchor stability", below) — so the ~1.31x/23.6% figures should be read as accurate in direction but not pinned to three significant figures. The anchor re-measurement in that section puts *clean-window* single-shot variability at `workers=1` at about ±3% (five spread-out BCF `workers=1` single-shots landed in 25.193s–26.829s), and the new-code side of this A/B, 26.103s, sits comfortably inside that cluster. Under a conservative ±10% band the worst-case reading (old run 10% low, new run 10% high: 30.761s vs. 28.713s) still gives ~1.07x — above 1x, but only just. The honest statement is therefore: in a clean measurement window the *direction* of this result is secure; what a single-shot A/B cannot defend against is one side landing in a contaminated window, which is exactly what the min-of-3 estimator used for the scaling curves guards against and this A/B does not. Also worth flagging: the two binaries were run at the same `seed=42`, but `block_rng` changed from `(seed, block_idx)` to `(seed, block_idx, stream)` in an earlier task in this same body of work, so the old and new binaries generate different variant *content* at that seed — same distributions, same record and cell counts (verified: both `n_records_total()` values and both `cells` denominators are 80,000,000), but not byte-identical records. At 20k records this is very unlikely to move timing, but it is an unstated difference between the two sides worth naming for completeness.

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

Speedup = 3.48x for 8x the workers — but 8 is the wrong denominator. Per the
Hardware and allocation section above, this box's Slurm allocation is **4
physical cores** (`0-3`) plus their SMT siblings (`48-51`), not 8 independent
cores; `workers` defaults from `nproc`/`available_parallelism()`, both of
which count logical CPUs and cannot see that distinction. Against the correct
denominator of 4, 3.48x is **87% parallel efficiency**, not 43.5%. The full
six-point sweep in Scaling curves (below) gives a cleaner reading of the same
quantity from one consistent series rather than this single two-point A/B:
BCF-unpinned reaches `S(4)=3.07` (**77% efficient**) by `workers=4` and keeps
climbing through the SMT range to `S(8)=3.52` (**88% efficient**) at
`workers=8`. VCF-unpinned (no bgzf pool) climbs monotonically to `S(4)=3.07`
(77%) and `S(8)=3.36` (84%); BCF-pinned (`taskset -c 0-3`) does not climb
monotonically — it dips slightly at `workers=6` before recovering — but
covers a similar high-70s-to-low-80s-percent efficiency range from
`workers=4` onward. (The 3.48x two-point reading above and the six-point
BCF-unpinned curve's `S(8)=3.52` agree closely, which they did not before the
`workers=1` anchor was re-measured — see "Measurement noise and anchor
stability", below.) Wall clock falls sharply from `workers=1` to `workers=4`
in all three curves (25.193s→8.205s BCF-unpinned, 23.996s→7.565s BCF-pinned,
18.529s→6.026s VCF-unpinned), confirming the serial `O(cells)` stage
Hypothesis 2 targeted is gone as a *dominant* cost. What remains is a
smaller, genuine sub-linear component, not the near-50%-efficiency shortfall
first reported.

One caveat applies to every `S(w)` figure in this document, including the
77%/88% above: the `workers=1` baseline is not itself single-threaded for the
BCF and VCF.gz paths. `workers` sizes both the rayon pool
(`rayon::ThreadPoolBuilder::new().num_threads(self.workers.get())`,
`src/bulk/mod.rs`) and the bgzf multithreaded writer's compression pool
(`bgzf::io::multithreaded_writer::Builder::default().set_worker_count(workers)`,
`src/bulk/writer.rs:149`) from the same value, so `workers=1` still runs a
second bgzf compression thread overlapping with generation — a two-thread
run is standing in as the "1" in `S(w) = min_s(1) / min_s(w)`. A genuinely
single-threaded BCF baseline would be **slower**, not faster: with no second
thread, DEFLATE would run serialized after each block's generation instead of
overlapping with it. Because the recorded baseline already has that overlap,
it is smaller than a true 1-thread time, and dividing by a too-small baseline
*understates* `S(w)`. So every BCF/VCF.gz `S(w)` in this document is a lower
bound on true parallel speedup rather than an overstatement — the conclusion
is unchanged, but note the mechanism runs through the baseline being fast, not
slow.

Uncompressed VCF is the exception: with no bgzf pool to size, its `workers=1`
is genuinely single-threaded, so its `S(w)` column needs no such correction.
Its `min_s(1)` of 18.529s is lower than BCF's 25.193s, but that is **not**
evidence for or against the overlap mechanism — overlap is a speedup
mechanism, and the gap is simply DEFLATE work that the VCF path never
performs. The two `workers=1` numbers are also not strictly like-for-like:
same `BulkSpec`, same seed, same 80M cells, but a different encoder path and a
materially different output file, so read the 18.529 vs 25.193 difference as
"compressed output costs more to produce", nothing finer.

The two candidates originally named here — allocator arena-lock contention
and rayon+bgzf thread oversubscription — were both tested by direct
measurement rather than left as profile-bucket speculation, and a third
candidate (a residual serial barrier in the per-chunk write loop) was found
and tested as well:

- **Allocator.** Task 4 ran two interventions rather than another profile
  (Allocator interventions, below). Forcing a single arena
  (`MALLOC_ARENA_MAX=1`) made `workers=8` 9.16x slower, so per-thread arenas
  are load-bearing — but that is exactly glibc's default behavior, already in
  effect under every measurement in this document. Swapping to mimalloc made
  every worker count 1.63x–1.73x faster in absolute wall clock, but
  mimalloc's own `S(4)`/`S(8)` (3.26/3.52) land level with glibc's
  (3.07/3.52) — and, decisively, mimalloc's *own* `S(4)` of 3.26 is itself
  only 81% efficient against 4 physical cores, so even an allocator built to
  avoid arena-lock contention entirely does not scale linearly on this
  workload. Verdict: allocator cost is large and real, but it is per-thread
  overhead, not the cause of the sub-linear curve.
- **Thread oversubscription.** The VCF ablation settles this directly: VCF
  has no bgzf pool at all, so at `workers=4` it is 4 rayon threads on 4
  physical cores with nothing else competing for them — and it still only
  reaches `S(4)=3.07`, 77% efficient (Scaling curves, below). Oversubscription
  cannot explain a shortfall that persists with no second thread pool to
  oversubscribe with — and after the `workers=1` anchor re-measurement the two
  unpinned curves' `S(4)` values are indistinguishable (VCF 3.075 with no bgzf
  pool, BCF 3.070 with one), so removing the pool entirely buys no scaling at
  all at the physical-core count. A genuine oversubscription signature *is* visible, but
  only in the BCF-pinned curve at worker counts above the physical core count
  (`S(6)=3.10` dipping below `S(4)=3.17` before recovering to `S(8)=3.29`) —
  real, but too small and too narrowly scoped to account for the residual
  seen at `workers=4` in any of the three curves.
- **Per-chunk serial write barrier.** Task 3 instrumented the boundary
  directly (Serial-fraction measurement, below) instead of reasoning from the
  Amdahl-implied figure alone: `stream_contigs` collects each chunk's encoded
  records with `par_iter().collect()`, which drains the rayon pool fully
  before a serial `write_encoded` loop runs. Measured `serial_frac` is
  0.029–0.031, well under the ~0.100 Amdahl's law implies is needed to explain
  the measured `S(4)` (0.1009 BCF-unpinned, 0.1003 VCF-unpinned). The barrier
  costs about 3%, not the ~10% the curve would require.

**Verdict: held, but the residual is real and, after this task, genuinely
unexplained.** All three candidates examined — allocator contention, bgzf
thread oversubscription, and the chunk-write serial barrier — were tested by
direct measurement or intervention, not left as speculation, and none of them
accounts for the gap between the measured curves and ideal linear scaling
against 4 physical cores. That gap is smaller than the original
43.5%-efficiency framing suggested (77%–88% across the BCF-unpinned curve from
`w=4` to `w=8`, and 77% at `w=4` on the VCF-unpinned curve), but it has
not been explained, and this document should say so plainly rather than
retire the question behind a plausible-sounding cause that was never tested.
Closing it is future work, not a conclusion this investigation reached:
candidates that remain **untested** are memory-bandwidth saturation,
per-thread allocation overhead that is inherent to the workload rather than
contended (distinct from the arena-contention question Task 4 already
answered), the separate `compute_layouts` parallel pass, and SMT siblings on
this box contributing measurably less than a full physical core each.

### Allocator interventions

Task 3 instrumented `stream_contigs`'s chunk barrier directly and exonerated
it: measured `serial_frac` (2.9%-3.1%) sits well under the Amdahl-implied
serial fraction (10.1% BCF, 10.0% VCF) needed to explain the measured `S(4)`.
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
— a 9.16x degradation against a widest-measured noise band of 61% (see
"Measurement noise and anchor stability", below), so the effect is more than
an order of magnitude past anything noise could produce and is unambiguously
real. (The default arm here, 8.467s, reads 18%
above Task 2's originally recorded BCF-unpinned `workers=8` `min_s` of 7.153s
from a separate session. That 18% is a *cross-session* shift in a min-of-3
estimator, not a within-cell spread: within Task 2's own `workers=8` cell the
three reps spanned only 9% (7.153–7.802). The two figures are measuring
different things and both are real — the honest band for comparing a `min_s`
from one session against a `min_s` from another is therefore ~20%, wider than
any single cell's internal spread. It changes nothing about the Step 1
reading: the arena effect is ~9x, dwarfing either baseline by close to an
order of magnitude.)

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

| workers | mimalloc min_s | glibc min_s (BCF-unpinned) | speedup |
|---|---|---|---|
| 1 | 15.416 | 25.193 | 1.63x |
| 4 | 4.736 | 8.205 | 1.73x |
| 8 | 4.375 | 7.153 | 1.63x |

The glibc column is the BCF-unpinned series used throughout this document, not
a fresh remeasurement in this session. Its `workers=1` entry is the
re-measured anchor (25.193s, min over 8 samples — see "Measurement noise and
anchor stability", below), which replaced the originally recorded 27.577s;
`workers=4` and `workers=8` are Task 2's originals. Per the
baseline-discrepancy note under Step 1 above, a same-session default-arm rerun
of `workers=8` read 18% higher (8.467s vs 7.153s), which is the honest size of
the session-to-session noise band this comparison is subject to; Task 2's
series is used here because the brief mandates it as the like-for-like
baseline, not because `w=4`/`w=8` were reverified this session. Note the
asymmetry this creates and read the comparison accordingly: glibc's
`workers=1` figure rests on 8 samples while mimalloc's rests on 3, and a
min-over-more-samples is biased low, so glibc's `w=1` number is the better
sampled — and the more pessimistic — of the two.

mimalloc's own scaling, `S(w) = mimalloc_min_s(1) / mimalloc_min_s(w)`,
against the same BCF-unpinned series used throughout this document:

| workers | mimalloc `S(w)` | glibc `S(w)` | mimalloc vs glibc | mimalloc efficiency vs 4 cores |
|---|---|---|---|---|
| 4 | 3.26 (`15.416/4.736`) | 3.07 (`25.193/8.205`) | 6.0% higher | 81% |
| 8 | 3.52 (`15.416/4.375`) | 3.52 (`25.193/7.153`) | 0.05% higher | 88% |

Two separate findings here. mimalloc makes every worker count **uniformly
faster in absolute terms** — 1.63x to 1.73x faster wall clock at `w=1`, `4`,
and `8` alike, with no trend toward a *larger* speedup at higher worker
counts. That is consistent with the profile's ~47% allocator self-time
finding: allocation cost is a large, real, per-call overhead, and a faster
allocator removes a roughly constant fraction of it independent of thread
count — exactly the "swap that makes everything uniformly faster without
changing the shape of the curve" case, not one that changes scaling.

The scaling comparison itself is a dead heat. mimalloc's `S(8)` equals
glibc's to within 0.05% (3.5237 vs 3.5220 — a difference far below any noise
band this document measures), and its `S(4)` is 6.0% higher, which is inside
the noise band *and* is exactly the direction the sampling asymmetry noted
above would produce on its own (glibc's `w=1` anchor is a min over 8 samples,
mimalloc's over 3). So the raw numbers support neither "mimalloc raised
`S(w)`" nor "mimalloc lowered `S(w)`" — they support "no measurable
difference."

The decisive observation is not the comparison at all, but mimalloc's own
efficiency column: **3.26 at `w=4` is 81% of linear, and 3.52 at `w=8` is 88%
of 4 cores.** An allocator with thread-local heaps and no shared arena lock —
a design that cannot suffer arena-lock contention — still falls well short of
linear scaling on this workload. Whatever caps `S(w)` here is therefore not
arena-lock contention, and that conclusion does not depend on the
glibc-vs-mimalloc delta being resolvable above the noise. One caveat on the
comparison is addressed next for completeness.

**Caveat — the Amdahl-mechanical artifact, and why it no longer bites.** If
some cost is roughly fixed in absolute seconds regardless of allocator (Task
3's chunk-barrier serial drain, `finish_and_index`, the summary JSON write,
process startup), then shrinking the parallel portion — exactly what a faster
allocator does — makes that fixed cost a *larger share* of a smaller total,
which mechanically depresses the faster allocator's `S(w)` with no change in
contention whatsoever. Under the originally recorded `workers=1` anchor
(27.577s) this mattered a great deal, because mimalloc's `S(8)` appeared 8.6%
*below* glibc's and the artifact was the leading candidate explanation for
that gap. With the re-measured anchor (25.193s) there is no such gap left to
explain: the two `S(8)` values are 3.5237 and 3.5220. Solving for the fixed
cost `F` that would equalize the two allocators' fixed-cost-corrected scaling,

```
(mimalloc_min_s(1) - F) / (mimalloc_min_s(8) - F)
  = (glibc_min_s(1) - F) / (glibc_min_s(8) - F)

(15.416 - F) / (4.375 - F) = (25.193 - F) / (7.153 - F)
=> F ≈ -0.007s
```

i.e. the equalizing `F` is indistinguishable from zero — the raw numbers are
already equal without any correction. Applying Task 3's actually measured
fixed cost anyway moves things marginally in mimalloc's *favour* rather than
against it: at `F = 0.229s`, glibc's corrected `S(8)` rises `3.52 → 3.61` and
mimalloc's rises `3.52 → 3.66`; at `F = 0.313s` (the worst case including the
untimed residual — see the Serial-fraction measurement section's
reconciliation), `3.52 → 3.64` and `3.52 → 3.72`. So the artifact is real in
principle, but it cannot be masking relieved contention here, because there is
no residual deficit in mimalloc's `S(8)` for it to mask; if anything the
correction says mimalloc scales fractionally better, by an amount far inside
the noise. The supportable reading is the neutral one: **on this evidence
mimalloc neither raised nor lowered `S(w)` measurably**, and the load-bearing
finding is instead that mimalloc's own curve is itself sub-linear (81% at
`w=4`, 88% at `w=8`).

**Verdict.** Combining both interventions: allocator cost is large and real
— Step 2 confirms this causally, not just as a profile bucket, since a
different allocator saves 1.6x-1.7x of wall time at every worker count — but
it is **not the primary cause of the sub-linear scaling curves** in the
Scaling curves section below. Step 1 shows glibc's per-thread arenas are
load-bearing scaffolding that, if removed, produces catastrophic (9.16x)
contention; but that scaffolding is exactly what the default configuration
under which every scaling curve in this document was measured already has,
and Step 2 shows a qualitatively different, more contention-resistant
allocator does not measurably unlock additional scaling headroom on top of
it — and, more directly, that this contention-immune allocator's *own* curve
is sub-linear to essentially the same degree (81% at `w=4`, 88% at `w=8`),
which no amount of relieved-but-masked contention can explain away. At this
thread count, allocator behavior reads as **mostly pure per-thread overhead,
not primarily a scaling limiter** — the concern belongs to issue #26
(allocation count/size), not to this investigation's search for what caps
`S(w)` below linear. The gap between measured `S(4)`/`S(8)` and the ideal
4x/8x — Amdahl-implied serial fractions of ~10% (10.1% BCF, 10.0% VCF)
against a measured chunk-barrier serial fraction of only ~3% (Task 3) —
**remains open** after this task, with
both originally-proposed candidates (chunk-barrier serial drain, allocator
contention) now tested and exonerated — the allocator one on two independent
grounds: swapping allocators does not move `S(w)` measurably, and the
contention-immune allocator's own curve is sub-linear by the same margin.

### Scaling curves

Fixed workload throughout (2000 samples × 20000 records, 80M cells, `seed=42`), swept over `VCFIXTURE_BENCH_WORKERS ∈ {1,2,3,4,6,8}`, 3 reps per cell (`VCFIXTURE_BENCH_REPS=3`), reporting `min_s`/`med_s`/`max_s` per cell. `S(w) = min_s(1) / min_s(w)`; efficiency is `S(w) / min(w, 4)` — against **4 physical cores**, per the hardware section above (so `w=6` and `w=8` both divide by 4, not by `w`). The BCF-unpinned `workers=1` anchor was subsequently re-measured and revised downward; see "Measurement noise and anchor stability" immediately below for the samples, the decision rule, and which curves it changed.

**BCF, unpinned** (`Cpus_allowed_list: 0-3,48-51`, all 8 logical CPUs available):

| workers | min_s | med_s | max_s | S(w) | S(w)/min(w,4) |
|---|---|---|---|---|---|
| 1 † | 25.193 | 26.829 | 36.312 | 1.00 | 1.00 |
| 2 | 13.309 | 14.042 | 18.509 | 1.89 | 0.95 |
| 3 | 10.614 | 12.355 | 12.526 | 2.37 | 0.79 |
| 4 | 8.205 | 8.305 | 8.374 | 3.07 | 0.77 |
| 6 | 7.580 | 7.650 | 8.024 | 3.32 | 0.83 |
| 8 | 7.153 | 7.363 | 7.802 | 3.52 | 0.88 |

† This cell pools **8 samples**, not 3: the original 3-rep cell
(27.577/31.642/36.312) plus the five spread-out single-shots from the anchor
re-measurement below, which found the original cell had landed in a
contaminated window. `min_s` is the min over all eight; `med_s` is the same
upper-middle order statistic the harness reports, applied to the pooled eight.
Every other cell in every curve still rests on 3 reps — see "Measurement noise
and anchor stability" below for the decision rule and for why this asymmetry
is disclosed rather than removed.

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

Two of the three curves rise monotonically from `w=1` through `w=8`; the BCF-pinned curve does not. BCF-unpinned reaches 3.07x by `w=4` (77% efficient against 4 physical cores) and keeps climbing more slowly through the SMT range to 3.52x at `w=8` (88% efficient) — consistent with the SMT siblings contributing real, if diminished, throughput on this workload rather than pure oversubscription noise. VCF-unpinned also climbs monotonically throughout, from 3.07x at `w=4` to 3.36x at `w=8`. BCF-pinned tracks BCF-unpinned reasonably closely through `w=4` (their `min_s` values differ by 8% at `w=4`, inside the cross-session `min_s` band characterized below) but shows the oversubscription signature the other two curves don't: `S(6)=3.10` dips slightly below `S(4)=3.17` before recovering to `S(8)=3.29` — running 6 or 8 rayon+bgzf threads on 4 pinned physical cores produces mild, non-monotonic thrashing rather than a clean plateau, expected once the process is confined to fewer cores than it requests threads for.

VCF's *relative* speedup — `S(w)` is normalized against each curve's own `workers=1` baseline, and the three baselines differ (25.193s, 23.996s, 18.529s), so an `S` comparison across curves compares two differently-normalized quantities and cannot be read as an absolute throughput comparison — is at or below BCF-**unpinned**'s from `w=4` onward: a statistical tie at `w=4` (3.075 vs 3.070, a 0.1% difference), then 3.14 vs 3.32 at `w=6` and 3.36 vs 3.52 at `w=8`. In absolute wall clock VCF is of course faster everywhere (6.026s vs 8.205s at `w=4`), because it does no DEFLATE work; that is a different statement from scaling better. VCF's peak (3.36x at `w=8`) is not materially higher than BCF-unpinned's peak (3.52x) — this is the comparison Gate B's "VCF reaches a materially higher peak speedup than BCF" clause evaluates, using the brief's default BCF-unpinned series, so Gate B does not fire on this evidence. VCF is *not*, however, uniformly below BCF-**pinned**: it trails at `w=4` (3.07 vs 3.17) but overtakes it at `w=6` (3.14 vs 3.10) and `w=8` (3.36 vs 3.29) — VCF's peak is in fact ~2% above BCF-pinned's peak (3.36 vs 3.29). A 2% crossover is far inside every noise band this document measures, so it reads as directionally suggestive rather than conclusive proof that uncompressed rayon-only writes outrun bgzf-pinned writes past `w=4` — but it is exactly the comparison Gate B's peak-speedup clause is asking about, so it is worth flagging even though it does not change the verdict (the gate's default series is BCF-unpinned, against which VCF's peak is lower, not higher). VCF's efficiency plateaus at 77%–84% from `w=4` onward (77% at `w=4`, 78% at `w=6`, 84% at `w=8`) rather than climbing toward 100%, indicating a real, if modest, sub-linear component in rayon-only fan-out scaling that is not explained by the bgzf pool at all — and the BCF-unpinned curve now plateaus over a similar 77%–88% range *with* the pool present, coinciding exactly at `w=4`, which is the cleanest single statement of how little the pool matters to scaling here.

### Measurement noise and anchor stability

This document previously carried three mutually inconsistent noise models — a
"~10%" figure quoted in several places, an "~18%" figure quoted in one, and
the scaling tables themselves, which imply neither. This section replaces all
of them with one model derived from the recorded data.

**Per-cell spread.** Every scaling-curve cell above is 3 reps taken back to
back, so `max_s/min_s - 1` is a direct read of within-cell, within-session
run-to-run spread. (The BCF-unpinned `workers=1` cell is shown here as its
*original* 3-rep cell, before the pooling described below, so that all 18
entries are like-for-like.)

| workers | BCF-unpinned | BCF-pinned | VCF-unpinned |
|---|---|---|---|
| 1 | 32% ‡ | 8% | 3% |
| 2 | 39% | 1% | 61% |
| 3 | 18% | 21% | 9% |
| 4 | 2% | 7% | 3% |
| 6 | 6% | 2% | 3% |
| 8 | 9% | 2% | 8% |

‡ the original 3-rep cell, before the anchor re-measurement below.

The structure is unmistakable and it is not "~10%" anywhere: **every spread of
18% or more occurs at `w ≤ 3`, and every cell at `w ≥ 4` falls between 2% and
9%.** Median spread is 18% for `w ≤ 3` and 3% for `w ≥ 4`. Two plausible
mechanisms compound, both pointing the same way. First, exposure: a `w=1` run
occupies the allocation for ~25s against ~7s at `w=8`, a 3.5x wider window in
which another tenant of this 96-CPU node, or the node's own background load,
can collide with it. Second, and more important, redundancy: when a run has
eight rayon threads and one is descheduled, the other seven keep working and
work-stealing rebalances around it, so an interruption costs a fraction of the
run; when a run has one thread, that thread's wall clock *is* the run's wall
clock and any stall lands on the total at full weight. Long, undiversified
runs are noisy; short, many-threaded ones are not.

**Per-run spread is not estimator error.** The two must not be conflated. The
tables report `min_s`, and a minimum over 3 draws is far more stable than any
single draw — that is precisely why `min_s` was chosen. The 61% spread in
VCF's `w=2` cell, for instance, comes from one 17.061s outlier against a
10.615s minimum; the minimum is what `S(w)` uses and the outlier does not
touch it. What the min-of-3 estimator does *not* protect against is all three
reps landing inside the same contaminated window, which is what happened at
BCF-unpinned `w=1` and is the subject of the re-measurement below.

**The ~18% figure, reconciled.** The Allocator interventions section observed
a same-configuration `workers=8` `min_s` of 8.467s against Task 2's 7.153s, an
18% difference. That is a *cross-session* shift in the estimator, a different
quantity from within-cell spread (9% for that same cell). Both are real. The
resulting rule of thumb, which the rest of this document now uses: **within a
session at `w ≥ 4`, single runs vary by 2%–9% and `min_s` is stable to a few
percent; across sessions, a `min_s`-to-`min_s` comparison carries about ±20%;
at `w ≤ 3` single runs can vary by tens of percent and no single-shot
comparison at those worker counts should be trusted to better than that.**

**Anchor re-measurement.** Because `S(w) = min_s(1) / min_s(w)`, every speedup
in this document divides by a `workers=1` cell — and BCF-unpinned's was the
noisiest `workers=1` cell in the study (32% spread, and it implied `S(2)=2.07`, i.e. 104%
efficiency: super-linear speedup, which is the standard red flag for an
inflated baseline). Two independent `workers=1` BCF measurements elsewhere in
this document (26.103s in the Hypothesis 1/2 sections, 25.485s in the
*instrumented* Serial-fraction run — i.e. carrying extra overhead and still
faster) both sat well below it. So the anchor was stress-tested directly.

Method: the release binary was rebuilt and verified free of instrumentation
(`strings ./target/release/examples/bulk_bench | grep -c '\[instr\]'` returned
0), then **five independent passes** were
run, each measuring the three `workers=1` configurations once each
(`VCFIXTURE_BENCH_REPS=1`) in the order BCF-unpinned, BCF-pinned
(`taskset -c 0-3`), VCF-unpinned. Spreading five single samples across ~7
minutes, rather than taking `REPS=5` on one cell, is deliberate: a single
contaminated window can swallow a whole `REPS=n` cell but cannot swallow five
passes separated in time. Individual times, in seconds:

| pass | BCF-unpinned | BCF-pinned | VCF-unpinned |
|---|---|---|---|
| 1 | 26.061 | 26.162 | 20.796 |
| 2 | 26.250 | 26.079 | 19.292 |
| 3 | 25.193 | 25.860 | 20.561 |
| 4 | 25.620 | 25.258 | 19.220 |
| 5 | 26.829 | 25.755 | 19.637 |

Decision rule, fixed in advance: for each curve take the minimum over *all*
samples now available (the five new ones plus the original cell's three) —
the same min estimator the study already uses, given more data. If that new
minimum is more than 5% below the recorded `min_s(1)`, adopt it and re-derive
the curve; otherwise keep the recorded anchor and record the stability check.

| curve | recorded `min_s(1)` | min over 5 new | min over all 8 | below recorded | branch |
|---|---|---|---|---|---|
| BCF-unpinned | 27.577 | 25.193 | **25.193** | 8.6% | **adopt** |
| BCF-pinned | 23.996 | 25.258 | **23.996** | 0.0% | keep |
| VCF-unpinned | 18.529 | 19.220 | **18.529** | 0.0% | keep |

So only BCF-unpinned changed: its anchor moved 27.577s → **25.193s**, and its
`S(w)` and efficiency columns above are computed against 25.193s. BCF-pinned
and VCF-unpinned are unchanged — for both, five fresh spread-out samples
failed to beat the recorded minimum at all, which is the stability check
passing.

Three things this makes visible, all worth recording:

- **The original `w=1` cell was a contaminated window, not a slow
  configuration.** The five new BCF-unpinned samples span 25.193–26.829, a 6.5%
  spread; the original three span 27.577–36.312. Six of the eight samples now
  lie within 9.5% of each other and the two outliers are both from the original
  triple. That is the signature of three consecutive reps hitting interference,
  exactly the failure mode `REPS=n` cannot detect and time-separated passes can.
- **The implausible readings are gone.** `S(2)` was 2.07 — super-linear, 104%
  efficient at two workers; it is now 1.89 (95%). And
  BCF-pinned's `w=1` was 13% *faster* than BCF-unpinned's, which would have
  meant restricting a 1-worker run to 4 of 8 logical CPUs sped it up; the two
  are now 4.8% apart, with the fresh samples showing them essentially level
  (unpinned 25.193–26.829, pinned 25.258–26.162).
- **BCF-pinned's anchor is now the under-sampled one.** It rests on 3 samples
  where BCF-unpinned's rests on 8. A min over more samples is biased low, so
  it is BCF-unpinned's own anchor that is pulled down the most — meaning
  BCF-unpinned's `S(w)` column is, if anything, the more *deflated* of the
  three, not the least: its published speedups are the most conservative in
  the document, the same direction noted above for glibc vs. mimalloc.
  Separately, and worth recording on its own: 23.996s also sits below all
  five of BCF-pinned's fresh samples, which speaks to 23.996s being a
  reasonable minimum for that curve, not to which curve's `S(w)` is biased
  relative to the other. Between the two effects, the net relative direction
  between BCF-pinned's and BCF-unpinned's `S(w)` columns is not determinable
  from this data. This asymmetry — anchors sampled 8-and-3 deep while every
  other cell has 3 — is disclosed rather than removed; closing it would mean
  re-running the whole sweep, which would answer no open question in this
  document.

None of this changes any verdict. It does move the BCF-unpinned headline
efficiencies down (84% → 77% at `w=4`, 96% → 88% at `w=8`), which makes the
residual sub-linearity *larger* than previously published, not smaller — the
correction is in the conservative direction.

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
**single-shot** measurements, not averaged over repetitions — read them at the
per-worker-count bands characterized in "Measurement noise and anchor
stability" above, which is tens of percent for the `workers=1` run and 2%–9%
for the `workers=8` one):

```
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_WORKERS=1 ./target/release/examples/bulk_bench
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

Results:

| workers | parallel (t_par) | serial (t_ser) | `serial_frac` |
|---|---|---|---|
| 1 | 24.670s | 0.742s | 0.0292 |
| 8 |  7.135s | 0.229s | 0.0310 |

The `workers=1` value is the control: with only one worker there is no rayon
pool to park at the barrier, so `serial_frac=0.0292` at `workers=1` bounds the
drain loop's cost in the absence of any barrier contention (the next paragraph
shows part of even that is bgzf back-pressure rather than the loop's own
memcpy). Going from
`workers=1` to `workers=8`, `serial_frac` barely moves — 0.0292 → 0.0310, a
0.18-point rise — which is consistent with no material barrier-contention
effect at 8 workers beyond the loop's fixed cost.

One feature of the table deserves naming because it cuts *for* the
exoneration. In absolute seconds `t_ser` falls 0.742s → 0.229s from
`workers=1` to `workers=8` — 3.2x faster. A genuinely serial memcpy cannot get
3.2x faster by adding workers, so a large part of what the serial timer
captures is evidently not fixed serial work at all but bgzf back-pressure: at
`workers=1` the single compression thread cannot drain the staging buffer as
fast as `write_encoded` fills it, and the drain loop blocks; at `workers=8` the
compression pool keeps up and the same loop is nearly free. The genuinely fixed
component is therefore *at most* 0.229s, not 0.742s — which makes the barrier
an even smaller candidate for the residual than the `serial_frac` figures
alone suggest.

This also puts Task 4's Amdahl-mechanical caveat, which models a fixed cost
`F` constant across worker counts, in tension with these very numbers. That
tension does not change its conclusion: recomputing that caveat's `S(8)`
comparison with a *variable* fixed cost — 0.742s at `w=1`, 0.229s at `w=8` —
gives glibc 3.531 and mimalloc 3.539, a gap of −0.008 against a raw gap of
−0.002. Both are indistinguishable from zero, exactly as the constant-`F`
treatment found, so the caveat's reading survives its own assumption being
imperfect.

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
2-decimal display (`25.193 / 8.205 = 3.070`, `18.529 / 6.026 = 3.075`; the
table rounds both to `3.07` — same numbers, tighter precision here). The
instrumented runs above use the default BCF format, so the like-for-like
comparison is the **BCF-unpinned** curve, whose `S(4) = 3.070` (min_s 8.205s
against the re-measured 25.193s anchor at `workers=1`, from the Scaling curves
table above):

```
s = (4 - 3.070) / (3.070 * (4 - 1)) = 0.930 / 9.211 ≈ 0.1009
```

i.e. Amdahl's law says a 10.09% serial fraction would fully explain
BCF-unpinned's `S(4) = 3.070`. For reference, since Gate C fired on the VCF
series, the same arithmetic against VCF-unpinned's `S(4) = 3.075` (min_s 6.026s
vs 18.529s):

```
s = (4 - 3.075) / (3.075 * (4 - 1)) = 0.925 / 9.225 ≈ 0.1003
```

i.e. VCF-unpinned's scaling would need a 10.03% serial fraction to be fully
explained by Amdahl's law alone. The two implied fractions are effectively
identical — 10.09% with a bgzf compression pool competing for cores, 10.03%
with no such pool at all. That coincidence is itself informative: whatever
serial-or-serializing component Amdahl's law is inferring here is the same
size with and without the pool, which independently corroborates the
oversubscription exoneration recorded under Hypothesis 2.

**Reading.** Amdahl's `s` is defined as the serial fraction of the
*sequential* run, so the like-for-like measured quantity is the `workers=1`
control: **0.0292 (2.9%) against BCF-unpinned's implied 0.1009 (10.1%)** —
under a third, and likewise under a third of VCF-unpinned's implied 0.1003.
The `workers=8` figure (0.031, 3.1%) tells the same story against the same
implied fractions, and its near-identity with the `workers=1` control is what
shows the barrier is not *becoming* more expensive under contention: the loop's
baseline cost accounts for nearly all of it, leaving at most ~0.2 points
attributable to contention at 8 workers. Both instrumented values sit
comfortably under the brief's 5% exoneration threshold. **The chunk barrier's
serial drain is exonerated: it is not the source of the sub-linear scaling
gate C flagged.** Where the gap does lie is not settled by this document. Both
of Hypothesis 2's originally named candidates were subsequently tested and
neither survived — the allocator by Task 4's two interventions (and by the
observation that mimalloc's own curve is sub-linear), thread oversubscription
by Task 2's uncompressed-VCF ablation, which has no bgzf pool at all and still
falls short. See the Hypothesis 2 verdict above for the resulting list of
candidates that remain genuinely untested.

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

## Allocation reduction (#26) and atomic destination writes (#27)

Measured on branch `fix-26-27-alloc-partial-output`, which stacks on this
document's branch. Two independent changes, measured separately so each can be
judged on its own:

- **Part A** — a per-thread scratch `RecordBuf` reused across records, taking
  steady-state allocations from four per sample per record to one.
- **Part B** — mimalloc as the binaries' global allocator.

Part C (atomic destination writes, #27) is a correctness fix with no
measurable performance component and is not benchmarked here.

### Protocol

Three binaries, all built `--release`:

```bash
# baseline: branch point 501ce4d, no source change, glibc
cargo build --release --features bulk --example bulk_bench

# Part A only: source change, glibc
cargo build --release --no-default-features --features bulk --example bulk_bench

# Parts A+B: source change, mimalloc
cargo build --release --features bulk --example bulk_bench
```

`bulk_bench` now stamps `alloc={mimalloc,system}` into its header line, so a
sweep cannot be misattributed to the wrong binary after the fact.

Worker counts were **round-robined across passes** — all three binaries at
`w=1`, then all three at `w=2`, and so on, five times — rather than running
each cell's reps back to back. Running reps consecutively lets one noise burst
contaminate all of them identically, which defeats a min-of-N estimator; that
is the failure that inflated this document's earlier BCF-unpinned anchor by
8.6%. Round-robining also means the five `workers=1` readings are separated by
a full sweep of other cells, so the anchor-stability check is built into the
protocol rather than bolted on afterwards.

Workload: 2000 samples x 20000 records = 80M cells, BCF, `Payload::GtOnly`,
seed 42, unpinned, one rep per invocation. Machine: `carter-cn-03`,
`Cpus_allowed_list: 0-3,48-51` — **4 physical cores** plus SMT siblings.
Efficiency below is `S(w) / min(w, 4)` against those 4 physical cores, not the
8 logical CPUs `nproc` reports.

### Raw readings

Seconds, in pass order. Spread is `(max - min) / min`.

| binary | w | pass 1 | pass 2 | pass 3 | pass 4 | pass 5 | min | spread |
|---|---|---|---|---|---|---|---|---|
| baseline | 1 | 22.152 | 21.478 | 21.221 | 21.378 | 21.225 | 21.221 | 4.4% |
| baseline | 2 | 11.370 | 11.275 | 11.311 | 11.228 | 11.335 | 11.228 | 1.3% |
| baseline | 3 | 9.272 | 9.226 | 8.374 | 8.338 | 8.174 | 8.174 | 13.4% |
| baseline | 4 | 6.617 | 7.072 | 6.654 | 6.605 | 6.785 | 6.605 | 7.1% |
| baseline | 6 | 6.902 | 6.836 | 6.514 | 6.525 | 6.678 | 6.514 | 6.0% |
| baseline | 8 | 6.451 | 6.607 | 6.752 | 6.491 | 6.433 | 6.433 | 5.0% |
| Part A | 1 | 12.490 | 12.470 | 12.594 | 12.530 | 12.514 | 12.470 | 1.0% |
| Part A | 2 | 6.622 | 6.699 | 6.640 | 6.665 | 6.618 | 6.618 | 1.2% |
| Part A | 3 | 5.190 | 5.078 | 4.784 | 6.130 | 5.681 | 4.784 | 28.1% |
| Part A | 4 | 3.829 | 3.813 | 3.823 | 3.910 | 3.846 | 3.813 | 2.5% |
| Part A | 6 | 3.866 | 3.828 | 3.542 | 3.770 | 3.726 | 3.542 | 9.1% |
| Part A | 8 | 3.528 | 3.434 | 3.495 | 3.474 | 3.434 | 3.434 | 2.7% |
| A+B | 1 | 8.063 | 8.045 | 8.125 | 8.067 | 8.113 | 8.045 | 1.0% |
| A+B | 2 | 4.383 | 4.451 | 4.407 | 4.370 | 4.443 | 4.370 | 1.9% |
| A+B | 3 | 3.569 | 3.383 | 3.198 | 3.595 | 3.353 | 3.198 | 12.4% |
| A+B | 4 | 2.610 | 2.633 | 2.718 | 2.567 | 2.835 | 2.567 | 10.4% |
| A+B | 6 | 2.393 | 2.428 | 2.353 | 2.333 | 2.396 | 2.333 | 4.1% |
| A+B | 8 | 2.178 | 2.174 | 2.258 | 2.189 | 2.134 | 2.134 | 5.8% |

**Anchor stability.** All three `workers=1` anchors are stable, so none needed
the re-measurement this document's earlier curve did. Part A and A+B span 1.0%
across five time-separated passes. The baseline spans 4.4%, but that is driven
entirely by a high first pass (22.152); the remaining four cluster within 1.2%
and three of them agree with the adopted minimum to within 0.02%.

**Noise.** `w=3` is the noisiest cell for all three binaries (13.4%, 28.1%,
12.4%) while every other cell is 1–10%. The same worker count being worst
across three independently built binaries suggests something systematic rather
than a stray burst — plausibly an interaction between 3 workers, the
`2 * workers` chunk size, and 4 physical cores — but this was not investigated
and is recorded as an observation only. It does not affect the conclusions,
which rest on `w=1` and `w>=4`.

### Speedup

Absolute wall clock at matched worker counts, from the minima above:

| w | baseline | Part A | vs base | A+B | vs base | B's marginal |
|---|---|---|---|---|---|---|
| 1 | 21.221 | 12.470 | **1.70x** | 8.045 | **2.64x** | 1.55x |
| 2 | 11.228 | 6.618 | **1.70x** | 4.370 | **2.57x** | 1.51x |
| 3 | 8.174 | 4.784 | **1.71x** | 3.198 | **2.56x** | 1.50x |
| 4 | 6.605 | 3.813 | **1.73x** | 2.567 | **2.57x** | 1.49x |
| 6 | 6.514 | 3.542 | **1.84x** | 2.333 | **2.79x** | 1.52x |
| 8 | 6.433 | 3.434 | **1.87x** | 2.134 | **3.01x** | 1.61x |

**Part A pays for itself independently of the allocator swap**, which is the
question this document committed to answering separately: 1.70x–1.87x on its
own, with no dependence on Part B. The two compose to 2.56x–3.01x.

**Part B's marginal contribution is 1.49x–1.61x**, below the 1.63x–1.73x
measured for mimalloc against *unmodified* code in the "Allocator
interventions" section above. That direction is expected rather than
contradictory: Part A removes three of four per-sample allocations, so there is
less allocator work left for a faster allocator to accelerate. The two
measurements are consistent, not in tension.

### Scaling

| w | baseline S(w) | eff | Part A S(w) | eff | A+B S(w) | eff |
|---|---|---|---|---|---|---|
| 1 | 1.00 | 100% | 1.00 | 100% | 1.00 | 100% |
| 2 | 1.89 | 94.5% | 1.88 | 94.2% | 1.84 | 92.0% |
| 3 | 2.60 | 86.5% | 2.61 | 86.9% | 2.52 | 83.9% |
| 4 | 3.21 | 80.3% | 3.27 | 81.8% | 3.13 | 78.4% |
| 6 | 3.26 | 81.4% | 3.52 | 88.0% | 3.45 | 86.2% |
| 8 | 3.30 | 82.5% | 3.63 | 90.8% | 3.77 | 94.2% |

Scaling *shape* is essentially unchanged through `w=4` (80.3% / 81.8% / 78.4%),
which is the same sub-linear behaviour this document has been unable to
explain. The apparent improvement at `w>=6` is real but must be read carefully:
`S(w)` divides by each binary's *own* `min_s(1)`, and Part A speeds up the
1-worker case by 1.70x while speeding up the 8-worker case by 1.87x, so the
ratio moves without any claim that contention was relieved.

The `workers=1` caveat from the "Scaling curves" section still applies to all
three curves: for BCF, `workers` sizes the bgzf compression pool as well as the
rayon pool, so a 1-worker run is not single-threaded and dividing by it
**understates** true scaling.

**The residual sub-linear deficit remains unexplained.** This work changes the
absolute cost of bulk generation; it does not resolve the open scaling
question. The candidates listed earlier — memory-bandwidth saturation, inherent
per-thread allocation overhead, the separate `compute_layouts` parallel pass,
and SMT siblings contributing less than a full core — remain **untested**.

### Peak RSS

MB, maximum over passes:

| w | baseline | Part A | A+B |
|---|---|---|---|
| 1 | 13.1 | 13.7 | 28.9 |
| 2 | 23.3 | 22.1 | 50.8 |
| 4 | 41.7 | 43.8 | 68.6 |
| 8 | 77.4 | 77.9 | 108.3 |

**Part A is RSS-neutral** (77.9 vs 77.4 MB at `w=8`), which was not a foregone
conclusion: holding scratch buffers alive for a worker thread's lifetime rather
than freeing them per record could have inflated the high-water mark, and does
not.

**mimalloc costs roughly 40% more peak RSS** (108.3 vs 77.4 MB at `w=8`),
consistent with per-thread heaps that are not returned to the OS as eagerly.
That is the price of Part B's 1.5x, and it is why the `mimalloc` feature has an
opt-out.

### Comparability caveat

The baseline binary here is built from commit `501ce4d` — the same code the
"Scaling curves" section measured — yet reads 21.221s at `w=1` against that
section's 25.193s, and 6.433s at `w=8` against 7.153s. Both are ~16% and ~10%
lower.

Absolute numbers from different measurement sessions on this cluster are
therefore **not comparable**, and none of the ratios in this section are
derived by mixing the two. Every figure above comes from binaries measured
against each other within a single round-robin sweep. It is also possible the
earlier re-anchored 25.193s was still noise-inflated; that cannot be settled
retrospectively and does not affect any conclusion here, since the earlier
section's claims are likewise internally consistent.
