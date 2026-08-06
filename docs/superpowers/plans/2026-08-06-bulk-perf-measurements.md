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

### Hypothesis 1 — "removing the redundant generation pass roughly halves generation CPU"

**Measured directly with an isolated A/B comparison**, using the `VCFIXTURE_BENCH_WORKERS` knob added in this task to factor parallelism out entirely. A throwaway `git worktree` was checked out at the pre-change baseline commit `9e64a0c` (which predates the `bulk_bench` harness), the current `examples/bulk_bench.rs` was copied into it unmodified (it uses only public `BulkSpec`/`Payload`/`Profile`/`Size` API that is identical at `9e64a0c`), built `--release`, and run once at `WORKERS=1` for the same workload used in the Hypothesis 2 scaling check:

```
OLD code (9e64a0c), WORKERS=1: 34.179s  (4.272e-7 s/cell)
NEW code (HEAD),    WORKERS=1: 26.103s  (3.263e-7 s/cell)
```

At 1 worker there is no parallelism benefit to factor out — this isolates whatever the pipeline-shape change (gap-only span pass + moving encode/summary into a single interleaved pass) is worth on its own, end to end. Result: **1.31x speedup, a 23.6% reduction in end-to-end single-threaded wall clock** (26.103s / 34.179s = 0.7637; 34.179 − 26.103 = 8.076s saved).

This A/B is one run per side, not averaged, on top of the ~10% run-to-run spread already established above for this exact workload — so the ~1.31x/23.6% figures should be read as accurate in direction but not pinned to three significant figures. Even at a worst-case reading (old run 5% low, new run 5% high: 32.470s vs. 27.408s) the ratio is still ~1.18x, safely above 1x — the *direction* of the result is not noise, only its precise magnitude is uncertain by a few points. Also worth flagging: the two binaries were run at the same `seed=42`, but `block_rng` changed from `(seed, block_idx)` to `(seed, block_idx, stream)` in an earlier task in this same body of work, so the old and new binaries generate different variant *content* at that seed — same distributions, same record and cell counts (verified: both `n_records_total()` values and both `cells` denominators are 80,000,000), but not byte-identical records. At 20k records this is very unlikely to move timing, but it is an unstated difference between the two sides worth naming for completeness.

**Read against what H1 actually claims — generation CPU, not end-to-end wall clock:** H1 says removing the redundant pass roughly halves *generation CPU*, not that it halves total single-threaded time. At `9e64a0c`, the spans loop calls `generate_contig` to lay out the file, and the write loop calls it again to produce records — two identical full-genotype generation passes. Deleting one of two identical passes halves generation CPU **by construction**; that mechanism is not in question. The 8.076s this A/B measures is consistent with exactly that: if generation CPU was halved, and generation was roughly a quarter of the old single-threaded total (8.076s / 34.179s ≈ 24%, close to the observed 23.6% end-to-end reduction, since removing *one whole* generation pass out of two removes very close to *all* the redundant-pass cost in one step), then an end-to-end reduction of ~24% — not ~50% — is exactly what a reader should expect. Encoding, header/CSI work, and I/O make up the other ~76% of single-threaded cost and are untouched by this specific mechanism (Hypothesis 2 addresses the effect of restructuring *those* stages, not H1).

Verdict: **held, as a claim about generation CPU specifically — not as a claim about end-to-end wall clock.** The mechanism H1 describes (deleting one of two identical generation passes) did occur and does halve generation CPU by construction. What this A/B rules out is reading H1 as implying a 2x *end-to-end* single-threaded speedup: generation is only a modest fraction (roughly a quarter, by this measurement) of total single-threaded cost, so halving it nets a real but smaller ~1.31x end-to-end improvement. Do not read the ~1.31x figure as H1 failing, and do not read it as licensing a claim of a 2x end-to-end win either — both would misstate what was measured.

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
