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

### Hypothesis 1 — "removing the redundant generation pass roughly halves generation CPU"

**Not directly measurable from wall-clock alone**, and the honest answer is: the harness has no separate timer for "generation" vs "encode," so a clean before/after CPU split for generation specifically was not obtained. What is measurable:

- The profile (see below) shows a single `gen_record` symbol at 2.49% self time, with no duplicate/second generation call site visible in the hot path — consistent with the gap-only span pass having removed the second full-record generation pass, since a leftover redundant pass would show up as a second `gen_record`-shaped consumer of CPU.
- The overall 2.5–3.1x wall-clock speedup is larger than a simple 2x (halved generation CPU alone would predict) or a "removed one of two redundant passes" story on its own, because Hypothesis 2 (parallel encode) is also contributing on top of it. The two effects are not separable with this harness.

Verdict: **plausible but not confirmed** by this measurement. The absence of a duplicate generation call in the profile is consistent with the pass having been removed, but no isolated "generation-only CPU" number exists to confirm "roughly halves." Do not read the 2.5–3.1x combined speedup as proof of the halving claim specifically.

### Hypothesis 2 — "moving encode/summary into the fan-out removes the serial O(cells) stage, so wall clock scales with `--threads`"

Measured directly at fixed workload (2000 samples × 20000 records, 80M cells):

```
VCFIXTURE_BENCH_WORKERS=1: 26.103s  (3.263e-7 s/cell)
VCFIXTURE_BENCH_WORKERS=8:  7.502s  (9.377e-8 s/cell)
```

Speedup = 3.48x for 8x the workers (43.5% parallel efficiency). Wall clock is **not flat** — it falls substantially, confirming the serial `O(cells)` stage is no longer the dominant cost. But scaling is clearly sub-linear, not the full 8x a purely-parallel fan-out with a cheap serial tail would predict.

Verdict: **held, partially.** The serial `O(cells)` encode stage that dominated before is gone — s/cell nearly triples in efficiency from 1→8 workers, it does not stay flat. But claiming full linear scaling would be dishonest: at ~44% parallel efficiency, something is still capping speedup below 8x. The profile below points at glibc allocator contention (`malloc`/`free`/`_int_malloc`/`_int_free`/`malloc_consolidate` together are ~47% of self time) as the likely cause — many small per-sample allocations (one `String` for GT, one `Vec<Option<Value>>` per sample) contending on the glibc allocator's arena lock under 8 concurrent workers — rather than a leftover serial pipeline stage. This is consistent with an Amdahl-style residual (implied serial fraction ≈ 18.6% from the 3.48x/8-worker numbers), but the source of that residual looks like allocator contention, not a serial `O(cells)` stage in the pipeline itself.

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

glibc allocator-family symbols (`_int_malloc` + `malloc` + `_int_free` + `malloc_consolidate` + `cfree` + `unlink_chunk`) sum to **~47.4% of self time**. This is now the single largest cost bucket by a wide margin, and it is not attributable to any one Rust function — it is the aggregate cost of many small per-sample heap allocations.

**Deferred question — `to_record_buf`'s per-sample `Vec<Option<Value>>`/`String` GT clone vs. noodles' BCF genotype encoder:**

Summing self time by side of the boundary (excluding the shared allocator cost, which both sides drive):

- vcfixture-side value construction (`SampleStats::new`, `core::fmt::write`, `<i8 as Display>::fmt`, `pad_integral`(+`write_prefix`), `SampleStats::value_for`, `Vec<Option<Value>>::from_iter`, `String::write_str`, the `to_record_buf` closure fold, `String::clone`): **~20.1%** of self time.
- noodles BCF genotype encoder (`encode_genotype_str`(+`::encode`), `write_genotype_values`, `write_samples`, `Sample::get_index`, the `Value` `From` conversion): **~11.8%** of self time.

`to_record_buf`'s own per-sample value construction is the larger of the two (~20% vs. ~12%), so the deferred question resolves in favor of **`to_record_buf`, not noodles' encoder**, as the dominant non-allocator cost. But the more important finding is that neither of these is the actual dominant cost: **glibc allocator overhead (~47%) exceeds both combined (~32%)**. The `RandomState`/`Hasher` hashing (~5.2% combined) sits on the encoder side of the boundary (`Sample::get_index` performing a keyed lookup) and was counted there.

If this is pursued further, the lever is allocation count, not formatting or encoding logic: `SampleStats::new` allocates a `String` per sample per record for GT (with `write!` formatting each allele digit), and `to_record_buf` allocates a `Vec<Option<Value>>` per sample per record. For `Payload::GtOnly` that is 2 allocations × 40M samples-per-record-pairs in the (2000, 20000) profiled workload. **Not optimized in this task per the brief** — recorded as a finding only. Filing a follow-up issue is left to the user since this measurement task should not itself decide to open one without being asked; the ~20%+~47% shares above are the numbers a follow-up issue would cite.

### Stop condition

The Phase-0 target was hit — a 2.5–3.1x wall-clock speedup with wall clock now falling substantially (not flat) as worker count rises clears the "benchmark shows it wins" bar this task exists to check, and the next hot spot (glibc allocator overhead, ~47% of self time) is a new, separate optimization opportunity rather than unfinished work from this change.
