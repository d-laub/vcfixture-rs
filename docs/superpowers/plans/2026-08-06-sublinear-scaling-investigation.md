# Sub-Linear Scaling Investigation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace PR #28's "sub-linear scaling is not fully explained" caveat with a measured account, and correct the two documents that publish the wrong 43.5%-parallel-efficiency figure.

**Architecture:** Extend the existing `bulk_bench` example with repetition and output-format knobs, collect three scaling curves on the reference workload, then run conditional deeper probes only where the curves leave a gap. Measurement work only — no `src/` change is committed.

**Tech Stack:** Rust 2021, `cargo run --release --example bulk_bench --features bulk`, `taskset`, glibc `MALLOC_ARENA_MAX`, `gh pr edit`.

## Global Constraints

- **Benchmark tasks MUST NOT run in parallel with each other or with any other
  CPU-consuming work.** The cgroup holds 4 physical cores; two concurrent
  benchmarks contend for them and invalidate both results. This plan is
  strictly sequential. Do not use `dispatching-parallel-agents` for Tasks 2–4.
- **Hardware ground truth:** `Cpus_allowed_list: 0-3,48-51` on `carter-cn-03`
  = 4 physical cores (0–3) plus their SMT siblings (48–51). 2 threads/core,
  Intel Xeon E5-4650 v3. `nproc` reports 8 (logical), which is what
  `std::thread::available_parallelism()` returns.
- **Reference workload, unchanged across every measurement:**
  `VCFIXTURE_BENCH_SAMPLES=2000`, `VCFIXTURE_BENCH_RECORDS=20000` (80M cells),
  `contigs = ["chr1", "chr2"]`, `Payload::GtOnly`, `seed = 42`.
- **Never detach a process** — no `nohup`, `setsid`, `disown`, or trailing `&`.
  An NFS-blocked process on this cluster is unkillable and drained a node on
  2026-07-29. Verify no bench process survives before reporting a task done.
- **`TMPDIR=$CLAUDE_JOB_DIR/tmp`** (`/carter/users/dlaub/.claude/jobs/968fbc40/tmp`)
  for every bench run. `bulk_bench` writes output via `env::temp_dir()`, and
  parallel background jobs share `/tmp`.
- **Never hand-edit `CHANGELOG.md` or the version in `Cargo.toml`** — `cz bump`
  generates both in CI.
- **Do not bypass prek hooks** with `--no-verify`. Hooks are already installed
  at `/carter/users/dlaub/projects/vcfixture-rs/.git/hooks/pre-commit`.
- **Raw benchmark output goes to `$CLAUDE_JOB_DIR/tmp`**, never under
  `~/.claude` directly and never onto NFS bulk paths.

---

### Task 1: Bench harness — repetition and format knobs

**Files:**
- Modify: `examples/bulk_bench.rs`

**Interfaces:**
- Consumes: `vcfixture::bulk::{BulkSpec, Format, Payload, Profile, Size}`.
  `Format` is `pub use writer::Format` from `src/bulk/mod.rs:44`, with variants
  `Bcf`, `VcfGz`, `Vcf`. `BulkSpec::format(f: Format) -> BulkSpec` at
  `src/bulk/mod.rs:428`.
- Produces: two new env knobs consumed by Tasks 2–4 —
  `VCFIXTURE_BENCH_REPS` (positive integer, default `1`) and
  `VCFIXTURE_BENCH_FORMAT` (one of `bcf` | `vcf` | `vcf.gz`, default `bcf`).
  A new output column layout: `samples records cells reps min_s med_s max_s
  s/cell peakRSS_MB`, where `s/cell` is derived from `min_s`.

- [ ] **Step 1: Add the two knob readers**

In `examples/bulk_bench.rs`, change the import line to add `Format`:

```rust
use vcfixture::bulk::{BulkSpec, Format, Payload, Profile, Size};
```

Then add these two functions immediately after the existing `workers()`
function:

```rust
/// Repetitions per sweep cell. Benchmarks on this box show ~10% run-to-run
/// spread, so a single shot cannot distinguish a real 5% effect from noise.
/// Reported as min/median/max; `s/cell` uses the min, the standard robust
/// estimator for "how fast can this machine do it" (noise only ever adds
/// time).
fn reps() -> usize {
    match env::var("VCFIXTURE_BENCH_REPS") {
        Ok(v) => {
            let n: usize = v.trim().parse().expect("reps must be a positive integer");
            assert!(n > 0, "reps must be a positive integer");
            n
        }
        Err(_) => 1,
    }
}

/// Output format, with the file extension to use for the bench output path.
/// `Vcf` is the ablation that matters: it is the only variant with no bgzf
/// compression pool, so it measures rayon scaling with nothing else competing
/// for cores.
fn format() -> (Format, &'static str) {
    match env::var("VCFIXTURE_BENCH_FORMAT").as_deref() {
        Ok("bcf") | Err(_) => (Format::Bcf, "bcf"),
        Ok("vcf") => (Format::Vcf, "vcf"),
        Ok("vcf.gz") => (Format::VcfGz, "vcf.gz"),
        Ok(other) => panic!("VCFIXTURE_BENCH_FORMAT must be bcf, vcf, or vcf.gz; got {other:?}"),
    }
}
```

- [ ] **Step 2: Rewrite `main`'s body to use them**

Replace the whole of `fn main()` with:

```rust
fn main() {
    let dir = env::temp_dir().join("vcfixture_bulk_bench");
    std::fs::create_dir_all(&dir).expect("create bench output dir");

    let samples = sweep("VCFIXTURE_BENCH_SAMPLES", &[500, 2_000, 8_000]);
    let records = sweep("VCFIXTURE_BENCH_RECORDS", &[5_000, 20_000]);

    let workers = workers();
    let reps = reps();
    let (format, ext) = format();
    println!("workers={workers} reps={reps} format={ext}");
    println!(
        "{:>8} {:>9} {:>12} {:>5} {:>9} {:>9} {:>9} {:>12} {:>10}",
        "samples", "records", "cells", "reps", "min_s", "med_s", "max_s", "s/cell", "peakRSS_MB"
    );

    for &n_samples in &samples {
        for &n_records in &records {
            let path = dir.join(format!("bench_{n_samples}_{n_records}.{ext}"));

            let mut times: Vec<f64> = Vec::with_capacity(reps);
            let mut n_records_total = 0u64;

            for _ in 0..reps {
                let profile = Profile::builtin("germline-1kgp").expect("built-in profile loads");

                let t0 = Instant::now();
                let summary = BulkSpec::new(profile)
                    .samples(n_samples as usize)
                    .contigs(["chr1", "chr2"])
                    .size(Size::Records(n_records))
                    .payload(Payload::GtOnly)
                    .seed(42)
                    .format(format)
                    .workers(workers)
                    .write(&path)
                    .expect("bulk generation succeeds");
                times.push(t0.elapsed().as_secs_f64());
                n_records_total = summary.n_records_total();

                // Clean between reps: every rep must write a fresh file, and
                // the index/summary siblings must not accumulate. The library
                // appends `.csi` to the full path (only for `Bcf`; see
                // `BulkWriter::finish_and_index`), not to the stem.
                let _ = std::fs::remove_file(&path);
                let _ = std::fs::remove_file(format!("{}.csi", path.display()));
                let _ = std::fs::remove_file(format!("{}.summary.json", path.display()));
            }

            times.sort_by(|a, b| a.partial_cmp(b).expect("elapsed times are finite"));
            let min = times[0];
            let med = times[times.len() / 2];
            let max = times[times.len() - 1];

            // Ploidy 2 is the germline-1kgp profile's dialed value; the cell
            // count is what cost is linear in.
            let cells = n_records_total * n_samples * 2;
            println!(
                "{n_samples:>8} {n_records:>9} {cells:>12} {reps:>5} {min:>9.3} {med:>9.3} \
                 {max:>9.3} {:>12.3e} {:>10.1}",
                min / cells as f64,
                peak_rss_kib() as f64 / 1024.0,
            );
        }
    }
}
```

- [ ] **Step 3: Update the module doc comment**

Replace the doc-comment block at the top of the file (lines 1–15, ending with
the `VCFIXTURE_BENCH_WORKERS` example) with:

```rust
//! Sweep bulk generation over cohort width and record count, reporting
//! seconds per cell so numbers line up directly with the ladder in issue
//! #22.
//!
//! Run with:
//!   cargo run --release --example bulk_bench --features bulk
//!
//! Override the sweep with `VCFIXTURE_BENCH_SAMPLES` /
//! `VCFIXTURE_BENCH_RECORDS` (comma-separated), e.g.
//!   VCFIXTURE_BENCH_SAMPLES=4000 VCFIXTURE_BENCH_RECORDS=250000 \
//!     cargo run --release --example bulk_bench --features bulk
//!
//! Override the worker count (passed to both the rayon pool and the bgzf
//! writer) with `VCFIXTURE_BENCH_WORKERS`, e.g. to measure scaling:
//!   VCFIXTURE_BENCH_WORKERS=1 cargo run --release --example bulk_bench --features bulk
//!
//! `VCFIXTURE_BENCH_REPS=3` repeats each sweep cell and reports min/median/max
//! seconds; `s/cell` is computed from the min. `VCFIXTURE_BENCH_FORMAT=vcf`
//! writes uncompressed VCF, which removes the bgzf compression pool entirely —
//! the ablation that separates rayon scaling from writer-pool contention.
```

- [ ] **Step 4: Build and verify the binary compiles clean**

Run:

```bash
cargo build --release --example bulk_bench --features bulk
```

Expected: compiles with no warnings.

- [ ] **Step 5: Verify the knobs work on a tiny workload**

Run each of these three commands (small sweep, seconds not minutes):

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=50 VCFIXTURE_BENCH_RECORDS=200 VCFIXTURE_BENCH_REPS=3 ./target/release/examples/bulk_bench
```

Expected: header line reads `workers=8 reps=3 format=bcf`; one data row with
`reps` column showing `3` and `min_s <= med_s <= max_s`.

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=50 VCFIXTURE_BENCH_RECORDS=200 VCFIXTURE_BENCH_FORMAT=vcf ./target/release/examples/bulk_bench
```

Expected: header reads `format=vcf`; the run succeeds (this proves the `Vcf`
path works end to end, since `finish_and_index` skips indexing for it).

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_FORMAT=bogus ./target/release/examples/bulk_bench
```

Expected: panics with `VCFIXTURE_BENCH_FORMAT must be bcf, vcf, or vcf.gz; got "bogus"`.

- [ ] **Step 6: Verify no bench output was left behind**

Run:

```bash
ls -la $CLAUDE_JOB_DIR/tmp/vcfixture_bulk_bench
```

Expected: the directory is empty (every file removed by the per-rep cleanup).
If files remain, the cleanup paths are wrong — fix before committing.

- [ ] **Step 7: Confirm the committed tree is still green**

Run:

```bash
cargo fmt --check
cargo clippy --all-features --all-targets -- -D warnings
```

Expected: both clean. (`cargo test` is not required here — this task touches
only an example binary, and the full suite runs in Task 5.)

- [ ] **Step 8: Commit**

```bash
git add examples/bulk_bench.rs
git commit -m "perf(bulk): add repetition and format knobs to the bench harness (#22)

VCFIXTURE_BENCH_REPS repeats each sweep cell and reports min/median/max,
so the ~10% run-to-run spread is quantified per point instead of assumed.
VCFIXTURE_BENCH_FORMAT selects the output format; the uncompressed vcf
setting removes the bgzf compression pool, isolating rayon scaling from
writer-pool contention.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Collect three scaling curves

**Files:**
- Create: `$CLAUDE_JOB_DIR/tmp/curve-bcf-unpinned.txt` (scratch, not committed)
- Create: `$CLAUDE_JOB_DIR/tmp/curve-bcf-pinned.txt` (scratch, not committed)
- Create: `$CLAUDE_JOB_DIR/tmp/curve-vcf-unpinned.txt` (scratch, not committed)
- Modify: `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md`

**Interfaces:**
- Consumes: `VCFIXTURE_BENCH_REPS` and `VCFIXTURE_BENCH_FORMAT` from Task 1.
- Produces: three speedup series `S(w) = t_min(1) / t_min(w)` for
  `w ∈ {1,2,3,4,6,8}`, one per variant, recorded as tables in the measurements
  doc. Tasks 3 and 4 gate on the BCF-unpinned and VCF-unpinned series.

**Note on scope:** the design doc made the VCF ablation conditional on Stage 1
showing turnover. It is unconditional here, and that is deliberate: at
`workers=4` the BCF path already runs 4 rayon + 4 bgzf threads on 4 physical
cores, so *every* point on the BCF curve is potentially oversubscribed, not
just the high ones. Without the VCF curve there is no clean read on rayon-only
scaling at any worker count. The ablation costs ~4 minutes.

- [ ] **Step 1: Record the hardware ground truth**

Run each command and capture the output verbatim for the doc section:

```bash
nproc
```
```bash
nproc --all
```
```bash
grep Cpus_allowed_list /proc/self/status
```
```bash
cat /sys/devices/system/cpu/cpu0/topology/thread_siblings_list
```
```bash
lscpu | rg -i 'thread|core|socket|model name'
```

Expected: `8`; `96`; `0-3,48-51`; `0,48`; 2 threads/core, 12 cores/socket, 4
sockets, Xeon E5-4650 v3. If any value differs from this, **stop and report** —
the allocation changed and every number in this plan needs re-derivation.

- [ ] **Step 2: Write the sweep script**

Create `$CLAUDE_JOB_DIR/tmp/sweep.sh`:

```bash
#!/bin/bash
# usage: sweep.sh <outfile> <format> [taskset-cpulist]
set -euo pipefail
out="$1"; fmt="$2"; pin="${3:-}"
: > "$out"
for w in 1 2 3 4 6 8; do
  if [ -n "$pin" ]; then
    runner="taskset -c $pin"
  else
    runner=""
  fi
  TMPDIR="$CLAUDE_JOB_DIR/tmp" \
  VCFIXTURE_BENCH_SAMPLES=2000 \
  VCFIXTURE_BENCH_RECORDS=20000 \
  VCFIXTURE_BENCH_REPS=3 \
  VCFIXTURE_BENCH_FORMAT="$fmt" \
  VCFIXTURE_BENCH_WORKERS="$w" \
    $runner ./target/release/examples/bulk_bench >> "$out" 2>&1
done
cat "$out"
```

Then `chmod +x $CLAUDE_JOB_DIR/tmp/sweep.sh`.

- [ ] **Step 3: Run the BCF unpinned sweep**

Run (foreground, ~4 minutes — do NOT background it):

```bash
$CLAUDE_JOB_DIR/tmp/sweep.sh $CLAUDE_JOB_DIR/tmp/curve-bcf-unpinned.txt bcf
```

Expected: six `workers=N reps=3 format=bcf` blocks, each with one data row.
`min_s` at `workers=1` should be near 26s, matching the recorded baseline.

- [ ] **Step 4: Run the BCF pinned-to-physical-cores sweep**

Run (foreground, ~5 minutes):

```bash
$CLAUDE_JOB_DIR/tmp/sweep.sh $CLAUDE_JOB_DIR/tmp/curve-bcf-pinned.txt bcf 0-3
```

This restricts the process to the 4 physical cores, excluding the SMT
siblings. Comparing its `workers=4` against the unpinned `workers=8` prices
SMT's actual contribution on this workload instead of assuming a textbook
figure.

- [ ] **Step 5: Run the VCF (no bgzf pool) unpinned sweep**

Run (foreground, ~4 minutes):

```bash
$CLAUDE_JOB_DIR/tmp/sweep.sh $CLAUDE_JOB_DIR/tmp/curve-vcf-unpinned.txt vcf
```

- [ ] **Step 6: Verify no bench process survived**

Run:

```bash
pgrep -a bulk_bench
```

Expected: no output (exit code 1). If anything is listed, kill it and
investigate before continuing — a surviving process poisons later
measurements and risks a node drain.

- [ ] **Step 7: Compute the speedup series**

For each of the three files, extract `min_s` per worker count and compute
`S(w) = min_s(1) / min_s(w)`. Also compute parallel efficiency against **4
physical cores**, i.e. `S(w) / min(w, 4)`.

Record as three tables with these exact columns:
`workers | min_s | med_s | max_s | S(w) | S(w)/min(w,4)`.

- [ ] **Step 8: Write the measurements-doc sections**

In `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md`, add a new
`## Hardware and allocation` section **before** the existing hypothesis
sections, containing the Step 1 output and this statement:

> Every measurement in this document ran inside a Slurm allocation of
> `Cpus_allowed_list: 0-3,48-51` on `carter-cn-03` — **4 physical cores** plus
> their SMT siblings, not 8 cores. `nproc` and
> `std::thread::available_parallelism()` both report 8 because they count
> logical CPUs, which is why `workers` defaulted to 8. Parallel efficiency must
> be computed against 4, not 8.

Then add a `### Scaling curves` subsection under Hypothesis 2 with the three
tables from Step 7 and a one-paragraph reading of the curve shapes.

Do **not** yet rewrite the Hypothesis 2 verdict — that happens in Task 5, once
any conditional probes have run.

- [ ] **Step 9: Commit**

```bash
git add docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md
git commit -m "docs(bulk): record the true core allocation and three scaling curves (#22)

The recorded 43.5% parallel efficiency used a denominator of 8 cores. The
allocation is Cpus_allowed_list 0-3,48-51 — 4 physical cores plus SMT
siblings — so the same 3.48x is ~87% efficiency against physical cores.

Adds curves at workers 1,2,3,4,6,8 (3 reps each) for BCF unpinned, BCF
pinned to physical cores, and uncompressed VCF, which has no bgzf pool.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

- [ ] **Step 10: Evaluate the gates and report them**

Using the BCF-unpinned series unless stated otherwise, evaluate all three
gates and state each verdict explicitly in the task report:

- **Gate A — question resolved.** `S(4) >= 3.2` **and** `S(8) >= S(4)`.
  Core count and SMT explain the curve; the residual serial fraction is small.
- **Gate B — writer-pool oversubscription.** `S(6) < S(4)` or `S(8) < S(4)` by
  more than 5%, **or** the VCF series reaches a materially higher peak speedup
  than the BCF series.
- **Gate C — real serial stage or contention.** `S(4) < 3.2` in the **VCF**
  series (the variant with no bgzf pool, so rayon scaling is unconfounded).

**Task 3 runs only if Gate C fires. Task 4 runs only if Gate C fires and Task
3 fails to account for the gap.** If only Gate A fires, skip directly to Task
5. If Gate B fires, note it for Task 5's write-up; it needs no extra probe
because the VCF curve already measures it.

The controller records skipped tasks in the ledger together with the numbers
that justified skipping.

---

### Task 3 (conditional — Gate C only): Measure the serial fraction directly

**Files:**
- Modify (temporarily, then revert): `src/bulk/mod.rs:735-791`
- Create: `$CLAUDE_JOB_DIR/tmp/instr.patch` (scratch, not committed)
- Modify: `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md`

**Interfaces:**
- Consumes: nothing from Task 1 beyond the built binary.
- Produces: a measured serial fraction (a single `f64` printed to stderr per
  run), to be compared against the Amdahl-implied fraction recomputed for 4
  cores.

**Why this task exists:** `stream_contigs` is a barrier per chunk.
`chunk.par_iter().map_init(...).collect()` fully drains before the serial
`for item in encoded` loop starts, so the O(bytes) `write_encoded` memcpy into
the bgzf staging buffer runs with the entire rayon pool parked. That is an
O(cells) serial stage — the shape Amdahl's law points at, and a candidate
neither originally-recorded hypothesis considered.

- [ ] **Step 1: Apply the instrumentation patch**

In `src/bulk/mod.rs`, inside `stream_contigs`, immediately after the line
`let chunk_blocks = 2 * self.workers.get();` add:

```rust
        let mut t_par = std::time::Duration::ZERO;
        let mut t_ser = std::time::Duration::ZERO;
```

Wrap the parallel region: change

```rust
                let encoded: Vec<Result<(Vec<u8>, BlockSummary), BulkError>> = pool.install(|| {
```

so it is preceded by `let t_par0 = std::time::Instant::now();`, and add
`t_par += t_par0.elapsed();` on the line immediately after the closing
`});` of that `pool.install` call.

Wrap the serial drain: change

```rust
                for item in encoded {
                    let (bytes, bs) = item?;
                    writer.write_encoded(&bytes)?;
                    summary.merge_block(id, &bs);
                }
```

to

```rust
                let t_ser0 = std::time::Instant::now();
                for item in encoded {
                    let (bytes, bs) = item?;
                    writer.write_encoded(&bytes)?;
                    summary.merge_block(id, &bs);
                }
                t_ser += t_ser0.elapsed();
```

Finally, immediately before the closing `Ok(summary)` of `stream_contigs`, add:

```rust
        let p = t_par.as_secs_f64();
        let s = t_ser.as_secs_f64();
        eprintln!(
            "[instr] parallel={p:.3}s serial={s:.3}s serial_frac={:.4}",
            s / (p + s)
        );
```

- [ ] **Step 2: Save the patch to scratch for reproducibility**

```bash
git diff src/bulk/mod.rs > $CLAUDE_JOB_DIR/tmp/instr.patch
```

- [ ] **Step 3: Build and run at the two informative worker counts**

```bash
cargo build --release --example bulk_bench --features bulk
```

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_WORKERS=1 ./target/release/examples/bulk_bench
```

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

Expected: each prints one `[instr] parallel=… serial=… serial_frac=…` line to
stderr. Record both `serial_frac` values.

Reading: if `serial_frac` at `workers=8` is close to the Amdahl-implied serial
fraction recomputed for 4 cores, the chunk barrier's drain is the explanation.
If `serial_frac` is small (say under 5%) the barrier is exonerated and the gap
lies elsewhere — proceed to Task 4.

- [ ] **Step 4: Revert the patch and prove the tree is clean**

```bash
git checkout -- src/bulk/mod.rs
```

```bash
git status --short
```

Expected: no output. The instrumentation must not reach any commit.

- [ ] **Step 5: Record the finding**

Add a `### Serial-fraction measurement` subsection to
`docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md` under
Hypothesis 2, stating the two `serial_frac` values, the reading above, and
quoting the contents of `$CLAUDE_JOB_DIR/tmp/instr.patch` in a fenced diff
block so the number is reproducible without the patch being committed.

- [ ] **Step 6: Commit**

```bash
git add docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md
git commit -m "docs(bulk): measure the chunk-barrier serial fraction directly (#22)

stream_contigs barriers per chunk: par_iter().collect() drains fully
before the serial write_encoded loop, so an O(bytes) memcpy runs with the
rayon pool parked. Measured rather than inferred from Amdahl's law. The
instrumentation patch is quoted in the doc, not committed.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4 (conditional — only if Task 3 leaves the gap unexplained): Allocator interventions

**Files:**
- Modify (temporarily, then revert): `Cargo.toml`, `examples/bulk_bench.rs`
- Modify: `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md`

**Interfaces:**
- Consumes: the built `bulk_bench` binary and the Task 2 curves.
- Produces: two causal readings on whether glibc's allocator limits scaling.

**Why these tests and not the originally-proposed one:** comparing glibc
self-time share between a 1-worker and an 8-worker profile is correlational —
a higher share at 8 workers is equally consistent with "more allocation work"
and with "contention". Both probes below are interventions.

- [ ] **Step 1: Arena-count intervention**

Run the reference workload at `workers=8` with glibc forced to a single arena,
then with the default:

```bash
MALLOC_ARENA_MAX=1 TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

Reading: if forcing one arena sharply degrades `min_s`, per-thread arenas are
load-bearing and allocator contention is real but *already mitigated* by
glibc. If the two are within noise, the allocator is not a scaling limiter at
this thread count — it is pure per-thread overhead, which is issue #26's
concern, not this investigation's.

- [ ] **Step 2: Allocator-swap intervention**

Add to `Cargo.toml` under `[dev-dependencies]`:

```toml
mimalloc = { version = "0.1", default-features = false }
```

Add to the top of `examples/bulk_bench.rs`, after the `use` block:

```rust
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
```

Then:

```bash
cargo build --release --example bulk_bench --features bulk
```

**If the build fails because crates.io is unreachable** (this is a cluster node
and network egress is not guaranteed), do not work around it: revert both
edits, record in the doc that the mimalloc arm was not run and why, and rely on
Step 1's result. That is an honest gap, not a failure.

If it builds, run the sweep at `workers ∈ {1, 4, 8}` and compare `S(4)` and
`S(8)` against the Task 2 BCF-unpinned series:

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=1 ./target/release/examples/bulk_bench
```

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=4 ./target/release/examples/bulk_bench
```

```bash
TMPDIR=$CLAUDE_JOB_DIR/tmp VCFIXTURE_BENCH_SAMPLES=2000 VCFIXTURE_BENCH_RECORDS=20000 VCFIXTURE_BENCH_REPS=3 VCFIXTURE_BENCH_WORKERS=8 ./target/release/examples/bulk_bench
```

- [ ] **Step 3: Revert both edits and prove the tree is clean**

```bash
git checkout -- Cargo.toml examples/bulk_bench.rs
```

```bash
git status --short
```

Expected: no output. Note that `Cargo.lock` may have changed; if
`git status --short` lists it, revert that too with
`git checkout -- Cargo.lock`.

- [ ] **Step 4: Record the findings and update issue #26**

Add an `### Allocator interventions` subsection to the measurements doc with
both readings. Then comment on issue #26 with whichever of these applies:

```bash
gh issue comment 26 --body "<the measured reading, stating whether the allocator limits scaling or is pure per-thread overhead>"
```

- [ ] **Step 5: Commit**

```bash
git add docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md
git commit -m "docs(bulk): test the allocator hypothesis by intervention (#22)

MALLOC_ARENA_MAX=1 versus default, and a mimalloc global-allocator swap,
both at the reference workload. Replaces the originally-proposed profile
self-time comparison, which was correlational.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Correct the published claims

**Files:**
- Modify: `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md`
- Modify: PR #28 body (via `gh pr edit`)

**Interfaces:**
- Consumes: the curves from Task 2 and whichever of Tasks 3–4 ran.
- Produces: the final, corrected public account. No further tasks depend on it.

- [ ] **Step 1: Rewrite the Hypothesis 2 verdict**

In `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md`, replace the
existing Hypothesis 2 verdict paragraphs — the ones asserting "43.5% parallel
efficiency", "implied serial fraction ≈ 18.6%", and "Both 'allocator arena-lock
contention' and 'rayon+bgzf thread oversubscription' are plausible, untested
explanations" — with the measured account. It must state:

- the corrected denominator (4 physical cores) and the corrected efficiency;
- what the three curves show;
- that the `workers=1` baseline was itself not single-threaded, since
  `workers=1` also sizes the bgzf pool to 1, overlapping compression with
  generation and thereby *understating* the measured speedup;
- whatever Tasks 3–4 established, or an explicit statement of what remains open
  and why, if a gap survives.

Leave the speedup-vs-baseline tables untouched: baseline and post-change runs
shared the same allocation, so those comparisons were never affected.

- [ ] **Step 2: Run the full suite on the final tree**

```bash
cargo test --all-features
```

Expected: 153 tests pass, 0 failures.

```bash
cargo fmt --check
```

```bash
cargo clippy --all-features --all-targets -- -D warnings
```

Expected: both clean.

- [ ] **Step 3: Prove no instrumentation survived**

```bash
git diff origin/main -- src/bulk/mod.rs Cargo.toml | rg -c 'instr|mimalloc|t_par|t_ser'
```

Expected: exit code 1 with no output (zero matches). If anything matches, a
temporary patch leaked into a commit — remove it before proceeding.

- [ ] **Step 4: Commit the doc**

```bash
git add docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md
git commit -m "docs(bulk): resolve the sub-linear scaling question (#22)

Replaces the 43.5%-efficiency claim and its two untested explanations with
a measured account against the true 4-physical-core allocation.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

- [ ] **Step 5: Push and correct the PR body**

```bash
git push
```

Then rewrite the PR #28 body's "Sub-linear scaling is not fully explained"
bullet to state the corrected finding. Write the new body to
`$CLAUDE_JOB_DIR/tmp/pr-body.md` first, then:

```bash
gh pr edit 28 --body-file $CLAUDE_JOB_DIR/tmp/pr-body.md
```

- [ ] **Step 6: Final process hygiene check**

```bash
pgrep -a bulk_bench
```

Expected: no output. Nothing this plan started may still be running.

---

## Self-Review

**Spec coverage.** Stage 0 → Task 2 Step 1/8. Stage 1 → Task 2 Steps 3–4, 7.
Stage 2 → Task 2 Step 5 (promoted to unconditional, justified in that task's
scope note). Stage 3 → Task 3. Stage 4 → Task 4. Harness changes → Task 1.
Deliverables 1–2 → Task 5; deliverable 3 → Task 1; deliverable 4 → Task 4
Step 4. Testing section → Task 1 Step 7 and Task 5 Steps 2–3.

**Placeholder scan.** No TBDs. Every code step carries the literal code; every
run step carries the literal command and its expected output. The one
deliberately open element is the *content* of the Task 5 rewrite, which cannot
be pre-written because it reports measurements that do not exist yet — Step 1
therefore specifies the four claims it must contain rather than the prose.

**Type consistency.** `Format` is imported in Task 1 Step 1 and used in Task 1
Step 2's `.format(format)`. `reps()` returns `usize`, consumed as
`Vec::with_capacity(reps)` and `0..reps`. `format()` returns
`(Format, &'static str)`, destructured as `let (format, ext)`.
`summary.n_records_total()` returns `u64`, multiplied by `n_samples` (`u64`
from `sweep`) — matching the original code's arithmetic. `t_par`/`t_ser` are
`Duration`, read via `.as_secs_f64()`.
