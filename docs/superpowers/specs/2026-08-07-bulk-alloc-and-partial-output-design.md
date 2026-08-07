# Bulk generation: cut per-sample allocation, stop leaving partial output

Closes #26 and #27. Stacks on PR #28 (branch `worktree-fix-22-bulk-perf`).

## Why

Issue #26 measured glibc allocator work at ~47% of self time in the bulk
generation profile — a larger share than `vcfixture`'s own logic (~20%) and
noodles' BCF encoder (~16%) combined. PR #28's follow-up investigation then
established, by intervention rather than by profile, that this is **per-thread
allocation overhead, not arena contention**: swapping to mimalloc made every
worker count 1.63x–1.73x faster in absolute wall clock while leaving the
scaling curve's shape unchanged.

That splits the problem cleanly in two, and both halves are worth taking:

- Each allocation is more expensive than it needs to be. An allocator swap
  fixes that without touching a line of generation code.
- There are far more allocations than the work requires. Only source changes
  fix that, and they help every consumer regardless of allocator.

Issue #27 is unrelated to performance but lives in the same call path: a
mid-stream encode or write failure leaves a truncated, un-indexed file at the
caller's destination path. The repository already contains the machinery to
prevent this and uses it for exactly one of the four `Size` variants.

## Scope

Three parts, independent enough to land and revert separately.

| Part | Issue | Surface |
|---|---|---|
| A — scratch-buffer reuse | #26 | `src/bulk/generate.rs`, `src/bulk/mod.rs` |
| B — mimalloc for binaries | #26 | `Cargo.toml`, `src/bin/*.rs`, `examples/bulk_bench.rs` |
| C — atomic destination | #27 | `src/bulk/mod.rs` |

Out of scope: writing a bespoke BCF record encoder to bypass `RecordBuf`
entirely. That would remove the last remaining per-sample allocation (noodles'
internal `Vec<i8>` in `encode_genotype`), but it means reimplementing encoder
internals against a moving upstream API, and Part A already recovers three of
the four allocations at a fraction of the risk.

## Part A — reuse a per-thread scratch record

### The allocation inventory

For `Payload::GtOnly`, each sample of each record currently costs four heap
allocations:

1. `String::with_capacity(alleles.len() * 2)` in `SampleStats::new`
   (`generate.rs:297`), filled digit by digit with `write!(gt, "{a}")`.
2. `self.gt.clone()` in `SampleStats::value_for` (`generate.rs:329`). The
   original `SampleStats` is dropped immediately afterward, so this copy is
   pure waste.
3. The inner `Vec<Option<Value>>` produced by the per-sample `.collect()`
   (`generate.rs:379`).
4. noodles' `Vec<i8>` inside `encode_genotype_str`, built while *reparsing*
   the string that (1) just formatted.

At the reference workload — 2000 samples x 20000 records — that is roughly
160 million allocations. Payloads with array-valued keys are worse:
`Payload::Mutect2` adds five more `Vec`s per sample.

Items 1–3 are recoverable. Item 4 is inside noodles and stays.

### The mechanism

noodles exposes `impl From<Samples> for (Keys, Vec<Vec<Option<Value>>>)`. A
`RecordBuf`'s sample block can therefore be destructured into its parts,
refilled, and reassembled with `Samples::new` — with the outer `Vec`, every
inner `Vec<Option<Value>>`, and every `Value`'s own backing buffer keeping
their capacity across records. `Genotype` implements `AsMut<Vec<Allele>>` and
each `Array` variant owns its `Vec` directly, so both are clearable in place.

A `RecordScratch` type in `generate.rs` owns one `RecordBuf` and is carried
across records. `stream_contigs` already uses `rayon`'s `map_init`
(`mod.rs:741`), whose init closure is the natural home for one scratch per
worker thread. `to_record_buf(&GenRecord, &Payload, bool) -> RecordBuf` is
replaced by `RecordScratch::fill(&mut self, &GenRecord, bool) -> &RecordBuf`.

The refill rule is uniform across every FORMAT key: match on the slot's
existing `Value` variant; if it matches the variant this key needs, clear its
buffer and refill; otherwise assign a fresh `Value`. After the first record of
a block, every slot matches, so steady-state allocation for items 1–3 is zero.

Record-level fields (`chrom`, `ref_`, `alts`, `pos`) are set by plain
assignment and cloning. They cost about three allocations per *record*, some
60k over the reference workload against ~160M for the sample path — chasing
them would add code for no measurable gain.

### GT stays a string — attempted and rejected

**Original design, since abandoned.** `SampleStats.gt: String` was to become a
structured `Genotype` yielding `Value::Genotype` rather than `Value::String`,
which would additionally have removed our `write!("{a}")` integer formatting
from the ~20% bucket and noodles' `encode_genotype_str` reparse from the ~11%
bucket, since the BCF encoder dispatches `Value::Genotype` straight to
`encode_genotype`.

**Why it was rejected.** It changes the *text* VCF output. `build_header`
declares `VCFv4.5`, so noodles' text writer takes its VCF-4.4-and-later
genotype branch (`io::writer::record::samples::sample::value::genotype`),
which writes a phasing separator before **every** allele including position 0.
A diploid `0|0` renders as `/0|0`. The leading indicator is 4.4-conformant in
principle, but it is not what this crate emitted before and not what consumers
of these fixtures expect.

This was caught by the byte-equality gate, in exactly the shape the gate was
built to detect: **all four BCF goldens passed while all eight VCF and VCF.gz
goldens failed.** That split is itself the evidence that the two BCF encoders
agree — the byte-compatibility reasoning about the phase bit on position 0 was
correct — and that the problem is confined to the text writer.

**What ships instead.** GT remains a `Value::String` whose buffer is cleared
and refilled rather than reallocated. This preserves the entire allocation
win, which is what issue #26 is actually about; only the CPU saving from
skipping the format/reparse round-trip is forgone. A single-digit fast path
(`gt.push((b'0' + a) as char)` for `a < 10`) recovers the `core::fmt` /
`pad_integral` portion of that without touching the output.

The allocation count is unchanged by this decision: four to one either way.

### Expected effect

Steady-state allocations per sample per record fall from four to one, and the
format/reparse round-trip disappears. No wall-clock prediction is made here;
the measurement plan below is what settles it.

## Part B — mimalloc as the binaries' global allocator

Add `mimalloc` behind a `mimalloc` feature, **on by default**, and install it
as `#[global_allocator]` in `src/bin/vcfixture.rs`, `src/bin/validate_profile.rs`,
and `examples/bulk_bench.rs`.

Never in `src/lib.rs`. A library that sets a global allocator imposes it on
every dependent binary in the graph, which is not a library's decision to
make. Placing it in the binaries confines it to the surface where the cost was
measured, and `--no-default-features` restores a pure-Rust build for anyone
who does not want a `cc` toolchain dependency.

This is orthogonal to Part A: Part A removes allocations, Part B makes the
survivors cheaper. Neither subsumes the other, and the measurement plan reads
them separately.

## Part C — never leave a partial file at the destination

### Current behaviour

`BulkSpec::write` creates the destination through `BulkWriter::create`
(`mod.rs:529-535`), streams into it, and only then calls `finish_and_index`
(`mod.rs:537`). If `stream_contigs` returns `Err`, the `?` propagates and
`finish_and_index` never runs. The caller does observe the error — this is not
a silent-success bug — but a file remains at `path` with no BGZF EOF block, no
`.csi`, and no `.summary.json`, and nothing removes it.

`Size::Target` does not have this problem. It already generates into a
`NamedTempFile` via `write_to_temp` and moves it into place with
`promote_temp` (`mod.rs:516`), and `write_to_temp` is documented as byte-exact
against `write`'s own path.

### Change

Route every `Size` variant through `write_to_temp` + `promote_temp`. The
destination is written by a single atomic rename after the output is complete
and indexed; a failure anywhere upstream leaves the destination untouched.
This also removes the duplicated create/stream/finish sequence — the two paths
become one.

**Temp files must be created in the destination's parent directory**, via
`NamedTempFile::new_in`, not `NamedTempFile::new`. The default lands in
`TMPDIR`, which is routinely a `tmpfs` or a different filesystem from the
output. `persist` would then fall back to `std::fs::copy` of a bulk-scale
file — doubling I/O, and potentially failing outright on a `tmpfs` too small
to hold it. Creating the temp alongside the destination makes promotion a
same-filesystem rename in every case. This applies to the `Size::Target`
calibration probes as well, which are the same order of size as the output.

Cleanup of the `.csi` companion — which `NamedTempFile`'s `Drop` does not know
about — currently exists as two duplicated best-effort blocks, one in
`resolve_target_counts` (`mod.rs:913-917`) and one in `measured_bytes`
(`mod.rs:985-990`), the second explicitly commented as "mirroring" the first.
Both move next to the temp creation so every temp is cleaned the same way and
a third call site cannot forget. Routing `write` through the same path adds a
third temp that would otherwise need its own copy of this block.

### Deliberate non-goals

Concurrent writers racing to the same destination are not addressed. Rename is
atomic, so the destination is never torn, but last-writer-wins is unchanged
and no locking is introduced.

## Testing

### Byte-equality gate (blocking)

Parts A and C must not change one byte of output. Before touching any source,
capture SHA-256 of the generated artifacts at a fixed seed for **all three
formats** — BCF, VCF, and VCF.gz — and for **each `Payload` preset**, since
`GtOnly` exercises only the genotype path while `Gatk` and `Mutect2` exercise
the array-valued refill branches. Compare after each part.

Covering all three formats is what tests the text writer's rendering of
`Value::Genotype`, which is the one risk in Part A that reading the code
cannot settle. A mismatch here is a design failure, not a test to adjust.

### Behavioural tests

- **#27 regression test.** Induce a mid-stream write failure and assert that
  nothing exists at the destination path afterwards — no output file, no
  `.csi`, no `.summary.json`. Without this the fix is unfalsifiable, since the
  success path looks identical either way.
- **Cross-filesystem promotion.** Assert the temp is created beside the
  destination, so the copy fallback is not silently reintroduced.
- The existing `tests/bulk.rs` suite must pass unchanged. Any test requiring
  amendment is a signal that behaviour changed, and must be justified rather
  than updated.

### Measurement

Perf claims follow the methodology PR #28 established, and its cautions are
binding:

- Efficiency is computed against **4 physical cores**, not the 8 logical CPUs
  `nproc` reports.
- **Round-robin the worker counts across passes** rather than running all reps
  of a cell back to back. Running reps consecutively lets one noise burst
  contaminate all of them identically, which is what defeats a min-of-N
  estimator — the failure that inflated PR #28's anchor by 8.6%.
- **Re-check the `workers=1` anchor** with time-separated passes before
  dividing anything by it.
- Benchmarks run alone. No parallel CPU work.

Report Part A and Part B as separate deltas against the same baseline, plus
the combination. Part B's contribution is already known to be 1.63x–1.73x;
if the combined figure does not exceed it, Part A did not pay for itself and
that should be reported plainly rather than absorbed into a joint number.

## Risks

| Risk | Handling |
|---|---|
| ~~`Value::Genotype` renders differently in text VCF~~ | **Occurred.** Caught by the gate: 4 BCF goldens passed, 8 text goldens failed. `Value::Genotype` abandoned; see "GT stays a string" above |
| Scratch reuse leaks state between records (stale slot from a previous record) | Uniform clear-then-refill for every key every record; byte-equality gate across all payload presets |
| `NamedTempFile::new_in` needs the destination's parent to exist | `BulkWriter::create` already required a writable destination, so this is not a new precondition; surface a clear error if the parent is missing |
| mimalloc's `cc` dependency breaks a consumer's build | Feature is opt-out via `--no-default-features`; library target never sets a global allocator |
| Part A's win is smaller than the diff's cost | Measured and reported separately from Part B, so it can be judged on its own and reverted independently |

## Constraints carried from PR #28

- Never hand-edit `CHANGELOG.md` or the version in `Cargo.toml`; `cz bump`
  generates both in CI.
- Do not bypass prek hooks with `--no-verify`.
- Benchmarks must not run concurrently with other CPU-consuming work.
