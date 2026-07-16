# Task 4 report: streaming writer + counting writer + second-pass index

## What was built

- `src/bulk/writer.rs` (new): `CountingWriter<W>`, `Format`, `BulkWriter`.
- `src/bulk/mod.rs`: added `pub mod writer;`.
- `Cargo.toml`: added `tempfile = "3"` to `[dev-dependencies]`.

Interfaces match the brief exactly, plus one deliberate addition:

```rust
pub struct CountingWriter<W> { .. }
impl<W: Write> CountingWriter<W> {
    pub fn new(inner: W) -> (CountingWriter<W>, Arc<AtomicU64>);
}

pub enum Format { Bcf, VcfGz, Vcf }

pub struct BulkWriter { .. }
impl BulkWriter {
    pub fn create(path: &Path, format: Format, header: &vcf::Header,
                  compression_level: u8, workers: NonZero<usize>) -> Result<BulkWriter, BulkError>;
    pub fn write(&mut self, header: &vcf::Header, record: &RecordBuf) -> Result<(), BulkError>;
    pub fn flush(&mut self) -> Result<(), BulkError>;               // ADDED, not in the brief
    pub fn compressed_bytes(&self) -> u64;
    pub fn finish_and_index(self, path: &Path) -> Result<(), BulkError>;
}
```

## TDD sequence

1. Wrote the brief's test verbatim (fixed one bug in it — see below), added `pub mod writer;`
   and `tempfile` dev-dep. `cargo test --features bulk writer` failed as expected:
   `CountingWriter`/`BulkWriter`/`Format` not found (4 compile errors).
2. Implemented per the brief's reference code, essentially unchanged (see "API verification"
   below — the brief's noodles calls all compiled as written). Re-ran: 1/2 tests passed,
   `writes_a_readable_indexed_bcf` failed on `assert!(w.compressed_bytes() > 0, ...)` — count
   was 0. Diagnosed (see below), fixed by adding `BulkWriter::flush` + a bounded poll in the
   test. Final run: both tests pass, and ran 5x in a row with no flakiness.
3. Added a third test, `output_is_byte_identical_regardless_of_worker_count`, to positively
   verify the "determinism regardless of thread count" global constraint (writes 500 records
   through `BulkWriter` with `workers=1` vs `workers=8`, asserts the resulting `.bcf` bytes are
   identical). Passes.

Final `cargo test --features bulk bulk::writer` output:

```
running 3 tests
test bulk::writer::tests::counting_writer_counts_bytes_through ... ok
test bulk::writer::tests::writes_a_readable_indexed_bcf ... ok
test bulk::writer::tests::output_is_byte_identical_regardless_of_worker_count ... ok

test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out
```

Full suite (`pixi run test` = `cargo test --all-features`): 73 lib unit tests + 5 roundtrip +
2 snapshot + 1 doctest, all green. `pixi run fmt` and `pixi run clippy` (`-D warnings`) clean.
Confirmed the default (no-features) build still excludes `noodles-bcf` (`cargo tree -e normal`).

## API verification against the real 0.81/0.45/0.53/0.83 crates

Before implementing, read the actual crate source under
`~/.cargo/registry/src/index.crates.io-*/noodles-{bcf,bgzf,csi,vcf}-<version>/src/` rather than
trusting the brief's uncompiled reference. Verified line-for-line:

- `bgzf::io::MultithreadedWriter::{new, with_worker_count, finish}` — exist as written.
  `finish(&mut self) -> io::Result<W>` flushes, joins the deflater and writer threads, appends
  the BGZF EOF block, and transitions internal state to `Done`.
- **`Drop for MultithreadedWriter` calls `finish()`** (`multithreaded_writer.rs:152-161`,
  `if !matches!(self.state, State::Done) { let _ = self.finish(); }`). So the brief's note "if
  it does not flush on drop, call `finish()` explicitly" turned out unnecessary — it does flush
  on drop in 0.45. `finish_and_index` relies on this: `drop(self.sink)` is sufficient, and the
  read-back-3-records test genuinely exercises and proves this (not skipped).
- `bcf::fs::index<P: AsRef<Path>>(src: P) -> io::Result<csi::Index>` — exists exactly as the
  brief wrote it (`src/fs/index.rs`), including the doc comment referencing `csi::fs::write`.
- `csi::fs::write<P: AsRef<Path>>(dst: P, index: &Index) -> io::Result<()>` — exists exactly as
  written.
- `impl<W> From<W> for bcf::io::Writer<W>` — unconstrained, exists as written, so
  `bcf::io::Writer::from(multithreaded_writer)` works.
- `vcf::io::Writer::<W>::new(inner: W) -> Self` — generic over any `W`, exists as written.
- `vcf::variant::io::Write` trait (`write_variant_header`/`write_variant_record`) — implemented
  by both `bcf::io::Writer<W: Write>` and `vcf::io::Writer<W: Write>`, exists as written.
- `noodles_bcf::io::reader::Builder::default().build_from_path(path)` — exists as written.

**Net result: the brief's reference `BulkWriter`/`CountingWriter` implementation compiled
against the pinned versions with zero API adaptation.** No noodles call needed to change.

## Where I did deviate from the brief, and why

1. **Test bug fix**: `header()` used
   `.add_contig("chr1", vcf::header::record::value::map::Contig::default())`, but
   `Header::builder().add_contig` takes `Map<Contig>`, not bare `Contig` — a real type mismatch
   in the brief's reference test. Fixed to
   `vcf::header::record::value::Map::<vcf::header::record::value::map::Contig>::new()`.

2. **`compression_level` is wired up, not ignored.** The brief's reference code took
   `_compression_level: u8` and discarded it. I used
   `bgzf::io::writer::CompressionLevel::try_from(compression_level)` and passed it through
   `bgzf::io::multithreaded_writer::Builder::default().set_worker_count(workers)
   .set_compression_level(level).build_from_writer(counting)` (the brief's
   `MultithreadedWriter::with_worker_count` always uses the default level 6 — there is no
   builder-free way to set both worker count and compression level, so I went through
   `multithreaded_writer::Builder` directly, which exists and is exactly what
   `with_worker_count` uses internally). An invalid level (out of 0..=9, since this build has no
   `libdeflate` feature so `flate2::Compression` caps at 9) now returns
   `BulkError::Invalid(..)` instead of silently being dropped on the floor — a public numeric
   knob that's silently ignored is exactly the kind of surprising API the coding principles
   ("make invalid states unrepresentable", fail fast) push against.

3. **Added `pub fn flush(&mut self) -> Result<(), BulkError>` to `BulkWriter`** — not in the
   brief's interface list. This was forced by a real behavioral finding, detailed next.

## The `compressed_bytes() > 0` investigation (the one substantive deviation)

The brief's test writes only 3 tiny records then immediately asserts
`w.compressed_bytes() > 0`, before calling `finish_and_index`. With the brief's code verbatim,
this failed: count was 0.

Root cause (confirmed by reading `multithreaded_writer.rs` and by an instrumented run with
`eprintln!` + `thread::sleep(200ms)`, which showed the count become 251 after the sleep):
`MultithreadedWriter::write()` only dispatches a block to the compression/writer threads once
its ~64 KiB uncompressed staging buffer is full (`if !self.has_remaining() { self.flush()?; }`),
and even an explicit `flush()` call only *enqueues* the block onto `crossbeam_channel`s — actual
compression + write to the inner `CountingWriter` happens asynchronously on background threads.
3 small BCF records (~50-100 bytes total) never fill the buffer, so nothing is ever dispatched,
and even if it were, dispatch completion is inherently async relative to the caller.

I considered and rejected making `write()` auto-flush after every record: that would turn every
single record into its own BGZF block (each with ~18-28 bytes of gzip/BGZF header overhead) and
its own compression job, which at "benchmark scale" (hundreds of thousands of records for a
100 MB target) would (a) bloat output size significantly by never letting DEFLATE exploit
cross-record redundancy within a block, and (b) add per-record thread-handoff overhead — directly
contradicting the design spec's stated goal that multithreaded batched compression is what the
throughput target depends on.

Instead: added `BulkWriter::flush()`, matching the design spec's own language ("the generator
polls between record blocks") — the not-yet-built size-targeting loop (a sibling task) is meant
to call `flush()` once per *batch* of records before checking `compressed_bytes()`, not after
every record. Docstrings on both `write()` and `flush()` explain this trade-off explicitly so
whoever builds the generator loop doesn't call `flush()` per-record by mistake.

Because dispatch is still asynchronous even after `flush()` returns, the test itself now calls
`w.flush()` once after the write loop, then polls `compressed_bytes()` in a tight loop with a
2-second bounded timeout (checking every 1ms) rather than assuming instantaneous consistency —
this is the honest way to test an async pipeline, and it ran clean across 5 consecutive runs. In
real usage at scale the async lag is irrelevant: by the time a caller has written enough records
to be worth polling, the 64 KiB buffer will already have auto-flushed multiple times.

## Files touched

- `src/bulk/writer.rs` (new) — `/carter/users/dlaub/projects/vcfixture-rs/.claude/worktrees/bulk-task-4/src/bulk/writer.rs`
- `src/bulk/mod.rs` — added `pub mod writer;`
- `Cargo.toml` — added `tempfile = "3"` dev-dependency
- `Cargo.lock` — updated by cargo for the new dev-dependency

## Fixes: review round 1

Covering test file: `src/bulk/writer.rs`, `#[cfg(test)] mod tests`.

### Important 1 — near-vacuous determinism test

`output_is_byte_identical_regardless_of_worker_count` previously wrote 500 tiny
records (~19,099 uncompressed bytes total per the reviewer's measurement) — well
under `MultithreadedWriter`'s `MAX_BUF_SIZE` (~65,498 uncompressed bytes/block),
and never called `flush()`, so the whole run became a single bgzf block dispatched
only on drop. With one block in flight, no worker count could produce reordering,
so the test could not fail even if a real bug existed.

Fix: the test now writes 1,000 records with a 400-byte padded ALT allele, and
independently measures the *true* uncompressed record payload by encoding the
same records into an in-memory, uncompressed `bcf::io::Writer<Vec<u8>>` (so the
size check doesn't depend on how well bgzf's DEFLATE happens to compress this
repetitive test data). It asserts that payload exceeds `3 * MAX_BUF_SIZE`
(196,494 bytes) before comparing outputs, so a future shrink of the payload
fails loudly instead of silently going vacuous again. Measured payload in this
run: **440,099 bytes (~6.72x MAX_BUF_SIZE)** — comfortably several block
boundaries. The test now compares workers=1 vs 4 vs 16 (previously 1 vs 8) and
asserts byte-for-byte equality, per the brief.

This exercises real concurrency and it passes: the underlying dispatch mechanism
is genuinely order-preserving, confirming the reviewer's independent 50k-record
finding.

### Important 2 — timing-dependent poll in `writes_a_readable_indexed_bcf`

Removed the bounded 2s poll loop entirely (no more `Instant::now()` / `sleep`
race). `finish_and_index` drops `self.sink`, whose `MultithreadedWriter::drop`
calls `finish()`, which synchronously joins the deflater/writer threads before
returning — so dispatch is guaranteed complete once `finish_and_index` returns.
The test now clones the private `count: Arc<AtomicU64>` field (accessible from
the same module via `use super::*`) before calling `finish_and_index`, and
asserts `count.load(Ordering::Relaxed) > 0` afterward — fully deterministic, no
timing dependency. Also added a `std::fs::metadata(&path)?.len() > 0` assertion
for the same reason. No separate poll-during-write test was added since none of
the deterministic alternatives available without a production API change would
avoid reintroducing a race; the design's poll-based `Size::Target` loop (Task 7)
is not blocked by this since `flush()` itself is unchanged and still exercised
in this test.

### Minor 3 — `compression_level` validated even for `Vcf`

Doc-only fix (preferred option in the review): `BulkWriter::create`'s doc
comment now states that `compression_level` is validated unconditionally, even
for `Vcf` where it goes unused, so code and docs agree. No behavior change.

### Commands run

```
cd /carter/users/dlaub/projects/vcfixture-rs/.claude/worktrees/bulk-task-4
cargo build --features bulk --tests
cargo test --features bulk writer -- --nocapture
pixi run fmt
pixi run clippy
pixi run test
```

`cargo build --features bulk --tests`: clean, no warnings.

`pixi run fmt`: clean (no diff after running; formatting was already applied by
the same command during development).

`pixi run clippy` (`cargo clippy --all-features -- -D warnings`): clean.

`pixi run test` (`cargo test --all-features`, full suite): 74 lib unit tests
(including all 3 in `bulk::writer::tests`) + 5 roundtrip + 2 snapshot + 1
doctest, all green.

### 10 consecutive repeat runs of the covering tests (flakiness evidence)

Command: `cargo test --features bulk writer` (run 10 times in a row, back to
back, no changes between runs):

```
=== run 1 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.17s
=== run 2 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.17s
=== run 3 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.18s
=== run 4 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.16s
=== run 5 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.17s
=== run 6 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.16s
=== run 7 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.16s
=== run 8 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.17s
=== run 9 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.16s
=== run 10 ===
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 71 filtered out; finished in 0.17s
```

10/10 green, all 3 tests passing each time (`counting_writer_counts_bytes_through`,
`writes_a_readable_indexed_bcf`, `output_is_byte_identical_regardless_of_worker_count`),
each full run completing in ~0.16-0.18s. The beefed-up determinism test (1,000
records x3 worker counts x1 in-memory measurement pass) added negligible wall
time versus the prior vacuous version — no perceptible slowdown at this scale.

Note: the reviewer separately confirmed at 50k padded records (workers=1 vs 16,
3 runs) that the underlying mechanism is genuinely order-preserving; this
fix makes our committed test actually exercise that property instead of
trivially passing.
