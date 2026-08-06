# Bulk generation: parallel encode and a positions-only span pass (issue #22)

**Status:** approved, ready for implementation
**Issue:** [#22](https://github.com/d-laub/vcfixture-rs/issues/22)
**Baseline:** v0.4.0 (`9e64a0c`)
**Breaking:** yes — generated output for a given seed changes, and
`summary.json`'s `genotype_checksum` changes meaning. Target v0.5.0.

## Problem

Generating a 6-corpus benchmark ladder (up to 32,000 samples × 1M variants)
took ~4.5 h on a dedicated 48-core allocation while `vcfixture` sat at **120 %
CPU** — ~46 cores idle throughout. Cost is dead linear in
`n_records × n_samples` at ~2.6e-7 s/cell across a 16× range of cell counts.

`bulk` is not missing parallelism. `generate_contig` already fans out over
`BLOCK_SIZE`-record blocks with rayon, `workers` already defaults to
`available_parallelism()`, and `BulkWriter` already wraps a multithreaded bgzf
writer. Both dominant costs are simply outside that fan-out.

### Finding 1 — every record is generated twice

`BulkSpec::write` is a two-pass structure: a span pass that generates each
contig only to compute `contig_span(&recs)`, then a write pass that regenerates
it. But `contig_span` reads only `r.g.pos`. The span pass therefore performs
the full `O(n_records × n_samples)` genotype synthesis — 1.6e10 genotype draws
for the 32,000 × 500,000 corpus — and discards everything except one `u64` per
contig.

The two-pass design is well motivated (it bounds peak memory to one contig
rather than the whole corpus) and is not what should change. What should change
is that pass 1 does the expensive half of the work at all.

**This is worse than the issue states.** The same span-pass-then-write-loop body
appears *three* times — `BulkSpec::write` (`mod.rs:431-468`),
`BulkSpec::write_to_temp` (`:680-710`), and
`BulkSpec::measure_compressed_bytes` (`:771-…`) — and `Size::Target`'s
two-point calibration calls `measure_compressed_bytes` twice
(`:622-623`) before calling `write_to_temp` once (`:639`). A `--target-size`
run therefore pays the double generation **three times**: six full corpus
generations to emit one file.

**Correction to the issue's premise.** The issue proposes a positions-only mode
on the grounds that "the positions are already a pure function of
`(seed, block_idx)`, so the determinism guarantee is unaffected" and that this
costs "no change to output bytes". That is not true of the current RNG layout:
`samplers.gap(&mut rng)` draws from the *same* block stream that `gen_record`
then advances by `n_samples` missingness draws plus the AC-placement draws
(`mod.rs:545-560`). A pass that skips genotype synthesis desynchronises the
stream, so every subsequent position differs. Positions can only be made
independent of `n_samples` by giving them their own stream — which does change
output bytes. That break is accepted here (see Compatibility).

### Finding 2 — per-record encoding is serial and scales with cohort width

The write loop is sequential per record:

```rust
for r in &recs {
    let buf = to_record_buf(&r.g, self.payload.clone(), r.phased);
    writer.write(&header, &buf)?;
    summary.observe(id, r.g.pos, r.g.class, &r.g.gts);
}
```

All three statements are `O(n_samples)` per record, on one thread:
`to_record_buf` materialises a `Vec<Option<Value>>` per sample (with a `String`
GT clone each), `writer.write` encodes all of them, and `summary.observe` folds
every allele into the checksum. bgzf *compression* downstream is multithreaded
and record *synthesis* upstream is block-parallel, but this stage — the one
that scales with cohort width — sits between them on a single core. That is
both the 120 % CPU and the flat s/cell.

### Finding 3 (not in the issue) — peak memory is a whole contig of records

`generate_contig` returns `Vec<Rec>` for an entire contig. At 32,000 samples
and ~45,000 records/contig that is ~2.9 GB of `gts` alone, which matches the
5.8 GB RSS observed. Fixing Finding 2 by streaming blocks fixes this too, and
must, or holding a contig's worth of *encoded* records would be worse.

### Already fixed upstream, not re-litigated here

v0.3.0 landed `perf(bulk): build GT into a reused String, not Vec<String>+join`
and `perf(bulk): rejection-sample sparse alt placement`. `SampleStats::new` now
writes GT into one `String`, and sparse AC placement no longer materialises a
`Vec<usize>` of non-missing indices. Neither is re-addressed.

## Non-goals

- No change to any sampling distribution, to realism, or to the fitted-profile
  schema. Record *content* changes only as a consequence of the RNG-stream
  split, never because a sampler changed.
- No hand-written BCF encoder. Records continue to go through
  `RecordBuf` + noodles' encoder; only *where* that runs changes.
- No further constant-factor work on `to_record_buf` (its per-sample `Vec` and
  `String` clone) unless a profile taken *after* this change shows it hot.
  Phase 4 decides that, not this spec.

## Design

### 1. Two RNG streams per block

`block_rng` gains an explicit stream domain rather than callers salting the
seed themselves:

```rust
pub enum Stream { Position, Content }
pub fn block_rng(seed: u64, block_idx: u64, stream: Stream) -> ChaCha8Rng;
```

The domain is folded into the splitmix64 finalizer alongside `block_idx` with
its own constant, so `(seed, block_idx, Position)` and
`(seed, block_idx, Content)` cannot alias for any `seed`/`block_idx` — an
explicit domain rather than a `seed ^ SALT` trick, which would only make
aliasing improbable rather than impossible.

A block then draws gaps from `Stream::Position` and everything else (class,
site, AC, missingness, placement, phasing) from `Stream::Content`. A block's
positions become a pure function of `(seed, block_idx, count)`, independent of
`n_samples` and of the payload — which is what the issue assumed was already
true.

### 2. `ContigLayout`: spans without generating genotypes

```rust
struct ContigLayout {
    block_spans: Vec<u64>,   // per-block local span (sum of that block's gaps)
}
```

`block_spans[i]` is the sum of block `i`'s gap draws. Contig span is their sum;
block `i`'s absolute position offset is their exclusive prefix sum. Both fall
out of one `Vec<u64>`.

Layouts for **all** contigs are computed in a single flat `par_iter` over
`(contig_idx, block_idx)` pairs that draws only from `Stream::Position` —
`O(total_records)` RNG draws with no per-sample work at all. For the 1M-record
ladder that is ~1M gap draws, milliseconds.

Memory is `n_blocks × 8` bytes per contig (≤ ~4 MB per contig even at the
smallest block size in §4), so all layouts are held at once without a bound
worth stating.

The redundant generation pass disappears entirely: **one** generation pass per
output file instead of two, and one per calibration probe instead of two.

### 3. `stream_contigs`: one implementation of the write loop

The three duplicated span-pass-plus-write-loop bodies collapse into one
helper:

```rust
fn stream_contigs(
    &self, pool: &rayon::ThreadPool, samplers: &Samplers, fitted: &Fitted,
    counts: &[u64], layouts: &[ContigLayout], header: &vcf::Header,
    writer: &mut BulkWriter,
) -> Result<Summary, BulkError>;
```

`write` and `write_to_temp` both call it. `measure_compressed_bytes` becomes
`write_to_temp(...)` discarding all but the byte count — three near-identical
bodies become one, which is the point of doing this refactor here rather than
duplicating the new block pipeline three times.

The two behavioural differences between the merged functions must be preserved
deliberately: `measure_compressed_bytes` builds no `Summary` (it now will, at
O(1) merges per block — negligible), and it best-effort removes the
`<tmp>.csi` companion that `NamedTempFile`'s `Drop` does not know about. That
cleanup moves to the calibration call site, which already performs exactly the
same cleanup for its own corrective rounds (`mod.rs:648-658`) — so the merge
removes a third copy of that too, rather than leaking temp files.

Per contig, blocks are processed in **chunks of `2 × workers`**, collected in
index order and drained:

```rust
for chunk in (0..n_blocks).collect::<Vec<_>>().chunks(chunk_blocks) {
    let out: Vec<(Vec<u8>, BlockSummary)> = pool.install(|| {
        chunk.par_iter().map(|&b| encode_block(b, ...)).collect()
    });
    for (bytes, bs) in out {
        writer.write_encoded(&bytes)?;
        summary.merge_block(id, bs);
    }
}
```

`.collect()` on a rayon parallel iterator preserves index order regardless of
which thread computed which item — the same property the existing block design
already relies on — so record order is unchanged. The serial half is now a
`write_all` memcpy into the bgzf staging buffer plus an O(1) summary merge.

`encode_block` runs entirely in the worker: draw gaps from `Stream::Position`,
add the block's absolute offset, synthesise records from `Stream::Content`,
`to_record_buf`, encode into a per-block byte buffer, and fold the block's own
`BlockSummary`. `self.payload` is borrowed, not cloned per record.

Encoding runs against an in-memory buffer. A header-*less* writer does not
work: `noodles_bcf::io::Writer` builds its `StringMaps` inside `write_header`
and keeps it in a private field, so `write_variant_record` on a fresh writer
fails with `chromosome not in string map`. (Empirically confirmed, not
assumed.) The working construction is a **per-worker** writer that wrote its
header once and rewinds per block:

```rust
// rayon `map_init`: once per worker thread, not once per block
let mut blk = bcf::io::Writer::from(Vec::new());
blk.write_header(&header)?;                 // populates StringMaps
let header_len = blk.get_ref().len();
// ... per block:
blk.get_mut().truncate(header_len);         // rewind, keep the string map
for r in &recs { blk.write_variant_record(&header, &buf)?; }
let bytes = &blk.get_ref()[header_len..];   // this block's records only
```

Header formatting therefore happens once per worker (it is not free: at 32,000
samples the header text is ~200 KB of sample names), and the record buffer is
reused across that worker's blocks. `vcf::io::Writer` takes the same shape for
`VcfGz`/`Vcf`.

`BulkWriter` gains `write_encoded(&[u8])`, which reaches the underlying sink
via `noodles_bcf::io::Writer::get_mut` (present in 0.81; `vcf::io::Writer` and
the plain-`Vcf` `File` sink get the same treatment).

**Risk:** the byte-identity claim is the load-bearing assumption of this whole
design. A scratch test has already confirmed it for both BCF and VCF text at
small scale; the plan promotes that into a permanent test against the real
generator — see Testing.

### 4. Block size scales with cohort width

`BLOCK_SIZE` (a flat 500 records) becomes a function of how wide the cohort is,
because a block's memory is cells, not records:

```rust
const TARGET_CELLS_PER_BLOCK: u64 = 4_000_000;
const MAX_BLOCK_RECORDS: u64 = 500;          // today's BLOCK_SIZE

fn block_records(n_samples: usize, ploidy: u8) -> u64 {
    (TARGET_CELLS_PER_BLOCK / (n_samples as u64 * ploidy as u64))
        .clamp(1, MAX_BLOCK_RECORDS)
}
```

At 32,000 diploid samples that is 62 records/block (~4 M cells, ~4 MB encoded
gt-only); at ≤ 4,000 samples it saturates at today's 500, so small-cohort
behaviour keeps its current shape. In-flight memory is then
`2 × workers × block bytes` ≈ 380 MB at 48 workers and 32,000 samples, against
5.8 GB today.

`TARGET_CELLS_PER_BLOCK` is a starting value to be confirmed in Phase 4 — small
enough to bound memory, large enough that per-block overhead (one `Vec`, one
writer, one RNG init) is negligible against 4 M cells of work.

Block count per contig rises as block size falls, which brings
`CONTIG_BLOCK_STRIDE` (1e6) within reach: at 1 record/block a contig may not
exceed 1e6 records without colliding into the next contig's block-index range.
Today that collision is **silent** — it would reuse another contig's streams.
It becomes an explicit `BulkError::TooManyBlocks { contig, n_blocks, stride }`.

### 5. Mergeable summary

`Summary::observe` (per record, serial) is replaced by per-block accumulation
plus an O(1) merge:

```rust
struct BlockSummary {
    n_records: u64, pos_min: u64, pos_max: u64,
    n_alleles_total: u64, n_alleles_nonref: u64,
    class_counts: [u64; VariantClass::COUNT],   // indexed by class, not by String
    checksum: u64,
}
```

Counts add, `pos_min`/`pos_max` min/max. `class_counts` becomes a fixed array
indexed by the `VariantClass` enum instead of a `BTreeMap<String, u64>` probed
by name once per record; `Summary`'s public serde shape keeps its
`BTreeMap<String, u64>`, so `summary.json`'s schema is unchanged.

`genotype_checksum` is order-sensitive FNV-1a and does not merge. It becomes
**hierarchical**:

- per block: FNV-1a over that block's allele bytes in record-then-slot order,
  from `FNV_OFFSET`;
- per file: FNV-1a over each block's 8 checksum bytes, in block order, from
  `FNV_OFFSET`.

This stays order-sensitive both within and across blocks (swapping two blocks
changes the folded sequence), and stays independent of `workers` and of the
chunk size, because block boundaries depend only on `block_records` and
`n_records` — never on how many blocks are in flight. The checksum *value*
changes, which is covered by the output-bytes break.

## Compatibility

**Breaking.** Same seed + same spec no longer reproduces v0.4.0's bytes, because
positions now come from their own stream and block boundaries now depend on
cohort width. `summary.json`'s `genotype_checksum` changes meaning as well.

What is preserved, and remains tested: same seed + same spec produces
byte-identical output regardless of `--threads`, and now also regardless of
chunking. That is the guarantee the book documents
(`docs/book/src/bulk-generation.md:99-101`); its wording needs a note that
byte-stability holds within a major version, not across one.

Requires: CHANGELOG `BREAKING CHANGE` entry, a v0.5.0 bump, and a book note.
Existing corpora must be regenerated once — which the same change makes
substantially cheaper.

## Testing

Correctness gates, all of which must pass before any tuning is kept:

1. **Serial-reference oracle.** A test-only reference path with the same block
   decomposition but `.map()` instead of `par_iter`, encoding inline. The
   optimized path must produce a **byte-identical file** and an **equal
   `Summary`**. This is what verifies the load-bearing header-less-encoding
   assumption in §3, and it verifies parallelism and chunking rather than the
   block math.
2. **Layout/realization agreement.** Spans now come from a pass that never
   generates genotypes, so a divergence would silently declare a wrong
   `##contig` length. Assert, across several seeds and sizes, that each
   contig's `ContigLayout` span equals the maximum position actually written,
   and that positions are strictly increasing across block boundaries.
3. **Thread- and chunk-independence.** Extend the existing
   `same_seed_gives_byte_identical_output_across_thread_counts` to also vary
   the chunk size, with the existing payload-size assertion kept so the test
   cannot be silently shrunk below a bgzf block boundary.
4. **Summary fold laws.** Merge is associative and identity-respecting;
   checksum detects a dropped record, a reordering within a block, and a
   reordering of blocks.
5. **Stream separation.** Positions for a given `(seed, contig, n_records)` are
   identical across different `n_samples` and different payloads.
6. **Stride overflow.** A contig exceeding `CONTIG_BLOCK_STRIDE` blocks errors
   rather than silently colliding.
7. Existing distribution/realism tests (SFS, missingness, class mix, phasing)
   must still pass unchanged — they are what confirms the RNG split did not
   perturb any sampler.

## Measurement plan

Per `performant-py-rust`; the harness lands **before** any optimization.

- **Phase 0 target.** Evidence: 120 % CPU on a 48-core allocation, s/cell flat
  (2.55–2.62e-7) across a 16× range of cell counts — serial per-cell work, not
  the block-parallel path. Target: no `O(cells)` stage outside the fan-out, so
  utilisation scales with `--threads`; and total generation work halved by
  removing the redundant pass.
- **Phase 1 dimensions.**

  | dimension | typical | max | grows? | notes |
  |-----------|---------|-----|--------|-------|
  | `n_samples` | 4,000 | 32,000 | yes, per study | sets per-record encode cost |
  | `n_records` | 250,000 | 1,000,000 | yes | sets block count |
  | cells (`n_records × n_samples × ploidy`) | 2e9 | 1.6e10 | **dominates** | cost is linear in this |
  | `n_contigs` | 22 | 22 | fixed | outer loop only |
  | payload keys | 1 (`gt-only`) | 7 (`mutect2`) | fixed per run | multiplies encode cost |
  | `workers` | 48 | 48 | machine | currently ~1.2 effective |

  Bound: CPU-bound, with a serial section. Lever: move the serial section into
  the existing data-parallel fan-out; do not add threads.
- **Phase 3 harness.** `examples/bulk_bench.rs`, sweeping `n_samples` and
  `n_records` independently and reporting **s/cell** and peak RSS, so numbers
  line up directly against the ladder in issue #22. Baseline recorded before
  any change.

  Note: this dev box has **8 cores**, so local numbers show the serial-section
  fix at 8× headroom, not 48×. Confirming the 48-core figure is a re-run of the
  ladder on the original allocation — the issue author has offered to do this
  and the corpora are content-cached.
- **Phase 4 loop.** Profile (`samply`) → one hypothesis → one change → re-run
  oracle *and* benchmark → keep or revert. Stop at the Amdahl ceiling or when
  added complexity outweighs the remaining gain, and state which.

## Expected outcome

Two hypotheses, to be confirmed by the sweep, not asserted:

1. Removing the redundant pass halves generation CPU (and cuts a
   `--target-size` run from six corpus generations to three).
2. Moving encode/summary into the fan-out removes the only `O(cells)` serial
   stage, so wall clock becomes bounded by generation + bgzf compression + I/O
   rather than by one core.

Peak RSS should fall from ~5.8 GB to a few hundred MB, bounded by
`2 × workers × block bytes` rather than by the largest contig.
