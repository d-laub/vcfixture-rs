# Bulk allocation reduction and atomic destination writes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut steady-state per-sample allocations in bulk generation from four to one (#26), make each surviving allocation cheaper via mimalloc in the binaries (#26), and stop leaving truncated output at the destination path on mid-stream failure (#27).

**Architecture:** A per-thread `RecordScratch` carries one `RecordBuf` across records, refilling its sample buffers in place instead of reallocating them; GT becomes a structured `Value::Genotype` so both our integer formatting and noodles' string reparse disappear. mimalloc is installed as `#[global_allocator]` in the binary and example targets only, never the library. Every `Size` variant routes through the existing `write_to_temp`/`promote_temp` pair, with temps created beside the destination so promotion is always a same-filesystem rename.

**Tech Stack:** Rust 2021 (rust-version 1.86), noodles-vcf 0.83, noodles-bcf 0.81, rayon, tempfile, insta (dev), mimalloc 0.1 (new, optional).

**Spec:** `docs/superpowers/specs/2026-08-07-bulk-alloc-and-partial-output-design.md`

## Global Constraints

- Parts A and C must not change one byte of generated output. Task 1's golden gate is the arbiter; a mismatch is a design failure, not a test to update.
- `#[global_allocator]` must never appear in `src/lib.rs` or any module it compiles. Binaries and examples only.
- Never hand-edit `CHANGELOG.md` or the `version` field in `Cargo.toml`. `cz bump` generates both in CI.
- Do not bypass prek hooks with `--no-verify`. If `prek` hooks are not installed, run `prek install` before the first commit.
- Benchmarks must not run concurrently with each other or with any other CPU-consuming work. The Slurm allocation is **4 physical cores** (`0-3`) plus SMT siblings (`48-51`); `nproc` reports 8 logical CPUs and is the wrong denominator.
- Never detach a process from the session — no `nohup`, `setsid`, `disown`, or trailing `&`. Long benchmark runs use the harness's background mode so they stay tracked and reaped.
- Write all scratch and benchmark output to `$CLAUDE_JOB_DIR/tmp`, never under `~/.claude` or any other NFS path.
- Every task ends with `cargo test --all-features` green (153 tests at branch point) plus `cargo fmt --check` and `cargo clippy --all-features -- -D warnings` clean.

## Task Dependency Graph

```
Task 1 (golden gate) ──┬── Task 3 (Part A: scratch reuse) ──┐
                       └── Task 4 (Part C: atomic dest)  ───┼── Task 5 (measure + docs)
Task 2 (Part B: mimalloc) ─────────────────────────────────┘
```

**Parallelism:** Tasks 1 and 2 touch disjoint files and may run concurrently. After Task 1 lands, Tasks 3 and 4 may run concurrently — they both modify `src/bulk/mod.rs` but own disjoint regions (Task 3 owns the `map_init` block at `mod.rs:741-783`; Task 4 owns `write` at `mod.rs:505-545` and the temp helpers at `mod.rs:860-1015`). Use `superpowers:dispatching-parallel-agents` with `superpowers:subagent-driven-development`. Implementation subagents use Sonnet or weaker; reserve stronger models for fix rounds where an implementer critically failed.

---

### Task 1: Golden byte-equality gate

The blocking prerequisite for Tasks 3 and 4. Neither may start until this is committed, because without a pinned golden the "output is unchanged" claim is unfalsifiable.

**Files:**
- Create: `tests/bulk_golden.rs`
- Create: `tests/snapshots/` (insta writes these; commit them)

Touches no shared file, so it cannot conflict with Task 2.

**Interfaces:**
- Consumes: nothing.
- Produces: committed insta snapshots that Tasks 3 and 4 verify against. No Rust API.

**Why FNV-1a and not sha2:** this gate detects accidental change, not adversarial collision, so a 5-line inline hash is sufficient and avoids adding a dependency to `Cargo.toml` — which is what keeps this task parallel with Task 2. The snapshot pins byte length alongside the digest, so a change must collide *and* preserve length to slip through.

- [ ] **Step 1: Write the golden test**

Create `tests/bulk_golden.rs`:

```rust
//! Golden byte-equality gate for bulk generation.
//!
//! Pins a digest of the generated artifact for every (format, payload)
//! combination. Refactors that are supposed to preserve output — the
//! scratch-buffer reuse of #26 and the temp-then-promote change of #27 —
//! must leave every one of these snapshots untouched.
//!
//! A snapshot change is a design failure, not a test to update. If noodles
//! is upgraded and the encoding legitimately changes, update these with an
//! explicit note in the commit message saying which upstream change caused
//! it.

#![cfg(feature = "bulk")]

use std::num::NonZero;
use std::path::Path;

use vcfixture::bulk::{BulkSpec, Format, Payload, Profile, Size};

/// Cohort width. `BulkSpec::block_records(8, 2)` is 500 records at this
/// width, and `RECORDS_PER_CONTIG` below is a multiple of it, so every run
/// spans several blocks and actually exercises the `map_init` path that
/// Task 3 rewrites. A single-block run would be structurally incapable of
/// catching a scratch-reuse bug that leaks state between blocks.
const SAMPLES: usize = 8;
const RECORDS_PER_CONTIG: u64 = 1200;

/// FNV-1a, 64-bit.
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

fn digest(path: &Path) -> String {
    let bytes = std::fs::read(path).expect("output file must exist");
    format!("{} bytes, fnv1a64={:016x}", bytes.len(), fnv1a64(&bytes))
}

fn generate(format: Format, payload: Payload, name: &str) -> String {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join(name);
    BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
        .samples(SAMPLES)
        .contigs(["chr1", "chr2"])
        .payload(payload)
        .format(format)
        .size(Size::RecordsPerContig(RECORDS_PER_CONTIG))
        .seed(42)
        .workers(NonZero::new(2).unwrap())
        .write(&path)
        .unwrap();
    digest(&path)
}

macro_rules! golden {
    ($test_name:ident, $format:expr, $payload:expr, $file:expr) => {
        #[test]
        fn $test_name() {
            insta::assert_snapshot!(generate($format, $payload, $file));
        }
    };
}

// BCF exercises `encode_genotype` / `encode_genotype_str`; VCF and VcfGz
// exercise the *text* writer's rendering of the genotype value, which is
// the one Part A risk that reading the encoder source cannot settle.
golden!(bcf_gt_only, Format::Bcf, Payload::GtOnly, "a.bcf");
golden!(bcf_gt_vaf, Format::Bcf, Payload::GtVaf, "a.bcf");
golden!(bcf_gatk, Format::Bcf, Payload::Gatk, "a.bcf");
golden!(bcf_mutect2, Format::Bcf, Payload::Mutect2, "a.bcf");

golden!(vcf_gt_only, Format::Vcf, Payload::GtOnly, "a.vcf");
golden!(vcf_gt_vaf, Format::Vcf, Payload::GtVaf, "a.vcf");
golden!(vcf_gatk, Format::Vcf, Payload::Gatk, "a.vcf");
golden!(vcf_mutect2, Format::Vcf, Payload::Mutect2, "a.vcf");

golden!(vcfgz_gt_only, Format::VcfGz, Payload::GtOnly, "a.vcf.gz");
golden!(vcfgz_gt_vaf, Format::VcfGz, Payload::GtVaf, "a.vcf.gz");
golden!(vcfgz_gatk, Format::VcfGz, Payload::Gatk, "a.vcf.gz");
golden!(vcfgz_mutect2, Format::VcfGz, Payload::Mutect2, "a.vcf.gz");
```

- [ ] **Step 2: Run and accept the snapshots**

```bash
cargo test --all-features --test bulk_golden 2>&1 | tail -20
```

Expected: 12 tests fail on first run with "snapshot assertion failed" / new snapshot. Then:

```bash
cargo insta accept --workspace
```

If `cargo insta` is not installed, rename each `tests/snapshots/*.snap.new` to `*.snap` by hand instead — do not `cargo install` inside this session.

- [ ] **Step 3: Verify the gate is real, not vacuous**

A golden that passes no matter what is worse than no golden. Prove it bites:

```bash
# Temporarily perturb the generated output.
sed -i 's/\.seed(42)/.seed(43)/' tests/bulk_golden.rs
cargo test --all-features --test bulk_golden 2>&1 | tail -20
```

Expected: **all 12 tests FAIL.** If any test passes with a different seed, that test is not actually reading the generated bytes — fix it before proceeding. Then revert:

```bash
sed -i 's/\.seed(43)/.seed(42)/' tests/bulk_golden.rs
cargo test --all-features --test bulk_golden 2>&1 | tail -5
```

Expected: 12 passed.

- [ ] **Step 4: Confirm determinism across repeat runs**

```bash
cargo test --all-features --test bulk_golden 2>&1 | tail -5
cargo test --all-features --test bulk_golden 2>&1 | tail -5
```

Expected: 12 passed both times. A flaky golden would make Tasks 3 and 4 unreviewable.

- [ ] **Step 5: Full suite and lints**

```bash
cargo test --all-features 2>&1 | tail -20
cargo fmt --check && cargo clippy --all-features -- -D warnings
```

Expected: 153 pre-existing tests plus 12 new = 165 passed, 0 failed. Lints clean.

- [ ] **Step 6: Commit**

```bash
git add tests/bulk_golden.rs tests/snapshots
git commit -m "test(bulk): pin golden digests for every format and payload

Byte-equality gate for the #26 and #27 refactors, which must both be
output-preserving. Covers all three formats so the text VCF writer's
rendering of the genotype value is pinned, not just BCF's encoder.

Refs #26, #27"
```

---

### Task 2: mimalloc as the binaries' global allocator (Part B)

Independent of every other task. May run concurrently with Task 1.

**Files:**
- Modify: `Cargo.toml:13-31` (dependency), `Cargo.toml:39-43` (features)
- Modify: `src/bin/vcfixture.rs` (top of file)
- Modify: `src/bin/validate_profile.rs` (top of file)
- Modify: `examples/bulk_bench.rs` (top of file)
- Create: `docs/book/src/` note — see Step 5

**Interfaces:**
- Consumes: nothing.
- Produces: a `mimalloc` cargo feature, on by default. No Rust API.

**Note on the default-on cost.** Adding `mimalloc` to `default` puts the crate in the dependency graph for anyone depending on `vcfixture` without `default-features = false`, so library consumers do pay a compile-time `cc` cost even though the library never installs the allocator. That is the accepted trade for the CLI getting the measured 1.63x–1.73x automatically. The opt-out is one line and must be documented in Step 5.

- [ ] **Step 1: Add the dependency and feature**

In `Cargo.toml`, add to `[dependencies]`:

```toml
mimalloc = { version = "0.1", default-features = false, optional = true }
```

and change `[features]`:

```toml
[features]
default = ["mimalloc"]
proptest = ["dep:proptest"]
bulk = ["dep:noodles-bcf", "dep:serde", "dep:serde_json", "dep:rayon", "dep:tempfile"]
cli = ["bulk", "dep:clap"]
mimalloc = ["dep:mimalloc"]
```

- [ ] **Step 2: Install it in the binaries and the bench example**

Add to the top of `src/bin/vcfixture.rs`, `src/bin/validate_profile.rs`, and `examples/bulk_bench.rs`, immediately after the module doc comment and before the `use` statements:

```rust
// mimalloc, not glibc malloc: bulk generation is allocation-dominated
// (~47% of profile self time, issue #26) and swapping the allocator
// measured 1.63x-1.73x lower wall clock at every worker count. Installed
// here rather than in the library because a library that sets a global
// allocator imposes it on every dependent binary in the graph, which is
// not a library's decision to make.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
```

- [ ] **Step 3: Verify it is absent from the library**

```bash
rg -n "global_allocator" src/ examples/
```

Expected: exactly three hits — `src/bin/vcfixture.rs`, `src/bin/validate_profile.rs`, `examples/bulk_bench.rs`. **Zero hits under `src/bulk/`, `src/lib.rs`, or any other library module.** If a hit appears anywhere else, remove it.

- [ ] **Step 4: Verify both feature states build**

```bash
cargo build --all-features 2>&1 | tail -5
cargo build --no-default-features --features cli 2>&1 | tail -5
cargo test --all-features 2>&1 | tail -20
cargo fmt --check && cargo clippy --all-features -- -D warnings
```

Expected: both builds succeed, tests pass, lints clean. The `--no-default-features` build is the load-bearing check — it proves the opt-out actually works and that no code references `mimalloc` unconditionally.

- [ ] **Step 5: Document the opt-out**

Find the installation or build section of the user-facing docs:

```bash
rg -ln "no-default-features|cargo install|\[features\]" docs/book/src/ README.md
```

Add to whichever file documents features or installation:

```markdown
### Allocator

The `vcfixture` binary uses [mimalloc](https://github.com/microsoft/mimalloc)
as its global allocator. Bulk generation is allocation-dominated, and the swap
measures 1.63x–1.73x lower wall clock at every worker count.

This is the `mimalloc` feature, on by default. It pulls a C dependency built
via `cc`. To build without it:

```bash
cargo build --no-default-features --features cli
```

The library never installs a global allocator, whatever the feature state —
that choice belongs to the binary at the top of the dependency graph.
```

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml Cargo.lock src/bin examples/bulk_bench.rs docs README.md
git commit -m "perf(bulk): use mimalloc in the binaries and bench harness

Bulk generation is allocation-dominated (~47% of profile self time),
and PR #28 measured mimalloc at 1.63x-1.73x lower wall clock at every
worker count. Installed in the [[bin]] and example targets only; a
library must not impose a global allocator on its dependents.

On by default via the mimalloc feature; opt out with
--no-default-features.

Refs #26"
```

---

### Task 3: Per-thread scratch record (Part A)

Depends on Task 1. Owns `src/bulk/generate.rs` and the `map_init` block at `src/bulk/mod.rs:741-783`. Does **not** touch `write` or the temp helpers — those belong to Task 4.

**Files:**
- Modify: `src/bulk/generate.rs:287-388` (`SampleStats`, `to_record_buf`)
- Modify: `src/bulk/generate.rs:555-565` (the existing unit test that calls `to_record_buf`)
- Modify: `src/bulk/mod.rs:36` (import), `src/bulk/mod.rs:741-776` (`map_init` init closure and the `to_record_buf` call site)
- Modify: `src/bulk/mod.rs:1267` (the comment referencing `to_record_buf`'s key list)

**Interfaces:**
- Consumes: Task 1's golden snapshots.
- Produces:
  - `pub struct RecordScratch` in `src/bulk/generate.rs`
  - `pub fn RecordScratch::new(payload: &Payload) -> RecordScratch`
  - `pub fn RecordScratch::fill(&mut self, r: &GenRecord, phased: bool) -> &RecordBuf`
  - `pub fn to_record_buf` is **removed**; its two call sites move to `RecordScratch`.

- [ ] **Step 1: Write the failing tests**

Append to the `mod tests` block in `src/bulk/generate.rs`:

```rust
/// The scratch must produce a record indistinguishable from a fresh one,
/// *and* keep doing so after being reused. A scratch that leaked state
/// between records — a stale slot, a genotype vector that was appended to
/// instead of cleared — would pass a single-record test and corrupt every
/// record after the first.
#[test]
fn scratch_reuse_matches_a_fresh_scratch() {
    use crate::bulk::profile::Payload;

    let (p, s) = fixture();
    let mut rng = block_rng(7, 0, Stream::Content);

    let records: Vec<GenRecord> = (0..8)
        .map(|i| gen_record(&mut rng, &s, "chr1", 100 * (i + 1), 4, 2, &p.fitted))
        .collect();

    for payload in [
        Payload::GtOnly,
        Payload::GtVaf,
        Payload::Gatk,
        Payload::Mutect2,
    ] {
        let mut reused = RecordScratch::new(&payload);
        for (i, r) in records.iter().enumerate() {
            let phased = i % 2 == 0;
            let from_reused = reused.fill(r, phased).clone();
            let from_fresh = RecordScratch::new(&payload).fill(r, phased).clone();
            assert_eq!(
                from_reused, from_fresh,
                "record {i} of {payload:?} differs after scratch reuse"
            );
        }
    }
}

/// Ploidy varies between records in a real run only via the profile, but
/// the scratch must still shrink and grow its per-sample genotype buffers
/// correctly. A `clear()` that was actually a `truncate` to the wrong
/// length would show up here and nowhere else.
#[test]
fn scratch_handles_changing_sample_count() {
    use crate::bulk::profile::Payload;

    let (p, s) = fixture();
    let mut rng = block_rng(11, 0, Stream::Content);
    let mut scratch = RecordScratch::new(&Payload::GtOnly);

    for n_samples in [6usize, 2, 9, 1, 6] {
        let r = gen_record(&mut rng, &s, "chr1", 500, n_samples, 2, &p.fitted);
        let got = scratch.fill(&r, false);
        assert_eq!(
            got.samples().values().count(),
            n_samples,
            "scratch must resize to {n_samples} samples"
        );
        let fresh = RecordScratch::new(&Payload::GtOnly).fill(&r, false).clone();
        assert_eq!(*got, fresh, "resized scratch must match a fresh one");
    }
}
```

- [ ] **Step 2: Run to verify they fail**

```bash
cargo test --all-features --lib bulk::generate 2>&1 | tail -20
```

Expected: FAIL — `cannot find type RecordScratch in this scope`.

- [ ] **Step 3: Rewrite `SampleStats` to refill in place**

In `src/bulk/generate.rs`, add the imports:

```rust
use noodles_vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles_vcf::variant::record_buf::samples::sample::value::genotype::Allele;
use noodles_vcf::variant::record_buf::samples::sample::value::{Array, Genotype};
use noodles_vcf::variant::record_buf::samples::{Keys, Samples};
```

Replace `SampleStats`' `gt: String` field and its `new` / `value_for` methods with:

```rust
/// Per-sample FORMAT values derived cheaply and deterministically from one
/// sample's allele calls. Realism of the non-GT values is out of scope —
/// only their presence, type, and cardinality affect the benchmark (per the
/// design spec); do not add fields beyond what each [`Payload`] preset asks
/// for.
///
/// No longer carries a formatted `GT` string: the genotype is written
/// straight into the destination slot as a structured [`Genotype`] by
/// [`SampleStats::refill`], which lets the slot's backing `Vec<Allele>` be
/// reused across records and removes both this crate's integer formatting
/// and noodles' string reparse from the encode path (issue #26).
struct SampleStats {
    dp: i32,
    ad: [i32; 2],
    vaf: f32,
}

impl SampleStats {
    fn new(alleles: &[i8]) -> SampleStats {
        let n_ref = alleles.iter().filter(|&&a| a == 0).count() as i32;
        let n_alt = alleles.iter().filter(|&&a| a == 1).count() as i32;
        let dp = alleles.iter().filter(|&&a| a != -1).count() as i32;
        let vaf = if dp > 0 {
            n_alt as f32 / dp as f32
        } else {
            0.0
        };

        SampleStats {
            dp,
            ad: [n_ref, n_alt],
            vaf,
        }
    }

    /// Writes this sample's value for `key` into `slot`, reusing `slot`'s
    /// existing heap buffer when the variant already matches. Every key is
    /// written unconditionally on every record, so no slot can carry a
    /// stale value forward.
    fn refill(&self, key: &str, alleles: &[i8], phased: bool, slot: &mut Option<Value>) {
        match key {
            "GT" => refill_genotype(slot, alleles, phased),
            "DP" => *slot = Some(Value::Integer(self.dp)),
            "GQ" => *slot = Some(Value::Integer(99)),
            "VAF" | "AF" => *slot = Some(Value::Float(self.vaf)),
            "AD" => refill_int_array(slot, &[self.ad[0], self.ad[1]]),
            "PL" => refill_int_array(slot, &[0, 30, 60]),
            "F1R2" => refill_int_array(slot, &[self.ad[0] / 2, self.ad[1] / 2]),
            "F2R1" => refill_int_array(
                slot,
                &[self.ad[0] - self.ad[0] / 2, self.ad[1] - self.ad[1] / 2],
            ),
            "SB" => refill_int_array(slot, &[0, 0, 0, 0]),
            other => unreachable!("unhandled FORMAT key {other}: not in any Payload preset"),
        }
    }
}

/// Rewrites `slot` as a genotype over `alleles`, reusing the existing
/// `Vec<Allele>` when there is one.
///
/// # Phase bit
///
/// Position 0 is always [`Phasing::Unphased`], regardless of `phased`.
/// This is not cosmetic: it is what keeps output byte-identical to the
/// string path this replaces. `encode_genotype_str` initialises its
/// `last_phasing` to `"/"` and so never sets the phase bit on position 0,
/// whatever separator follows. Note that noodles' genotype *parser* is not
/// symmetric with its encoder here — it maps the single-allele string
/// `"0"` to `Phasing::Phased` — so the parser's behaviour must not be used
/// to infer what the encoder expects.
fn refill_genotype(slot: &mut Option<Value>, alleles: &[i8], phased: bool) {
    if !matches!(slot, Some(Value::Genotype(_))) {
        *slot = Some(Value::Genotype(Genotype::default()));
    }
    let Some(Value::Genotype(g)) = slot else {
        unreachable!("slot was just set to a Genotype")
    };

    let v = g.as_mut();
    v.clear();
    for (i, &a) in alleles.iter().enumerate() {
        let position = if a < 0 { None } else { Some(a as usize) };
        let phasing = if i > 0 && phased {
            Phasing::Phased
        } else {
            Phasing::Unphased
        };
        v.push(Allele::new(position, phasing));
    }
}

/// Rewrites `slot` as an integer array, reusing the existing `Vec` when
/// there is one.
fn refill_int_array(slot: &mut Option<Value>, values: &[i32]) {
    if !matches!(slot, Some(Value::Array(Array::Integer(_)))) {
        *slot = Some(Value::Array(Array::Integer(Vec::new())));
    }
    let Some(Value::Array(Array::Integer(v))) = slot else {
        unreachable!("slot was just set to an integer array")
    };

    v.clear();
    v.extend(values.iter().copied().map(Some));
}
```

- [ ] **Step 4: Replace `to_record_buf` with `RecordScratch`**

Delete `pub fn to_record_buf` entirely and add in its place:

```rust
/// The FORMAT keys each [`Payload`] preset renders, in order.
///
/// Must stay in sync with the `##FORMAT` header lines `BulkSpec::build_header`
/// emits (`src/bulk/mod.rs`) — a key rendered here with no header line
/// produces an unreadable file.
fn payload_keys(payload: &Payload) -> &'static [&'static str] {
    match payload {
        Payload::GtOnly => &["GT"],
        Payload::GtVaf => &["GT", "VAF"],
        Payload::Gatk => &["GT", "AD", "DP", "GQ", "PL"],
        Payload::Mutect2 => &["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
    }
}

/// A reusable [`RecordBuf`] plus the machinery to refill it from a
/// [`GenRecord`] without reallocating its per-sample buffers.
///
/// # Why this exists
///
/// The obvious shape — build and return a fresh `RecordBuf` per record —
/// costs four heap allocations per sample per record for `Payload::GtOnly`:
/// a formatted `GT` string, a clone of it, the per-sample
/// `Vec<Option<Value>>`, and noodles' own `Vec<i8>` while reparsing the
/// string. At the reference benchmark workload (2000 samples x 20000
/// records) that is roughly 160 million allocations, and allocator work
/// measured ~47% of profile self time (issue #26).
///
/// Holding one `RecordBuf` across records and refilling it in place removes
/// three of those four. The fourth is inside noodles and stays. Record-level
/// fields (`chrom`, `ref_`, `alts`) are still cloned per record — that is
/// about three allocations per *record* against ~160M on the sample path,
/// so reusing them would add code for no measurable gain.
///
/// # Reuse safety
///
/// Every FORMAT key of every sample is written on every [`RecordScratch::fill`]
/// call, and the per-sample slot vector is resized to the record's sample
/// count, so no value can survive from a previous record. The
/// `scratch_reuse_matches_a_fresh_scratch` test pins this against a fresh
/// scratch for all four payload presets.
pub struct RecordScratch {
    buf: RecordBuf,
    key_names: &'static [&'static str],
}

impl RecordScratch {
    /// Builds a scratch record for `payload`. The [`Keys`] set is
    /// constructed once here and then moved in and out of the record's
    /// [`Samples`] on each fill, never rebuilt or cloned.
    pub fn new(payload: &Payload) -> RecordScratch {
        let key_names = payload_keys(payload);
        let keys: Keys = key_names.iter().map(|k| k.to_string()).collect();
        let mut buf = RecordBuf::default();
        *buf.samples_mut() = Samples::new(keys, Vec::new());
        RecordScratch { buf, key_names }
    }

    /// Refills this scratch from `r` and returns it, ready to encode.
    ///
    /// The returned reference borrows the scratch, so the caller must finish
    /// with it before the next `fill`. That is exactly the encode-then-next
    /// shape `BulkSpec::stream_contigs` already has.
    pub fn fill(&mut self, r: &GenRecord, phased: bool) -> &RecordBuf {
        let ploidy = r.ploidy as usize;
        // `GenRecord` is a flat `pub` struct, so nothing prevents a caller
        // from constructing one with `ploidy: 0` or a `gts.len()` that isn't
        // a multiple of `ploidy`. `checked_div`'s `unwrap_or(0)` below
        // silently turns that into a zero-sample (or truncated) record
        // rather than failing, so assert the invariant explicitly first --
        // in debug/test builds this fails fast instead of silently
        // mis-encoding.
        debug_assert!(
            ploidy > 0 && r.gts.len() % ploidy == 0,
            "ploidy must be > 0 and evenly divide gts.len() (ploidy={ploidy}, gts.len()={})",
            r.gts.len()
        );
        let n_samples = r.gts.len().checked_div(ploidy).unwrap_or(0);

        // Take the sample block apart so the outer `Vec`, every per-sample
        // `Vec<Option<Value>>`, and every `Value`'s own buffer keep their
        // capacity. `Samples` exposes no `values_mut`, so this `From` impl
        // is the only way to reach them.
        let samples = std::mem::take(self.buf.samples_mut());
        let (keys, mut values) = <(Keys, Vec<Vec<Option<Value>>>)>::from(samples);

        values.resize_with(n_samples, Vec::new);
        for (i, slots) in values.iter_mut().enumerate() {
            let alleles = &r.gts[i * ploidy..(i + 1) * ploidy];
            let stats = SampleStats::new(alleles);
            slots.resize_with(self.key_names.len(), || None);
            for (slot, &k) in slots.iter_mut().zip(self.key_names) {
                stats.refill(k, alleles, phased, slot);
            }
        }

        *self.buf.samples_mut() = Samples::new(keys, values);

        // Record-level fields: cloned, not reused. See the type's docs.
        self.buf.reference_sequence_name_mut().clear();
        self.buf.reference_sequence_name_mut().push_str(&r.chrom);
        *self.buf.variant_start_mut() =
            Some(Position::try_from(r.pos as usize).expect("pos must be >= 1"));
        self.buf.reference_bases_mut().clear();
        self.buf.reference_bases_mut().push_str(&r.ref_);
        let alts = self.buf.alternate_bases_mut().as_mut();
        alts.clear();
        alts.extend(r.alts.iter().cloned());

        &self.buf
    }
}
```

Remove the now-unused `AlternateBases`, `Position` builder, `Samples`, and `Keys` imports only if the compiler says they are unused — `Position` is still needed above.

- [ ] **Step 5: Update the existing unit test at `generate.rs:555-565`**

Replace the `to_record_buf` call. Bind the scratch to a local rather than
chaining off a temporary — the following lines read `buf.samples().keys()`,
and a temporary scratch would not outlive the borrow:

```rust
let mut scratch = RecordScratch::new(&payload);
let buf = scratch.fill(&r, true);
```

- [ ] **Step 6: Update the call site in `stream_contigs`**

In `src/bulk/mod.rs:36`, change the import:

```rust
use generate::{block_rng, gen_record, RecordScratch, Stream};
```

In the `map_init` at `mod.rs:741`, change the init closure to build a scratch alongside the encoder:

```rust
.map_init(
    || {
        let enc = BlockEncoder::new(self.format, header).expect(
            "encoding a header into a Vec<u8> cannot fail; the same \
             construction was validated before the loop",
        );
        // One scratch record per rayon worker, reused for every record
        // that worker encodes. See `RecordScratch` for why.
        (enc, RecordScratch::new(&self.payload))
    },
    |(enc, scratch), &b| {
```

and the record call at `mod.rs:775`:

```rust
let buf = scratch.fill(&g, phased);
enc.push(header, buf)?;
```

Note `buf` is now a `&RecordBuf`, so `enc.push(header, buf)` replaces `enc.push(header, &buf)`.

- [ ] **Step 7: Update the stale comment at `mod.rs:1267`**

It currently says "Must stay in sync with `generate::to_record_buf`'s own (private) `key_names`". Change `to_record_buf` to `payload_keys`.

- [ ] **Step 8: Run the new tests**

```bash
cargo test --all-features --lib bulk::generate 2>&1 | tail -20
```

Expected: PASS, including `scratch_reuse_matches_a_fresh_scratch` and `scratch_handles_changing_sample_count`.

- [ ] **Step 9: Run the golden gate — the real verdict**

```bash
cargo test --all-features --test bulk_golden 2>&1 | tail -20
```

Expected: **12 passed, 0 failed, and no `.snap.new` files created.**

```bash
ls tests/snapshots/*.snap.new 2>/dev/null && echo "REGRESSION" || echo "clean"
```

Expected: `clean`. If any snapshot changed, **stop**. Do not accept the new snapshot. A changed digest means Part A altered the output, which the spec forbids. The likeliest cause is the phase bit on allele position 0 — re-read `refill_genotype`'s doc comment. The second likeliest is the text VCF writer rendering `Value::Genotype` differently from `Value::String`, which would show as VCF and VcfGz failing while all four BCF goldens pass. Report that shape explicitly if you see it; it changes the design, not the test.

- [ ] **Step 10: Full suite and lints**

```bash
cargo test --all-features 2>&1 | tail -20
cargo fmt --check && cargo clippy --all-features -- -D warnings
```

Expected: 165 + 2 new = 167 passed, 0 failed. Lints clean.

- [ ] **Step 11: Commit**

```bash
git add src/bulk/generate.rs src/bulk/mod.rs
git commit -m "perf(bulk): reuse a per-thread scratch record across records

Holding one RecordBuf per rayon worker and refilling its sample buffers
in place takes steady-state allocations from four per sample per record
to one. GT becomes a structured Value::Genotype, which also drops this
crate's integer formatting and noodles' string reparse from the encode
path.

Output is byte-identical: the golden digests for all three formats and
all four payload presets are unchanged.

Refs #26"
```

---

### Task 4: Write through a temp file beside the destination (Part C)

Depends on Task 1. Owns `write` (`mod.rs:505-545`) and the temp helpers (`mod.rs:860-1015`). Does **not** touch the `map_init` block — that belongs to Task 3.

**Files:**
- Modify: `src/bulk/mod.rs:505-545` (`BulkSpec::write`)
- Modify: `src/bulk/mod.rs:860-920` (`resolve_target_counts`' temp handling)
- Modify: `src/bulk/mod.rs:940-995` (`write_to_temp`, `measured_bytes`)
- Modify: `tests/bulk.rs` (new regression test)

**Interfaces:**
- Consumes: Task 1's golden snapshots.
- Produces:
  - `fn BulkSpec::write_to_temp(&self, pool, samplers, fitted, per_contig_count, dir: &Path) -> Result<(NamedTempFile, u64, Summary), BulkError>` — gains a `dir` parameter.
  - `fn BulkSpec::measured_bytes(&self, pool, samplers, fitted, per_contig_count, dir: &Path) -> Result<u64, BulkError>` — gains a `dir` parameter.
  - `fn BulkSpec::resolve_target_counts(...)` — gains a `dir: &Path` parameter.
  - `fn temp_dir_for(dest: &Path) -> &Path` — free function.
  - `fn cleanup_csi(tmp_path: &Path, format: Format)` — free function, replaces two duplicated inline blocks.
  - `fn write_summary_json(dest: &Path, summary: &Summary) -> Result<(), BulkError>` — free function, replaces two duplicated inline blocks.

- [ ] **Step 1: Write the failing regression test**

Append to `tests/bulk.rs`:

```rust
/// The #27 regression guard. A mid-stream failure must leave nothing at
/// the caller's destination — not a truncated BCF, not a stray `.csi`,
/// not a `.summary.json`.
///
/// The failure is induced by making the destination's parent directory
/// read-only *after* the spec is built, so `BulkWriter::create` (or the
/// temp creation that replaces it) fails. This is deliberately a failure
/// the old code would also have caught; what it pins is the *absence of
/// debris*, which the old code did not guarantee.
#[test]
fn failed_write_leaves_nothing_at_the_destination() {
    let dir = tempfile::tempdir().unwrap();
    let out_dir = dir.path().join("out");
    std::fs::create_dir(&out_dir).unwrap();
    let path = out_dir.join("a.bcf");

    // Make the output directory unwritable so generation cannot complete.
    let mut perms = std::fs::metadata(&out_dir).unwrap().permissions();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        perms.set_mode(0o500);
    }
    std::fs::set_permissions(&out_dir, perms).unwrap();

    let result = spec().size(Size::RecordsPerContig(600)).write(&path);

    // Restore permissions so the assertions below (and tempdir cleanup)
    // can read the directory.
    let mut perms = std::fs::metadata(&out_dir).unwrap().permissions();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        perms.set_mode(0o700);
    }
    std::fs::set_permissions(&out_dir, perms).unwrap();

    assert!(result.is_err(), "write into an unwritable dir must fail");
    assert!(!path.exists(), "no output file may remain at the destination");
    assert!(
        !path.with_extension("bcf.csi").exists(),
        "no .csi may remain at the destination"
    );
    let leftovers: Vec<_> = std::fs::read_dir(&out_dir)
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .collect();
    assert!(
        leftovers.is_empty(),
        "destination directory must be empty, found {leftovers:?}"
    );
}

/// Temps must be created beside the destination, not in `TMPDIR`.
///
/// `TMPDIR` is routinely a `tmpfs` or a different filesystem from the
/// output. If promotion falls back to `std::fs::copy`, bulk-scale writes
/// double their I/O and can fail outright on a `tmpfs` too small to hold
/// the output. Pointing `TMPDIR` at a path that does not exist proves the
/// write never goes there.
#[test]
fn temp_files_are_created_beside_the_destination() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");

    let bogus = dir.path().join("this-tmpdir-does-not-exist");
    temp_env::with_var("TMPDIR", Some(&bogus), || {
        spec()
            .size(Size::RecordsPerContig(600))
            .write(&path)
            .expect("write must not depend on TMPDIR");
    });

    assert!(path.exists(), "output must land at the destination");
}
```

**Note:** if `temp_env` is not already a dev-dependency, do not add one. Replace the `temp_env::with_var` block with a direct `std::env::set_var`/`remove_var` pair guarded by a `#[serial]`-style comment, or simpler — drop this second test and instead assert the property directly in a unit test inside `src/bulk/mod.rs` that calls `temp_dir_for` and checks it returns the destination's parent. Prefer the unit test; it is deterministic and needs no environment mutation.

- [ ] **Step 2: Run to verify they fail**

```bash
cargo test --all-features --test bulk failed_write_leaves_nothing 2>&1 | tail -20
```

Expected: FAIL — a truncated `a.bcf` remains in the destination directory.

- [ ] **Step 3: Add the three extracted helpers**

Add near the other free functions at the bottom of `src/bulk/mod.rs`:

```rust
/// The directory temp files are created in for a given destination.
///
/// Deliberately not `TMPDIR`: `NamedTempFile::new` puts temps there, which
/// is routinely a `tmpfs` or a different filesystem from the output. A
/// cross-filesystem `persist` falls back to `std::fs::copy`, which for
/// bulk-scale output doubles the I/O and can fail outright when `TMPDIR`
/// is too small to hold a second copy. Creating the temp beside the
/// destination makes promotion a same-filesystem rename in every case.
fn temp_dir_for(dest: &Path) -> &Path {
    match dest.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        // A bare filename like `out.bcf` has a parent of `""`.
        _ => Path::new("."),
    }
}

/// Best-effort removal of a temp file's `.csi` companion, which
/// `NamedTempFile`'s `Drop` does not know about.
fn cleanup_csi(tmp_path: &Path, format: Format) {
    if matches!(format, Format::Bcf) {
        let mut csi_path = tmp_path.as_os_str().to_os_string();
        csi_path.push(".csi");
        let _ = std::fs::remove_file(csi_path);
    }
}

/// Writes `<dest>.summary.json`.
fn write_summary_json(dest: &Path, summary: &Summary) -> Result<(), BulkError> {
    let json = summary.to_json()?;
    let mut summary_path = dest.as_os_str().to_os_string();
    summary_path.push(".summary.json");
    std::fs::write(&summary_path, json)?;
    Ok(())
}
```

- [ ] **Step 4: Thread `dir` through the temp helpers**

In `write_to_temp`, add a `dir: &Path` parameter and change the temp construction:

```rust
let tmp = tempfile::NamedTempFile::new_in(dir)?;
```

In `measured_bytes`, add a `dir: &Path` parameter, forward it to `write_to_temp`, and replace the inline `.csi` block with `cleanup_csi(tmp.path(), self.format);`.

In `resolve_target_counts`, add a `dir: &Path` parameter, forward it to every `write_to_temp` / `measured_bytes` call, and replace its inline `.csi` block (`mod.rs:913-917`) with `cleanup_csi(tmp.path(), self.format);`.

- [ ] **Step 5: Route every `Size` variant through temp-then-promote**

Replace the body of `BulkSpec::write` from the `counts` match through the end with:

```rust
let tmp_dir = temp_dir_for(path);

// Every `Size` variant generates into a temp file beside the
// destination and is promoted by a single rename once the output is
// complete and indexed. Writing straight to `path` would leave a
// truncated, un-indexed file there if encoding or I/O failed
// mid-stream (issue #27); with this shape a failure leaves the
// destination untouched.
let counts: Vec<u64> = match &self.size {
    Size::RecordsPerContig(n) => vec![*n; self.contig_ids.len()],
    Size::Records(total) => distribute_by_n_variants(fitted, &self.contig_ids, *total),
    Size::PerContig(map) => per_contig_counts(map, &self.contig_ids)?,
    Size::Target(target_bytes) => {
        let (_counts, tmp, _bytes, summary) =
            self.resolve_target_counts(&pool, &samplers, fitted, *target_bytes, tmp_dir)?;
        Self::promote_temp(tmp, path, self.format)?;
        write_summary_json(path, &summary)?;
        return Ok(summary);
    }
};

let (tmp, _bytes, summary) =
    self.write_to_temp(&pool, &samplers, fitted, &counts, tmp_dir)?;
Self::promote_temp(tmp, path, self.format)?;
write_summary_json(path, &summary)?;

Ok(summary)
```

The `compute_layouts` / `build_header` / `BulkWriter::create` / `stream_contigs` / `finish_and_index` sequence that used to live inline here is deleted — `write_to_temp` already performs exactly it, and its doc comment already states it is byte-exact against this path.

- [ ] **Step 6: Run the regression test**

```bash
cargo test --all-features --test bulk failed_write_leaves_nothing 2>&1 | tail -20
```

Expected: PASS.

- [ ] **Step 7: Run the golden gate**

```bash
cargo test --all-features --test bulk_golden 2>&1 | tail -20
ls tests/snapshots/*.snap.new 2>/dev/null && echo "REGRESSION" || echo "clean"
```

Expected: 12 passed, `clean`. Byte-for-byte identical output is the whole claim of `write_to_temp` being byte-exact; if a snapshot moves, that claim was false and must be investigated rather than papered over.

- [ ] **Step 8: Full suite and lints**

```bash
cargo test --all-features 2>&1 | tail -20
cargo fmt --check && cargo clippy --all-features -- -D warnings
```

Expected: all green. `target_size_is_byte_identical_across_runs` and `same_seed_gives_byte_identical_output_across_thread_counts` must still pass untouched.

- [ ] **Step 9: Commit**

```bash
git add src/bulk/mod.rs tests/bulk.rs
git commit -m "fix(bulk): never leave a partial file at the destination

Every Size variant now generates into a temp file beside the
destination and is promoted by a single rename once the output is
complete and indexed, so a mid-stream encode or write failure leaves
the destination untouched. Size::Target already worked this way; this
extends it to the rest and deletes the duplicated
create/stream/finish sequence.

Temps are created in the destination's directory rather than TMPDIR so
promotion is always a same-filesystem rename, never a copy of a
bulk-scale file.

Closes #27"
```

---

### Task 5: Measure and publish

Depends on Tasks 2, 3, and 4. Sequential — benchmarks must not share the machine.

**Files:**
- Modify: `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md` (append a section)
- Comment on issues #26 and #27

**Interfaces:**
- Consumes: all prior tasks.
- Produces: measured deltas; no code.

- [ ] **Step 1: Confirm the machine is quiet**

```bash
cat /proc/self/status | rg Cpus_allowed_list
pgrep -a bulk_bench || echo "no bench running"
uptime
```

Expected: `Cpus_allowed_list: 0-3,48-51` (4 physical cores + SMT siblings), no bench running, load average near zero. **Do not proceed if anything else is running.**

- [ ] **Step 2: Build the three binaries to compare**

Three configurations, each a separate build:

```bash
# glibc, all source changes (Parts A + C)
cargo build --release --example bulk_bench --features bulk --no-default-features
cp target/release/examples/bulk_bench "$CLAUDE_JOB_DIR/tmp/bench-glibc"

# mimalloc, all source changes (Parts A + B + C)
cargo build --release --example bulk_bench --features bulk
cp target/release/examples/bulk_bench "$CLAUDE_JOB_DIR/tmp/bench-mimalloc"
```

For the baseline, build from the branch point:

```bash
git stash push -u -m "task5-bench-wip"
git switch --detach 501ce4d
cargo build --release --example bulk_bench --features bulk
cp target/release/examples/bulk_bench "$CLAUDE_JOB_DIR/tmp/bench-baseline"
git switch -
git stash list --format='%H %gs' | rg task5-bench-wip
```

Restore with `git stash apply <sha>` (never bare `git stash pop` — the stash stack is shared with other worktrees and sessions), then drop the entry by re-finding its `stash@{n}` by tag.

Verify no stale instrumentation survived into any binary:

```bash
for b in baseline glibc mimalloc; do
  printf '%s: ' "$b"
  strings "$CLAUDE_JOB_DIR/tmp/bench-$b" | rg -c '\[instr\]' || echo 0
done
```

Expected: `0` for all three.

- [ ] **Step 3: Sweep, round-robin across worker counts**

Reps of a single cell must **not** run back to back. Running them consecutively lets one noise burst contaminate all reps identically, which defeats the min-of-N estimator — the exact failure that inflated PR #28's anchor by 8.6%. Round-robin instead:

```bash
cd "$CLAUDE_JOB_DIR/tmp"
for pass in 1 2 3 4 5; do
  for w in 1 2 3 4 6 8; do
    for b in baseline glibc mimalloc; do
      VCFIXTURE_BENCH_REPS=1 VCFIXTURE_BENCH_FORMAT=bcf \
        ./bench-$b $w >> "sweep-$b.txt" 2>&1
    done
  done
done
```

Run this with the harness's background mode. Never with `nohup`, `setsid`, `disown`, or a trailing `&`.

- [ ] **Step 4: Re-check the `workers=1` anchor**

Before dividing anything by `min_s(1)`, take five further time-separated single-shot passes at `w=1` for each binary and confirm the minimum is stable. If any anchor is more than ~5% above the minimum of those passes, adopt the lower value and say so explicitly — this is the check whose absence produced PR #28's inflated figure.

- [ ] **Step 5: Compute and report**

For each binary compute `S(w) = min_s(1) / min_s(w)` and efficiency `S(w) / min(w, 4)` — **against 4 physical cores, not 8 logical CPUs.**

Report three separate deltas against the baseline, never a single joint number:

- **Part A alone** (`bench-glibc` vs `bench-baseline`): the allocation-count reduction.
- **Part B alone**: already known at 1.63x–1.73x from PR #28. Confirm it reproduces.
- **A + B combined** (`bench-mimalloc` vs `bench-baseline`).

If the combined figure does not exceed Part B alone, **Part A did not pay for itself.** Report that plainly rather than absorbing it into a joint number, and recommend reverting Task 3 — the spec commits to judging it on its own.

Also report peak RSS for each configuration. Scratch reuse holds buffers alive for a thread's lifetime rather than freeing them per record, so a modest RSS increase is expected and worth stating; a large one is a finding.

- [ ] **Step 6: Append to the measurements doc**

Add a `## Allocation reduction (#26)` section to `docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md` with the exact commands, the raw per-pass readings (not just the minima), the anchor re-check from Step 4, and the three deltas. State the sample count behind every figure.

Note: sections in that document are ordered by commit, not by investigation order. Do not silently reorder existing ones.

- [ ] **Step 7: Comment on the issues**

On #26: the measured effect of each part, and whether the allocation-count reduction paid for itself independently of the allocator swap. Do not claim the allocator question is closed — PR #28 established that a residual sub-linear scaling deficit remains genuinely unexplained, with memory-bandwidth saturation, inherent per-thread allocation overhead, the separate `compute_layouts` pass, and SMT siblings contributing less than a full core all still **untested**. This task changes the absolute cost, not that open question.

On #27: state that the fix is behavioural and pinned by `failed_write_leaves_nothing_at_the_destination`.

- [ ] **Step 8: Verify nothing is still running**

```bash
pgrep -a bulk_bench && echo "STILL RUNNING - kill it" || echo "clean"
```

An unreaped child at session teardown blocks on NFS in uninterruptible D state and drains the Slurm node. This check is not optional.

- [ ] **Step 9: Commit and open the stacked PR**

```bash
git add docs/superpowers/plans/2026-08-06-bulk-perf-measurements.md
git commit -m "docs(bulk): measured effect of allocation reduction and mimalloc

Refs #26"
git push -u origin fix-26-27-alloc-partial-output
gh pr create --base worktree-fix-22-bulk-perf \
  --title "perf(bulk): cut per-sample allocations; write atomically to the destination" \
  --body "..."
```

**`--base worktree-fix-22-bulk-perf`, not `main`.** This PR stacks on #28 and must not be retargeted until #28 merges.

---

## Self-Review

**Spec coverage.** Part A → Task 3. Part B → Task 2. Part C → Task 4. Byte-equality gate → Task 1 (all 3 formats x 4 payloads, as the spec requires). #27 regression test → Task 4 Step 1. Cross-filesystem promotion → Task 4 Steps 3-4 plus the `temp_dir_for` unit test. Measurement methodology (4 physical cores, round-robin, anchor re-check, benchmarks alone) → Task 5 Steps 1, 3, 4, 5, 8. Separate reporting of Parts A and B → Task 5 Step 5. Every risk in the spec's risk table maps to a step.

**Placeholder scan.** No TBD/TODO. Every code step carries real code. The one `--body "..."` in Task 5 Step 9 is a PR description written at execution time from the measured results, which cannot be known now; the surrounding steps specify what must be in it.

**Type consistency.** `RecordScratch::new(&Payload) -> RecordScratch` and `RecordScratch::fill(&mut self, &GenRecord, bool) -> &RecordBuf` are used identically in Tasks 3 Steps 4, 5, 6 and in both new unit tests. `payload_keys` is the name used in both the function definition and the `mod.rs:1267` comment fix. `temp_dir_for`, `cleanup_csi`, and `write_summary_json` are declared in Task 4 Step 3 and used with matching signatures in Steps 4 and 5. `write_to_temp`, `measured_bytes`, and `resolve_target_counts` all gain the same trailing `dir: &Path` parameter.

**Known soft spot.** Task 4 Step 1's second test depends on `temp_env`, which is not a dev-dependency; the step explicitly directs the implementer to the unit-test alternative instead of adding one. Flagged rather than left to discovery.
