# BulkError::Invalid Split — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the catch-all `BulkError::Invalid(String)` with per-class
variants, so the `"invalid profile: "` prefix stops misdescribing spec,
argument-parsing, runtime, and profile-loading failures.

**Architecture:** Additive-then-subtractive. Task 1 *adds* the new variants
while leaving `Invalid` in place, so the tree compiles and tests pass after
every task. Tasks 2–4 re-route call sites file-group by file-group. Task 5
deletes the now-unused `Invalid` and runs the end-to-end gates. This ordering
is what makes each task independently testable and committable — a "define the
new enum and delete the old variant in one step" ordering would leave the crate
uncompilable until every call site had been touched.

**Tech Stack:** Rust 2021 (rust-version 1.86), `thiserror` 2, `clap` 4.6,
`rayon` 1, `cargo test` / `clippy` / `fmt`, `prek` git hooks, `commitizen`
(`cz`) for versioning.

## Global Constraints

- **Spec:** `docs/superpowers/specs/2026-08-05-bulk-error-split-design.md`.
  Read it before starting.
- **Feature flags:** `BulkError` lives behind the `bulk` feature; the CLI
  behind `cli` (which implies `bulk`). Always build and test with
  `--all-features`, or nothing in scope compiles.
- **Error-message text is preserved verbatim** except where prefix removal
  requires a change. Do not reword messages while moving them.
- **`PayloadPloidy` keeps the `invalid profile:` prefix** — it *is* a profile
  error. Every other new variant drops the prefix.
- **Structure where a caller can branch; string where a caller can only
  display.** Do not invent extra variants for the 13 profile-content messages
  that stay behind `InvalidProfile(String)`.
- **No behaviour changes.** This plan re-routes existing errors. It adds and
  removes no validation checks. If you find a check you think is wrong, leave
  it and report it — do not fix it here.
- **`CHANGELOG.md` is generated**, not hand-edited: `.cz.toml` sets
  `update_changelog_on_bump = true`, so `cz bump` writes it from commit
  messages. The breaking change is declared via the commit message in Task 5
  (`refactor(bulk)!:` subject + `BREAKING CHANGE:` footer), not by editing the
  file. Do not touch `CHANGELOG.md`.
- **Commit hooks:** `prek` hooks are already installed in this worktree and run
  `cargo fmt`, `cargo clippy`, `cargo check`, and a commitizen message check on
  every commit. A commit that fails them is a task that is not done.
- **Gates for every task:** `cargo test --all-features` passes, and
  `cargo clippy --all-features --all-targets -- -D warnings` is clean.

## Parallelization

Tasks 2, 3, and 4 are mutually independent — they touch disjoint `src/` files
and each is a complete red→green→commit cycle on its own. They all depend on
Task 1; Task 5 depends on all of them.

```
Task 1  ──┬──  Task 2 (profile.rs, sample.rs)  ──┬──  Task 5
          ├──  Task 3 (mod.rs)                  ──┤
          └──  Task 4 (writer.rs, bin)          ──┘
```

**Recommendation: run them sequentially anyway.** The whole change is ~25
mechanical edits behind one enum; dispatching three parallel agents in separate
worktrees costs more in merge handling than it saves, and Tasks 2 and 3 both
touch `tests/bulk.rs` (at lines ~481 and ~457 respectively — separate hunks
that git will usually auto-merge, but "usually" is the problem). If you do
parallelize, use `superpowers:dispatching-parallel-agents` with
`superpowers:subagent-driven-development`, give each agent its own worktree,
and merge Task 3 last since it has the largest `tests/bulk.rs` footprint.

## File Structure

| File | Role in this change |
|---|---|
| `src/bulk/mod.rs` | Owns `BulkError`. Task 1 (enum), Task 3 (6 call sites + `parse_size`), Task 5 (delete `Invalid`) |
| `src/bulk/profile.rs` | 13 profile-content sites → `InvalidProfile`; 1 site → `PayloadPloidy` |
| `src/bulk/sample.rs` | 1 profile-content site (`titv`) → `InvalidProfile` |
| `src/bulk/writer.rs` | 1 site → `CompressionLevel` |
| `src/bin/vcfixture.rs` | 1 site → `ProfileLoad`; 1 site deleted in favour of a clap value parser |
| `tests/bulk.rs` | Two over-broad assertions tightened; new spec-variant tests |
| `tests/cli.rs` | `parse_size` assertion tightened; `--threads 0` usage-error test |

---

### Task 1: Add the new variants

Adds every new variant to `BulkError` while leaving `Invalid` in place. Nothing
is routed yet, so the crate still compiles and every existing test still
passes. The deliverable is the enum itself plus `Display` tests pinning each
message.

**Files:**
- Modify: `src/bulk/mod.rs:44-55` (the `BulkError` enum)
- Test: `src/bulk/mod.rs` (the existing `#[cfg(test)] mod tests` at line ~939)

**Interfaces:**
- Consumes: nothing.
- Produces: the variants every later task routes to —
  - `BulkError::InvalidProfile(String)`
  - `BulkError::PayloadPloidy { payload: Payload, ploidy: u8 }`
  - `BulkError::NoContigs`
  - `BulkError::NoSamples`
  - `BulkError::DuplicateContig(String)`
  - `BulkError::BadSize(String)`
  - `BulkError::CompressionLevel(String)`
  - `BulkError::ProfileLoad { path: String, source: std::io::Error }`
  - `BulkError::WorkerPool(rayon::ThreadPoolBuildError)` (has `#[from]`)
  - `BulkError::TargetNotReached { target_bytes: u64, corrections: usize }`

- [ ] **Step 1: Write the failing test**

Append to the existing `mod tests` block at the bottom of `src/bulk/mod.rs`
(inside it, after the last existing test):

```rust
    /// The whole point of the split: only profile-content errors carry the
    /// `invalid profile:` prefix. Every other class names its own problem.
    #[test]
    fn error_messages_are_classified() {
        // Profile content keeps the prefix -- these messages name a profile
        // field, not themselves, so the prefix is load-bearing.
        assert_eq!(
            BulkError::InvalidProfile("ploidy must be >= 1".into()).to_string(),
            "invalid profile: ploidy must be >= 1"
        );
        assert!(BulkError::PayloadPloidy {
            payload: Payload::Gatk,
            ploidy: 3,
        }
        .to_string()
        .starts_with("invalid profile: payload Gatk emits AD and/or PL"));

        // Everything else must NOT claim the profile is at fault.
        for e in [
            BulkError::NoContigs,
            BulkError::NoSamples,
            BulkError::DuplicateContig("chr1".into()),
            BulkError::BadSize("banana".into()),
            BulkError::CompressionLevel("level 99 out of range".into()),
            BulkError::TargetNotReached {
                target_bytes: 1024,
                corrections: 4,
            },
        ] {
            let msg = e.to_string();
            assert!(
                !msg.starts_with("invalid profile:"),
                "{e:?} must not be described as an invalid profile, got: {msg}"
            );
        }

        assert_eq!(BulkError::NoContigs.to_string(), "need >= 1 output contig");
        assert_eq!(BulkError::NoSamples.to_string(), "need >= 1 sample");
        assert!(BulkError::DuplicateContig("chr1".into())
            .to_string()
            .starts_with(r#"duplicate output contig name: "chr1""#));
        assert!(BulkError::BadSize("banana".into())
            .to_string()
            .starts_with(r#"bad size: "banana""#));
        assert_eq!(
            BulkError::TargetNotReached {
                target_bytes: 1024,
                corrections: 4,
            }
            .to_string(),
            "could not reach target size 1024 bytes within 4 corrective rounds"
        );
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib bulk::tests::error_messages_are_classified`
Expected: FAIL — compile error, `no variant named InvalidProfile found for enum BulkError` (and similar for each new variant).

- [ ] **Step 3: Write minimal implementation**

Replace the `BulkError` enum at `src/bulk/mod.rs:44-55` in full:

```rust
/// Errors from bulk generation.
///
/// Variants are grouped by *who* is at fault, because the message a user
/// sees depends on it: a malformed profile JSON, a spec the caller built
/// wrong, an unparseable argument, or a failure at generation time. The
/// `invalid profile:` prefix appears on exactly the first group -- those
/// messages name a profile field rather than describing themselves, so they
/// need the context; every other variant's message stands alone.
///
/// Profile-content failures share one `InvalidProfile(String)` rather than
/// getting a variant each: a caller cannot act differently on "histogram
/// edges must be increasing" than on "histogram weights must sum > 0", so
/// the extra variants would only ever reach `Display`.
#[derive(Debug, thiserror::Error)]
pub enum BulkError {
    #[error("unknown builtin profile: {0}")]
    UnknownProfile(String),

    // --- profile content --------------------------------------------------
    /// A fitted or dialed statistic that fails [`Profile::validate`].
    #[error("invalid profile: {0}")]
    InvalidProfile(String),
    /// The one profile-content failure a caller can fix without editing the
    /// profile JSON -- by choosing a payload that doesn't emit AD/PL.
    #[error(
        "invalid profile: payload {payload:?} emits AD and/or PL, which are \
         hard-coded for diploid (ploidy 2) calls, but ploidy is {ploidy}"
    )]
    PayloadPloidy { payload: Payload, ploidy: u8 },

    // --- spec / caller validation -----------------------------------------
    #[error("need >= 1 output contig")]
    NoContigs,
    #[error("need >= 1 sample")]
    NoSamples,
    #[error(
        "duplicate output contig name: {0:?} (each requested contig must be \
         unique; duplicates produce backwards positions and a CSI that \
         silently drops region-query hits)"
    )]
    DuplicateContig(String),

    // --- argument parsing -------------------------------------------------
    #[error("bad size: {0:?} (expected a byte count, optionally suffixed KB/MB/GB)")]
    BadSize(String),
    #[error("invalid compression level: {0}")]
    CompressionLevel(String),

    // --- profile loading --------------------------------------------------
    #[error("profile {path:?} is not a builtin name and could not be read as a file: {source}")]
    ProfileLoad {
        path: String,
        source: std::io::Error,
    },

    // --- runtime ----------------------------------------------------------
    #[error("failed to build worker pool: {0}")]
    WorkerPool(#[from] rayon::ThreadPoolBuildError),
    #[error("could not reach target size {target_bytes} bytes within {corrections} corrective rounds")]
    TargetNotReached { target_bytes: u64, corrections: usize },

    /// Deprecated catch-all, removed in Task 5 once every call site is
    /// routed. Do not add new uses.
    #[error("invalid profile: {0}")]
    Invalid(String),

    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}
```

Note `ProfileLoad` deliberately uses a plain `source:` field with **no**
`#[from]` and **no** `#[source]` attribute conflict: `BulkError` already has
`Io(#[from] std::io::Error)`, so a second `#[from] std::io::Error` would be a
duplicate-`From`-impl compile error. Naming the field `source` is enough for
`thiserror` to wire up `Error::source()`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib bulk::tests::error_messages_are_classified`
Expected: PASS

Then confirm nothing else regressed:

Run: `cargo test --all-features`
Expected: PASS (all existing tests still green — no call site has moved yet)

- [ ] **Step 5: Commit**

```bash
git add src/bulk/mod.rs
git commit -m "refactor(bulk): add per-class BulkError variants (#16)

Adds the variants that replace the BulkError::Invalid catch-all, leaving
Invalid in place so call sites can be re-routed incrementally."
```

---

### Task 2: Route profile-content errors

**Files:**
- Modify: `src/bulk/profile.rs` (14 `BulkError::Invalid` sites at lines ~130,
  133, 145, 149, 160, 163, 173, 178, 183, 186, 189, 213, 221)
- Modify: `src/bulk/sample.rs:103`
- Test: `tests/bulk.rs:481` (tighten existing assertion)

**Interfaces:**
- Consumes: `BulkError::InvalidProfile`, `BulkError::PayloadPloidy` from Task 1.
- Produces: nothing new. `Profile::validate`, `Histogram::validate`,
  `ClassMix::validate`, and `Samplers::new` keep their existing signatures
  (`-> Result<_, BulkError>`).

- [ ] **Step 1: Write the failing test**

In `tests/bulk.rs`, in `non_diploid_profile_rejects_payloads_declaring_pl_or_ad`
(around line 481), replace the assertion:

```rust
        assert!(
            matches!(result, Err(BulkError::Invalid(_))),
            "payload {payload:?} declares PL/AD, which are diploid-only; \
             ploidy 3 must be rejected, got: {result:?}"
        );
```

with:

```rust
        assert!(
            matches!(
                result,
                Err(BulkError::PayloadPloidy {
                    ploidy: 3,
                    payload: ref p,
                }) if *p == payload
            ),
            "payload {payload:?} declares PL/AD, which are diploid-only; \
             ploidy 3 must be rejected with PayloadPloidy naming the payload \
             and the offending ploidy, got: {result:?}"
        );
```

This is the point of the change: the old assertion passed if `validate()`
failed for *any* reason, so it never actually proved the PL/AD guard ran.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --test bulk non_diploid_profile_rejects_payloads_declaring_pl_or_ad`
Expected: FAIL — the assertion trips, because `validate()` still returns
`BulkError::Invalid(_)`.

- [ ] **Step 3: Write minimal implementation**

In `src/bulk/profile.rs`, replace the payload/ploidy guard (lines ~132-138):

```rust
        if self.dialed.payload.needs_diploid() && self.dialed.ploidy != 2 {
            return Err(BulkError::Invalid(format!(
                "payload {:?} emits AD and/or PL, which are hard-coded for \
                 diploid (ploidy 2) calls, but ploidy is {}",
                self.dialed.payload, self.dialed.ploidy
            )));
        }
```

with:

```rust
        if self.dialed.payload.needs_diploid() && self.dialed.ploidy != 2 {
            return Err(BulkError::PayloadPloidy {
                payload: self.dialed.payload.clone(),
                ploidy: self.dialed.ploidy,
            });
        }
```

(`Payload` derives `Clone` but not `Copy`, hence the `.clone()`.)

Then rename every *remaining* `BulkError::Invalid` in `src/bulk/profile.rs` and
`src/bulk/sample.rs` to `BulkError::InvalidProfile`, changing nothing else —
same message strings, same `.into()` / `format!` calls. There are 13 in
`profile.rs` and 1 in `sample.rs` (`"titv must be > 0"`).

```bash
sed -i 's/BulkError::Invalid(/BulkError::InvalidProfile(/g' src/bulk/profile.rs src/bulk/sample.rs
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --all-features`
Expected: PASS

Run: `cargo clippy --all-features --all-targets -- -D warnings`
Expected: clean

Confirm the routing is complete for these two files:

Run: `rg -n "BulkError::Invalid\(" src/bulk/profile.rs src/bulk/sample.rs`
Expected: no output

- [ ] **Step 5: Commit**

```bash
git add src/bulk/profile.rs src/bulk/sample.rs tests/bulk.rs
git commit -m "refactor(bulk): route profile-content errors to InvalidProfile (#16)

Also tightens the PL/AD ploidy test, which previously passed on any
validation failure rather than the guard it names."
```

---

### Task 3: Route spec, parsing, and runtime errors in `mod.rs`

**Files:**
- Modify: `src/bulk/mod.rs:278` (`NoContigs`), `:281` (`NoSamples`), `:294`
  (`DuplicateContig`), `:313` (`WorkerPool`), `:561` (`TargetNotReached`),
  `:729` (`BadSize` in `parse_size`)
- Test: `tests/bulk.rs:457` (tighten existing assertion), plus two new tests
- Test: `tests/cli.rs:5-12` (tighten `parse_size` assertion)

**Interfaces:**
- Consumes: `NoContigs`, `NoSamples`, `DuplicateContig`, `BadSize`,
  `WorkerPool`, `TargetNotReached` from Task 1.
- Produces: nothing new. `BulkSpec::write` and `parse_size` keep their
  signatures.

- [ ] **Step 1: Write the failing tests**

(a) In `tests/bulk.rs`, in `duplicate_contig_names_are_rejected` (~line 457),
replace:

```rust
    assert!(
        matches!(result, Err(BulkError::Invalid(_))),
        "duplicate output contig names must be rejected as invalid: {result:?}"
    );
```

with:

```rust
    assert!(
        matches!(&result, Err(BulkError::DuplicateContig(id)) if id == "chr1"),
        "duplicate output contig names must be rejected with DuplicateContig \
         naming the offending contig: {result:?}"
    );
```

(b) Append two new tests to the end of `tests/bulk.rs`:

```rust
/// An empty contig list and a zero sample count are caller mistakes, not
/// profile mistakes -- they must not be reported as an invalid profile.
#[test]
fn empty_spec_dimensions_are_rejected_as_spec_errors() {
    let dir = tempfile::tempdir().unwrap();

    let no_contigs = spec()
        .contigs(Vec::<String>::new())
        .size(Size::RecordsPerContig(10))
        .write(dir.path().join("a.bcf"));
    assert!(
        matches!(no_contigs, Err(BulkError::NoContigs)),
        "an empty contig list must be a spec error: {no_contigs:?}"
    );

    let no_samples = spec()
        .samples(0)
        .size(Size::RecordsPerContig(10))
        .write(dir.path().join("b.bcf"));
    assert!(
        matches!(no_samples, Err(BulkError::NoSamples)),
        "a zero sample count must be a spec error: {no_samples:?}"
    );

    for e in [no_contigs, no_samples] {
        let msg = e.unwrap_err().to_string();
        assert!(
            !msg.starts_with("invalid profile:"),
            "a spec error must not blame the profile, got: {msg}"
        );
    }
}
```

(c) In `tests/cli.rs`, replace `assert!(parse_size("banana").is_err());` with:

```rust
    let bad = parse_size("banana");
    assert!(
        matches!(&bad, Err(vcfixture::bulk::BulkError::BadSize(s)) if s == "banana"),
        "an unparseable size is an argument error, not an invalid profile: {bad:?}"
    );
    assert!(!bad.unwrap_err().to_string().starts_with("invalid profile:"));
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --all-features --test bulk duplicate_contig_names_are_rejected`
Expected: FAIL — assertion trips, still `Invalid(_)`.

Run: `cargo test --all-features --test bulk empty_spec_dimensions_are_rejected_as_spec_errors`
Expected: FAIL — assertion trips.

Run: `cargo test --all-features --test cli parses_a_size_with_units`
Expected: FAIL — assertion trips.

- [ ] **Step 3: Write minimal implementation**

In `src/bulk/mod.rs`, make six edits.

`:278` and `:281` —

```rust
        if self.contig_ids.is_empty() {
            return Err(BulkError::NoContigs);
        }
        if self.n_samples == 0 {
            return Err(BulkError::NoSamples);
        }
```

`:294` — inside the duplicate-name loop, replace the `format!` with:

```rust
                    return Err(BulkError::DuplicateContig(id.clone()));
```

`:313` — the worker pool. `WorkerPool` has `#[from]`, so the closure goes away:

```rust
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.workers.get())
            .build()?;
```

`:561` — the target-size give-up. Hoist `MAX_CORRECTIONS` so the error carries
the real bound instead of interpolating it into a string:

```rust
        Err(BulkError::TargetNotReached {
            target_bytes,
            corrections: MAX_CORRECTIONS,
        })
```

`MAX_CORRECTIONS` is declared as `const MAX_CORRECTIONS: usize = 4;` at
`src/bulk/mod.rs:536` — a function-body const, not a loop-scoped one, and the
line you are replacing already interpolates it. It is in scope; no move needed.

`:729` — `parse_size`:

```rust
        .map_err(|_| BulkError::BadSize(s.to_string()))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --all-features`
Expected: PASS

Run: `cargo clippy --all-features --all-targets -- -D warnings`
Expected: clean

Run: `rg -n "BulkError::Invalid\(" src/bulk/mod.rs`
Expected: no output

- [ ] **Step 5: Commit**

```bash
git add src/bulk/mod.rs tests/bulk.rs tests/cli.rs
git commit -m "refactor(bulk): route spec, parsing, and runtime errors (#16)

Spec, argument, and runtime failures no longer render as 'invalid
profile'. TargetNotReached carries MAX_CORRECTIONS instead of
interpolating it, so the message cannot drift from the loop bound. Also
tightens the duplicate-contig test, which previously passed on any
failure."
```

---

### Task 4: Route writer and CLI errors, and make `--threads` unrepresentably-zero

**Files:**
- Modify: `src/bulk/writer.rs:66`
- Modify: `src/bin/vcfixture.rs:72-73` (the `threads` arg), `:117-133`
  (`resolve_profile`), `:175-179` (the manual `NonZero` check)
- Test: `src/bulk/writer.rs` (existing `#[cfg(test)] mod tests` at line ~142)
- Test: `tests/cli.rs` (new test for the `--threads 0` usage error)

**Interfaces:**
- Consumes: `BulkError::CompressionLevel`, `BulkError::ProfileLoad` from Task 1.
- Produces: `fn parse_threads(s: &str) -> Result<NonZero<usize>, String>` in
  `src/bin/vcfixture.rs` (private to the binary; no other task uses it). The
  `Cmd::Bulk` field `threads` changes type from `Option<usize>` to
  `Option<NonZero<usize>>`.

- [ ] **Step 1: Write the failing tests**

(a) Append to `mod tests` in `src/bulk/writer.rs`:

```rust
    /// An out-of-range bgzf compression level is an argument error. It is
    /// the caller's number that is wrong, not the profile.
    #[test]
    fn out_of_range_compression_level_is_an_argument_error() {
        let dir = tempfile::tempdir().unwrap();
        // `header()` is the existing helper at the top of this test module.
        let result = BulkWriter::create(
            &dir.path().join("a.bcf"),
            Format::Bcf,
            &header(),
            99,
            NonZero::new(1).unwrap(),
        );
        assert!(
            matches!(result, Err(BulkError::CompressionLevel(_))),
            "compression level 99 is out of range and must be an argument \
             error, not an invalid profile"
        );
        assert!(!result
            .err()
            .unwrap()
            .to_string()
            .starts_with("invalid profile:"));
    }
```

If `BulkWriter::create` returns a non-`Debug` success type, bind with
`let Err(e) = result else { panic!("level 99 must be rejected") };` and assert
on `e` instead of using `matches!` on the whole `Result`.

(b) Append to `tests/cli.rs`:

```rust
/// `--threads 0` is rejected by clap as a usage error before any library
/// code runs, so it can never reach BulkError at all.
#[test]
fn zero_threads_is_a_clap_usage_error() {
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_vcfixture"))
        .args(["bulk", "--threads", "0", "-o", "unused.bcf"])
        .output()
        .expect("binary should run");
    assert!(!out.status.success(), "--threads 0 must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("must be >= 1"),
        "clap should explain the constraint, got: {stderr}"
    );
    assert!(
        !stderr.contains("invalid profile"),
        "a bad --threads value must not blame the profile, got: {stderr}"
    );
    assert!(
        !std::path::Path::new("unused.bcf").exists(),
        "nothing should be written for a usage error"
    );
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --all-features --lib bulk::writer::tests::out_of_range_compression_level_is_an_argument_error`
Expected: FAIL — still returns `Invalid(_)`.

Run: `cargo test --all-features --test cli zero_threads_is_a_clap_usage_error`
Expected: FAIL — stderr reads `error: invalid profile: --threads must be >= 1`.

- [ ] **Step 3: Write minimal implementation**

(a) `src/bulk/writer.rs:65-66`:

```rust
        let level = bgzf::io::writer::CompressionLevel::try_from(compression_level)
            .map_err(|e| BulkError::CompressionLevel(e.to_string()))?;
```

(b) `src/bin/vcfixture.rs`, `resolve_profile` (~line 123):

```rust
            let text = std::fs::read_to_string(name_or_path).map_err(|e| {
                BulkError::ProfileLoad {
                    path: name_or_path.to_string(),
                    source: e,
                }
            })?;
```

(c) `src/bin/vcfixture.rs`, the `threads` argument (~line 72). Clap 4.6 has no
built-in `NonZero` value parser (verified against the locked `clap_builder`
4.6.2 source), so supply one:

```rust
        /// Worker thread count for compression/generation. Defaults to all
        /// available cores.
        #[arg(long, value_parser = parse_threads)]
        threads: Option<NonZero<usize>>,
```

and add, next to `resolve_profile`:

```rust
/// Parses `--threads` straight into a `NonZero<usize>` so zero is rejected
/// by clap's own usage error and never becomes a library error variant.
fn parse_threads(s: &str) -> Result<NonZero<usize>, String> {
    let n: usize = s
        .parse()
        .map_err(|_| format!("{s:?} is not a thread count"))?;
    NonZero::new(n).ok_or_else(|| "must be >= 1".to_string())
}
```

(d) `src/bin/vcfixture.rs`, `run()` (~line 175) — the manual check disappears:

```rust
    if let Some(n) = threads {
        spec = spec.workers(n);
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --all-features`
Expected: PASS

Run: `cargo clippy --all-features --all-targets -- -D warnings`
Expected: clean

Run: `rg -n "BulkError::Invalid\(" src/`
Expected: no output — every call site in `src/` is now routed.

- [ ] **Step 5: Commit**

```bash
git add src/bulk/writer.rs src/bin/vcfixture.rs tests/cli.rs
git commit -m "refactor(bulk): route writer and CLI errors, parse --threads as NonZero (#16)

--threads is now parsed into NonZero<usize> by clap, so the zero case is
a usage error rather than a BulkError variant."
```

---

### Task 5: Delete `Invalid` and verify end to end

**Files:**
- Modify: `src/bulk/mod.rs` (remove the `Invalid` variant)
- Verify: `src/bin/validate_profile.rs` (must not match on `Invalid`)

**Interfaces:**
- Consumes: everything from Tasks 2–4.
- Produces: the final public `BulkError`. `Invalid` no longer exists.

- [ ] **Step 1: Confirm nothing still references `Invalid`**

Run: `rg -n "BulkError::Invalid" src/ tests/ examples/ docs/ benches/ 2>/dev/null`
Expected: matches only inside
`docs/superpowers/specs/2026-08-05-bulk-error-split-design.md` and this plan
(both quote the old code deliberately). Any hit in `src/`, `tests/`,
`examples/`, or `benches/` means an earlier task is incomplete — go finish it
before continuing.

Also check the second binary specifically, since it postdates the issue:

Run: `rg -n "BulkError" src/bin/validate_profile.rs`
Expected: matches only in a `use` statement or a `Result<_, BulkError>` return
type — no `match` arm on `Invalid`. If it does match on `Invalid`, route it
like Task 2 (it validates profiles) and note the deviation in your report.

- [ ] **Step 2: Delete the variant**

Remove these four lines from the `BulkError` enum in `src/bulk/mod.rs`:

```rust
    /// Deprecated catch-all, removed in Task 5 once every call site is
    /// routed. Do not add new uses.
    #[error("invalid profile: {0}")]
    Invalid(String),
```

- [ ] **Step 3: Verify the tree still builds and every gate passes**

Run: `cargo test --all-features`
Expected: PASS

Run: `cargo clippy --all-features --all-targets -- -D warnings`
Expected: clean

Run: `cargo fmt --check`
Expected: no output

Run: `cargo doc --all-features --no-deps`
Expected: builds without warnings (the enum's new doc comment references
`Profile::validate`, which must resolve)

- [ ] **Step 4: Reproduce the issue's own repro and confirm it is fixed**

Run:

```bash
cargo run --all-features --bin vcfixture -- bulk --contigs chr1,chr1 \
  --records-per-contig 10 -o /tmp/issue16.bcf; echo "exit=$?"
```

Expected: exits non-zero, and stderr reads

```
error: duplicate output contig name: "chr1" (each requested contig must be unique; ...)
```

**not** `error: invalid profile: duplicate output contig name: ...`. Paste the
actual output into your task report — this is the user-visible symptom the
whole change exists to fix, and it is the one thing no unit test observes end
to end.

Then confirm the profile path still reads correctly:

```bash
cargo run --all-features --bin vcfixture -- bulk --profile /nonexistent.json \
  -o /tmp/issue16.bcf; echo "exit=$?"
```

Expected: `error: profile "/nonexistent.json" is not a builtin name and could
not be read as a file: No such file or directory (os error 2)`

- [ ] **Step 5: Commit**

The `!` and the `BREAKING CHANGE:` footer are what make `cz bump` classify this
correctly and write the `CHANGELOG.md` entry. Do not edit `CHANGELOG.md` by
hand.

```bash
git add src/bulk/mod.rs
git commit -m "refactor(bulk)!: remove the BulkError::Invalid catch-all (#16)

Every call site now routes to a variant that describes its own failure
class, so spec, argument, runtime, and profile-loading errors are no
longer rendered as 'invalid profile: ...'.

Closes #16

BREAKING CHANGE: BulkError::Invalid is removed. Profile-validation
failures are now BulkError::InvalidProfile (same message, same prefix);
other failures use NoContigs, NoSamples, DuplicateContig, PayloadPloidy,
BadSize, CompressionLevel, ProfileLoad, WorkerPool, or TargetNotReached."
```

- [ ] **Step 6: Push and open a PR**

```bash
git push -u origin worktree-bulk-error-split
gh pr create --fill --draft
```

---

## Self-Review

**Spec coverage** — every section of
`docs/superpowers/specs/2026-08-05-bulk-error-split-design.md` maps to a task:

| Spec requirement | Task |
|---|---|
| New enum shape, all 10 new variants | 1 |
| `InvalidProfile` keeps prefix, 13 messages behind it | 2 |
| `PayloadPloidy` structured, keeps prefix | 2 |
| `NoContigs` / `NoSamples` / `DuplicateContig` | 3 |
| `BadSize`, `CompressionLevel` | 3, 4 |
| `ProfileLoad` | 4 |
| `WorkerPool` via `#[from]`, `TargetNotReached` carries `MAX_CORRECTIONS` | 3 |
| `--threads` becomes a clap `NonZero` parser | 4 |
| Tighten the two over-broad assertions | 2, 3 |
| Tests for new spec/parsing variants | 1, 3, 4 |
| Profile failures still carry the prefix | 1 (unit), 2 (integration) |
| CLI smoke check | 5 |
| `Invalid` removed, no deprecation shim | 5 |
| Verify `validate_profile.rs` | 5 |
| Breaking change declared | 5 (commit footer) |

**Deviation from the spec, deliberate:** the spec says `CHANGELOG.md` gets a
hand-written `Refactor` entry. It does not — `.cz.toml` sets
`update_changelog_on_bump = true`, so `cz bump` generates the file from commit
messages, and a hand-edit would be overwritten. The breaking change is declared
in Task 5's commit footer instead. This is recorded in Global Constraints.

**Placeholder scan:** no TBDs; every code step carries the actual code; the two
"if it doesn't compile, do this instead" notes (the `MAX_CORRECTIONS` scope
check in Task 3, the `Debug`-bound fallback in Task 4) are explicit
alternatives with concrete instructions, not deferred decisions.

**Type consistency:** `InvalidProfile(String)`, `PayloadPloidy { payload:
Payload, ploidy: u8 }`, `DuplicateContig(String)`, `BadSize(String)`,
`CompressionLevel(String)`, `ProfileLoad { path: String, source:
std::io::Error }`, `TargetNotReached { target_bytes: u64, corrections: usize }`
are spelled identically in Task 1's definition and in every consuming task.
`ploidy: u8` matches `Dialed::ploidy` (`src/bulk/profile.rs:83`);
`corrections: usize` matches `const MAX_CORRECTIONS: usize`; `target_bytes:
u64` matches `Size::Target(u64)`. `parse_threads` returns
`NonZero<usize>`, matching `BulkSpec::workers`' parameter type.
