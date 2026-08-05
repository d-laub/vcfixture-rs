# Split `BulkError::Invalid` into per-class variants (issue #16)

**Date:** 2026-08-05
**Status:** approved, ready for planning
**Scope:** one public error enum (`src/bulk/mod.rs`), its 25 call sites, and
the tests that match on it. Breaking change to a `0.x` public API.

## Context

`BulkError::Invalid` is declared as:

```rust
#[error("invalid profile: {0}")]
Invalid(String),
```

but it is the catch-all for 25 call sites spanning four unrelated failure
classes. For three of the four the `invalid profile:` prefix is a lie. A caller
who passes a duplicate contig name is told their *profile* is invalid:

```
$ vcfixture bulk --contigs chr1,chr1 -o out.bcf
error: invalid profile: duplicate output contig name: "chr1" (each requested
contig must be unique; ...)
```

There is no one-line fix. Dropping the prefix to `"{0}"` de-contextualises the
profile-validation messages (`"ploidy must be >= 1"`, `"histogram needs >= 2
edges"`), which are not self-describing and rely on it.

A second, unreported problem shares the same root. Because every validation
failure collapses to one variant, tests can only assert
`matches!(result, Err(BulkError::Invalid(_)))` — which passes if the call
failed for *any* reason. `duplicate_contig_names_are_rejected`
(`tests/bulk.rs:457`) would pass even if `write()` had failed on an unrelated
profile check, never exercising the guard it names. Splitting the variant fixes
the messages and makes those assertions precise.

### Call-site inventory (at `origin/main`, v0.3.0)

| Class | Sites | Files |
|---|---|---|
| Profile content | 14 | `profile.rs` (13), `sample.rs` (1) |
| Spec / caller | 3 | `mod.rs` |
| Argument parsing | 3 | `mod.rs` (`parse_size`), `writer.rs`, `bin/vcfixture.rs` |
| Runtime | 2 | `mod.rs` |
| Profile loading | 1 | `bin/vcfixture.rs` |
| Test assertions | 2 | `tests/bulk.rs` |

Note: the PL/AD-vs-ploidy guard moved from `BulkSpec::write` into
`Profile::validate` in the v0.2/v0.3 work (`payload` now lives in `Dialed`, and
the check runs at parse time). It is therefore a **profile-content** error on
this base, not a spec error — issue #16's inventory predates that move.
Issue #16 also refers to `Size::PerContig`'s two name-set errors; that variant
was never added — #15 was fixed by changing the split key to `n_variants`
(`distribute_by_n_variants`), so there is no such class to route.

## Design

### The new enum

```rust
#[derive(Debug, thiserror::Error)]
pub enum BulkError {
    #[error("unknown builtin profile: {0}")]
    UnknownProfile(String),

    // --- profile content --------------------------------------------------
    /// Fitted/dialed statistics that fail `Profile::validate`. The prefix is
    /// load-bearing: these messages name a profile field, not themselves.
    #[error("invalid profile: {0}")]
    InvalidProfile(String),

    /// Split out of `InvalidProfile` because it is the one profile-content
    /// failure a caller can act on without editing the profile JSON — by
    /// choosing a different payload.
    #[error(
        "invalid profile: payload {payload:?} emits AD and/or PL, which are \
         hard-coded for diploid (ploidy 2) calls, but ploidy is {ploidy}"
    )]
    PayloadPloidy { payload: Payload, ploidy: u8 },

    // --- spec / caller validation ----------------------------------------
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
    #[error(
        "profile {path:?} is not a builtin name and could not be read as a \
         file: {source}"
    )]
    ProfileLoad { path: String, source: std::io::Error },

    // --- runtime ----------------------------------------------------------
    #[error("failed to build worker pool: {0}")]
    WorkerPool(#[from] rayon::ThreadPoolBuildError),
    #[error(
        "could not reach target size {target_bytes} bytes within \
         {corrections} corrective rounds"
    )]
    TargetNotReached { target_bytes: u64, corrections: usize },

    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}
```

`Invalid` is removed outright.

### The rule that decides the split

**Structure where a caller can branch; string where a caller can only
display.**

Thirteen of the fourteen profile-content messages stay behind a single
`InvalidProfile(String)`.
Nobody can usefully do something different for `"histogram edges must be
increasing"` than for `"histogram weights must sum > 0"` — both mean *your
profile JSON is malformed, here is where*. Giving each its own variant would
add fifteen variants that only ever reach `Display`, and would make every
future profile check a public API addition. The one exception is
`PayloadPloidy`: unlike the others it is not a malformed profile but a
mismatched *choice*, actionable by picking another payload, and it is the
subject of an existing test that deserves a precise assertion.

Everything outside profile content gets a named variant. There are only eight,
each is a distinct failure a caller might handle, and each message is
self-describing without a prefix.

This matches the crate's existing house style: `BuildError` in `src/error.rs`
is ~30 structured variants with named fields.

### Sites that disappear rather than move

- **`--threads must be >= 1`** (`bin/vcfixture.rs:177`). The manual
  `NonZero::new(n).ok_or_else(...)` becomes a clap `value_parser` that parses
  straight into `NonZero<usize>`, so the invalid state is unrepresentable and
  clap renders the usage error. Field type changes to
  `Option<NonZero<usize>>`.
- **`MAX_CORRECTIONS` interpolation** (`mod.rs:561`). The constant moves into
  `TargetNotReached.corrections`, so the message cannot drift from the loop
  bound.

### Error-message text

Messages are preserved verbatim except where the prefix removal requires a
change. `PayloadPloidy` keeps its `invalid profile:` prefix (it *is* a profile
error). Spec, parsing, runtime, and loading variants lose the prefix and keep
the rest of their wording, so existing message-substring assertions elsewhere
continue to match.

## Testing

- **Tighten the two existing assertions.** `tests/bulk.rs:457` →
  `Err(BulkError::DuplicateContig(_))`; `tests/bulk.rs:481` →
  `Err(BulkError::PayloadPloidy { .. })`. These are the assertions that
  currently pass on any failure.
- **One test per new spec/parsing variant not already covered**: `NoContigs`,
  `NoSamples`, `BadSize` (a few malformed inputs plus the valid `KB/MB/GB`
  suffixes as a regression guard on `parse_size`).
- **One test that profile-content failure still surfaces as
  `InvalidProfile`** and still carries the `invalid profile:` prefix in its
  `Display` output — this is the property the whole split exists to protect.
- **CLI smoke check**: run `vcfixture bulk --contigs chr1,chr1` and confirm
  stderr no longer begins `invalid profile:`. This reproduces the issue's own
  repro (adapted — the issue's `--records-for` flag does not exist).
- **Gates**: `cargo test --all-features`, `cargo clippy --all-features -- -D
  warnings`, `cargo fmt --check`.

## Compatibility

Breaking change to a public enum. `BulkError` is behind the `bulk` feature and
the crate is `0.3.0`, so it ships as a minor bump under Cargo's `0.x` rules
with no deprecation shim — a shim is impossible for a removed enum variant
anyway, and any `match` on `BulkError` is non-exhaustive-safe only if the
caller has a wildcard arm.

`CHANGELOG.md` gets a `Refactor` entry flagging the break, in the existing
commitizen style (`.cz.toml`).

Also verify `src/bin/validate_profile.rs` (added in v0.2) does not pattern-match
`BulkError::Invalid`; a grep at spec time shows it does not, but the
implementation must confirm after the rename.

## Out of scope

- Converting the remaining 13 profile-content messages into structured
  variants.
- Introducing nested sub-enums (`ProfileError`, `SpecError`). The crate's
  existing `BuildError` is flat; a flat `BulkError` stays consistent and keeps
  `?` working without extra `From` impls.
- Any change to validation *behaviour*. This spec re-routes existing errors; it
  adds and removes no checks.
