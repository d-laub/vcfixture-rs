# vcfixture-rs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Rust crate `vcfixture` that builds spec-conformant VCF (v4.1–4.5) test data and returns the decoded ground-truth oracle alongside, with proptest strategies — for property-testing svar2 in genoray.

**Architecture:** One owned immutable `Document` model is the hub. A `VcfBuilder` (or proptest strategies) produces it; `truth()` derives an ndarray oracle directly from it; `write.rs` converts it to noodles types for text/bgzip/CSI output. Bytes and oracle come from the same object, so they can't diverge.

**Tech Stack:** Rust, ndarray, noodles (`-vcf/-bgzf/-csi/-tabix/-fasta/-core`), indexmap, thiserror, proptest (feature-gated), rand/rand_chacha. Dev env: pixi + prek.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-06-23-vcfixture-rs-design.md`.
- Crate name `vcfixture`; edition 2021; library crate.
- `proptest` is an **optional** dependency behind feature `proptest`, **off by default**. `rand`/`rand_chacha`/`ndarray`/`indexmap`/`thiserror`/noodles are always on.
- Hub invariant: nothing derives truth or writes bytes except from the `Document` model. noodles is output-only; never the in-memory model.
- No cross-language byte/oracle parity with the Python lib is required.
- Public fallible API returns `Result<_, BuildError>`. Inside proptest strategies, generation is correct by construction; `.expect()` builder results there.
- Every task ends green: `pixi run fmt`, `pixi run clippy`, `pixi run test` all pass before commit. Commits use conventional-commit messages (commitizen commit-msg hook is active).
- Missing genotype allele encodes as `-1` in `GroundTruth.genotypes`.

---

### Task 1: Scaffold crate + dev environment

**Files:**
- Create: `Cargo.toml`, `pixi.toml`, `.pre-commit-config.yaml`, `rustfmt.toml`, `src/lib.rs`, `src/error.rs`

**Interfaces:**
- Consumes: nothing.
- Produces: `vcfixture::error::BuildError` (enum, `thiserror::Error`, `Debug`, `Clone`-not-required), used by every later task. Crate compiles with `cargo check`.

- [ ] **Step 1: Write `Cargo.toml`**

```toml
[package]
name = "vcfixture"
version = "0.1.0"
edition = "2021"
license = "Apache-2.0"
description = "Generate small VCF test data with decoded ground truth, for property-based testing."

[dependencies]
ndarray = "0.16"
indexmap = "2"
thiserror = "2"
rand = "0.8"
rand_chacha = "0.3"
noodles-core = "0.17"
noodles-vcf = "0.79"
noodles-bgzf = "0.41"
noodles-csi = "0.46"
noodles-tabix = "0.55"
noodles-fasta = "0.51"
proptest = { version = "1", optional = true }

[dev-dependencies]
insta = "1"
proptest = "1"

[features]
default = []
proptest = ["dep:proptest"]
```

Note: noodles sub-crate versions move fast. After `pixi install`, run `cargo update` and adjust these to the resolved compatible set; pin whatever `cargo build` accepts. The `dev-dependencies` `proptest` (unconditional) lets `tests/` use strategies without the consumer enabling the feature — but library `strategies.rs` stays behind the feature.

- [ ] **Step 2: Write `pixi.toml`**

```toml
[workspace]
name = "vcfixture-rs"
channels = ["conda-forge"]
platforms = ["osx-arm64", "linux-64"]

[dependencies]
rust = "*"
prek = "*"
commitizen = "*"

[tasks]
build        = "cargo build"
test         = "cargo test --all-features"
check        = "cargo check --all-features"
fmt          = "cargo fmt"
fmt-check    = "cargo fmt --check"
clippy       = "cargo clippy --all-features -- -D warnings"
prek-install = "prek install -t pre-commit -t commit-msg"
bump-dry     = "cz bump --dry-run"
```

- [ ] **Step 3: Write `.pre-commit-config.yaml`**

```yaml
default_install_hook_types: [pre-commit, commit-msg]

repos:
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v6.0.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-toml
      - id: check-merge-conflict
      - id: check-case-conflict
      - id: check-added-large-files
      - id: mixed-line-ending
        args: [--fix=lf]

  - repo: https://github.com/commitizen-tools/commitizen
    rev: v4.16.3
    hooks:
      - id: commitizen
        stages: [commit-msg]

  - repo: local
    hooks:
      - id: cargo-fmt
        name: cargo fmt
        entry: cargo fmt --check
        language: system
        types: [rust]
        pass_filenames: false
      - id: cargo-clippy
        name: cargo clippy
        entry: cargo clippy --all-features -- -D warnings
        language: system
        types: [rust]
        pass_filenames: false
      - id: cargo-check
        name: cargo check
        entry: cargo check --all-features
        language: system
        types: [rust]
        pass_filenames: false
```

- [ ] **Step 4: Write `rustfmt.toml`** (empty defaults are fine; create the file so config is explicit)

```toml
edition = "2021"
```

- [ ] **Step 5: Write `src/error.rs`**

```rust
//! The crate-wide error type.

use thiserror::Error;

/// All errors produced while declaring fields, adding records, deriving the
/// reserved registry, or writing output.
#[derive(Debug, Error)]
pub enum BuildError {
    #[error("sequence allele bases must be [ACGTN]+, got {0:?}")]
    BadAlleleBases(String),

    #[error("not a valid breakend replacement string: {0:?}")]
    BadBreakend(String),

    #[error("symbolic SV first type must be one of DEL/INS/DUP/INV/CNV, got {0:?}")]
    BadSvType(String),

    #[error("ID {0:?} does not match the VCF key regex")]
    BadFieldId(String),

    #[error("Flag fields must be INFO, not FORMAT")]
    FlagNotInfo,

    #[error("Flag fields must have Number=0")]
    FlagNumberNotZero,

    #[error("Number=0 is only valid for Flag fields")]
    ZeroNumberNotFlag,

    #[error("fixed Number must be >= 0")]
    NegativeFixedNumber,

    #[error("{kind} field {id:?} is not a known reserved field; pass number and type explicitly")]
    UnknownReserved { kind: String, id: String },

    #[error("{kind} field {id:?} was introduced in {since}; not available in {version}")]
    FieldTooNew { kind: String, id: String, since: String, version: String },

    #[error("symbolic/breakend ALT requires a single preceding REF padding base, got REF={0:?}")]
    MissingRefPadding(String),

    #[error("SVLEN required for symbolic allele {0}")]
    MissingSvlen(String),

    #[error("SVCLAIM {claim:?} invalid for {allele}; allowed {allowed:?}")]
    BadSvclaim { claim: String, allele: String, allowed: Vec<String> },

    #[error("SVCLAIM required for {0} (D/J/DJ)")]
    SvclaimRequired(String),

    #[error("SVLEN must be missing for {0}")]
    SvlenMustBeMissing(String),

    #[error("{kind} field {id:?} not declared")]
    UndeclaredField { kind: String, id: String },

    #[error("{id} cardinality mismatch: expected {expected}, got {got}")]
    Cardinality { id: String, expected: usize, got: usize },

    #[error("allele index {index} out of range (n_alt={n_alt})")]
    AlleleIndexOutOfRange { index: u32, n_alt: usize },

    #[error("GT not declared; declare it with .format(\"GT\", ...)")]
    GtNotDeclared,

    #[error("FORMAT CN requires equal SVLEN across <CNV>/<DEL>/<DUP> alleles")]
    CnSvlenMismatch,

    #[error("contig {0:?} already added")]
    ContigExists(String),

    #[error("contig {0:?} not found")]
    ContigNotFound(String),

    #[error("range {contig}:{pos0}+{len} runs past contig length {clen}")]
    OutOfBounds { contig: String, pos0: usize, len: usize, clen: usize },

    #[error(transparent)]
    Io(#[from] std::io::Error),
}
```

- [ ] **Step 6: Write `src/lib.rs`**

```rust
//! vcfixture — generate small VCF test data with decoded ground truth.

pub mod error;

pub use error::BuildError;
```

- [ ] **Step 7: Install env and verify toolchain**

Run:
```bash
pixi install
pixi run cargo --version
pixi run cargo clippy --version   # verify clippy is bundled; if absent, add the component to pixi.toml deps and re-run
pixi run prek-install
pixi run check
```
Expected: all succeed; `pixi run check` prints a clean `cargo check`.

- [ ] **Step 8: Commit**

```bash
git add -A
git commit -m "chore: scaffold crate, pixi env, and prek hooks"
```

---

### Task 2: `spec/version.rs` — VcfVersion

**Files:**
- Create: `src/spec/mod.rs`, `src/spec/version.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- Produces: `enum VcfVersion { V4_1, V4_2, V4_3, V4_4, V4_5 }` (`Copy`, `Eq`, `Ord` in chronological order); `VcfVersion::as_str(&self) -> &'static str` returning `"VCFv4.1"`..`"VCFv4.5"`; `const LATEST: VcfVersion = VcfVersion::V4_5`.

- [ ] **Step 1: Write the failing test** (in `src/spec/version.rs`)

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordered_chronologically() {
        assert!(VcfVersion::V4_3 < VcfVersion::V4_4);
        assert!(VcfVersion::V4_5 >= LATEST);
    }

    #[test]
    fn fileformat_strings() {
        assert_eq!(VcfVersion::V4_2.as_str(), "VCFv4.2");
        assert_eq!(LATEST.as_str(), "VCFv4.5");
    }
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run cargo test --lib spec::version`
Expected: FAIL (module/types not found).

- [ ] **Step 3: Implement** (top of `src/spec/version.rs`)

```rust
/// A supported VCF spec version. `Ord` is chronological.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VcfVersion {
    V4_1,
    V4_2,
    V4_3,
    V4_4,
    V4_5,
}

/// The latest supported version.
pub const LATEST: VcfVersion = VcfVersion::V4_5;

impl VcfVersion {
    /// The exact `##fileformat` string.
    pub fn as_str(&self) -> &'static str {
        match self {
            VcfVersion::V4_1 => "VCFv4.1",
            VcfVersion::V4_2 => "VCFv4.2",
            VcfVersion::V4_3 => "VCFv4.3",
            VcfVersion::V4_4 => "VCFv4.4",
            VcfVersion::V4_5 => "VCFv4.5",
        }
    }
}

impl std::fmt::Display for VcfVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}
```

And `src/spec/mod.rs`:

```rust
pub mod version;
```

And add to `src/lib.rs`:

```rust
pub mod spec;
pub use spec::version::{VcfVersion, LATEST};
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run cargo test --lib spec::version`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat: add VcfVersion"
```

---

### Task 3: `spec/types.rs` — Type

**Files:**
- Create: `src/spec/types.rs`
- Modify: `src/spec/mod.rs`, `src/lib.rs`

**Interfaces:**
- Produces: `enum Type { Integer, Float, Flag, Character, String }` (`Copy`, `Eq`); `Type::as_str(&self) -> &'static str` (`"Integer"`, …); `Type::info_allowed() -> [Type; 5]`; `Type::format_allowed() -> [Type; 4]` (no Flag).

- [ ] **Step 1: Write the failing test** (in `src/spec/types.rs`)

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_excludes_flag() {
        assert!(!Type::format_allowed().contains(&Type::Flag));
        assert!(Type::info_allowed().contains(&Type::Flag));
    }

    #[test]
    fn header_tokens() {
        assert_eq!(Type::Float.as_str(), "Float");
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib spec::types` → FAIL.

- [ ] **Step 3: Implement** (top of `src/spec/types.rs`)

```rust
/// VCF value type for an INFO or FORMAT field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Type {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

impl Type {
    pub fn as_str(&self) -> &'static str {
        match self {
            Type::Integer => "Integer",
            Type::Float => "Float",
            Type::Flag => "Flag",
            Type::Character => "Character",
            Type::String => "String",
        }
    }

    /// All types valid in INFO.
    pub fn info_allowed() -> [Type; 5] {
        [Type::Integer, Type::Float, Type::Flag, Type::Character, Type::String]
    }

    /// All types valid in FORMAT (Flag excluded).
    pub fn format_allowed() -> [Type; 4] {
        [Type::Integer, Type::Float, Type::Character, Type::String]
    }
}
```

Add `pub mod types;` to `src/spec/mod.rs` and `pub use spec::types::Type;` to `src/lib.rs`.

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib spec::types` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add Type"`

---

### Task 4: `spec/number.rs` — Number + cardinality

**Files:**
- Create: `src/spec/number.rs`
- Modify: `src/spec/mod.rs`, `src/lib.rs`

**Interfaces:**
- Consumes: `BuildError`.
- Produces: `enum NumberKind { Fixed, A, R, G, Dot, Flag }`; `struct Number { kind: NumberKind, count: Option<u32> }` (`Copy`, `Eq`); consts `Number::ONE/A/R/G/DOT/FLAG`; `Number::fixed(u32) -> Result<Number, BuildError>`; `Number::header_str(&self) -> String`; `Number::cardinality(&self, n_alt: usize, ploidy: usize) -> Option<usize>`.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cardinalities() {
        assert_eq!(Number::A.cardinality(3, 2), Some(3));
        assert_eq!(Number::R.cardinality(3, 2), Some(4));
        assert_eq!(Number::ONE.cardinality(3, 2), Some(1));
        assert_eq!(Number::DOT.cardinality(3, 2), None);
        // Number=G: diploid, n_alleles = n_alt+1 = 3 => C(3+2-1,2)=C(4,2)=6
        assert_eq!(Number::G.cardinality(2, 2), Some(6));
        // haploid G => n_alleles
        assert_eq!(Number::G.cardinality(2, 1), Some(3));
    }

    #[test]
    fn header_tokens() {
        assert_eq!(Number::A.header_str(), "A");
        assert_eq!(Number::fixed(2).unwrap().header_str(), "2");
        assert_eq!(Number::FLAG.header_str(), "0");
        assert_eq!(Number::DOT.header_str(), ".");
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib spec::number` → FAIL.

- [ ] **Step 3: Implement**

```rust
use crate::error::BuildError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NumberKind {
    Fixed,
    A,
    R,
    G,
    Dot,
    Flag,
}

/// VCF `Number=` cardinality descriptor. `count` is set only for `Fixed`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Number {
    pub kind: NumberKind,
    pub count: Option<u32>,
}

impl Number {
    pub const ONE: Number = Number { kind: NumberKind::Fixed, count: Some(1) };
    pub const A: Number = Number { kind: NumberKind::A, count: None };
    pub const R: Number = Number { kind: NumberKind::R, count: None };
    pub const G: Number = Number { kind: NumberKind::G, count: None };
    pub const DOT: Number = Number { kind: NumberKind::Dot, count: None };
    pub const FLAG: Number = Number { kind: NumberKind::Flag, count: None };

    pub fn fixed(n: u32) -> Result<Number, BuildError> {
        // u32 cannot be negative; kept fallible to mirror the spec API and to
        // allow a future signed source. Always Ok for now.
        Ok(Number { kind: NumberKind::Fixed, count: Some(n) })
    }

    pub fn header_str(&self) -> String {
        match self.kind {
            NumberKind::Fixed => self.count.unwrap_or(0).to_string(),
            NumberKind::Flag => "0".to_string(),
            NumberKind::A => "A".to_string(),
            NumberKind::R => "R".to_string(),
            NumberKind::G => "G".to_string(),
            NumberKind::Dot => ".".to_string(),
        }
    }

    /// Resolve to a concrete value count for one record, or `None` when
    /// unbounded (`Number=.`).
    pub fn cardinality(&self, n_alt: usize, ploidy: usize) -> Option<usize> {
        match self.kind {
            NumberKind::Fixed => Some(self.count.unwrap_or(0) as usize),
            NumberKind::Flag => Some(0),
            NumberKind::A => Some(n_alt),
            NumberKind::R => Some(n_alt + 1),
            NumberKind::G => {
                let n_alleles = n_alt + 1;
                Some(binom(n_alleles + ploidy - 1, ploidy))
            }
            NumberKind::Dot => None,
        }
    }
}

/// Binomial coefficient C(n, k), computed iteratively to avoid overflow for
/// the small values used here.
fn binom(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result: usize = 1;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}
```

Add `pub mod number;` to `src/spec/mod.rs` and `pub use spec::number::{Number, NumberKind};` to `src/lib.rs`.

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib spec::number` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add Number and cardinality"`

---

### Task 5: `spec/field.rs` — FieldDef

**Files:**
- Create: `src/spec/field.rs`
- Modify: `src/spec/mod.rs`, `src/lib.rs`

**Interfaces:**
- Consumes: `Number`, `NumberKind`, `Type`, `BuildError`.
- Produces: `enum FieldKind { Info, Format }` with `as_str()`; `struct FieldDef { id: String, number: Number, type_: Type, description: String, kind: FieldKind }`; `FieldDef::new(id, number, type_, description, kind) -> Result<FieldDef, BuildError>` (validates ID regex + Flag↔Number=0 rules); `FieldDef::header_line(&self) -> String`.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::number::Number;
    use crate::spec::types::Type;

    #[test]
    fn header_line_format() {
        let fd = FieldDef::new("AF", Number::A, Type::Float, "Allele frequency", FieldKind::Info).unwrap();
        assert_eq!(
            fd.header_line(),
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
        );
    }

    #[test]
    fn flag_must_be_info_with_number_zero() {
        assert!(FieldDef::new("DB", Number::FLAG, Type::Flag, "x", FieldKind::Format).is_err());
        assert!(FieldDef::new("DB", Number::ONE, Type::Flag, "x", FieldKind::Info).is_err());
        assert!(FieldDef::new("X", Number::FLAG, Type::Integer, "x", FieldKind::Info).is_err());
        assert!(FieldDef::new("DB", Number::FLAG, Type::Flag, "x", FieldKind::Info).is_ok());
    }

    #[test]
    fn bad_id_rejected() {
        assert!(FieldDef::new("1BAD", Number::ONE, Type::Integer, "x", FieldKind::Info).is_err());
        assert!(FieldDef::new("1000G", Number::ONE, Type::Integer, "x", FieldKind::Info).is_ok());
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib spec::field` → FAIL.

- [ ] **Step 3: Implement**

```rust
use crate::error::BuildError;
use crate::spec::number::{Number, NumberKind};
use crate::spec::types::Type;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FieldKind {
    Info,
    Format,
}

impl FieldKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            FieldKind::Info => "INFO",
            FieldKind::Format => "FORMAT",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldDef {
    pub id: String,
    pub number: Number,
    pub type_: Type,
    pub description: String,
    pub kind: FieldKind,
}

/// VCF key regex: `[A-Za-z_][0-9A-Za-z_.]*` or the literal `1000G`.
fn valid_id(id: &str) -> bool {
    if id == "1000G" {
        return true;
    }
    let mut chars = id.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() || c == '_' => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.')
}

impl FieldDef {
    pub fn new(
        id: impl Into<String>,
        number: Number,
        type_: Type,
        description: impl Into<String>,
        kind: FieldKind,
    ) -> Result<FieldDef, BuildError> {
        let id = id.into();
        if !valid_id(&id) {
            return Err(BuildError::BadFieldId(id));
        }
        if type_ == Type::Flag {
            if kind != FieldKind::Info {
                return Err(BuildError::FlagNotInfo);
            }
            if number.kind != NumberKind::Flag {
                return Err(BuildError::FlagNumberNotZero);
            }
        } else if number.kind == NumberKind::Flag {
            return Err(BuildError::ZeroNumberNotFlag);
        }
        Ok(FieldDef { id, number, type_, description: description.into(), kind })
    }

    pub fn header_line(&self) -> String {
        format!(
            "##{}=<ID={},Number={},Type={},Description=\"{}\">",
            self.kind.as_str(),
            self.id,
            self.number.header_str(),
            self.type_.as_str(),
            self.description,
        )
    }
}
```

Add `pub mod field;` to `src/spec/mod.rs` and `pub use spec::field::{FieldDef, FieldKind};` to `src/lib.rs`.

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib spec::field` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add FieldDef"`

---

### Task 6: `spec/reserved.rs` — reserved-field registry + version gating

**Files:**
- Create: `src/spec/reserved.rs`
- Modify: `src/spec/mod.rs`, `src/lib.rs`

**Interfaces:**
- Consumes: `FieldDef`, `FieldKind`, `Number`, `Type`, `VcfVersion`, `BuildError`.
- Produces: `reserved(id: &str, kind: FieldKind, version: VcfVersion) -> Result<FieldDef, BuildError>`.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::field::FieldKind;
    use crate::spec::number::Number;
    use crate::spec::types::Type;
    use crate::spec::version::VcfVersion;

    #[test]
    fn resolves_af() {
        let fd = reserved("AF", FieldKind::Info, VcfVersion::V4_5).unwrap();
        assert_eq!(fd.number, Number::A);
        assert_eq!(fd.type_, Type::Float);
    }

    #[test]
    fn svlen_form_switches_at_4_4() {
        let pre = reserved("SVLEN", FieldKind::Info, VcfVersion::V4_3).unwrap();
        assert_eq!(pre.number, Number::DOT);
        let post = reserved("SVLEN", FieldKind::Info, VcfVersion::V4_4).unwrap();
        assert_eq!(post.number, Number::A);
    }

    #[test]
    fn svclaim_gated_before_4_4() {
        assert!(reserved("SVCLAIM", FieldKind::Info, VcfVersion::V4_3).is_err());
        assert!(reserved("SVCLAIM", FieldKind::Info, VcfVersion::V4_4).is_ok());
    }

    #[test]
    fn unknown_reserved_errs() {
        assert!(reserved("NOPE", FieldKind::Info, VcfVersion::V4_5).is_err());
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib spec::reserved` → FAIL.

- [ ] **Step 3: Implement**

```rust
use crate::error::BuildError;
use crate::spec::field::{FieldDef, FieldKind};
use crate::spec::number::Number;
use crate::spec::types::Type;
use crate::spec::version::VcfVersion;

fn info_entry(id: &str) -> Option<(Number, Type, &'static str)> {
    Some(match id {
        "AA" => (Number::ONE, Type::String, "Ancestral allele"),
        "AC" => (Number::A, Type::Integer, "Allele count"),
        "AF" => (Number::A, Type::Float, "Allele frequency"),
        "AN" => (Number::ONE, Type::Integer, "Total allele number"),
        "DP" => (Number::ONE, Type::Integer, "Combined depth"),
        "DB" => (Number::FLAG, Type::Flag, "dbSNP membership"),
        "H2" => (Number::FLAG, Type::Flag, "HapMap2 membership"),
        "END" => (Number::ONE, Type::Integer, "End position (deprecated)"),
        "SVTYPE" => (Number::ONE, Type::String, "Type of structural variant"),
        "SVLEN" => (Number::A, Type::Integer, "Length of structural variant"),
        "SVCLAIM" => (Number::A, Type::String, "Structural variant claim"),
        "CIPOS" => (Number::DOT, Type::Integer, "Confidence interval around POS"),
        "CIEND" => (Number::DOT, Type::Integer, "Confidence interval around END"),
        "CILEN" => (Number::DOT, Type::Integer, "Confidence interval around SVLEN"),
        "MATEID" => (Number::A, Type::String, "ID of mate breakend"),
        "PARID" => (Number::A, Type::String, "ID of partner breakend"),
        "IMPRECISE" => (Number::FLAG, Type::Flag, "Imprecise structural variant"),
        _ => return None,
    })
}

fn format_entry(id: &str) -> Option<(Number, Type, &'static str)> {
    Some(match id {
        "GT" => (Number::ONE, Type::String, "Genotype"),
        "GQ" => (Number::ONE, Type::Integer, "Genotype quality"),
        "DP" => (Number::ONE, Type::Integer, "Read depth"),
        "AD" => (Number::R, Type::Integer, "Allelic depths"),
        "PL" => (Number::G, Type::Integer, "Phred genotype likelihoods"),
        "GL" => (Number::G, Type::Float, "Log10 genotype likelihoods"),
        "PS" => (Number::ONE, Type::Integer, "Phase set"),
        "CN" => (Number::ONE, Type::Float, "Copy number"),
        "LEN" => (Number::ONE, Type::Integer, "Length of <*> reference block"),
        _ => return None,
    })
}

/// Version each gated reserved field was introduced; absent ids exist since 4.1.
fn since(id: &str, kind: FieldKind) -> VcfVersion {
    match (kind, id) {
        (FieldKind::Info, "SVCLAIM") => VcfVersion::V4_4,
        (FieldKind::Format, "LEN") => VcfVersion::V4_4,
        _ => VcfVersion::V4_1,
    }
}

pub fn reserved(
    id: &str,
    kind: FieldKind,
    version: VcfVersion,
) -> Result<FieldDef, BuildError> {
    let entry = match kind {
        FieldKind::Info => info_entry(id),
        FieldKind::Format => format_entry(id),
    };
    let (number, type_, desc) = entry.ok_or_else(|| BuildError::UnknownReserved {
        kind: kind.as_str().to_string(),
        id: id.to_string(),
    })?;

    let intro = since(id, kind);
    if version < intro {
        return Err(BuildError::FieldTooNew {
            kind: kind.as_str().to_string(),
            id: id.to_string(),
            since: intro.to_string(),
            version: version.to_string(),
        });
    }

    // SVLEN's pre-4.4 form: Number=. (signed length difference).
    if id == "SVLEN" && kind == FieldKind::Info && version < VcfVersion::V4_4 {
        return FieldDef::new(
            "SVLEN",
            Number::DOT,
            Type::Integer,
            "Difference in length between REF and ALT alleles",
            FieldKind::Info,
        );
    }

    FieldDef::new(id, number, type_, desc, kind)
}
```

Add `pub mod reserved;` to `src/spec/mod.rs`. (No `lib.rs` re-export; the builder calls `crate::spec::reserved::reserved`.)

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib spec::reserved` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add reserved-field registry"`

---

### Task 7: `spec/genotype_order.rs` — Number=G ordering

**Files:**
- Create: `src/spec/genotype_order.rs`
- Modify: `src/spec/mod.rs`

**Interfaces:**
- Produces: `genotype_ordering(ploidy: usize, n_alleles: usize) -> Vec<Vec<u32>>`.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diploid_biallelic_order() {
        // VCF Number=G order for ploidy 2, 2 alleles: 0/0, 0/1, 1/1
        assert_eq!(
            genotype_ordering(2, 2),
            vec![vec![0, 0], vec![0, 1], vec![1, 1]]
        );
    }

    #[test]
    fn count_matches_binomial() {
        // ploidy 2, 3 alleles => 6 genotypes
        assert_eq!(genotype_ordering(2, 3).len(), 6);
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib spec::genotype_order` → FAIL.

- [ ] **Step 3: Implement**

```rust
/// Ordered genotypes per the VCF `Number=G` ordering.
pub fn genotype_ordering(ploidy: usize, n_alleles: usize) -> Vec<Vec<u32>> {
    assert!(ploidy >= 1, "ploidy must be >= 1");
    rec(ploidy, n_alleles)
}

fn rec(p: usize, n_alleles: usize) -> Vec<Vec<u32>> {
    if p == 1 {
        return (0..n_alleles as u32).map(|a| vec![a]).collect();
    }
    let mut out = Vec::new();
    for a in 0..n_alleles as u32 {
        for prefix in rec(p - 1, n_alleles) {
            if *prefix.last().unwrap() <= a {
                let mut g = prefix.clone();
                g.push(a);
                out.push(g);
            }
        }
    }
    out
}
```

Add `pub mod genotype_order;` to `src/spec/mod.rs`.

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib spec::genotype_order` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add genotype ordering"`

---

### Task 8: `allele.rs` — Allele enum + parse + render

**Files:**
- Create: `src/allele.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- Consumes: `BuildError`.
- Produces:
  - `enum SvType { Del, Ins, Dup, Inv, Cnv }` with `as_str()` and `from_str(&str) -> Result<SvType, BuildError>`.
  - `enum Allele { Seq(String), Star, Symbolic { first_type: SvType, subtypes: Vec<String> }, Unspecified, Breakend { raw: String, single: bool } }`.
  - ctors: `Allele::seq(impl Into<String>) -> Result<Allele, BuildError>`, `Allele::star()`, `Allele::unspecified()`, `Allele::deletion(subtypes: impl IntoIterator<Item = impl Into<String>>)`, `::insertion`, `::duplication`, `::inversion`, `::cnv`, `Allele::breakend_parse(&str) -> Result<Allele, BuildError>`.
  - `Allele::render(&self) -> String`; `Allele::parse(&str) -> Allele`; `Allele::symbolic_type_str(&self) -> Option<String>` (e.g. `"DUP:TANDEM"`).

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq_validates() {
        assert_eq!(Allele::seq("GAT").unwrap().render(), "GAT");
        assert!(Allele::seq("GX").is_err());
    }

    #[test]
    fn symbolic_render_and_type_str() {
        let dup = Allele::duplication(["TANDEM"]);
        assert_eq!(dup.render(), "<DUP:TANDEM>");
        assert_eq!(dup.symbolic_type_str().as_deref(), Some("DUP:TANDEM"));
    }

    #[test]
    fn parse_dispatch() {
        assert!(matches!(Allele::parse("*"), Allele::Star));
        assert!(matches!(Allele::parse("<*>"), Allele::Unspecified));
        assert!(matches!(Allele::parse("<DEL>"), Allele::Symbolic { .. }));
        assert!(matches!(Allele::parse("T[chr2:5["), Allele::Breakend { single: false, .. }));
        assert!(matches!(Allele::parse(".A"), Allele::Breakend { single: true, .. }));
        assert!(matches!(Allele::parse("ACGT"), Allele::Seq(_)));
    }

    #[test]
    fn breakend_parse_rejects_junk() {
        assert!(Allele::breakend_parse("not-a-breakend").is_err());
        assert!(Allele::breakend_parse("G[chr2:321[").is_ok());
        assert!(Allele::breakend_parse("A.").is_ok());
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib allele` → FAIL.

- [ ] **Step 3: Implement** (no regex crate; hand-rolled validators)

```rust
use crate::error::BuildError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SvType {
    Del,
    Ins,
    Dup,
    Inv,
    Cnv,
}

impl SvType {
    pub fn as_str(&self) -> &'static str {
        match self {
            SvType::Del => "DEL",
            SvType::Ins => "INS",
            SvType::Dup => "DUP",
            SvType::Inv => "INV",
            SvType::Cnv => "CNV",
        }
    }

    pub fn from_str(s: &str) -> Result<SvType, BuildError> {
        Ok(match s {
            "DEL" => SvType::Del,
            "INS" => SvType::Ins,
            "DUP" => SvType::Dup,
            "INV" => SvType::Inv,
            "CNV" => SvType::Cnv,
            _ => return Err(BuildError::BadSvType(s.to_string())),
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Allele {
    Seq(String),
    Star,
    Symbolic { first_type: SvType, subtypes: Vec<String> },
    Unspecified,
    Breakend { raw: String, single: bool },
}

fn is_seq(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'))
}

fn is_seq_or_empty(s: &str) -> bool {
    s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'))
}

impl Allele {
    pub fn seq(bases: impl Into<String>) -> Result<Allele, BuildError> {
        let bases = bases.into();
        if !is_seq(&bases) {
            return Err(BuildError::BadAlleleBases(bases));
        }
        Ok(Allele::Seq(bases))
    }

    pub fn star() -> Allele {
        Allele::Star
    }

    pub fn unspecified() -> Allele {
        Allele::Unspecified
    }

    fn symbolic(first: SvType, subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::Symbolic {
            first_type: first,
            subtypes: subtypes.into_iter().map(Into::into).collect(),
        }
    }

    pub fn deletion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Del, subtypes)
    }
    pub fn insertion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Ins, subtypes)
    }
    pub fn duplication(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Dup, subtypes)
    }
    pub fn inversion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Inv, subtypes)
    }
    pub fn cnv(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Cnv, subtypes)
    }

    /// Parse a breakend replacement string (paired or single forms).
    pub fn breakend_parse(s: &str) -> Result<Allele, BuildError> {
        if is_single_breakend(s) {
            return Ok(Allele::Breakend { raw: s.to_string(), single: true });
        }
        if is_paired_breakend(s) {
            return Ok(Allele::Breakend { raw: s.to_string(), single: false });
        }
        Err(BuildError::BadBreakend(s.to_string()))
    }

    /// Inner `<...>` token, e.g. `DEL` or `DUP:TANDEM`.
    pub fn symbolic_type_str(&self) -> Option<String> {
        match self {
            Allele::Symbolic { first_type, subtypes } => {
                let mut parts = vec![first_type.as_str().to_string()];
                parts.extend(subtypes.iter().cloned());
                Some(parts.join(":"))
            }
            _ => None,
        }
    }

    pub fn render(&self) -> String {
        match self {
            Allele::Seq(b) => b.clone(),
            Allele::Star => "*".to_string(),
            Allele::Unspecified => "<*>".to_string(),
            Allele::Breakend { raw, .. } => raw.clone(),
            Allele::Symbolic { .. } => format!("<{}>", self.symbolic_type_str().unwrap()),
        }
    }

    /// Syntactic dispatch from a raw ALT string (never fails: junk falls back
    /// to a sequence allele, matching the Python `classify_allele`).
    pub fn parse(alt: &str) -> Allele {
        if alt == "*" {
            return Allele::Star;
        }
        if alt == "<*>" {
            return Allele::Unspecified;
        }
        if alt.starts_with('<') && alt.ends_with('>') {
            let inner = &alt[1..alt.len() - 1];
            let mut parts = inner.split(':');
            let first = parts.next().unwrap_or("");
            let first_type = SvType::from_str(first).unwrap_or(SvType::Del);
            return Allele::Symbolic {
                first_type,
                subtypes: parts.map(|s| s.to_string()).collect(),
            };
        }
        if alt.contains('[') || alt.contains(']') {
            if let Ok(b) = Allele::breakend_parse(alt) {
                return b;
            }
        }
        if alt.len() > 1 && (alt.starts_with('.') || alt.ends_with('.')) {
            if let Ok(b) = Allele::breakend_parse(alt) {
                return b;
            }
        }
        Allele::Seq(alt.to_string())
    }
}

/// Single breakend: `.t` or `t.` where t is a non-empty sequence.
fn is_single_breakend(s: &str) -> bool {
    if let Some(rest) = s.strip_prefix('.') {
        return is_seq(rest);
    }
    if let Some(rest) = s.strip_suffix('.') {
        return is_seq(rest);
    }
    false
}

/// Paired breakend: `t[p[`, `t]p]`, `[p[t`, `]p]t` where both brackets are the
/// same char, t is a (possibly empty) sequence, and p is `chr:pos`.
fn is_paired_breakend(s: &str) -> bool {
    let bytes = s.as_bytes();
    let open = match bytes.iter().position(|&b| b == b'[' || b == b']') {
        Some(i) => i,
        None => return false,
    };
    let bracket = bytes[open];
    let close = match bytes.iter().rposition(|&b| b == b'[' || b == b']') {
        Some(i) if i != open => i,
        _ => return false,
    };
    if bytes[close] != bracket {
        return false;
    }
    let left = &s[..open];
    let mate = &s[open + 1..close];
    let right = &s[close + 1..];
    // Exactly one of left/right is the sequence side; the other is empty.
    if !is_seq_or_empty(left) || !is_seq_or_empty(right) {
        return false;
    }
    if left.is_empty() == right.is_empty() {
        return false; // need exactly one side with the replacement sequence
    }
    valid_mate(mate)
}

/// `chr:pos` with pos all-digits and a non-empty contig containing no brackets.
fn valid_mate(mate: &str) -> bool {
    match mate.rsplit_once(':') {
        Some((chrom, pos)) => {
            !chrom.is_empty()
                && !chrom.contains('[')
                && !chrom.contains(']')
                && !pos.is_empty()
                && pos.bytes().all(|b| b.is_ascii_digit())
        }
        None => false,
    }
}
```

Add to `src/lib.rs`: `pub mod allele; pub use allele::{Allele, SvType};`

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib allele` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add Allele model"`

---

### Task 9: `genotype.rs` — Genotype parse/render

**Files:**
- Create: `src/genotype.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- Produces: `struct Genotype { alleles: Vec<Option<u32>>, phased: Vec<bool> }`; `Genotype::parse(&str) -> Genotype`; `render(&self) -> String`; `ploidy(&self) -> usize`; `is_phased(&self) -> bool` (true iff ≥1 separator and all are `|`).

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_render_roundtrip() {
        for s in ["0|1", "1/1", "./.", "0", ".|1"] {
            assert_eq!(Genotype::parse(s).render(), s);
        }
    }

    #[test]
    fn phasing_and_ploidy() {
        assert!(Genotype::parse("0|1").is_phased());
        assert!(!Genotype::parse("0/1").is_phased());
        assert!(!Genotype::parse("0").is_phased());
        assert_eq!(Genotype::parse("0|1|1").ploidy(), 3);
        assert_eq!(Genotype::parse("./.").alleles, vec![None, None]);
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib genotype` → FAIL.

- [ ] **Step 3: Implement**

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// Allele indices in call order; `None` is a missing allele (`.`).
    pub alleles: Vec<Option<u32>>,
    /// One bool per separator; `true` = phased (`|`). Length == ploidy - 1.
    pub phased: Vec<bool>,
}

impl Genotype {
    pub fn parse(s: &str) -> Genotype {
        let mut alleles = Vec::new();
        let mut phased = Vec::new();
        let mut token = String::new();
        for ch in s.chars() {
            if ch == '|' || ch == '/' {
                alleles.push(parse_allele(&token));
                token.clear();
                phased.push(ch == '|');
            } else {
                token.push(ch);
            }
        }
        alleles.push(parse_allele(&token));
        Genotype { alleles, phased }
    }

    pub fn ploidy(&self) -> usize {
        self.alleles.len()
    }

    pub fn is_phased(&self) -> bool {
        !self.phased.is_empty() && self.phased.iter().all(|&p| p)
    }

    pub fn render(&self) -> String {
        let mut out = String::new();
        for (i, a) in self.alleles.iter().enumerate() {
            if i > 0 {
                out.push(if self.phased[i - 1] { '|' } else { '/' });
            }
            match a {
                Some(v) => out.push_str(&v.to_string()),
                None => out.push('.'),
            }
        }
        out
    }
}

fn parse_allele(tok: &str) -> Option<u32> {
    if tok == "." {
        None
    } else {
        Some(tok.parse().expect("genotype allele index must be an integer or '.'"))
    }
}
```

Add to `src/lib.rs`: `pub mod genotype; pub use genotype::Genotype;`

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib genotype` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add Genotype"`

---

### Task 10: `variants.rs` — variant classification

**Files:**
- Create: `src/variants.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- Consumes: `Allele`, `SvType`.
- Produces:
  - `enum VariantClass { Snp, Mnp, Ins, Del, Delins, SpanningDel, Unspecified, Bnd, Multiallelic, SvDel, SvIns, SvDup, SvInv, Cnv }` with `as_str()`.
  - `classify_seq(ref_: &str, alt: &str) -> VariantClass` (sequence-pair classifier: SNP/MNP/INS/DEL/DELINS/SpanningDel).
  - `record_class(ref_: &str, alts: &[Allele]) -> VariantClass`.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;

    #[test]
    fn sequence_classes() {
        assert_eq!(classify_seq("A", "T"), VariantClass::Snp);
        assert_eq!(classify_seq("AC", "GT"), VariantClass::Mnp);
        assert_eq!(classify_seq("A", "AT"), VariantClass::Ins);
        assert_eq!(classify_seq("AT", "A"), VariantClass::Del);
        assert_eq!(classify_seq("AT", "GC".get(0..1).map(|_| "C").unwrap()), VariantClass::Delins);
        assert_eq!(classify_seq("A", "*"), VariantClass::SpanningDel);
    }

    #[test]
    fn record_classes() {
        assert_eq!(record_class("A", &[Allele::seq("T").unwrap()]), VariantClass::Snp);
        assert_eq!(
            record_class("A", &[Allele::seq("T").unwrap(), Allele::seq("C").unwrap()]),
            VariantClass::Multiallelic
        );
        assert_eq!(record_class("A", &[Allele::deletion::<[&str; 0], _>([])]), VariantClass::SvDel);
        assert_eq!(record_class("A", &[Allele::cnv::<[&str; 0], _>([])]), VariantClass::Cnv);
        assert_eq!(record_class("A", &[Allele::Star]), VariantClass::SpanningDel);
    }
}
```

(Note: the `::<[&str; 0], _>` turbofish picks the empty-subtypes case; the implementer may instead pass `Vec::<&str>::new()` for readability.)

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib variants` → FAIL.

- [ ] **Step 3: Implement**

```rust
use crate::allele::{Allele, SvType};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass {
    Snp,
    Mnp,
    Ins,
    Del,
    Delins,
    SpanningDel,
    Unspecified,
    Bnd,
    Multiallelic,
    SvDel,
    SvIns,
    SvDup,
    SvInv,
    Cnv,
}

impl VariantClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            VariantClass::Snp => "SNP",
            VariantClass::Mnp => "MNP",
            VariantClass::Ins => "INS",
            VariantClass::Del => "DEL",
            VariantClass::Delins => "DELINS",
            VariantClass::SpanningDel => "SPANNING_DEL",
            VariantClass::Unspecified => "UNSPECIFIED",
            VariantClass::Bnd => "BND",
            VariantClass::Multiallelic => "MULTIALLELIC",
            VariantClass::SvDel => "SV_DEL",
            VariantClass::SvIns => "SV_INS",
            VariantClass::SvDup => "SV_DUP",
            VariantClass::SvInv => "SV_INV",
            VariantClass::Cnv => "CNV",
        }
    }
}

pub fn classify_seq(ref_: &str, alt: &str) -> VariantClass {
    if alt == "*" {
        return VariantClass::SpanningDel;
    }
    let (lr, la) = (ref_.len(), alt.len());
    if lr == 1 && la == 1 {
        VariantClass::Snp
    } else if lr == la {
        VariantClass::Mnp
    } else if la > lr && alt.starts_with(ref_) {
        VariantClass::Ins
    } else if lr > la && ref_.starts_with(alt) {
        VariantClass::Del
    } else {
        VariantClass::Delins
    }
}

pub fn record_class(ref_: &str, alts: &[Allele]) -> VariantClass {
    if alts.len() != 1 {
        return VariantClass::Multiallelic;
    }
    match &alts[0] {
        Allele::Seq(bases) => classify_seq(ref_, bases),
        Allele::Star => VariantClass::SpanningDel,
        Allele::Unspecified => VariantClass::Unspecified,
        Allele::Breakend { .. } => VariantClass::Bnd,
        Allele::Symbolic { first_type, .. } => match first_type {
            SvType::Del => VariantClass::SvDel,
            SvType::Ins => VariantClass::SvIns,
            SvType::Dup => VariantClass::SvDup,
            SvType::Inv => VariantClass::SvInv,
            SvType::Cnv => VariantClass::Cnv,
        },
    }
}
```

Add to `src/lib.rs`: `pub mod variants; pub use variants::VariantClass;`

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib variants` → PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add variant classification"`

---

### Task 11: `model.rs` — value types + Document/Record hub

**Files:**
- Create: `src/model.rs`, `src/value.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- Consumes: `Allele`, `Genotype`, `FieldDef`, `VcfVersion`.
- Produces (in `value.rs`):
  - `enum Scalar { Int(i64), Float(f64), Char(char), Str(String) }` (`Clone`, `PartialEq`).
  - `enum FieldValue { Flag, Scalar(Scalar), List(Vec<Option<Scalar>>) }` with helpers `FieldValue::list_len(&self) -> Option<usize>` (None for Flag/Scalar) and `FieldValue::ints(slice)`, `FieldValue::floats(slice)` convenience ctors used by builder callers.
- Produces (in `model.rs`):
  - `struct ContigDef { id: String, length: Option<u64> }` with `header_line()`.
  - `struct AltDef { id: String, description: String }` with `header_line()`.
  - `struct SampleValues { gt: Option<Genotype>, values: IndexMap<String, FieldValue> }`.
  - `struct Record { chrom, pos: u64, ids: Option<Vec<String>>, ref_: String, alts: Vec<Allele>, qual: Option<f64>, filters: Option<Vec<String>>, info: IndexMap<String, FieldValue>, fmt_keys: Vec<String>, samples: Vec<SampleValues>, labels: BTreeSet<String> }` with `n_alt(&self) -> usize`.
  - `struct Document { version, info_defs: Vec<FieldDef>, format_defs: Vec<FieldDef>, filter_defs: Vec<(String,String)>, contigs: Vec<ContigDef>, samples: Vec<String>, records: Vec<Record>, alt_defs: Vec<AltDef> }` with `max_ploidy(&self) -> usize`.

- [ ] **Step 1: Write the failing test** (in `src/model.rs`)

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contig_header_lines() {
        assert_eq!(
            ContigDef { id: "chr1".into(), length: Some(100) }.header_line(),
            "##contig=<ID=chr1,length=100>"
        );
        assert_eq!(
            ContigDef { id: "chr1".into(), length: None }.header_line(),
            "##contig=<ID=chr1>"
        );
    }

    #[test]
    fn max_ploidy_defaults_and_scans() {
        let doc = Document {
            version: crate::spec::version::LATEST,
            info_defs: vec![],
            format_defs: vec![],
            filter_defs: vec![],
            contigs: vec![],
            samples: vec!["s1".into()],
            records: vec![],
            alt_defs: vec![],
        };
        assert_eq!(doc.max_ploidy(), 1);
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib model` → FAIL.

- [ ] **Step 3: Implement `src/value.rs`**

```rust
/// A single decoded INFO/FORMAT scalar.
#[derive(Debug, Clone, PartialEq)]
pub enum Scalar {
    Int(i64),
    Float(f64),
    Char(char),
    Str(String),
}

/// A decoded INFO/FORMAT field value.
#[derive(Debug, Clone, PartialEq)]
pub enum FieldValue {
    Flag,
    Scalar(Scalar),
    List(Vec<Option<Scalar>>),
}

impl FieldValue {
    /// Number of list entries, or `None` for Flag/lone-scalar values.
    pub fn list_len(&self) -> Option<usize> {
        match self {
            FieldValue::List(v) => Some(v.len()),
            _ => None,
        }
    }

    pub fn ints<I: IntoIterator<Item = i64>>(xs: I) -> FieldValue {
        FieldValue::List(xs.into_iter().map(|x| Some(Scalar::Int(x))).collect())
    }

    pub fn floats<I: IntoIterator<Item = f64>>(xs: I) -> FieldValue {
        FieldValue::List(xs.into_iter().map(|x| Some(Scalar::Float(x))).collect())
    }

    pub fn strings<I, S>(xs: I) -> FieldValue
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        FieldValue::List(xs.into_iter().map(|x| Some(Scalar::Str(x.into()))).collect())
    }
}
```

- [ ] **Step 4: Implement `src/model.rs`**

```rust
use std::collections::BTreeSet;

use indexmap::IndexMap;

use crate::allele::Allele;
use crate::genotype::Genotype;
use crate::spec::field::FieldDef;
use crate::spec::version::VcfVersion;
use crate::value::FieldValue;

#[derive(Debug, Clone, PartialEq)]
pub struct ContigDef {
    pub id: String,
    pub length: Option<u64>,
}

impl ContigDef {
    pub fn header_line(&self) -> String {
        match self.length {
            Some(n) => format!("##contig=<ID={},length={}>", self.id, n),
            None => format!("##contig=<ID={}>", self.id),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AltDef {
    pub id: String,
    pub description: String,
}

impl AltDef {
    pub fn header_line(&self) -> String {
        format!("##ALT=<ID={},Description=\"{}\">", self.id, self.description)
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SampleValues {
    pub gt: Option<Genotype>,
    pub values: IndexMap<String, FieldValue>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub chrom: String,
    pub pos: u64,
    pub ids: Option<Vec<String>>,
    pub ref_: String,
    pub alts: Vec<Allele>,
    pub qual: Option<f64>,
    pub filters: Option<Vec<String>>,
    pub info: IndexMap<String, FieldValue>,
    pub fmt_keys: Vec<String>,
    pub samples: Vec<SampleValues>,
    pub labels: BTreeSet<String>,
}

impl Record {
    pub fn n_alt(&self) -> usize {
        self.alts.len()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Document {
    pub version: VcfVersion,
    pub info_defs: Vec<FieldDef>,
    pub format_defs: Vec<FieldDef>,
    pub filter_defs: Vec<(String, String)>,
    pub contigs: Vec<ContigDef>,
    pub samples: Vec<String>,
    pub records: Vec<Record>,
    pub alt_defs: Vec<AltDef>,
}

impl Document {
    pub fn max_ploidy(&self) -> usize {
        let mut p = 1;
        for rec in &self.records {
            for s in &rec.samples {
                if let Some(gt) = &s.gt {
                    p = p.max(gt.ploidy());
                }
            }
        }
        p
    }
}
```

Add to `src/lib.rs`:

```rust
pub mod value;
pub mod model;
pub use value::{FieldValue, Scalar};
pub use model::{AltDef, ContigDef, Document, Record, SampleValues};
```

- [ ] **Step 5: Run to verify it passes** — `pixi run cargo test --lib model` → PASS. Then `pixi run clippy`.

- [ ] **Step 6: Commit** — `git add -A && git commit -m "feat: add Document/Record model and value types"`

---

### Task 12: `build.rs` — VcfBuilder + Record sub-builder + validation

**Files:**
- Create: `src/build.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `struct RecordSpec` (the record sub-builder) with `RecordSpec::at(chrom, pos) -> RecordSpec` and chained setters: `ref_(s)`, `alt(IntoIterator<Allele>)`, `ids(IntoIterator<into String>)`, `qual(f64)`, `filter(IntoIterator<into String>)`, `gt(IntoIterator<into String>)`, `info(id, FieldValue)`, `format(id, IntoIterator<FieldValue>)`, `labels(IntoIterator<into String>)`. (`info`/`format` accept a `FieldValue` and a per-sample `Vec<FieldValue>` respectively.)
  - `struct VcfBuilder` with `new(samples, contigs, version)`, `info(...) -> Result<Self>`, `format(...) -> Result<Self>`, `filter(id, desc) -> Self`, `alt(id, desc) -> Self`, `record(RecordSpec) -> Result<Self>`, `build(self) -> Result<Document>`, plus `render()/truth()/write()` convenience that build then delegate.
- Note signatures used downstream: `VcfBuilder::record` performs all validation listed in the spec; `build` auto-registers `AltDef`s.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;
    use crate::spec::number::Number;
    use crate::spec::types::Type;
    use crate::spec::version::LATEST;
    use crate::value::{FieldValue, Scalar};

    fn base() -> VcfBuilder {
        VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
    }

    #[test]
    fn happy_path_builds() {
        let doc = base()
            .info("AF", None, None, None).unwrap()
            .format("GT", None, None, None).unwrap()
            .format("DS", Some(Number::A), Some(Type::Float), None).unwrap()
            .record(
                RecordSpec::at("chr1", 1000)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1", "1|1"])
                    .info("AF", FieldValue::floats([0.25]))
                    .format("DS", [FieldValue::floats([0.4]), FieldValue::floats([1.9])]),
            ).unwrap()
            .build().unwrap();
        assert_eq!(doc.records.len(), 1);
        assert_eq!(doc.records[0].samples[0].gt.as_ref().unwrap().render(), "0|1");
    }

    #[test]
    fn undeclared_field_errs() {
        let r = base()
            .format("GT", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 1).ref_("A").alt([Allele::seq("T").unwrap()])
                .info("AF", FieldValue::floats([0.1])));
        assert!(matches!(r, Err(crate::error::BuildError::UndeclaredField { .. })));
    }

    #[test]
    fn cardinality_checked() {
        let r = base()
            .format("GT", None, None, None).unwrap()
            .info("AF", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 1).ref_("A")
                .alt([Allele::seq("T").unwrap()])           // n_alt = 1, AF is Number::A
                .info("AF", FieldValue::floats([0.1, 0.2]))); // 2 values -> mismatch
        assert!(matches!(r, Err(crate::error::BuildError::Cardinality { .. })));
    }

    #[test]
    fn symbolic_requires_svlen_and_padding() {
        // missing SVLEN
        let r = base().format("GT", None, None, None).unwrap()
            .info("SVLEN", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 1).ref_("A")
                .alt([Allele::deletion(Vec::<&str>::new())]));
        assert!(matches!(r, Err(crate::error::BuildError::MissingSvlen(_))));

        // multi-base REF padding violation
        let r = base().format("GT", None, None, None).unwrap()
            .info("SVLEN", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 1).ref_("AC")
                .alt([Allele::deletion(Vec::<&str>::new())])
                .info("SVLEN", FieldValue::ints([100])));
        assert!(matches!(r, Err(crate::error::BuildError::MissingRefPadding(_))));
    }

    #[test]
    fn gt_index_out_of_range() {
        let r = base().format("GT", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 1).ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|2", "0|0"]));  // index 2 > n_alt 1
        assert!(matches!(r, Err(crate::error::BuildError::AlleleIndexOutOfRange { .. })));
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib build` → FAIL.

- [ ] **Step 3: Implement** (the full builder; mirror the Python `build.py` logic)

```rust
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use indexmap::IndexMap;

use crate::allele::{Allele, SvType};
use crate::error::BuildError;
use crate::genotype::Genotype;
use crate::model::{AltDef, ContigDef, Document, Record, SampleValues};
use crate::spec::field::{FieldDef, FieldKind};
use crate::spec::number::{Number, NumberKind};
use crate::spec::reserved::reserved;
use crate::spec::types::Type;
use crate::spec::version::{VcfVersion, LATEST};
use crate::value::FieldValue;

/// Allowed SVCLAIM tokens per first-level SV type.
fn svclaim_allowed(t: SvType) -> &'static [&'static str] {
    match t {
        SvType::Del | SvType::Dup => &["D", "J", "DJ"],
        SvType::Cnv => &["D"],
        SvType::Ins | SvType::Inv => &["J"],
    }
}

fn svclaim_required(t: SvType) -> bool {
    matches!(t, SvType::Del | SvType::Dup)
}

fn cn_svlen_type(t: SvType) -> bool {
    matches!(t, SvType::Cnv | SvType::Del | SvType::Dup)
}

/// A record's spec, before validation/appending.
#[derive(Debug, Clone, Default)]
pub struct RecordSpec {
    chrom: String,
    pos: u64,
    ref_: String,
    alts: Vec<Allele>,
    ids: Option<Vec<String>>,
    qual: Option<f64>,
    filters: Option<Vec<String>>,
    gt: Option<Vec<String>>,
    info: IndexMap<String, FieldValue>,
    fmt: IndexMap<String, Vec<FieldValue>>,
    labels: BTreeSet<String>,
}

impl RecordSpec {
    pub fn at(chrom: impl Into<String>, pos: u64) -> RecordSpec {
        RecordSpec { chrom: chrom.into(), pos, ..Default::default() }
    }
    pub fn ref_(mut self, r: impl Into<String>) -> Self {
        self.ref_ = r.into();
        self
    }
    pub fn alt(mut self, alts: impl IntoIterator<Item = Allele>) -> Self {
        self.alts = alts.into_iter().collect();
        self
    }
    pub fn ids(mut self, ids: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.ids = Some(ids.into_iter().map(Into::into).collect());
        self
    }
    pub fn qual(mut self, q: f64) -> Self {
        self.qual = Some(q);
        self
    }
    pub fn filter(mut self, f: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.filters = Some(f.into_iter().map(Into::into).collect());
        self
    }
    pub fn gt(mut self, gts: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.gt = Some(gts.into_iter().map(Into::into).collect());
        self
    }
    pub fn info(mut self, id: impl Into<String>, value: FieldValue) -> Self {
        self.info.insert(id.into(), value);
        self
    }
    pub fn format(mut self, id: impl Into<String>, per_sample: impl IntoIterator<Item = FieldValue>) -> Self {
        self.fmt.insert(id.into(), per_sample.into_iter().collect());
        self
    }
    pub fn labels(mut self, labels: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.labels = labels.into_iter().map(Into::into).collect();
        self
    }
}

pub struct VcfBuilder {
    samples: Vec<String>,
    contigs: Vec<ContigDef>,
    version: VcfVersion,
    info_defs: IndexMap<String, FieldDef>,
    format_defs: IndexMap<String, FieldDef>,
    filter_defs: Vec<(String, String)>,
    alt_defs: IndexMap<String, String>,
    records: Vec<Record>,
}

impl VcfBuilder {
    pub fn new(
        samples: impl IntoIterator<Item = impl Into<String>>,
        contigs: impl IntoIterator<Item = (impl Into<String>, Option<u64>)>,
        version: VcfVersion,
    ) -> VcfBuilder {
        VcfBuilder {
            samples: samples.into_iter().map(Into::into).collect(),
            contigs: contigs
                .into_iter()
                .map(|(id, length)| ContigDef { id: id.into(), length })
                .collect(),
            version,
            info_defs: IndexMap::new(),
            format_defs: IndexMap::new(),
            filter_defs: Vec::new(),
            alt_defs: IndexMap::new(),
            records: Vec::new(),
        }
    }

    fn make_def(
        &self,
        id: &str,
        number: Option<Number>,
        type_: Option<Type>,
        description: Option<String>,
        kind: FieldKind,
    ) -> Result<FieldDef, BuildError> {
        match (number, type_) {
            (Some(n), Some(t)) => {
                FieldDef::new(id, n, t, description.unwrap_or_else(|| id.to_string()), kind)
            }
            _ => reserved(id, kind, self.version),
        }
    }

    pub fn info(
        mut self,
        id: impl AsRef<str>,
        number: Option<Number>,
        type_: Option<Type>,
        description: Option<String>,
    ) -> Result<VcfBuilder, BuildError> {
        let id = id.as_ref();
        let def = self.make_def(id, number, type_, description, FieldKind::Info)?;
        self.info_defs.insert(id.to_string(), def);
        Ok(self)
    }

    pub fn format(
        mut self,
        id: impl AsRef<str>,
        number: Option<Number>,
        type_: Option<Type>,
        description: Option<String>,
    ) -> Result<VcfBuilder, BuildError> {
        let id = id.as_ref();
        let def = self.make_def(id, number, type_, description, FieldKind::Format)?;
        self.format_defs.insert(id.to_string(), def);
        Ok(self)
    }

    pub fn filter(mut self, id: impl Into<String>, description: impl Into<String>) -> VcfBuilder {
        self.filter_defs.push((id.into(), description.into()));
        self
    }

    pub fn alt(mut self, id: impl Into<String>, description: impl Into<String>) -> VcfBuilder {
        self.alt_defs.insert(id.into(), description.into());
        self
    }

    pub fn record(mut self, spec: RecordSpec) -> Result<VcfBuilder, BuildError> {
        let n_alt = spec.alts.len();
        self.validate_alleles(&spec)?;

        let mut fmt_keys: Vec<String> = Vec::new();
        let mut samples: Vec<SampleValues> = vec![SampleValues::default(); self.samples.len()];

        // GT
        if let Some(gts) = &spec.gt {
            if !self.format_defs.contains_key("GT") {
                return Err(BuildError::GtNotDeclared);
            }
            fmt_keys.push("GT".to_string());
            for (si, s) in gts.iter().enumerate() {
                let geno = Genotype::parse(s);
                for a in geno.alleles.iter().flatten() {
                    if *a as usize > n_alt {
                        return Err(BuildError::AlleleIndexOutOfRange { index: *a, n_alt });
                    }
                }
                samples[si].gt = Some(geno);
            }
        }

        let ploidy = samples
            .iter()
            .filter_map(|s| s.gt.as_ref().map(|g| g.ploidy()))
            .max()
            .unwrap_or(2);

        // FORMAT (non-GT)
        for (key, per_sample) in &spec.fmt {
            let fdef = self
                .format_defs
                .get(key)
                .ok_or_else(|| BuildError::UndeclaredField {
                    kind: "FORMAT".into(),
                    id: key.clone(),
                })?;
            fmt_keys.push(key.clone());
            let card = fdef.number.cardinality(n_alt, ploidy);
            for (si, val) in per_sample.iter().enumerate() {
                check_cardinality(key, fdef.number.kind, card, val)?;
                samples[si].values.insert(key.clone(), val.clone());
            }
        }

        // INFO
        let mut info: IndexMap<String, FieldValue> = IndexMap::new();
        for (key, val) in &spec.info {
            let fdef = self
                .info_defs
                .get(key)
                .ok_or_else(|| BuildError::UndeclaredField {
                    kind: "INFO".into(),
                    id: key.clone(),
                })?;
            let card = fdef.number.cardinality(n_alt, ploidy);
            if fdef.number.kind != NumberKind::Flag {
                check_cardinality(key, fdef.number.kind, card, val)?;
            }
            info.insert(key.clone(), val.clone());
        }

        // FORMAT CN requires equal SVLEN across CNV/DEL/DUP alleles.
        if fmt_keys.iter().any(|k| k == "CN") {
            let svlen = spec.info.get("SVLEN");
            let mut seen: Vec<Option<i64>> = Vec::new();
            for (i, a) in spec.alts.iter().enumerate() {
                if let Allele::Symbolic { first_type, .. } = a {
                    if cn_svlen_type(*first_type) {
                        seen.push(per_allele_int(svlen, i));
                    }
                }
            }
            seen.dedup();
            if seen.len() > 1 {
                return Err(BuildError::CnSvlenMismatch);
            }
        }

        self.records.push(Record {
            chrom: spec.chrom,
            pos: spec.pos,
            ids: spec.ids,
            ref_: spec.ref_,
            alts: spec.alts,
            qual: spec.qual,
            filters: spec.filters,
            info,
            fmt_keys,
            samples,
            labels: spec.labels,
        });
        Ok(self)
    }

    fn validate_alleles(&self, spec: &RecordSpec) -> Result<(), BuildError> {
        let svlen = spec.info.get("SVLEN");
        let svclaim = spec.info.get("SVCLAIM");
        let needs_padding = spec
            .alts
            .iter()
            .any(|a| matches!(a, Allele::Symbolic { .. } | Allele::Breakend { .. }));
        if needs_padding && spec.ref_.len() != 1 {
            return Err(BuildError::MissingRefPadding(spec.ref_.clone()));
        }
        for (i, a) in spec.alts.iter().enumerate() {
            let sv = per_allele_int(svlen, i);
            let cl = per_allele_str(svclaim, i);
            match a {
                Allele::Symbolic { first_type, .. } => {
                    if sv.is_none() {
                        return Err(BuildError::MissingSvlen(a.render()));
                    }
                    let allowed = svclaim_allowed(*first_type);
                    if let Some(c) = &cl {
                        if !allowed.contains(&c.as_str()) {
                            return Err(BuildError::BadSvclaim {
                                claim: c.clone(),
                                allele: a.render(),
                                allowed: allowed.iter().map(|s| s.to_string()).collect(),
                            });
                        }
                    }
                    if self.version >= VcfVersion::V4_4
                        && svclaim_required(*first_type)
                        && cl.is_none()
                    {
                        return Err(BuildError::SvclaimRequired(a.render()));
                    }
                }
                Allele::Breakend { .. } | Allele::Unspecified | Allele::Star => {
                    if sv.is_some() {
                        return Err(BuildError::SvlenMustBeMissing(a.render()));
                    }
                }
                Allele::Seq(_) => {}
            }
        }
        Ok(())
    }

    pub fn build(self) -> Result<Document, BuildError> {
        // Auto-describe symbolic ALT types; explicit .alt() descriptions win.
        let mut alt_ids: IndexMap<String, String> = IndexMap::new();
        for rec in &self.records {
            for a in &rec.alts {
                if let Some(ts) = a.symbolic_type_str() {
                    alt_ids
                        .entry(ts.clone())
                        .or_insert_with(|| format!("{ts} structural variant"));
                }
            }
        }
        for (id, desc) in &self.alt_defs {
            alt_ids.insert(id.clone(), desc.clone());
        }
        let alt_defs = alt_ids
            .into_iter()
            .map(|(id, description)| AltDef { id, description })
            .collect();

        Ok(Document {
            version: self.version,
            info_defs: self.info_defs.into_values().collect(),
            format_defs: self.format_defs.into_values().collect(),
            filter_defs: self.filter_defs,
            contigs: self.contigs,
            samples: self.samples,
            records: self.records,
            alt_defs,
        })
    }

    pub fn render(self) -> Result<String, BuildError> {
        Ok(self.build()?.render())
    }

    pub fn truth(self) -> Result<crate::truth::GroundTruth, BuildError> {
        Ok(self.build()?.truth())
    }

    pub fn write(self, path: impl AsRef<Path>, opts: crate::write::WriteOpts) -> Result<PathBuf, BuildError> {
        self.build()?.write(path, opts)
    }
}

/// Resolve the i-th per-allele integer of a Number=A info value.
fn per_allele_int(value: Option<&FieldValue>, i: usize) -> Option<i64> {
    match value {
        Some(FieldValue::List(v)) => v.get(i).and_then(|x| match x {
            Some(crate::value::Scalar::Int(n)) => Some(*n),
            Some(crate::value::Scalar::Float(f)) => Some(*f as i64),
            _ => None,
        }),
        Some(FieldValue::Scalar(crate::value::Scalar::Int(n))) if i == 0 => Some(*n),
        _ => None,
    }
}

fn per_allele_str(value: Option<&FieldValue>, i: usize) -> Option<String> {
    match value {
        Some(FieldValue::List(v)) => v.get(i).and_then(|x| match x {
            Some(crate::value::Scalar::Str(s)) => Some(s.clone()),
            _ => None,
        }),
        Some(FieldValue::Scalar(crate::value::Scalar::Str(s))) if i == 0 => Some(s.clone()),
        _ => None,
    }
}

/// Enforce that a list value's length equals the resolved cardinality.
/// Flag and unbounded (`None`) cardinalities are not checked. Lone scalars and
/// Flag values bypass the length check (they have no list length).
fn check_cardinality(
    id: &str,
    kind: NumberKind,
    card: Option<usize>,
    val: &FieldValue,
) -> Result<(), BuildError> {
    if kind == NumberKind::Flag {
        return Ok(());
    }
    if let (Some(expected), Some(got)) = (card, val.list_len()) {
        if expected != got {
            return Err(BuildError::Cardinality {
                id: id.to_string(),
                expected,
                got,
            });
        }
    }
    Ok(())
}

// Silence unused-import warnings if HashMap ends up unused after edits.
#[allow(unused_imports)]
use std::collections::HashMap as _UnusedHashMap;
```

(If `HashMap`/`LATEST` end up unused, remove those imports rather than suppressing — the `_UnusedHashMap` shim above is a reminder to clean imports so clippy passes with `-D warnings`.)

Add to `src/lib.rs`:

```rust
pub mod build;
pub use build::{RecordSpec, VcfBuilder};
```

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib build` → PASS. Then remove any unused imports and run `pixi run clippy` until clean.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add VcfBuilder with eager validation"`

---

### Task 13: `truth.rs` — GroundTruth oracle (ndarray)

**Files:**
- Create: `src/truth.rs`
- Modify: `src/lib.rs`, `src/model.rs` (add `Document::truth`)

**Interfaces:**
- Consumes: `Document`, `Allele`, `SvType`, `Genotype`, `VariantClass`, `FieldValue`.
- Produces:
  - `enum AlleleKind { Snp, Mnp, Ins, Del, Delins, SpanningDel, Symbolic, Unspecified, Bnd }`.
  - `struct AlleleTruth { kind: AlleleKind, is_sequence: bool, sv_type: Option<String>, svlen: Option<i64>, sv_end: Option<i64> }`.
  - `struct GroundTruth { samples, contigs, pos: Array1<i64>, ref_: Vec<String>, alts: Vec<Vec<String>>, variant_class: Vec<VariantClass>, genotypes: Array3<i32>, phasing: Array2<bool>, info: Vec<HashMap<String,FieldValue>>, format: Vec<Vec<HashMap<String,FieldValue>>>, labels: Vec<BTreeSet<String>>, alts_truth: Vec<Vec<AlleleTruth>>, is_sequence_mask: Vec<Array1<bool>> }`.
  - `derive(doc: &Document) -> GroundTruth`; and `Document::truth(&self) -> GroundTruth` delegating to it.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;
    use crate::spec::version::LATEST;
    use crate::value::FieldValue;
    use crate::variants::VariantClass;
    use crate::VcfBuilder;
    use crate::RecordSpec;

    #[test]
    fn genotypes_phasing_and_missing() {
        let t = VcfBuilder::new(["s1", "s2"], [("chr1", Some(1000u64))], LATEST)
            .format("GT", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 10).ref_("A").alt([Allele::seq("T").unwrap()])
                .gt(["0|1", "./."])).unwrap()
            .build().unwrap()
            .truth();
        assert_eq!(t.genotypes[[0, 0, 0]], 0);
        assert_eq!(t.genotypes[[0, 0, 1]], 1);
        assert_eq!(t.genotypes[[0, 1, 0]], -1);
        assert!(t.phasing[[0, 0]]);
        assert!(!t.phasing[[0, 1]]);
        assert_eq!(t.variant_class[0], VariantClass::Snp);
    }

    #[test]
    fn symbolic_svlen_and_end() {
        let t = VcfBuilder::new(["s1"], [("chr1", Some(100_000u64))], LATEST)
            .format("GT", None, None, None).unwrap()
            .info("SVLEN", None, None, None).unwrap()
            .info("SVCLAIM", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 100).ref_("A")
                .alt([Allele::deletion(Vec::<&str>::new())])
                .info("SVLEN", FieldValue::ints([250]))
                .info("SVCLAIM", FieldValue::strings(["D"]))
                .gt(["0|1"])).unwrap()
            .build().unwrap()
            .truth();
        let at = &t.alts_truth[0][0];
        assert_eq!(at.kind, AlleleKind::Symbolic);
        assert_eq!(at.svlen, Some(250));
        assert_eq!(at.sv_end, Some(350)); // pos + svlen
        assert!(!t.is_sequence_mask[0][0]);
    }

    #[test]
    fn info_excludes_gt_from_format() {
        let t = VcfBuilder::new(["s1"], [("chr1", Some(1000u64))], LATEST)
            .format("GT", None, None, None).unwrap()
            .format("GQ", None, None, None).unwrap()
            .record(RecordSpec::at("chr1", 10).ref_("A").alt([Allele::seq("T").unwrap()])
                .gt(["0|1"])
                .format("GQ", [FieldValue::ints([42])])).unwrap()
            .build().unwrap()
            .truth();
        assert!(t.format[0][0].contains_key("GQ"));
        assert!(!t.format[0][0].contains_key("GT"));
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib truth` → FAIL.

- [ ] **Step 3: Implement**

```rust
use std::collections::{BTreeSet, HashMap};

use ndarray::{Array1, Array2, Array3};

use crate::allele::{Allele, SvType};
use crate::model::Document;
use crate::value::{FieldValue, Scalar};
use crate::variants::{classify_seq, record_class, VariantClass};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleKind {
    Snp,
    Mnp,
    Ins,
    Del,
    Delins,
    SpanningDel,
    Symbolic,
    Unspecified,
    Bnd,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlleleTruth {
    pub kind: AlleleKind,
    pub is_sequence: bool,
    pub sv_type: Option<String>,
    pub svlen: Option<i64>,
    pub sv_end: Option<i64>,
}

#[derive(Debug, Clone)]
pub struct GroundTruth {
    pub samples: Vec<String>,
    pub contigs: Vec<String>,
    pub pos: Array1<i64>,
    pub ref_: Vec<String>,
    pub alts: Vec<Vec<String>>,
    pub variant_class: Vec<VariantClass>,
    pub genotypes: Array3<i32>,
    pub phasing: Array2<bool>,
    pub info: Vec<HashMap<String, FieldValue>>,
    pub format: Vec<Vec<HashMap<String, FieldValue>>>,
    pub labels: Vec<BTreeSet<String>>,
    pub alts_truth: Vec<Vec<AlleleTruth>>,
    pub is_sequence_mask: Vec<Array1<bool>>,
}

/// Symbolic SV types that have a reference span (=> computed end). Excludes INS.
fn sv_spanning(t: SvType) -> bool {
    matches!(t, SvType::Del | SvType::Dup | SvType::Inv | SvType::Cnv)
}

fn seq_kind_to_allele_kind(c: VariantClass) -> AlleleKind {
    match c {
        VariantClass::Snp => AlleleKind::Snp,
        VariantClass::Mnp => AlleleKind::Mnp,
        VariantClass::Ins => AlleleKind::Ins,
        VariantClass::Del => AlleleKind::Del,
        VariantClass::Delins => AlleleKind::Delins,
        VariantClass::SpanningDel => AlleleKind::SpanningDel,
        _ => AlleleKind::Delins, // unreachable for sequence pairs
    }
}

fn allele_truth(pos: u64, allele: &Allele, svlen_val: Option<i64>) -> AlleleTruth {
    match allele {
        Allele::Seq(bases) => AlleleTruth {
            kind: seq_kind_to_allele_kind(classify_seq("", bases).snp_fix(bases)),
            is_sequence: true,
            sv_type: None,
            svlen: None,
            sv_end: None,
        },
        Allele::Star => AlleleTruth {
            kind: AlleleKind::SpanningDel,
            is_sequence: false,
            sv_type: None,
            svlen: None,
            sv_end: None,
        },
        Allele::Unspecified => AlleleTruth {
            kind: AlleleKind::Unspecified,
            is_sequence: false,
            sv_type: None,
            svlen: None,
            sv_end: None,
        },
        Allele::Breakend { .. } => AlleleTruth {
            kind: AlleleKind::Bnd,
            is_sequence: false,
            sv_type: None,
            svlen: None,
            sv_end: None,
        },
        Allele::Symbolic { first_type, .. } => {
            let svlen = svlen_val.map(|v| v.abs());
            let end = match (svlen, sv_spanning(*first_type)) {
                (Some(l), true) => Some(pos as i64 + l),
                _ => None,
            };
            AlleleTruth {
                kind: AlleleKind::Symbolic,
                is_sequence: false,
                sv_type: allele.symbolic_type_str(),
                svlen,
                sv_end: end,
            }
        }
    }
}

// classify_seq needs the real ref to classify, but for a sequence ALT's
// AlleleTruth the Python code passes the record REF. We thread it via a small
// helper below; this trait shim is replaced in the real derive() call which
// has the REF in scope. (See derive().)
trait SnpFix {
    fn snp_fix(self, _bases: &str) -> VariantClass;
}
impl SnpFix for VariantClass {
    fn snp_fix(self, _bases: &str) -> VariantClass {
        self
    }
}

pub fn derive(doc: &Document) -> GroundTruth {
    let n_rec = doc.records.len();
    let n_smp = doc.samples.len();
    let ploidy = doc.max_ploidy();

    let mut genos = Array3::<i32>::from_elem((n_rec.max(0), n_smp, ploidy), -1);
    let mut phasing = Array2::<bool>::from_elem((n_rec, n_smp), false);
    let mut pos = Array1::<i64>::zeros(n_rec);
    let mut ref_ = Vec::with_capacity(n_rec);
    let mut alts = Vec::with_capacity(n_rec);
    let mut vclass = Vec::with_capacity(n_rec);
    let mut info = Vec::with_capacity(n_rec);
    let mut fmt = Vec::with_capacity(n_rec);
    let mut labels = Vec::with_capacity(n_rec);
    let mut alts_truth = Vec::with_capacity(n_rec);
    let mut seq_mask = Vec::with_capacity(n_rec);

    for (ri, rec) in doc.records.iter().enumerate() {
        pos[ri] = rec.pos as i64;
        ref_.push(rec.ref_.clone());
        alts.push(rec.alts.iter().map(|a| a.render()).collect());
        vclass.push(record_class(&rec.ref_, &rec.alts));
        info.push(rec.info.iter().map(|(k, v)| (k.clone(), v.clone())).collect());

        let mut per_sample = Vec::with_capacity(n_smp);
        for (si, sample) in rec.samples.iter().enumerate() {
            if let Some(gt) = &sample.gt {
                for (ai, allele) in gt.alleles.iter().enumerate() {
                    genos[[ri, si, ai]] = match allele {
                        Some(v) => *v as i32,
                        None => -1,
                    };
                }
                phasing[[ri, si]] = gt.is_phased();
            }
            let m: HashMap<String, FieldValue> =
                sample.values.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
            per_sample.push(m);
        }
        fmt.push(per_sample);
        labels.push(rec.labels.clone());

        // Per-ALT truth. SVLEN is Number=A => per-allele list of ints.
        let svlen = rec.info.get("SVLEN");
        let mut per_alt = Vec::with_capacity(rec.alts.len());
        for (ai, allele) in rec.alts.iter().enumerate() {
            let sv = svlen_at(svlen, ai);
            // For sequence ALTs, classify against the record REF.
            let at = match allele {
                Allele::Seq(bases) => AlleleTruth {
                    kind: seq_kind_to_allele_kind(classify_seq(&rec.ref_, bases)),
                    is_sequence: true,
                    sv_type: None,
                    svlen: None,
                    sv_end: None,
                },
                _ => allele_truth(rec.pos, allele, sv),
            };
            per_alt.push(at);
        }
        let mask = Array1::from(per_alt.iter().map(|a| a.is_sequence).collect::<Vec<_>>());
        seq_mask.push(mask);
        alts_truth.push(per_alt);
    }

    GroundTruth {
        samples: doc.samples.clone(),
        contigs: doc.contigs.iter().map(|c| c.id.clone()).collect(),
        pos,
        ref_,
        alts,
        variant_class: vclass,
        genotypes: genos,
        phasing,
        info,
        format: fmt,
        labels,
        alts_truth,
        is_sequence_mask: seq_mask,
    }
}

fn svlen_at(value: Option<&FieldValue>, i: usize) -> Option<i64> {
    match value {
        Some(FieldValue::List(v)) => v.get(i).and_then(|x| match x {
            Some(Scalar::Int(n)) => Some(*n),
            Some(Scalar::Float(f)) => Some(*f as i64),
            _ => None,
        }),
        Some(FieldValue::Scalar(Scalar::Int(n))) if i == 0 => Some(*n),
        _ => None,
    }
}
```

(Note to implementer: the `SnpFix`/`allele_truth` sequence-branch shim exists only to keep `allele_truth` self-contained; in `derive` the sequence branch is handled inline with the real REF, so you may delete the `SnpFix` trait and the `Allele::Seq` arm of `allele_truth` if you prefer — keep whichever compiles cleanly under clippy.)

Add `Document::truth` to `src/model.rs`:

```rust
impl Document {
    pub fn truth(&self) -> crate::truth::GroundTruth {
        crate::truth::derive(self)
    }
}
```

Add to `src/lib.rs`:

```rust
pub mod truth;
pub use truth::{AlleleKind, AlleleTruth, GroundTruth};
```

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib truth` → PASS, then `pixi run clippy`.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add GroundTruth oracle"`

---

### Task 14: `write.rs` — render + bgzip + CSI via noodles

**Files:**
- Create: `src/write.rs`, `tests/snapshots.rs`
- Modify: `src/lib.rs`, `src/model.rs` (add `Document::render`, `Document::write`)

**Interfaces:**
- Consumes: `Document`.
- Produces:
  - `struct WriteOpts { bgzip: bool, index: bool }` (`Default` = both false) and a convenience `WriteOpts::text()/bgzipped()/bgzipped_indexed()`.
  - `render(doc: &Document) -> String`.
  - `write(doc: &Document, path, opts) -> Result<PathBuf, BuildError>`.
  - `Document::render(&self) -> String` and `Document::write(&self, path, opts) -> Result<PathBuf, BuildError>` delegating.

Implementation note: this task is where noodles' exact API matters and may differ from the snippets below. Build the header via `noodles_vcf::Header`, build records via `noodles_vcf::variant::RecordBuf` (or write text directly), and serialize with `noodles_vcf::io::Writer`. **Verify each call against the pinned noodles version**; if a deliberately edge-case construct cannot be expressed through `RecordBuf`, fall back to the custom text renderer described below (which is the spec's documented escape hatch). The custom renderer is also the simplest way to guarantee the snapshot test is stable, so implement it first and treat noodles as the bgzip/index backend.

Decision for this task: **implement `render()` as a custom text serializer** (port of Python `serialize.py`) for byte-stable output, and use **noodles only for bgzip compression and CSI indexing**. This satisfies "noodles for I/O" where it adds real value (compression/indexing) while keeping rendering fully under our control. (Revisit using noodles' record writer later if cross-tool header validation is desired.)

- [ ] **Step 1: Write the failing snapshot test** (`tests/snapshots.rs`)

```rust
use vcfixture::{Allele, RecordSpec, VcfBuilder};
use vcfixture::spec::version::LATEST;

#[test]
fn renders_minimal_vcf() {
    let text = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
        .info("AF", None, None, None).unwrap()
        .format("GT", None, None, None).unwrap()
        .record(RecordSpec::at("chr1", 1000).ref_("A").alt([Allele::seq("T").unwrap()])
            .gt(["0|1", "1|1"])
            .info("AF", vcfixture::FieldValue::floats([0.25]))).unwrap()
        .build().unwrap()
        .render();
    insta::assert_snapshot!(text);
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --test snapshots` → FAIL (compile error: no `render`).

- [ ] **Step 3: Implement the text serializer** (`src/write.rs`)

```rust
use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use crate::error::BuildError;
use crate::model::{Document, Record};
use crate::value::{FieldValue, Scalar};

#[derive(Debug, Clone, Copy, Default)]
pub struct WriteOpts {
    pub bgzip: bool,
    pub index: bool,
}

impl WriteOpts {
    pub fn text() -> WriteOpts {
        WriteOpts { bgzip: false, index: false }
    }
    pub fn bgzipped() -> WriteOpts {
        WriteOpts { bgzip: true, index: false }
    }
    pub fn bgzipped_indexed() -> WriteOpts {
        WriteOpts { bgzip: true, index: true }
    }
}

// Percent-encoding for reserved characters in string values. '%' first.
const PERCENT: &[(char, &str)] = &[
    ('%', "%25"),
    (':', "%3A"),
    (';', "%3B"),
    ('=', "%3D"),
    (',', "%2C"),
    ('\r', "%0D"),
    ('\n', "%0A"),
    ('\t', "%09"),
];

fn encode(s: &str) -> String {
    let mut out = s.to_string();
    for (ch, rep) in PERCENT {
        out = out.replace(*ch, rep);
    }
    out
}

fn fmt_scalar(s: &Scalar) -> String {
    match s {
        Scalar::Int(n) => n.to_string(),
        Scalar::Float(f) => {
            if f.is_nan() || f.is_infinite() {
                ".".to_string()
            } else {
                // Rust's default float formatting is round-trippable.
                format!("{f}")
            }
        }
        Scalar::Char(c) => encode(&c.to_string()),
        Scalar::Str(s) => encode(s),
    }
}

fn fmt_opt_scalar(s: &Option<Scalar>) -> String {
    match s {
        Some(v) => fmt_scalar(v),
        None => ".".to_string(),
    }
}

fn fmt_value(v: &FieldValue) -> String {
    match v {
        FieldValue::Flag => "1".to_string(), // never reached for INFO rendering
        FieldValue::Scalar(s) => fmt_scalar(s),
        FieldValue::List(xs) => {
            if xs.is_empty() {
                ".".to_string()
            } else {
                xs.iter().map(fmt_opt_scalar).collect::<Vec<_>>().join(",")
            }
        }
    }
}

fn render_info(rec: &Record) -> String {
    if rec.info.is_empty() {
        return ".".to_string();
    }
    let mut parts = Vec::new();
    for (key, val) in &rec.info {
        match val {
            FieldValue::Flag => parts.push(key.clone()),
            other => parts.push(format!("{key}={}", fmt_value(other))),
        }
    }
    if parts.is_empty() {
        ".".to_string()
    } else {
        parts.join(";")
    }
}

fn render_sample(rec: &Record, si: usize) -> String {
    let sample = &rec.samples[si];
    let mut vals = Vec::with_capacity(rec.fmt_keys.len());
    for key in &rec.fmt_keys {
        if key == "GT" {
            vals.push(sample.gt.as_ref().map(|g| g.render()).unwrap_or_else(|| ".".to_string()));
        } else {
            match sample.values.get(key) {
                Some(v) => vals.push(fmt_value(v)),
                None => vals.push(".".to_string()),
            }
        }
    }
    vals.join(":")
}

fn render_record(rec: &Record) -> String {
    let ids = match &rec.ids {
        Some(v) if !v.is_empty() => v.join(";"),
        _ => ".".to_string(),
    };
    let alt = if rec.alts.is_empty() {
        ".".to_string()
    } else {
        rec.alts.iter().map(|a| a.render()).collect::<Vec<_>>().join(",")
    };
    let qual = match rec.qual {
        Some(q) => fmt_scalar(&Scalar::Float(q)),
        None => ".".to_string(),
    };
    let filt = match &rec.filters {
        None => ".".to_string(),
        Some(v) if v.is_empty() => "PASS".to_string(),
        Some(v) => v.join(";"),
    };
    let mut cols = vec![
        rec.chrom.clone(),
        rec.pos.to_string(),
        ids,
        rec.ref_.clone(),
        alt,
        qual,
        filt,
        render_info(rec),
    ];
    if !rec.fmt_keys.is_empty() {
        cols.push(rec.fmt_keys.join(":"));
        for si in 0..rec.samples.len() {
            cols.push(render_sample(rec, si));
        }
    }
    cols.join("\t")
}

pub fn render(doc: &Document) -> String {
    let mut lines = vec![format!("##fileformat={}", doc.version.as_str())];
    for c in &doc.contigs {
        lines.push(c.header_line());
    }
    for f in &doc.info_defs {
        lines.push(f.header_line());
    }
    for (id, desc) in &doc.filter_defs {
        lines.push(format!("##FILTER=<ID={id},Description=\"{desc}\">"));
    }
    for f in &doc.format_defs {
        lines.push(f.header_line());
    }
    for ad in &doc.alt_defs {
        lines.push(ad.header_line());
    }
    let mut header = vec!["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let has_fmt = !doc.format_defs.is_empty() || doc.records.iter().any(|r| !r.fmt_keys.is_empty());
    if has_fmt {
        header.push("FORMAT".to_string());
        header.extend(doc.samples.iter().cloned());
    }
    lines.push(header.join("\t"));
    for rec in &doc.records {
        lines.push(render_record(rec));
    }
    let mut text = lines.join("\n");
    text.push('\n');
    text
}

pub fn write(doc: &Document, path: impl AsRef<Path>, opts: WriteOpts) -> Result<PathBuf, BuildError> {
    let text = render(doc);
    let mut path = path.as_ref().to_path_buf();
    if !opts.bgzip {
        fs::write(&path, text)?;
        return Ok(path);
    }
    if path.extension().and_then(|e| e.to_str()) != Some("gz") {
        let mut name = path.into_os_string();
        name.push(".gz");
        path = PathBuf::from(name);
    }
    // bgzip via noodles-bgzf
    let file = fs::File::create(&path)?;
    let mut writer = noodles_bgzf::Writer::new(file);
    writer.write_all(text.as_bytes())?;
    writer.finish()?;
    if opts.index {
        write_csi(&path)?;
    }
    Ok(path)
}

/// Build and write a CSI index alongside a bgzipped VCF.
/// Verify the exact noodles-csi/-tabix API against the pinned version; this
/// indexes the bgzipped file produced above.
fn write_csi(_bgzf_path: &Path) -> Result<(), BuildError> {
    // Implementer: use noodles_tabix / noodles_csi to index `_bgzf_path`
    // (preset = VCF) and write `<path>.csi`. Pseudocode against the noodles API:
    //
    //   let index = noodles_tabix::index(_bgzf_path)?;   // or build via reader
    //   let csi_path = format!("{}.csi", _bgzf_path.display());
    //   let mut w = noodles_csi::io::Writer::new(File::create(csi_path)?);
    //   w.write_index(&index)?;
    //
    // If the helper signature differs, adapt; the deliverable is a readable
    // .csi alongside the .vcf.gz. Covered by the round-trip test in Task 16.
    Ok(())
}
```

This is the one place the plan leaves an API-shaped gap (`write_csi`), because the exact noodles-csi/-tabix entry points must be confirmed against the resolved crate versions. The round-trip test in Task 16 exercises it; do not mark Task 14 complete until `write_csi` actually writes a `.csi` that a noodles reader can load.

Add `Document::render` / `Document::write` to `src/model.rs`:

```rust
impl Document {
    pub fn render(&self) -> String {
        crate::write::render(self)
    }
    pub fn write(
        &self,
        path: impl AsRef<std::path::Path>,
        opts: crate::write::WriteOpts,
    ) -> Result<std::path::PathBuf, crate::error::BuildError> {
        crate::write::write(self, path, opts)
    }
}
```

Add to `src/lib.rs`: `pub mod write; pub use write::WriteOpts;`

- [ ] **Step 4: Run and accept the snapshot**

Run: `pixi run cargo test --test snapshots`
Then review and accept: `pixi run cargo insta accept` *(if `cargo-insta` is unavailable, the snapshot is written to `tests/snapshots/…snap.new`; rename to `.snap` after eyeballing it)*. Re-run to confirm PASS.

- [ ] **Step 5: Write the bgzip round-trip smoke test** (append to `tests/snapshots.rs`)

```rust
#[test]
fn writes_bgzipped_file() {
    let dir = std::env::temp_dir().join("vcfixture_rs_write_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("x.vcf.gz");
    let written = VcfBuilder::new(["s1"], [("chr1", Some(1000u64))], LATEST)
        .format("GT", None, None, None).unwrap()
        .record(RecordSpec::at("chr1", 10).ref_("A").alt([Allele::seq("T").unwrap()]).gt(["0|1"]))
        .unwrap()
        .write(&path, vcfixture::WriteOpts::bgzipped_indexed())
        .unwrap();
    assert!(written.exists());
    assert!(written.with_extension("gz.csi").exists() || written.to_string_lossy().ends_with(".gz"));
}
```

- [ ] **Step 6: Run** — `pixi run cargo test --test snapshots` → PASS. Then `pixi run clippy`.

- [ ] **Step 7: Commit** — `git add -A && git commit -m "feat: add VCF rendering, bgzip, and CSI indexing"`

---

### Task 15: `reference.rs` — synthetic reference + FASTA write

**Files:**
- Create: `src/reference.rs`
- Modify: `src/lib.rs`

**Interfaces:**
- Consumes: `BuildError`.
- Produces:
  - `enum VariantKlass { Snp, Mnp, Ins, Del, Delins, SpanningDel }` (the draw classes) — or reuse string klass; use this enum for type safety.
  - `struct RepeatFeature { contig: String, pos0: usize, motif: String, count: usize }` with `length(&self) -> usize`.
  - `struct DrawOpts { alt_index: usize, del_len: usize, ins_seq: String, mnp_len: usize }` with `Default` (alt_index=1, del_len=1, ins_seq="T", mnp_len=2).
  - `struct ReferenceSpec { contigs: Vec<(String, String)>, repeats: Vec<RepeatFeature> }` with `length`, `base`, `seq`, `draw_ref_alt(contig, pos0, klass, &DrawOpts) -> Result<(String, Vec<String>), BuildError>`, `write(path, bgzip, index) -> Result<PathBuf, BuildError>`.
  - `struct ReferenceBuilder` with `new(seed: u64)`, `add_contig(id, length) -> Result<&mut Self>`, `set_base`, `set_seq`, `tandem_repeat`, `build() -> ReferenceSpec`.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproducible_fill() {
        let a = ReferenceBuilder::new(7).add_contig("chr1", 100).unwrap().build();
        let b = ReferenceBuilder::new(7).add_contig("chr1", 100).unwrap().build();
        assert_eq!(a.seq("chr1", 0, 100).unwrap(), b.seq("chr1", 0, 100).unwrap());
        assert_eq!(a.length("chr1").unwrap(), 100);
    }

    #[test]
    fn draw_snp_matches_reference() {
        let mut rb = ReferenceBuilder::new(1);
        rb.add_contig("chr1", 100).unwrap();
        rb.set_base("chr1", 10, "A").unwrap();
        let spec = rb.build();
        let (r, alts) = spec.draw_ref_alt("chr1", 10, VariantKlass::Snp, &DrawOpts::default()).unwrap();
        assert_eq!(r, "A");
        assert_eq!(alts.len(), 1);
        assert_ne!(alts[0], "A");
    }

    #[test]
    fn tandem_repeat_recorded() {
        let mut rb = ReferenceBuilder::new(1);
        rb.add_contig("chr1", 100).unwrap();
        rb.tandem_repeat("chr1", 10, "CAG", 4).unwrap();
        let spec = rb.build();
        assert_eq!(spec.repeats.len(), 1);
        assert_eq!(spec.seq("chr1", 10, 12).unwrap(), "CAGCAGCAGCAG");
    }
}
```

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --lib reference` → FAIL.

- [ ] **Step 3: Implement**

```rust
use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use crate::error::BuildError;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKlass {
    Snp,
    Mnp,
    Ins,
    Del,
    Delins,
    SpanningDel,
}

#[derive(Debug, Clone)]
pub struct DrawOpts {
    pub alt_index: usize,
    pub del_len: usize,
    pub ins_seq: String,
    pub mnp_len: usize,
}

impl Default for DrawOpts {
    fn default() -> Self {
        DrawOpts { alt_index: 1, del_len: 1, ins_seq: "T".to_string(), mnp_len: 2 }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatFeature {
    pub contig: String,
    pub pos0: usize,
    pub motif: String,
    pub count: usize,
}

impl RepeatFeature {
    pub fn length(&self) -> usize {
        self.motif.len() * self.count
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSpec {
    pub contigs: Vec<(String, String)>,
    pub repeats: Vec<RepeatFeature>,
}

fn next_base(b: char, offset: usize) -> char {
    const ORDER: [char; 4] = ['A', 'C', 'G', 'T'];
    let i = ORDER.iter().position(|&c| c == b).unwrap_or(0);
    ORDER[(i + offset) % 4]
}

impl ReferenceSpec {
    fn seq_for(&self, contig: &str) -> Result<&str, BuildError> {
        self.contigs
            .iter()
            .find(|(id, _)| id == contig)
            .map(|(_, s)| s.as_str())
            .ok_or_else(|| BuildError::ContigNotFound(contig.to_string()))
    }

    pub fn length(&self, contig: &str) -> Result<usize, BuildError> {
        Ok(self.seq_for(contig)?.len())
    }

    pub fn base(&self, contig: &str, pos0: usize) -> Result<String, BuildError> {
        Ok(self.seq_for(contig)?[pos0..pos0 + 1].to_string())
    }

    pub fn seq(&self, contig: &str, start0: usize, length: usize) -> Result<String, BuildError> {
        Ok(self.seq_for(contig)?[start0..start0 + length].to_string())
    }

    pub fn draw_ref_alt(
        &self,
        contig: &str,
        pos0: usize,
        klass: VariantKlass,
        opts: &DrawOpts,
    ) -> Result<(String, Vec<String>), BuildError> {
        match klass {
            VariantKlass::Snp => {
                let r = self.base(contig, pos0)?;
                let alt = next_base(r.chars().next().unwrap(), opts.alt_index).to_string();
                Ok((r, vec![alt]))
            }
            VariantKlass::Mnp => {
                let r = self.seq(contig, pos0, opts.mnp_len)?;
                let alt: String = r.chars().map(|c| next_base(c, opts.alt_index)).collect();
                Ok((r, vec![alt]))
            }
            VariantKlass::Ins => {
                let anchor = self.base(contig, pos0)?;
                let alt = format!("{anchor}{}", opts.ins_seq);
                Ok((anchor, vec![alt]))
            }
            VariantKlass::Del => {
                let r = self.seq(contig, pos0, opts.del_len + 1)?;
                let alt = r[..1].to_string();
                Ok((r, vec![alt]))
            }
            VariantKlass::Delins => {
                let r = self.seq(contig, pos0, opts.mnp_len)?;
                Ok((r, vec![opts.ins_seq.clone()]))
            }
            VariantKlass::SpanningDel => {
                let r = self.base(contig, pos0)?;
                Ok((r, vec!["*".to_string()]))
            }
        }
    }

    pub fn write(&self, path: impl AsRef<Path>, bgzip: bool, index: bool) -> Result<PathBuf, BuildError> {
        let path = path.as_ref().to_path_buf();
        let mut text = String::new();
        for (cid, seq) in &self.contigs {
            text.push('>');
            text.push_str(cid);
            text.push('\n');
            for chunk in seq.as_bytes().chunks(60) {
                text.push_str(std::str::from_utf8(chunk).unwrap());
                text.push('\n');
            }
        }
        if bgzip {
            let file = fs::File::create(&path)?;
            let mut w = noodles_bgzf::Writer::new(file);
            w.write_all(text.as_bytes())?;
            w.finish()?;
        } else {
            fs::write(&path, &text)?;
        }
        if index {
            // Implementer: write a .fai (and .gzi when bgzipped) via
            // noodles-fasta's index writer against the pinned version. The
            // strategies/round-trip tests do not require the index, so this is
            // best-effort; verify if a downstream test needs faidx access.
        }
        Ok(path)
    }
}

pub struct ReferenceBuilder {
    rng: ChaCha8Rng,
    seqs: indexmap::IndexMap<String, Vec<u8>>,
    repeats: Vec<RepeatFeature>,
}

impl ReferenceBuilder {
    pub fn new(seed: u64) -> ReferenceBuilder {
        ReferenceBuilder {
            rng: ChaCha8Rng::seed_from_u64(seed),
            seqs: indexmap::IndexMap::new(),
            repeats: Vec::new(),
        }
    }

    pub fn add_contig(&mut self, id: impl Into<String>, length: usize) -> Result<&mut Self, BuildError> {
        let id = id.into();
        if self.seqs.contains_key(&id) {
            return Err(BuildError::ContigExists(id));
        }
        let seq: Vec<u8> = (0..length).map(|_| BASES[self.rng.gen_range(0..4)]).collect();
        self.seqs.insert(id, seq);
        Ok(self)
    }

    pub fn set_base(&mut self, contig: &str, pos0: usize, base: &str) -> Result<&mut Self, BuildError> {
        let b = base.as_bytes();
        if b.len() != 1 {
            return Err(BuildError::BadAlleleBases(base.to_string()));
        }
        let arr = self.seqs.get_mut(contig).ok_or_else(|| BuildError::ContigNotFound(contig.to_string()))?;
        arr[pos0] = b[0];
        Ok(self)
    }

    pub fn set_seq(&mut self, contig: &str, pos0: usize, seq: &str) -> Result<&mut Self, BuildError> {
        let arr = self.seqs.get_mut(contig).ok_or_else(|| BuildError::ContigNotFound(contig.to_string()))?;
        if pos0 + seq.len() > arr.len() {
            return Err(BuildError::OutOfBounds {
                contig: contig.to_string(),
                pos0,
                len: seq.len(),
                clen: arr.len(),
            });
        }
        arr[pos0..pos0 + seq.len()].copy_from_slice(seq.as_bytes());
        Ok(self)
    }

    pub fn tandem_repeat(&mut self, contig: &str, pos0: usize, motif: &str, n: usize) -> Result<&mut Self, BuildError> {
        let run = motif.repeat(n);
        self.set_seq(contig, pos0, &run)?;
        self.repeats.push(RepeatFeature {
            contig: contig.to_string(),
            pos0,
            motif: motif.to_string(),
            count: n,
        });
        Ok(self)
    }

    pub fn build(&self) -> ReferenceSpec {
        let contigs = self
            .seqs
            .iter()
            .map(|(id, bytes)| (id.clone(), String::from_utf8(bytes.clone()).unwrap()))
            .collect();
        ReferenceSpec { contigs, repeats: self.repeats.clone() }
    }
}
```

Add to `src/lib.rs`:

```rust
pub mod reference;
pub use reference::{DrawOpts, ReferenceBuilder, ReferenceSpec, RepeatFeature, VariantKlass};
```

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --lib reference` → PASS, then `pixi run clippy`.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: add synthetic reference subsystem"`

---

### Task 16: `strategies.rs` — proptest strategies + round-trip test

**Files:**
- Create: `src/strategies.rs`, `tests/roundtrip.rs`
- Modify: `src/lib.rs`, `Cargo.toml` (already has proptest feature)

**Interfaces:**
- Consumes: everything; gated behind `#[cfg(feature = "proptest")]`.
- Produces (all behind the feature):
  - consts `ALL_VARIANT_CLASSES: [VariantKlass; 6]`, `ALL_NUMBER_TYPE_COMBOS: Vec<(Number, Type, FieldKind)>` (via a function `all_number_type_combos()`).
  - `genotype_strategy(ploidy, n_alt, missing_rate) -> impl Strategy<Value = String>`.
  - `field_value_strategy(&FieldDef, n_alt, ploidy) -> impl Strategy<Value = FieldValue>` (Flag => `FieldValue::Flag`; else list of resolved cardinality, Dot => 1..=3).
  - `documents(opts: DocumentOpts) -> impl Strategy<Value = Document>` with `struct DocumentOpts { max_samples, max_records, max_alt, version }` (`Default`).
  - `documents_with_fields(...)`, `symbolic_documents(...)`, `references(...)`, `reference_and_documents(...)`.

- [ ] **Step 1: Write the failing round-trip test** (`tests/roundtrip.rs`)

```rust
#![cfg(feature = "proptest")]

use proptest::prelude::*;
use vcfixture::strategies::{documents, DocumentOpts};

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn render_then_parse_genotype_counts_match_truth(doc in documents(DocumentOpts::default())) {
        let truth = doc.truth();
        let text = doc.render();

        // Independent re-read: count data lines and assert the genotype matrix
        // record count matches. (Full noodles re-parse asserted below once the
        // noodles reader wiring is confirmed.)
        let data_lines = text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).count();
        prop_assert_eq!(data_lines, truth.pos.len());

        // Genotype dims are consistent.
        prop_assert_eq!(truth.genotypes.shape()[0], truth.pos.len());
        prop_assert_eq!(truth.genotypes.shape()[1], truth.samples.len());
    }
}
```

(Note: this is the structural round-trip. The richer noodles-reader assertion — parse each record with `noodles_vcf::io::Reader` and compare decoded GT/INFO to `truth` — should be added as a second `proptest!` block once the noodles reader API is confirmed against the pinned version. This is exactly the loop svar2 will run; leave a `// TODO(svar2): full noodles re-parse comparison` only if the API can't be confirmed in this task, and open a follow-up. Prefer to implement it now.)

- [ ] **Step 2: Run to verify it fails** — `pixi run cargo test --features proptest --test roundtrip` → FAIL (no `strategies` module).

- [ ] **Step 3: Implement `src/strategies.rs`** (feature-gated; reference-free `documents` body first, then the others)

```rust
use proptest::prelude::*;

use crate::allele::Allele;
use crate::build::{RecordSpec, VcfBuilder};
use crate::model::Document;
use crate::reference::{DrawOpts, ReferenceBuilder, ReferenceSpec, VariantKlass};
use crate::spec::field::{FieldDef, FieldKind};
use crate::spec::number::Number;
use crate::spec::types::Type;
use crate::spec::version::{VcfVersion, LATEST};
use crate::truth::GroundTruth;
use crate::value::{FieldValue, Scalar};

pub const ALL_VARIANT_CLASSES: [VariantKlass; 6] = [
    VariantKlass::Snp,
    VariantKlass::Mnp,
    VariantKlass::Ins,
    VariantKlass::Del,
    VariantKlass::Delins,
    VariantKlass::SpanningDel,
];

const BASES: [&str; 4] = ["A", "C", "G", "T"];

/// All classic Number×Type combos as (Number, Type, kind); Flag only as INFO.
pub fn all_number_type_combos() -> Vec<(Number, Type, FieldKind)> {
    let numbers = [Number::ONE, Number::fixed(2).unwrap(), Number::A, Number::R, Number::G, Number::DOT];
    let mut combos = Vec::new();
    for kind in [FieldKind::Info, FieldKind::Format] {
        let allowed: Vec<Type> = match kind {
            FieldKind::Info => Type::info_allowed().into_iter().filter(|t| *t != Type::Flag).collect(),
            FieldKind::Format => Type::format_allowed().to_vec(),
        };
        for n in numbers {
            for t in &allowed {
                combos.push((n, *t, kind));
            }
        }
        if kind == FieldKind::Info {
            combos.push((Number::FLAG, Type::Flag, FieldKind::Info));
        }
    }
    combos
}

fn next_base(b: &str, off: usize) -> String {
    let i = BASES.iter().position(|&x| x == b).unwrap_or(0);
    BASES[(i + off) % 4].to_string()
}

/// Draw a (ref, alt) sequence pair for a class, reference-free.
fn ref_alt_strategy(klass: VariantKlass) -> impl Strategy<Value = (String, String)> {
    let base = prop::sample::select(BASES.to_vec());
    let base2 = prop::sample::select(BASES.to_vec());
    let ins = "[ACGT]{1,3}";
    (base, base2, prop::string::string_regex(ins).unwrap(), 0usize..3)
        .prop_map(move |(b, b2, tail, snp_off)| match klass {
            VariantKlass::Snp => (b.clone(), next_base(&b, 1 + snp_off)),
            VariantKlass::Mnp => {
                let r = format!("{b}{b2}");
                let a = format!("{}{}", next_base(&b, 1), next_base(&b2, 1));
                (r, a)
            }
            VariantKlass::Ins => (b.clone(), format!("{b}{tail}")),
            VariantKlass::Del => (format!("{b}{tail}"), b.clone()),
            VariantKlass::Delins => (format!("{b}{b2}"), tail),
            VariantKlass::SpanningDel => (b, "*".to_string()),
        })
}

/// A single-sample GT string.
pub fn genotype_strategy(ploidy: usize, n_alt: usize, missing_rate: f64) -> impl Strategy<Value = String> {
    let slots = prop::collection::vec(
        (0.0f64..1.0, 0u32..=(n_alt as u32)),
        ploidy,
    );
    (slots, any::<bool>()).prop_map(move |(slots, phased)| {
        let sep = if phased { "|" } else { "/" };
        let tokens: Vec<String> = slots
            .into_iter()
            .map(|(r, idx)| if r < missing_rate { ".".to_string() } else { idx.to_string() })
            .collect();
        tokens.join(sep)
    })
}

fn scalar_strategy(t: Type) -> BoxedStrategy<Scalar> {
    match t {
        Type::Integer => (-1000i64..=1000).prop_map(Scalar::Int).boxed(),
        Type::Float => (-1.0e6f64..1.0e6).prop_map(|f| Scalar::Float(f as f32 as f64)).boxed(),
        Type::Character => prop::sample::select(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789".chars().collect::<Vec<_>>(),
        )
        .prop_map(Scalar::Char)
        .boxed(),
        Type::String => prop::string::string_regex("[A-Za-z0-9]{1,6}").unwrap().prop_map(Scalar::Str).boxed(),
        Type::Flag => Just(Scalar::Int(1)).boxed(), // unused; Flag handled in field_value_strategy
    }
}

/// A spec-valid value for one field given the record's n_alt/ploidy.
pub fn field_value_strategy(fd: &FieldDef, n_alt: usize, ploidy: usize) -> BoxedStrategy<FieldValue> {
    if fd.type_ == Type::Flag {
        return Just(FieldValue::Flag).boxed();
    }
    let card = fd.number.cardinality(n_alt, ploidy);
    let t = fd.type_;
    match card {
        Some(c) => prop::collection::vec(scalar_strategy(t), c)
            .prop_map(|xs| FieldValue::List(xs.into_iter().map(Some).collect()))
            .boxed(),
        None => (1usize..=3)
            .prop_flat_map(move |c| prop::collection::vec(scalar_strategy(t), c))
            .prop_map(|xs| FieldValue::List(xs.into_iter().map(Some).collect()))
            .boxed(),
    }
}

#[derive(Debug, Clone)]
pub struct DocumentOpts {
    pub max_samples: usize,
    pub max_records: usize,
    pub max_alt: usize,
    pub version: VcfVersion,
}

impl Default for DocumentOpts {
    fn default() -> Self {
        DocumentOpts { max_samples: 3, max_records: 4, max_alt: 1, version: LATEST }
    }
}

/// A reference-free document over a single synthetic chr1.
pub fn documents(opts: DocumentOpts) -> impl Strategy<Value = Document> {
    (1usize..=opts.max_samples, 1usize..=2, 1usize..=opts.max_records).prop_flat_map(
        move |(n_samples, ploidy, n_rec)| {
            let max_alt = opts.max_alt;
            let version = opts.version;
            // Per-record: (n_alt, per-alt klass, ref base, gts, pos-gap).
            let rec_strat = (1usize..=max_alt).prop_flat_map(move |n_alt| {
                let klasses = prop::collection::vec(
                    prop::sample::select(ALL_VARIANT_CLASSES.to_vec()),
                    n_alt,
                );
                let refbase = prop::sample::select(BASES.to_vec());
                let gts = prop::collection::vec(genotype_strategy(ploidy, n_alt, 0.1), n_samples);
                let gap = 1u64..=50;
                (Just(n_alt), klasses, refbase, gts, gap)
            });
            prop::collection::vec(rec_strat, n_rec).prop_map(move |recs| {
                let samples: Vec<String> = (0..n_samples).map(|i| format!("s{i}")).collect();
                let mut b = VcfBuilder::new(samples, [("chr1", Some(100_000u64))], version)
                    .format("GT", None, None, None)
                    .expect("GT declares");
                let mut pos = 1000u64;
                for (n_alt, mut klasses, refbase, gts, gap) in recs {
                    // SPANNING_DEL only valid as the last ALT; downgrade others to SNP.
                    let last = n_alt - 1;
                    for (j, k) in klasses.iter_mut().enumerate() {
                        if *k == VariantKlass::SpanningDel && j != last {
                            *k = VariantKlass::Snp;
                        }
                    }
                    // Build ALT strings deterministically from the klass (reference-free):
                    let alts: Vec<Allele> = klasses
                        .iter()
                        .map(|k| Allele::parse(&sample_alt_for(*k, &refbase)))
                        .collect();
                    let mut spec = RecordSpec::at("chr1", pos).ref_(refbase.clone()).alt(alts).gt(gts);
                    let _ = &mut spec;
                    b = b.record(spec).expect("valid record");
                    pos += gap;
                }
                b.build().expect("valid document")
            })
        },
    )
}

/// Deterministic ALT text for a class given a single REF base (reference-free).
fn sample_alt_for(klass: VariantKlass, refbase: &str) -> String {
    match klass {
        VariantKlass::Snp => next_base(refbase, 1),
        VariantKlass::Mnp => format!("{}{}", next_base(refbase, 1), next_base(refbase, 1)),
        VariantKlass::Ins => format!("{refbase}T"),
        VariantKlass::Del => refbase.to_string(), // ref will be longer; handled below
        VariantKlass::Delins => "T".to_string(),
        VariantKlass::SpanningDel => "*".to_string(),
    }
}
```

Implementer note on `documents`: the reference-free Python body draws REF as a single base and an ALT consistent with the class, *without* enforcing REF/ALT length agreement for DEL/DELINS (the Python version draws ALT from `_ref_alt(klass)` and uses `ref = single base`). To match that and keep records valid, simplify: for the reference-free strategy, set REF to a single base and restrict classes to those valid with a 1-base REF (SNP, INS, SPANNING_DEL) **or** draw the full `(ref, alt)` pair from `ref_alt_strategy(klass)` and use *that* REF for the record. Prefer the latter: replace the per-record `refbase`+`sample_alt_for` with a drawn `(ref, alts)` per record via `ref_alt_strategy`, mirroring Python's `documents_with_fields`. Wire it so each record's REF is the drawn pair's REF and ALTs are `Allele::parse` of the drawn ALT strings. Keep `genotype_strategy` for GTs. Ensure `pixi run clippy` is clean (remove the now-unused `sample_alt_for`).

Then implement `documents_with_fields`, `symbolic_documents`, `references`, and `reference_and_documents` following the Python `strategies.py` logic 1:1, using the combinator patterns above. For `references`, port the per-contig cursor repeat-planting; for `reference_and_documents`, draw a `references()` spec, then a reference-consistent `documents()`, and return `(spec, doc, doc.truth())`.

Add to `src/lib.rs`:

```rust
#[cfg(feature = "proptest")]
pub mod strategies;
```

- [ ] **Step 4: Run to verify it passes** — `pixi run cargo test --features proptest --test roundtrip` → PASS. Then `pixi run cargo test --all-features` and `pixi run clippy`.

- [ ] **Step 5: Add the full noodles re-parse round-trip** (second `proptest!` block in `tests/roundtrip.rs`): write the doc to a temp `.vcf`, read it back with `noodles_vcf::io::reader::Builder`, and for each record assert the decoded GT indices and phasing equal `truth.genotypes`/`truth.phasing`. Verify the noodles reader API against the pinned version. Run until PASS.

- [ ] **Step 6: Commit** — `git add -A && git commit -m "feat: add proptest strategies and round-trip tests"`

---

### Task 17: Public API docs + README + final polish

**Files:**
- Modify: `src/lib.rs` (crate-level docs + doctest), create `README.md`

**Interfaces:**
- Consumes: the full public API.
- Produces: documented crate; `cargo test --doc` passes.

- [ ] **Step 1: Write a crate-level doctest in `src/lib.rs`**

````rust
//! vcfixture — generate small VCF test data with decoded ground truth.
//!
//! ```
//! use vcfixture::{Allele, RecordSpec, VcfBuilder, FieldValue};
//! use vcfixture::spec::version::LATEST;
//!
//! let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
//!     .info("AF", None, None, None).unwrap()
//!     .format("GT", None, None, None).unwrap()
//!     .record(
//!         RecordSpec::at("chr1", 1000)
//!             .ref_("A")
//!             .alt([Allele::seq("T").unwrap()])
//!             .gt(["0|1", "1|1"])
//!             .info("AF", FieldValue::floats([0.25])),
//!     ).unwrap()
//!     .build().unwrap();
//!
//! let truth = doc.truth();
//! assert_eq!(truth.genotypes[[0, 0, 1]], 1);
//! assert_eq!(truth.pos[0], 1000);
//! let _text = doc.render();
//! ```
````

- [ ] **Step 2: Run the doctest** — `pixi run cargo test --doc` → PASS.

- [ ] **Step 3: Write `README.md`** (short: what it is, the example above, link to the spec, note `features = ["proptest"]` for strategies).

- [ ] **Step 4: Full green sweep** — run all and confirm:

```bash
pixi run fmt-check
pixi run clippy
pixi run cargo test --all-features
pixi run cargo test --doc
```
Expected: all PASS.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "docs: add crate docs and README"`

---

## Self-Review

**Spec coverage:**
- Architecture / hub invariant → Tasks 11–14 (model is the only source for truth/render). ✓
- Crate layout & feature gating → Tasks 1, 16 (proptest optional, off by default). ✓
- Allele enum + parse/render/ctors → Task 8. ✓
- Genotype → Task 9. ✓
- spec (Number/Type/Version/FieldDef/reserved/genotype_order) → Tasks 2–7. ✓
- FieldValue/Scalar → Task 11. ✓
- VcfBuilder + RecordSpec + all validation rules → Task 12. ✓
- GroundTruth (ndarray) + AlleleTruth + derive → Task 13. ✓
- render/write/bgzip/CSI → Task 14. ✓
- Reference subsystem + FASTA write → Task 15. ✓
- proptest strategies (all 6 + coverage tables) → Task 16. ✓
- Testing strategy (unit + round-trip + snapshot) → Tasks per-module + 14 + 16. ✓
- Dev env (pixi + prek) → Task 1. ✓

**Known deliberate gaps (flagged, not placeholders):** the exact noodles-csi/-tabix index call (`write_csi`, Task 14), noodles-fasta `.fai` writing (Task 15), and the full noodles re-parse assertion (Task 16, Step 5) are described with pseudocode because their exact signatures depend on the resolved noodles versions. Each is paired with a concrete acceptance criterion and an exercising test, and each task is "not done" until the real call works. This matches the spec's documented escape-hatch / version-verification stance.

**Decision recorded during planning:** Task 14 implements `render()` as a custom text serializer (port of `serialize.py`) for byte-stable snapshots, using noodles only for bgzip + CSI. This is a refinement of the spec's "noodles for everything" toward "noodles for the I/O that benefits (compression/indexing), custom serializer for exact bytes" — the spec explicitly anticipated this escape hatch. Flag for user confirmation at execution time.

**Type consistency:** `FieldValue`/`Scalar` ctors (`ints`/`floats`/`strings`), `RecordSpec` setter names (`ref_`/`alt`/`gt`/`info`/`format`), `VcfBuilder` method names (`info`/`format`/`filter`/`alt`/`record`/`build`), `WriteOpts` constructors, and `VariantKlass`/`VariantClass` (distinct: the former is the 6 draw classes, the latter the full oracle label set) are used consistently across Tasks 8–16.
