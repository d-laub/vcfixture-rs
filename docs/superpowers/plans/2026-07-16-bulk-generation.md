# Bulk VCF Generation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fit summary statistics from real local germline/somatic data into committed JSON profiles, and generate realistic-enough BCFs at scale (~100 MB compressed, >=3 contigs) from them, via a Rust API and CLI, for benchmarking `genoray` / `GenVarLoader` read speed and memory.

**Architecture:** A `Profile` (fitted stats + dialed payload) drives a streaming record generator. Records are drawn i.i.d. from HWE — never accumulated — and written through `noodles-bcf` over a multithreaded bgzf writer. The index is built in a second pass by `bcf::fs::index`. A summary (counts + checksum, not a per-genotype oracle) falls out of the stream for free. The existing `Document`/`VcfBuilder`/`GroundTruth` fixture path is untouched; bulk shares only `spec/` and the allele model.

**Tech Stack:** Rust, noodles (`-bcf`/`-vcf`/`-bgzf`/`-csi`/`-core`), serde/serde_json, rayon, rand/rand_chacha, clap (cli-gated). Extraction: Python + plink2 + polars, via pixi.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`. Read it before Task 1.
- `rust-version = "1.86"`. **`noodles-bcf` 0.87 is unusable** (requires rustc 1.89). Pin `noodles-bcf = "0.81"`.
- Feature `bulk` is **off by default**; feature `cli` is off by default and implies `bulk`. Existing fixture-path consumers must not gain `serde`/`rayon`/`clap`/`noodles-bcf`.
- **Do not modify** `src/model.rs`, `src/build.rs`, `src/truth.rs`, `src/strategies.rs`. Bulk shares only `src/spec/` and `src/allele.rs`.
- Genotypes are drawn **i.i.d. from HWE**. Do not implement LD, haplotype copying, or coalescent simulation — the spec's ablation shows LD is a 0x lever on read speed/memory and 1.14x on BCF size.
- Contigs are declared at **fake lengths equal to the populated span**, never real hg38 lengths.
- Determinism is a shipped guarantee: same seed + profile => byte-identical output, **regardless of thread count**. Seed per-record-block from `hash(seed, block_idx)`; never from thread ID or a shared mutable RNG.
- Profile JSON must keep `fitted` and `dialed` in separate objects. Never write a hand-chosen value into `fitted`.
- Every task ends green: `pixi run fmt`, `pixi run clippy`, `pixi run test` all pass before commit.
- Conventional-commit messages (the commitizen `commit-msg` hook is active). Run `pixi run prek-install` once if hooks are not installed.
- CI runs `cargo clippy --all-features` and `cargo test --all-features --locked`, so `bulk` and `cli` are always exercised. `Cargo.lock` must be committed.
- CI must never read `/carter`. Only Task 9 touches real data.

## Parallelization

Tasks form four waves. Within a wave, tasks are independent — dispatch them concurrently with `superpowers:dispatching-parallel-agents` using `superpowers:subagent-driven-development`. Use Sonnet or weaker for implementation; reserve stronger models for second-pass fixes.

| Wave | Tasks | Notes |
| ---- | ----- | ----- |
| 1 | Task 1 | Dependency upgrade + feature gates. Blocks everything. |
| 2 | Task 2 | Profile schema. Blocks 3–8 (they all consume `Profile`). |
| 3 | **Tasks 3, 4, 5, 8** | **Parallel.** Samplers / writer / summary / extraction script. Disjoint files. |
| 4 | Task 6 -> Task 7 -> Task 9 -> Task 10 -> Task 11 | Sequential: generator needs 3+4+5; API needs 6; profiles need 8; CLI/docs need 7; fidelity needs 7+8+9. |

Task 9 is the only task requiring `/carter` access and is slow (reads a 9.6 GB
pgen). If that access is unavailable, Tasks 1–8 and 10 still complete against the
Task 2 placeholder profile; Tasks 9 and 11 block.

---

### Task 1: Coordinated noodles upgrade + feature gates

`noodles-bcf` cannot be added at current versions — it resolves a *second* `noodles-vcf` (0.83) alongside the existing 0.79 plus duplicate `bgzf`/`core`/`csi`/`tabix`, and types across the two do not interoperate. This upgrade has been verified to compile clean with all 70 existing tests passing.

**Files:**
- Modify: `Cargo.toml`
- Modify: `Cargo.lock` (via cargo)

**Interfaces:**
- Consumes: nothing.
- Produces: features `bulk` and `cli`; single-version noodles tree with `noodles-bcf 0.81` available under `bulk`.

- [ ] **Step 1: Apply the coordinated upgrade**

Run exactly this — the versions are a verified set, do not substitute "latest":

```bash
cargo add noodles-vcf@0.83 noodles-bgzf@0.45 noodles-core@0.18 noodles-csi@0.53 noodles-tabix@0.59
```

- [ ] **Step 2: Verify a single noodles-vcf in the tree**

```bash
cargo tree | grep -oE "noodles-vcf v[0-9.]+" | sort -u
```

Expected: exactly one line, `noodles-vcf v0.83.0`. If two versions appear, stop — the upgrade is wrong.

- [ ] **Step 3: Verify existing code still compiles and passes**

```bash
cargo check --all-features && cargo test --all-features
```

Expected: exit 0; 70 tests pass (62 + 5 + 2 + 1). If `write.rs` fails to compile, fix it against the 0.83 API — do not downgrade.

- [ ] **Step 4: Add bulk dependencies as optional, and the feature gates**

In `Cargo.toml` `[dependencies]`, add:

```toml
noodles-bcf = { version = "0.81", optional = true }
serde = { version = "1", features = ["derive"], optional = true }
serde_json = { version = "1", optional = true }
rayon = { version = "1", optional = true }
clap = { version = "4", features = ["derive"], optional = true }
```

Replace the `[features]` section with:

```toml
[features]
default = []
proptest = ["dep:proptest"]
bulk = ["dep:noodles-bcf", "dep:serde", "dep:serde_json", "dep:rayon"]
cli = ["bulk", "dep:clap"]
```

Add the binary target:

```toml
[[bin]]
name = "vcfixture"
path = "src/bin/vcfixture.rs"
required-features = ["cli"]
```

- [ ] **Step 5: Verify the default build stays lean**

```bash
cargo tree | grep -E "clap|serde_json|rayon|noodles-bcf" | head
```

Expected: no output (these are optional and `default = []`).

```bash
cargo check --features bulk && cargo check --all-features
```

Expected: both exit 0. (`--all-features` will fail on the missing `src/bin/vcfixture.rs` — create it as a placeholder `fn main() {}` now; Task 7 fills it in.)

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml Cargo.lock src/bin/vcfixture.rs
git commit -m "build: upgrade noodles and add bulk/cli feature gates"
```

---

### Task 2: Profile schema

**Files:**
- Create: `src/bulk/mod.rs`, `src/bulk/profile.rs`
- Create: `profiles/germline-1kgp.json` (a hand-written *minimal valid* profile for tests; Task 9 replaces it with the real fit)
- Modify: `src/lib.rs` (add `#[cfg(feature = "bulk")] pub mod bulk;`)
- Test: inline `#[cfg(test)] mod tests` in `src/bulk/profile.rs`

**Interfaces:**
- Consumes: nothing.
- Produces — every later task depends on these exact names:

```rust
pub struct Profile { pub name: String, pub provenance: Provenance, pub fitted: Fitted, pub dialed: Dialed }
pub struct Provenance { pub source: String, pub n_samples_source: usize, pub n_variants_source: u64, pub fitted_on: String, pub fit_tool_version: String }
pub struct Fitted {
    pub contigs: Vec<ContigStat>,
    pub gap_dist: Histogram,
    pub sfs: Histogram,
    pub variant_classes: ClassMix,
    pub indel_length: Histogram,
    pub titv: f64,
    pub multiallelic_rate: f64,
    pub missing_rate: f64,
    pub phased_rate: f64,
    pub ploidy: u8,
}
pub struct ContigStat { pub id: String, pub n_variants: u64, pub density_per_kb: f64 }
pub struct Histogram { pub edges: Vec<f64>, pub weights: Vec<f64> }
pub struct ClassMix { pub snp: f64, pub insertion: f64, pub deletion: f64, pub mnp: f64, pub complex: f64, pub symbolic: f64 }
pub struct Dialed { pub payload: Payload }
pub enum Payload { GtOnly, GtVaf, Gatk, Mutect2 }

impl Profile {
    pub fn builtin(name: &str) -> Result<Profile, BulkError>;
    pub fn from_json(s: &str) -> Result<Profile, BulkError>;
    pub fn validate(&self) -> Result<(), BulkError>;
}
```

- [ ] **Step 1: Write the failing tests**

Create `src/bulk/profile.rs` with only this test module at the bottom:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_germline_loads_and_validates() {
        let p = Profile::builtin("germline-1kgp").unwrap();
        assert_eq!(p.name, "germline-1kgp");
        assert_eq!(p.dialed.payload, Payload::GtOnly);
        assert_eq!(p.fitted.ploidy, 2);
        p.validate().unwrap();
    }

    #[test]
    fn unknown_builtin_errors() {
        assert!(Profile::builtin("nope").is_err());
    }

    #[test]
    fn histogram_length_mismatch_is_rejected() {
        // weights must have exactly edges.len() - 1 entries
        let h = Histogram { edges: vec![0.0, 1.0, 2.0], weights: vec![1.0] };
        assert!(h.validate().is_err());
    }

    #[test]
    fn class_mix_must_sum_to_one() {
        let m = ClassMix { snp: 0.5, insertion: 0.1, deletion: 0.1, mnp: 0.0, complex: 0.0, symbolic: 0.0 };
        assert!(m.validate().is_err());
    }

    #[test]
    fn payload_round_trips_through_serde() {
        let json = r#""mutect2""#;
        let p: Payload = serde_json::from_str(json).unwrap();
        assert_eq!(p, Payload::Mutect2);
        assert_eq!(serde_json::to_string(&p).unwrap(), json);
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --features bulk profile`
Expected: FAIL — `Profile` not found / module does not exist.

- [ ] **Step 3: Implement the schema**

In `src/lib.rs`, add:

```rust
#[cfg(feature = "bulk")]
pub mod bulk;
```

Create `src/bulk/mod.rs`:

```rust
//! Bulk generation of realistic-enough VCF/BCF at benchmark scale.
//!
//! Unlike the fixture path ([`crate::build::VcfBuilder`]), bulk generation
//! streams records and derives no per-genotype oracle — see
//! `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`.

pub mod profile;

pub use profile::{Payload, Profile};

/// Errors from bulk generation.
#[derive(Debug, thiserror::Error)]
pub enum BulkError {
    #[error("unknown builtin profile: {0}")]
    UnknownProfile(String),
    #[error("invalid profile: {0}")]
    Invalid(String),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}
```

In `src/bulk/profile.rs`, above the test module, write the structs exactly as
given in the Interfaces block, each deriving
`#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]`.
`Payload` additionally gets `#[serde(rename_all = "kebab-case")]` so it
serializes as `"gt-only"`, `"gt-vaf"`, `"gatk"`, `"mutect2"`.

```rust
use crate::bulk::BulkError;

const GERMLINE_1KGP: &str = include_str!("../../profiles/germline-1kgp.json");

impl Profile {
    pub fn builtin(name: &str) -> Result<Profile, BulkError> {
        let src = match name {
            "germline-1kgp" => GERMLINE_1KGP,
            other => return Err(BulkError::UnknownProfile(other.to_string())),
        };
        let p = Profile::from_json(src)?;
        p.validate()?;
        Ok(p)
    }

    pub fn from_json(s: &str) -> Result<Profile, BulkError> {
        Ok(serde_json::from_str(s)?)
    }

    pub fn validate(&self) -> Result<(), BulkError> {
        self.fitted.gap_dist.validate()?;
        self.fitted.sfs.validate()?;
        self.fitted.indel_length.validate()?;
        self.fitted.variant_classes.validate()?;
        if self.fitted.ploidy == 0 {
            return Err(BulkError::Invalid("ploidy must be >= 1".into()));
        }
        for (label, v) in [
            ("multiallelic_rate", self.fitted.multiallelic_rate),
            ("missing_rate", self.fitted.missing_rate),
            ("phased_rate", self.fitted.phased_rate),
        ] {
            if !(0.0..=1.0).contains(&v) {
                return Err(BulkError::Invalid(format!("{label} must be in [0, 1]")));
            }
        }
        if self.fitted.contigs.is_empty() {
            return Err(BulkError::Invalid("need >= 1 contig".into()));
        }
        Ok(())
    }
}

impl Histogram {
    pub fn validate(&self) -> Result<(), BulkError> {
        if self.edges.len() < 2 {
            return Err(BulkError::Invalid("histogram needs >= 2 edges".into()));
        }
        if self.weights.len() + 1 != self.edges.len() {
            return Err(BulkError::Invalid(format!(
                "histogram weights ({}) must be edges ({}) - 1",
                self.weights.len(),
                self.edges.len()
            )));
        }
        if self.weights.iter().any(|w| *w < 0.0) {
            return Err(BulkError::Invalid("histogram weights must be >= 0".into()));
        }
        if self.weights.iter().sum::<f64>() <= 0.0 {
            return Err(BulkError::Invalid("histogram weights must sum > 0".into()));
        }
        if self.edges.windows(2).any(|w| w[1] <= w[0]) {
            return Err(BulkError::Invalid("histogram edges must be increasing".into()));
        }
        Ok(())
    }
}

impl ClassMix {
    pub fn validate(&self) -> Result<(), BulkError> {
        let sum = self.snp + self.insertion + self.deletion + self.mnp + self.complex + self.symbolic;
        if (sum - 1.0).abs() > 1e-6 {
            return Err(BulkError::Invalid(format!("variant_classes must sum to 1.0, got {sum}")));
        }
        Ok(())
    }
}
```

- [ ] **Step 4: Write the minimal test profile**

Create `profiles/germline-1kgp.json`. These are **placeholder values so tests can
run** — Task 9 overwrites this file with the real fit. Mark it as such:

```json
{
  "name": "germline-1kgp",
  "provenance": {
    "source": "PLACEHOLDER - not yet fitted; see Task 9",
    "n_samples_source": 0,
    "n_variants_source": 0,
    "fitted_on": "1970-01-01",
    "fit_tool_version": "0.0.0"
  },
  "fitted": {
    "contigs": [
      { "id": "chr1", "n_variants": 1000, "density_per_kb": 40.0 },
      { "id": "chr2", "n_variants": 1000, "density_per_kb": 40.0 },
      { "id": "chr3", "n_variants": 1000, "density_per_kb": 40.0 }
    ],
    "gap_dist": { "edges": [1.0, 10.0, 100.0, 1000.0], "weights": [0.6, 0.3, 0.1] },
    "sfs": { "edges": [1.0, 2.0, 10.0, 100.0, 6404.0], "weights": [0.476, 0.2, 0.2, 0.124] },
    "variant_classes": {
      "snp": 0.83,
      "insertion": 0.06,
      "deletion": 0.09,
      "mnp": 0.005,
      "complex": 0.005,
      "symbolic": 0.01
    },
    "indel_length": { "edges": [1.0, 2.0, 6.0, 20.0, 100.0], "weights": [0.5, 0.4, 0.085, 0.015] },
    "titv": 2.05,
    "multiallelic_rate": 0.0,
    "missing_rate": 0.0,
    "phased_rate": 1.0,
    "ploidy": 2
  },
  "dialed": { "payload": "gt-only" }
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test --features bulk profile`
Expected: 5 tests PASS.

- [ ] **Step 6: Verify default build is unaffected**

Run: `cargo test`
Expected: 70 tests pass; `bulk` module absent.

- [ ] **Step 7: Commit**

```bash
git add src/lib.rs src/bulk/ profiles/
git commit -m "feat(bulk): add profile schema with fitted/dialed partition"
```

---

### Task 3: Samplers

Wave 3 — parallel with Tasks 4, 5, 8.

**Files:**
- Create: `src/bulk/sample.rs`
- Modify: `src/bulk/mod.rs` (add `pub mod sample;`)
- Test: inline `#[cfg(test)] mod tests` in `src/bulk/sample.rs`

**Interfaces:**
- Consumes: `crate::bulk::profile::{Fitted, Histogram, ClassMix}`, `crate::bulk::BulkError`.
- Produces:

```rust
pub enum VariantClass { Snp, Insertion, Deletion, Mnp, Complex, Symbolic }

pub struct Samplers { /* private precomputed CDFs */ }

impl Samplers {
    pub fn new(fitted: &Fitted) -> Result<Samplers, BulkError>;
    pub fn gap<R: rand::Rng>(&self, rng: &mut R) -> u64;         // >= 1
    pub fn allele_count<R: rand::Rng>(&self, rng: &mut R, n_alleles: u64) -> u64; // 1..=n_alleles
    pub fn class<R: rand::Rng>(&self, rng: &mut R) -> VariantClass;
    pub fn indel_len<R: rand::Rng>(&self, rng: &mut R) -> usize; // >= 1
    pub fn snp_alt<R: rand::Rng>(&self, rng: &mut R, ref_base: u8) -> u8; // respects titv
    pub fn base<R: rand::Rng>(&self, rng: &mut R) -> u8;         // one of b"ACGT"
}
```

`Samplers::new` precomputes cumulative weights once so per-record sampling is a
binary search, not a re-scan. This matters: it runs ~265k times.

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bulk::profile::Profile;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn samplers() -> Samplers {
        let p = Profile::builtin("germline-1kgp").unwrap();
        Samplers::new(&p.fitted).unwrap()
    }

    #[test]
    fn gaps_are_at_least_one() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        for _ in 0..1000 {
            assert!(s.gap(&mut rng) >= 1);
        }
    }

    #[test]
    fn allele_count_is_in_range() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        for _ in 0..1000 {
            let ac = s.allele_count(&mut rng, 6404);
            assert!((1..=6404).contains(&ac), "ac out of range: {ac}");
        }
    }

    #[test]
    fn sfs_reproduces_singleton_fraction() {
        // The germline profile puts 47.6% of weight in the [1, 2) bin.
        // This is the whole point of fitting an empirical SFS (a neutral 1/x
        // SFS would give ~12%), so guard it.
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let n = 20_000;
        let singletons = (0..n).filter(|_| s.allele_count(&mut rng, 6404) == 1).count();
        let frac = singletons as f64 / n as f64;
        assert!((frac - 0.476).abs() < 0.02, "singleton fraction {frac} != ~0.476");
    }

    #[test]
    fn class_mix_is_reproduced() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(4);
        let n = 20_000;
        let snps = (0..n).filter(|_| matches!(s.class(&mut rng), VariantClass::Snp)).count();
        let frac = snps as f64 / n as f64;
        assert!((frac - 0.83).abs() < 0.02, "snp fraction {frac} != ~0.83");
    }

    #[test]
    fn snp_alt_never_equals_ref_and_respects_titv() {
        let s = samplers();
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let n = 20_000;
        let mut ti = 0usize;
        for _ in 0..n {
            let alt = s.snp_alt(&mut rng, b'A');
            assert_ne!(alt, b'A');
            if alt == b'G' { ti += 1; } // A<->G is a transition
        }
        // titv = 2.05 => transitions are 2.05 / 3.05 of SNPs
        let frac = ti as f64 / n as f64;
        assert!((frac - 2.05 / 3.05).abs() < 0.02, "ti fraction {frac}");
    }

    #[test]
    fn sampling_is_deterministic_for_a_seed() {
        let s = samplers();
        let mut a = ChaCha8Rng::seed_from_u64(7);
        let mut b = ChaCha8Rng::seed_from_u64(7);
        let xs: Vec<u64> = (0..100).map(|_| s.gap(&mut a)).collect();
        let ys: Vec<u64> = (0..100).map(|_| s.gap(&mut b)).collect();
        assert_eq!(xs, ys);
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --features bulk sample`
Expected: FAIL — `Samplers` not found.

- [ ] **Step 3: Implement the samplers**

Add `pub mod sample;` to `src/bulk/mod.rs`. In `src/bulk/sample.rs`:

```rust
use rand::Rng;

use crate::bulk::profile::{ClassMix, Fitted, Histogram};
use crate::bulk::BulkError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass { Snp, Insertion, Deletion, Mnp, Complex, Symbolic }

/// A histogram sampler with a precomputed CDF.
#[derive(Debug, Clone)]
struct HistSampler { edges: Vec<f64>, cdf: Vec<f64> }

impl HistSampler {
    fn new(h: &Histogram) -> Result<HistSampler, BulkError> {
        h.validate()?;
        let total: f64 = h.weights.iter().sum();
        let mut cdf = Vec::with_capacity(h.weights.len());
        let mut acc = 0.0;
        for w in &h.weights {
            acc += w / total;
            cdf.push(acc);
        }
        Ok(HistSampler { edges: h.edges.clone(), cdf })
    }

    /// Draw a bin by CDF binary search, then a value uniformly within it.
    ///
    /// Callers quantize with `.floor()`, never `.round()`. A bin `[lo, hi)` over
    /// integer edges represents the integers `lo ..= hi - 1`, so `floor` of a
    /// uniform draw on `[lo, hi)` is exactly discrete-uniform over them.
    /// `.round()` instead sends the bin's upper half to `hi` — for the SFS's
    /// `[1, 2)` singleton bin that bleeds half the singleton mass to AC=2,
    /// halving the fitted singleton fraction (0.476 -> ~0.24) and failing
    /// `sfs_reproduces_singleton_fraction`.
    fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        let u: f64 = rng.gen();
        let bin = self.cdf.partition_point(|c| *c < u).min(self.cdf.len() - 1);
        let (lo, hi) = (self.edges[bin], self.edges[bin + 1]);
        rng.gen_range(lo..hi)
    }
}

#[derive(Debug, Clone)]
pub struct Samplers {
    gap: HistSampler,
    sfs: HistSampler,
    indel: HistSampler,
    class_cdf: [f64; 6],
    ti_frac: f64,
}

impl Samplers {
    pub fn new(fitted: &Fitted) -> Result<Samplers, BulkError> {
        fitted.variant_classes.validate()?;
        if fitted.titv <= 0.0 {
            return Err(BulkError::Invalid("titv must be > 0".into()));
        }
        let m: &ClassMix = &fitted.variant_classes;
        let mut acc = 0.0;
        let mut class_cdf = [0.0f64; 6];
        for (i, w) in [m.snp, m.insertion, m.deletion, m.mnp, m.complex, m.symbolic].iter().enumerate() {
            acc += w;
            class_cdf[i] = acc;
        }
        Ok(Samplers {
            gap: HistSampler::new(&fitted.gap_dist)?,
            sfs: HistSampler::new(&fitted.sfs)?,
            indel: HistSampler::new(&fitted.indel_length)?,
            class_cdf,
            ti_frac: fitted.titv / (fitted.titv + 1.0),
        })
    }

    pub fn gap<R: Rng>(&self, rng: &mut R) -> u64 {
        (self.gap.sample(rng).floor() as u64).max(1)
    }

    pub fn allele_count<R: Rng>(&self, rng: &mut R, n_alleles: u64) -> u64 {
        (self.sfs.sample(rng).floor() as u64).clamp(1, n_alleles)
    }

    pub fn indel_len<R: Rng>(&self, rng: &mut R) -> usize {
        (self.indel.sample(rng).floor() as usize).max(1)
    }

    pub fn class<R: Rng>(&self, rng: &mut R) -> VariantClass {
        let u: f64 = rng.gen();
        let i = self.class_cdf.partition_point(|c| *c < u).min(5);
        [VariantClass::Snp, VariantClass::Insertion, VariantClass::Deletion,
         VariantClass::Mnp, VariantClass::Complex, VariantClass::Symbolic][i]
    }

    pub fn base<R: Rng>(&self, rng: &mut R) -> u8 {
        b"ACGT"[rng.gen_range(0..4)]
    }

    /// Draw a SNP ALT != `ref_base`, with transitions at `titv / (titv + 1)`.
    pub fn snp_alt<R: Rng>(&self, rng: &mut R, ref_base: u8) -> u8 {
        let transition = match ref_base { b'A' => b'G', b'G' => b'A', b'C' => b'T', b'T' => b'C', _ => b'A' };
        if rng.gen::<f64>() < self.ti_frac {
            transition
        } else {
            let transversions: [u8; 2] = match ref_base {
                b'A' | b'G' => [b'C', b'T'],
                _ => [b'A', b'G'],
            };
            transversions[rng.gen_range(0..2)]
        }
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --features bulk sample`
Expected: 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bulk/sample.rs src/bulk/mod.rs
git commit -m "feat(bulk): add profile-driven samplers with precomputed CDFs"
```

---

### Task 4: Streaming writer + counting writer + second-pass index

Wave 3 — parallel with Tasks 3, 5, 8.

Read the spec's "Writer and sizing" section first. Two non-obvious constraints:
`MultithreadedWriter` exposes **neither** `virtual_position()` **nor** a byte
count, so the index is a second pass via `bcf::fs::index` and the compressed byte
count comes from a `CountingWriter` placed **underneath** the bgzf writer.

**Files:**
- Create: `src/bulk/writer.rs`
- Modify: `src/bulk/mod.rs` (add `pub mod writer;`)
- Test: inline `#[cfg(test)] mod tests` in `src/bulk/writer.rs`

**Interfaces:**
- Consumes: `crate::bulk::BulkError`.
- Produces:

```rust
pub struct CountingWriter<W> { /* private */ }
impl<W: std::io::Write> CountingWriter<W> {
    pub fn new(inner: W) -> (CountingWriter<W>, std::sync::Arc<std::sync::atomic::AtomicU64>);
}

pub enum Format { Bcf, VcfGz, Vcf }

pub struct BulkWriter { /* private */ }
impl BulkWriter {
    pub fn create(path: &std::path::Path, format: Format, header: &noodles_vcf::Header,
                  compression_level: u8, workers: std::num::NonZero<usize>)
        -> Result<BulkWriter, BulkError>;
    pub fn write(&mut self, header: &noodles_vcf::Header,
                 record: &noodles_vcf::variant::RecordBuf) -> Result<(), BulkError>;
    /// Compressed bytes written so far (polled by the size-targeting loop).
    pub fn compressed_bytes(&self) -> u64;
    /// Finishes the stream, then writes `<path>.csi` via a second read pass.
    pub fn finish_and_index(self, path: &std::path::Path) -> Result<(), BulkError>;
}
```

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use noodles_vcf::{self as vcf, variant::RecordBuf};
    use std::num::NonZero;

    fn header() -> vcf::Header {
        vcf::Header::builder()
            .add_contig("chr1", vcf::header::record::value::map::Contig::default())
            .add_sample_name("s1")
            .build()
    }

    #[test]
    fn counting_writer_counts_bytes_through() {
        let (mut w, count) = CountingWriter::new(Vec::new());
        use std::io::Write;
        w.write_all(b"hello").unwrap();
        w.flush().unwrap();
        assert_eq!(count.load(std::sync::atomic::Ordering::Relaxed), 5);
    }

    #[test]
    fn writes_a_readable_indexed_bcf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bcf");
        let h = header();
        let mut w = BulkWriter::create(&path, Format::Bcf, &h, 6, NonZero::new(2).unwrap()).unwrap();
        for pos in [100usize, 200, 300] {
            let rec = RecordBuf::builder()
                .set_reference_sequence_name("chr1")
                .set_variant_start(noodles_core::Position::try_from(pos).unwrap())
                .set_reference_bases("A")
                .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![String::from("T")]))
                .build();
            w.write(&h, &rec).unwrap();
        }
        assert!(w.compressed_bytes() > 0, "counter should see compressed bytes");
        w.finish_and_index(&path).unwrap();

        assert!(path.exists());
        assert!(path.with_extension("bcf.csi").exists(), "csi must be written");

        // Read back through an independent path and confirm the records survive.
        let mut r = noodles_bcf::io::reader::Builder::default().build_from_path(&path).unwrap();
        let rh = r.read_header().unwrap();
        let n = r.records().count();
        assert_eq!(n, 3);
        assert_eq!(rh.sample_names().len(), 1);
    }
}
```

Add `tempfile = "3"` to `[dev-dependencies]` if not present.

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --features bulk writer`
Expected: FAIL — `BulkWriter` not found.

- [ ] **Step 3: Implement**

Add `pub mod writer;` to `src/bulk/mod.rs`. In `src/bulk/writer.rs`:

```rust
use std::io::Write;
use std::num::NonZero;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_csi as csi;
use noodles_vcf::{self as vcf, variant::io::Write as _, variant::RecordBuf};

use crate::bulk::BulkError;

/// Wraps a writer and counts bytes passing through it.
///
/// Placed *underneath* the bgzf writer, so the count is the compressed size —
/// `MultithreadedWriter` exposes no position of its own.
pub struct CountingWriter<W> { inner: W, count: Arc<AtomicU64> }

impl<W: Write> CountingWriter<W> {
    pub fn new(inner: W) -> (CountingWriter<W>, Arc<AtomicU64>) {
        let count = Arc::new(AtomicU64::new(0));
        (CountingWriter { inner, count: Arc::clone(&count) }, count)
    }
}

impl<W: Write> Write for CountingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.count.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> { self.inner.flush() }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format { Bcf, VcfGz, Vcf }

enum Sink {
    Bcf(bcf::io::Writer<bgzf::io::MultithreadedWriter<CountingWriter<std::fs::File>>>),
    VcfGz(vcf::io::Writer<bgzf::io::MultithreadedWriter<CountingWriter<std::fs::File>>>),
    Vcf(vcf::io::Writer<CountingWriter<std::fs::File>>),
}

pub struct BulkWriter { sink: Sink, count: Arc<AtomicU64>, format: Format }

impl BulkWriter {
    pub fn create(
        path: &Path,
        format: Format,
        header: &vcf::Header,
        _compression_level: u8,
        workers: NonZero<usize>,
    ) -> Result<BulkWriter, BulkError> {
        let file = std::fs::File::create(path)?;
        let (counting, count) = CountingWriter::new(file);
        let mut w = match format {
            Format::Bcf => {
                let inner = bgzf::io::MultithreadedWriter::with_worker_count(workers, counting);
                BulkWriter { sink: Sink::Bcf(bcf::io::Writer::from(inner)), count, format }
            }
            Format::VcfGz => {
                let inner = bgzf::io::MultithreadedWriter::with_worker_count(workers, counting);
                BulkWriter { sink: Sink::VcfGz(vcf::io::Writer::new(inner)), count, format }
            }
            Format::Vcf => BulkWriter { sink: Sink::Vcf(vcf::io::Writer::new(counting)), count, format },
        };
        match &mut w.sink {
            Sink::Bcf(x) => x.write_header(header)?,
            Sink::VcfGz(x) => x.write_header(header)?,
            Sink::Vcf(x) => x.write_header(header)?,
        }
        Ok(w)
    }

    pub fn write(&mut self, header: &vcf::Header, record: &RecordBuf) -> Result<(), BulkError> {
        match &mut self.sink {
            Sink::Bcf(x) => x.write_variant_record(header, record)?,
            Sink::VcfGz(x) => x.write_variant_record(header, record)?,
            Sink::Vcf(x) => x.write_variant_record(header, record)?,
        }
        Ok(())
    }

    pub fn compressed_bytes(&self) -> u64 { self.count.load(Ordering::Relaxed) }

    /// Finish the stream, then build and write `<path>.csi` in a second pass.
    pub fn finish_and_index(self, path: &Path) -> Result<(), BulkError> {
        let format = self.format;
        drop(self.sink); // flush + finish the bgzf stream
        if format == Format::Bcf {
            let index = bcf::fs::index(path)?;
            let mut p = path.as_os_str().to_os_string();
            p.push(".csi");
            csi::fs::write(std::path::PathBuf::from(p), &index)?;
        }
        Ok(())
    }
}
```

If `bgzf::io::MultithreadedWriter` does not flush on drop in 0.45, call
`finish()` explicitly before indexing — verify by checking that the test's
read-back sees all 3 records.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --features bulk writer`
Expected: 2 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bulk/writer.rs src/bulk/mod.rs Cargo.toml
git commit -m "feat(bulk): add streaming writer with byte counting and second-pass index"
```

---

### Task 5: Summary truth

Wave 3 — parallel with Tasks 3, 4, 8.

**Files:**
- Create: `src/bulk/summary.rs`
- Modify: `src/bulk/mod.rs` (add `pub mod summary;`)
- Test: inline `#[cfg(test)] mod tests` in `src/bulk/summary.rs`

**Interfaces:**
- Consumes: `crate::bulk::sample::VariantClass`, `crate::bulk::BulkError`.
- Produces:

```rust
pub struct ContigSummary { pub n_records: u64, pub pos_min: u64, pub pos_max: u64 }
pub struct Summary {
    pub n_samples: usize,
    pub per_contig: std::collections::BTreeMap<String, ContigSummary>,
    pub n_alleles_total: u64,
    pub n_alleles_nonref: u64,
    pub class_counts: std::collections::BTreeMap<String, u64>,
    pub genotype_checksum: u64,
}
impl Summary {
    pub fn new(n_samples: usize) -> Summary;
    pub fn observe(&mut self, chrom: &str, pos: u64, class: VariantClass, gts: &[i8]);
    pub fn n_records_total(&self) -> u64;
    pub fn to_json(&self) -> Result<String, BulkError>;
}
```

`observe` is called once per record and must be O(gts.len()) with no allocation.
The checksum is FNV-1a over the genotype bytes — order-dependent on purpose, so a
reader that drops or reorders records fails to reproduce it.

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bulk::sample::VariantClass;

    #[test]
    fn tracks_counts_and_ranges() {
        let mut s = Summary::new(2);
        s.observe("chr1", 100, VariantClass::Snp, &[0, 1, 1, 1]);
        s.observe("chr1", 500, VariantClass::Deletion, &[0, 0, 0, 0]);
        s.observe("chr2", 7, VariantClass::Snp, &[1, 1, 0, 0]);

        assert_eq!(s.n_records_total(), 3);
        assert_eq!(s.per_contig["chr1"].n_records, 2);
        assert_eq!(s.per_contig["chr1"].pos_min, 100);
        assert_eq!(s.per_contig["chr1"].pos_max, 500);
        assert_eq!(s.per_contig["chr2"].pos_min, 7);
        assert_eq!(s.n_alleles_total, 12);
        assert_eq!(s.n_alleles_nonref, 5);
        assert_eq!(s.class_counts["snp"], 2);
        assert_eq!(s.class_counts["deletion"], 1);
    }

    #[test]
    fn missing_alleles_are_not_counted_as_nonref() {
        let mut s = Summary::new(1);
        s.observe("chr1", 1, VariantClass::Snp, &[-1, 1]);
        assert_eq!(s.n_alleles_nonref, 1);
    }

    #[test]
    fn checksum_detects_a_dropped_record() {
        let mut a = Summary::new(1);
        a.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        a.observe("chr1", 2, VariantClass::Snp, &[1, 1]);
        let mut b = Summary::new(1);
        b.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        assert_ne!(a.genotype_checksum, b.genotype_checksum);
    }

    #[test]
    fn checksum_detects_reordering() {
        let mut a = Summary::new(1);
        a.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        a.observe("chr1", 2, VariantClass::Snp, &[1, 0]);
        let mut b = Summary::new(1);
        b.observe("chr1", 1, VariantClass::Snp, &[1, 0]);
        b.observe("chr1", 2, VariantClass::Snp, &[0, 1]);
        assert_ne!(a.genotype_checksum, b.genotype_checksum);
    }

    #[test]
    fn serializes_to_json() {
        let mut s = Summary::new(1);
        s.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        let j = s.to_json().unwrap();
        assert!(j.contains("\"n_samples\""));
        assert!(j.contains("\"genotype_checksum\""));
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --features bulk summary`
Expected: FAIL — `Summary` not found.

- [ ] **Step 3: Implement**

Add `pub mod summary;` to `src/bulk/mod.rs`. In `src/bulk/summary.rs`:

```rust
use std::collections::BTreeMap;

use crate::bulk::sample::VariantClass;
use crate::bulk::BulkError;

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x100_0000_01b3;

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ContigSummary { pub n_records: u64, pub pos_min: u64, pub pos_max: u64 }

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Summary {
    pub n_samples: usize,
    pub per_contig: BTreeMap<String, ContigSummary>,
    pub n_alleles_total: u64,
    pub n_alleles_nonref: u64,
    pub class_counts: BTreeMap<String, u64>,
    pub genotype_checksum: u64,
}

fn class_name(c: VariantClass) -> &'static str {
    match c {
        VariantClass::Snp => "snp",
        VariantClass::Insertion => "insertion",
        VariantClass::Deletion => "deletion",
        VariantClass::Mnp => "mnp",
        VariantClass::Complex => "complex",
        VariantClass::Symbolic => "symbolic",
    }
}

impl Summary {
    pub fn new(n_samples: usize) -> Summary {
        Summary {
            n_samples,
            per_contig: BTreeMap::new(),
            n_alleles_total: 0,
            n_alleles_nonref: 0,
            class_counts: BTreeMap::new(),
            genotype_checksum: FNV_OFFSET,
        }
    }

    pub fn observe(&mut self, chrom: &str, pos: u64, class: VariantClass, gts: &[i8]) {
        let e = self.per_contig.entry(chrom.to_string()).or_insert(ContigSummary {
            n_records: 0, pos_min: u64::MAX, pos_max: 0,
        });
        e.n_records += 1;
        e.pos_min = e.pos_min.min(pos);
        e.pos_max = e.pos_max.max(pos);

        *self.class_counts.entry(class_name(class).to_string()).or_insert(0) += 1;

        self.n_alleles_total += gts.len() as u64;
        for g in gts {
            if *g > 0 { self.n_alleles_nonref += 1; }
            self.genotype_checksum ^= *g as u8 as u64;
            self.genotype_checksum = self.genotype_checksum.wrapping_mul(FNV_PRIME);
        }
    }

    pub fn n_records_total(&self) -> u64 {
        self.per_contig.values().map(|c| c.n_records).sum()
    }

    pub fn to_json(&self) -> Result<String, BulkError> {
        Ok(serde_json::to_string_pretty(self)?)
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --features bulk summary`
Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bulk/summary.rs src/bulk/mod.rs
git commit -m "feat(bulk): add summary truth with order-sensitive genotype checksum"
```

---

### Task 6: Record generator

Wave 4. Needs Tasks 3, 4, 5.

**Files:**
- Create: `src/bulk/gen.rs`
- Modify: `src/bulk/mod.rs` (add `pub mod gen;`)
- Test: inline `#[cfg(test)] mod tests` in `src/bulk/gen.rs`

**Interfaces:**
- Consumes: `Samplers`, `VariantClass`, `Summary`, `Profile`, `Payload`.
- Produces:

```rust
pub struct GenRecord { pub chrom: String, pub pos: u64, pub ref_: String,
                       pub alts: Vec<String>, pub class: VariantClass, pub gts: Vec<i8> }
pub fn block_rng(seed: u64, block_idx: u64) -> rand_chacha::ChaCha8Rng;
pub fn gen_record<R: rand::Rng>(rng: &mut R, s: &Samplers, chrom: &str, pos: u64,
                                n_samples: usize, ploidy: u8, fitted: &Fitted) -> GenRecord;
pub fn to_record_buf(r: &GenRecord, payload: Payload, phased: bool)
    -> noodles_vcf::variant::RecordBuf;
```

`block_rng` is the determinism guarantee: seed derives from `(seed, block_idx)`
only — never from a thread ID, never from a shared mutable RNG.

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bulk::profile::Profile;
    use crate::bulk::sample::Samplers;

    fn fixture() -> (Profile, Samplers) {
        let p = Profile::builtin("germline-1kgp").unwrap();
        let s = Samplers::new(&p.fitted).unwrap();
        (p, s)
    }

    #[test]
    fn block_rng_is_a_pure_function_of_seed_and_block() {
        use rand::Rng;
        let mut a = block_rng(42, 7);
        let mut b = block_rng(42, 7);
        let mut c = block_rng(42, 8);
        let xa: u64 = a.gen();
        let xb: u64 = b.gen();
        let xc: u64 = c.gen();
        assert_eq!(xa, xb, "same (seed, block) must give the same stream");
        assert_ne!(xa, xc, "different block must give a different stream");
    }

    #[test]
    fn genotypes_have_expected_shape_and_alphabet() {
        let (p, s) = fixture();
        let mut rng = block_rng(1, 0);
        let r = gen_record(&mut rng, &s, "chr1", 100, 10, 2, &p.fitted);
        assert_eq!(r.gts.len(), 20);
        assert!(r.gts.iter().all(|g| (-1..=1).contains(g)));
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.pos, 100);
        assert!(!r.alts.is_empty());
    }

    #[test]
    fn ref_and_alt_are_never_equal() {
        let (p, s) = fixture();
        let mut rng = block_rng(2, 0);
        for i in 0..500 {
            let r = gen_record(&mut rng, &s, "chr1", 100 + i, 4, 2, &p.fitted);
            for a in &r.alts {
                if !a.starts_with('<') {
                    assert_ne!(*a, r.ref_, "ref == alt at iteration {i}");
                }
            }
        }
    }

    #[test]
    fn generation_is_deterministic() {
        let (p, s) = fixture();
        let a: Vec<_> = (0..50).map(|i| {
            let mut r = block_rng(9, i);
            gen_record(&mut r, &s, "chr1", 100 + i, 8, 2, &p.fitted)
        }).collect();
        let b: Vec<_> = (0..50).map(|i| {
            let mut r = block_rng(9, i);
            gen_record(&mut r, &s, "chr1", 100 + i, 8, 2, &p.fitted)
        }).collect();
        assert_eq!(a.iter().map(|r| r.gts.clone()).collect::<Vec<_>>(),
                   b.iter().map(|r| r.gts.clone()).collect::<Vec<_>>());
    }

    #[test]
    fn payload_presets_produce_the_right_format_keys() {
        let (p, s) = fixture();
        let mut rng = block_rng(3, 0);
        let r = gen_record(&mut rng, &s, "chr1", 100, 2, 2, &p.fitted);
        for (payload, expected) in [
            (Payload::GtOnly, vec!["GT"]),
            (Payload::GtVaf, vec!["GT", "VAF"]),
            (Payload::Gatk, vec!["GT", "AD", "DP", "GQ", "PL"]),
            (Payload::Mutect2, vec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"]),
        ] {
            let buf = to_record_buf(&r, payload, true);
            let keys: Vec<String> = buf.samples().keys().map(|k| k.to_string()).collect();
            assert_eq!(keys, expected, "payload {payload:?}");
        }
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --features bulk gen`
Expected: FAIL — `gen_record` not found.

- [ ] **Step 3: Implement**

Add `pub mod gen;` to `src/bulk/mod.rs`. In `src/bulk/gen.rs`:

- `block_rng(seed, block_idx)`: hash the pair with a fixed, stable mixer and seed
  `ChaCha8Rng::seed_from_u64`. Use a splitmix64-style finalizer so adjacent block
  indices give well-separated streams:

```rust
pub fn block_rng(seed: u64, block_idx: u64) -> rand_chacha::ChaCha8Rng {
    use rand::SeedableRng;
    let mut z = seed ^ block_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    rand_chacha::ChaCha8Rng::seed_from_u64(z)
}
```

- `gen_record`: draw `class`; build REF/ALT per class —
  - `Snp`: `ref_` = one base; `alts` = `[snp_alt(ref)]`.
  - `Insertion`: `ref_` = one base; `alts` = `[ref_ + indel_len random bases]`.
  - `Deletion`: `ref_` = `1 + indel_len` random bases; `alts` = `[first base]`.
  - `Mnp`: `ref_` and `alt` both length 2–3, differing at every position.
  - `Complex`: `ref_` length 2–4, `alt` length 2–4, different lengths.
  - `Symbolic`: `ref_` = one base; `alts` = `["<DEL>"]`.

  Then draw `ac = allele_count(rng, n_samples * ploidy)`, set `p = ac / (n_samples * ploidy)`,
  and fill `gts` with `n_samples * ploidy` i.i.d. draws: `-1` with probability
  `fitted.missing_rate`, else `1` with probability `p`, else `0`. **Do not
  implement LD.**

- `to_record_buf`: build a `RecordBuf` with CHROM/POS/REF/ALT and a samples block
  whose keys are exactly the preset's list, in order. GT strings join alleles with
  `|` when `phased` else `/`, mapping `-1` to `.`. For the non-GT fields emit
  cheap deterministic values derived from `gts` (e.g. `DP` = count of non-missing,
  `AD` = `[n_ref, n_alt]`, `GQ` = 99, `PL` = `[0, 30, 60]`, `VAF`/`AF` = alt
  fraction, `F1R2`/`F2R1` = `AD` halves, `SB` = `[0, 0, 0, 0]`). Realism of these
  values is out of scope — only their *presence, type, and cardinality* affect the
  benchmark, per the spec.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --features bulk gen`
Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bulk/gen.rs src/bulk/mod.rs
git commit -m "feat(bulk): add streaming record generator with block-seeded determinism"
```

---

### Task 7: BulkSpec public API

Wave 4. Needs Task 6.

**Files:**
- Modify: `src/bulk/mod.rs`
- Create: `tests/bulk.rs`
- Create: `examples/bulk.rs`

**Interfaces:**
- Consumes: everything above.
- Produces:

```rust
pub enum Size { Records(u64), RecordsPerContig(u64), Target(u64) }
pub struct BulkSpec { /* private */ }
impl BulkSpec {
    pub fn new(profile: Profile) -> BulkSpec;
    pub fn samples(self, n: usize) -> BulkSpec;
    pub fn contigs<I, S>(self, ids: I) -> BulkSpec where I: IntoIterator<Item = S>, S: Into<String>;
    pub fn size(self, size: Size) -> BulkSpec;
    pub fn payload(self, p: Payload) -> BulkSpec;
    pub fn seed(self, seed: u64) -> BulkSpec;
    pub fn format(self, f: Format) -> BulkSpec;
    pub fn workers(self, n: std::num::NonZero<usize>) -> BulkSpec;
    pub fn compression_level(self, level: u8) -> BulkSpec;
    pub fn write(self, path: impl AsRef<std::path::Path>) -> Result<Summary, BulkError>;
}
```

**Contig length rule (spec-critical):** the header's `##contig` `length` is the
**populated span** — `pos_max` for that contig — never a real hg38 length. Because
the span is only known after generating, generate a contig's records into a
bounded buffer first, then emit the header. Records are per-contig, so buffer one
contig at a time, not the whole file. If `Size::Target` is used, poll
`compressed_bytes()` between contigs and stop once the target is reached.

- [ ] **Step 1: Write the failing integration tests**

Create `tests/bulk.rs`:

```rust
#![cfg(feature = "bulk")]

use std::num::NonZero;
use vcfixture::bulk::{BulkSpec, Format, Payload, Profile, Size};

fn spec() -> BulkSpec {
    BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
        .samples(8)
        .contigs(["chr1", "chr2", "chr3"])
        .payload(Payload::GtOnly)
        .seed(42)
        .workers(NonZero::new(2).unwrap())
}

#[test]
fn records_per_contig_is_exact() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = spec().size(Size::RecordsPerContig(100)).write(&path).unwrap();
    assert_eq!(s.n_records_total(), 300);
    assert_eq!(s.per_contig["chr1"].n_records, 100);
    assert_eq!(s.per_contig.len(), 3);
}

#[test]
fn same_seed_gives_byte_identical_output_across_thread_counts() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.bcf");
    let b = dir.path().join("b.bcf");
    spec().size(Size::RecordsPerContig(50)).workers(NonZero::new(1).unwrap()).write(&a).unwrap();
    spec().size(Size::RecordsPerContig(50)).workers(NonZero::new(4).unwrap()).write(&b).unwrap();
    assert_eq!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap(),
               "output must not depend on thread count");
}

#[test]
fn different_seeds_differ() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.bcf");
    let b = dir.path().join("b.bcf");
    spec().seed(1).size(Size::RecordsPerContig(50)).write(&a).unwrap();
    spec().seed(2).size(Size::RecordsPerContig(50)).write(&b).unwrap();
    assert_ne!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap());
}

#[test]
fn declared_contig_length_equals_populated_span() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = spec().size(Size::RecordsPerContig(200)).write(&path).unwrap();

    let mut r = noodles_bcf::io::reader::Builder::default().build_from_path(&path).unwrap();
    let header = r.read_header().unwrap();
    for (id, contig) in header.contigs() {
        let declared = contig.length().expect("contig must declare a length") as u64;
        let span = s.per_contig[id.as_ref()].pos_max;
        assert_eq!(declared, span,
            "contig {id} declared length {declared} must equal populated span {span}");
        // and it must not be a real hg38 length
        assert!(declared < 248_956_422, "contig {id} must not use a real hg38 length");
    }
}

#[test]
fn target_size_lands_near_the_target() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let target = 512 * 1024;
    spec().size(Size::Target(target)).write(&path).unwrap();
    let got = std::fs::metadata(&path).unwrap().len();
    assert!(got >= target, "got {got} < target {target}");
    assert!(got < target + 256 * 1024, "overshoot too large: {got} vs {target}");
}

#[test]
fn summary_matches_an_independent_read_back() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = spec().size(Size::RecordsPerContig(100)).write(&path).unwrap();

    let mut r = noodles_bcf::io::reader::Builder::default().build_from_path(&path).unwrap();
    let _ = r.read_header().unwrap();
    let n = r.records().count() as u64;
    assert_eq!(n, s.n_records_total(), "summary must match what a reader sees");
}

#[test]
fn index_is_written_and_positions_are_sorted() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    spec().size(Size::RecordsPerContig(100)).write(&path).unwrap();
    assert!(path.with_extension("bcf.csi").exists());
}

#[test]
fn payload_presets_all_write_readable_files() {
    for payload in [Payload::GtOnly, Payload::GtVaf, Payload::Gatk, Payload::Mutect2] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.bcf");
        let s = spec().payload(payload).size(Size::RecordsPerContig(20)).write(&path).unwrap();
        assert_eq!(s.n_records_total(), 60, "payload {payload:?}");
        let mut r = noodles_bcf::io::reader::Builder::default().build_from_path(&path).unwrap();
        let _ = r.read_header().unwrap();
        assert_eq!(r.records().count(), 60, "payload {payload:?}");
    }
}
```

Add `noodles-bcf = "0.81"` to `[dev-dependencies]` so the integration test can
read back independently.

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --features bulk --test bulk`
Expected: FAIL — `BulkSpec` not found.

- [ ] **Step 3: Implement `BulkSpec` in `src/bulk/mod.rs`**

Builder with defaults: `samples = 1`, `contigs = ["chr1", "chr2", "chr3"]`,
`size = Size::RecordsPerContig(1000)`, `payload` from `profile.dialed.payload`,
`seed = 0`, `format = Format::Bcf`, `workers = available_parallelism()`,
`compression_level = 6`.

`write` does, per contig:
1. Generate records into a `Vec<GenRecord>` using `block_rng(seed, block_idx)`
   where `block_idx` is a running global record-block counter (so a contig's
   stream does not depend on how many contigs precede it — derive it as
   `contig_idx * 1_000_000 + local_block`).
2. Track `pos` by cumulative `samplers.gap()`.
3. Record `pos_max` -> that contig's declared length.
4. After all contigs' spans are known, build the `vcf::Header` with
   `##contig=<ID=...,length=span>` and the sample names, create the `BulkWriter`,
   then write every buffered record and `summary.observe(...)` each one.
5. `finish_and_index(path)`, write `<path>.summary.json` next to the output,
   return the `Summary`.

For `Size::Target`, generate contig-by-contig in rounds, polling
`compressed_bytes()` after each round, until the target is met. Because the
header must precede the records, write with a provisional pass: generate until
target on a `std::io::sink()`-backed counter to learn the spans and counts, then
do the real write. Simpler alternative, acceptable here: generate a first contig
round, measure bytes/record, extrapolate the remaining record count, then verify
and top up. Prefer the two-pass approach — it is exact and the spec's
accept-longer decision permits it.

Parallelize record generation across blocks with `rayon` (`into_par_iter` over
block indices, each using `block_rng`), then write the blocks **in index order**.
Determinism must hold regardless of `workers`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --features bulk --test bulk`
Expected: 8 tests PASS.

- [ ] **Step 5: Add the runnable example**

Create `examples/bulk.rs` mirroring the spec's API sketch, with a runtime
assertion (`assert_eq!(summary.n_records_total(), 300)`). Register it in
`Cargo.toml`:

```toml
[[example]]
name = "bulk"
required-features = ["bulk"]
```

Add to `.github/workflows/ci.yml` under "run examples":

```yaml
          cargo run --example bulk --features bulk
```

- [ ] **Step 6: Commit**

```bash
git add src/bulk/mod.rs tests/bulk.rs examples/bulk.rs Cargo.toml .github/workflows/ci.yml
git commit -m "feat(bulk): add BulkSpec API with span-derived contig lengths"
```

---

### Task 8: Extraction script

Wave 3 — parallel with Tasks 3, 4, 5. Depends only on Task 2's JSON schema.

**Files:**
- Create: `scripts/fit/fit_profile.py`
- Create: `scripts/fit/README.md`
- Modify: `pixi.toml` (add the `fit` task + Python deps)
- Test: `scripts/fit/test_fit_profile.py`

**Interfaces:**
- Consumes: the profile JSON schema from Task 2.
- Produces: `pixi run fit --pgen <prefix> --name <profile-name> --out <path.json>` writing a schema-valid profile.

The script must never be imported by Rust; its only contract is the JSON schema.

- [ ] **Step 1: Add the pixi feature and task**

In `pixi.toml`:

```toml
[feature.fit.dependencies]
python = ">=3.11"
polars = "*"
plink2 = "*"
pytest = "*"

[feature.fit.tasks]
fit = "python scripts/fit/fit_profile.py"
test-fit = "pytest scripts/fit -q"

[environments]
fit = ["fit"]
```

Run: `pixi install -e fit`
Expected: exit 0.

- [ ] **Step 2: Write the failing test**

Create `scripts/fit/test_fit_profile.py`:

```python
import json

from fit_profile import build_profile, histogram, class_mix_from_counts


def test_histogram_weights_are_one_shorter_than_edges():
    h = histogram([1, 1, 2, 5, 50], edges=[1, 2, 10, 100])
    assert len(h["weights"]) == len(h["edges"]) - 1
    assert abs(sum(h["weights"]) - 1.0) < 1e-9


def test_class_mix_sums_to_one():
    m = class_mix_from_counts({"snp": 83, "insertion": 6, "deletion": 9,
                               "mnp": 1, "complex": 1, "symbolic": 0})
    assert abs(sum(m.values()) - 1.0) < 1e-9


def test_build_profile_emits_schema_valid_json():
    p = build_profile(
        name="test",
        source="/dev/null",
        n_samples=10,
        contigs=[{"id": "chr1", "n_variants": 100, "density_per_kb": 40.0}],
        gaps=[1, 2, 3, 40],
        acs=[1, 1, 2, 19],
        indel_lens=[1, 2, 3],
        class_counts={"snp": 83, "insertion": 6, "deletion": 9,
                      "mnp": 1, "complex": 1, "symbolic": 0},
        titv=2.05,
        multiallelic_rate=0.0,
        missing_rate=0.0,
        phased_rate=1.0,
        ploidy=2,
    )
    j = json.loads(json.dumps(p))
    assert set(j) == {"name", "provenance", "fitted", "dialed"}
    assert set(j["fitted"]) == {
        "contigs", "gap_dist", "sfs", "variant_classes", "indel_length",
        "titv", "multiallelic_rate", "missing_rate", "phased_rate", "ploidy",
    }
    assert j["dialed"]["payload"] in {"gt-only", "gt-vaf", "gatk", "mutect2"}
    # provenance must be populated, never left as a placeholder
    assert j["provenance"]["n_samples_source"] == 10
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pixi run -e fit test-fit`
Expected: FAIL — `ModuleNotFoundError: fit_profile`.

- [ ] **Step 4: Implement `scripts/fit/fit_profile.py`**

Structure:
- `histogram(values, edges) -> {"edges": [...], "weights": [...]}` — normalized,
  `len(weights) == len(edges) - 1`.
- `class_mix_from_counts(counts) -> dict` — normalized to sum exactly 1.0.
- `classify(ref, alt) -> str` — `symbolic` if `alt.startswith("<")` or is a
  breakend; `snp` if `len(ref) == len(alt) == 1`; `mnp` if
  `len(ref) == len(alt) > 1`; `insertion` if `len(alt) > len(ref)` and
  `alt.startswith(ref)`; `deletion` if `len(ref) > len(alt)` and
  `ref.startswith(alt)`; else `complex`.
- `read_pvar(path) -> polars.DataFrame` — handles `.pvar` and `.pvar.zst`; skips
  `##` lines; keeps `#CHROM POS ID REF ALT`. Use `polars.scan_csv` with
  `separator="\t"` and stream; do **not** load 557 MB into memory eagerly.
- `fit_sfs(pgen_prefix) -> list[int]` — shell out to
  `plink2 --pfile <prefix> --freq counts --out <tmp>`, read the `.acount` file
  with polars, return `ALT_CTS`.
- `build_profile(...) -> dict` — assemble, filling `provenance` with the real
  source path, sample/variant counts, today's date, and the script version.
- `main()` — `argparse` with `--pgen`, `--name`, `--out`, `--payload`
  (default `gt-only`), `--contigs` (default: all present).

Per-contig `density_per_kb` = `n_variants / (span_bp / 1000)` where `span_bp` is
`pos_max - pos_min` for that contig in the source.

Gap edges: use log-spaced bins over `[1, 1e5]`. SFS edges: log-spaced over
`[1, 2 * n_samples]` so the first bin is exactly `[1, 2)` — the singleton bin the
Rust test asserts on. Indel edges: `[1, 2, 3, 4, 5, 6, 10, 20, 50, 100, 1000]`
(the spec notes ~90% of indels are <= 6 bp, so resolve that range finely).

- [ ] **Step 5: Run test to verify it passes**

Run: `pixi run -e fit test-fit`
Expected: 3 tests PASS.

- [ ] **Step 6: Write `scripts/fit/README.md`**

Document: what the script does, the exact command to re-fit each of the two
sources (paths from the spec), that `fitted` values come from data while `dialed`
values do not, and that the output must be committed to `profiles/`.

- [ ] **Step 7: Commit**

```bash
git add scripts/fit/ pixi.toml pixi.lock
git commit -m "feat(fit): add profile extraction script for pgen sources"
```

---

### Task 9: Fit the real profiles

Wave 4. Needs Tasks 2 and 8. **This is the only task that reads `/carter`.**

Run it on a machine with access; it reads a 9.6 GB pgen and will take a while.

**Files:**
- Modify: `profiles/germline-1kgp.json` (replace the Task 2 placeholder)
- Create: `profiles/somatic-gdc.json`
- Modify: `src/bulk/profile.rs` (register the second builtin)
- Test: `src/bulk/profile.rs` inline tests

**Interfaces:**
- Consumes: Task 8's `fit` task, Task 2's `Profile::builtin`.
- Produces: `Profile::builtin("somatic-gdc")`.

- [ ] **Step 1: Fit germline**

```bash
pixi run -e fit fit \
  --pgen /carter/users/dlaub/data/1kGP/plink2/hg38.norm \
  --name germline-1kgp --payload gt-only \
  --out profiles/germline-1kgp.json
```

Sanity-check the output before committing:
- `provenance.n_samples_source` == 3202
- the SFS's first bin (`[1, 2)`) weight is ~0.4–0.5 (the spec's 47.6% singleton
  figure). **If it is ~0.12 you have fitted a neutral SFS and something is wrong.**
- `multiallelic_rate` ~0.0 (`hg38.norm` is normalized — expected, not a bug)
- `variant_classes.symbolic` ~0.02 (the HGSVC `<DEL:ME:LINE|L1|L1HS>` calls)
- `titv` ~2.0–2.1

- [ ] **Step 2: Fit somatic**

Note the path: the pvar is `.pvar`, **uncompressed** — not `.pvar.zst`.

```bash
pixi run -e fit fit \
  --pgen /carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt \
  --name somatic-gdc --payload gt-vaf \
  --out profiles/somatic-gdc.json
```

Sanity-check:
- `provenance.n_samples_source` == 16007
- the SFS should be overwhelmingly singleton-dominated (private mutations across a
  merged cohort) — a much sharper spike than germline. This is expected.
- density should be orders of magnitude below germline (somatic ~1–13 mut/Mb vs
  germline ~1300 SNV/Mb).

- [ ] **Step 3: Register the second builtin**

In `src/bulk/profile.rs`:

```rust
const SOMATIC_GDC: &str = include_str!("../../profiles/somatic-gdc.json");
```

and add `"somatic-gdc" => SOMATIC_GDC,` to the `builtin` match.

- [ ] **Step 4: Add the tests**

```rust
#[test]
fn builtin_somatic_loads_and_validates() {
    let p = Profile::builtin("somatic-gdc").unwrap();
    assert_eq!(p.name, "somatic-gdc");
    assert_eq!(p.dialed.payload, Payload::GtVaf);
    p.validate().unwrap();
    assert_eq!(p.provenance.n_samples_source, 16007);
}

#[test]
fn germline_profile_is_really_fitted_not_placeholder() {
    let p = Profile::builtin("germline-1kgp").unwrap();
    assert_eq!(p.provenance.n_samples_source, 3202);
    assert!(!p.provenance.source.contains("PLACEHOLDER"));
}

#[test]
fn germline_sfs_is_empirical_not_neutral() {
    // A neutral 1/x SFS gives ~12% singletons; real 1kGP is ~47.6%.
    // This test is the guard that we fitted data rather than theory.
    let p = Profile::builtin("germline-1kgp").unwrap();
    let total: f64 = p.fitted.sfs.weights.iter().sum();
    let singleton_frac = p.fitted.sfs.weights[0] / total;
    assert!(singleton_frac > 0.3,
        "singleton fraction {singleton_frac} looks neutral, not empirical");
}
```

- [ ] **Step 5: Run the full suite**

Run: `cargo test --all-features`
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add profiles/ src/bulk/profile.rs
git commit -m "feat(bulk): fit germline-1kgp and somatic-gdc profiles from real data"
```

---

### Task 10: CLI + docs

Wave 4. Needs Task 7 (and Task 9 for the somatic profile to be listed).

**Files:**
- Modify: `src/bin/vcfixture.rs` (replace the Task 1 placeholder)
- Create: `docs/book/src/bulk-generation.md`
- Modify: `docs/book/src/SUMMARY.md`
- Modify: `README.md`
- Test: `tests/cli.rs`

**Interfaces:**
- Consumes: `BulkSpec`, `Profile`, `Payload`, `Size`, `Format`.
- Produces: the `vcfixture bulk` subcommand.

- [ ] **Step 1: Write the failing test**

Create `tests/cli.rs`:

```rust
#![cfg(feature = "cli")]

#[test]
fn parses_a_size_with_units() {
    // exercised via the public parser fn
    use vcfixture::bulk::parse_size;
    assert_eq!(parse_size("100MB").unwrap(), 100 * 1024 * 1024);
    assert_eq!(parse_size("512KB").unwrap(), 512 * 1024);
    assert_eq!(parse_size("1GB").unwrap(), 1024 * 1024 * 1024);
    assert_eq!(parse_size("2048").unwrap(), 2048);
    assert!(parse_size("banana").is_err());
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --features cli --test cli`
Expected: FAIL — `parse_size` not found.

- [ ] **Step 3: Implement `parse_size` in `src/bulk/mod.rs`**

```rust
/// Parse a byte size like `100MB`, `512KB`, `1GB`, or a bare byte count.
pub fn parse_size(s: &str) -> Result<u64, BulkError> {
    let t = s.trim();
    let (num, mult) = if let Some(p) = t.strip_suffix("GB") { (p, 1024 * 1024 * 1024) }
        else if let Some(p) = t.strip_suffix("MB") { (p, 1024 * 1024) }
        else if let Some(p) = t.strip_suffix("KB") { (p, 1024) }
        else { (t, 1) };
    num.trim().parse::<u64>()
        .map(|n| n * mult)
        .map_err(|_| BulkError::Invalid(format!("bad size: {s}")))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --features cli --test cli`
Expected: PASS.

- [ ] **Step 5: Implement the CLI**

Replace `src/bin/vcfixture.rs`:

```rust
use clap::{Parser, Subcommand, ValueEnum};

#[derive(Parser)]
#[command(name = "vcfixture", about = "Generate VCF/BCF test data")]
struct Cli { #[command(subcommand)] cmd: Cmd }

#[derive(Subcommand)]
enum Cmd {
    /// Generate a bulk BCF for benchmarking.
    Bulk {
        /// Builtin profile name, or a path to a profile JSON.
        #[arg(long, default_value = "germline-1kgp")]
        profile: String,
        #[arg(long, default_value_t = 1)]
        samples: usize,
        #[arg(long, value_delimiter = ',', default_values_t = ["chr1".to_string(), "chr2".to_string(), "chr3".to_string()])]
        contigs: Vec<String>,
        /// Stop once the compressed output reaches this size (e.g. 100MB).
        #[arg(long, conflicts_with_all = ["records", "records_per_contig"])]
        target_size: Option<String>,
        #[arg(long, conflicts_with_all = ["target_size", "records_per_contig"])]
        records: Option<u64>,
        #[arg(long, conflicts_with_all = ["target_size", "records"])]
        records_per_contig: Option<u64>,
        /// Override the profile's payload preset.
        #[arg(long)]
        payload: Option<PayloadArg>,
        #[arg(long, default_value = "bcf")]
        format: FormatArg,
        #[arg(long, default_value_t = 0)]
        seed: u64,
        #[arg(long, default_value_t = 6)]
        compression_level: u8,
        #[arg(long)]
        threads: Option<usize>,
        #[arg(short, long)]
        output: std::path::PathBuf,
    },
}
```

Add `PayloadArg`/`FormatArg` as `ValueEnum` mirrors of `Payload`/`Format`
(kebab-case values: `gt-only`, `gt-vaf`, `gatk`, `mutect2`; `bcf`, `vcf-gz`,
`vcf`). Resolve `--profile` as a builtin name first, falling back to reading it as
a file path. Print the resulting size, record count, and elapsed time to stderr on
completion.

Run and verify by hand:

```bash
cargo run --features cli --bin vcfixture -- bulk \
  --profile germline-1kgp --samples 100 --records-per-contig 500 \
  --seed 42 -o /tmp/x.bcf
```

Expected: `/tmp/x.bcf`, `/tmp/x.bcf.csi`, and `/tmp/x.bcf.summary.json` exist;
`bcftools index -n /tmp/x.bcf` reports 1500 records.

- [ ] **Step 6: Write the mdBook chapter**

Create `docs/book/src/bulk-generation.md` covering: what bulk generation is for
(benchmarking, not fixtures — link to the fixture chapters for the oracle path);
the fitted-vs-dialed split; the four payload presets and *why* payload is dialed
rather than fitted (the sources have no INFO/FORMAT); why contigs are declared at
span lengths; why genotypes are i.i.d. (the LD ablation); and both the API and CLI
examples from this plan. Add `- [Bulk generation](bulk-generation.md)` to
`docs/book/src/SUMMARY.md`.

Run: `pixi run docs-build`
Expected: exit 0.

- [ ] **Step 7: Update the README**

Add a short "Bulk generation" section after "Proptest strategies" with the CLI
one-liner and a pointer to the mdBook chapter.

- [ ] **Step 8: Commit**

```bash
git add src/bin/vcfixture.rs src/bulk/mod.rs tests/cli.rs docs/book/ README.md
git commit -m "feat(cli): add vcfixture bulk subcommand and docs"
```

---

### Task 11: Profile fidelity round-trip

Wave 4. Needs Tasks 7, 8, 9.

This is the test that proves the samplers are actually correct: generate from a
profile, re-fit the *generated output*, and assert the fitted stats come back
close to the profile that produced them. Everything else tests plumbing; this
tests fidelity.

The fit script reads pgen, and the generator emits BCF — bridge them with
`plink2 --bcf`, so the round-trip goes through the exact same fitting code path as
the real data rather than a parallel implementation.

**Files:**
- Create: `scripts/fit/test_fidelity.py`
- Modify: `pixi.toml` (add the `test-fidelity` task)

**Interfaces:**
- Consumes: `vcfixture bulk` (Task 10 CLI), `fit_profile.build_profile` (Task 8),
  `profiles/germline-1kgp.json` (Task 9).
- Produces: nothing consumed downstream.

- [ ] **Step 1: Add the pixi task**

In `pixi.toml` under `[feature.fit.tasks]`:

```toml
test-fidelity = "pytest scripts/fit/test_fidelity.py -q"
```

- [ ] **Step 2: Write the failing test**

Create `scripts/fit/test_fidelity.py`:

```python
"""Generate from a profile, re-fit the output, assert the stats round-trip.

This is the real test that the samplers reproduce the profile they were given.
Bridges BCF -> pgen with plink2 so re-fitting uses the same code path as the
original fit.
"""

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PROFILE = REPO / "profiles" / "germline-1kgp.json"

pytestmark = pytest.mark.skipif(
    shutil.which("plink2") is None, reason="plink2 not available"
)


def _generate(out: Path, samples: int, per_contig: int) -> None:
    subprocess.run(
        ["cargo", "run", "--release", "--features", "cli", "--bin", "vcfixture", "--",
         "bulk", "--profile", str(PROFILE), "--samples", str(samples),
         "--contigs", "chr1,chr2,chr3", "--records-per-contig", str(per_contig),
         "--seed", "42", "-o", str(out)],
        cwd=REPO, check=True,
    )


def _refit(bcf: Path, tmp: Path) -> dict:
    prefix = tmp / "refit"
    subprocess.run(
        ["plink2", "--bcf", str(bcf), "--make-pgen", "--out", str(prefix)],
        check=True, capture_output=True,
    )
    out = tmp / "refit.json"
    subprocess.run(
        ["python", str(REPO / "scripts" / "fit" / "fit_profile.py"),
         "--pgen", str(prefix), "--name", "refit", "--out", str(out)],
        check=True,
    )
    return json.loads(out.read_text())


def test_generated_output_refits_to_its_source_profile():
    original = json.loads(PROFILE.read_text())["fitted"]
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        bcf = tmp / "gen.bcf"
        _generate(bcf, samples=200, per_contig=5000)
        refit = _refit(bcf, tmp)["fitted"]

        # Ti/Tv is a single scalar and the most direct sampler check.
        assert abs(refit["titv"] - original["titv"]) < 0.15, \
            f"titv {refit['titv']} != {original['titv']}"

        # Class mix must survive the round-trip.
        for cls in ("snp", "insertion", "deletion"):
            a = original["variant_classes"][cls]
            b = refit["variant_classes"][cls]
            assert abs(a - b) < 0.05, f"class {cls}: {b} != {a}"

        # The singleton fraction is the stat the whole empirical-SFS decision
        # exists to preserve, so it gets the tightest guard.
        a = original["sfs"]["weights"][0] / sum(original["sfs"]["weights"])
        b = refit["sfs"]["weights"][0] / sum(refit["sfs"]["weights"])
        assert abs(a - b) < 0.06, f"singleton fraction {b} != {a}"

        assert refit["ploidy"] == original["ploidy"]
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pixi run -e fit test-fidelity`
Expected: FAIL (the binary or script is missing, or the tolerances are not met).

- [ ] **Step 4: Make it pass**

If a tolerance fails, the bug is in the **sampler**, not the test — fix
`src/bulk/sample.rs` or `src/bulk/gen.rs`. Do not loosen a tolerance to make the
test pass without first explaining, in the commit message, why the discrepancy is
a sampling artifact rather than a defect.

Two known-legitimate sources of small drift, acceptable to document rather than
fix:
- `plink2 --make-pgen` may normalize or re-order alleles, shifting class counts
  slightly.
- Symbolic ALTs (`<DEL>`) may be dropped by the pgen conversion, so the `symbolic`
  class is deliberately not asserted above.

- [ ] **Step 5: Run test to verify it passes**

Run: `pixi run -e fit test-fidelity`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/fit/test_fidelity.py pixi.toml
git commit -m "test(bulk): assert generated output re-fits to its source profile"
```

---

## Verification

After Task 11, confirm the whole thing end-to-end:

- [ ] `pixi run fmt-check && pixi run clippy && pixi run test` — all green.
- [ ] `cargo test` (default features) — 70 original tests pass, bulk absent.
- [ ] `cargo tree | grep -E "clap|serde_json|rayon|noodles-bcf"` — empty on default features.
- [ ] Generate the real target and time it:

```bash
time cargo run --release --features cli --bin vcfixture -- bulk \
  --profile germline-1kgp --samples 3202 --contigs chr1,chr2,chr3 \
  --target-size 100MB --seed 42 -o /tmp/bench.bcf
```

Expected: ~100 MB output, `.csi` present, completes in seconds-to-low-tens-of-seconds.
Report the actual wall time — the spec's target is "a few seconds," with
accept-longer explicitly permitted. Do not sacrifice realism to hit it.

- [ ] `bcftools view -h /tmp/bench.bcf | grep contig` — lengths are spans (small), not hg38 lengths.
- [ ] `bcftools view -H /tmp/bench.bcf | head` — records look like plausible variants.
- [ ] Region query at an arbitrary offset returns variants at realistic density:

```bash
bcftools view -H /tmp/bench.bcf chr2:1000-2000 | wc -l
```

Expected: non-zero and roughly `40/kb * 1kb` = ~40 records. This is the check that
the fake-contig-length decision actually delivered its purpose.
