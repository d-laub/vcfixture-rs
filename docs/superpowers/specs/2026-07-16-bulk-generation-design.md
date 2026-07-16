# Bulk VCF generation for benchmarking

## Problem

Downstream libraries (`genoray`, `GenVarLoader`) need sizable, realistic variant
files to benchmark read speed and memory against. Today there is no way to make
them: `vcfixture`'s `Document`/`VcfBuilder`/`GroundTruth` path is built for small,
precise fixtures with a decoded oracle, and cannot scale — a `Document` holds
every `Record` in memory with a `Vec<SampleValues>` of `IndexMap`s per record, and
`render()` returns the whole file as one `String`.

Hand-cut slices of real data are the status quo alternative. They are unshareable
(controlled access), unparameterizable (you get the samples and density the slice
happens to have), and cannot sweep an axis.

This adds a **bulk generation** path: fit summary statistics from real local data
once, commit the resulting profile, and generate realistic-enough BCFs at scale
from it. "Realistic enough" is scoped precisely: **fidelity on the statistics that
demonstrably move read speed, memory, and compression**, and nothing else.

## Goals

- Generate ~100 MB (compressed) BCFs with >=3 contigs. Target a few seconds;
  accept longer if realism would otherwise suffer.
- Rust API + CLI.
- Statistics fitted from real local data, not invented.
- Re-fitting new data is a first-class, repeatable workflow.
- Reproducible: same seed and profile produce byte-identical output.

## Non-goals

- A per-genotype ground-truth oracle at bulk scale (see "Summary truth").
- Biological realism: no coalescent, no mutation-rate model, no phylogeny.
- LD / haplotype correlation structure (see "Evidence" — not worth the budget).
- Replacing or modifying the existing fixture path.

## Evidence

Two findings from a literature review plus a controlled ablation drive this
design. Both invert what one might assume from the genotype-compression
literature.

### LD is a ~0x lever on read speed and memory

Ablation: msprime coalescent, 1000 diploids, 5 Mb, 23,747 sites. Sample labels
were permuted independently at each variant — this preserves every site's allele
count exactly, so SFS, genotype sparsity, variant count, and raw byte size are all
held fixed; **only LD is destroyed**. Both raw VCFs were byte-identical at
95,652,697 bytes.

| Metric                           | LD (real) | LD destroyed | ratio     |
| -------------------------------- | --------- | ------------ | --------- |
| bgzip                            | 3.787 MB  | 3.957 MB     | 1.045x    |
| **BCF**                          | 2.817 MB  | 3.204 MB     | **1.14x** |
| pgen                             | 2.319 MB  | 2.507 MB     | 1.08x     |
| xz -6                            | 1.447 MB  | 2.965 MB     | 2.05x     |
| decode `bcftools view -H` (.bcf) | 0.77 s    | 0.78 s       | **1.00**  |
| peak RSS (.bcf)                  | 8,636 KB  | 8,636 KB     | **1.000** |

LD only pays where the compression window can span rows — `xz` sees it, BCF's
bgzip window (32 KB) essentially cannot. **Consequence: i.i.d. genotype draws do
not produce unrealistically incompressible BCF.** Li-Stephens haplotype copying,
block resampling, and coalescent simulation are all budget spent on a 1.14x
compression effect and a 0x speed effect. We draw i.i.d. from HWE.

Caveat: this holds for parse-bound readers. PBWT-decoding formats (savvy, xSI,
GTShark) may couple decode cost to compressed size, where LD could leak into read
speed. Untested — flagged as an open question, not a blocker, since our targets
(genoray, GVL, bcftools) are parse-bound.

### Per-record FORMAT payload is the dominant lever

Same sites, same genotypes, same LD; only adding GATK-style `GT:AD:DP:GQ:PL`:

| Metric        | GT only | GT:AD:DP:GQ:PL | ratio    |
| ------------- | ------- | -------------- | -------- |
| raw           | 95.7 MB | 643.0 MB       | 6.7x     |
| **BCF**       | 2.82 MB | 185.2 MB       | **66x**  |
| decode (.bcf) | 0.76 s  | 4.07 s         | **5.4x** |

The vcf-zarr authors concede the same point in print: their GT-only simulations
are "something of a best-case scenario for specialised genotype compression
methods" (Czech et al., GigaScience 2025).

**But this stat cannot be fitted from our sources.** 1kGP `hg38.norm.pvar` is
`#CHROM POS ID REF ALT` with no INFO column, and pgen stores only genotypes. The
somatic BCF has zero `##INFO` lines and exactly two FORMAT fields (`GT`, `VAF`).
Therefore payload is a **dialed** axis, explicitly marked as such, never presented
as measured.

### The SFS must be empirical, not neutral

A neutral constant-Ne coalescent yields a 12.3% singleton fraction, matching the
Watterson 1/x prediction. Real 1kGP high-coverage is **47.6%** — a ~4x
under-count, driven by recent human population growth. Rare-variant sparsity is
the entire premise of pgen/savvy/spVCF, so a neutral-SFS generator would
systematically flatter sparse formats. Fitting the empirical SFS from the local
pgen is the single highest-value thing the extraction step does.

### Prior art and the gap

`bcftools +simulate` does not exist. bcftools ships the *estimation* half
(`+counts`, `+af-dist`, `+allele-length`, `+indel-stats`) with no emitter.
RAREsim2 and HAPNEST genuinely close the estimate->emit loop but fit only AFS
(+LD/kinship), are germline/biobank-shaped only, and do not target VCF/BCF output.
Every genotype-format paper (savvy, xSI, GTC, GTShark, GSC, SeqArray, BGT, GQT,
vcf-zarr) benchmarks on subsampled real data or GT-only msprime output. Nothing
models FORMAT payload — the stat that actually governs read speed.

## Design

### Module layout

Bulk shares only `spec/` (Number/Type/reserved) and the allele model with the
existing path. It does not touch `Document`, `VcfBuilder`, `truth.rs`, or
`strategies.rs`.

```
src/bulk/
  mod.rs       - public API: BulkSpec, generate()
  profile.rs   - Profile struct + serde; built-in profiles via include_str!
  sample.rs    - samplers: SFS, gap, variant class, indel length
  gen.rs       - streaming record generator
  writer.rs    - streaming BCF/VCF writer + on-the-fly index
  summary.rs   - summary truth
src/bin/vcfixture.rs   - CLI (clap)
profiles/{germline-1kgp,somatic-gdc}.json
scripts/fit/           - Python extraction (pixi task)
```

Feature gates keep the crate light for existing fixture consumers:

- `bulk` (default off) — the bulk module; pulls `serde`, `serde_json`,
  `noodles-bcf`, PRNG.
- `cli` (default off, implies `bulk`) — the binary; pulls `clap`.

### Profile schema

One JSON file per profile, partitioned so a measured stat can never be mistaken
for a chosen one.

**The values below are illustrative schema examples, not fitted results.** Real
values are produced by the extraction step; nothing here should be copied into a
committed profile by hand.

```json
{
  "name": "germline-1kgp",
  "provenance": {
    "source": "/carter/users/dlaub/data/1kGP/plink2/hg38.norm.{pgen,pvar.zst}",
    "n_samples_source": 3202,
    "n_variants_source": 125000000,
    "fitted_on": "2026-07-16",
    "fit_tool_version": "0.2.0"
  },
  "fitted": {
    "contigs": [{ "id": "chr1", "n_variants": 10000000, "density_per_kb": 40.1 }],
    "gap_dist": { "edges": [], "weights": [] },
    "sfs": { "ac_edges": [], "weights": [] },
    "variant_classes": {
      "snp": 0.83,
      "insertion": 0.06,
      "deletion": 0.09,
      "mnp": 0.005,
      "complex": 0.005,
      "symbolic": 0.02
    },
    "indel_length": { "edges": [], "weights": [] },
    "titv": 2.05,
    "multiallelic_rate": 0.0,
    "missing_rate": 0.0,
    "phased_rate": 1.0,
    "ploidy": 2
  },
  "dialed": {
    "payload": "gt-only"
  }
}
```

Notes:

- **`gap_dist` over mean density.** One histogram of inter-variant gaps captures
  clustering, and sampling it (draw gaps, cumsum) yields sorted positions for
  free. Strictly cheaper and strictly more faithful than a scalar density.
- **`sfs` is an allele-count histogram**, not a frequency model. This is the fix
  for the singleton under-count. For the somatic profile it will be
  singleton-dominated (private mutations across a 16,007-sample merged cohort) —
  which falls out of the empirical fit automatically rather than needing a
  separate somatic model.
- **`multiallelic_rate` is expected ~0 for germline.** `hg38.norm` is normalized,
  so multiallelics are split. This legitimately differs from the published 6.6%
  figure for 1kGP high-coverage, and is captured rather than assumed.
- Symbolic alleles are ~2% of the 1kGP pvar (`<DEL:ME:LINE|L1|L1HS>` HGSVC calls);
  `vcfixture` already models these.

### Payload presets

Four, all `dialed`:

| preset    | FORMAT                     | rationale                                |
| --------- | -------------------------- | ---------------------------------------- |
| `gt-only` | `GT`                       | 1kGP-faithful; germline default          |
| `gt-vaf`  | `GT:VAF`                   | somatic-source-faithful; somatic default |
| `gatk`    | `GT:AD:DP:GQ:PL`           | sweeps the 5.4x read-speed axis          |
| `mutect2` | `GT:AD:AF:DP:F1R2:F2R1:SB` | matches the real upstream somatic caller |

Each profile defaults to whatever its source actually carried; the others are
available to sweep the dominant axis deliberately.

### Generation

Records stream into the writer; nothing accumulates. Per record: draw a gap from
`gap_dist` and advance POS; draw a class from `variant_classes`; draw REF/ALT
(indel length from `indel_length`, SNP ALT respecting `titv`); draw an allele
count from `sfs`; draw genotypes i.i.d. from HWE at the implied frequency; apply
`missing_rate` and `phased_rate`; emit payload per the preset.

**Contigs are declared at fake lengths equal to the populated span.** A 100 MB
BCF at 3202 samples is ~265k records; 1kGP density is ~40 variants/kb, so those
records span only a few Mb. Declaring real hg38 lengths and populating a prefix
would make any query outside the prefix return nothing — the exact pathological
case a benchmark trips over. Declaring `length` = generated span means density is
realistic *everywhere in the declared contig*, so **arbitrary region queries always
yield realistic-density data**. The span is an outcome of (records, fitted
density), not a knob.

**Determinism under parallelism.** Each record block seeds its own PRNG from
`hash(seed, record_idx)`, so output is byte-identical regardless of thread count
or scheduling. Generation parallelizes with rayon over record blocks.

### Writer and sizing

Streaming BCF + CSI. The index is built from live bgzf virtual offsets as records
are written, replacing the existing `write_csi` approach (re-render the document,
replay it through an in-memory bgzf writer) which is O(size) twice and cannot
survive bulk scale. `noodles-bgzf`'s `MultithreadedWriter` does the compression.

Output defaults to `.bcf` + `.csi`. `--format` also accepts `vcf.gz` and `vcf`.

Two sizing modes:

- `--target-size 100MB` — stream until the compressed offset reaches the target.
  Exact, not estimated: compressed bytes/record depends on allele frequency, so it
  cannot be predicted in closed form, but it can be *observed* while writing.
- `--records N` / `--records-per-contig N` — emit exact counts.

`--samples` is always explicit. Size targeting solves for records with samples
pinned.

### Summary truth

No per-genotype oracle: a `GroundTruth` genotype array for 3202 samples x 265k
records x ploidy 2 is ~1.3 GB as `i32`, larger than the BCF it describes.

Instead, computed for free while streaming and written to `<out>.summary.json` (or
returned from the API):

- per-contig record counts and POS ranges
- total and non-ref allele counts
- variant-class counts
- a checksum over the genotype stream

This guards the failure mode that makes benchmarks lie: a library that silently
drops records looks *fast*. It keeps faith with the library's premise (never
assert against hand-coded literals) at ~0 cost.

### Extraction script

`scripts/fit/`, Python, run as a `pixi` task:

```bash
pixi run fit --pgen /path/to/data --out profiles/name.json
```

`plink2 --freq` / `--geno-counts` for the SFS; polars over the pvar for classes,
gaps, and Ti/Tv. Writes the profile JSON with a populated `provenance` block.
Re-fitting new data is the same one command.

Sources for the two committed profiles:

- germline: `/carter/users/dlaub/data/1kGP/plink2/hg38.norm.{pgen,psam,pvar.zst}`
  (3202 samples)
- somatic:
  `/carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt.{pgen,psam,pvar}`
  (16,007 samples; note `.pvar` is uncompressed, not `.pvar.zst`)

CI never reads `/carter` — it consumes the committed profiles.

## API sketch

```rust
use vcfixture::bulk::{BulkSpec, Profile, Payload, Size};

let profile = Profile::builtin("germline-1kgp")?;
let summary = BulkSpec::new(profile)
    .samples(3202)
    .contigs(["chr1", "chr2", "chr3"])
    .size(Size::Target(100 * 1024 * 1024))
    .payload(Payload::GtOnly)
    .seed(42)
    .write("bench.bcf")?;
```

```bash
vcfixture bulk --profile germline-1kgp --samples 3202 \
  --contigs chr1,chr2,chr3 --target-size 100MB --seed 42 -o bench.bcf
```

## Risks

**Throughput is the engineering constraint.** 100 MB compressed BCF is ~1.7 GB raw
(~1.7e9 genotype bytes at 3202 samples). A few seconds requires ~500 MB/s of bgzf
compression — `MultithreadedWriter` plus rayon over record blocks, on ~10+ cores.
Per explicit decision: aim for a few seconds, but **accept longer generation time
rather than compromise realism**. Compression level is the available lever if
needed and is exposed as `--compression-level`.

**Profiles can go stale.** They are committed artifacts fitted from data that may
change. Mitigated by `provenance` (source path, counts, date, tool version) and by
making re-fit a single command.

**PBWT-format read speed under destroyed LD is untested.** If a future target is
savvy/xSI rather than genoray/GVL, the LD conclusion needs re-validation.

## Testing

- Round-trip: generate small, read back with `noodles`, assert the summary matches.
- Determinism: same seed => byte-identical output; thread count does not change
  output.
- Profile fidelity: generate from a profile, re-fit the output, assert fitted stats
  land within tolerance of the input profile. This closes the loop and is the real
  test that the samplers are correct.
- Sizing: `--target-size` lands within a bgzf block of target; `--records` exact.
- Contig spans: declared length equals populated span; a region query at any offset
  in the declared contig returns variants at fitted density.
- Existing fixture-path tests must be unaffected (bulk is feature-gated).
