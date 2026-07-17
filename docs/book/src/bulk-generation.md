# Bulk generation

`vcfixture::bulk` (feature `bulk`, plus `cli` for the `vcfixture bulk`
subcommand) generates large, realistic-enough BCF/VCF files for
**benchmarking** — measuring a reader's speed, memory, or compression at
scale. It is a different tool from the rest of this guide: `VcfBuilder` and
`Document` build small, exact fixtures with a decoded per-genotype
[ground-truth oracle](ground-truth.md), and hold every record in memory to do
it. Bulk generation streams records and keeps only summary statistics, so it
scales to files far larger than would fit as a `Document` — but it gives up
the oracle to do so.

Use `VcfBuilder` when a test needs to assert against exact expected values.
Use `vcfixture::bulk` when a benchmark needs a big, plausible file and does
not care what's in any individual record — only that the file behaves like
real data at the sizes and shapes that matter.

## Fitted vs. dialed

Every [`Profile`](https://docs.rs/vcfixture) is split into two parts that are
never allowed to blur together:

- **`fitted`** — statistics estimated from a real cohort (inter-variant gap
  distribution, site-frequency spectrum, variant-class mix, indel lengths,
  Ti/Tv, multiallelic/missing/phased rates, ploidy). Every number here traces
  back to a fit against real data; nothing is hand-picked.
- **`dialed`** — a payload preset chosen explicitly by the profile (or
  overridden by the caller), independent of any fit.

The payload is dialed, not fitted, for a concrete reason: the sources this
crate fits from are genotype-only. The 1000 Genomes `pvar` used for
`germline-1kgp` has no INFO column and no per-sample FORMAT beyond the
genotype call itself, so there is nothing to fit a `GT:AD:DP:GQ:PL`-style
payload *from*. Presenting a dialed value as fitted would misrepresent it as
measured when it is actually chosen — so this guide (and the crate's docs)
never do that.

## Payload presets

Four presets, selectable via `Payload` (Rust) or `--payload` (CLI,
kebab-case):

| preset    | FORMAT                     | why it exists                                    |
| --------- | --------------------------- | ------------------------------------------------- |
| `gt-only` | `GT`                        | matches the 1kGP source; the germline default      |
| `gt-vaf`  | `GT:VAF`                    | matches a genotype+VAF somatic source               |
| `gatk`    | `GT:AD:DP:GQ:PL`            | sweeps per-record FORMAT payload, the dominant lever on read speed and compressed size |
| `mutect2` | `GT:AD:AF:DP:F1R2:F2R1:SB`  | matches a real upstream somatic variant caller      |

Payload size matters more than almost anything else you can dial: in a
controlled ablation, adding `GT:AD:DP:GQ:PL` to otherwise-identical sites
and genotypes made the compressed BCF **~66x** larger and decode time
**~5.4x** slower. If a benchmark cares about payload-heavy callers (GATK,
Mutect2), pick the matching preset explicitly — don't assume `gt-only` is
representative.

## Why contigs are declared at fake lengths

The output header's `##contig` `length` is **never** a real chromosome
length. It is set to the *populated span* of whatever was actually generated
for that contig — the position of the last record written, not (say) hg38's
`chr1` length of ~248 Mb.

This is deliberate. A benchmark-scale BCF only populates a prefix of a real
chromosome (a 100 MB compressed BCF at 1kGP density is on the order of
265k records spanning a few Mb, not the whole chromosome). Declaring the real
hg38 length over that sparse prefix would mean most region queries a
benchmark issues — anything outside the populated prefix — return nothing,
which is exactly the pathological case a benchmark must not trip over.
Declaring `length` as the generated span instead means density is realistic
*everywhere in the declared contig*, so an arbitrary region query always
lands on realistic data.

## Why genotypes are drawn i.i.d.

Genotypes are drawn independently per site from Hardy-Weinberg equilibrium at
each site's sampled allele count — there is no linkage disequilibrium (LD) or
haplotype structure between sites. Within a site, the sampled allele count
`ac` is placed exactly: `ac` alt alleles are assigned uniformly at random
among that site's non-missing genotype slots (sampling without replacement),
rather than drawing each slot independently at the implied frequency
`ac / n_alleles` — the independent draw would re-randomise the realised
allele count away from `ac` and destroy the fitted site-frequency spectrum on
any single record.

This follows directly from a controlled ablation: permuting sample labels
independently at every variant (destroying LD while holding site-level
statistics — allele counts, SFS, variant count, raw byte size — exactly
fixed) changed BCF size by only **1.14x**, `bcftools view -H` decode time by
**1.00x**, and peak RSS by **1.000x**. LD only pays off where a compressor's
window can span multiple rows (`xz -6` saw 2.05x); BCF's 32 KB bgzip window
essentially cannot. Since the crate's benchmarking targets (genoray, GVL,
bcftools) are parse-bound readers over BCF, simulating LD — coalescent
simulation, Li-Stephens haplotype copying, block resampling — would spend
real budget on a ~0x lever on the metrics this tool exists to support.

## Determinism

Same seed, profile, and spec always produce byte-identical output, regardless
of worker/thread count. Each block of records seeds its own PRNG from a pure
function of `(seed, block index)`, so results assemble the same way no matter
which thread computed which block or in what order.

## Sizing

Three ways to say how much to generate:

- `Size::RecordsPerContig(n)` / `--records-per-contig N` — exactly `n`
  records for each requested contig.
- `Size::Records(n)` / `--records N` — exactly `n` records total, split
  across contigs proportional to each contig's fitted density.
- `Size::Target(bytes)` / `--target-size 100MB` — generate until the
  compressed output reaches *at least* the target, then stop. This
  overshoots rather than undershoots, and the overshoot is a small
  percentage of the target rather than a fixed byte budget (observed on the
  order of ~9% in practice) — don't rely on the output landing within any
  absolute byte window of the target.

## API example

```rust
{{#rustdoc_include ../../../examples/bulk.rs:bulk}}
```

## CLI example

```bash
vcfixture bulk --profile germline-1kgp --samples 3202 \
  --contigs chr1,chr2,chr3 --target-size 100MB --seed 42 -o bench.bcf
```

`--profile` accepts either a builtin profile name or a path to a profile JSON
file, tried in that order. Three builtin profiles are shipped:

- `germline-1kgp` — fitted from the 1kGP phased panel.
- `germline-1kgp-unphased` — fitted from the same cohort's raw unphased
  GATK sites-only callset.
- `somatic-gdc` — fitted from a merged 16,007-sample GDC cohort.

Two germline profiles exist because phased and unphased data can't be one
file: phasing drops unphaseable singletons, so the phased panel has no
singletons, while the unphased callset has no phase. Pass a path (e.g. one
produced by `pixi run fit`, or any JSON matching the `Profile` schema) to use
anything else.

Output is `<output>` plus a `.csi` index and a `<output>.summary.json` — a
summary truth (per-contig record counts and position ranges, allele counts,
variant-class counts, and a genotype checksum) computed for free while
streaming, in place of a per-genotype oracle that would be larger than the
file it describes at this scale. See `Summary` in the API reference for its
full shape.
