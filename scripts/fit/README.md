# `scripts/fit/`

`fit_profile.py` turns a real cohort into a `vcfixture` bulk-generation
profile JSON (`profiles/*.json`), the file `Profile::builtin` embeds via
`include_str!` in `src/bulk/profile.rs`. It supports two mutually exclusive
input modes:

- **`--pgen`** — a plink2 fileset with genotypes. `n_samples`, `sfs`,
  `missing_rate`, and `phased_rate` are all fitted from it.
- **`--sites-vcf`** — a sites-only VCF/BCF with no genotype columns at all
  (only `INFO/AC`/`INFO/AN`), e.g. the 1000 Genomes raw GATK callset -- the
  only available source of a realistic 1kGP allele-frequency spectrum (the
  phased panel already fit has 0% singletons, since phasing drops
  unphaseable singletons). `n_samples` and `phased_rate` cannot be derived
  from a sites-only file and must be passed explicitly via `--n-samples`
  and `--phased-rate`.

## What it does

The profile schema splits every profile into two parts, and this script only ever
writes the first one:

- **`fitted`** — statistics measured from the source cohort: per-contig
  variant density, the inter-variant gap distribution, the site-frequency
  spectrum (allele-count histogram), variant-class mix, indel-length
  distribution, Ti/Tv, multiallelic rate, missing rate, phased rate, and
  ploidy. Every value here is computed from the source (`--pgen` or
  `--sites-vcf`) -- with the one exception of `phased_rate` under
  `--sites-vcf`, which is passed verbatim via `--phased-rate` since a
  sites-only source has no genotypes to fit it from. The script never
  hand-picks any other `fitted` value.
- **`dialed`** — generation choices independent of any fit. Today that's just
  `payload` (the FORMAT-field preset: `gt-only`, `gt-vaf`, `gatk`, or
  `mutect2`), which cannot be inferred from either source (neither carries
  a per-record FORMAT payload) and is instead passed via `--payload`.

`provenance` (source path, sample/variant counts, fit date, tool version) is
always populated from the real run -- never left as a placeholder.

Internally, for `--pgen`: `polars.scan_csv` (lazy, streaming) reads the pvar
for contigs, gaps, variant classification, indel lengths, Ti/Tv, and
multiallelic rate; `plink2 --freq counts` supplies the site-frequency
spectrum; `plink2 --missing` supplies the missing-call rate; a bounded VCF
export (`--phase-sample-mb`, default 1 Mb from the first fitted contig)
supplies an estimate of the phased-call rate, since pgen has no direct
phase-fraction report.

For `--sites-vcf`: `bcftools query -i 'FILTER="PASS"'` reads
`CHROM/POS/REF/ALT/AC/AN` (dropping VQSR-tranche/`LowQual` records, which
make up ~14-16% of the raw 1kGP callset) into the same lazy pipeline that
computes contigs/gaps/classification/indel-lengths/Ti-Tv/multiallelic-rate
for `--pgen` (`compute_pvar_stats` is shared, unduplicated, between both
paths). The site-frequency spectrum is built directly from `INFO/AC`
(multiallelics are already split into one row per ALT allele upstream, so
no explode step is needed); the missing rate is derived from `INFO/AN`
relative to `2 * --n-samples`.

**Multiallelic records.** plink2 pvar retains multiallelic sites natively
(it does NOT auto-split them the way `bcftools norm -m-` would), so `ALT`
can be a comma-joined list like `"G,T"`. Two different units of observation
apply and must not be conflated:

- per-**record**: `contigs` (density), `gap_dist`, and `multiallelic_rate`
  count one pvar row as one observation, regardless of how many ALT alleles
  it lists.
- per-**allele**: `variant_classes`, `indel_length`, and `titv` each count
  one observation per REF/ALT pair -- a site with 2 ALT alleles
  (`_explode_alleles`) contributes 2 independent observations to each of
  these three, computed by `compute_pvar_stats`. `fit_sfs` mirrors this on
  the plink2 side: `plink2 --freq counts` keeps `ALT_CTS` comma-joined and
  positionally aligned with `ALT` (e.g. `ALT="G,T"`, `ALT_CTS="1,1"`), and
  `_split_alt_cts` splits that into one allele-count observation per
  `sfs` input.

A biallelic-only pvar is unaffected either way -- splitting a single-allele
ALT is a no-op.

**Out-of-range histogram values.** `histogram()` drops values outside the
edge range (as `numpy.histogram` already does) but warns whenever *any*
values are dropped, not just when *all* of them are -- e.g. `gap_dist`
values beyond the 1e5 bp cap (real inter-variant gaps over centromeres/
telomeres) would otherwise vanish from the fitted distribution with no
diagnostic output.

The **singleton SFS bin is exactly `[1, 2)`** by construction (`_sfs_edges`
doubles from 1 upward): real 1kGP high-coverage data is ~47.6% singletons
vs. ~12.3% for a neutral coalescent SFS, and getting that first bin wrong
would silently produce a neutral-looking spectrum instead of the empirical
one the whole point of this script is to capture.

This script is never imported by Rust. Its only contract with the rest of
the crate is the JSON schema in `src/bulk/profile.rs`.

## Usage

```bash
# --pgen: genotyped plink2 fileset
pixi run -e fit fit -- \
    --pgen <plink2 prefix> \
    --name <profile-name> \
    --out profiles/<profile-name>.json \
    [--payload gt-only|gt-vaf|gatk|mutect2] \
    [--contigs chr1 chr2 ...] \
    [--ploidy 2]

# --sites-vcf: sites-only VCF/BCF, no genotypes
pixi run -e fit fit -- \
    --sites-vcf <path/to/sites.vcf.gz> \
    --name <profile-name> \
    --out profiles/<profile-name>.json \
    --n-samples <cohort size> \
    --phased-rate <fixed rate in [0, 1]> \
    [--payload gt-only|gt-vaf|gatk|mutect2] \
    [--contigs chr1 chr2 ...] \
    [--ploidy 2]
```

`--pgen` and `--sites-vcf` are mutually exclusive; exactly one is required.

`--pgen` is a plink2 fileset *prefix*: the script expects
`<prefix>.pgen`, `<prefix>.psam`, and either `<prefix>.pvar.zst` or
`<prefix>.pvar`.

`--sites-vcf` requires `--n-samples` (a sites-only VCF has no sample columns
to derive the cohort size from) and `--phased-rate` (no genotypes to count
phased/unphased calls from); both are rejected if passed with `--pgen`,
which fits both automatically instead.

## Re-fitting the two committed profiles

```bash
# germline-1kgp (3,202 samples; gt-only is 1kGP-faithful)
pixi run -e fit fit -- \
    --pgen /carter/users/dlaub/data/1kGP/plink2/hg38.norm \
    --name germline-1kgp \
    --out profiles/germline-1kgp.json \
    --payload gt-only

# somatic-gdc (16,007 samples; note the source .pvar is uncompressed,
# not .pvar.zst -- the prefix is the same either way)
pixi run -e fit fit -- \
    --pgen /carter/shared/data/gdc/somatic/wgs_DR45/results/gdc_wgs_DR45.gt \
    --name somatic-gdc \
    --out profiles/somatic-gdc.json \
    --payload gt-vaf
```

Re-fitting new data is the same one command with a different `--pgen` prefix.
The output must be committed to `profiles/` -- `Profile::builtin` embeds it
at compile time via `include_str!`, so the crate only ever reads the
committed JSON, never `/carter` at build or run time.

## Tests

```bash
pixi run -e fit test-fit
```

The test suite (`test_fit_profile.py`) exercises only the small synthetic
inputs it constructs inline -- it never touches `/carter` and requires no
network access. Most tests need only `polars`/`numpy` (no `plink2`/
`bcftools`); a handful of end-to-end tests build a tiny synthetic
multiallelic VCF, convert it with `plink2 --make-pgen`, and run it through
the full extraction pipeline (the gold-standard regression coverage for the
multiallelic crash this script used to hit) -- those are marked
`@pytest.mark.skipif(not PLINK2_AVAILABLE, ...)` and skip cleanly wherever
`plink2` is not on `PATH`, such as CI. Similarly, `read_sites_vcf` and
`--sites-vcf` end-to-end tests build a small bgzipped+tabix-indexed sites
VCF with `bcftools view -Oz` + `bcftools index -t` and are marked
`@pytest.mark.skipif(not BCFTOOLS_AVAILABLE, ...)`.
