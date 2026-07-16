# `scripts/fit/`

`fit_profile.py` turns a real plink2 pgen/pvar cohort into a `vcfixture` bulk-generation
profile JSON (`profiles/*.json`), the file `Profile::builtin` embeds via
`include_str!` in `src/bulk/profile.rs`.

## What it does

The profile schema splits every profile into two parts, and this script only ever
writes the first one:

- **`fitted`** — statistics measured from the source cohort: per-contig
  variant density, the inter-variant gap distribution, the site-frequency
  spectrum (allele-count histogram), variant-class mix, indel-length
  distribution, Ti/Tv, multiallelic rate, missing rate, phased rate, and
  ploidy. Every value here is computed from `--pgen`; the script never
  hand-picks a `fitted` value.
- **`dialed`** — generation choices independent of any fit. Today that's just
  `payload` (the FORMAT-field preset: `gt-only`, `gt-vaf`, `gatk`, or
  `mutect2`), which cannot be inferred from a pgen (pgen stores genotypes
  only, no per-record FORMAT payload) and is instead passed via `--payload`.

`provenance` (source path, sample/variant counts, fit date, tool version) is
always populated from the real run -- never left as a placeholder.

Internally: `polars.scan_csv` (lazy, streaming) reads the pvar for contigs,
gaps, variant classification, indel lengths, Ti/Tv, and multiallelic rate;
`plink2 --freq counts` supplies the site-frequency spectrum; `plink2
--missing` supplies the missing-call rate; a bounded VCF export
(`--phase-sample-mb`, default 1 Mb from the first fitted contig) supplies an
estimate of the phased-call rate, since pgen has no direct phase-fraction
report.

The **singleton SFS bin is exactly `[1, 2)`** by construction (`_sfs_edges`
doubles from 1 upward): real 1kGP high-coverage data is ~47.6% singletons
vs. ~12.3% for a neutral coalescent SFS, and getting that first bin wrong
would silently produce a neutral-looking spectrum instead of the empirical
one the whole point of this script is to capture.

This script is never imported by Rust. Its only contract with the rest of
the crate is the JSON schema in `src/bulk/profile.rs`.

## Usage

```bash
pixi run -e fit fit -- \
    --pgen <plink2 prefix> \
    --name <profile-name> \
    --out profiles/<profile-name>.json \
    [--payload gt-only|gt-vaf|gatk|mutect2] \
    [--contigs chr1 chr2 ...] \
    [--ploidy 2]
```

`--pgen` is a plink2 fileset *prefix*: the script expects
`<prefix>.pgen`, `<prefix>.psam`, and either `<prefix>.pvar.zst` or
`<prefix>.pvar`.

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
inputs it constructs inline -- it never touches `/carter`, requires no
network access, and does not require `plink2` (the plink2-shelling
functions, `fit_sfs`/`fit_missing_rate`/`fit_phased_rate`, are exercised only
by `main()`, which is not covered by CI).
