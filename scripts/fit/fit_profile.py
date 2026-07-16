#!/usr/bin/env python
"""Fit a vcfixture bulk-generation ``Profile`` JSON from real plink2 pgen/pvar data.

This script is the only place that turns a real cohort (1kGP high-coverage
pgen, GDC somatic pgen, ...) into the committed profile JSON consumed by the
Rust bulk generator (see ``src/bulk/profile.rs``). It is never imported by
Rust and has no contract other than that JSON schema.

The profile schema deliberately separates ``fitted`` (statistics measured
from ``--pgen``) from ``dialed`` (generation choices, e.g. FORMAT payload,
that are picked by the user and never claimed to be measured). Every value
this script writes under ``fitted`` is derived from the source data passed
in; nothing here should ever be a hand-picked literal.

Usage
-----
    pixi run -e fit fit -- \\
        --pgen /path/to/prefix --name germline-1kgp --out profiles/germline-1kgp.json

See ``scripts/fit/README.md`` for the exact commands used to fit the two
committed profiles.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import polars as pl

__version__ = "0.1.0"

# Must match the field names of `ClassMix` in src/bulk/profile.rs exactly.
CLASS_NAMES = ("snp", "insertion", "deletion", "mnp", "complex", "symbolic")

# ~90% of indels are <= 6 bp (see the design spec), so resolve that range
# finely and taper off for the long tail.
INDEL_EDGES = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 20.0, 50.0, 100.0, 1000.0]

# Log-spaced bins over [1, 1e5] for inter-variant gaps (bp).
GAP_LOW, GAP_HIGH, GAP_N_BINS = 1.0, 1e5, 24

TRANSITION_PAIRS = frozenset({"AG", "GA", "CT", "TC"})


# --------------------------------------------------------------------------
# Schema-facing helpers (unit tested directly, see test_fit_profile.py)
# --------------------------------------------------------------------------


def histogram(
    values: Sequence[float] | np.ndarray | pl.Series, edges: Sequence[float]
) -> dict[str, list[float]]:
    """Bin `values` into `edges` and normalize to weights summing to 1.0.

    Returns ``{"edges": [...], "weights": [...]}`` with
    ``len(weights) == len(edges) - 1``, matching the Rust ``Histogram``
    schema in ``src/bulk/profile.rs``.

    Values outside ``[edges[0], edges[-1]]`` are dropped from the
    normalization (as `numpy.histogram` already does implicitly) but, unlike
    a bare `numpy.histogram` call, this warns whenever *any* values are
    dropped -- not just when *all* of them are -- since a silently-dropped
    tail (e.g. inter-variant gaps over centromeres/telomeres exceeding the
    `gap_dist` cap) would otherwise skew the fitted distribution with zero
    diagnostic output.
    """
    edges = [float(e) for e in edges]
    if len(edges) < 2:
        raise ValueError("histogram needs >= 2 edges")
    if any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError("histogram edges must be strictly increasing")
    arr = np.asarray(values, dtype=np.float64)
    counts, _ = np.histogram(arr, bins=edges)
    n_in_range = int(counts.sum())
    n_total = int(arr.size)
    n_dropped = n_total - n_in_range
    if n_dropped > 0:
        frac = n_dropped / n_total
        warnings.warn(
            f"histogram: dropped {n_dropped}/{n_total} values ({frac:.1%}) "
            f"outside the edge range [{edges[0]}, {edges[-1]}]",
            stacklevel=2,
        )
    if n_in_range <= 0:
        raise ValueError("no values fell within the histogram edges")
    weights = (counts / n_in_range).tolist()
    return {"edges": edges, "weights": weights}


def class_mix_from_counts(counts: Mapping[str, int]) -> dict[str, float]:
    """Normalize variant-class counts to frequencies summing to 1.0.

    `counts` must have exactly the six `CLASS_NAMES` keys, matching the
    Rust `ClassMix` struct fields.
    """
    missing = set(CLASS_NAMES) - set(counts)
    if missing:
        raise ValueError(f"missing variant classes: {sorted(missing)}")
    total = sum(counts[name] for name in CLASS_NAMES)
    if total <= 0:
        raise ValueError("variant class counts must sum > 0")
    return {name: counts[name] / total for name in CLASS_NAMES}


def classify(ref: str, alt: str) -> str:
    """Classify a single REF/ALT pair into one of `CLASS_NAMES`.

    This is the scalar reference implementation. It assumes `alt` is a
    *single* allele -- a comma-joined multiallelic ALT (e.g. "G,T") must be
    split into one call per allele by the caller; passing the raw joined
    string here silently misclassifies (e.g. `classify("A", "A,G")` reads
    as an insertion, not two independent alleles). The extraction pipeline
    handles this via `_explode_alleles` before `_classify_expr` ever runs.

    The extraction pipeline uses `_classify_expr`, a vectorized polars
    expression that computes the same classification over an entire pvar
    column at once -- looping this function over 10s-100s of millions of
    pvar rows in Python would be far too slow.
    """
    if alt.startswith("<") or "[" in alt or "]" in alt:
        return "symbolic"
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(ref) == len(alt) and len(ref) > 1:
        return "mnp"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "insertion"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "deletion"
    return "complex"


def _classify_expr(ref: pl.Expr, alt: pl.Expr) -> pl.Expr:
    """Vectorized equivalent of `classify`, for use over a full pvar column.

    Like `classify`, this assumes `alt` is a single allele per row -- run it
    only after `_explode_alleles` has split any comma-joined multiallelic
    ALT into one row per allele.
    """
    is_symbolic = alt.str.starts_with("<") | alt.str.contains(r"[\[\]]", literal=False)
    ref_len = ref.str.len_chars()
    alt_len = alt.str.len_chars()
    is_snp = (ref_len == 1) & (alt_len == 1)
    is_mnp = (ref_len == alt_len) & (ref_len > 1)
    is_insertion = (alt_len > ref_len) & alt.str.starts_with(ref)
    is_deletion = (ref_len > alt_len) & ref.str.starts_with(alt)
    return (
        pl.when(is_symbolic)
        .then(pl.lit("symbolic"))
        .when(is_snp)
        .then(pl.lit("snp"))
        .when(is_mnp)
        .then(pl.lit("mnp"))
        .when(is_insertion)
        .then(pl.lit("insertion"))
        .when(is_deletion)
        .then(pl.lit("deletion"))
        .otherwise(pl.lit("complex"))
    )


def read_pvar(path: str | Path) -> pl.DataFrame:
    """Read a `.pvar` or `.pvar.zst` file, keeping only #CHROM/POS/ID/REF/ALT.

    A 1kGP-scale pvar can be 500+ MB of text. This scans lazily
    (`polars.scan_csv`) and collects with the streaming engine so the file is
    never fully materialized as an eager, unstreamed read; `.zst` decoding is
    handled natively by polars based on the file extension.

    ``ALT`` may be a comma-joined multiallelic list (e.g. "G,T") exactly as
    plink2 pvar stores it natively -- pgen does NOT auto-split multiallelic
    records. Callers that need a single REF/ALT pair per row (classification,
    indel length, Ti/Tv) must run `_explode_alleles` first; callers that
    operate per-*record* (contig density, inter-variant gaps,
    `multiallelic_rate`) should use the frame as returned here.
    """
    # CHROM must stay a string even for numeric-looking contig names like
    # "1" or "22" -- ContigStat.id is a String in src/bulk/profile.rs, and
    # letting polars infer it as an integer breaks every downstream
    # plink2 --chr argument (and str/int comparisons) built from it.
    lf = pl.scan_csv(
        str(path),
        separator="\t",
        comment_prefix="##",
        schema_overrides={"#CHROM": pl.Utf8, "ID": pl.Utf8, "REF": pl.Utf8, "ALT": pl.Utf8},
    ).rename({"#CHROM": "CHROM"})
    return lf.select(["CHROM", "POS", "ID", "REF", "ALT"]).collect(engine="streaming")


def build_profile(
    *,
    name: str,
    source: str,
    n_samples: int,
    contigs: list[dict],
    gaps: Sequence[float],
    acs: Sequence[float],
    indel_lens: Sequence[float],
    class_counts: Mapping[str, int],
    titv: float,
    multiallelic_rate: float,
    missing_rate: float,
    phased_rate: float,
    ploidy: int,
    payload: str = "gt-only",
    n_variants_source: int | None = None,
) -> dict:
    """Assemble a schema-valid Profile dict (see `src/bulk/profile.rs`).

    Every field under "fitted" is derived from `gaps`/`acs`/`indel_lens`/
    `class_counts`/`contigs`/the scalar rates passed in by the caller (which
    `main()` computes from real pgen/pvar data) -- never hand-picked here.
    `dialed.payload` is the one deliberate exception: it is a generation
    choice, not a fitted statistic.
    """
    if n_variants_source is None:
        n_variants_source = sum(c["n_variants"] for c in contigs)

    fitted = {
        "contigs": contigs,
        "gap_dist": histogram(gaps, _gap_edges()),
        "sfs": histogram(acs, _sfs_edges(n_samples)),
        "variant_classes": class_mix_from_counts(class_counts),
        "indel_length": histogram(indel_lens, INDEL_EDGES),
        "titv": titv,
        "multiallelic_rate": multiallelic_rate,
        "missing_rate": missing_rate,
        "phased_rate": phased_rate,
        "ploidy": ploidy,
    }
    return {
        "name": name,
        "provenance": {
            "source": source,
            "n_samples_source": n_samples,
            "n_variants_source": n_variants_source,
            "fitted_on": _dt.date.today().isoformat(),
            "fit_tool_version": __version__,
        },
        "fitted": fitted,
        "dialed": {"payload": payload},
    }


# --------------------------------------------------------------------------
# Edge construction
# --------------------------------------------------------------------------


def _log_spaced_edges(low: float, high: float, n_bins: int) -> list[float]:
    if high <= low:
        raise ValueError(f"high ({high}) must be > low ({low})")
    return np.geomspace(low, high, n_bins + 1).tolist()


def _gap_edges() -> list[float]:
    """Log-spaced edges over [1, 1e5] for the inter-variant gap distribution."""
    return _log_spaced_edges(GAP_LOW, GAP_HIGH, GAP_N_BINS)


def _sfs_edges(n_samples: int) -> list[float]:
    """Log-spaced (base-2) edges over [1, 2*n_samples], first bin exactly [1, 2).

    The singleton bin [1, 2) is the statistically load-bearing one: real
    1kGP high-coverage data is ~47.6% singletons vs. ~12.3% for a neutral
    constant-Ne coalescent SFS (see
    docs/superpowers/specs/2026-07-16-bulk-generation-design.md). A base-2
    doubling sequence guarantees the first two edges are exactly 1.0 and 2.0
    for any n_samples, unlike a generic `geomspace` whose step size depends
    on the chosen bin count.
    """
    high = 2 * n_samples
    if high < 2:
        raise ValueError("need at least 1 sample to build an SFS")
    edges = [1.0]
    e = 2.0
    while e < high:
        edges.append(e)
        e *= 2.0
    edges.append(float(high))
    return edges


# --------------------------------------------------------------------------
# Extraction from a pvar DataFrame (vectorized polars, no Python row loops)
# --------------------------------------------------------------------------


def _explode_alleles(df: pl.DataFrame) -> pl.DataFrame:
    """Split multiallelic ALT into one row per REF/ALT allele pair.

    plink2 pvar retains multiallelic records natively -- it does NOT
    auto-split them the way `bcftools norm -m-` would -- so `ALT` can be a
    comma-joined list (e.g. "G,T"). `classify()` / `_classify_expr` are only
    defined for a single REF/ALT pair, so this must run before any
    per-allele statistic (variant class, indel length, Ti/Tv) is computed.

    This is a per-*allele* view: a site with 2 ALT alleles becomes 2 rows,
    each with the shared CHROM/POS/ID/REF and its own single ALT. A
    biallelic record (the common case) splits into a 1-element list and
    explodes back to the same single row, so this is a no-op for it.
    Callers that need per-*record* semantics instead (contig density,
    inter-variant gaps, `multiallelic_rate`) must use the un-exploded frame.

    Vectorized via `Series.str.split` + `DataFrame.explode` -- no Python
    loop over rows, so this streams fine over a 500+ MB pvar.
    """
    return df.with_columns(pl.col("ALT").str.split(",")).explode("ALT", empty_as_null=False)


def _classify_df(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(_classify_expr(pl.col("REF"), pl.col("ALT")).alias("class"))


def _contig_stats(df: pl.DataFrame) -> list[dict]:
    """Per-contig record count and density. `df` is per-RECORD (not exploded)."""
    agg = (
        df.group_by("CHROM", maintain_order=True)
        .agg(
            n_variants=pl.len(),
            pos_min=pl.col("POS").min(),
            pos_max=pl.col("POS").max(),
        )
        .with_columns(span_bp=(pl.col("pos_max") - pl.col("pos_min")).clip(lower_bound=1))
        .with_columns(density_per_kb=pl.col("n_variants") / (pl.col("span_bp") / 1000.0))
        .sort("CHROM")
    )
    return [
        {
            "id": row["CHROM"],
            "n_variants": row["n_variants"],
            "density_per_kb": row["density_per_kb"],
        }
        for row in agg.iter_rows(named=True)
    ]


def _gaps(df: pl.DataFrame) -> pl.Series:
    """Inter-variant gaps (bp) within each contig, sorted-position diffs.

    `df` is per-RECORD (not exploded): a multiallelic site is still one
    position, so exploding first would fabricate spurious zero-length gaps
    between the alleles of the same site.
    """
    return (
        df.sort(["CHROM", "POS"])
        .select(pl.col("POS").diff().over("CHROM").alias("gap"))
        .drop_nulls()
        .filter(pl.col("gap") > 0)["gap"]
    )


def _class_counts(alleles: pl.DataFrame) -> dict[str, int]:
    """Tally classes over a per-ALLELE frame (post `_explode_alleles` + `_classify_df`)."""
    counts = dict.fromkeys(CLASS_NAMES, 0)
    tally = alleles.group_by("class").agg(pl.len().alias("n"))
    for row in tally.iter_rows(named=True):
        counts[row["class"]] = row["n"]
    return counts


def _indel_lengths(alleles: pl.DataFrame) -> pl.Series:
    """Indel lengths over a per-ALLELE frame (post `_explode_alleles` + `_classify_df`)."""
    indels = alleles.filter(pl.col("class").is_in(["insertion", "deletion"]))
    # str.len_chars() is UInt32; subtracting two UInt32 columns wraps around
    # on underflow instead of going negative, so cast to a signed type first.
    alt_len = pl.col("ALT").str.len_chars().cast(pl.Int64)
    ref_len = pl.col("REF").str.len_chars().cast(pl.Int64)
    return indels.select((alt_len - ref_len).abs().alias("len"))["len"]


def _multiallelic_rate(df: pl.DataFrame) -> float:
    """Fraction of pvar RECORDS whose ALT field lists more than one allele.

    Must run on the un-exploded, per-record frame: this counts sites, not
    alleles, so a triallelic site still contributes exactly 1 to both the
    numerator and denominator, not 2.
    """
    n = df.height
    if n == 0:
        return 0.0
    n_multi = df.select(pl.col("ALT").str.contains(",", literal=True).sum()).item()
    return n_multi / n


def _titv(alleles: pl.DataFrame) -> float:
    """Transition/transversion ratio over SNP alleles (post `_explode_alleles` + `_classify_df`)."""
    snps = alleles.filter(pl.col("class") == "snp")
    if snps.height == 0:
        raise ValueError("no SNPs found; cannot compute Ti/Tv")
    pair = pl.concat_str([pl.col("REF"), pl.col("ALT")])
    n_ts = snps.select(pair.is_in(TRANSITION_PAIRS).sum()).item()
    n_tv = snps.height - n_ts
    if n_tv == 0:
        raise ValueError("no transversions found; cannot compute a finite Ti/Tv ratio")
    return n_ts / n_tv


def compute_pvar_stats(df: pl.DataFrame) -> dict:
    """Compute every `fitted` statistic derived purely from the pvar frame.

    `df` is the raw, per-record frame from `read_pvar` (optionally
    contig-filtered) -- ALT may still be comma-joined for multiallelic
    sites. Two different units of observation are in play here and must
    not be conflated:

    - per-RECORD: `contigs` (density), `gap_dist` inputs, and
      `multiallelic_rate` all count one pvar row as one observation,
      regardless of how many ALT alleles it lists.
    - per-ALLELE: `variant_classes`, `indel_length`, and `titv` each count
      one observation per REF/ALT pair -- a site with 2 ALT alleles
      contributes 2 independent class/indel-length/Ti-Tv observations.
      `classify()` is only defined for a single allele, so `ALT` is split
      via `_explode_alleles` before any of these three run.

    A biallelic-only pvar is unaffected either way: `_explode_alleles` is a
    no-op per row when ALT has no comma.
    """
    contigs = _contig_stats(df)
    gaps = _gaps(df)
    multiallelic_rate = _multiallelic_rate(df)

    alleles = _classify_df(_explode_alleles(df))
    class_counts = _class_counts(alleles)
    indel_lens = _indel_lengths(alleles)
    titv = _titv(alleles)

    return {
        "contigs": contigs,
        "gaps": gaps,
        "multiallelic_rate": multiallelic_rate,
        "class_counts": class_counts,
        "indel_lens": indel_lens,
        "titv": titv,
    }


# --------------------------------------------------------------------------
# plink2 subprocess helpers
# --------------------------------------------------------------------------


def _run_plink2(args: list[str]) -> None:
    result = subprocess.run(["plink2", *args], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"plink2 {' '.join(args)} failed (exit {result.returncode}):\n"
            f"{result.stdout}\n{result.stderr}"
        )


def _split_alt_cts(alt_cts: Sequence[str] | pl.Series) -> list[float]:
    """Parse plink2 `--freq counts` ALT_CTS into one count per ALT allele.

    For a multiallelic site, plink2's `.acount` keeps `ALT` and `ALT_CTS`
    comma-joined and positionally aligned (e.g. `ALT="G,T"`,
    `ALT_CTS="1,1"`): the i-th count belongs to the i-th ALT allele. This
    splits that alignment into one float observation per allele, matching
    the per-allele split `_explode_alleles` does for classification -- so
    each ALT allele contributes its own count to the `sfs` histogram
    instead of the whole comma-joined string being handed to
    `np.asarray(..., dtype=float)`, which raises on a string like "1,1".

    A biallelic row (`ALT_CTS` with no comma) splits into a 1-element list
    and is unaffected. Vectorized via polars `Series.str.split` + explode --
    no Python loop, though `.acount` files are one row per site so this
    would be cheap either way.
    """
    s = pl.Series("ALT_CTS", list(alt_cts), dtype=pl.Utf8)
    return s.str.split(",").explode(empty_as_null=False).cast(pl.Float64).to_list()


def fit_sfs(pgen_prefix: str | Path, contigs: Iterable[str] | None = None) -> list[float]:
    """Shell out to `plink2 --freq counts` and return one ALT_CTS per allele.

    ALT_CTS is the observed non-reference allele count per site -- exactly
    the site-frequency-spectrum input the `sfs` histogram is fit from. For
    multiallelic sites plink2 keeps ALT_CTS comma-joined (aligned with the
    equally comma-joined ALT column); `_split_alt_cts` splits that into one
    count per allele.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "sfs")
        args = ["--pfile", str(pgen_prefix), "--freq", "counts", "--out", out_prefix]
        if contigs:
            args += ["--chr", ",".join(contigs)]
        _run_plink2(args)
        df = pl.read_csv(
            f"{out_prefix}.acount", separator="\t", schema_overrides={"ALT_CTS": pl.Utf8}
        )
        return _split_alt_cts(df["ALT_CTS"])


def fit_missing_rate(pgen_prefix: str | Path, contigs: Iterable[str] | None = None) -> float:
    """Shell out to `plink2 --missing` and return the global hardcall missing rate."""
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "miss")
        args = [
            "--pfile",
            str(pgen_prefix),
            "--missing",
            "variant-only",
            "--out",
            out_prefix,
        ]
        if contigs:
            args += ["--chr", ",".join(contigs)]
        _run_plink2(args)
        df = pl.read_csv(f"{out_prefix}.vmiss", separator="\t")
        n_missing = df["MISSING_CT"].sum()
        n_obs = df["OBS_CT"].sum()
        return n_missing / n_obs if n_obs else 0.0


def fit_phased_rate(
    pgen_prefix: str | Path,
    contig: str,
    pos_min: int,
    window_bp: int = 1_000_000,
) -> float:
    """Estimate the fraction of genotype calls that are phased.

    pgen has no direct "phased fraction" report, so this exports a bounded
    window (`window_bp` starting at `pos_min` on `contig`) to VCF and counts
    the phased ('|') vs unphased ('/') GT separators plink2 actually wrote,
    across every sample x variant call in that window.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "phase")
        _run_plink2(
            [
                "--pfile",
                str(pgen_prefix),
                "--chr",
                contig,
                "--from-bp",
                str(pos_min),
                "--to-bp",
                str(pos_min + window_bp),
                "--export",
                "vcf",
                "--out",
                out_prefix,
            ]
        )
        vcf_path = Path(f"{out_prefix}.vcf")
        n_phased = 0
        n_total = 0
        with open(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                for gt_field in fields[9:]:
                    gt = gt_field.split(":", 1)[0]
                    if "|" in gt:
                        n_phased += 1
                    n_total += 1
        return n_phased / n_total if n_total else 1.0


def _read_n_samples(psam_path: Path) -> int:
    n = 0
    with open(psam_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            n += 1
    return n


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Fit a vcfixture bulk-generation profile JSON from plink2 pgen/pvar data."
    )
    parser.add_argument(
        "--pgen",
        required=True,
        help="plink2 fileset prefix (expects <prefix>.pgen/.psam/.pvar[.zst])",
    )
    parser.add_argument("--name", required=True, help="profile name, e.g. germline-1kgp")
    parser.add_argument("--out", required=True, type=Path, help="output profile JSON path")
    parser.add_argument(
        "--payload",
        default="gt-only",
        choices=["gt-only", "gt-vaf", "gatk", "mutect2"],
        help="dialed FORMAT preset (not fitted from data; see src/bulk/profile.rs::Payload)",
    )
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=None,
        help="restrict fitting to these contig IDs (default: all present in the source)",
    )
    parser.add_argument(
        "--ploidy", type=int, default=2, help="sample ploidy recorded in the profile"
    )
    parser.add_argument(
        "--phase-sample-mb",
        type=float,
        default=1.0,
        help="window size (Mb), from the first fitted contig's start, used to estimate phased_rate",
    )
    args = parser.parse_args(argv)

    prefix = str(args.pgen)
    psam_path = Path(f"{prefix}.psam")
    pvar_path = Path(f"{prefix}.pvar.zst")
    if not pvar_path.exists():
        pvar_path = Path(f"{prefix}.pvar")
    if not pvar_path.exists():
        raise FileNotFoundError(f"no .pvar or .pvar.zst found for prefix {prefix}")

    n_samples = _read_n_samples(psam_path)

    df = read_pvar(pvar_path)
    if args.contigs:
        df = df.filter(pl.col("CHROM").is_in(args.contigs))

    stats = compute_pvar_stats(df)
    contigs = stats["contigs"]
    if not contigs:
        raise ValueError("no variants found for the requested contig(s)")
    contig_ids = [c["id"] for c in contigs]

    gaps = stats["gaps"]
    class_counts = stats["class_counts"]
    indel_lens = stats["indel_lens"]
    titv = stats["titv"]
    multiallelic_rate = stats["multiallelic_rate"]

    chr_filter = contig_ids if args.contigs else None
    acs = fit_sfs(args.pgen, contigs=chr_filter)
    missing_rate = fit_missing_rate(args.pgen, contigs=chr_filter)

    first_contig = contigs[0]
    pos_min = int(df.filter(pl.col("CHROM") == first_contig["id"])["POS"].min())
    phased_rate = fit_phased_rate(
        args.pgen,
        contig=first_contig["id"],
        pos_min=pos_min,
        window_bp=int(args.phase_sample_mb * 1_000_000),
    )

    profile = build_profile(
        name=args.name,
        source=prefix,
        n_samples=n_samples,
        contigs=contigs,
        gaps=gaps,
        acs=acs,
        indel_lens=indel_lens,
        class_counts=class_counts,
        titv=titv,
        multiallelic_rate=multiallelic_rate,
        missing_rate=missing_rate,
        phased_rate=phased_rate,
        ploidy=args.ploidy,
        payload=args.payload,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(profile, indent=2) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
