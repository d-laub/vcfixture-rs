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
    """
    edges = [float(e) for e in edges]
    if len(edges) < 2:
        raise ValueError("histogram needs >= 2 edges")
    if any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError("histogram edges must be strictly increasing")
    counts, _ = np.histogram(np.asarray(values, dtype=np.float64), bins=edges)
    total = counts.sum()
    if total <= 0:
        raise ValueError("no values fell within the histogram edges")
    weights = (counts / total).tolist()
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

    This is the scalar reference implementation. The extraction pipeline
    uses `_classify_expr`, a vectorized polars expression that computes the
    same classification over an entire pvar column at once -- looping this
    function over 10s-100s of millions of pvar rows in Python would be far
    too slow.
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
    """Vectorized equivalent of `classify`, for use over a full pvar column."""
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


def _classify_df(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(_classify_expr(pl.col("REF"), pl.col("ALT")).alias("class"))


def _contig_stats(df: pl.DataFrame) -> list[dict]:
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
    """Inter-variant gaps (bp) within each contig, sorted-position diffs."""
    return (
        df.sort(["CHROM", "POS"])
        .select(pl.col("POS").diff().over("CHROM").alias("gap"))
        .drop_nulls()
        .filter(pl.col("gap") > 0)["gap"]
    )


def _class_counts(df: pl.DataFrame) -> dict[str, int]:
    counts = dict.fromkeys(CLASS_NAMES, 0)
    tally = df.group_by("class").agg(pl.len().alias("n"))
    for row in tally.iter_rows(named=True):
        counts[row["class"]] = row["n"]
    return counts


def _indel_lengths(df: pl.DataFrame) -> pl.Series:
    indels = df.filter(pl.col("class").is_in(["insertion", "deletion"]))
    # str.len_chars() is UInt32; subtracting two UInt32 columns wraps around
    # on underflow instead of going negative, so cast to a signed type first.
    alt_len = pl.col("ALT").str.len_chars().cast(pl.Int64)
    ref_len = pl.col("REF").str.len_chars().cast(pl.Int64)
    return indels.select((alt_len - ref_len).abs().alias("len"))["len"]


def _multiallelic_rate(df: pl.DataFrame) -> float:
    """Fraction of pvar rows whose ALT field lists more than one allele."""
    n = df.height
    if n == 0:
        return 0.0
    n_multi = df.select(pl.col("ALT").str.contains(",", literal=True).sum()).item()
    return n_multi / n


def _titv(df: pl.DataFrame) -> float:
    """Transition/transversion ratio over SNP rows (`_classify_df` must have run)."""
    snps = df.filter(pl.col("class") == "snp")
    if snps.height == 0:
        raise ValueError("no SNPs found; cannot compute Ti/Tv")
    pair = pl.concat_str([pl.col("REF"), pl.col("ALT")])
    n_ts = snps.select(pair.is_in(TRANSITION_PAIRS).sum()).item()
    n_tv = snps.height - n_ts
    if n_tv == 0:
        raise ValueError("no transversions found; cannot compute a finite Ti/Tv ratio")
    return n_ts / n_tv


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


def fit_sfs(pgen_prefix: str | Path, contigs: Iterable[str] | None = None) -> list[int]:
    """Shell out to `plink2 --freq counts` and return the ALT_CTS column.

    ALT_CTS is the observed non-reference allele count per site -- exactly
    the site-frequency-spectrum input the `sfs` histogram is fit from.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "sfs")
        args = ["--pfile", str(pgen_prefix), "--freq", "counts", "--out", out_prefix]
        if contigs:
            args += ["--chr", ",".join(contigs)]
        _run_plink2(args)
        df = pl.read_csv(f"{out_prefix}.acount", separator="\t")
        return df["ALT_CTS"].to_list()


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
    df = _classify_df(df)

    contigs = _contig_stats(df)
    if not contigs:
        raise ValueError("no variants found for the requested contig(s)")
    contig_ids = [c["id"] for c in contigs]

    gaps = _gaps(df)
    class_counts = _class_counts(df)
    indel_lens = _indel_lengths(df)
    titv = _titv(df)
    multiallelic_rate = _multiallelic_rate(df)

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
