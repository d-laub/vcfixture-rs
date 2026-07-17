#!/usr/bin/env python
"""Fit a vcfixture bulk-generation ``Profile`` JSON from real cohort data.

This script is the only place that turns a real cohort (1kGP high-coverage
pgen, GDC somatic pgen, a sites-only VCF such as the 1000 Genomes raw GATK
callset, ...) into the committed profile JSON consumed by the Rust bulk
generator (see ``src/bulk/profile.rs``). It is never imported by Rust and
has no contract other than that JSON schema.

Two mutually exclusive input modes are supported: ``--pgen`` (a plink2
fileset with genotypes, from which ``n_samples``, ``sfs``, ``missing_rate``,
and ``phased_rate`` are all fitted) and ``--sites-vcf`` (a sites-only
VCF/BCF with no genotype columns at all, only ``INFO/AC``/``INFO/AN`` --
``n_samples`` and ``phased_rate`` cannot be derived from it and must be
passed explicitly via ``--n-samples``/``--phased-rate``).

The profile schema deliberately separates ``fitted`` (statistics measured
from the source) from ``dialed`` (generation choices, e.g. FORMAT payload,
that are picked by the user and never claimed to be measured). Every value
this script writes under ``fitted`` is derived from the source data passed
in (or, for ``phased_rate`` under ``--sites-vcf``, passed explicitly since
no genotypes exist to fit it from); nothing here should ever be a
hand-picked literal presented as measured.

Every statistic in this module is computed as a *lazy* polars aggregation
whose collected result is small (a handful of rows, a histogram, a scalar)
-- never the full pvar/`.acount`/sites-VCF-TSV frame. Real cohorts are 75M
(1kGP germline) to 350M (somatic) pvar rows, or tens of millions of
sites-only VCF records per chromosome, with variable-length REF/ALT
strings; materializing that as an eager `pl.DataFrame` or a Python list is
what used to OOM-kill this script well before it finished (see git history
/ the design doc for the measured 20+ GB RSS). Keep every new statistic
lazy end-to-end: build a `pl.LazyFrame` pipeline, and only `.collect()` a
bounded-size result.

Usage
-----
    pixi run -e fit fit -- \\
        --pgen /path/to/prefix --name germline-1kgp --out profiles/germline-1kgp.json

    pixi run -e fit fit -- \\
        --sites-vcf /path/to/sites.vcf.gz --name germline-1kgp-sites \\
        --out profiles/germline-1kgp-sites.json \\
        --n-samples 3202 --phased-rate 0.999

See ``scripts/fit/README.md`` for the exact commands used to fit the two
committed profiles.
"""

from __future__ import annotations

import argparse
import atexit
import datetime as _dt
import json
import os
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import polars as pl

__version__ = "0.1.0"

# Enforced against the Rust enums by test_fit_profile.py
CLASS_NAMES = ("snp", "insertion", "deletion", "mnp", "complex", "symbolic")

# Enforced against the Rust enums by test_fit_profile.py
_PAYLOAD_CHOICES = ("gt-only", "gt-vaf", "gatk", "mutect2")


def _payload_choices() -> tuple[str, ...]:
    return _PAYLOAD_CHOICES


# ~90% of indels are <= 6 bp (see the design spec), so resolve that range
# finely and taper off for the long tail.
INDEL_EDGES = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 20.0, 50.0, 100.0, 1000.0]

# Log-spaced bins over [1, 1e5] for inter-variant gaps (bp).
GAP_LOW, GAP_HIGH, GAP_N_BINS = 1.0, 1e5, 24


# --------------------------------------------------------------------------
# Schema-facing helpers (unit tested directly, see test_fit_profile.py)
# --------------------------------------------------------------------------


def _validate_edges(edges: Sequence[float]) -> list[float]:
    edges = [float(e) for e in edges]
    if len(edges) < 2:
        raise ValueError("histogram needs >= 2 edges")
    if any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError("histogram edges must be strictly increasing")
    return edges


def _finalize_histogram(
    counts: Sequence[int], edges: Sequence[float], n_dropped: int
) -> dict[str, list[float]]:
    """Normalize per-bin `counts` (length ``len(edges) - 1``) into weights.

    Shared tail end of both the small-scale `histogram()` (which bins an
    already-materialized array with `numpy.histogram`) and the large-scale
    `histogram_lazy()` path (which bins hundreds of millions of rows via a
    lazy polars `group_by` and only ever collects the tiny per-bin counts) --
    the warning/normalization semantics must be identical either way.
    """
    edges = _validate_edges(edges)
    n_in_range = int(sum(counts))
    n_total = n_in_range + int(n_dropped)
    if n_dropped > 0:
        frac = n_dropped / n_total
        warnings.warn(
            f"histogram: dropped {n_dropped}/{n_total} values ({frac:.1%}) "
            f"outside the edge range [{edges[0]}, {edges[-1]}]",
            stacklevel=3,
        )
    if n_in_range <= 0:
        raise ValueError("no values fell within the histogram edges")
    weights = [c / n_in_range for c in counts]
    return {"edges": edges, "weights": weights}


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

    This bins an already-materialized `values` array with `numpy.histogram`
    -- fine for small inputs (tests, small cohorts) but NOT for the
    real-scale pvar-derived statistics (75M-350M rows), which go through
    `histogram_lazy` instead so the values are never collected into a single
    in-memory array.
    """
    edges = _validate_edges(edges)
    arr = np.asarray(values, dtype=np.float64)
    counts, _ = np.histogram(arr, bins=edges)
    n_dropped = int(arr.size) - int(counts.sum())
    return _finalize_histogram(counts.tolist(), edges, n_dropped)


def _bucket_index_expr(value: pl.Expr, edges: Sequence[float]) -> pl.Expr:
    """Vectorized equivalent of `numpy.histogram`'s bin assignment.

    Returns an `Int64` expression giving the 0-indexed bin each row's
    `value` falls into, or `null` if it falls outside ``[edges[0],
    edges[-1]]``. Bins are half-open ``[edges[i], edges[i+1])`` except the
    last, which is closed on both ends -- matching `numpy.histogram`'s
    convention exactly (a value equal to the final edge lands in the last
    bin, not out of range).

    This is the piece that lets histograms be computed over a `LazyFrame`
    with hundreds of millions of rows: the number of edges is small and
    static (tens, not millions), so the `when/then` chain built here is
    cheap regardless of how many rows the expression is evaluated over, and
    evaluating it never requires holding more than one row in memory at a
    time -- unlike `numpy.histogram`, which needs the full array.
    """
    edges = _validate_edges(edges)
    n_bins = len(edges) - 1
    if n_bins == 1:
        # single bin is closed on both ends, matching numpy.histogram
        return (
            pl.when((value >= edges[0]) & (value <= edges[1]))
            .then(pl.lit(0, dtype=pl.Int64))
            .otherwise(None)
        )
    expr = pl.when((value >= edges[0]) & (value < edges[1])).then(pl.lit(0, dtype=pl.Int64))
    for i in range(1, n_bins - 1):
        expr = expr.when((value >= edges[i]) & (value < edges[i + 1])).then(
            pl.lit(i, dtype=pl.Int64)
        )
    expr = expr.when((value >= edges[-2]) & (value <= edges[-1])).then(
        pl.lit(n_bins - 1, dtype=pl.Int64)
    )
    return expr.otherwise(None)


def histogram_lazy(lf: pl.LazyFrame, value: pl.Expr, edges: Sequence[float]) -> pl.LazyFrame:
    """Lazily bucketize `value` (evaluated over `lf`) into `edges` bins.

    Returns a `LazyFrame` with one row per populated bin index (`_bin`,
    `n`), plus possibly one row with `_bin == null` holding the count of
    out-of-range values. Collecting this yields at most ``len(edges)`` rows
    -- a handful to a few dozen -- no matter how many rows `lf` has, because
    the `group_by` reduces before anything is materialized. Pair with
    `_finalize_histogram_from_binned` to get the same
    ``{"edges": ..., "weights": ...}`` shape `histogram()` produces.
    """
    edges = _validate_edges(edges)
    return (
        lf.select(_bucket_index_expr(value, edges).alias("_bin"))
        .group_by("_bin")
        .agg(pl.len().alias("n"))
    )


def _finalize_histogram_from_binned(
    binned: pl.DataFrame, edges: Sequence[float]
) -> dict[str, list[float]]:
    """Turn a collected `histogram_lazy` result into the schema-facing dict."""
    edges = _validate_edges(edges)
    n_bins = len(edges) - 1
    counts = [0] * n_bins
    n_dropped = 0
    for b, n in binned.iter_rows():
        if b is None:
            n_dropped = int(n)
        else:
            counts[int(b)] = int(n)
    return _finalize_histogram(counts, edges, n_dropped)


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


def _count_leading_meta_lines(path: str | Path) -> int:
    """Count leading ``##`` metadata lines so scan_csv can `skip_lines` past
    them instead of `comment_prefix='##'` (which forces polars to materialise
    the whole file -- measured 6.2 GB / 3x slower on a 5.9 GB somatic pvar).

    Handles the pvar's optional .zst compression: the meta block is a small
    text prefix, so only the first chunk needs decompressing. The `.zst`
    branch requires the optional `zstandard` package (not a hard dependency
    of the `fit` env); it is only imported when a `.zst` path is passed.
    """
    p = str(path)
    if p.endswith(".zst"):
        import zstandard  # in the `fit` env

        n = 0
        with open(p, "rb") as raw:
            dctx = zstandard.ZstdDecompressor()
            with dctx.stream_reader(raw) as r:
                buf = r.read(1 << 16).split(b"\n")
                for line in buf:
                    if line.startswith(b"##"):
                        n += 1
                    else:
                        break
        return n
    n = 0
    with open(p, "rb") as fh:
        for line in fh:
            if line.startswith(b"##"):
                n += 1
            else:
                break
    return n


def read_pvar(path: str | Path) -> pl.LazyFrame:
    """Lazily scan a `.pvar` or `.pvar.zst` file, keeping only #CHROM/POS/ID/REF/ALT.

    A 1kGP-scale pvar can be 500+ MB of compressed text (75M+ rows); the
    somatic cohort's is 5.9 GB uncompressed (348M rows) with REF/ALT strings
    up to 60+ characters. This returns an *unmaterialized* `pl.LazyFrame` --
    callers must keep every downstream statistic expressed as a lazy
    aggregation (`compute_pvar_stats` does this) and only collect small,
    bounded results. Calling `.collect()` on the frame this returns without
    first reducing it (e.g. via `group_by`/aggregation) will eagerly pull
    every row into memory and reproduce the OOM this function exists to
    avoid -- a streaming `.collect(engine="streaming")` streams the
    *computation*, not the *result*: the result is still a fully
    materialized `DataFrame` holding every row. `.zst` decoding is handled
    natively by polars based on the file extension either way.

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
    n_meta = _count_leading_meta_lines(path)
    lf = pl.scan_csv(
        str(path),
        separator="\t",
        skip_lines=n_meta,
        schema_overrides={"#CHROM": pl.Utf8, "ID": pl.Utf8, "REF": pl.Utf8, "ALT": pl.Utf8},
    ).rename({"#CHROM": "CHROM"})
    # ID is carried only to keep read_pvar/read_sites_vcf schemas identical
    # for compute_pvar_stats; nothing reads it, and in this file it is "."
    # (REF/ALT carry the wide strings), so dropping it saves ~0 memory --
    # keep it for schema parity with read_sites_vcf's fabricated ID.
    return lf.select(["CHROM", "POS", "ID", "REF", "ALT"])


def read_sites_vcf(path: str | Path, contigs: Iterable[str] | None = None) -> pl.LazyFrame:
    """Lazily read a sites-only VCF/BCF's PASS records via `bcftools query`.

    This is the sites-only-VCF analogue of `read_pvar`, for cohorts with no
    genotype columns at all -- e.g. the 1000 Genomes raw GATK callset, the
    only available source of a realistic 1kGP allele-frequency spectrum
    (the phased panel already fit has 0% singletons, since phasing drops
    unphaseable singletons). It carries `INFO/AC` and `INFO/AN` instead of
    per-sample calls. It is also VQSR-tranched: ~14-16% of raw records carry
    a `VQSRTrancheSNP99.80to100.00` / `VQSRTrancheINDEL99.00to100.00` /
    `LowQual` FILTER, so the `bcftools query -i 'FILTER="PASS"'` below is the
    cleaning step here, not an optimization -- those records must never
    reach the profile.

    `bcftools query` streams its TSV output directly to a file on disk (via
    `_run_bcftools`'s `stdout_path`), never through a materialized Python
    string, and this function returns an *unmaterialized* `pl.LazyFrame`
    scanning that file -- exactly the memory contract `read_pvar`'s
    docstring documents at length, for the same reason: a 1kGP-chromosome
    sites-only VCF is ~337 MB with tens of millions of PASS records, and a
    prior version of this script was OOM-killed at 20.8 GB for eagerly
    materializing data at this scale. Callers must keep every downstream
    statistic expressed as a lazy aggregation and only collect small,
    bounded results -- never call `.collect()` on the frame returned here
    without first reducing it.

    Returns columns `CHROM (Utf8), POS (Int64), ID (Utf8), REF (Utf8), ALT
    (Utf8), AC (Int64), AN (Int64)`. `ID` is always the literal `"."` --
    sites-only callsets carry no meaningful per-record ID, but `pvar`'s
    schema (and `compute_pvar_stats`, which this frame is fed into via
    `.select(["CHROM", "POS", "ID", "REF", "ALT"])`) has one, so a
    placeholder keeps the two readers' output interchangeable there.
    Multiallelic sites are assumed already split into separate biallelic
    rows upstream (true of the 1kGP raw callset), each with its own scalar
    `AC` and a shared `AN` -- unlike `read_pvar`'s comma-joined `ALT`, no
    `_explode_alleles` step is needed or possible here.

    The TSV lives in a directory made with `tempfile.mkdtemp()`, not a
    `tempfile.TemporaryDirectory()` context manager, deliberately: this
    function hands back an *unresolved* `LazyFrame` that still needs to scan
    that file whenever a caller eventually collects it, so the directory
    must outlive this function's return -- a context manager would delete
    the TSV out from under the returned frame the moment `read_sites_vcf`
    returns. Instead, cleanup is registered with `atexit`: the directory is
    removed when the process exits, which is exactly the returned
    `LazyFrame`'s lifetime (this script runs once per `fit` invocation, not
    as a long-lived process), so nothing needs to collect the frame early
    just to let the directory be deleted, and a real run's hundreds-of-MB
    TSV no longer lingers in the system temp dir waiting for the OS to age
    it out.
    """
    d = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, d, ignore_errors=True)
    tsv_path = Path(d) / "sites.tsv"
    args = ["query", "-i", 'FILTER="PASS"', "-f", "%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n"]
    if contigs:
        args += ["-r", ",".join(contigs)]
    args.append(str(path))
    _run_bcftools(args, stdout_path=tsv_path)

    # CHROM must stay a string for the same reason as read_pvar's #CHROM:
    # numeric-looking contig ids like "1"/"22" must not be inferred as
    # integers -- ContigStat.id is a String in src/bulk/profile.rs.
    lf = pl.scan_csv(
        str(tsv_path),
        separator="\t",
        has_header=False,
        new_columns=["CHROM", "POS", "REF", "ALT", "AC", "AN"],
        schema_overrides={
            "CHROM": pl.Utf8,
            "POS": pl.Int64,
            "REF": pl.Utf8,
            "ALT": pl.Utf8,
            "AC": pl.Int64,
            "AN": pl.Int64,
        },
    )
    return lf.with_columns(pl.lit(".").alias("ID")).select(
        ["CHROM", "POS", "ID", "REF", "ALT", "AC", "AN"]
    )


def build_profile(
    *,
    name: str,
    source: str,
    n_samples: int,
    contigs: list[dict],
    gap_dist: dict,
    sfs: dict,
    indel_length: dict,
    class_counts: Mapping[str, int],
    titv: float,
    multiallelic_rate: float,
    missing_rate: float,
    phased_rate: float,
    ploidy: int,
    supplied: list[str],
    payload: str = "gt-only",
    n_variants_source: int | None = None,
) -> dict:
    """Assemble a schema-valid Profile dict (see `src/bulk/profile.rs`).

    `gap_dist`, `sfs`, and `indel_length` are already-finalized histogram
    dicts (``{"edges": [...], "weights": [...]}``, as produced by
    `histogram`/`histogram_lazy` + `_finalize_histogram*`) rather than raw
    value sequences: at real scale those sequences would be 75M-350M
    elements, so the histogram must already have been reduced by the caller
    (`main()`, via lazy polars aggregations) before it ever reaches this
    function. Every field under "fitted" is otherwise derived from
    `class_counts`/`contigs`/the scalar rates passed in by the caller --
    never hand-picked here. `dialed.payload` and `dialed.ploidy` are the
    deliberate exceptions: they are generation choices, not fitted
    statistics.

    `supplied` names every field the caller passed in rather than measured
    from the source data (e.g. `ploidy` always; `phased_rate`/`n_samples`
    too when fitting from a sites-only VCF, since neither is derivable from
    one) -- it makes which values in "fitted"/"dialed" are hand-supplied,
    not measured, auditable from the JSON alone. No default: callers must
    say explicitly.
    """
    if n_variants_source is None:
        n_variants_source = sum(c["n_variants"] for c in contigs)

    fitted = {
        "contigs": contigs,
        "gap_dist": gap_dist,
        "sfs": sfs,
        "variant_classes": class_mix_from_counts(class_counts),
        "indel_length": indel_length,
        "titv": titv,
        "multiallelic_rate": multiallelic_rate,
        "missing_rate": missing_rate,
        "phased_rate": phased_rate,
    }
    return {
        "name": name,
        "provenance": {
            "source": source,
            "n_samples_source": n_samples,
            "n_variants_source": n_variants_source,
            "fitted_on": _dt.date.today().isoformat(),
            "fit_tool_version": __version__,
            "supplied": sorted(supplied),
        },
        "fitted": fitted,
        "dialed": {"payload": payload, "ploidy": ploidy},
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
# Extraction from a pvar LazyFrame (vectorized polars, no Python row loops)
# --------------------------------------------------------------------------


def _explode_alleles(lf: pl.LazyFrame) -> pl.LazyFrame:
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

    Vectorized via `Expr.str.split` + `LazyFrame.explode` -- no Python loop
    over rows, and this stays lazy: nothing is materialized until a
    downstream `.collect()` on a reduced (e.g. `group_by`) result.
    """
    return lf.with_columns(pl.col("ALT").str.split(",")).explode("ALT", empty_as_null=False)


def _classify_df(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.with_columns(_classify_expr(pl.col("REF"), pl.col("ALT")).alias("class"))


def _contig_stats_lazy(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Per-contig record count, span, and position bounds. `lf` is per-RECORD (not exploded).

    Collecting this yields one row per contig (tens, not millions) --
    `pos_min` is kept alongside the public `id`/`n_variants`/`density_per_kb`
    fields purely so `main()` can look up the first fitted contig's start
    position for `fit_phased_rate` without an extra full-file scan; it is
    dropped before the profile's `contigs` list (see `_contig_rows_from_df`)
    since it isn't part of the `ContigStat` schema in `src/bulk/profile.rs`.
    """
    return (
        lf.group_by("CHROM")
        .agg(
            n_variants=pl.len(),
            pos_min=pl.col("POS").min(),
            pos_max=pl.col("POS").max(),
        )
        .with_columns(span_bp=(pl.col("pos_max") - pl.col("pos_min")).clip(lower_bound=1))
        .with_columns(density_per_kb=pl.col("n_variants") / (pl.col("span_bp") / 1000.0))
        .sort("CHROM")
    )


def _contig_rows_from_df(df: pl.DataFrame) -> list[dict]:
    return [
        {
            "id": row["CHROM"],
            "n_variants": row["n_variants"],
            "density_per_kb": row["density_per_kb"],
        }
        for row in df.iter_rows(named=True)
    ]


def _contig_pos_min_from_df(df: pl.DataFrame) -> dict[str, int]:
    return {row["CHROM"]: row["pos_min"] for row in df.iter_rows(named=True)}


def _gap_bins_lazy(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Inter-variant gap (bp) histogram bin counts, within each contig.

    `lf` is per-RECORD (not exploded): a multiallelic site is still one
    position, so exploding first would fabricate spurious zero-length gaps
    between the alleles of the same site.

    `lf` is assumed coordinate-sorted within each contig (plink2 emits pvar
    sorted by CHROM, POS; guarded by `assert_pvar_sorted`). Gaps are a
    straight `POS.diff()` masked to same-contig adjacent rows (`CHROM ==
    CHROM.shift(1)`), NOT `sort().diff().over("CHROM")`: the sort is a
    full-frame pipeline breaker and `.over()` is a non-streaming window,
    which together cost ~20 GB at genome scale (348M rows). shift+mask
    streams in ~220 MB and is bit-identical on sorted input.

    Returns the `histogram_lazy` bin-count frame directly (bounded by
    `GAP_N_BINS`), never the underlying gaps themselves -- a 1kGP-scale
    contig can have tens of millions of gaps, which is exactly the array
    `read_pvar`'s docstring warns against materializing.
    """
    gaps = (
        lf.select(
            pl.col("POS").diff().alias("gap"),
            (pl.col("CHROM") == pl.col("CHROM").shift(1)).alias("same_contig"),
        )
        .filter(
            pl.col("same_contig")
            & pl.col("gap").is_not_null()
            & (pl.col("gap") > 0)
        )
        .select("gap")
    )
    return histogram_lazy(gaps, pl.col("gap"), _gap_edges())


def assert_pvar_sorted(lf: pl.LazyFrame) -> None:
    """Fail if any within-contig POS is out of order (the precondition
    `_gap_bins_lazy` relies on after dropping its sort). Streams: one boolean
    reduction, bounded memory.

    Caveat: this checks only descending POS between *adjacent same-contig*
    rows. It does not detect a contig split into non-contiguous blocks
    (interleaved contigs), which would also break the sort-free gap path --
    plink2 always emits contig-grouped pvar, so that case does not arise
    for our inputs."""
    bad = (
        lf.select(
            (pl.col("CHROM") == pl.col("CHROM").shift(1)).alias("same"),
            (pl.col("POS") < pl.col("POS").shift(1)).alias("descending"),
        )
        .select((pl.col("same") & pl.col("descending")).sum().alias("n_desc"))
        .collect(engine="streaming")
        .item()
    )
    if bad:
        raise ValueError(
            f"pvar is not sorted within contigs ({bad} descending steps); "
            "gap fitting requires coordinate-sorted input"
        )


def _class_counts_from_df(df: pl.DataFrame) -> dict[str, int]:
    counts = dict.fromkeys(CLASS_NAMES, 0)
    for row in df.iter_rows(named=True):
        counts[row["class"]] = row["n"]
    return counts


def _indel_bins_lazy(alleles: pl.LazyFrame) -> pl.LazyFrame:
    """Indel-length histogram bin counts over a per-ALLELE frame (post `_explode_alleles` + `_classify_df`)."""
    indels = alleles.filter(pl.col("class").is_in(["insertion", "deletion"]))
    # str.len_chars() is UInt32; subtracting two UInt32 columns wraps around
    # on underflow instead of going negative, so cast to a signed type first.
    alt_len = pl.col("ALT").str.len_chars().cast(pl.Int64)
    ref_len = pl.col("REF").str.len_chars().cast(pl.Int64)
    length = (alt_len - ref_len).abs()
    return histogram_lazy(indels, length, INDEL_EDGES)


def _multiallelic_rate_lazy(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Record count and multiallelic-record count (both scalars), per-RECORD frame.

    Must run on the un-exploded, per-record frame: this counts sites, not
    alleles, so a triallelic site still contributes exactly 1 to both the
    numerator and denominator, not 2.
    """
    return lf.select(
        n=pl.len(),
        n_multi=pl.col("ALT").str.contains(",", literal=True).sum(),
    )


def _multiallelic_rate_from_row(n: int, n_multi: int) -> float:
    return (n_multi / n) if n else 0.0


def _titv_lazy(alleles: pl.LazyFrame) -> pl.LazyFrame:
    """SNP count and transition count (both scalars) over a per-ALLELE frame.

    Transitions are the four purine<->purine / pyrimidine<->pyrimidine pairs,
    expressed as direct (REF,ALT) comparisons rather than
    `concat_str([REF,ALT]).is_in(TRANSITION_PAIRS)`: measured, the `is_in`
    against a literal list costs ~10 GB extra at 348M rows (16.7 GB vs the
    6.4 GB scan floor) for a bit-identical result.
    """
    snps = alleles.filter(pl.col("class") == "snp")
    r, a = pl.col("REF"), pl.col("ALT")
    is_ts = (
        ((r == "A") & (a == "G"))
        | ((r == "G") & (a == "A"))
        | ((r == "C") & (a == "T"))
        | ((r == "T") & (a == "C"))
    )
    return snps.select(n_snps=pl.len(), n_ts=is_ts.sum())


def _titv_from_row(n_snps: int, n_ts: int) -> float:
    if n_snps == 0:
        raise ValueError("no SNPs found; cannot compute Ti/Tv")
    n_tv = n_snps - n_ts
    if n_tv == 0:
        raise ValueError("no transversions found; cannot compute a finite Ti/Tv ratio")
    return n_ts / n_tv


def compute_pvar_stats(lf: pl.LazyFrame) -> dict:
    """Compute every `fitted` statistic derived purely from the pvar frame.

    `lf` is the raw, per-record `pl.LazyFrame` from `read_pvar` (optionally
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

    Every one of these is built as a lazy query whose *collected* output is
    small (contig count, histogram bin counts, or a one-row scalar frame) --
    never the full per-record or per-allele frame. All six queries below
    share two common subplans (the raw `lf` scan and the exploded `alleles`
    scan), but they are collected one at a time rather than via
    `pl.collect_all`: `collect_all`'s common-subplan elimination *caches*
    the shared 348M-row scan/explode output so it can be reused across the
    six plans, costing +20 GB (26.4 GB vs 6.4 GB) for no real gain, since a
    re-scan from warm page cache is both leaner and faster than holding a
    materialised copy around.
    """
    assert_pvar_sorted(lf)

    contig_lf = _contig_stats_lazy(lf)
    gap_bins_lf = _gap_bins_lazy(lf)
    multi_lf = _multiallelic_rate_lazy(lf)

    alleles = _classify_df(_explode_alleles(lf))
    class_lf = alleles.group_by("class").agg(pl.len().alias("n"))
    indel_bins_lf = _indel_bins_lazy(alleles)
    titv_lf = _titv_lazy(alleles)

    # Collect each plan on its own rather than pl.collect_all: collect_all's
    # common-subplan elimination *caches* the shared 348M-row scan/explode
    # output, costing +20 GB (26.4 GB vs 6.4 GB) and buying nothing -- a
    # re-scan from warm page cache is both leaner and faster than holding a
    # materialised copy. Measured on the 348M-row somatic pvar.
    contig_df = contig_lf.collect(engine="streaming")
    gap_bins_df = gap_bins_lf.collect(engine="streaming")
    multi_df = multi_lf.collect(engine="streaming")
    class_df = class_lf.collect(engine="streaming")
    indel_bins_df = indel_bins_lf.collect(engine="streaming")
    titv_df = titv_lf.collect(engine="streaming")

    n, n_multi = multi_df.row(0)
    n_snps, n_ts = titv_df.row(0)

    return {
        "contigs": _contig_rows_from_df(contig_df),
        "contig_pos_min": _contig_pos_min_from_df(contig_df),
        "gap_dist": _finalize_histogram_from_binned(gap_bins_df, _gap_edges()),
        "multiallelic_rate": _multiallelic_rate_from_row(n, n_multi),
        "class_counts": _class_counts_from_df(class_df),
        "indel_length": _finalize_histogram_from_binned(indel_bins_df, INDEL_EDGES),
        "titv": _titv_from_row(n_snps, n_ts),
    }


# --------------------------------------------------------------------------
# plink2 subprocess helpers
# --------------------------------------------------------------------------

# plink2 sizes its working-memory ("bigstack") off the *node's* physical RAM
# as reported by the OS, not off the cgroup it is actually confined to -- on
# a shared SLURM node with e.g. a 32 GiB per-job cgroup but ~1 TB of node
# RAM, plink2 auto-detects ~half the node total (~500 GB) and plans huge
# genotype blocks accordingly. Running `--freq counts` genome-wide (348M
# variants x 16k samples) it then grew resident memory past the cgroup limit
# (~23 GiB) and was OOM-killed. Every plink2 call here is block-processed
# (`--freq`/`--missing`/`--export vcf`), so `--memory` only changes the block
# size, never the output -- capping it keeps resident memory bounded.
#
# The cap must be large enough, not just small: plink2's bigstack is a
# sparsely-resident bump allocator. For 348M variants it *reserves* ~23 GiB
# up front for the (worst-case-sized) variant index, but only ~15 GiB ever
# becomes resident. Caps below that reservation fail early with plink2's own
# "Out of memory" -- measured on this fileset: 8000/16000/20000/26000 MiB all
# fail (26000 leaves < 2.8 GiB, too little for the ~2.8 GiB pgen block),
# while 30000 MiB succeeds with plink2 RSS ~15 GiB and total job-cgroup anon
# RSS ~23 GiB (agent + plink2), ~11 GiB under the 32 GiB limit. The 30000 MiB
# is a *virtual* reservation that only materializes what the data needs, so
# it is harmless on smaller inputs. The polars stages elsewhere in this
# script (~19 GiB peak python RSS on the genome-wide fit) run before/between
# the plink2 calls, not concurrently with plink2's bigstack.
#
# The default is tuned for a 32 GiB cgroup; override via the
# VCFIXTURE_PLINK2_MEMORY_MB env var for a different job allocation.
_PLINK2_MEMORY_MB = int(os.environ.get("VCFIXTURE_PLINK2_MEMORY_MB", "30000"))


def _run_plink2(args: list[str]) -> None:
    cmd = ["plink2", *args, "--memory", str(_PLINK2_MEMORY_MB)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"plink2 {' '.join(cmd[1:])} failed (exit {result.returncode}):\n"
            f"{result.stdout}\n{result.stderr}"
        )


def _run_bcftools(args: list[str], stdout_path: str | Path) -> None:
    """Shell out to `bcftools`, mirroring `_run_plink2`'s error-handling style.

    bcftools' stdout is streamed directly to `stdout_path` by the subprocess
    itself (`stdout=<file handle>`), never captured through a Python string
    first -- required (not optional) for exactly the reason `read_pvar`'s
    docstring explains at length: a 1kGP-scale sites-only VCF's PASS records
    can be tens of millions of rows, and `subprocess.run(capture_output=True)`
    would materialize the whole TSV as one Python string before polars ever
    saw it. Every caller streams to a file, so there is no in-memory-capture
    mode to keep around.
    """
    with open(stdout_path, "wb") as out_fh:
        result = subprocess.run(["bcftools", *args], stdout=out_fh, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"bcftools {' '.join(args)} failed (exit {result.returncode}):\n"
            f"{result.stderr.decode(errors='replace')}"
        )


def _pfile_args(pgen_prefix: str | Path, vzs: bool) -> list[str]:
    """Build the `--pfile` argument list, adding plink2's `vzs` marker when needed.

    plink2 does NOT auto-detect a `.pvar.zst` the way `read_pvar` does: if
    only the `.zst`-compressed pvar exists (no plain `.pvar` alongside it,
    which is how the committed 1kGP/somatic filesets are laid out), a bare
    `--pfile <prefix>` fails with "Failed to open <prefix>.pvar : No such
    file or directory" -- `vzs` must be passed as a second positional token
    after the prefix to tell plink2 to look for `.pvar.zst` instead.
    """
    args = ["--pfile", str(pgen_prefix)]
    if vzs:
        args.append("vzs")
    return args


def _sfs_from_ac(lf: pl.LazyFrame, ac_col: str, n_samples: int) -> dict:
    """Reduce a LazyFrame's per-allele allele-count column to the finalized `sfs` histogram.

    The shared tail of both `fit_sfs` (plink2 `.acount`'s split `ALT_CTS`)
    and `fit_sfs_from_sites_vcf` (a sites-only VCF's `INFO/AC`, already one
    scalar count per allele since multiallelics are pre-split upstream) --
    both need only the same `_sfs_edges` + `histogram_lazy` reduction, so
    this factors it out rather than duplicating it. `lf` must already carry
    one row per allele observation in `ac_col`; the caller is responsible
    for whatever split/explode got it there. Stays lazy end-to-end except
    for the final `.collect()`, whose result is bounded by the number of
    histogram bins, not the number of input rows.
    """
    edges = _sfs_edges(n_samples)
    binned = histogram_lazy(lf, pl.col(ac_col), edges).collect(engine="streaming")
    return _finalize_histogram_from_binned(binned, edges)


def fit_sfs(
    pgen_prefix: str | Path,
    n_samples: int,
    contigs: Iterable[str] | None = None,
    vzs: bool = False,
) -> dict:
    """Shell out to `plink2 --freq counts` and return the fitted `sfs` histogram.

    ALT_CTS is the observed non-reference allele count per site -- exactly
    the site-frequency-spectrum input the `sfs` histogram is fit from. For
    multiallelic sites plink2 keeps ALT_CTS comma-joined (aligned with the
    equally comma-joined ALT column); this splits that into one count per
    allele, lazily via `pl.scan_csv` + `histogram_lazy`, and returns the
    finalized histogram directly -- the `.acount` file has one row per pvar
    record (up to ~350M for the somatic cohort), so collecting every count
    into a Python list first would materialize hundreds of millions of
    floats.

    The temporary directory holding `.acount` is deleted once this function
    returns, so the full lazy-scan-to-histogram pipeline (including the
    final `.collect()`, which only ever produces a few dozen bin-count rows)
    must run inside the `with tempfile.TemporaryDirectory()` block below --
    returning an uncollected `LazyFrame` scanning that file would go stale.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "sfs")
        # --nonfounders: plink2's --freq counts defaults to founders-only,
        # which would silently compute ALT_CTS over fewer than the
        # `n_samples` used to build `_sfs_edges` (a real cohort has related
        # individuals / nonfounders present, e.g. 1kGP is 2583 founders out
        # of 3202 total samples) -- without this, plink2 refuses to run at
        # all ("--freq counts specified, but with neither --ac-founders nor
        # --nonfounders; and nonfounders are present").
        args = _pfile_args(pgen_prefix, vzs) + [
            "--freq",
            "counts",
            "--nonfounders",
            "--out",
            out_prefix,
        ]
        if contigs:
            args += ["--chr", ",".join(contigs)]
        _run_plink2(args)
        lf = pl.scan_csv(
            f"{out_prefix}.acount", separator="\t", schema_overrides={"ALT_CTS": pl.Utf8}
        )
        alt_cts = (
            lf.select(pl.col("ALT_CTS").str.split(","))
            .explode("ALT_CTS", empty_as_null=False)
            .select(pl.col("ALT_CTS").cast(pl.Float64))
        )
        return _sfs_from_ac(alt_cts, "ALT_CTS", n_samples)


def fit_missing_rate(
    pgen_prefix: str | Path, contigs: Iterable[str] | None = None, vzs: bool = False
) -> float:
    """Shell out to `plink2 --missing` and return the global hardcall missing rate.

    `.vmiss` also has one row per pvar record, so this scans it lazily and
    sums the two count columns via a single reducing `.select()` -- the
    collected result is one row -- rather than reading the whole file
    eagerly with `pl.read_csv`.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "miss")
        args = _pfile_args(pgen_prefix, vzs) + [
            "--missing",
            "variant-only",
            "--out",
            out_prefix,
        ]
        if contigs:
            args += ["--chr", ",".join(contigs)]
        _run_plink2(args)
        totals = (
            pl.scan_csv(f"{out_prefix}.vmiss", separator="\t")
            .select(n_missing=pl.col("MISSING_CT").sum(), n_obs=pl.col("OBS_CT").sum())
            .collect(engine="streaming")
        )
        n_missing, n_obs = totals.row(0)
        return n_missing / n_obs if n_obs else 0.0


def fit_phased_rate(
    pgen_prefix: str | Path,
    contig: str,
    pos_min: int,
    window_bp: int = 1_000_000,
    vzs: bool = False,
) -> float:
    """Estimate the fraction of genotype calls that are phased.

    pgen has no direct "phased fraction" report, so this exports a bounded
    window (`window_bp` starting at `pos_min` on `contig`) to VCF and counts
    the phased ('|') vs unphased ('/') GT separators plink2 actually wrote,
    across every sample x variant call in that window. The window (default
    1 Mb) is bounded regardless of cohort size, so the line-by-line VCF
    parse below stays small and isn't a materialization risk the way an
    unbounded pvar/`.acount` scan would be.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "phase")
        _run_plink2(
            _pfile_args(pgen_prefix, vzs)
            + [
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
        if n_total == 0:
            raise ValueError(
                "no genotype calls found in the phase-sampling window; "
                "cannot estimate phased_rate"
            )
        return n_phased / n_total


# --------------------------------------------------------------------------
# Sites-only VCF (bcftools, no genotypes): the two genotype-derived stats
# that a sites-only source can still support from INFO/AC and INFO/AN.
# `phased_rate` has no such substitute -- there are no genotype calls at
# all to inspect -- so it is taken verbatim from `--phased-rate` in `main`.
# --------------------------------------------------------------------------


def fit_sfs_from_sites_vcf(lf: pl.LazyFrame, n_samples: int) -> dict:
    """Build the fitted `sfs` histogram from a `read_sites_vcf` frame's `AC` column.

    Unlike the pgen path (`fit_sfs`), a sites-only VCF's `INFO/AC` is
    already one scalar count per allele -- multiallelic sites are pre-split
    upstream (true of the 1kGP raw callset) -- so no split/explode step is
    needed before handing it to the shared `_sfs_from_ac` reduction. Uses
    the same `_sfs_edges(n_samples)` edges as the pgen path, so profiles
    fitted from either source stay directly comparable.
    """
    return _sfs_from_ac(lf.select(pl.col("AC").cast(pl.Float64)), "AC", n_samples)


def fit_missing_rate_from_sites_vcf(lf: pl.LazyFrame, n_samples: int) -> float:
    """Estimate the missing rate from a `read_sites_vcf` frame's `AN` column.

    A sites-only VCF has no per-sample hardcalls to count missing directly
    the way `fit_missing_rate`'s plink2 `.vmiss` does, but `INFO/AN` -- the
    number of alleles actually called at each site -- gives the same rate
    indirectly: at full ploidy-2 calling every site would have
    `AN == 2 * n_samples`, so the shortfall from that maximum, averaged
    over all sites, is the missing rate. Computed as a single lazy
    `.mean()` aggregation collecting one scalar, never the full `AN`
    column. Clamped to `[0, 1]`, since a mean `AN` above `2 * n_samples`
    (e.g. from a mismatched `--n-samples`) would otherwise yield a negative
    rate.
    """
    row = lf.select(mean_an=pl.col("AN").mean()).collect(engine="streaming").row(0)
    mean_an = row[0]
    if mean_an is None:
        return 0.0
    rate = 1.0 - (mean_an / (2 * n_samples))
    return min(1.0, max(0.0, rate))


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


def _fit_from_pgen(args: argparse.Namespace) -> dict:
    """Fit a profile from a plink2 pgen/pvar/psam fileset (the original, genotyped path)."""
    prefix = str(args.pgen)
    psam_path = Path(f"{prefix}.psam")
    pvar_path = Path(f"{prefix}.pvar.zst")
    vzs = pvar_path.exists()
    if not vzs:
        pvar_path = Path(f"{prefix}.pvar")
    if not pvar_path.exists():
        raise FileNotFoundError(f"no .pvar or .pvar.zst found for prefix {prefix}")

    n_samples = _read_n_samples(psam_path)

    lf = read_pvar(pvar_path)
    if args.contigs:
        lf = lf.filter(pl.col("CHROM").is_in(args.contigs))

    stats = compute_pvar_stats(lf)
    contigs = stats["contigs"]
    if not contigs:
        raise ValueError("no variants found for the requested contig(s)")
    contig_ids = [c["id"] for c in contigs]

    chr_filter = contig_ids if args.contigs else None
    sfs = fit_sfs(args.pgen, n_samples, contigs=chr_filter, vzs=vzs)
    missing_rate = fit_missing_rate(args.pgen, contigs=chr_filter, vzs=vzs)

    first_contig = contigs[0]
    pos_min = int(stats["contig_pos_min"][first_contig["id"]])
    phased_rate = fit_phased_rate(
        args.pgen,
        contig=first_contig["id"],
        pos_min=pos_min,
        window_bp=int(args.phase_sample_mb * 1_000_000),
        vzs=vzs,
    )

    return build_profile(
        name=args.name,
        source=prefix,
        n_samples=n_samples,
        contigs=contigs,
        gap_dist=stats["gap_dist"],
        sfs=sfs,
        indel_length=stats["indel_length"],
        class_counts=stats["class_counts"],
        titv=stats["titv"],
        multiallelic_rate=stats["multiallelic_rate"],
        missing_rate=missing_rate,
        phased_rate=phased_rate,
        ploidy=args.ploidy,
        payload=args.payload,
        supplied=["ploidy"],
    )


def _fit_from_sites_vcf(args: argparse.Namespace) -> dict:
    """Fit a profile from a sites-only VCF: no genotypes, so `--n-samples` and
    `--phased-rate` (validated as required in `main` before this runs) stand
    in for what `_fit_from_pgen` derives from the pgen fileset itself.
    """
    lf = read_sites_vcf(args.sites_vcf, contigs=args.contigs)

    # compute_pvar_stats is the shared stats core (contigs, gap_dist,
    # variant_classes, indel_length, titv, multiallelic_rate) -- it only
    # needs the CHROM/POS/ID/REF/ALT columns `read_pvar` also produces.
    stats = compute_pvar_stats(lf.select(["CHROM", "POS", "ID", "REF", "ALT"]))
    contigs = stats["contigs"]
    if not contigs:
        raise ValueError("no variants found for the requested contig(s)")

    sfs = fit_sfs_from_sites_vcf(lf, args.n_samples)
    missing_rate = fit_missing_rate_from_sites_vcf(lf, args.n_samples)

    return build_profile(
        name=args.name,
        source=str(args.sites_vcf),
        n_samples=args.n_samples,
        contigs=contigs,
        gap_dist=stats["gap_dist"],
        sfs=sfs,
        indel_length=stats["indel_length"],
        class_counts=stats["class_counts"],
        titv=stats["titv"],
        multiallelic_rate=stats["multiallelic_rate"],
        missing_rate=missing_rate,
        phased_rate=args.phased_rate,
        ploidy=args.ploidy,
        payload=args.payload,
        supplied=["ploidy", "phased_rate", "n_samples"],
    )


def _validate_with_rust(path: Path) -> None:
    """Self-check a freshly-written profile against `Profile::validate`.

    Runs ``cargo run -q --features bulk --bin validate-profile -- <path>``
    (from the repo root, so it works regardless of the caller's cwd) --
    the exact same `Profile::from_json` + `Profile::validate` a profile
    goes through once the crate embeds it via ``include_str!``. This turns
    a bad fit (e.g. `ploidy == 0`, a rate like `missing_rate` outside
    [0, 1], a NaN histogram bin) into an immediate failure here instead of
    a much-later failure inside the Rust crate.

    Behavior:
    - If ``cargo`` is not found on ``PATH`` (e.g. a Python-only sandbox
      with no Rust toolchain), prints a warning and returns without
      running anything -- this check is a courtesy, not a hard dependency
      of fitting a profile.
    - If the validator exits non-zero, raises ``SystemExit`` carrying its
      stderr as the message, which aborts `main()` with exit status 1 and
      prints that stderr for the caller to see.
    - On success (exit 0), returns silently.
    """
    if shutil.which("cargo") is None:
        warnings.warn(
            f"cargo not found on PATH; skipping Rust validation of {path}",
            stacklevel=2,
        )
        return
    repo_root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        ["cargo", "run", "-q", "--features", "bulk", "--bin", "validate-profile", "--", str(path)],
        cwd=repo_root,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise SystemExit(f"validate-profile rejected {path}:\n{result.stderr}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Fit a vcfixture bulk-generation profile JSON from a real cohort."
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--pgen",
        help="plink2 fileset prefix (expects <prefix>.pgen/.psam/.pvar[.zst])",
    )
    source_group.add_argument(
        "--sites-vcf",
        help=(
            "sites-only VCF/BCF path (no genotype columns; requires "
            "--n-samples and --phased-rate, since neither is derivable "
            "from a sites-only file)"
        ),
    )
    parser.add_argument("--name", required=True, help="profile name, e.g. germline-1kgp")
    parser.add_argument("--out", required=True, type=Path, help="output profile JSON path")
    parser.add_argument(
        "--payload",
        default="gt-only",
        choices=_PAYLOAD_CHOICES,
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
        default=None,
        help=(
            "window size (Mb, default 1.0), from the first fitted contig's "
            "start, used to estimate phased_rate (--pgen only)"
        ),
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=None,
        help="cohort size (required with --sites-vcf; derived from --pgen's .psam otherwise)",
    )
    parser.add_argument(
        "--phased-rate",
        type=float,
        default=None,
        help=(
            "fixed phased_rate in [0, 1] (required with --sites-vcf, which "
            "has no genotypes to count phased/unphased calls from; fitted "
            "automatically from --pgen otherwise)"
        ),
    )
    args = parser.parse_args(argv)

    if args.sites_vcf:
        if args.n_samples is None:
            parser.error("--sites-vcf requires --n-samples (not derivable from a sites-only VCF)")
        if args.phased_rate is None:
            parser.error(
                "--sites-vcf requires --phased-rate "
                "(no genotypes to count phased/unphased calls from)"
            )
        if not (0.0 <= args.phased_rate <= 1.0):
            parser.error(f"--phased-rate must be in [0, 1], got {args.phased_rate}")
        if args.phase_sample_mb is not None:
            parser.error(
                "--phase-sample-mb is only valid with --pgen (it only affects "
                "the phased_rate window sampled from genotypes; --sites-vcf's "
                "phased_rate comes from --phased-rate instead)"
            )
    else:
        if args.n_samples is not None:
            parser.error(
                "--n-samples is only valid with --sites-vcf (--pgen derives it from .psam)"
            )
        if args.phased_rate is not None:
            parser.error(
                "--phased-rate is only valid with --sites-vcf (--pgen fits it from genotypes)"
            )
        if args.phase_sample_mb is None:
            args.phase_sample_mb = 1.0

    profile = _fit_from_sites_vcf(args) if args.sites_vcf else _fit_from_pgen(args)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(profile, indent=2) + "\n")
    print(f"wrote {args.out}")
    _validate_with_rust(args.out)


if __name__ == "__main__":
    main()
