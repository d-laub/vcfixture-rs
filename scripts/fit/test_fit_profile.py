import json
import re
import shutil
import subprocess
import warnings
from pathlib import Path

import polars as pl
import pytest

from fit_profile import (
    CLASS_NAMES,
    INDEL_EDGES,
    build_profile,
    classify,
    _classify_df,
    _explode_alleles,
    _gap_bins_lazy,
    _gap_edges,
    _multiallelic_rate_from_row,
    _multiallelic_rate_lazy,
    _payload_choices,
    _sfs_edges,
    _titv_lazy,
    assert_pvar_sorted,
    class_mix_from_counts,
    compute_pvar_stats,
    fit_missing_rate_from_sites_vcf,
    fit_sfs,
    fit_sfs_from_sites_vcf,
    histogram,
    histogram_lazy,
    main,
    read_pvar,
    read_sites_vcf,
)

PLINK2_AVAILABLE = shutil.which("plink2") is not None
BCFTOOLS_AVAILABLE = shutil.which("bcftools") is not None
CARGO_AVAILABLE = shutil.which("cargo") is not None

_PROFILE_RS = Path(__file__).resolve().parents[2] / "src" / "bulk" / "profile.rs"


def _rust_classmix_fields() -> list[str]:
    src = _PROFILE_RS.read_text()
    block = re.search(r"pub struct ClassMix \{(.*?)\}", src, re.S).group(1)
    return re.findall(r"pub (\w+): f64", block)


def _rust_payload_variants() -> list[str]:
    src = _PROFILE_RS.read_text()
    block = re.search(r"pub enum Payload \{(.*?)\}", src, re.S).group(1)
    camel = re.findall(r"\b([A-Z]\w+),", block)
    return [re.sub(r"(?<!^)(?=[A-Z])", "-", c).lower() for c in camel]


def test_class_names_match_rust_classmix():
    assert list(CLASS_NAMES) == _rust_classmix_fields()


def test_payload_choices_match_rust_enum():
    assert sorted(_payload_choices()) == sorted(_rust_payload_variants())


def test_histogram_weights_are_one_shorter_than_edges():
    h = histogram([1, 1, 2, 5, 50], edges=[1, 2, 10, 100])
    assert len(h["weights"]) == len(h["edges"]) - 1
    assert abs(sum(h["weights"]) - 1.0) < 1e-9


def test_class_mix_sums_to_one():
    m = class_mix_from_counts({"snp": 83, "insertion": 6, "deletion": 9,
                               "mnp": 1, "complex": 1, "symbolic": 0})
    assert abs(sum(m.values()) - 1.0) < 1e-9


def test_build_profile_emits_schema_valid_json():
    # gap_dist/sfs/indel_length are already-finalized histogram dicts at
    # this point -- build_profile no longer computes histograms itself
    # (that would require materializing the full gaps/acs/indel_lens
    # sequences, which is exactly what OOMs at real cohort scale). They're
    # built here with the small-scale `histogram()` helper purely because
    # this test's inputs are tiny literals.
    p = build_profile(
        name="test",
        source="/dev/null",
        n_samples=10,
        contigs=[{"id": "chr1", "n_variants": 100, "density_per_kb": 40.0}],
        gap_dist=histogram([1, 2, 3, 40], _gap_edges()),
        sfs=histogram([1, 1, 2, 19], _sfs_edges(10)),
        indel_length=histogram([1, 2, 3], INDEL_EDGES),
        class_counts={"snp": 83, "insertion": 6, "deletion": 9,
                      "mnp": 1, "complex": 1, "symbolic": 0},
        titv=2.05,
        multiallelic_rate=0.0,
        missing_rate=0.0,
        phased_rate=1.0,
        ploidy=2,
        supplied=["ploidy"],
    )
    j = json.loads(json.dumps(p))
    assert set(j) == {"name", "provenance", "fitted", "dialed"}
    assert set(j["fitted"]) == {
        "contigs", "gap_dist", "sfs", "variant_classes", "indel_length",
        "titv", "multiallelic_rate", "missing_rate", "phased_rate",
    }
    assert j["dialed"]["payload"] in {"gt-only", "gt-vaf", "gatk", "mutect2"}
    assert j["dialed"]["ploidy"] == 2
    # provenance must be populated, never left as a placeholder
    assert j["provenance"]["n_samples_source"] == 10


def test_build_profile_records_supplied_fields():
    prof = build_profile(
        name="t", source="x", n_samples=10,
        contigs=[{"id": "chr1", "n_variants": 100, "density_per_kb": 40.0}],
        gap_dist={"edges": [1.0, 2.0], "weights": [1.0]},
        sfs={"edges": [1.0, 2.0], "weights": [1.0]},
        indel_length={"edges": [1.0, 2.0], "weights": [1.0]},
        class_counts={n: 1 for n in CLASS_NAMES},
        titv=2.0, multiallelic_rate=0.1, missing_rate=0.0,
        phased_rate=1.0, ploidy=2, supplied=["ploidy", "phased_rate"],
    )
    assert prof["provenance"]["supplied"] == ["phased_rate", "ploidy"]
    assert "ploidy" not in prof["fitted"]
    assert prof["dialed"]["ploidy"] == 2


# --------------------------------------------------------------------------
# histogram(): out-of-range values must warn, not vanish silently
# --------------------------------------------------------------------------


def test_histogram_warns_when_some_values_are_out_of_range():
    # Regression test for the "Important" review finding: histogram() only
    # used to raise if ALL values fell outside edges. A partially
    # out-of-range input (e.g. gap_dist values beyond the 1e5 bp cap) must
    # at least emit a diagnostic instead of silently vanishing from the
    # normalization.
    with pytest.warns(UserWarning, match=r"dropped 1/6.*16\.7%"):
        h = histogram([1, 1, 2, 5, 50, 99999], edges=[1, 2, 10, 100])
    assert h["weights"] == pytest.approx([0.4, 0.4, 0.2])
    assert sum(h["weights"]) == pytest.approx(1.0)


def test_histogram_raises_when_all_values_are_out_of_range():
    with pytest.raises(ValueError, match="no values fell within"):
        histogram([1000, 2000], edges=[1, 2, 10, 100])


def test_histogram_does_not_warn_when_all_values_in_range():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        histogram([1, 1, 2, 5], edges=[1, 2, 10, 100])


# --------------------------------------------------------------------------
# Multiallelic ALT: per-allele split (the "Critical" finding)
# --------------------------------------------------------------------------


def test_fit_sfs_regression_comma_joined_alt_cts_no_longer_crashes(monkeypatch, tmp_path):
    # Simulate what fit_sfs() sees after plink2 --freq counts on a
    # multiallelic pgen, without needing plink2 itself: a comma-joined
    # ALT_CTS column aligned with a comma-joined ALT column.
    acount = tmp_path / "sfs.acount"
    acount.write_text(
        "#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n"
        "1\trs1\tA\tG,T\t2,1\t6\n"
        "1\trs2\tC\tT\t3\t6\n"
    )

    def fake_run_plink2(args):
        return None

    monkeypatch.setattr("fit_profile._run_plink2", fake_run_plink2)
    monkeypatch.setattr("tempfile.TemporaryDirectory", lambda: _FixedTmpDir(tmp_path))
    # n_samples=3 -> _sfs_edges(3) == [1.0, 2.0, 4.0, 6.0]. The 3 ALT_CTS
    # observations after splitting "2,1" -> 2.0, 1.0 are [2.0, 1.0, 3.0]:
    # 1.0 lands in [1,2) (1 obs), 2.0 and 3.0 both land in [2,4) (2 obs).
    # fit_sfs used to crash trying to float() the raw "2,1" string directly;
    # now it returns the finalized histogram, never a materialized list.
    sfs = fit_sfs("unused-prefix", n_samples=3)
    assert sfs["edges"] == [1.0, 2.0, 4.0, 6.0]
    assert sfs["weights"] == pytest.approx([1 / 3, 2 / 3, 0.0])


class _FixedTmpDir:
    def __init__(self, path):
        self._path = str(path)

    def __enter__(self):
        return self._path

    def __exit__(self, *exc):
        return False


def test_multiallelic_alt_corrupts_scalar_classify_but_not_the_split_pipeline():
    # Documents the silent-corruption half of the crash bug from the code
    # review: classify()/_classify_expr are only defined for a single
    # REF/ALT pair. Feeding them a raw, un-split multiallelic ALT silently
    # misclassifies -- exactly the reviewer's repro.
    assert classify("A", "A,G") == "insertion"  # wrong: two SNP alleles misread as an insertion

    df = pl.DataFrame({
        "CHROM": ["1", "1"],
        "POS": [100, 200],
        "ID": ["rs1", "rs2"],
        "REF": ["A", "C"],
        "ALT": ["G,T", "T"],
    })
    alleles = _classify_df(_explode_alleles(df.lazy())).collect()
    # 3 alleles total: rs1/G, rs1/T, rs2/T -- all real SNPs once split.
    assert alleles["class"].to_list() == ["snp", "snp", "snp"]


def test_explode_alleles_splits_multiallelic_and_is_a_no_op_for_biallelic():
    df = pl.DataFrame({
        "CHROM": ["1", "1"],
        "POS": [100, 200],
        "ID": ["rs1", "rs2"],
        "REF": ["A", "AT"],
        "ALT": ["G,T", "A"],
    })
    alleles = _explode_alleles(df.lazy()).collect()
    assert alleles.height == 3
    assert alleles["ALT"].to_list() == ["G", "T", "A"]
    # shared per-record fields are broadcast onto every allele row
    assert alleles.filter(pl.col("ID") == "rs1")["REF"].to_list() == ["A", "A"]


def test_multiallelic_rate_counts_records_not_alleles():
    # A triallelic site (2 ALT alleles) must contribute exactly 1 to both
    # the numerator and denominator of multiallelic_rate, not 2 -- this
    # must be computed on the un-exploded, per-record frame.
    df = pl.DataFrame({
        "CHROM": ["1", "1", "1"],
        "POS": [100, 200, 300],
        "ID": ["rs1", "rs2", "rs3"],
        "REF": ["A", "C", "G"],
        "ALT": ["G,T,C", "T", "A"],
    })
    n, n_multi = _multiallelic_rate_lazy(df.lazy()).collect().row(0)
    assert _multiallelic_rate_from_row(n, n_multi) == pytest.approx(1 / 3)


# --------------------------------------------------------------------------
# _titv_lazy: direct (REF,ALT) comparisons must match the
# concat_str().is_in(TRANSITION_PAIRS) reference behavior.
# --------------------------------------------------------------------------


def test_titv_direct_matches_is_in_reference():
    alleles = pl.LazyFrame({
        "class": ["snp","snp","snp","snp","insertion"],
        "REF":   ["A","G","C","A","A"],
        "ALT":   ["G","A","T","C","AT"],  # A>G, G>A, C>T ts; A>C tv
    })
    got = _titv_lazy(alleles).collect().row(0)  # (n_snps, n_ts)
    assert got == (4, 3)


# --------------------------------------------------------------------------
# _gap_bins_lazy: shift+mask must match sort+window on sorted input, and
# assert_pvar_sorted must guard the precondition that makes them equal.
# --------------------------------------------------------------------------


def test_gap_bins_matches_sorted_window_reference():
    lf = pl.LazyFrame({
        "CHROM": ["1","1","1","2","2"],
        "POS":   [100, 250, 251, 5, 30],
        "ID": ["."]*5, "REF": ["A"]*5, "ALT": ["T"]*5,
    })
    got = _gap_bins_lazy(lf).collect()
    ref_gaps = (lf.sort(["CHROM","POS"])
        .select(pl.col("POS").diff().over("CHROM").alias("gap"))
        .filter(pl.col("gap").is_not_null() & (pl.col("gap") > 0)))
    ref = histogram_lazy(ref_gaps, pl.col("gap"), _gap_edges()).collect()
    # group_by("_bin") has no defined row order (confirmed non-deterministic
    # across repeated runs even between two logically-identical plans), so
    # sort by the bin index before comparing -- the invariant under test is
    # bit-identical bin *counts*, not group_by iteration order.
    assert got.sort("_bin").equals(ref.sort("_bin"))


def test_assert_pvar_sorted_rejects_unsorted():
    lf = pl.LazyFrame({"CHROM": ["1","1"], "POS": [200, 100],
                       "ID": [".","."], "REF": ["A","A"], "ALT": ["T","T"]})
    with pytest.raises(ValueError, match="not sorted"):
        assert_pvar_sorted(lf)


def test_read_pvar_skips_meta_and_reads_all_rows(tmp_path):
    p = tmp_path / "t.pvar"
    p.write_text(
        "##fileformat=PVARv1.0\n"
        "##contig=<ID=1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\t.\tA\tG\n"
        "1\t200\t.\tC\tT\n"
        "2\t50\t.\tG\tA\n"
    )
    lf = read_pvar(p)
    df = lf.collect()
    assert df.height == 3
    assert df["CHROM"].to_list() == ["1", "1", "2"]
    assert df["POS"].to_list() == [100, 200, 50]


def test_compute_pvar_stats_gives_two_class_and_two_indel_observations_per_multiallelic_site():
    # The brief's stated invariant: a multiallelic site with 2 ALTs
    # contributes 2 class observations and (if applicable) 2 indel-length
    # observations, while still counting as ONE record for
    # multiallelic_rate. A second, plain SNP record is included purely so
    # `_titv` (which compute_pvar_stats also computes) has a transition and
    # a transversion to work with -- unrelated to what this test checks.
    df = pl.DataFrame({
        "CHROM": ["1", "1"],
        "POS": [100, 200],
        "ID": ["rs1", "rs2"],
        "REF": ["A", "A"],
        "ALT": ["AG,ATT", "C"],  # rs1: two insertions of length 1 and 2; rs2: A>C transversion
    })
    stats = compute_pvar_stats(df.lazy())
    assert stats["class_counts"]["insertion"] == 2
    assert stats["class_counts"]["snp"] == 1
    # indel_length is now a finalized histogram (never a materialized list
    # of lengths): lengths 1 and 2 land in INDEL_EDGES bins [1,2) and [2,3)
    # respectively, one observation each.
    indel_hist = stats["indel_length"]
    assert indel_hist["edges"] == INDEL_EDGES
    assert indel_hist["weights"][0] == pytest.approx(0.5)  # length 1 -> [1,2)
    assert indel_hist["weights"][1] == pytest.approx(0.5)  # length 2 -> [2,3)
    assert sum(indel_hist["weights"]) == pytest.approx(1.0)
    # 1 multiallelic record (rs1) out of 2 total records.
    assert stats["multiallelic_rate"] == pytest.approx(0.5)
    # Regression guard for the six-plan collection (D4): both records are on
    # contig "1", so exactly one contig; rs2 (A>C) is the only SNP and it's
    # a transversion, so titv = n_ts / n_tv = 0 / 1 = 0.0.
    assert len(stats["contigs"]) == 1
    assert stats["titv"] == pytest.approx(0.0)
    # One inter-record gap (POS 100 -> 200, both on contig "1"), landing in
    # a single gap_dist bin -- weights still sum to 1.0.
    assert sum(stats["gap_dist"]["weights"]) == pytest.approx(1.0)


# --------------------------------------------------------------------------
# SFS singleton bin invariant: must stay exactly [1, 2)
# --------------------------------------------------------------------------


@pytest.mark.parametrize("n", [1, 2, 3, 10, 3202, 16007])
def test_sfs_edges_first_bin_is_exactly_one_two(n):
    edges = _sfs_edges(n)
    assert edges[0] == 1.0
    assert edges[1] == 2.0


# --------------------------------------------------------------------------
# End-to-end: reviewer's gold-standard repro (skips gracefully without plink2)
# --------------------------------------------------------------------------


def _write_synthetic_vcf(path):
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
        "1\t100\trs1\tA\tG,T\t.\t.\t.\tGT\t1/2\t0/0\t0|1\n"
        "1\t200\trs2\tC\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n"
        "1\t300\trs3\tAT\tA\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n"
    )


def _replica_validate_profile(p: dict) -> None:
    """Replicates the invariants `Profile::validate` enforces in Rust
    (src/bulk/profile.rs), without depending on the Rust crate itself:
    strictly-increasing finite histogram edges, non-negative finite
    weights summing > 0, a ClassMix summing to 1.0 within 1e-6, and
    every rate in [0, 1].
    """
    for hist_name in ("gap_dist", "sfs", "indel_length"):
        h = p["fitted"][hist_name]
        assert len(h["weights"]) == len(h["edges"]) - 1
        assert all(b > a for a, b in zip(h["edges"], h["edges"][1:]))
        assert all(w >= 0 for w in h["weights"])
        assert sum(h["weights"]) > 0
    classmix = p["fitted"]["variant_classes"]
    assert abs(sum(classmix.values()) - 1.0) < 1e-6
    for rate in ("multiallelic_rate", "missing_rate", "phased_rate"):
        assert 0.0 <= p["fitted"][rate] <= 1.0
    assert p["dialed"]["ploidy"] >= 1
    assert len(p["fitted"]["contigs"]) >= 1


@pytest.mark.skipif(not PLINK2_AVAILABLE, reason="plink2 not installed")
def test_end_to_end_multiallelic_pgen_does_not_crash_and_has_per_allele_semantics(tmp_path):
    # This is the reviewer's exact repro: a pvar with a genuine
    # multiallelic record (REF=A ALT=G,T), converted with
    # `plink2 --make-pgen` (which retains it natively, unlike bcftools
    # norm -m-), run all the way through the extraction script.
    vcf_path = tmp_path / "test.vcf"
    _write_synthetic_vcf(vcf_path)
    prefix = str(tmp_path / "prefix")
    subprocess.run(
        ["plink2", "--vcf", str(vcf_path), "--make-pgen", "--out", prefix],
        check=True,
        capture_output=True,
    )

    out_path = tmp_path / "profile.json"
    # main() previously raised ValueError: could not convert string to
    # float: '2,1' inside fit_sfs() -> histogram() on this exact input.
    main(["--pgen", prefix, "--name", "synthtest", "--out", str(out_path)])

    p = json.loads(out_path.read_text())
    _replica_validate_profile(p)

    # Per-allele semantics: rs1 contributes 2 SNP observations (A/G, A/T),
    # rs2 contributes 1 SNP observation (C/T), rs3 contributes 1 deletion.
    assert p["fitted"]["variant_classes"]["snp"] == pytest.approx(0.75)
    assert p["fitted"]["variant_classes"]["deletion"] == pytest.approx(0.25)
    # multiallelic_rate is per-RECORD: 1 multiallelic site out of 3 records.
    assert p["fitted"]["multiallelic_rate"] == pytest.approx(1 / 3)
    # sfs singleton bin edges are untouched by any of this.
    assert p["fitted"]["sfs"]["edges"][0] == 1.0
    assert p["fitted"]["sfs"]["edges"][1] == 2.0


@pytest.mark.skipif(not PLINK2_AVAILABLE, reason="plink2 not installed")
def test_fit_sfs_preserves_positional_alignment_between_alt_and_alt_cts(tmp_path):
    # A stronger check than symmetric counts: rs1's two ALT alleles have
    # *different* counts (G=2, T=1), so a bug that mixed up which count
    # belongs to which allele -- or just summed/averaged them -- would be
    # caught here, not just a crash.
    vcf_path = tmp_path / "test.vcf"
    _write_synthetic_vcf(vcf_path)
    prefix = str(tmp_path / "prefix")
    subprocess.run(
        ["plink2", "--vcf", str(vcf_path), "--make-pgen", "--out", prefix],
        check=True,
        capture_output=True,
    )
    # n_samples=3 -> _sfs_edges(3) == [1.0, 2.0, 4.0, 6.0]. The 4 per-allele
    # observations are rs1/G=2.0, rs1/T=1.0, rs2=3.0, rs3=3.0: if G and T's
    # counts were summed/averaged instead of kept independent, or
    # misaligned with the wrong allele, this bin split would come out wrong.
    sfs = fit_sfs(prefix, n_samples=3)
    assert sfs["edges"] == [1.0, 2.0, 4.0, 6.0]
    assert sfs["weights"] == pytest.approx([0.25, 0.75, 0.0])


# --------------------------------------------------------------------------
# read_sites_vcf(): the sites-only-VCF input path (task 9b)
# --------------------------------------------------------------------------


def _write_sites_vcf(tmp_path):
    """Small bgzipped+tabix-indexed sites-only VCF fixture with a PASS/VQSR mix.

    Two contigs (numeric-looking ids, to exercise the CHROM-stays-string
    contract) each carry one PASS record and one non-PASS record that must
    be dropped -- one VQSR-tranche, one LowQual, matching the two FILTER
    values actually seen on the real 1kGP raw callset.
    """
    vcf_path = tmp_path / "sites.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=21>\n"
        "##contig=<ID=22>\n"
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100\n"
        "21\t200\t.\tC\tT\t.\tVQSRTrancheSNP99.80to100.00\tAC=5;AN=100\n"
        "22\t150\t.\tG\tA\t.\tPASS\tAC=50;AN=100\n"
        "22\t250\t.\tT\tC\t.\tLowQual\tAC=2;AN=100\n"
    )
    gz_path = tmp_path / "sites.vcf.gz"
    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", str(gz_path), str(vcf_path)],
        check=True,
        capture_output=True,
    )
    subprocess.run(["bcftools", "index", "-t", str(gz_path)], check=True, capture_output=True)
    return gz_path


@pytest.mark.skipif(not BCFTOOLS_AVAILABLE, reason="bcftools not installed")
def test_read_sites_vcf_drops_non_pass_records(tmp_path):
    vcf_path = _write_sites_vcf(tmp_path)
    df = read_sites_vcf(vcf_path).collect()
    # Only the two PASS records (21:100, 22:150) survive; the
    # VQSRTranche/LowQual records must not reach the profile.
    assert sorted(df["POS"].to_list()) == [100, 150]
    assert set(df["AC"].to_list()) == {1, 50}


@pytest.mark.skipif(not BCFTOOLS_AVAILABLE, reason="bcftools not installed")
def test_read_sites_vcf_keeps_numeric_looking_contig_ids_as_strings(tmp_path):
    vcf_path = _write_sites_vcf(tmp_path)
    df = read_sites_vcf(vcf_path).collect()
    assert df["CHROM"].dtype == pl.Utf8
    assert set(df["CHROM"].to_list()) == {"21", "22"}


@pytest.mark.skipif(not BCFTOOLS_AVAILABLE, reason="bcftools not installed")
def test_read_sites_vcf_returns_a_lazyframe(tmp_path):
    # Guards the memory contract documented in read_sites_vcf's docstring:
    # callers must be able to reduce this before collecting, at real
    # (tens-of-millions-of-rows) scale.
    vcf_path = _write_sites_vcf(tmp_path)
    lf = read_sites_vcf(vcf_path)
    assert isinstance(lf, pl.LazyFrame)


# --------------------------------------------------------------------------
# fit_missing_rate_from_sites_vcf(): AN -> missing_rate arithmetic
# --------------------------------------------------------------------------


def test_fit_missing_rate_from_sites_vcf_full_an_gives_zero():
    # AN == 2 * n_samples everywhere -> no missingness.
    lf = pl.DataFrame({"AN": [20, 20, 20]}).lazy()
    assert fit_missing_rate_from_sites_vcf(lf, n_samples=10) == pytest.approx(0.0)


def test_fit_missing_rate_from_sites_vcf_half_an_gives_half():
    # AN == n_samples everywhere (half of 2 * n_samples) -> 50% missing.
    lf = pl.DataFrame({"AN": [10, 10, 10]}).lazy()
    assert fit_missing_rate_from_sites_vcf(lf, n_samples=10) == pytest.approx(0.5)


def test_fit_sfs_from_sites_vcf_uses_the_shared_sfs_edges():
    # n_samples=3 -> _sfs_edges(3) == [1.0, 2.0, 4.0, 6.0]. AC values
    # [1, 1, 2, 3]: the two 1s land in [1, 2), the 2 and 3 land in [2, 4).
    lf = pl.DataFrame({"AC": [1, 1, 2, 3]}).lazy()
    sfs = fit_sfs_from_sites_vcf(lf, n_samples=3)
    assert sfs["edges"] == _sfs_edges(3)
    assert sfs["weights"] == pytest.approx([0.5, 0.5, 0.0])


# --------------------------------------------------------------------------
# CLI validation: --pgen / --sites-vcf are mutually exclusive, and
# --n-samples / --phased-rate are required with (and only with) --sites-vcf.
# --------------------------------------------------------------------------


def test_main_rejects_pgen_and_sites_vcf_together(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "--pgen", "unused-prefix",
                "--sites-vcf", "unused.vcf.gz",
                "--name", "x",
                "--out", str(tmp_path / "out.json"),
            ]
        )


def test_main_rejects_sites_vcf_without_n_samples(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "--sites-vcf", "unused.vcf.gz",
                "--phased-rate", "0.5",
                "--name", "x",
                "--out", str(tmp_path / "out.json"),
            ]
        )


def test_main_rejects_sites_vcf_without_phased_rate(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "--sites-vcf", "unused.vcf.gz",
                "--n-samples", "100",
                "--name", "x",
                "--out", str(tmp_path / "out.json"),
            ]
        )


def test_main_rejects_out_of_range_phased_rate_with_sites_vcf(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "--sites-vcf", "unused.vcf.gz",
                "--n-samples", "100",
                "--phased-rate", "1.5",
                "--name", "x",
                "--out", str(tmp_path / "out.json"),
            ]
        )


def test_main_rejects_n_samples_with_pgen(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "--pgen", "unused-prefix",
                "--n-samples", "100",
                "--name", "x",
                "--out", str(tmp_path / "out.json"),
            ]
        )


def test_main_rejects_phased_rate_with_pgen(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "--pgen", "unused-prefix",
                "--phased-rate", "0.5",
                "--name", "x",
                "--out", str(tmp_path / "out.json"),
            ]
        )


def test_main_rejects_phase_sample_mb_with_sites_vcf(tmp_path):
    # --phase-sample-mb only affects the pgen path's phased_rate window
    # sampling; under --sites-vcf, phased_rate comes verbatim from
    # --phased-rate, so a hand-supplied --phase-sample-mb would be silently
    # meaningless -- reject it like --n-samples/--phased-rate are rejected
    # under --pgen.
    with pytest.raises(SystemExit):
        main(
            [
                "--sites-vcf", "unused.vcf.gz",
                "--n-samples", "100",
                "--phased-rate", "0.5",
                "--phase-sample-mb", "2.0",
                "--name", "x",
                "--out", str(tmp_path / "out.json"),
            ]
        )


# --------------------------------------------------------------------------
# End-to-end: --sites-vcf through main() produces a schema-valid profile
# --------------------------------------------------------------------------


def _write_sites_vcf_for_e2e(tmp_path):
    """A single-contig fixture shaped so `compute_pvar_stats` can run to
    completion on the PASS-only records: three well-separated PASS
    positions (nonzero gaps for `gap_dist`) with one transition and one
    transversion (both a Ti and a Tv observation for `titv`) plus one
    deletion (a non-empty `indel_length` histogram), and one VQSR-tranche
    and one LowQual record that must be dropped.
    """
    vcf_path = tmp_path / "sites_e2e.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=21>\n"
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100\n"          # SNP transition
        "21\t150\t.\tC\tT\t.\tVQSRTrancheSNP99.80to100.00\tAC=5;AN=100\n"
        "21\t200\t.\tG\tT\t.\tPASS\tAC=50;AN=100\n"         # SNP transversion
        "21\t250\t.\tT\tC\t.\tLowQual\tAC=2;AN=100\n"
        "21\t300\t.\tAT\tA\t.\tPASS\tAC=3;AN=100\n"         # deletion, length 1
    )
    gz_path = tmp_path / "sites_e2e.vcf.gz"
    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", str(gz_path), str(vcf_path)],
        check=True,
        capture_output=True,
    )
    subprocess.run(["bcftools", "index", "-t", str(gz_path)], check=True, capture_output=True)
    return gz_path


@pytest.mark.skipif(not BCFTOOLS_AVAILABLE, reason="bcftools not installed")
def test_end_to_end_sites_vcf_produces_a_valid_profile(tmp_path):
    vcf_path = _write_sites_vcf_for_e2e(tmp_path)
    out_path = tmp_path / "profile.json"
    main(
        [
            "--sites-vcf", str(vcf_path),
            "--name", "sitestest",
            "--out", str(out_path),
            "--n-samples", "50",
            "--phased-rate", "0.9",
        ]
    )

    p = json.loads(out_path.read_text())
    _replica_validate_profile(p)
    assert p["provenance"]["source"] == str(vcf_path)
    assert p["provenance"]["n_samples_source"] == 50
    # 3 PASS records survive (POS 100, 200, 300) out of 5 raw records.
    assert p["provenance"]["n_variants_source"] == 3
    # phased_rate is taken verbatim from --phased-rate: no genotypes exist
    # to fit it from.
    assert p["fitted"]["phased_rate"] == pytest.approx(0.9)
    # AN == 100 == 2 * n_samples everywhere -> no missingness.
    assert p["fitted"]["missing_rate"] == pytest.approx(0.0)


# --------------------------------------------------------------------------
# validate-profile binary: CI gate for freshly-written profiles
# --------------------------------------------------------------------------


@pytest.mark.skipif(not CARGO_AVAILABLE, reason="cargo not installed")
def test_validate_profile_binary_rejects_nan(tmp_path):
    prof = build_profile(
        name="bad", source="x", n_samples=10,
        contigs=[{"id": "chr1", "n_variants": 100, "density_per_kb": 40.0}],
        gap_dist={"edges": [1.0, 2.0], "weights": [1.0]},
        sfs={"edges": [1.0, 2.0], "weights": [1.0]},
        indel_length={"edges": [1.0, 2.0], "weights": [1.0]},
        class_counts={n: 1 for n in CLASS_NAMES},
        titv=2.0, multiallelic_rate=0.1, missing_rate=0.0,
        phased_rate=1.0, ploidy=2, supplied=["ploidy"],
    )
    prof["fitted"]["variant_classes"]["snp"] = float("nan")  # poison
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(prof))
    r = subprocess.run(
        ["cargo", "run", "--quiet", "--features", "bulk",
         "--bin", "validate-profile", "--", str(p)],
        capture_output=True, text=True,
    )
    assert r.returncode != 0, r.stdout + r.stderr
