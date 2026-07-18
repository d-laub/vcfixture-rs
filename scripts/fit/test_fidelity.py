"""Generate from a profile, re-fit the output, assert the stats round-trip.

This is the real test that the samplers reproduce the profile they were given.
Bridges BCF -> pgen with plink2 so re-fitting uses the same code path as the
original fit.

Profile choice: `germline-1kgp` is fitted from the 1kGP *phased* panel, and
phasing is precisely what removes unphaseable singletons -- its singleton
fraction is 0.0, which would make the SFS assertion below pass vacuously.
`germline-1kgp-unphased` is fitted from the same cohort's raw unphased
callset and has a real, sharp singleton fraction (0.3579), so it is the
profile that actually exercises the SFS sampler.

Sample count: generating at anything less than the profile's native
`n_samples_source` (3202) would put fewer than `2 * samples` alleles in
existence, forcing every SFS bin above that ceiling to be clamped/rescaled --
which would distort the very singleton fraction under test. Generating at
3202 removes that confound.
"""

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PROFILE = REPO / "profiles" / "germline-1kgp-unphased.json"

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
    original_full = json.loads(PROFILE.read_text())
    original = original_full["fitted"]
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        bcf = tmp / "gen.bcf"
        _generate(bcf, samples=3202, per_contig=2000)
        refit_full = _refit(bcf, tmp)
        refit = refit_full["fitted"]

        # Ti/Tv is a single scalar and the most direct sampler check.
        assert abs(refit["titv"] - original["titv"]) < 0.15, \
            f"titv {refit['titv']} != {original['titv']}"

        # Class mix must survive the round-trip.
        for cls in ("snp", "insertion", "deletion"):
            a = original["variant_classes"][cls]
            b = refit["variant_classes"][cls]
            assert abs(a - b) < 0.05, f"class {cls}: {b} != {a}"

        # The singleton fraction is the stat the whole empirical-SFS decision
        # exists to preserve, so it gets the tightest guard. This profile's
        # singleton fraction (0.3579) is a real, sharp, easy-to-break value --
        # unlike the phased germline-1kgp profile, where it is 0.0 and this
        # assertion would pass vacuously.
        a = original["sfs"]["weights"][0] / sum(original["sfs"]["weights"])
        b = refit["sfs"]["weights"][0] / sum(refit["sfs"]["weights"])
        assert abs(a - b) < 0.06, f"singleton fraction {b} != {a}"

        assert refit_full["dialed"]["ploidy"] == original_full["dialed"]["ploidy"]
