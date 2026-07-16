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
