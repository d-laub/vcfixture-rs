#![cfg(feature = "bulk")]

use std::num::NonZero;
use vcfixture::bulk::{BulkSpec, Payload, Profile, Size};

fn spec() -> BulkSpec {
    BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
        .samples(8)
        .contigs(["chr1", "chr2", "chr3"])
        .payload(Payload::GtOnly)
        .seed(42)
        .workers(NonZero::new(2).unwrap())
}

#[test]
fn records_per_contig_is_exact() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = spec()
        .size(Size::RecordsPerContig(100))
        .write(&path)
        .unwrap();
    assert_eq!(s.n_records_total(), 300);
    assert_eq!(s.per_contig["chr1"].n_records, 100);
    assert_eq!(s.per_contig.len(), 3);
}

#[test]
fn same_seed_gives_byte_identical_output_across_thread_counts() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.bcf");
    let b = dir.path().join("b.bcf");
    spec()
        .size(Size::RecordsPerContig(50))
        .workers(NonZero::new(1).unwrap())
        .write(&a)
        .unwrap();
    spec()
        .size(Size::RecordsPerContig(50))
        .workers(NonZero::new(4).unwrap())
        .write(&b)
        .unwrap();
    assert_eq!(
        std::fs::read(&a).unwrap(),
        std::fs::read(&b).unwrap(),
        "output must not depend on thread count"
    );
}

#[test]
fn different_seeds_differ() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.bcf");
    let b = dir.path().join("b.bcf");
    spec()
        .seed(1)
        .size(Size::RecordsPerContig(50))
        .write(&a)
        .unwrap();
    spec()
        .seed(2)
        .size(Size::RecordsPerContig(50))
        .write(&b)
        .unwrap();
    assert_ne!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap());
}

#[test]
fn declared_contig_length_equals_populated_span() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = spec()
        .size(Size::RecordsPerContig(200))
        .write(&path)
        .unwrap();

    let mut r = noodles_bcf::io::reader::Builder::default()
        .build_from_path(&path)
        .unwrap();
    let header = r.read_header().unwrap();
    for (id, contig) in header.contigs() {
        let declared = contig.length().expect("contig must declare a length") as u64;
        let span = s.per_contig[id.as_str()].pos_max;
        assert_eq!(
            declared, span,
            "contig {id} declared length {declared} must equal populated span {span}"
        );
        // and it must not be a real hg38 length
        assert!(
            declared < 248_956_422,
            "contig {id} must not use a real hg38 length"
        );
    }
}

#[test]
fn target_size_lands_near_the_target() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let target = 512 * 1024;
    spec().size(Size::Target(target)).write(&path).unwrap();
    let got = std::fs::metadata(&path).unwrap().len();
    assert!(got >= target, "got {got} < target {target}");
    assert!(
        got < target + 256 * 1024,
        "overshoot too large: {got} vs {target}"
    );
}

#[test]
fn summary_matches_an_independent_read_back() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = spec()
        .size(Size::RecordsPerContig(100))
        .write(&path)
        .unwrap();

    let mut r = noodles_bcf::io::reader::Builder::default()
        .build_from_path(&path)
        .unwrap();
    let _ = r.read_header().unwrap();
    let n = r.records().count() as u64;
    assert_eq!(
        n,
        s.n_records_total(),
        "summary must match what a reader sees"
    );
}

#[test]
fn index_is_written_and_positions_are_sorted() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    spec()
        .size(Size::RecordsPerContig(100))
        .write(&path)
        .unwrap();
    assert!(path.with_extension("bcf.csi").exists());
}

#[test]
fn payload_presets_all_write_readable_files() {
    for payload in [
        Payload::GtOnly,
        Payload::GtVaf,
        Payload::Gatk,
        Payload::Mutect2,
    ] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.bcf");
        let s = spec()
            .payload(payload.clone())
            .size(Size::RecordsPerContig(20))
            .write(&path)
            .unwrap();
        assert_eq!(s.n_records_total(), 60, "payload {payload:?}");
        let mut r = noodles_bcf::io::reader::Builder::default()
            .build_from_path(&path)
            .unwrap();
        let _ = r.read_header().unwrap();
        assert_eq!(r.records().count(), 60, "payload {payload:?}");
    }
}

/// `Size::Records(total)` is not per-contig-exact like `RecordsPerContig` —
/// it splits `total` across contigs proportional to fitted density. This is
/// not covered by the brief's 8 tests but is part of the required public
/// `Size` interface, so it gets its own coverage: the split must sum to
/// exactly `total` and must be non-uniform for a profile with non-uniform
/// per-contig density (the placeholder profile's 3 contigs are all
/// identical, so a custom profile is needed to actually exercise weighting).
#[test]
fn records_total_splits_by_density_and_sums_exactly() {
    let profile = Profile::from_json(NONUNIFORM_DENSITY_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1", "chr2", "chr3"])
        .seed(7)
        .size(Size::Records(300))
        .write(&path)
        .unwrap();

    assert_eq!(
        s.n_records_total(),
        300,
        "must split to exactly the requested total"
    );
    // chr1's density (80) is 8x chr3's (10) in the fitted profile, so it
    // must get noticeably more records, not an even 100/100/100 split.
    assert!(
        s.per_contig["chr1"].n_records > s.per_contig["chr3"].n_records * 2,
        "higher-density contig must get proportionally more records: {:?}",
        s.per_contig
    );
}

/// The real-data finding this task was built around: committed profiles are
/// fit from pvar files with *bare* contig ids (`"1"`, `"2"`, `"3"`), not
/// `chr1`-style ids, while output contig names are chosen independently by
/// the caller. Generation must not fail, error, or silently zero out a
/// contig's stats just because the requested output name doesn't textually
/// match any `fitted.contigs[].id` -- it must fall back positionally (see
/// `resolve_contig_stat`'s doc comment in `src/bulk/mod.rs`).
#[test]
fn bare_id_profile_resolves_positionally_for_chr_prefixed_output_names() {
    let profile = Profile::from_json(NONUNIFORM_DENSITY_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    // The profile's fitted contig ids are bare ("1", "2", "3"); the
    // requested output names are chr-prefixed and in the same order.
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1", "chr2", "chr3"])
        .seed(7)
        .size(Size::Records(300))
        .write(&path)
        .unwrap();

    assert_eq!(s.per_contig.len(), 3);
    assert!(s.per_contig.contains_key("chr1"));
    // Positional resolution (index 0 -> fitted contig 0, id "1", density 80)
    // must still pick up the highest-density contig's weight for "chr1".
    assert!(
        s.per_contig["chr1"].n_records > s.per_contig["chr3"].n_records * 2,
        "chr1 (positionally resolved to fitted id \"1\", density 80) must \
         outweigh chr3 (fitted id \"3\", density 10): {:?}",
        s.per_contig
    );
}

/// Same profile as the two tests above: fitted contig ids are bare
/// (`"1"`/`"2"`/`"3"`), with 8x/4x/1x relative density, so the weighted
/// split has an unambiguous, checkable direction.
const NONUNIFORM_DENSITY_PROFILE: &str = r#"
{
  "name": "nonuniform-density-test",
  "provenance": {
    "source": "test fixture",
    "n_samples_source": 0,
    "n_variants_source": 0,
    "fitted_on": "1970-01-01",
    "fit_tool_version": "0.0.0"
  },
  "fitted": {
    "contigs": [
      { "id": "1", "n_variants": 8000, "density_per_kb": 80.0 },
      { "id": "2", "n_variants": 4000, "density_per_kb": 40.0 },
      { "id": "3", "n_variants": 1000, "density_per_kb": 10.0 }
    ],
    "gap_dist": { "edges": [1.0, 10.0, 100.0, 1000.0], "weights": [0.6, 0.3, 0.1] },
    "sfs": { "edges": [1.0, 2.0, 10.0, 100.0, 6404.0], "weights": [0.476, 0.2, 0.2, 0.124] },
    "variant_classes": {
      "snp": 0.83,
      "insertion": 0.06,
      "deletion": 0.09,
      "mnp": 0.005,
      "complex": 0.005,
      "symbolic": 0.01
    },
    "indel_length": { "edges": [1.0, 2.0, 6.0, 20.0, 100.0], "weights": [0.5, 0.4, 0.085, 0.015] },
    "titv": 2.05,
    "multiallelic_rate": 0.0,
    "missing_rate": 0.0,
    "phased_rate": 1.0,
    "ploidy": 2
  },
  "dialed": { "payload": "gt-only" }
}
"#;
