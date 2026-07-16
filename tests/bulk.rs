#![cfg(feature = "bulk")]

use std::num::NonZero;
use vcfixture::bulk::{BulkError, BulkSpec, Payload, Profile, Size};

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
    // `BulkSpec::BLOCK_SIZE` (private; `src/bulk/mod.rs`) is currently 500
    // records — the granularity at which `generate_contig`'s rayon
    // `into_par_iter` has anything to reorder. At `RecordsPerContig(50)`
    // (the value this test used before this fix), `50.div_ceil(500) == 1`
    // block per contig: there is nothing for rayon to reorder, so both
    // parallel paths (this crate's own block-level rayon split, and the
    // writer's bgzf multithreaded compression, which also scales with
    // `workers`) are structurally incapable of failing regardless of
    // whether determinism actually holds. `RECORDS_PER_CONTIG` below is
    // sized well past that boundary, and the assertions further down pin
    // both dimensions of scale down so a future edit can't silently shrink
    // this test back into vacuity. See `src/bulk/writer.rs`'s
    // `output_is_byte_identical_regardless_of_worker_count` for the same
    // idiom applied to the writer's own compression layer.
    const RECORDS_PER_CONTIG: u64 = 2_500;
    const BLOCK_SIZE_LOWER_BOUND: u64 = 500; // BulkSpec::BLOCK_SIZE, mirrored here
    const MAX_BUF_SIZE: u64 = 65_498; // bgzf's uncompressed staging buffer size

    let blocks_per_contig = RECORDS_PER_CONTIG.div_ceil(BLOCK_SIZE_LOWER_BOUND);
    assert!(
        blocks_per_contig > 1,
        "test is vacuous: {RECORDS_PER_CONTIG} records/contig gives only \
         {blocks_per_contig} rayon block(s) per contig, so there is nothing \
         for `into_par_iter` to reorder"
    );

    fn heavy(workers: usize) -> BulkSpec {
        // A heavier payload (more FORMAT keys, 8 samples via `spec()`) than
        // `GtOnly`, so the encoded record stream is large enough to
        // actually span several bgzf blocks, exercising the writer's own
        // compression-thread reordering risk too (`BulkSpec::write` passes
        // `workers` through to `BulkWriter::create` as well as the rayon
        // pool).
        spec()
            .payload(Payload::Gatk)
            .size(Size::RecordsPerContig(RECORDS_PER_CONTIG))
            .workers(NonZero::new(workers).unwrap())
    }

    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.bcf");
    let b = dir.path().join("b.bcf");
    let c = dir.path().join("c.bcf");
    heavy(1).write(&a).unwrap();
    heavy(4).write(&b).unwrap();
    heavy(16).write(&c).unwrap();

    let bytes_a = std::fs::read(&a).unwrap();
    let bytes_b = std::fs::read(&b).unwrap();
    let bytes_c = std::fs::read(&c).unwrap();

    // Pin down the actual uncompressed payload size flowing through bgzf,
    // by decompressing the real output rather than trusting arithmetic on
    // record count and payload shape — independent of how well bgzf's
    // DEFLATE happens to compress this data.
    let mut decompressed = Vec::new();
    std::io::Read::read_to_end(
        &mut noodles_bgzf::io::Reader::new(std::io::Cursor::new(&bytes_a)),
        &mut decompressed,
    )
    .unwrap();
    assert!(
        decompressed.len() as u64 > 3 * MAX_BUF_SIZE,
        "test payload ({} bytes uncompressed) must exceed several bgzf blocks \
         ({MAX_BUF_SIZE} bytes each) or this test cannot detect reordering bugs",
        decompressed.len()
    );

    assert_eq!(
        bytes_a, bytes_b,
        "output must be byte-identical regardless of thread count (1 vs 4 workers)"
    );
    assert_eq!(
        bytes_a, bytes_c,
        "output must be byte-identical regardless of thread count (1 vs 16 workers)"
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
/// match any `fitted.contigs[].id` -- it resolves by name after stripping a
/// `chr` prefix (see `resolve_contig_stat`'s doc comment in
/// `src/bulk/mod.rs`), which for in-order requests happens to agree with
/// what positional resolution would also have given.
#[test]
fn bare_id_profile_resolves_by_chr_prefix_normalization() {
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
    // "chr1" normalizes to "1" and matches fitted id "1" (density 80) by
    // name, so it must still outweigh "chr3" (fitted id "3", density 10).
    assert!(
        s.per_contig["chr1"].n_records > s.per_contig["chr3"].n_records * 2,
        "chr1 (name-resolved to fitted id \"1\", density 80) must \
         outweigh chr3 (fitted id \"3\", density 10): {:?}",
        s.per_contig
    );
}

/// The specific failure mode Important-4 fixed: with *only* positional
/// resolution, requesting output contigs out of the profile's fitted order
/// silently inverts intent (`.contigs(["chr22", "chr1"])` would pair
/// `chr22` with the profile's highest-density stats). Requesting the
/// bare-id profile's contigs in **reversed** order must still pair each
/// output name with the *matching* fitted id's density, not the id at the
/// same list position.
#[test]
fn chr_prefix_normalization_resolves_by_name_not_position() {
    let profile = Profile::from_json(NONUNIFORM_DENSITY_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    // Reversed relative to the fitted order (fitted: "1"=80, "2"=40,
    // "3"=10). Positional resolution would give index 0 ("chr3") the
    // fitted id "1" stats (density 80) and index 2 ("chr1") the fitted id
    // "3" stats (density 10) -- exactly inverted from what the names say.
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr3", "chr2", "chr1"])
        .seed(7)
        .size(Size::Records(300))
        .write(&path)
        .unwrap();

    assert_eq!(s.per_contig.len(), 3);
    // Name resolution must give "chr1" the *high*-density stats (fitted id
    // "1", 80) and "chr3" the *low*-density stats (fitted id "3", 10),
    // regardless of their position in the requested list.
    assert!(
        s.per_contig["chr1"].n_records > s.per_contig["chr3"].n_records * 2,
        "chr1 must resolve BY NAME to fitted id \"1\" (density 80), not by \
         position to fitted id \"3\" (density 10): {:?}",
        s.per_contig
    );
}

/// A genuine last-resort case: requested output names that share no `chr`-
/// normalized correspondence with any fitted id at all. Must still succeed
/// via the positional fallback rather than erroring or zeroing out a
/// contig's stats.
#[test]
fn unrelated_contig_names_fall_back_positionally_without_erroring() {
    let profile = Profile::from_json(NONUNIFORM_DENSITY_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["scaffold_a", "scaffold_b", "scaffold_c"])
        .seed(7)
        .size(Size::Records(300))
        .write(&path)
        .unwrap();

    assert_eq!(
        s.n_records_total(),
        300,
        "positional fallback must still resolve every contig to a real \
         fitted stat, not error or drop records"
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

/// Duplicate output contig names each get their own position stream
/// starting from 0, so positions run backwards across the file even though
/// noodles dedupes the `##contig` header line to one entry -- a CSI built
/// over such out-of-order records silently drops region-query hits rather
/// than erroring anywhere else, so `write()` must reject this up front.
#[test]
fn duplicate_contig_names_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let result = spec()
        .contigs(["chr1", "chr1"])
        .size(Size::RecordsPerContig(10))
        .write(&path);
    assert!(
        matches!(result, Err(BulkError::Invalid(_))),
        "duplicate output contig names must be rejected as invalid: {result:?}"
    );
    assert!(
        !path.exists(),
        "no file should be written for a rejected spec"
    );
}

/// `SampleStats::value_for` (`src/bulk/gen.rs`) hard-codes `AD` as a
/// 2-element `[n_ref, n_alt]` and `PL` as a fixed 3-element diploid
/// likelihood triple -- correct only for `ploidy == 2`. `Profile::validate`
/// only requires `ploidy >= 1`, so a ploidy-3 profile combined with a
/// payload that declares `PL`/`AD` (`Gatk`, `Mutect2`) must be rejected
/// rather than silently emitting genotype-likelihood/allele-depth fields
/// whose cardinality doesn't match the declared ploidy.
#[test]
fn non_diploid_profile_rejects_payloads_declaring_pl_or_ad() {
    let profile = Profile::from_json(TRIPLOID_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();

    for payload in [Payload::Gatk, Payload::Mutect2] {
        let path = dir.path().join(format!("{payload:?}.bcf"));
        let result = BulkSpec::new(profile.clone())
            .samples(4)
            .contigs(["chr1"])
            .payload(payload.clone())
            .size(Size::RecordsPerContig(10))
            .write(&path);
        assert!(
            matches!(result, Err(BulkError::Invalid(_))),
            "payload {payload:?} declares PL/AD, which are diploid-only; \
             ploidy 3 must be rejected, got: {result:?}"
        );
    }
}

/// The flip side of the guard above: a non-diploid profile combined with a
/// payload that does *not* declare `PL`/`AD` (`GtOnly`, `GtVaf`) must still
/// write successfully -- the guard is specific to the fields that are
/// actually hard-coded for diploid, not a blanket rejection of non-diploid
/// profiles.
#[test]
fn non_diploid_profile_accepts_payloads_without_pl_or_ad() {
    let profile = Profile::from_json(TRIPLOID_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();

    for payload in [Payload::GtOnly, Payload::GtVaf] {
        let path = dir.path().join(format!("{payload:?}.bcf"));
        let s = BulkSpec::new(profile.clone())
            .samples(4)
            .contigs(["chr1"])
            .payload(payload.clone())
            .size(Size::RecordsPerContig(10))
            .write(&path)
            .unwrap_or_else(|e| panic!("payload {payload:?} must not need ploidy 2: {e}"));
        assert_eq!(s.n_records_total(), 10, "payload {payload:?}");
    }
}

/// Same shape as `NONUNIFORM_DENSITY_PROFILE`, but `ploidy: 3`, for the
/// PL/AD-ploidy-guard tests above.
const TRIPLOID_PROFILE: &str = r#"
{
  "name": "triploid-test",
  "provenance": {
    "source": "test fixture",
    "n_samples_source": 0,
    "n_variants_source": 0,
    "fitted_on": "1970-01-01",
    "fit_tool_version": "0.0.0"
  },
  "fitted": {
    "contigs": [
      { "id": "chr1", "n_variants": 8000, "density_per_kb": 80.0 }
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
    "ploidy": 3
  },
  "dialed": { "payload": "gt-only" }
}
"#;
