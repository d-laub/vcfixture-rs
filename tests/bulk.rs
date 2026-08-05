#![cfg(feature = "bulk")]

use std::collections::BTreeMap;
use std::num::NonZero;
use vcfixture::bulk::{BulkError, BulkSpec, Format, Payload, Profile, Size};

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
    // `BulkSpec::BLOCK_SIZE` (`src/bulk/mod.rs`) is currently 500 records —
    // the granularity at which `generate_contig`'s rayon `into_par_iter` has
    // anything to reorder. At `RecordsPerContig(50)` (the value this test
    // used before this fix), `50.div_ceil(500) == 1` block per contig: there
    // is nothing for rayon to reorder, so both parallel paths (this crate's
    // own block-level rayon split, and the writer's bgzf multithreaded
    // compression, which also scales with `workers`) are structurally
    // incapable of failing regardless of whether determinism actually
    // holds. `RECORDS_PER_CONTIG` below is sized well past that boundary,
    // and the assertions further down pin both dimensions of scale down so
    // a future edit can't silently shrink this test back into vacuity.
    //
    // This guard reads `BulkSpec::BLOCK_SIZE` directly rather than mirroring
    // its value as a local literal: a mirrored literal would silently stop
    // catching a vacuous test the moment the real constant changed (it
    // already regressed twice in this branch), since the guard would keep
    // comparing against its own stale copy instead of the source of truth.
    // See `src/bulk/writer.rs`'s `output_is_byte_identical_regardless_of_
    // worker_count` for the same idiom applied to the writer's own
    // compression layer.
    const RECORDS_PER_CONTIG: u64 = 2_500;
    const MAX_BUF_SIZE: u64 = 65_498; // bgzf's uncompressed staging buffer size

    let blocks_per_contig = RECORDS_PER_CONTIG.div_ceil(BulkSpec::BLOCK_SIZE);
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
    // `Size::Target`'s overshoot is a *proportional* margin (the 15% top-up
    // margin in `BulkSpec::resolve_target_counts`, observed on the order of
    // ~9.4% in practice -- a 4 MB target overshot by +377,520 bytes, i.e.
    // ~9.0%), not a fixed byte budget. An absolute cap (e.g. `target + 256
    // KiB`) only looks like it bounds the real behavior at whatever target
    // size makes the two coincide -- 256 KiB happens to be 50% of the 512
    // KiB target this test used to use, which is why that version of this
    // assertion passed without actually checking proportionality. Bound it
    // as a percentage of `target` instead, comfortably above the ~9.4%
    // observed overshoot so this doesn't flake, but still tight enough to
    // catch a regression back toward unbounded overshoot.
    const MAX_OVERSHOOT_FRACTION: f64 = 0.25;

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let target = 512 * 1024;
    spec().size(Size::Target(target)).write(&path).unwrap();
    let got = std::fs::metadata(&path).unwrap().len();
    assert!(got >= target, "got {got} < target {target}");
    let max_allowed = target + (target as f64 * MAX_OVERSHOOT_FRACTION) as u64;
    assert!(
        got < max_allowed,
        "overshoot too large: {got} vs target {target} (max allowed {max_allowed}, \
         i.e. {MAX_OVERSHOOT_FRACTION} of target)"
    );
}

/// Guards the calibrate+promote rewrite of `resolve_target_counts`: two runs
/// of the same seed and target must produce byte-identical output. This is
/// weaker than asserting the resolved per-contig counts directly (no new
/// `Summary` accessor needed for that), but it is sufficient to catch any
/// nondeterminism the rewrite could introduce -- e.g. the corrective rounds
/// depending on wall-clock-observed byte counts rather than purely on
/// `(seed, contig, count)`.
#[test]
fn target_size_is_byte_identical_across_runs() {
    let dir = tempfile::tempdir().unwrap();
    let profile = Profile::from_json(NONUNIFORM_DENSITY_PROFILE).unwrap();
    let run = |name: &str| {
        let out = dir.path().join(name);
        BulkSpec::new(profile.clone())
            .samples(4)
            .contigs(["1", "2", "3"])
            .format(Format::VcfGz)
            .size(Size::Target(512 * 1024))
            .seed(99)
            .write(&out)
            .unwrap();
        std::fs::read(&out).unwrap()
    };
    let a = run("a.vcf.gz");
    let b = run("b.vcf.gz");
    assert!(a.len() as u64 >= 512 * 1024, "must reach target");
    assert_eq!(a, b, "calibrate+promote must be byte-identical run to run");
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
/// it splits `total` across contigs proportional to each contig's fitted
/// `n_variants`. This covers the exactness half of that contract: the split
/// must sum to exactly `total` and must be non-uniform for a profile with
/// non-uniform per-contig statistics (the placeholder profile's 3 contigs
/// are all identical, so a custom profile is needed to exercise weighting
/// at all). `records_total_splits_by_variants_not_density` covers *which*
/// statistic does the weighting -- this fixture cannot, since its
/// `n_variants` and `density_per_kb` are both 8:4:1.
#[test]
fn records_total_sums_to_exactly_the_requested_total() {
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
    // chr1's fitted n_variants (8000) is 8x chr3's (1000), so it must get
    // noticeably more records, not an even 100/100/100 split.
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
    // "chr1" normalizes to "1" and matches fitted id "1" (n_variants 8000)
    // by name, so it must still outweigh "chr3" (fitted id "3",
    // n_variants 1000).
    assert!(
        s.per_contig["chr1"].n_records > s.per_contig["chr3"].n_records * 2,
        "chr1 (name-resolved to fitted id \"1\", n_variants 8000) must \
         outweigh chr3 (fitted id \"3\", n_variants 1000): {:?}",
        s.per_contig
    );
}

/// The specific failure mode Important-4 fixed: with *only* positional
/// resolution, requesting output contigs out of the profile's fitted order
/// silently inverts intent (`.contigs(["chr22", "chr1"])` would pair
/// `chr22` with the profile's highest-count stats). Requesting the
/// bare-id profile's contigs in **reversed** order must still pair each
/// output name with the *matching* fitted id's variant count, not the id at
/// the same list position.
#[test]
fn chr_prefix_normalization_resolves_by_name_not_position() {
    let profile = Profile::from_json(NONUNIFORM_DENSITY_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    // Reversed relative to the fitted order (fitted n_variants: "1"=8000,
    // "2"=4000, "3"=1000). Positional resolution would give index 0
    // ("chr3") the fitted id "1" stats (n_variants 8000) and index 2
    // ("chr1") the fitted id "3" stats (n_variants 1000) -- exactly
    // inverted from what the names say.
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr3", "chr2", "chr1"])
        .seed(7)
        .size(Size::Records(300))
        .write(&path)
        .unwrap();

    assert_eq!(s.per_contig.len(), 3);
    // Name resolution must give "chr1" the *high*-count stats (fitted id
    // "1", n_variants 8000) and "chr3" the *low*-count stats (fitted id
    // "3", n_variants 1000), regardless of their position in the requested
    // list.
    assert!(
        s.per_contig["chr1"].n_records > s.per_contig["chr3"].n_records * 2,
        "chr1 must resolve BY NAME to fitted id \"1\" (n_variants 8000), \
         not by position to fitted id \"3\" (n_variants 1000): {:?}",
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

/// The regression guard issue #15 asks for. Densities are identical across
/// this profile's three contigs, so a density-weighted split gives a dead
/// even 300/300/300; the shipped `n_variants` weight (6:2:1) gives
/// 600/200/100. 900 divides evenly by the 6+2+1 weight sum, so this can
/// assert exact counts without depending on largest-remainder tie-breaking.
#[test]
fn records_total_splits_by_variants_not_density() {
    let profile = Profile::from_json(SKEWED_VARIANTS_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1", "chr2", "chr3"])
        .seed(7)
        .size(Size::Records(900))
        .write(&path)
        .unwrap();

    assert_eq!(s.n_records_total(), 900, "must sum to exactly the request");
    let got = (
        s.per_contig["chr1"].n_records,
        s.per_contig["chr2"].n_records,
        s.per_contig["chr3"].n_records,
    );
    assert_eq!(
        got,
        (600, 200, 100),
        "must follow fitted n_variants (6:2:1), not fitted density (all \
         25.0/kb, which would give 300/300/300): {:?}",
        s.per_contig
    );
}

/// `Size::Target` resolves its per-contig counts through the same helper as
/// `Size::Records`, but via a different call path (`resolve_target_counts`).
/// Cover it separately so a change to one weight site cannot pass while the
/// other regresses.
#[test]
fn target_size_split_also_follows_variants() {
    let profile = Profile::from_json(SKEWED_VARIANTS_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1", "chr2", "chr3"])
        .seed(7)
        .size(Size::Target(512 * 1024))
        .write(&path)
        .unwrap();

    let (c1, c3) = (
        s.per_contig["chr1"].n_records,
        s.per_contig["chr3"].n_records,
    );
    // True ratio is 6:1. Assert only a direction and a loose magnitude --
    // the exact totals depend on how many calibration rounds `Size::Target`
    // needed, which is not what this test is pinning down.
    assert!(
        c1 > c3 * 3,
        "Size::Target must weight by n_variants (6:1 here), not by the \
         uniform density: chr1={c1} chr3={c3} {:?}",
        s.per_contig
    );
}

/// The core `Size::PerContig` contract: each requested contig gets exactly
/// the count it was given, with no profile-derived reweighting. A count of
/// 0 is legal and means "generate nothing here" -- such a contig never
/// reaches `Summary::observe`, so it has no `per_contig` entry at all and
/// must be read through `get`, not indexed.
#[test]
fn per_contig_gives_exact_requested_counts() {
    let profile = Profile::from_json(SKEWED_VARIANTS_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let counts = BTreeMap::from([
        ("chr1".to_string(), 37u64),
        ("chr2".to_string(), 250u64),
        ("chr3".to_string(), 0u64),
    ]);
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1", "chr2", "chr3"])
        .seed(7)
        .size(Size::PerContig(counts))
        .write(&path)
        .unwrap();

    assert_eq!(s.per_contig["chr1"].n_records, 37);
    assert_eq!(s.per_contig["chr2"].n_records, 250);
    assert_eq!(
        s.per_contig.get("chr3").map_or(0, |c| c.n_records),
        0,
        "a 0-count contig must produce no records: {:?}",
        s.per_contig
    );
    assert_eq!(s.n_records_total(), 287);
}

/// A requested contig with no entry in the map would otherwise generate
/// silently empty output. Reject it, and name it.
#[test]
fn per_contig_missing_contig_is_an_error() {
    let profile = Profile::from_json(SKEWED_VARIANTS_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let counts = BTreeMap::from([("chr1".to_string(), 10u64)]);
    let err = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1", "chr2"])
        .size(Size::PerContig(counts))
        .write(&path)
        .unwrap_err();

    let msg = err.to_string();
    assert!(matches!(err, BulkError::Invalid(_)), "{msg}");
    assert!(
        msg.contains("chr2"),
        "error must name the contig with no count: {msg}"
    );
}

/// The near-miss that pins down exact name matching: `"1"` is what the
/// profile calls this contig, and `resolve_contig_stat` *would* normalize
/// `"chr1"` onto it -- but `Size::PerContig` keys are matched exactly
/// against the requested output names, so `"1"` is an unknown key, not a
/// silently-accepted alias.
#[test]
fn per_contig_unknown_key_is_an_error() {
    let profile = Profile::from_json(SKEWED_VARIANTS_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let counts = BTreeMap::from([("chr1".to_string(), 10u64), ("1".to_string(), 10u64)]);
    let err = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1"])
        .size(Size::PerContig(counts))
        .write(&path)
        .unwrap_err();

    let msg = err.to_string();
    assert!(matches!(err, BulkError::Invalid(_)), "{msg}");
    assert!(
        msg.contains("\"1\""),
        "error must name the unrequested key: {msg}"
    );
}

/// Explicit counts must beat the profile's fitted shape outright, not be
/// blended with it: this map inverts SKEWED_VARIANTS_PROFILE's 6:2:1 skew,
/// giving the *least*-variant contig the most records.
#[test]
fn per_contig_ignores_profile_shape() {
    let profile = Profile::from_json(SKEWED_VARIANTS_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.bcf");
    let counts = BTreeMap::from([
        ("chr1".to_string(), 50u64),
        ("chr2".to_string(), 100u64),
        ("chr3".to_string(), 300u64),
    ]);
    let s = BulkSpec::new(profile)
        .samples(4)
        .contigs(["chr1", "chr2", "chr3"])
        .seed(7)
        .size(Size::PerContig(counts))
        .write(&path)
        .unwrap();

    assert_eq!(s.per_contig["chr1"].n_records, 50);
    assert_eq!(s.per_contig["chr2"].n_records, 100);
    assert_eq!(s.per_contig["chr3"].n_records, 300);
}

/// The fixture issue #15 needed and `NONUNIFORM_DENSITY_PROFILE` cannot
/// provide. Its two per-contig statistics deliberately **disagree**:
/// `density_per_kb` is identical (25.0) across all three contigs while
/// `n_variants` is 6:2:1. A split weighted by density is therefore dead
/// even, and a split weighted by variant counts is 6:2:1 -- so a test
/// against this profile can tell the two apart.
///
/// `NONUNIFORM_DENSITY_PROFILE` cannot: it declares both statistics in the
/// same 8:4:1 ratio, so every test using it passes under either weight.
/// That is exactly why the density-weighted split shipped unnoticed until
/// issue #15 measured a real corpus, and why this fixture exists.
const SKEWED_VARIANTS_PROFILE: &str = r#"
{
  "name": "skewed-variants-test",
  "provenance": {
    "source": "test fixture",
    "n_samples_source": 0,
    "n_variants_source": 0,
    "fitted_on": "1970-01-01",
    "fit_tool_version": "0.0.0"
  },
  "fitted": {
    "contigs": [
      { "id": "1", "n_variants": 6000, "density_per_kb": 25.0 },
      { "id": "2", "n_variants": 2000, "density_per_kb": 25.0 },
      { "id": "3", "n_variants": 1000, "density_per_kb": 25.0 }
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
    "phased_rate": 1.0
  },
  "dialed": { "payload": "gt-only", "ploidy": 2 }
}
"#;

/// Same profile as the two tests above: fitted contig ids are bare
/// (`"1"`/`"2"`/`"3"`), with 8x/4x/1x relative n_variants *and* density, so
/// the weighted split has an unambiguous direction under either statistic
/// -- see `SKEWED_VARIANTS_PROFILE` for a fixture that distinguishes them.
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
    "phased_rate": 1.0
  },
  "dialed": { "payload": "gt-only", "ploidy": 2 }
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
        matches!(&result, Err(BulkError::DuplicateContig(id)) if id == "chr1"),
        "duplicate output contig names must be rejected with DuplicateContig \
         naming the offending contig: {result:?}"
    );
    assert!(
        !path.exists(),
        "no file should be written for a rejected spec"
    );
}

/// `SampleStats::value_for` (`src/bulk/generate.rs`) hard-codes `AD` as a
/// 2-element `[n_ref, n_alt]` and `PL` as a fixed 3-element diploid
/// likelihood triple -- correct only for `ploidy == 2`. A ploidy-3 profile
/// combined with a payload that declares `PL`/`AD` (`Gatk`, `Mutect2`) must
/// be rejected by `Profile::validate` at parse time, rather than silently
/// emitting genotype-likelihood/allele-depth fields whose cardinality
/// doesn't match the declared ploidy.
#[test]
fn non_diploid_profile_rejects_payloads_declaring_pl_or_ad() {
    let mut profile = Profile::from_json(TRIPLOID_PROFILE).unwrap();

    for payload in [Payload::Gatk, Payload::Mutect2] {
        profile.dialed.payload = payload.clone();
        let result = profile.validate();
        assert!(
            matches!(
                result,
                Err(BulkError::PayloadPloidy {
                    ploidy: 3,
                    payload: ref p,
                }) if *p == payload
            ),
            "payload {payload:?} declares PL/AD, which are diploid-only; \
             ploidy 3 must be rejected with PayloadPloidy naming the payload \
             and the offending ploidy, got: {result:?}"
        );
    }
}

/// The flip side of the guard above: a non-diploid profile combined with a
/// payload that does *not* declare `PL`/`AD` (`GtOnly`, `GtVaf`) must still
/// pass `validate()` and write successfully -- the guard is specific to the
/// fields that are actually hard-coded for diploid, not a blanket rejection
/// of non-diploid profiles.
#[test]
fn non_diploid_profile_accepts_payloads_without_pl_or_ad() {
    let mut profile = Profile::from_json(TRIPLOID_PROFILE).unwrap();
    let dir = tempfile::tempdir().unwrap();

    for payload in [Payload::GtOnly, Payload::GtVaf] {
        profile.dialed.payload = payload.clone();
        profile
            .validate()
            .unwrap_or_else(|e| panic!("payload {payload:?} must not need ploidy 2: {e}"));

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
    "phased_rate": 1.0
  },
  "dialed": { "payload": "gt-only", "ploidy": 3 }
}
"#;

/// An empty contig list and a zero sample count are caller mistakes, not
/// profile mistakes -- they must not be reported as an invalid profile.
#[test]
fn empty_spec_dimensions_are_rejected_as_spec_errors() {
    let dir = tempfile::tempdir().unwrap();

    let no_contigs = spec()
        .contigs(Vec::<String>::new())
        .size(Size::RecordsPerContig(10))
        .write(dir.path().join("a.bcf"));
    assert!(
        matches!(no_contigs, Err(BulkError::NoContigs)),
        "an empty contig list must be a spec error: {no_contigs:?}"
    );

    let no_samples = spec()
        .samples(0)
        .size(Size::RecordsPerContig(10))
        .write(dir.path().join("b.bcf"));
    assert!(
        matches!(no_samples, Err(BulkError::NoSamples)),
        "a zero sample count must be a spec error: {no_samples:?}"
    );

    for e in [no_contigs, no_samples] {
        let msg = e.unwrap_err().to_string();
        assert!(
            !msg.starts_with("invalid profile:"),
            "a spec error must not blame the profile, got: {msg}"
        );
    }
}
