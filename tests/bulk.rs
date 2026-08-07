#![cfg(feature = "bulk")]

use std::collections::BTreeMap;
use std::num::NonZero;
use vcfixture::bulk::{BulkError, BulkSpec, Format, Payload, Profile, Size};

/// Single source of truth for the cohort width `spec()` builds with. Read by
/// `same_seed_gives_byte_identical_output_across_thread_counts`'s vacuity
/// guard via `BulkSpec::block_records(SAMPLES, 2)` -- a mirrored literal
/// there has already silently regressed that guard twice (see its doc
/// comment), so the guard must derive from this constant instead of
/// repeating it.
const SAMPLES: usize = 8;

fn spec() -> BulkSpec {
    BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
        .samples(SAMPLES)
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
    // `BulkSpec::block_records(SAMPLES, 2)` (`src/bulk/mod.rs`) is currently
    // 500 records at this cohort width (8 samples, ploidy 2) — the
    // granularity at which `BulkSpec::stream_contigs`' rayon `par_iter` has
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
    // This guard reads `BulkSpec::block_records` directly rather than
    // mirroring its value as a local literal: a mirrored literal would
    // silently stop catching a vacuous test the moment the real sizing
    // changed (it already regressed twice in this branch when the value was
    // a flat constant), since the guard would keep comparing against its
    // own stale copy instead of the source of truth. Block size is now
    // computed from cohort width rather than constant, so the guard passes
    // `SAMPLES` (the same constant `spec()` builds with) and ploidy 2 (the
    // germline-1kgp profile's dialed ploidy -- a real fitted profile, not a
    // placeholder; see `germline_profile_is_really_fitted_not_placeholder`
    // in `src/bulk/profile.rs`) rather than reading a bare constant. See
    // `src/bulk/writer.rs`'s
    // `output_is_byte_identical_regardless_of_worker_count` for the same
    // idiom applied to the writer's own compression layer.
    const RECORDS_PER_CONTIG: u64 = 2_500;
    const MAX_BUF_SIZE: u64 = 65_498; // bgzf's uncompressed staging buffer size

    let blocks_per_contig = RECORDS_PER_CONTIG.div_ceil(BulkSpec::block_records(SAMPLES, 2));
    assert!(
        blocks_per_contig > 1,
        "test is vacuous: {RECORDS_PER_CONTIG} records/contig gives only \
         {blocks_per_contig} rayon block(s) per contig, so there is nothing \
         for `par_iter` to reorder"
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
    // `Size::Target`'s overshoot is a *proportional* margin, not a fixed
    // byte budget: each corrective round in `BulkSpec::resolve_target_counts`
    // refits a local slope `k_eff` from the last two rounds' actual
    // records/bytes, tops up by `shortfall / k_eff` plus a 2% margin, and
    // doubles the top-up on a round that made no byte progress. Observed on
    // the order of ~9.4% in practice -- a 4 MB target overshot by +377,520
    // bytes, i.e. ~9.0%. An absolute cap (e.g. `target + 256
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
/// 0 is legal and means "generate nothing here" -- such a contig produces
/// no non-empty `BlockSummary` to merge, so it has no `per_contig` entry at
/// all and must be read through `get`, not indexed.
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
    assert!(
        matches!(&err, BulkError::PerContigMissing(missing) if missing.iter().any(|c| c == "chr2")),
        "must name the contig with no count: {msg}"
    );
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
    assert!(
        matches!(&err, BulkError::PerContigUnknown(unknown) if unknown.iter().any(|c| c == "1")),
        "must name the unrequested key: {msg}"
    );
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

#[test]
fn layout_span_equals_the_realized_max_position() {
    // Spans are now computed by a pass that never generates genotypes, so
    // a divergence would silently declare a wrong ##contig length. Check
    // the declared length against what was actually written, across seeds,
    // cohort widths, and record counts.
    //
    // Every case here is sized to at least 3 blocks per contig at
    // `BulkSpec::block_records(samples, 2)` (500 records/block at these
    // cohort widths) -- with only 1 block per contig, every block's
    // absolute offset is 0 and this assertion degenerates into a pure
    // identity that cannot catch anything. `(1, 2, 1500)` and `(99, 300,
    // 2000)` land on an exact multiple of 500 (every block, including the
    // last, is full-width); `(7, 37, 1300)` leaves a partial final block
    // (500, 500, 300) -- both shapes are covered because `ContigLayout::
    // block_len`'s `.min(n_records - start)` only has anything to clamp on
    // the partial-tail case.
    //
    // NOTE ON STRENGTH: this test is deliberately *weak*, and must not be
    // mistaken for evidence that the block pipeline is correct. The write
    // pass takes `offsets[b]` and `block_len(b)` from the same
    // `ContigLayout` whose `block_spans` `span()` sums, so an error in an
    // interior block's span cancels exactly on both sides and this
    // assertion still passes. What it does pin is the *plumbing*: that the
    // header is built from the layout the records were actually written
    // under, and that the final block's span is not dropped or
    // double-counted. The genuinely independent cross-check -- the full
    // position vector against a gap walk recomputed from `block_rng` +
    // `Samplers::gap`, never consulting `ContigLayout` -- is
    // `realized_positions_match_an_independently_recomputed_gap_walk`
    // below.
    for (seed, samples, records) in [(1u64, 2usize, 1500u64), (7, 37, 1300), (99, 300, 2000)] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bcf");
        let summary = BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
            .samples(samples)
            .contigs(["chr1", "chr2"])
            .size(Size::RecordsPerContig(records))
            .payload(Payload::GtOnly)
            .seed(seed)
            .workers(NonZero::new(3).unwrap())
            .write(&path)
            .unwrap();

        let mut r = noodles_bcf::io::reader::Builder::default()
            .build_from_path(&path)
            .unwrap();
        let header = r.read_header().unwrap();

        for (id, contig) in header.contigs() {
            let declared = contig.length().expect("contig length is declared") as u64;
            let observed = summary.per_contig[id.as_str()].pos_max;
            assert_eq!(
                declared, observed,
                "declared length must equal the realized max position \
                 (seed={seed}, samples={samples}, records={records}, contig={id})"
            );
        }
    }
}

#[test]
fn positions_do_not_depend_on_payload_and_are_stable_across_widths_sharing_a_block_size() {
    // Positions are independent of payload unconditionally: `Payload` only
    // selects which FORMAT keys `to_record_buf` renders, and never touches
    // `Stream::Position`'s gap draws, which come from `Samplers::gap` alone.
    //
    // Positions are NOT independent of cohort width in general. Cell-based
    // block sizing (`BulkSpec::block_records`) deliberately makes block
    // boundaries a function of `n_samples` once a cohort is wide enough to
    // shrink below `MAX_BLOCK_RECORDS` (diploid, over ~4,000 samples) --
    // each block reseeds its `Stream::Position` RNG from `block_idx` alone,
    // so a different partition draws a different sequence and DOES move
    // positions. This test can only compare cohort widths that happen to
    // land on the same `block_records`, which the assertion below checks
    // and documents explicitly rather than relying on it silently.
    assert_eq!(
        BulkSpec::block_records(2, 2),
        BulkSpec::block_records(64, 2),
        "this test compares cohort widths 2 and 64 only because they \
         currently share a block size (both saturate MAX_BLOCK_RECORDS at \
         these small cell counts); if TARGET_CELLS_PER_BLOCK or \
         MAX_BLOCK_RECORDS ever change so these two widths land in \
         different block sizes, the width comparison below stops being \
         valid and this must fail loudly here rather than silently \
         comparing incomparable partitions"
    );

    fn positions(samples: usize, payload: Payload) -> Vec<u64> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bcf");
        BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
            .samples(samples)
            .contigs(["chr1"])
            .size(Size::RecordsPerContig(300))
            .payload(payload)
            .seed(5)
            .workers(NonZero::new(2).unwrap())
            .write(&path)
            .unwrap();

        let mut r = noodles_bcf::io::reader::Builder::default()
            .build_from_path(&path)
            .unwrap();
        let _ = r.read_header().unwrap();
        r.records()
            .map(|rec| usize::from(rec.unwrap().variant_start().unwrap().unwrap()) as u64)
            .collect()
    }

    let a = positions(2, Payload::GtOnly);
    assert_eq!(
        a,
        positions(64, Payload::GtOnly),
        "cohort width moved positions, despite both widths sharing a block size"
    );
    assert_eq!(a, positions(2, Payload::Gatk), "payload moved positions");
}

#[test]
fn output_is_byte_identical_across_thread_counts_and_chunkings() {
    // The oracle for the block pipeline: whatever the worker count (which
    // sets both the rayon pool size and the in-flight chunk width), the
    // bytes and the Summary must be identical. Worker count 1 is the
    // serial reference — one block encoded and written at a time.
    //
    // The spec dimensions are constants so that the vacuity guards below
    // derive from the very values `run` builds with, rather than mirroring
    // them -- a mirrored literal stops catching a vacuous test the moment
    // the real sizing changes, which has already silently regressed a
    // guard twice in this branch.
    const N_SAMPLES: usize = 40;
    const RECORDS_PER_CONTIG: u64 = 1_500;
    /// Worker counts under test, ascending. Index 0 is the serial
    /// reference; it is also the count with the *narrowest* chunk
    /// (`chunk_blocks = 2 * workers`), so it is what the chunk-boundary
    /// guard below is computed from.
    const WORKERS: [usize; 4] = [1, 2, 5, 16];

    fn run(workers: usize) -> (Vec<u8>, vcfixture::bulk::Summary) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bcf");
        let summary = BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
            .samples(N_SAMPLES)
            .contigs(["chr1", "chr2"])
            .size(Size::RecordsPerContig(RECORDS_PER_CONTIG))
            .payload(Payload::Gatk)
            .seed(21)
            .workers(NonZero::new(workers).unwrap())
            .write(&path)
            .unwrap();
        (std::fs::read(&path).unwrap(), summary)
    }

    // Vacuity guard 1 -- RAYON BLOCKS. This test's whole subject is the
    // block fan-out, and with a single block per contig `par_iter` has
    // nothing to compute out of order and `.collect()` nothing to reassemble:
    // the test would pass no matter how badly ordering was broken. Ploidy 2
    // is the germline-1kgp profile's dialed ploidy (see
    // `germline_profile_is_really_fitted_not_placeholder` in
    // `src/bulk/profile.rs`). Same idiom as
    // `same_seed_gives_byte_identical_output_across_thread_counts` above.
    let blocks_per_contig = RECORDS_PER_CONTIG.div_ceil(BulkSpec::block_records(N_SAMPLES, 2));
    assert!(
        blocks_per_contig > 1,
        "test is vacuous: {RECORDS_PER_CONTIG} records/contig at {N_SAMPLES} \
         samples gives only {blocks_per_contig} rayon block(s) per contig, so \
         there is nothing for the fan-out to reorder"
    );

    // Vacuity guard 2 -- CHUNK BOUNDARIES. The other half of what this test
    // claims: that `chunk_blocks = 2 * workers` bounds only how many blocks
    // are in flight, never the bytes. If every worker count swallowed a
    // whole contig in one chunk, no contig would ever be split across a
    // chunk boundary and the chunking half would go untested. Only the
    // smallest worker count needs to split (a wider chunk at a larger
    // worker count is the degenerate case this is comparing *against*), so
    // the guard is over `WORKERS[0]` alone.
    let narrowest_chunk_blocks = 2 * WORKERS[0] as u64;
    assert!(
        narrowest_chunk_blocks < blocks_per_contig,
        "test is vacuous for chunking: at {} workers a chunk holds \
         {narrowest_chunk_blocks} blocks and a contig is only \
         {blocks_per_contig} blocks, so no contig is ever split across a \
         chunk boundary",
        WORKERS[0]
    );

    let (bytes1, sum1) = run(WORKERS[0]);
    for workers in WORKERS.into_iter().skip(1) {
        let (bytes, sum) = run(workers);
        assert_eq!(
            bytes1, bytes,
            "output must be byte-identical at {workers} workers"
        );
        assert_eq!(
            sum1.genotype_checksum, sum.genotype_checksum,
            "checksum must be identical at {workers} workers"
        );
        assert_eq!(sum1.per_contig, sum.per_contig);
        assert_eq!(sum1.class_counts, sum.class_counts);
        assert_eq!(sum1.n_alleles_nonref, sum.n_alleles_nonref);
    }

    // Guard against vacuity: the payload must span several bgzf blocks, or
    // neither the rayon block fan-out nor the writer's own compression
    // pool would have anything to reorder.
    //
    // The threshold is bgzf's *uncompressed* staging buffer (`MAX_BUF_SIZE`,
    // ~65,498 bytes): a block is dispatched once that many uncompressed
    // bytes have accumulated. So the quantity to compare against it is the
    // decompressed payload, not the compressed file length -- at level 6 on
    // this repetitive data the file is ~57 KB while the payload behind it
    // is megabytes, and comparing the compressed length would have made
    // this guard fail on a perfectly non-vacuous test. Same idiom as
    // `same_seed_gives_byte_identical_output_across_thread_counts` above
    // and `writer::tests::output_is_byte_identical_regardless_of_worker_count`.
    const MAX_BUF_SIZE: usize = 65_498;
    let mut decompressed = Vec::new();
    std::io::Read::read_to_end(
        &mut noodles_bgzf::io::Reader::new(std::io::Cursor::new(&bytes1)),
        &mut decompressed,
    )
    .unwrap();
    assert!(
        decompressed.len() > 3 * MAX_BUF_SIZE,
        "test payload ({} bytes uncompressed, {} compressed) must exceed \
         several bgzf blocks ({MAX_BUF_SIZE} bytes each)",
        decompressed.len(),
        bytes1.len()
    );
}

#[test]
fn declared_contig_lengths_are_never_zero() {
    // `##contig length` is a 1-based coordinate bound, so `length=0`
    // declares a contig that cannot hold even position 1 -- a malformed
    // header, and one that region-query tooling reads as an empty
    // reference sequence.
    //
    // This is reachable, not hypothetical: `ContigLayout::span()` sums an
    // empty block list to `0` for a zero-record contig, and
    // `Size::Records(total)` with `total` below the contig count really
    // does produce zero-record contigs (`distribute_by_n_variants` floors
    // the proportional split, then largest-remainder hands the single
    // record to exactly one contig). The pre-#22 `contig_span` hid this
    // behind an `unwrap_or(1)`; the floor now lives in
    // `BulkSpec::build_header`.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bcf");
    let summary = BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
        .samples(2)
        .contigs(["chr1", "chr2", "chr3"])
        .size(Size::Records(1))
        .payload(Payload::GtOnly)
        .seed(4)
        .workers(NonZero::new(2).unwrap())
        .write(&path)
        .unwrap();

    assert_eq!(summary.n_records_total(), 1);
    // Vacuity guard: if every contig somehow got a record, there is no
    // zero-span contig here and the floor is never exercised.
    assert!(
        summary.per_contig.len() < 3,
        "test is vacuous: no contig ended up with zero records ({:?})",
        summary.per_contig
    );

    let mut r = noodles_bcf::io::reader::Builder::default()
        .build_from_path(&path)
        .unwrap();
    let header = r.read_header().unwrap();
    assert_eq!(
        header.contigs().len(),
        3,
        "every requested contig must still be declared, even at zero records"
    );
    for (id, contig) in header.contigs() {
        let declared = contig.length().expect("contig length is declared");
        assert!(
            declared >= 1,
            "contig {id} declared length {declared}, but a ##contig length \
             must be >= 1 even for a contig with no records"
        );
    }
}

#[test]
fn realized_positions_match_an_independently_recomputed_gap_walk() {
    // The independent cross-check on the block pipeline.
    //
    // `layout_span_equals_the_realized_max_position` cannot catch an
    // interior block-offset error: the write pass reads `offsets[b]` from
    // the very `ContigLayout` whose `block_spans` `span()` sums, so both
    // sides of that assertion move together. This test never touches
    // `ContigLayout` at all. It rebuilds the expected position of *every*
    // record from the public primitives the design is specified in --
    // `block_rng(seed, block_idx, Stream::Position)` plus `Samplers::gap`,
    // with its own prefix sum over blocks -- and compares the full vector,
    // in order, against what was actually written to the file.
    //
    // Why that is not an identity: the expectation is derived from the
    // *specification* of the partition (block-index arithmetic over
    // `BulkSpec::block_records` and `BulkSpec::CONTIG_BLOCK_STRIDE`), not
    // from the artifact `compute_layouts` produced. A block seeded from the
    // wrong `block_idx`, a block handed the wrong record count, a wrong
    // per-block offset, a mis-sliced `block_spans` copy-back across
    // contigs, or a chunk written out of order all change interior
    // positions and fail here -- including the cases that leave the
    // contig's maximum position untouched (e.g. two interior blocks'
    // spans transposed in the prefix sum).
    use vcfixture::bulk::generate::{block_rng, Stream};
    use vcfixture::bulk::sample::Samplers;

    // (seed, samples, records_per_contig, contigs)
    let cases: [(u64, usize, u64, &[&str]); 2] = [
        // Narrow cohort: block_records saturates at MAX_BLOCK_RECORDS
        // (500), so 1300 records is 500/500/300 -- a partial final block --
        // across two contigs, which also exercises CONTIG_BLOCK_STRIDE.
        (3, 5, 1_300, &["chr1", "chr2"]),
        // Wide cohort: block_records is cell-bounded (4e6 / (8000*2) =
        // 250), so 650 records is 250/250/150 under a partition a flat-500
        // write pass would get wrong.
        (8, 8_000, 650, &["chr1"]),
    ];

    for (seed, samples, records, contigs) in cases {
        let profile = Profile::builtin("germline-1kgp").unwrap();
        let ploidy = profile.dialed.ploidy;
        let block_records = BulkSpec::block_records(samples, ploidy);
        let n_blocks = records.div_ceil(block_records);
        assert!(
            n_blocks >= 3,
            "test is vacuous at samples={samples}, records={records}: \
             {n_blocks} block(s) leaves no interior block whose offset \
             could be wrong"
        );

        // The expectation: gaps only, straight from the block streams.
        let samplers = Samplers::new(
            &profile.fitted,
            2 * profile.provenance.n_samples_source as u64,
        )
        .unwrap();
        let mut expected: BTreeMap<String, Vec<u64>> = BTreeMap::new();
        for (ci, id) in contigs.iter().enumerate() {
            let mut positions = Vec::with_capacity(records as usize);
            let mut offset = 0u64;
            for b in 0..n_blocks {
                let block_idx = ci as u64 * BulkSpec::CONTIG_BLOCK_STRIDE + b;
                let mut rng = block_rng(seed, block_idx, Stream::Position);
                let count = block_records.min(records - b * block_records);
                let mut local = 0u64;
                for _ in 0..count {
                    local += samplers.gap(&mut rng);
                    positions.push(offset + local);
                }
                offset += local;
            }
            expected.insert((*id).to_string(), positions);
        }

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bcf");
        BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
            .samples(samples)
            .contigs(contigs.iter().copied())
            .size(Size::RecordsPerContig(records))
            .payload(Payload::GtOnly)
            .seed(seed)
            .workers(NonZero::new(4).unwrap())
            .write(&path)
            .unwrap();

        let mut r = noodles_bcf::io::reader::Builder::default()
            .build_from_path(&path)
            .unwrap();
        let header = r.read_header().unwrap();
        let names: Vec<String> = header.contigs().keys().cloned().collect();

        let mut realized: BTreeMap<String, Vec<u64>> = BTreeMap::new();
        for rec in r.records() {
            let rec = rec.unwrap();
            let chrom = names[rec.reference_sequence_id().unwrap()].clone();
            let pos = usize::from(rec.variant_start().unwrap().unwrap()) as u64;
            realized.entry(chrom).or_default().push(pos);
        }

        assert_eq!(
            realized, expected,
            "realized positions must match the independently recomputed gap \
             walk (seed={seed}, samples={samples}, records={records})"
        );

        // The header's declared length must agree with that same
        // independent walk, not merely with whatever the write pass did.
        for (id, contig) in header.contigs() {
            let declared = contig.length().expect("contig length is declared") as u64;
            assert_eq!(
                declared,
                *expected[id.as_str()].last().unwrap(),
                "declared length for {id} must equal the independently \
                 recomputed final position (seed={seed}, samples={samples})"
            );
        }
    }
}

/// The #27 regression guard: a failed write must leave nothing behind at the
/// caller's destination — no truncated output, no stray `.csi`, no
/// `.summary.json`.
///
/// Before the temp-then-promote change, `BulkWriter::create` made the
/// destination file up front and a mid-stream failure propagated past
/// `finish_and_index`, leaving a headerless, un-indexed artifact sitting
/// where a valid one was expected. The caller did see the `Err`, so this is
/// not about error reporting — it is about the debris.
///
/// # Why this specific setup
///
/// The failure is induced by making the output *directory* read-only while
/// the destination file itself already exists and is writable. On Unix,
/// opening an existing file for writing needs write permission on the file,
/// not on its directory — so the old `BulkWriter::create(path, ..)` succeeded
/// and **truncated the destination**, then failed later when
/// `finish_and_index` tried to create the `.csi` entry. That is exactly the
/// issue #27 scenario: a retry into a path that already holds a good file
/// leaves a corrupt one there instead.
///
/// A read-only directory with no pre-existing destination would be a vacuous
/// test — the old code's `create` would fail up front and leave nothing
/// behind either. The pre-existing file is what makes this discriminate.
#[test]
#[cfg(unix)]
fn failed_write_leaves_the_destination_untouched() {
    use std::os::unix::fs::PermissionsExt;

    let dir = tempfile::tempdir().unwrap();
    let out_dir = dir.path().join("out");
    std::fs::create_dir(&out_dir).unwrap();
    let path = out_dir.join("a.bcf");

    const SENTINEL: &[u8] = b"a previously generated file that must survive a failed retry";
    std::fs::write(&path, SENTINEL).unwrap();

    let set_mode = |mode: u32| {
        let mut perms = std::fs::metadata(&out_dir).unwrap().permissions();
        perms.set_mode(mode);
        std::fs::set_permissions(&out_dir, perms).unwrap();
    };

    set_mode(0o500); // r-x: traversable, and existing files stay writable.
    let result = spec().size(Size::RecordsPerContig(600)).write(&path);
    set_mode(0o700); // Restore so the assertions and cleanup can proceed.

    assert!(
        result.is_err(),
        "writing into a read-only directory must fail"
    );
    assert_eq!(
        std::fs::read(&path).unwrap(),
        SENTINEL,
        "a failed write must not touch the destination file"
    );

    let mut leftovers: Vec<_> = std::fs::read_dir(&out_dir)
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .collect();
    leftovers.sort();
    assert_eq!(
        leftovers,
        vec![std::ffi::OsString::from("a.bcf")],
        "a failed write must leave no debris beside the destination"
    );
}
