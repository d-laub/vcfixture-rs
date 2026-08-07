//! Golden byte-equality gate for bulk generation.
//!
//! Pins a digest of the generated artifact for every (format, payload)
//! combination. Refactors that are supposed to preserve output — the
//! scratch-buffer reuse of #26 and the temp-then-promote change of #27 —
//! must leave every one of these snapshots untouched.
//!
//! A snapshot change is a design failure, not a test to update. If noodles
//! is upgraded and the encoding legitimately changes, update these with an
//! explicit note in the commit message saying which upstream change caused
//! it.
//!
//! This is deliberately separate from
//! `same_seed_gives_byte_identical_output_across_thread_counts` in
//! `tests/bulk.rs`: that test compares two runs of the *same* code to each
//! other, so it cannot notice a refactor that changes the output of both
//! runs identically. Only a committed golden catches that.

#![cfg(feature = "bulk")]

use std::num::NonZero;
use std::path::Path;

use vcfixture::bulk::{BulkSpec, Format, Payload, Profile, Size};

/// Cohort width. `BulkSpec::block_records(8, 2)` is 500 records at this
/// width, and `RECORDS_PER_CONTIG` below is a multiple of it, so every run
/// spans several blocks and actually exercises the `map_init` path that the
/// scratch-buffer change rewrites. A single-block run would be structurally
/// incapable of catching a scratch-reuse bug that leaks state between
/// blocks.
const SAMPLES: usize = 8;
const RECORDS_PER_CONTIG: u64 = 1200;

/// FNV-1a, 64-bit.
///
/// A cryptographic digest would be overkill: this gate detects accidental
/// change, not adversarial collision. Keeping it inline also keeps this
/// file off `Cargo.toml`, which is what lets it land in parallel with the
/// mimalloc change.
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Byte length alongside the digest, so a change would have to collide
/// *and* preserve length to slip through.
fn digest(path: &Path) -> String {
    let bytes = std::fs::read(path).expect("output file must exist");
    format!("{} bytes, fnv1a64={:016x}", bytes.len(), fnv1a64(&bytes))
}

fn generate(format: Format, payload: Payload, name: &str) -> String {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join(name);
    BulkSpec::new(Profile::builtin("germline-1kgp").unwrap())
        .samples(SAMPLES)
        .contigs(["chr1", "chr2"])
        .payload(payload)
        .format(format)
        .size(Size::RecordsPerContig(RECORDS_PER_CONTIG))
        .seed(42)
        .workers(NonZero::new(2).unwrap())
        .write(&path)
        .unwrap();
    digest(&path)
}

macro_rules! golden {
    ($test_name:ident, $format:expr, $payload:expr, $file:expr) => {
        #[test]
        fn $test_name() {
            insta::assert_snapshot!(generate($format, $payload, $file));
        }
    };
}

// BCF exercises noodles' `encode_genotype` / `encode_genotype_str`; VCF and
// VcfGz exercise the *text* writer's rendering of the genotype value, which
// is the one risk in the scratch-buffer change that reading the encoder
// source cannot settle.
golden!(bcf_gt_only, Format::Bcf, Payload::GtOnly, "a.bcf");
golden!(bcf_gt_vaf, Format::Bcf, Payload::GtVaf, "a.bcf");
golden!(bcf_gatk, Format::Bcf, Payload::Gatk, "a.bcf");
golden!(bcf_mutect2, Format::Bcf, Payload::Mutect2, "a.bcf");

golden!(vcf_gt_only, Format::Vcf, Payload::GtOnly, "a.vcf");
golden!(vcf_gt_vaf, Format::Vcf, Payload::GtVaf, "a.vcf");
golden!(vcf_gatk, Format::Vcf, Payload::Gatk, "a.vcf");
golden!(vcf_mutect2, Format::Vcf, Payload::Mutect2, "a.vcf");

golden!(vcfgz_gt_only, Format::VcfGz, Payload::GtOnly, "a.vcf.gz");
golden!(vcfgz_gt_vaf, Format::VcfGz, Payload::GtVaf, "a.vcf.gz");
golden!(vcfgz_gatk, Format::VcfGz, Payload::Gatk, "a.vcf.gz");
golden!(vcfgz_mutect2, Format::VcfGz, Payload::Mutect2, "a.vcf.gz");
