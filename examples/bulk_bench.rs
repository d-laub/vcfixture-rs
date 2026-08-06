//! Sweep bulk generation over cohort width and record count, reporting
//! seconds per cell so numbers line up directly with the ladder in issue
//! #22.
//!
//! Run with:
//!   cargo run --release --example bulk_bench --features bulk
//!
//! Override the sweep with `VCFIXTURE_BENCH_SAMPLES` /
//! `VCFIXTURE_BENCH_RECORDS` (comma-separated), e.g.
//!   VCFIXTURE_BENCH_SAMPLES=4000 VCFIXTURE_BENCH_RECORDS=250000 \
//!     cargo run --release --example bulk_bench --features bulk
//!
//! Override the worker count (passed to both the rayon pool and the bgzf
//! writer) with `VCFIXTURE_BENCH_WORKERS`, e.g. to measure scaling:
//!   VCFIXTURE_BENCH_WORKERS=1 cargo run --release --example bulk_bench --features bulk

use std::env;
use std::num::NonZeroUsize;
use std::time::Instant;

use vcfixture::bulk::{BulkSpec, Payload, Profile, Size};

/// Peak resident set size in KiB, from `/proc/self/status` (Linux only;
/// returns 0 elsewhere, so the column is informational, never asserted on).
fn peak_rss_kib() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse().ok())
        })
        .unwrap_or(0)
}

fn sweep(var: &str, default: &[u64]) -> Vec<u64> {
    match env::var(var) {
        Ok(v) => v
            .split(',')
            .map(|t| t.trim().parse().expect("sweep values must be integers"))
            .collect(),
        Err(_) => default.to_vec(),
    }
}

/// Worker count for both the rayon pool and the bgzf writer. Reads
/// `VCFIXTURE_BENCH_WORKERS` in the same style as the other sweep knobs so
/// scaling can be measured with `RAYON_NUM_THREADS` alone being insufficient
/// (the harness passes this value to `.workers(...)` directly, independent
/// of rayon's global thread pool).
fn workers() -> NonZeroUsize {
    match env::var("VCFIXTURE_BENCH_WORKERS") {
        Ok(v) => NonZeroUsize::new(
            v.trim()
                .parse()
                .expect("workers must be a positive integer"),
        )
        .expect("workers must be nonzero"),
        Err(_) => std::thread::available_parallelism().expect("parallelism"),
    }
}

fn main() {
    let dir = env::temp_dir().join("vcfixture_bulk_bench");
    std::fs::create_dir_all(&dir).expect("create bench output dir");

    let samples = sweep("VCFIXTURE_BENCH_SAMPLES", &[500, 2_000, 8_000]);
    let records = sweep("VCFIXTURE_BENCH_RECORDS", &[5_000, 20_000]);

    let workers = workers();
    println!("workers={workers}");
    println!(
        "{:>8} {:>9} {:>12} {:>10} {:>12} {:>10}",
        "samples", "records", "cells", "secs", "s/cell", "peakRSS_MB"
    );

    for &n_samples in &samples {
        for &n_records in &records {
            let profile = Profile::builtin("germline-1kgp").expect("built-in profile loads");
            let path = dir.join(format!("bench_{n_samples}_{n_records}.bcf"));

            let t0 = Instant::now();
            let summary = BulkSpec::new(profile)
                .samples(n_samples as usize)
                .contigs(["chr1", "chr2"])
                .size(Size::Records(n_records))
                .payload(Payload::GtOnly)
                .seed(42)
                .workers(workers)
                .write(&path)
                .expect("bulk generation succeeds");
            let secs = t0.elapsed().as_secs_f64();

            // Ploidy 2 is the germline-1kgp profile's dialed value; the
            // cell count is what cost is linear in.
            let cells = summary.n_records_total() * n_samples * 2;
            println!(
                "{n_samples:>8} {n_records:>9} {cells:>12} {secs:>10.3} {:>12.3e} {:>10.1}",
                secs / cells as f64,
                peak_rss_kib() as f64 / 1024.0,
            );

            let _ = std::fs::remove_file(&path);
            let _ = std::fs::remove_file(path.with_extension("bcf.csi"));
            let _ = std::fs::remove_file(format!("{}.summary.json", path.display()));
        }
    }
}
