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
//!
//! `VCFIXTURE_BENCH_REPS=3` repeats each sweep cell and reports min/median/max
//! seconds; `s/cell` is computed from the min. `VCFIXTURE_BENCH_FORMAT=vcf`
//! writes uncompressed VCF, which removes the bgzf compression pool entirely —
//! the ablation that separates rayon scaling from writer-pool contention.

// See the `#[global_allocator]` comment in `src/bin/vcfixture.rs`.
//
// This is on by default, so a plain `--features bulk` build of this harness
// measures *mimalloc*. To measure the glibc baseline, build with
// `--no-default-features --features bulk`. Reporting a number without
// stating which of the two produced it is how an allocator comparison
// becomes meaningless.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::env;
use std::num::NonZeroUsize;
use std::time::Instant;

use vcfixture::bulk::{BulkSpec, Format, Payload, Profile, Size};

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

/// Repetitions per sweep cell. Benchmarks on this box show ~10% run-to-run
/// spread, so a single shot cannot distinguish a real 5% effect from noise.
/// Reported as min/median/max; `s/cell` uses the min, the standard robust
/// estimator for "how fast can this machine do it" (noise only ever adds
/// time).
fn reps() -> usize {
    match env::var("VCFIXTURE_BENCH_REPS") {
        Ok(v) => {
            let n: usize = v.trim().parse().expect("reps must be a positive integer");
            assert!(n > 0, "reps must be a positive integer");
            n
        }
        Err(_) => 1,
    }
}

/// Output format, with the file extension to use for the bench output path.
/// `Vcf` is the ablation that matters: it is the only variant with no bgzf
/// compression pool, so it measures rayon scaling with nothing else competing
/// for cores.
fn format() -> (Format, &'static str) {
    match env::var("VCFIXTURE_BENCH_FORMAT").as_deref() {
        Ok("bcf") | Err(_) => (Format::Bcf, "bcf"),
        Ok("vcf") => (Format::Vcf, "vcf"),
        Ok("vcf.gz") => (Format::VcfGz, "vcf.gz"),
        Ok(other) => panic!("VCFIXTURE_BENCH_FORMAT must be bcf, vcf, or vcf.gz; got {other:?}"),
    }
}

fn main() {
    let dir = env::temp_dir().join("vcfixture_bulk_bench");
    std::fs::create_dir_all(&dir).expect("create bench output dir");

    let samples = sweep("VCFIXTURE_BENCH_SAMPLES", &[500, 2_000, 8_000]);
    let records = sweep("VCFIXTURE_BENCH_RECORDS", &[5_000, 20_000]);

    let workers = workers();
    let reps = reps();
    let (format, ext) = format();
    // Stamp the allocator into the output. Sweeps from the two builds are
    // otherwise indistinguishable once the numbers are pasted into a
    // document, and an allocator comparison that cannot say which binary
    // produced which column is not a comparison.
    let alloc = if cfg!(feature = "mimalloc") {
        "mimalloc"
    } else {
        "system"
    };
    println!("workers={workers} reps={reps} format={ext} alloc={alloc}");
    println!(
        "{:>8} {:>9} {:>12} {:>5} {:>9} {:>9} {:>9} {:>12} {:>10}",
        "samples", "records", "cells", "reps", "min_s", "med_s", "max_s", "s/cell", "peakRSS_MB"
    );

    for &n_samples in &samples {
        for &n_records in &records {
            let path = dir.join(format!("bench_{n_samples}_{n_records}.{ext}"));

            let mut times: Vec<f64> = Vec::with_capacity(reps);
            let mut n_records_total = 0u64;

            for _ in 0..reps {
                let profile = Profile::builtin("germline-1kgp").expect("built-in profile loads");

                let t0 = Instant::now();
                let summary = BulkSpec::new(profile)
                    .samples(n_samples as usize)
                    .contigs(["chr1", "chr2"])
                    .size(Size::Records(n_records))
                    .payload(Payload::GtOnly)
                    .seed(42)
                    .format(format)
                    .workers(workers)
                    .write(&path)
                    .expect("bulk generation succeeds");
                times.push(t0.elapsed().as_secs_f64());
                n_records_total = summary.n_records_total();

                // Clean between reps: every rep must write a fresh file, and
                // the index/summary siblings must not accumulate. The library
                // appends `.csi` to the full path (only for `Bcf`; see
                // `BulkWriter::finish_and_index`), not to the stem.
                let _ = std::fs::remove_file(&path);
                let _ = std::fs::remove_file(format!("{}.csi", path.display()));
                let _ = std::fs::remove_file(format!("{}.summary.json", path.display()));
            }

            times.sort_by(|a, b| a.partial_cmp(b).expect("elapsed times are finite"));
            let min = times[0];
            let med = times[times.len() / 2];
            let max = times[times.len() - 1];

            // Ploidy 2 is the germline-1kgp profile's dialed value; the cell
            // count is what cost is linear in.
            let cells = n_records_total * n_samples * 2;
            println!(
                "{n_samples:>8} {n_records:>9} {cells:>12} {reps:>5} {min:>9.3} {med:>9.3} \
                 {max:>9.3} {:>12.3e} {:>10.1}",
                min / cells as f64,
                peak_rss_kib() as f64 / 1024.0,
            );
        }
    }
}
