//! Bulk generation: build a small indexed BCF from the built-in germline
//! profile and confirm the summary matches what was requested.
//!
//! Run with: `cargo run --example bulk --features bulk`

use std::env;
use std::num::NonZero;

// ANCHOR: bulk
use vcfixture::bulk::{BulkSpec, Payload, Profile, Size};

fn main() {
    let dir = env::temp_dir().join("vcfixture_bulk_example");
    std::fs::create_dir_all(&dir).expect("create example output dir");
    let path = dir.join("bench.bcf");

    let profile = Profile::builtin("germline-1kgp").expect("built-in profile loads");

    let summary = BulkSpec::new(profile)
        .samples(8)
        .contigs(["chr1", "chr2", "chr3"])
        .size(Size::RecordsPerContig(100))
        .payload(Payload::GtOnly)
        .seed(42)
        .workers(NonZero::new(2).unwrap())
        .write(&path)
        .expect("bulk generation succeeds");

    assert_eq!(summary.n_records_total(), 300);
    println!(
        "wrote {} records ({} contigs) to {}",
        summary.n_records_total(),
        summary.per_contig.len(),
        path.display()
    );
    // ANCHOR_END: bulk

    // Clean up: this example's whole point is the runtime assertion above,
    // not the artifact it leaves behind.
    let _ = std::fs::remove_file(&path);
    let mut csi = path.clone().into_os_string();
    csi.push(".csi");
    let _ = std::fs::remove_file(csi);
    let mut summary_json = path.into_os_string();
    summary_json.push(".summary.json");
    let _ = std::fs::remove_file(summary_json);
}
