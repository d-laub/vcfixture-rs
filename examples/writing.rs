//! Writing a document to a bgzipped, CSI-indexed `.vcf.gz` file.
//!
//! Run with: `cargo run --example writing`

// ANCHOR: writing
use std::env;
use std::fs;

use vcfixture::spec::version::LATEST;
use vcfixture::write::WriteOpts;
use vcfixture::{Allele, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
        .info("AF")
        .format("GT")
        .record(
            RecordSpec::at("chr1", 1000)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|1", "1|1"])
                .info("AF", FieldValue::floats([0.25])),
        )
        .build()
        .expect("document is valid");

    // Render in-memory text whenever you just need a string:
    let _text = doc.render();

    // Or write a bgzipped + CSI-indexed file to disk. The `.gz` extension is
    // ensured for you; `write` returns the final path.
    let dir = env::temp_dir().join("vcfixture_example");
    fs::create_dir_all(&dir).expect("create temp dir");
    let out = doc
        .write(dir.join("out.vcf"), WriteOpts::bgzipped_indexed())
        .expect("write succeeds");

    assert!(out.exists());
    assert_eq!(out.extension().and_then(|e| e.to_str()), Some("gz"));
    // The CSI index sits next to the data file.
    let csi = out.with_extension("gz.csi");
    assert!(csi.exists());
    println!("wrote {} and {}", out.display(), csi.display());
}
// ANCHOR_END: writing
