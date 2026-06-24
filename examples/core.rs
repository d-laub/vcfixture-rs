//! Core workflow: build a tiny VCF document, read its ground-truth oracle, and
//! render it to text.
//!
//! Run with: `cargo run --example core`

// ANCHOR: build
use vcfixture::spec::version::LATEST;
use vcfixture::{Allele, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    // Two samples, one contig, the latest VCF version.
    let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
        .info("AF") // reserved INFO field: Number=A, Type=Float
        .format("GT") // reserved FORMAT field
        .record(
            RecordSpec::at("chr1", 1000)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|1", "1|1"])
                .info("AF", FieldValue::floats([0.25])),
        )
        .build()
        .expect("document is valid");
    // ANCHOR_END: build

    // ANCHOR: truth
    let truth = doc.truth();

    // `genotypes` is an [records, samples, ploidy] array of allele indices
    // (0 = REF, 1 = first ALT, -1 = missing/padding).
    assert_eq!(truth.genotypes[[0, 0, 0]], 0); // record 0, sample s1, allele 0 = REF
    assert_eq!(truth.genotypes[[0, 0, 1]], 1); // s1 allele 1 = first ALT
    assert_eq!(truth.genotypes[[0, 1, 1]], 1); // s2 allele 1 = first ALT

    // `phasing` is [records, samples]: both genotypes used '|'.
    assert!(truth.phasing[[0, 0]]);
    assert!(truth.phasing[[0, 1]]);

    // `pos` is 1-based, per the VCF spec.
    assert_eq!(truth.pos[0], 1000);
    // ANCHOR_END: truth

    // ANCHOR: render
    let text = doc.render();
    assert!(text.starts_with("##fileformat=VCFv"));
    assert!(text.contains("AF=0.25"));
    assert!(text.contains("0|1"));
    print!("{text}");
    // ANCHOR_END: render
}
