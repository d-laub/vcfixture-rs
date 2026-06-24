//! Symbolic structural variants, breakends, and per-allele truth.
//!
//! Symbolic/breakend ALTs require a single-base REF pad. Symbolic alleles
//! require SVLEN; DEL/DUP additionally require SVCLAIM at VCF >= 4.4. Breakends
//! must NOT carry SVLEN.
//!
//! Run with: `cargo run --example symbolic`

// ANCHOR: symbolic
use vcfixture::spec::version::LATEST;
use vcfixture::{Allele, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    let doc = VcfBuilder::new(["s1"], [("chr1", Some(1_000_000u64))], LATEST)
        .info("SVLEN")
        .info("SVCLAIM")
        .format("GT")
        // Symbolic deletion: single-base REF pad, SVLEN required, and at
        // VCF >= 4.4 a DEL also requires SVCLAIM.
        .record(
            RecordSpec::at("chr1", 5000)
                .ref_("A")
                .alt([Allele::deletion(Vec::<&str>::new())])
                .gt(["0|1"])
                .info("SVLEN", FieldValue::ints([100]))
                .info("SVCLAIM", FieldValue::strings(["D"])),
        )
        // Symbolic insertion: SVLEN required, no SVCLAIM requirement.
        .record(
            RecordSpec::at("chr1", 8000)
                .ref_("C")
                .alt([Allele::insertion(Vec::<&str>::new())])
                .gt(["1|1"])
                .info("SVLEN", FieldValue::ints([50])),
        )
        // Paired breakend: the raw replacement string carries the mate locus.
        // Breakends must NOT carry SVLEN.
        .record(
            RecordSpec::at("chr1", 9000)
                .ref_("G")
                .alt([Allele::breakend_parse("G]chr2:321]").unwrap()])
                .gt(["0|1"]),
        )
        .build()
        .expect("document is valid");

    let truth = doc.truth();
    // The deletion's per-allele truth is decoded for you.
    let del = &truth.alts_truth[0][0];
    assert_eq!(del.sv_type.as_deref(), Some("DEL"));
    assert_eq!(del.svlen, Some(100));
    assert!(!del.is_sequence);

    let text = doc.render();
    assert!(text.contains("<DEL>"));
    assert!(text.contains("<INS>"));
    assert!(text.contains("G]chr2:321]"));
    // Symbolic ALT types are auto-described in the header.
    assert!(text.contains("##ALT=<ID=DEL,"));
    print!("{text}");
}
// ANCHOR_END: symbolic
