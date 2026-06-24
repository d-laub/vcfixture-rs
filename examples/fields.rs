//! Declaring INFO/FORMAT fields and constructing typed values.
//!
//! Run with: `cargo run --example fields`

// ANCHOR: fields
use vcfixture::spec::number::Number;
use vcfixture::spec::types::Type;
use vcfixture::spec::version::LATEST;
use vcfixture::{Allele, Field, FieldValue, RecordSpec, VcfBuilder};

fn main() {
    let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
        // Reserved: looked up in the spec registry (AF => Number=A, Type=Float).
        .info("AF")
        // Reserved via the explicit constructor (identical to the &str form).
        .info(Field::reserved("DP"))
        // Typed: you choose Number and Type, plus an optional description.
        .info(Field::typed("AC", Number::A, Type::Integer).description("allele count"))
        // Flag: INFO-only, Number=0, Type=Flag.
        .info(Field::flag("SOMATIC"))
        .format("GT")
        // Per-allele FORMAT field (Number=A).
        .format(Field::typed("DS", Number::A, Type::Float))
        .record(
            RecordSpec::at("chr1", 2000)
                .ref_("G")
                .alt([Allele::seq("C").unwrap()])
                .gt(["0|0", "0|1"])
                .info("AF", FieldValue::floats([0.5])) // list of floats
                .info("DP", FieldValue::ints([42])) // single int (as a 1-elem list)
                .info("AC", FieldValue::ints([1]))
                .info("SOMATIC", FieldValue::Flag) // flag present
                // FORMAT DS: one value per sample.
                .format("DS", [FieldValue::floats([0.4]), FieldValue::floats([1.9])]),
        )
        .build()
        .expect("document is valid");

    // The decoded oracle exposes INFO and FORMAT per record/sample.
    let truth = doc.truth();
    assert_eq!(truth.info[0]["DP"], FieldValue::ints([42]));
    assert_eq!(truth.format[0][1]["DS"], FieldValue::floats([1.9]));

    let text = doc.render();
    // Typed and flag fields render deterministic headers.
    assert!(text.contains("##INFO=<ID=AC,Number=A,Type=Integer,"));
    assert!(text.contains("##INFO=<ID=SOMATIC,Number=0,Type=Flag,"));
    assert!(text.contains("##FORMAT=<ID=DS,Number=A,Type=Float,"));
    // INFO column joins fields with ';' in declaration order.
    assert!(text.contains("AF=0.5;DP=42;AC=1;SOMATIC"));
    print!("{text}");
}
// ANCHOR_END: fields
