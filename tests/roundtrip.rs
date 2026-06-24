#![cfg(feature = "proptest")]

use proptest::prelude::*;
use vcfixture::strategies::{
    documents, documents_with_fields, reference_and_documents, symbolic_documents, DocumentOpts,
};

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn render_then_parse_genotype_counts_match_truth(doc in documents(DocumentOpts::default())) {
        let truth = doc.truth();
        let text = doc.render();

        // Independent re-read: count data lines and assert the genotype matrix
        // record count matches.
        let data_lines = text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).count();
        prop_assert_eq!(data_lines, truth.pos.len());

        // Genotype dims are consistent.
        prop_assert_eq!(truth.genotypes.shape()[0], truth.pos.len());
        prop_assert_eq!(truth.genotypes.shape()[1], truth.samples.len());
    }

    #[test]
    fn noodles_reparse_gt_matches_truth(doc in documents(DocumentOpts::default())) {
        use noodles_vcf as vcf;
        use vcf::variant::record_buf::samples::sample::value::genotype::Genotype as NoodlesGenotype;

        let truth = doc.truth();
        let text = doc.render();

        let mut reader = vcf::io::Reader::new(text.as_bytes());
        let header = reader.read_header()
            .expect("noodles header parse");

        let n_samples = truth.samples.len();

        for (ri, result) in reader.record_bufs(&header).enumerate() {
            let rec_buf = result.expect("noodles record parse");
            let samples = rec_buf.samples();

            // Verify GT series.
            let gt_series = samples.select("GT").expect("GT column present");
            for si in 0..n_samples {
                let raw_value = gt_series.get(si).expect("sample index in-bounds");
                let Some(val) = raw_value else {
                    // Missing sample: truth should also be all -1.
                    let ploidy = truth.genotypes.shape()[2];
                    for ai in 0..ploidy {
                        prop_assert_eq!(truth.genotypes[[ri, si, ai]], -1);
                    }
                    continue;
                };
                use vcf::variant::record_buf::samples::sample::Value;
                let noodles_gt: NoodlesGenotype = match val {
                    Value::Genotype(g) => g.clone(),
                    Value::String(s) => s.parse::<NoodlesGenotype>()
                        .expect("parse GT string"),
                    other => panic!("unexpected GT value type: {:?}", other),
                };

                let alleles = noodles_gt.as_ref();
                let ploidy = truth.genotypes.shape()[2];
                for ai in 0..ploidy {
                    let truth_idx = truth.genotypes[[ri, si, ai]];
                    if ai < alleles.len() {
                        let noodles_idx = alleles[ai].position().map(|p| p as i32).unwrap_or(-1);
                        prop_assert_eq!(
                            noodles_idx,
                            truth_idx,
                            "GT mismatch at rec={} sample={} allele={}", ri, si, ai
                        );
                    } else {
                        // Extra ploidy slots (due to padding for max ploidy) should be -1.
                        prop_assert_eq!(
                            truth_idx,
                            -1,
                            "Expected padding -1 at rec={} sample={} allele={}", ri, si, ai
                        );
                    }
                }

                // Phasing: truth.phasing[ri, si] == all separators are '|'.
                // noodles: alleles[1..] each carry the phasing of their separator.
                let truth_phased = truth.phasing[[ri, si]];
                // Only check phasing when ploidy > 1 (haploid has no separator).
                if alleles.len() > 1 {
                    use vcf::variant::record::samples::series::value::genotype::Phasing;
                    let noodles_all_phased = alleles.iter().skip(1).all(|a| {
                        a.phasing() == Phasing::Phased
                    });
                    prop_assert_eq!(
                        noodles_all_phased,
                        truth_phased,
                        "Phasing mismatch at rec={} sample={}", ri, si
                    );
                }
            }
        }
    }

    /// `documents_with_fields` must be valid-by-construction: build/truth/render
    /// never panic and the rendered data-line count equals truth.pos.len().
    #[test]
    fn documents_with_fields_valid(doc in documents_with_fields(DocumentOpts::default())) {
        let truth = doc.truth();
        let text = doc.render();
        let data_lines = text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).count();
        prop_assert_eq!(data_lines, truth.pos.len());
        // The curated extra fields should be present in the rendered header.
        prop_assert!(text.contains("##INFO=<ID=AF,"));
        prop_assert!(text.contains("##FORMAT=<ID=GQ,"));
    }

    /// `reference_and_documents` must be valid-by-construction and the returned
    /// truth must match the document.
    #[test]
    fn reference_and_documents_valid(triple in reference_and_documents(DocumentOpts::default())) {
        let (ref_spec, doc, truth) = triple;
        prop_assert!(!ref_spec.contigs.is_empty());
        let text = doc.render();
        let data_lines = text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).count();
        prop_assert_eq!(data_lines, truth.pos.len());
        prop_assert_eq!(truth.genotypes.shape()[0], truth.pos.len());

        // Each record's REF must match the reference sequence at its position.
        for (ri, rec) in doc.records.iter().enumerate() {
            let pos0 = (truth.pos[ri] - 1) as usize;
            let ref_len = rec.ref_.len();
            let ref_from_spec = ref_spec
                .seq(&rec.chrom, pos0, ref_len)
                .expect("ref position in bounds");
            prop_assert_eq!(&rec.ref_, &ref_from_spec);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// `symbolic_documents` is the intricate SV path (REF padding / SVLEN /
    /// SVCLAIM). Lock in its valid-by-construction guarantee: build/truth/render
    /// never panic and the structural counts are consistent. Records render with
    /// symbolic ALTs (`<DEL>` etc.), so we only assert structural invariants and
    /// do not re-parse symbolic SVs through the noodles reader.
    #[test]
    fn symbolic_documents_valid(doc in symbolic_documents(DocumentOpts::default())) {
        let truth = doc.truth();
        let text = doc.render();
        let data_lines = text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).count();
        prop_assert_eq!(data_lines, truth.pos.len());
        prop_assert_eq!(truth.genotypes.shape()[0], truth.pos.len());
    }
}
