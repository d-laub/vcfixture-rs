use vcfixture::spec::version::LATEST;
use vcfixture::{Allele, FieldValue, RecordSpec, VcfBuilder};

#[test]
fn renders_minimal_vcf() {
    let text = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
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
        .unwrap()
        .render();
    insta::assert_snapshot!(text);
}

#[test]
fn writes_bgzipped_file() {
    let dir = std::env::temp_dir().join("vcfixture_rs_write_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("x.vcf.gz");
    let written = VcfBuilder::new(["s1"], [("chr1", Some(1000u64))], LATEST)
        .format("GT")
        .record(
            RecordSpec::at("chr1", 10)
                .ref_("A")
                .alt([Allele::seq("T").unwrap()])
                .gt(["0|1"]),
        )
        .write(&path, vcfixture::WriteOpts::bgzipped_indexed())
        .unwrap();
    assert!(written.exists(), "bgzipped file must exist");
    // Strengthened: the .csi index must exist
    let csi_path = {
        let mut p = written.clone().into_os_string();
        p.push(".csi");
        std::path::PathBuf::from(p)
    };
    assert!(
        csi_path.exists(),
        "CSI index must exist at {}",
        csi_path.display()
    );
    // Round-trip: load the CSI back with noodles_csi
    let _index = noodles_csi::fs::read(&csi_path).unwrap_or_else(|e| {
        panic!(
            "noodles_csi::fs::read failed on {}: {}",
            csi_path.display(),
            e
        )
    });
}
