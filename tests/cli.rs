#![cfg(feature = "cli")]

#[test]
fn parses_a_size_with_units() {
    // exercised via the public parser fn
    use vcfixture::bulk::parse_size;
    assert_eq!(parse_size("100MB").unwrap(), 100 * 1024 * 1024);
    assert_eq!(parse_size("512KB").unwrap(), 512 * 1024);
    assert_eq!(parse_size("1GB").unwrap(), 1024 * 1024 * 1024);
    assert_eq!(parse_size("2048").unwrap(), 2048);
    assert!(parse_size("banana").is_err());
}

#[test]
fn parses_records_for_entries() {
    use vcfixture::bulk::parse_records_for;

    let toks = ["chr1=100".to_string(), "chr2=250".to_string()];
    assert_eq!(
        parse_records_for(&toks).unwrap(),
        vec![("chr1".to_string(), 100), ("chr2".to_string(), 250)]
    );

    // Order is load-bearing: with no `--contigs`, these names *are* the
    // output contig order, so the parser must not sort or dedupe them into
    // a map on the caller's behalf.
    let rev = ["chr2=250".to_string(), "chr1=100".to_string()];
    assert_eq!(parse_records_for(&rev).unwrap()[0].0, "chr2");

    // Zero is a legal count ("generate nothing for this contig").
    assert_eq!(parse_records_for(&["chr1=0".to_string()]).unwrap()[0].1, 0);

    // Surrounding whitespace is tolerated.
    assert_eq!(
        parse_records_for(&[" chr1 = 100 ".to_string()]).unwrap(),
        vec![("chr1".to_string(), 100)]
    );
}

#[test]
fn rejects_malformed_records_for_entries() {
    use vcfixture::bulk::parse_records_for;

    let cases = [
        ("chr1", "no '=' separator"),
        ("=100", "empty contig name"),
        ("chr1=", "empty count"),
        ("chr1=banana", "non-numeric count"),
        ("chr1=-5", "negative count"),
    ];
    for (tok, why) in cases {
        assert!(
            parse_records_for(&[tok.to_string()]).is_err(),
            "must reject {tok:?} ({why})"
        );
    }

    // A duplicate name would silently drop one of the two counts once the
    // pairs are collected into a BTreeMap, so reject it at parse time.
    assert!(parse_records_for(&["chr1=1".to_string(), "chr1=2".to_string()]).is_err());
}
