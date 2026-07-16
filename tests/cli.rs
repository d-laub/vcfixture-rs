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
