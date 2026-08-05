#![cfg(feature = "cli")]

#[test]
fn parses_a_size_with_units() {
    // exercised via the public parser fn
    use vcfixture::bulk::parse_size;
    assert_eq!(parse_size("100MB").unwrap(), 100 * 1024 * 1024);
    assert_eq!(parse_size("512KB").unwrap(), 512 * 1024);
    assert_eq!(parse_size("1GB").unwrap(), 1024 * 1024 * 1024);
    assert_eq!(parse_size("2048").unwrap(), 2048);
    let bad = parse_size("banana");
    assert!(
        matches!(&bad, Err(vcfixture::bulk::BulkError::BadSize(s)) if s == "banana"),
        "an unparseable size is an argument error, not an invalid profile: {bad:?}"
    );
    assert!(!bad.unwrap_err().to_string().starts_with("invalid profile:"));
}

/// `--threads 0` is rejected by clap as a usage error before any library
/// code runs, so it can never reach BulkError at all.
#[test]
fn zero_threads_is_a_clap_usage_error() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("unused.bcf");
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_vcfixture"))
        .args(["bulk", "--threads", "0", "-o"])
        .arg(&out_path)
        .output()
        .expect("binary should run");
    assert!(!out.status.success(), "--threads 0 must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("must be >= 1"),
        "clap should explain the constraint, got: {stderr}"
    );
    assert!(
        !stderr.contains("invalid profile"),
        "a bad --threads value must not blame the profile, got: {stderr}"
    );
    assert!(
        !out_path.exists(),
        "nothing should be written for a usage error"
    );
}
