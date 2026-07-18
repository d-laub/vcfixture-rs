//! Validate a profile JSON file through the same `Profile::validate` the
//! embedded profiles pass, so a freshly-fitted profile fails here in CI
//! rather than later at `include_str!` time.
use std::process::ExitCode;
use vcfixture::bulk::Profile;

fn main() -> ExitCode {
    let path = match std::env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("usage: validate-profile <path.json>");
            return ExitCode::FAILURE;
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("{path}: read error: {e}");
            return ExitCode::FAILURE;
        }
    };
    match Profile::from_json(&text).and_then(|p| p.validate().map(|_| p)) {
        Ok(p) => {
            println!("{path}: OK ({} contigs)", p.fitted.contigs.len());
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("{path}: INVALID: {e}");
            ExitCode::FAILURE
        }
    }
}
