//! Fuzzing a VCF parser against the ground-truth oracle.
//!
//! Requires the `proptest` feature:
//!   cargo run --example proptest_fuzz --features proptest
//!
//! In your own test suite you would write this as a `proptest!` block (see the
//! note at the bottom). Here we drive the strategy directly so the example has
//! a runnable `main`.

// ANCHOR: proptest
use proptest::strategy::{Strategy, ValueTree};
use proptest::test_runner::TestRunner;

use vcfixture::strategies::{documents, DocumentOpts};

fn main() {
    let mut runner = TestRunner::default();
    let strategy = documents(DocumentOpts::default());

    for _ in 0..32 {
        // Draw one valid-by-construction document.
        let doc = strategy
            .new_tree(&mut runner)
            .expect("strategy produces a value")
            .current();

        // Derive the oracle and render the document.
        let truth = doc.truth();
        let text = doc.render();

        // A real parser test would parse `text` and compare against `truth`.
        // Here we assert the structural invariant the oracle guarantees: one
        // data line per record in the genotype matrix.
        let data_lines = text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .count();
        assert_eq!(data_lines, truth.pos.len());
        assert_eq!(truth.genotypes.shape()[0], truth.pos.len());
    }

    println!("checked 32 generated documents against their oracle");
}

// In your crate's tests, the idiomatic form is:
//
//   use proptest::prelude::*;
//   use vcfixture::strategies::{documents, DocumentOpts};
//
//   proptest! {
//       #[test]
//       fn parser_matches_oracle(doc in documents(DocumentOpts::default())) {
//           let truth = doc.truth();
//           let text = doc.render();
//           // parse `text` with your parser and prop_assert_eq! against `truth`.
//       }
//   }
// ANCHOR_END: proptest
