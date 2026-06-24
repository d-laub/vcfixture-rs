# vcfixture

Generate small VCF test data with a decoded ground-truth oracle, for property-testing VCF parsers.

## Overview

`vcfixture` lets you build a [`Document`] in code, render it to a VCF string (or file), and get
back a [`GroundTruth`] with numpy-style arrays of positions, genotypes, and per-allele metadata —
no hand-coded expected arrays needed.

## Example

```rust
use vcfixture::{Allele, RecordSpec, VcfBuilder, FieldValue};
use vcfixture::spec::version::LATEST;

let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
    .info("AF", None, None, None).unwrap()
    .format("GT", None, None, None).unwrap()
    .record(
        RecordSpec::at("chr1", 1000)
            .ref_("A")
            .alt([Allele::seq("T").unwrap()])
            .gt(["0|1", "1|1"])
            .info("AF", FieldValue::floats([0.25])),
    ).unwrap()
    .build().unwrap();

let truth = doc.truth();
assert_eq!(truth.genotypes[[0, 0, 1]], 1);
assert_eq!(truth.pos[0], 1000);
let _text = doc.render();
```

## Proptest strategies

Hypothesis-style strategies for fuzzing a VCF parser are available behind the `proptest` feature:

```toml
[dev-dependencies]
vcfixture = { version = "0.1", features = ["proptest"] }
```

## Design

See the [design spec](docs/superpowers/specs/2026-06-23-vcfixture-rs-design.md) for the full
architecture and implementation plan.
