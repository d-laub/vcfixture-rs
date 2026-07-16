# vcfixture

[![CI](https://github.com/d-laub/vcfixture-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/d-laub/vcfixture-rs/actions/workflows/ci.yml)
[![docs.rs](https://img.shields.io/docsrs/vcfixture)](https://docs.rs/vcfixture)
[![Guide](https://img.shields.io/badge/guide-d--laub.github.io-blue)](https://d-laub.github.io/vcfixture-rs/)

Generate small VCF (v4.x) test data with a decoded ground-truth oracle, for
property-testing VCF parsers. Build a `Document` in code, render it to VCF text
(or a bgzipped, indexed file), and get back a `GroundTruth` with arrays of
positions, genotypes, and per-allele metadata — no hand-coded expected arrays.

```rust
use vcfixture::{Allele, Field, RecordSpec, VcfBuilder, FieldValue};
use vcfixture::spec::number::Number;
use vcfixture::spec::types::Type;
use vcfixture::spec::version::LATEST;

let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
    .info("AF")
    .format("GT")
    .format(Field::typed("DS", Number::A, Type::Float))
    .record(
        RecordSpec::at("chr1", 1000)
            .ref_("A")
            .alt([Allele::seq("T").unwrap()])
            .gt(["0|1", "1|1"])
            .info("AF", FieldValue::floats([0.25])),
    )
    .build().unwrap();

let truth = doc.truth();
assert_eq!(truth.genotypes[[0, 0, 1]], 1);
assert_eq!(truth.pos[0], 1000);
let _text = doc.render();
```

## Examples

Runnable examples live in [`examples/`](examples/):

```bash
cargo run --example core       # build -> truth -> render
cargo run --example fields     # field declarations and typed values
cargo run --example symbolic   # symbolic SVs and breakends
cargo run --example writing    # write a bgzipped, indexed file
cargo run --example proptest_fuzz --features proptest
```

## Proptest strategies

Hypothesis-style strategies for fuzzing a VCF parser are available behind the
`proptest` feature:

```toml
[dev-dependencies]
vcfixture = { version = "0.1", features = ["proptest"] }
```

## Bulk generation

Behind the `cli` feature (which implies `bulk`), the `vcfixture bulk`
subcommand generates large, realistic-enough BCF/VCF files for benchmarking a
reader's speed, memory, or compression — not for exact-value fixtures. It
fits statistics from real data ahead of time into a committed `Profile`, then
streams records from it, so it scales far past what `VcfBuilder` can hold in
memory.

```bash
cargo install vcfixture --features cli
vcfixture bulk --profile germline-1kgp --samples 3202 \
  --contigs chr1,chr2,chr3 --target-size 100MB --seed 42 -o bench.bcf
```

See the [Bulk generation guide](https://d-laub.github.io/vcfixture-rs/bulk-generation.html)
for the fitted-vs-dialed split, the payload presets, and the API.

## Documentation

- [User guide](https://d-laub.github.io/vcfixture-rs/) (mdBook)
- [API reference](https://docs.rs/vcfixture) (docs.rs)
- [Design spec](docs/superpowers/specs/2026-06-23-vcfixture-rs-design.md)
