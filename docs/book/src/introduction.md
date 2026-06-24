# Introduction

`vcfixture` generates small, spec-conformant VCF (v4.x) test data and returns a
decoded **ground-truth oracle** alongside it. Parser tests assert against the
oracle instead of hand-coded expected arrays.

The primary consumer is property-based testing of VCF/SparseVar parsers. You
build a document in code, render it to text (or a bgzipped, indexed file), and
derive a `GroundTruth` of positions, genotypes, and per-allele metadata.

See the [API reference on docs.rs](https://docs.rs/vcfixture) for full type
signatures. Every code block in this guide is taken verbatim from a compiled
example in the crate's `examples/` directory.
