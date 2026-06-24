# Property testing

With the `proptest` feature enabled, `vcfixture::strategies` provides
valid-by-construction `Document` strategies for fuzzing a parser against the
oracle.

```rust
{{#rustdoc_include ../../../examples/proptest_fuzz.rs:proptest}}
```
