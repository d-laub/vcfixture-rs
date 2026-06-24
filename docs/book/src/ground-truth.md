# Ground truth

`Document::truth()` derives the oracle. `genotypes` is an
`[records, samples, ploidy]` array of allele indices (`-1` = missing/padding),
`phasing` is `[records, samples]`, and `pos` is 1-based.

```rust
{{#rustdoc_include ../../../examples/core.rs:truth}}
```

INFO and FORMAT are decoded per record (and per sample for FORMAT) into maps of
field id to `FieldValue`; per-allele structural metadata lives in `alts_truth`.
