# Building documents

`VcfBuilder` accumulates samples, contigs, field declarations, and records. It
is infallible until `build()`, which validates everything at once.

```rust
{{#rustdoc_include ../../../examples/core.rs:build}}
```

Validation is deferred: declaration order does not matter, and a single
`BuildError` (tagged with the offending record index) is returned from `build()`
if anything is inconsistent.
