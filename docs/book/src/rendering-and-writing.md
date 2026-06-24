# Rendering and writing

`Document::render()` returns VCF text. `Document::write()` writes a file,
optionally bgzipped and CSI-indexed via `WriteOpts`.

```rust
{{#rustdoc_include ../../../examples/writing.rs:writing}}
```
