# Fields and values

Declare INFO/FORMAT fields as reserved (looked up in the spec registry), typed
(you choose `Number` and `Type`), or flag (INFO-only). Values are built with
`FieldValue`.

```rust
{{#rustdoc_include ../../../examples/fields.rs:fields}}
```
