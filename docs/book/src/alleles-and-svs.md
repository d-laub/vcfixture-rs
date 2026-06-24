# Alleles and structural variants

ALTs are typed `Allele` values: sequence, symbolic SVs (`<DEL>`, `<INS>`, ...),
breakends, `*`, and `<*>`. Symbolic and breakend alleles require a single-base
REF pad; symbolic alleles require `SVLEN`, and DEL/DUP require `SVCLAIM` at
VCF ≥ 4.4.

```rust
{{#rustdoc_include ../../../examples/symbolic.rs:symbolic}}
```
