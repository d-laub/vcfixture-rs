## v0.1.0 (2026-06-23)

### Feat

- add Field declaration sub-spec for VcfBuilder
- make Genotype::parse and Allele::parse fallible (error on unsupported input)
- populate fields in documents_with_fields and make reference_and_documents reference-consistent
- add proptest strategies and round-trip tests
- add synthetic reference subsystem
- add VCF rendering, bgzip, and CSI indexing
- add GroundTruth oracle
- add VcfBuilder with eager validation
- add Document/Record model and value types
- add variant classification
- add Genotype
- add Allele model
- add genotype ordering
- add reserved-field registry
- add FieldDef
- add Number and cardinality
- add Type
- add VcfVersion

### Fix

- return SampleCountMismatch instead of panicking on over-long sample vectors
- return OutOfBounds instead of panicking in reference base/seq/set_base

### Refactor

- defer VcfBuilder validation to build()
