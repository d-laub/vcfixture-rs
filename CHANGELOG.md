## v0.2.0 (2026-07-17)

### Feat

- **bulk**: register germline-1kgp-unphased and somatic-gdc builtins
- **fit**: add sites-only VCF input path to fit_profile.py
- **cli**: add vcfixture bulk subcommand and docs
- **bulk**: add BulkSpec API with span-derived contig lengths
- **fit**: add profile extraction script for pgen sources
- **bulk**: add streaming record generator with block-seeded determinism
- **bulk**: add summary truth with order-sensitive genotype checksum
- **bulk**: add streaming writer with byte counting and second-pass index
- **bulk**: add profile-driven samplers with precomputed CDFs
- **bulk**: add profile schema with fitted/dialed partition

### Fix

- **release**: reference the GH_ACTIONS secret for the release PAT
- **release**: push bump commit via admin PAT to satisfy main ruleset
- **ci**: reclaim workspace ownership after commitizen bump
- **bulk**: bound Size::Target memory, rescale SFS to requested cohort
- **fit**: correct README invariant claim, document 3rd profile, scope test-fit
- **bulk**: place exact AC instead of i.i.d. Bernoulli draw
- **fit**: address code review findings on sites-vcf commit
- **fit**: stream pvar/acount/vmiss stats to bound fit_profile.py memory
- **bulk**: bound BulkSpec memory to one contig and harden its tests
- **fit**: handle multiallelic pvar records and warn on dropped histogram values
- **bulk**: reject infinite ClassMix components in validate
- **bulk**: reject NaN and infinite values in profile validators

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
