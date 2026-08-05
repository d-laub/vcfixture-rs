## v0.4.0 (2026-08-05)

### BREAKING CHANGE

- BulkError also gains BadRecordsFor(String),
PerContigMissing(Vec<String>), and PerContigUnknown(Vec<String>) for the
--records-for and Size::PerContig failures introduced in #17. Their
message text is unchanged; only the misleading "invalid profile: " prefix
is dropped. Downstream `match` arms on BulkError need a wildcard arm.
- BulkError::Invalid is removed. Profile-validation
failures are now BulkError::InvalidProfile (same message, same prefix);
other failures use NoContigs, NoSamples, DuplicateContig, PayloadPloidy,
BadSize, CompressionLevel, ProfileLoad, WorkerPool, or TargetNotReached.

### Feat

- **cli**: add --records-for for explicit per-contig counts
- **bulk**: add Size::PerContig for explicit per-contig counts

### Refactor

- **bulk**: route --records-for and PerContig errors (#16)
- **bulk**: remove the BulkError::Invalid catch-all (#16)
- **bulk**: route writer and CLI errors, parse --threads as NonZero (#16)
- **bulk**: route spec, parsing, and runtime errors (#16)
- **bulk**: route profile-content errors to InvalidProfile (#16)
- **bulk**: add per-class BulkError variants (#16)

## v0.3.0 (2026-07-18)

### Feat

- **bulk**: add provenance.supplied naming non-measured fields (#9)
- **release**: add publish_only dispatch input to recover failed publishes

### Fix

- **fit**: make plink2 memory cap overridable, correct comment (#7)
- **fit**: re-fit somatic-gdc genome-wide (#7)
- **fit**: cap plink2 --memory so it respects the cgroup (#7)
- **fit**: count edges[-1] in single-bin histograms (#10)
- **bulk**: clean up BCF .csi temp on corrective-round drop (#8)
- **bulk**: split records by n_variants, not density_per_kb (#10)
- **bulk**: reject non-diploid AD/PL payloads at validate time (#10)
- **release**: gitignore .release-notes.md so cargo publish sees a clean tree

### Refactor

- **bulk**: move ploidy from fitted to dialed (#9)
- **bulk**: rename module gen -> generate (#6)

### Perf

- **fit**: collect pvar-stat plans sequentially, not collect_all (#7)
- **fit**: compute Ti/Tv by direct base compares, not is_in (#7)
- **fit**: compute gaps by shift+mask, not sort+window (#7)
- **fit**: scan pvar with skip_lines, not comment_prefix (#7)
- **bulk**: calibrate byte target in 2 points, promote the file (#8)
- **bulk**: rejection-sample sparse alt placement (#10)
- **bulk**: build GT into a reused String, not Vec<String>+join (#10)

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
