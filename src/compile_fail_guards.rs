//! Compile-fail guards for the `#[non_exhaustive]` promise on this crate's
//! public error enums (#19).
//!
//! # Why every variant is listed
//!
//! Each block below matches **every** variant and omits only the wildcard
//! arm. That completeness is the whole point: with all variants covered, the
//! sole remaining reason the `match` can be rejected is `#[non_exhaustive]`
//! itself. A shorter block listing two or three variants would fail to
//! compile whether or not the attribute were present -- it would pass either
//! way and prove nothing. `compile_fail` only asserts *that* compilation
//! failed, so the block has to be built so that just one cause is possible.
//!
//! The blocks are doctests rather than `trybuild` fixtures (which could pin
//! the exact rustc note instead of relying on completeness) because
//! `trybuild` requires rustc 1.88 and this crate's MSRV is 1.86 -- a CI job
//! builds and tests on 1.86, and a dev-dependency would break it.
//!
//! # When you add a variant
//!
//! Add its arm to the matching block below. Until you do, that block fails
//! to compile for the ordinary "you forgot a variant" reason, `compile_fail`
//! keeps passing, and the guard quietly stops testing anything. This is a
//! known limitation of the doctest approach and the reason the note also
//! appears above each enum.
//!
//! Hidden from the rendered docs: it is a test, not documentation. The
//! caller-facing version of this rule is the `# Non-exhaustive` section on
//! [`crate::error::BuildError`] and [`crate::bulk::BulkError`].

/// [`crate::error::BuildError`] must reject an exhaustive downstream `match`.
///
/// ```compile_fail,E0004
/// use vcfixture::error::BuildError;
///
/// fn classify(e: &BuildError) -> u8 {
///     match e {
///         BuildError::BadAlleleBases(_) => 1,
///         BuildError::BadGenotype(_) => 2,
///         BuildError::BadBreakend(_) => 3,
///         BuildError::BadSvType(_) => 4,
///         BuildError::BadFieldId(_) => 5,
///         BuildError::FlagNotInfo => 6,
///         BuildError::FlagNumberNotZero => 7,
///         BuildError::ZeroNumberNotFlag => 8,
///         BuildError::NegativeFixedNumber => 9,
///         BuildError::UnknownReserved { .. } => 10,
///         BuildError::FieldTooNew { .. } => 11,
///         BuildError::MissingRefPadding(_) => 12,
///         BuildError::MissingSvlen(_) => 13,
///         BuildError::BadSvclaim { .. } => 14,
///         BuildError::SvclaimRequired(_) => 15,
///         BuildError::SvlenMustBeMissing(_) => 16,
///         BuildError::UndeclaredField { .. } => 17,
///         BuildError::Cardinality { .. } => 18,
///         BuildError::AlleleIndexOutOfRange { .. } => 19,
///         BuildError::SampleCountMismatch { .. } => 20,
///         BuildError::GtNotDeclared => 21,
///         BuildError::CnSvlenMismatch => 22,
///         BuildError::ContigExists(_) => 23,
///         BuildError::ContigNotFound(_) => 24,
///         BuildError::OutOfBounds { .. } => 25,
///         BuildError::InRecord { .. } => 26,
///         BuildError::Io(_) => 27,
///     }
/// }
/// ```
///
/// The same `match` with a trailing `_ => 0` arm compiles, which is what
/// makes the failure above attributable to the attribute rather than to the
/// arms:
///
/// ```
/// use vcfixture::error::BuildError;
///
/// fn classify(e: &BuildError) -> u8 {
///     match e {
///         BuildError::BadAlleleBases(_) => 1,
///         BuildError::BadGenotype(_) => 2,
///         BuildError::BadBreakend(_) => 3,
///         BuildError::BadSvType(_) => 4,
///         BuildError::BadFieldId(_) => 5,
///         BuildError::FlagNotInfo => 6,
///         BuildError::FlagNumberNotZero => 7,
///         BuildError::ZeroNumberNotFlag => 8,
///         BuildError::NegativeFixedNumber => 9,
///         BuildError::UnknownReserved { .. } => 10,
///         BuildError::FieldTooNew { .. } => 11,
///         BuildError::MissingRefPadding(_) => 12,
///         BuildError::MissingSvlen(_) => 13,
///         BuildError::BadSvclaim { .. } => 14,
///         BuildError::SvclaimRequired(_) => 15,
///         BuildError::SvlenMustBeMissing(_) => 16,
///         BuildError::UndeclaredField { .. } => 17,
///         BuildError::Cardinality { .. } => 18,
///         BuildError::AlleleIndexOutOfRange { .. } => 19,
///         BuildError::SampleCountMismatch { .. } => 20,
///         BuildError::GtNotDeclared => 21,
///         BuildError::CnSvlenMismatch => 22,
///         BuildError::ContigExists(_) => 23,
///         BuildError::ContigNotFound(_) => 24,
///         BuildError::OutOfBounds { .. } => 25,
///         BuildError::InRecord { .. } => 26,
///         BuildError::Io(_) => 27,
///         _ => 0,
///     }
/// }
/// ```
pub mod build_error {}

/// [`crate::bulk::BulkError`] must reject an exhaustive downstream `match`.
///
/// ```compile_fail,E0004
/// use vcfixture::bulk::BulkError;
///
/// fn exit_code(e: &BulkError) -> i32 {
///     match e {
///         BulkError::UnknownProfile(_) => 1,
///         BulkError::InvalidProfile(_) => 2,
///         BulkError::PayloadPloidy { .. } => 3,
///         BulkError::NoContigs => 4,
///         BulkError::NoSamples => 5,
///         BulkError::DuplicateContig(_) => 6,
///         BulkError::PerContigMissing(_) => 7,
///         BulkError::PerContigUnknown(_) => 8,
///         BulkError::BadSize(_) => 9,
///         BulkError::BadRecordsFor(_) => 10,
///         BulkError::CompressionLevel(_) => 11,
///         BulkError::ProfileLoad { .. } => 12,
///         BulkError::WorkerPool(_) => 13,
///         BulkError::TargetNotReached { .. } => 14,
///         BulkError::Json(_) => 15,
///         BulkError::Io(_) => 16,
///     }
/// }
/// ```
///
/// The same `match` with a trailing `_ => 0` arm compiles:
///
/// ```
/// use vcfixture::bulk::BulkError;
///
/// fn exit_code(e: &BulkError) -> i32 {
///     match e {
///         BulkError::UnknownProfile(_) => 1,
///         BulkError::InvalidProfile(_) => 2,
///         BulkError::PayloadPloidy { .. } => 3,
///         BulkError::NoContigs => 4,
///         BulkError::NoSamples => 5,
///         BulkError::DuplicateContig(_) => 6,
///         BulkError::PerContigMissing(_) => 7,
///         BulkError::PerContigUnknown(_) => 8,
///         BulkError::BadSize(_) => 9,
///         BulkError::BadRecordsFor(_) => 10,
///         BulkError::CompressionLevel(_) => 11,
///         BulkError::ProfileLoad { .. } => 12,
///         BulkError::WorkerPool(_) => 13,
///         BulkError::TargetNotReached { .. } => 14,
///         BulkError::Json(_) => 15,
///         BulkError::Io(_) => 16,
///         _ => 0,
///     }
/// }
/// ```
#[cfg(feature = "bulk")]
pub mod bulk_error {}
