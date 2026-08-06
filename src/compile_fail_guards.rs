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
//! keeps passing, and the guard quietly stops testing anything.
//!
//! Nothing about a doctest can notice that on its own, so the
//! `variant_lists_are_complete` tests at the bottom of this file do it
//! instead: they are ordinary in-crate exhaustive matches over the same
//! enums. `#[non_exhaustive]` does not restrict matching inside the defining
//! crate, so they compile today and stop compiling the moment a variant is
//! added -- which is the reminder to update the block above. The duplicated
//! arm lists are the price of the guards not being able to go stale in
//! silence.
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

/// In-crate exhaustive matches that fail to compile when a variant is added.
///
/// These do not test `#[non_exhaustive]` -- inside the defining crate the
/// attribute has no effect, so they compile whether or not it is present.
/// Their job is narrower: to make a stale guard block above impossible to
/// miss. Adding a variant breaks these loudly, at the file that says what
/// else needs updating.
#[cfg(test)]
mod variant_lists_are_complete {
    use crate::error::BuildError;

    /// If this stops compiling, a `BuildError` variant was added. Add its arm
    /// here *and* to the `build_error` guard block above -- the guard is only
    /// sound while it covers every variant.
    #[test]
    fn build_error() {
        fn classify(e: &BuildError) -> u8 {
            match e {
                BuildError::BadAlleleBases(_) => 1,
                BuildError::BadGenotype(_) => 2,
                BuildError::BadBreakend(_) => 3,
                BuildError::BadSvType(_) => 4,
                BuildError::BadFieldId(_) => 5,
                BuildError::FlagNotInfo => 6,
                BuildError::FlagNumberNotZero => 7,
                BuildError::ZeroNumberNotFlag => 8,
                BuildError::NegativeFixedNumber => 9,
                BuildError::UnknownReserved { .. } => 10,
                BuildError::FieldTooNew { .. } => 11,
                BuildError::MissingRefPadding(_) => 12,
                BuildError::MissingSvlen(_) => 13,
                BuildError::BadSvclaim { .. } => 14,
                BuildError::SvclaimRequired(_) => 15,
                BuildError::SvlenMustBeMissing(_) => 16,
                BuildError::UndeclaredField { .. } => 17,
                BuildError::Cardinality { .. } => 18,
                BuildError::AlleleIndexOutOfRange { .. } => 19,
                BuildError::SampleCountMismatch { .. } => 20,
                BuildError::GtNotDeclared => 21,
                BuildError::CnSvlenMismatch => 22,
                BuildError::ContigExists(_) => 23,
                BuildError::ContigNotFound(_) => 24,
                BuildError::OutOfBounds { .. } => 25,
                BuildError::InRecord { .. } => 26,
                BuildError::Io(_) => 27,
            }
        }
        assert_eq!(classify(&BuildError::GtNotDeclared), 21);
    }

    /// If this stops compiling, a `BulkError` variant was added. Add its arm
    /// here *and* to the `bulk_error` guard block above.
    #[cfg(feature = "bulk")]
    #[test]
    fn bulk_error() {
        use crate::bulk::BulkError;

        fn exit_code(e: &BulkError) -> i32 {
            match e {
                BulkError::UnknownProfile(_) => 1,
                BulkError::InvalidProfile(_) => 2,
                BulkError::PayloadPloidy { .. } => 3,
                BulkError::NoContigs => 4,
                BulkError::NoSamples => 5,
                BulkError::DuplicateContig(_) => 6,
                BulkError::PerContigMissing(_) => 7,
                BulkError::PerContigUnknown(_) => 8,
                BulkError::BadSize(_) => 9,
                BulkError::BadRecordsFor(_) => 10,
                BulkError::CompressionLevel(_) => 11,
                BulkError::ProfileLoad { .. } => 12,
                BulkError::WorkerPool(_) => 13,
                BulkError::TargetNotReached { .. } => 14,
                BulkError::Json(_) => 15,
                BulkError::Io(_) => 16,
            }
        }
        assert_eq!(exit_code(&BulkError::NoContigs), 4);
    }
}
