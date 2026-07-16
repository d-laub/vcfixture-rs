//! Bulk generation of realistic-enough VCF/BCF at benchmark scale.
//!
//! Unlike the fixture path ([`crate::build::VcfBuilder`]), bulk generation
//! streams records and derives no per-genotype oracle — see
//! `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`.

pub mod profile;
pub mod sample;
pub mod summary;
pub mod writer;

pub use profile::{Payload, Profile};

/// Errors from bulk generation.
#[derive(Debug, thiserror::Error)]
pub enum BulkError {
    #[error("unknown builtin profile: {0}")]
    UnknownProfile(String),
    #[error("invalid profile: {0}")]
    Invalid(String),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}
