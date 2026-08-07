//! Generate small VCF test data with a decoded ground-truth oracle.
//!
//! `vcfixture` builds a VCF [`Document`] in code, renders it to text (or a
//! bgzipped, indexed file), and derives a [`GroundTruth`] — arrays of
//! positions, genotypes, and per-allele metadata — so parser tests assert
//! against a known oracle instead of hand-coded literals.
//!
//! # Workflow
//!
//! 1. [`VcfBuilder`] accumulates samples, contigs, field declarations, and
//!    records. It is infallible until [`VcfBuilder::build`], which validates
//!    everything at once and returns a [`Document`] or a [`BuildError`].
//! 2. [`Document::render`] produces VCF text; [`Document::write`] writes a file;
//!    [`Document::truth`] derives the [`GroundTruth`] oracle.
//!
//! Property-test strategies for fuzzing a parser live in [`strategies`], behind
//! the `proptest` feature (off by default).
//!
//! # Example
//!
//! ```
//! use vcfixture::{Allele, Field, RecordSpec, VcfBuilder, FieldValue};
//! use vcfixture::spec::number::Number;
//! use vcfixture::spec::types::Type;
//! use vcfixture::spec::version::LATEST;
//!
//! let doc = VcfBuilder::new(["s1", "s2"], [("chr1", Some(100_000u64))], LATEST)
//!     .info("AF")
//!     .format("GT")
//!     .format(Field::typed("DS", Number::A, Type::Float))
//!     .record(
//!         RecordSpec::at("chr1", 1000)
//!             .ref_("A")
//!             .alt([Allele::seq("T").unwrap()])
//!             .gt(["0|1", "1|1"])
//!             .info("AF", FieldValue::floats([0.25])),
//!     )
//!     .build().unwrap();
//!
//! let truth = doc.truth();
//! assert_eq!(truth.genotypes[[0, 0, 1]], 1);
//! assert_eq!(truth.pos[0], 1000);
//! let _text = doc.render();
//! ```

#[cfg(feature = "proptest")]
pub mod strategies;

#[cfg(feature = "bulk")]
pub mod bulk;

pub mod allele;
pub mod build;
/// Compile-fail tests for the public error enums' `#[non_exhaustive]`
/// promise. Hidden: a test, not part of the API.
#[doc(hidden)]
pub mod compile_fail_guards;
pub mod error;
pub mod genotype;
pub mod model;
pub mod reference;
pub mod spec;
pub mod truth;
pub mod value;
pub mod variants;
pub mod write;

pub use allele::{Allele, SvType};
pub use build::{Field, RecordSpec, VcfBuilder};
pub use error::BuildError;
pub use genotype::Genotype;
pub use model::{AltDef, ContigDef, Document, Record, SampleValues};
pub use reference::{DrawOpts, ReferenceBuilder, ReferenceSpec, RepeatFeature, VariantKlass};
pub use spec::field::{FieldDef, FieldKind};
pub use spec::number::{Number, NumberKind};
pub use spec::types::Type;
pub use spec::version::{VcfVersion, LATEST};
pub use truth::{AlleleKind, AlleleTruth, GroundTruth};
pub use value::{FieldValue, Scalar};
pub use variants::VariantClass;
pub use write::WriteOpts;
