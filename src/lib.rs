//! vcfixture — generate small VCF test data with decoded ground truth.
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

pub mod allele;
pub mod build;
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
