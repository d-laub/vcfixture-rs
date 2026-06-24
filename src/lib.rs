//! vcfixture — generate small VCF test data with decoded ground truth.

pub mod allele;
pub mod build;
pub mod error;
pub mod genotype;
pub mod model;
pub mod spec;
pub mod truth;
pub mod value;
pub mod variants;

pub use allele::{Allele, SvType};
pub use build::{RecordSpec, VcfBuilder};
pub use error::BuildError;
pub use genotype::Genotype;
pub use model::{AltDef, ContigDef, Document, Record, SampleValues};
pub use spec::field::{FieldDef, FieldKind};
pub use spec::number::{Number, NumberKind};
pub use spec::types::Type;
pub use spec::version::{VcfVersion, LATEST};
pub use truth::{AlleleKind, AlleleTruth, GroundTruth};
pub use value::{FieldValue, Scalar};
pub use variants::VariantClass;
