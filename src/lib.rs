//! vcfixture — generate small VCF test data with decoded ground truth.

pub mod allele;
pub mod error;
pub mod genotype;
pub mod spec;
pub mod variants;

pub use allele::{Allele, SvType};
pub use error::BuildError;
pub use genotype::Genotype;
pub use spec::field::{FieldDef, FieldKind};
pub use spec::number::{Number, NumberKind};
pub use spec::types::Type;
pub use spec::version::{VcfVersion, LATEST};
pub use variants::VariantClass;
