//! vcfixture — generate small VCF test data with decoded ground truth.

pub mod error;
pub mod spec;

pub use error::BuildError;
pub use spec::types::Type;
pub use spec::version::{VcfVersion, LATEST};
