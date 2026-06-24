//! The crate-wide error type.

use thiserror::Error;

/// All errors produced while declaring fields, adding records, deriving the
/// reserved registry, or writing output.
#[derive(Debug, Error)]
pub enum BuildError {
    #[error("sequence allele bases must be [ACGTN]+, got {0:?}")]
    BadAlleleBases(String),

    #[error("not a valid genotype string: {0:?}")]
    BadGenotype(String),

    #[error("not a valid breakend replacement string: {0:?}")]
    BadBreakend(String),

    #[error("symbolic SV first type must be one of DEL/INS/DUP/INV/CNV, got {0:?}")]
    BadSvType(String),

    #[error("ID {0:?} does not match the VCF key regex")]
    BadFieldId(String),

    #[error("Flag fields must be INFO, not FORMAT")]
    FlagNotInfo,

    #[error("Flag fields must have Number=0")]
    FlagNumberNotZero,

    #[error("Number=0 is only valid for Flag fields")]
    ZeroNumberNotFlag,

    #[error("fixed Number must be >= 0")]
    NegativeFixedNumber,

    #[error("{kind} field {id:?} is not a known reserved field; pass number and type explicitly")]
    UnknownReserved { kind: String, id: String },

    #[error("{kind} field {id:?} was introduced in {since}; not available in {version}")]
    FieldTooNew {
        kind: String,
        id: String,
        since: String,
        version: String,
    },

    #[error("symbolic/breakend ALT requires a single preceding REF padding base, got REF={0:?}")]
    MissingRefPadding(String),

    #[error("SVLEN required for symbolic allele {0}")]
    MissingSvlen(String),

    #[error("SVCLAIM {claim:?} invalid for {allele}; allowed {allowed:?}")]
    BadSvclaim {
        claim: String,
        allele: String,
        allowed: Vec<String>,
    },

    #[error("SVCLAIM required for {0} (D/J/DJ)")]
    SvclaimRequired(String),

    #[error("SVLEN must be missing for {0}")]
    SvlenMustBeMissing(String),

    #[error("{kind} field {id:?} not declared")]
    UndeclaredField { kind: String, id: String },

    #[error("{id} cardinality mismatch: expected {expected}, got {got}")]
    Cardinality {
        id: String,
        expected: usize,
        got: usize,
    },

    #[error("allele index {index} out of range (n_alt={n_alt})")]
    AlleleIndexOutOfRange { index: u32, n_alt: usize },

    #[error("{kind} provides {got} per-sample values but {expected} samples are declared")]
    SampleCountMismatch {
        kind: String,
        expected: usize,
        got: usize,
    },

    #[error("GT not declared; declare it with .format(\"GT\", ...)")]
    GtNotDeclared,

    #[error("FORMAT CN requires equal SVLEN across <CNV>/<DEL>/<DUP> alleles")]
    CnSvlenMismatch,

    #[error("contig {0:?} already added")]
    ContigExists(String),

    #[error("contig {0:?} not found")]
    ContigNotFound(String),

    #[error("range {contig}:{pos0}+{len} runs past contig length {clen}")]
    OutOfBounds {
        contig: String,
        pos0: usize,
        len: usize,
        clen: usize,
    },

    #[error(transparent)]
    Io(#[from] std::io::Error),
}
