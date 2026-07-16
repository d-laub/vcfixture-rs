//! Minimal stub for Task 3's `sample` module.
//!
//! This file is intentionally bare — it exists only so `src/bulk/summary.rs`
//! (Task 5) can compile and be tested in isolation while Task 3 is developed
//! in parallel in a sibling worktree. The controller will replace this file
//! with Task 3's fuller version, which will include the fields/derives that
//! module actually needs plus the sampling logic. Do not add anything here.

/// The structural class of a variant record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass {
    Snp,
    Insertion,
    Deletion,
    Mnp,
    Complex,
    Symbolic,
}
