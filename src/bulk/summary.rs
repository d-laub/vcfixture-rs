use std::collections::BTreeMap;

use crate::bulk::sample::VariantClass;
use crate::bulk::BulkError;

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x100_0000_01b3;

/// Per-contig record count and observed position range.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ContigSummary {
    pub n_records: u64,
    pub pos_min: u64,
    pub pos_max: u64,
}

/// Counts and an order-sensitive checksum folded out of a bulk record stream.
///
/// Unlike [`crate::truth::GroundTruth`], this is not a per-genotype oracle —
/// it is cheap-to-derive summary truth that a consumer checks a read-back
/// against: same counts, same ranges, same [`Summary::genotype_checksum`].
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Summary {
    pub n_samples: usize,
    pub per_contig: BTreeMap<String, ContigSummary>,
    pub n_alleles_total: u64,
    pub n_alleles_nonref: u64,
    pub class_counts: BTreeMap<String, u64>,
    pub genotype_checksum: u64,
}

fn class_name(c: VariantClass) -> &'static str {
    match c {
        VariantClass::Snp => "snp",
        VariantClass::Insertion => "insertion",
        VariantClass::Deletion => "deletion",
        VariantClass::Mnp => "mnp",
        VariantClass::Complex => "complex",
        VariantClass::Symbolic => "symbolic",
    }
}

/// Number of [`VariantClass`] variants. A fixed-size array indexed by class
/// replaces a `BTreeMap<String, u64>` probed by name once per record, so a
/// block's class counts cost no allocation and merge by elementwise add.
pub const N_VARIANT_CLASSES: usize = 6;

fn class_index(c: VariantClass) -> usize {
    match c {
        VariantClass::Snp => 0,
        VariantClass::Insertion => 1,
        VariantClass::Deletion => 2,
        VariantClass::Mnp => 3,
        VariantClass::Complex => 4,
        VariantClass::Symbolic => 5,
    }
}

fn class_by_index(i: usize) -> VariantClass {
    match i {
        0 => VariantClass::Snp,
        1 => VariantClass::Insertion,
        2 => VariantClass::Deletion,
        3 => VariantClass::Mnp,
        4 => VariantClass::Complex,
        5 => VariantClass::Symbolic,
        _ => unreachable!("class index {i} is out of range 0..{N_VARIANT_CLASSES}"),
    }
}

/// One block's folded contribution to a [`Summary`].
///
/// Accumulated entirely inside a rayon block task (so the `O(n_samples)`
/// per-record fold runs in the fan-out, not on the serial writer thread),
/// then merged in O(1) by [`Summary::merge_block`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockSummary {
    pub n_records: u64,
    pub pos_min: u64,
    pub pos_max: u64,
    pub n_alleles_total: u64,
    pub n_alleles_nonref: u64,
    pub class_counts: [u64; N_VARIANT_CLASSES],
    /// FNV-1a over this block's allele bytes, in record-then-slot order.
    pub checksum: u64,
}

impl Default for BlockSummary {
    fn default() -> BlockSummary {
        BlockSummary::new()
    }
}

impl BlockSummary {
    pub fn new() -> BlockSummary {
        BlockSummary {
            n_records: 0,
            // `pos_min` starts at the maximum so the first `min` wins; an
            // empty block is never merged (see `Summary::merge_block`), so
            // this sentinel can never reach a `Summary`.
            pos_min: u64::MAX,
            pos_max: 0,
            n_alleles_total: 0,
            n_alleles_nonref: 0,
            class_counts: [0; N_VARIANT_CLASSES],
            checksum: FNV_OFFSET,
        }
    }

    /// Fold one record's genotypes into this block. `O(gts.len())`, no
    /// allocation.
    pub fn observe(&mut self, pos: u64, class: VariantClass, gts: &[i8]) {
        self.n_records += 1;
        self.pos_min = self.pos_min.min(pos);
        self.pos_max = self.pos_max.max(pos);
        self.class_counts[class_index(class)] += 1;
        self.n_alleles_total += gts.len() as u64;
        for &g in gts {
            if g > 0 {
                self.n_alleles_nonref += 1;
            }
            self.checksum ^= g as u8 as u64;
            self.checksum = self.checksum.wrapping_mul(FNV_PRIME);
        }
    }
}

impl Summary {
    pub fn new(n_samples: usize) -> Summary {
        Summary {
            n_samples,
            per_contig: BTreeMap::new(),
            n_alleles_total: 0,
            n_alleles_nonref: 0,
            class_counts: BTreeMap::new(),
            genotype_checksum: FNV_OFFSET,
        }
    }

    /// Merge one block's fold into the running summary, in block order.
    ///
    /// O(1) in the block's record and allele counts: the per-allele work
    /// already happened in [`BlockSummary::observe`] inside the worker.
    ///
    /// `genotype_checksum` is hierarchical — FNV-1a over each block's own
    /// checksum bytes, in the order blocks are merged. It therefore stays
    /// order-sensitive both within a block and across blocks, while being
    /// independent of worker count and of how many blocks are in flight
    /// (block boundaries depend only on record count and cohort width,
    /// never on `workers`).
    pub fn merge_block(&mut self, chrom: &str, b: &BlockSummary) {
        // An empty block must not perturb the checksum: folding its
        // untouched `FNV_OFFSET` would make the result depend on how many
        // empty blocks happened to exist.
        if b.n_records == 0 {
            return;
        }

        match self.per_contig.get_mut(chrom) {
            Some(e) => {
                e.n_records += b.n_records;
                e.pos_min = e.pos_min.min(b.pos_min);
                e.pos_max = e.pos_max.max(b.pos_max);
            }
            None => {
                self.per_contig.insert(
                    chrom.to_string(),
                    ContigSummary {
                        n_records: b.n_records,
                        pos_min: b.pos_min,
                        pos_max: b.pos_max,
                    },
                );
            }
        }

        for (i, &n) in b.class_counts.iter().enumerate() {
            if n == 0 {
                continue;
            }
            let name = class_name(class_by_index(i));
            match self.class_counts.get_mut(name) {
                Some(count) => *count += n,
                None => {
                    self.class_counts.insert(name.to_string(), n);
                }
            }
        }

        self.n_alleles_total += b.n_alleles_total;
        self.n_alleles_nonref += b.n_alleles_nonref;

        for byte in b.checksum.to_le_bytes() {
            self.genotype_checksum ^= byte as u64;
            self.genotype_checksum = self.genotype_checksum.wrapping_mul(FNV_PRIME);
        }
    }

    /// Temporary shim over [`Summary::merge_block`], kept only until the
    /// block pipeline lands (issue #22 task 6) and rewires the still-serial
    /// call sites in `src/bulk/mod.rs` to build and merge real
    /// [`BlockSummary`] blocks directly. Constructs a one-record block and
    /// merges it, so callers are unaffected in shape but pay the
    /// merge-per-record cost this task exists to move off the writer
    /// thread. `merge_block` is the real entry point; this wrapper is
    /// deleted once nothing calls it.
    pub fn observe(&mut self, chrom: &str, pos: u64, class: VariantClass, gts: &[i8]) {
        let mut b = BlockSummary::new();
        b.observe(pos, class, gts);
        self.merge_block(chrom, &b);
    }

    pub fn n_records_total(&self) -> u64 {
        self.per_contig.values().map(|c| c.n_records).sum()
    }

    pub fn to_json(&self) -> Result<String, BulkError> {
        Ok(serde_json::to_string_pretty(self)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bulk::sample::VariantClass;

    fn block(records: &[(u64, VariantClass, &[i8])]) -> BlockSummary {
        let mut b = BlockSummary::new();
        for (pos, class, gts) in records {
            b.observe(*pos, *class, gts);
        }
        b
    }

    #[test]
    fn tracks_counts_and_ranges() {
        let mut s = Summary::new(2);
        s.merge_block(
            "chr1",
            &block(&[
                (100, VariantClass::Snp, &[0, 1, 1, 1]),
                (500, VariantClass::Deletion, &[0, 0, 0, 0]),
            ]),
        );
        s.merge_block("chr2", &block(&[(7, VariantClass::Snp, &[1, 1, 0, 0])]));

        assert_eq!(s.n_records_total(), 3);
        assert_eq!(s.per_contig["chr1"].n_records, 2);
        assert_eq!(s.per_contig["chr1"].pos_min, 100);
        assert_eq!(s.per_contig["chr1"].pos_max, 500);
        assert_eq!(s.per_contig["chr2"].pos_min, 7);
        assert_eq!(s.n_alleles_total, 12);
        assert_eq!(s.n_alleles_nonref, 5);
        assert_eq!(s.class_counts["snp"], 2);
        assert_eq!(s.class_counts["deletion"], 1);
    }

    #[test]
    fn missing_alleles_are_not_counted_as_nonref() {
        let mut s = Summary::new(1);
        s.merge_block("chr1", &block(&[(1, VariantClass::Snp, &[-1, 1])]));
        assert_eq!(s.n_alleles_nonref, 1);
    }

    #[test]
    fn checksum_detects_a_dropped_record() {
        let mut a = Summary::new(1);
        a.merge_block(
            "chr1",
            &block(&[
                (1, VariantClass::Snp, &[0, 1]),
                (2, VariantClass::Snp, &[1, 1]),
            ]),
        );
        let mut b = Summary::new(1);
        b.merge_block("chr1", &block(&[(1, VariantClass::Snp, &[0, 1])]));
        assert_ne!(a.genotype_checksum, b.genotype_checksum);
    }

    #[test]
    fn checksum_detects_reordering_within_a_block() {
        let mut a = Summary::new(1);
        a.merge_block(
            "chr1",
            &block(&[
                (1, VariantClass::Snp, &[0, 1]),
                (2, VariantClass::Snp, &[1, 0]),
            ]),
        );
        let mut b = Summary::new(1);
        b.merge_block(
            "chr1",
            &block(&[
                (1, VariantClass::Snp, &[1, 0]),
                (2, VariantClass::Snp, &[0, 1]),
            ]),
        );
        assert_ne!(a.genotype_checksum, b.genotype_checksum);
    }

    #[test]
    fn checksum_detects_reordering_of_blocks() {
        let b0 = block(&[(1, VariantClass::Snp, &[0, 1])]);
        let b1 = block(&[(2, VariantClass::Snp, &[1, 1])]);

        let mut a = Summary::new(1);
        a.merge_block("chr1", &b0);
        a.merge_block("chr1", &b1);

        let mut b = Summary::new(1);
        b.merge_block("chr1", &b1);
        b.merge_block("chr1", &b0);

        assert_ne!(a.genotype_checksum, b.genotype_checksum);
    }

    #[test]
    fn block_decomposition_does_not_change_the_summary_except_the_checksum() {
        // Same records, one block vs two: every count and range must agree.
        // (The checksum is deliberately block-structured and is expected to
        // differ — that is what `checksum_detects_reordering_of_blocks` pins.)
        let recs: [(u64, VariantClass, &[i8]); 4] = [
            (1, VariantClass::Snp, &[0, 1]),
            (2, VariantClass::Deletion, &[1, 1]),
            (3, VariantClass::Snp, &[-1, -1]),
            (4, VariantClass::Mnp, &[0, 0]),
        ];

        let mut one = Summary::new(1);
        one.merge_block("chr1", &block(&recs));

        let mut two = Summary::new(1);
        two.merge_block("chr1", &block(&recs[..2]));
        two.merge_block("chr1", &block(&recs[2..]));

        assert_eq!(one.per_contig, two.per_contig);
        assert_eq!(one.n_alleles_total, two.n_alleles_total);
        assert_eq!(one.n_alleles_nonref, two.n_alleles_nonref);
        assert_eq!(one.class_counts, two.class_counts);
    }

    #[test]
    fn merging_an_empty_block_is_a_no_op() {
        let mut a = Summary::new(1);
        a.merge_block("chr1", &block(&[(1, VariantClass::Snp, &[0, 1])]));
        let before = a.clone();
        a.merge_block("chr1", &BlockSummary::new());
        assert_eq!(a, before);
    }

    #[test]
    fn serializes_to_json() {
        let mut s = Summary::new(1);
        s.merge_block("chr1", &block(&[(1, VariantClass::Snp, &[0, 1])]));
        let j = s.to_json().unwrap();
        assert!(j.contains("\"n_samples\""));
        assert!(j.contains("\"genotype_checksum\""));
    }
}
