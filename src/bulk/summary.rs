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

    /// Fold one record's genotypes into the running summary.
    ///
    /// Called once per record (~265k times for a benchmark-scale run), so
    /// this must be O(`gts.len()`) with no allocation in the common case.
    /// `BTreeMap::entry` always takes its key by value, which would force a
    /// `String` allocation on *every* call even when the contig/class entry
    /// already exists. Instead we probe with `get_mut` using a borrowed key
    /// (`&str` / `&'static str` — free, since `String: Borrow<str>`) and
    /// only allocate a `String` the first time a given contig or class is
    /// seen; every subsequent call for that key updates the existing entry
    /// in place with zero allocation. The genotype loop itself never
    /// allocates.
    pub fn observe(&mut self, chrom: &str, pos: u64, class: VariantClass, gts: &[i8]) {
        match self.per_contig.get_mut(chrom) {
            Some(e) => {
                e.n_records += 1;
                e.pos_min = e.pos_min.min(pos);
                e.pos_max = e.pos_max.max(pos);
            }
            None => {
                self.per_contig.insert(
                    chrom.to_string(),
                    ContigSummary {
                        n_records: 1,
                        pos_min: pos,
                        pos_max: pos,
                    },
                );
            }
        }

        let cname = class_name(class);
        match self.class_counts.get_mut(cname) {
            Some(count) => *count += 1,
            None => {
                self.class_counts.insert(cname.to_string(), 1);
            }
        }

        self.n_alleles_total += gts.len() as u64;
        for &g in gts {
            if g > 0 {
                self.n_alleles_nonref += 1;
            }
            self.genotype_checksum ^= g as u8 as u64;
            self.genotype_checksum = self.genotype_checksum.wrapping_mul(FNV_PRIME);
        }
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

    #[test]
    fn tracks_counts_and_ranges() {
        let mut s = Summary::new(2);
        s.observe("chr1", 100, VariantClass::Snp, &[0, 1, 1, 1]);
        s.observe("chr1", 500, VariantClass::Deletion, &[0, 0, 0, 0]);
        s.observe("chr2", 7, VariantClass::Snp, &[1, 1, 0, 0]);

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
        s.observe("chr1", 1, VariantClass::Snp, &[-1, 1]);
        assert_eq!(s.n_alleles_nonref, 1);
    }

    #[test]
    fn checksum_detects_a_dropped_record() {
        let mut a = Summary::new(1);
        a.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        a.observe("chr1", 2, VariantClass::Snp, &[1, 1]);
        let mut b = Summary::new(1);
        b.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        assert_ne!(a.genotype_checksum, b.genotype_checksum);
    }

    #[test]
    fn checksum_detects_reordering() {
        let mut a = Summary::new(1);
        a.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        a.observe("chr1", 2, VariantClass::Snp, &[1, 0]);
        let mut b = Summary::new(1);
        b.observe("chr1", 1, VariantClass::Snp, &[1, 0]);
        b.observe("chr1", 2, VariantClass::Snp, &[0, 1]);
        assert_ne!(a.genotype_checksum, b.genotype_checksum);
    }

    #[test]
    fn serializes_to_json() {
        let mut s = Summary::new(1);
        s.observe("chr1", 1, VariantClass::Snp, &[0, 1]);
        let j = s.to_json().unwrap();
        assert!(j.contains("\"n_samples\""));
        assert!(j.contains("\"genotype_checksum\""));
    }
}
