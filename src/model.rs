use std::collections::BTreeSet;

use indexmap::IndexMap;

use crate::allele::Allele;
use crate::genotype::Genotype;
use crate::spec::field::FieldDef;
use crate::spec::version::VcfVersion;
use crate::value::FieldValue;

#[derive(Debug, Clone, PartialEq)]
pub struct ContigDef {
    pub id: String,
    pub length: Option<u64>,
}

impl ContigDef {
    pub fn header_line(&self) -> String {
        match self.length {
            Some(n) => format!("##contig=<ID={},length={}>", self.id, n),
            None => format!("##contig=<ID={}>", self.id),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AltDef {
    pub id: String,
    pub description: String,
}

impl AltDef {
    pub fn header_line(&self) -> String {
        format!(
            "##ALT=<ID={},Description=\"{}\">",
            self.id, self.description
        )
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SampleValues {
    pub gt: Option<Genotype>,
    pub values: IndexMap<String, FieldValue>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub chrom: String,
    pub pos: u64,
    pub ids: Option<Vec<String>>,
    pub ref_: String,
    pub alts: Vec<Allele>,
    pub qual: Option<f64>,
    pub filters: Option<Vec<String>>,
    pub info: IndexMap<String, FieldValue>,
    pub fmt_keys: Vec<String>,
    pub samples: Vec<SampleValues>,
    pub labels: BTreeSet<String>,
}

impl Record {
    pub fn n_alt(&self) -> usize {
        self.alts.len()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Document {
    pub version: VcfVersion,
    pub info_defs: Vec<FieldDef>,
    pub format_defs: Vec<FieldDef>,
    pub filter_defs: Vec<(String, String)>,
    pub contigs: Vec<ContigDef>,
    pub samples: Vec<String>,
    pub records: Vec<Record>,
    pub alt_defs: Vec<AltDef>,
}

impl Document {
    pub fn max_ploidy(&self) -> usize {
        let mut p = 1;
        for rec in &self.records {
            for s in &rec.samples {
                if let Some(gt) = &s.gt {
                    p = p.max(gt.ploidy());
                }
            }
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contig_header_lines() {
        assert_eq!(
            ContigDef {
                id: "chr1".into(),
                length: Some(100)
            }
            .header_line(),
            "##contig=<ID=chr1,length=100>"
        );
        assert_eq!(
            ContigDef {
                id: "chr1".into(),
                length: None
            }
            .header_line(),
            "##contig=<ID=chr1>"
        );
    }

    #[test]
    fn max_ploidy_defaults_and_scans() {
        let doc = Document {
            version: crate::spec::version::LATEST,
            info_defs: vec![],
            format_defs: vec![],
            filter_defs: vec![],
            contigs: vec![],
            samples: vec!["s1".into()],
            records: vec![],
            alt_defs: vec![],
        };
        assert_eq!(doc.max_ploidy(), 1);
    }
}
