use crate::allele::{Allele, SvType};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass {
    Snp,
    Mnp,
    Ins,
    Del,
    Delins,
    SpanningDel,
    Unspecified,
    Bnd,
    Multiallelic,
    SvDel,
    SvIns,
    SvDup,
    SvInv,
    Cnv,
}

impl VariantClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            VariantClass::Snp => "SNP",
            VariantClass::Mnp => "MNP",
            VariantClass::Ins => "INS",
            VariantClass::Del => "DEL",
            VariantClass::Delins => "DELINS",
            VariantClass::SpanningDel => "SPANNING_DEL",
            VariantClass::Unspecified => "UNSPECIFIED",
            VariantClass::Bnd => "BND",
            VariantClass::Multiallelic => "MULTIALLELIC",
            VariantClass::SvDel => "SV_DEL",
            VariantClass::SvIns => "SV_INS",
            VariantClass::SvDup => "SV_DUP",
            VariantClass::SvInv => "SV_INV",
            VariantClass::Cnv => "CNV",
        }
    }
}

pub fn classify_seq(ref_: &str, alt: &str) -> VariantClass {
    if alt == "*" {
        return VariantClass::SpanningDel;
    }
    let (lr, la) = (ref_.len(), alt.len());
    if lr == 1 && la == 1 {
        VariantClass::Snp
    } else if lr == la {
        VariantClass::Mnp
    } else if la > lr && alt.starts_with(ref_) {
        VariantClass::Ins
    } else if lr > la && ref_.starts_with(alt) {
        VariantClass::Del
    } else {
        VariantClass::Delins
    }
}

pub fn record_class(ref_: &str, alts: &[Allele]) -> VariantClass {
    if alts.len() != 1 {
        return VariantClass::Multiallelic;
    }
    match &alts[0] {
        Allele::Seq(bases) => classify_seq(ref_, bases),
        Allele::Star => VariantClass::SpanningDel,
        Allele::Unspecified => VariantClass::Unspecified,
        Allele::Breakend { .. } => VariantClass::Bnd,
        Allele::Symbolic { first_type, .. } => match first_type {
            SvType::Del => VariantClass::SvDel,
            SvType::Ins => VariantClass::SvIns,
            SvType::Dup => VariantClass::SvDup,
            SvType::Inv => VariantClass::SvInv,
            SvType::Cnv => VariantClass::Cnv,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;

    #[test]
    fn sequence_classes() {
        assert_eq!(classify_seq("A", "T"), VariantClass::Snp);
        assert_eq!(classify_seq("AC", "GT"), VariantClass::Mnp);
        assert_eq!(classify_seq("A", "AT"), VariantClass::Ins);
        assert_eq!(classify_seq("AT", "A"), VariantClass::Del);
        assert_eq!(classify_seq("AT", "C"), VariantClass::Delins);
        assert_eq!(classify_seq("A", "*"), VariantClass::SpanningDel);
    }

    #[test]
    fn record_classes() {
        assert_eq!(
            record_class("A", &[Allele::seq("T").unwrap()]),
            VariantClass::Snp
        );
        assert_eq!(
            record_class("A", &[Allele::seq("T").unwrap(), Allele::seq("C").unwrap()]),
            VariantClass::Multiallelic
        );
        assert_eq!(
            record_class("A", &[Allele::deletion(Vec::<&str>::new())]),
            VariantClass::SvDel
        );
        assert_eq!(
            record_class("A", &[Allele::cnv(Vec::<&str>::new())]),
            VariantClass::Cnv
        );
        assert_eq!(
            record_class("A", &[Allele::Star]),
            VariantClass::SpanningDel
        );
    }
}
