use crate::error::BuildError;
use crate::spec::field::{FieldDef, FieldKind};
use crate::spec::number::Number;
use crate::spec::types::Type;
use crate::spec::version::VcfVersion;

fn info_entry(id: &str) -> Option<(Number, Type, &'static str)> {
    Some(match id {
        "AA" => (Number::ONE, Type::String, "Ancestral allele"),
        "AC" => (Number::A, Type::Integer, "Allele count"),
        "AF" => (Number::A, Type::Float, "Allele frequency"),
        "AN" => (Number::ONE, Type::Integer, "Total allele number"),
        "DP" => (Number::ONE, Type::Integer, "Combined depth"),
        "DB" => (Number::FLAG, Type::Flag, "dbSNP membership"),
        "H2" => (Number::FLAG, Type::Flag, "HapMap2 membership"),
        "END" => (Number::ONE, Type::Integer, "End position (deprecated)"),
        "SVTYPE" => (Number::ONE, Type::String, "Type of structural variant"),
        "SVLEN" => (Number::A, Type::Integer, "Length of structural variant"),
        "SVCLAIM" => (Number::A, Type::String, "Structural variant claim"),
        "CIPOS" => (Number::DOT, Type::Integer, "Confidence interval around POS"),
        "CIEND" => (Number::DOT, Type::Integer, "Confidence interval around END"),
        "CILEN" => (
            Number::DOT,
            Type::Integer,
            "Confidence interval around SVLEN",
        ),
        "MATEID" => (Number::A, Type::String, "ID of mate breakend"),
        "PARID" => (Number::A, Type::String, "ID of partner breakend"),
        "IMPRECISE" => (Number::FLAG, Type::Flag, "Imprecise structural variant"),
        _ => return None,
    })
}

fn format_entry(id: &str) -> Option<(Number, Type, &'static str)> {
    Some(match id {
        "GT" => (Number::ONE, Type::String, "Genotype"),
        "GQ" => (Number::ONE, Type::Integer, "Genotype quality"),
        "DP" => (Number::ONE, Type::Integer, "Read depth"),
        "AD" => (Number::R, Type::Integer, "Allelic depths"),
        "PL" => (Number::G, Type::Integer, "Phred genotype likelihoods"),
        "GL" => (Number::G, Type::Float, "Log10 genotype likelihoods"),
        "PS" => (Number::ONE, Type::Integer, "Phase set"),
        "CN" => (Number::ONE, Type::Float, "Copy number"),
        "LEN" => (Number::ONE, Type::Integer, "Length of <*> reference block"),
        _ => return None,
    })
}

/// Version each gated reserved field was introduced; absent ids exist since 4.1.
fn since(id: &str, kind: FieldKind) -> VcfVersion {
    match (kind, id) {
        (FieldKind::Info, "SVCLAIM") => VcfVersion::V4_4,
        (FieldKind::Format, "LEN") => VcfVersion::V4_4,
        _ => VcfVersion::V4_1,
    }
}

pub fn reserved(id: &str, kind: FieldKind, version: VcfVersion) -> Result<FieldDef, BuildError> {
    let entry = match kind {
        FieldKind::Info => info_entry(id),
        FieldKind::Format => format_entry(id),
    };
    let (number, type_, desc) = entry.ok_or_else(|| BuildError::UnknownReserved {
        kind: kind.as_str().to_string(),
        id: id.to_string(),
    })?;

    let intro = since(id, kind);
    if version < intro {
        return Err(BuildError::FieldTooNew {
            kind: kind.as_str().to_string(),
            id: id.to_string(),
            since: intro.to_string(),
            version: version.to_string(),
        });
    }

    // SVLEN's pre-4.4 form: Number=. (signed length difference).
    if id == "SVLEN" && kind == FieldKind::Info && version < VcfVersion::V4_4 {
        return FieldDef::new(
            "SVLEN",
            Number::DOT,
            Type::Integer,
            "Difference in length between REF and ALT alleles",
            FieldKind::Info,
        );
    }

    FieldDef::new(id, number, type_, desc, kind)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::field::FieldKind;
    use crate::spec::number::Number;
    use crate::spec::types::Type;
    use crate::spec::version::VcfVersion;

    #[test]
    fn resolves_af() {
        let fd = reserved("AF", FieldKind::Info, VcfVersion::V4_5).unwrap();
        assert_eq!(fd.number, Number::A);
        assert_eq!(fd.type_, Type::Float);
    }

    #[test]
    fn svlen_form_switches_at_4_4() {
        let pre = reserved("SVLEN", FieldKind::Info, VcfVersion::V4_3).unwrap();
        assert_eq!(pre.number, Number::DOT);
        let post = reserved("SVLEN", FieldKind::Info, VcfVersion::V4_4).unwrap();
        assert_eq!(post.number, Number::A);
    }

    #[test]
    fn svclaim_gated_before_4_4() {
        assert!(reserved("SVCLAIM", FieldKind::Info, VcfVersion::V4_3).is_err());
        assert!(reserved("SVCLAIM", FieldKind::Info, VcfVersion::V4_4).is_ok());
    }

    #[test]
    fn unknown_reserved_errs() {
        assert!(reserved("NOPE", FieldKind::Info, VcfVersion::V4_5).is_err());
    }
}
