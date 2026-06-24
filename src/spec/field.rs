use crate::error::BuildError;
use crate::spec::number::{Number, NumberKind};
use crate::spec::types::Type;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FieldKind {
    Info,
    Format,
}

impl FieldKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            FieldKind::Info => "INFO",
            FieldKind::Format => "FORMAT",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldDef {
    pub id: String,
    pub number: Number,
    pub type_: Type,
    pub description: String,
    pub kind: FieldKind,
}

/// VCF key regex: `[A-Za-z_][0-9A-Za-z_.]*` or the literal `1000G`.
fn valid_id(id: &str) -> bool {
    if id == "1000G" {
        return true;
    }
    let mut chars = id.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() || c == '_' => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.')
}

impl FieldDef {
    pub fn new(
        id: impl Into<String>,
        number: Number,
        type_: Type,
        description: impl Into<String>,
        kind: FieldKind,
    ) -> Result<FieldDef, BuildError> {
        let id = id.into();
        if !valid_id(&id) {
            return Err(BuildError::BadFieldId(id));
        }
        if type_ == Type::Flag {
            if kind != FieldKind::Info {
                return Err(BuildError::FlagNotInfo);
            }
            if number.kind != NumberKind::Flag {
                return Err(BuildError::FlagNumberNotZero);
            }
        } else if number.kind == NumberKind::Flag {
            return Err(BuildError::ZeroNumberNotFlag);
        }
        Ok(FieldDef {
            id,
            number,
            type_,
            description: description.into(),
            kind,
        })
    }

    pub fn header_line(&self) -> String {
        format!(
            "##{}=<ID={},Number={},Type={},Description=\"{}\">",
            self.kind.as_str(),
            self.id,
            self.number.header_str(),
            self.type_.as_str(),
            self.description,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::number::Number;
    use crate::spec::types::Type;

    #[test]
    fn header_line_format() {
        let fd = FieldDef::new(
            "AF",
            Number::A,
            Type::Float,
            "Allele frequency",
            FieldKind::Info,
        )
        .unwrap();
        assert_eq!(
            fd.header_line(),
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
        );
    }

    #[test]
    fn flag_must_be_info_with_number_zero() {
        assert!(matches!(
            FieldDef::new("DB", Number::FLAG, Type::Flag, "x", FieldKind::Format),
            Err(crate::error::BuildError::FlagNotInfo)
        ));
        assert!(matches!(
            FieldDef::new("DB", Number::ONE, Type::Flag, "x", FieldKind::Info),
            Err(crate::error::BuildError::FlagNumberNotZero)
        ));
        assert!(matches!(
            FieldDef::new("X", Number::FLAG, Type::Integer, "x", FieldKind::Info),
            Err(crate::error::BuildError::ZeroNumberNotFlag)
        ));
        assert!(FieldDef::new("DB", Number::FLAG, Type::Flag, "x", FieldKind::Info).is_ok());
    }

    #[test]
    fn bad_id_rejected() {
        assert!(FieldDef::new("1BAD", Number::ONE, Type::Integer, "x", FieldKind::Info).is_err());
        assert!(FieldDef::new("1000G", Number::ONE, Type::Integer, "x", FieldKind::Info).is_ok());
    }
}
