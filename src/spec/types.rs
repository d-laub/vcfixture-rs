/// VCF value type for an INFO or FORMAT field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Type {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

impl Type {
    pub fn as_str(&self) -> &'static str {
        match self {
            Type::Integer => "Integer",
            Type::Float => "Float",
            Type::Flag => "Flag",
            Type::Character => "Character",
            Type::String => "String",
        }
    }

    /// All types valid in INFO.
    pub fn info_allowed() -> [Type; 5] {
        [
            Type::Integer,
            Type::Float,
            Type::Flag,
            Type::Character,
            Type::String,
        ]
    }

    /// All types valid in FORMAT (Flag excluded).
    pub fn format_allowed() -> [Type; 4] {
        [Type::Integer, Type::Float, Type::Character, Type::String]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_excludes_flag() {
        assert!(!Type::format_allowed().contains(&Type::Flag));
        assert!(Type::info_allowed().contains(&Type::Flag));
    }

    #[test]
    fn header_tokens() {
        assert_eq!(Type::Float.as_str(), "Float");
    }
}
