/// A supported VCF spec version. `Ord` is chronological.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VcfVersion {
    V4_1,
    V4_2,
    V4_3,
    V4_4,
    V4_5,
}

/// The latest supported version.
pub const LATEST: VcfVersion = VcfVersion::V4_5;

impl VcfVersion {
    /// The exact `##fileformat` string.
    pub fn as_str(&self) -> &'static str {
        match self {
            VcfVersion::V4_1 => "VCFv4.1",
            VcfVersion::V4_2 => "VCFv4.2",
            VcfVersion::V4_3 => "VCFv4.3",
            VcfVersion::V4_4 => "VCFv4.4",
            VcfVersion::V4_5 => "VCFv4.5",
        }
    }
}

impl std::fmt::Display for VcfVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordered_chronologically() {
        assert!(VcfVersion::V4_3 < VcfVersion::V4_4);
        assert!(VcfVersion::V4_5 >= LATEST);
    }

    #[test]
    fn fileformat_strings() {
        assert_eq!(VcfVersion::V4_2.as_str(), "VCFv4.2");
        assert_eq!(LATEST.as_str(), "VCFv4.5");
    }
}
