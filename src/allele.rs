//! [`Allele`] — sequence, symbolic structural-variant, and breakend ALTs, with
//! constructors and VCF-string rendering/parsing.

use std::str::FromStr;

use crate::error::BuildError;

/// First-level symbolic SV type token (e.g. `DEL`, `DUP`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SvType {
    Del,
    Ins,
    Dup,
    Inv,
    Cnv,
}

impl SvType {
    /// Return the VCF symbolic ID string (e.g. `"DEL"`, `"DUP"`).
    pub fn as_str(&self) -> &'static str {
        match self {
            SvType::Del => "DEL",
            SvType::Ins => "INS",
            SvType::Dup => "DUP",
            SvType::Inv => "INV",
            SvType::Cnv => "CNV",
        }
    }
}

impl FromStr for SvType {
    type Err = BuildError;

    fn from_str(s: &str) -> Result<SvType, BuildError> {
        Ok(match s {
            "DEL" => SvType::Del,
            "INS" => SvType::Ins,
            "DUP" => SvType::Dup,
            "INV" => SvType::Inv,
            "CNV" => SvType::Cnv,
            _ => return Err(BuildError::BadSvType(s.to_string())),
        })
    }
}

/// A VCF ALT allele — sequence, symbolic SV, breakend, or special token.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Allele {
    Seq(String),
    Star,
    Symbolic {
        first_type: SvType,
        subtypes: Vec<String>,
    },
    Unspecified,
    Breakend {
        raw: String,
        single: bool,
    },
}

fn is_seq(s: &str) -> bool {
    !s.is_empty()
        && s.bytes().all(|b| {
            matches!(
                b,
                b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
            )
        })
}

fn is_seq_or_empty(s: &str) -> bool {
    s.bytes().all(|b| {
        matches!(
            b,
            b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
        )
    })
}

impl Allele {
    /// Construct a sequence allele, returning an error if `bases` contains non-IUPAC characters.
    pub fn seq(bases: impl Into<String>) -> Result<Allele, BuildError> {
        let bases = bases.into();
        if !is_seq(&bases) {
            return Err(BuildError::BadAlleleBases(bases));
        }
        Ok(Allele::Seq(bases))
    }

    /// Construct the spanning deletion token `*`.
    pub fn star() -> Allele {
        Allele::Star
    }

    /// Construct the unspecified allele token `<*>`.
    pub fn unspecified() -> Allele {
        Allele::Unspecified
    }

    fn symbolic(first: SvType, subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::Symbolic {
            first_type: first,
            subtypes: subtypes.into_iter().map(Into::into).collect(),
        }
    }

    /// Construct a `<DEL[:subtype]>` symbolic allele.
    pub fn deletion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Del, subtypes)
    }

    /// Construct a `<INS[:subtype]>` symbolic allele.
    pub fn insertion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Ins, subtypes)
    }

    /// Construct a `<DUP[:subtype]>` symbolic allele.
    pub fn duplication(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Dup, subtypes)
    }

    /// Construct a `<INV[:subtype]>` symbolic allele.
    pub fn inversion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Inv, subtypes)
    }

    /// Construct a `<CNV[:subtype]>` symbolic allele.
    pub fn cnv(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Cnv, subtypes)
    }

    /// Parse a breakend replacement string (paired or single forms).
    pub fn breakend_parse(s: &str) -> Result<Allele, BuildError> {
        if is_single_breakend(s) {
            return Ok(Allele::Breakend {
                raw: s.to_string(),
                single: true,
            });
        }
        if is_paired_breakend(s) {
            return Ok(Allele::Breakend {
                raw: s.to_string(),
                single: false,
            });
        }
        Err(BuildError::BadBreakend(s.to_string()))
    }

    /// Inner `<...>` token, e.g. `DEL` or `DUP:TANDEM`.
    pub fn symbolic_type_str(&self) -> Option<String> {
        match self {
            Allele::Symbolic {
                first_type,
                subtypes,
            } => {
                let mut parts = vec![first_type.as_str().to_string()];
                parts.extend(subtypes.iter().cloned());
                Some(parts.join(":"))
            }
            _ => None,
        }
    }

    /// Render the allele to its VCF string representation.
    pub fn render(&self) -> String {
        match self {
            Allele::Seq(b) => b.clone(),
            Allele::Star => "*".to_string(),
            Allele::Unspecified => "<*>".to_string(),
            Allele::Breakend { raw, .. } => raw.clone(),
            Allele::Symbolic { .. } => format!("<{}>", self.symbolic_type_str().unwrap()),
        }
    }

    /// Syntactic dispatch from a raw ALT string.
    ///
    /// Returns `Ok(Allele)` for well-formed inputs and `Err(BuildError)` for:
    /// - unknown symbolic SV types (e.g. `<UNKNOWN>`) → `BadSvType`
    /// - malformed breakend strings → `BadBreakend`
    /// - invalid sequence bases (including empty string) → `BadAlleleBases`
    pub fn parse(alt: &str) -> Result<Allele, BuildError> {
        if alt == "*" {
            return Ok(Allele::Star);
        }
        if alt == "<*>" {
            return Ok(Allele::Unspecified);
        }
        if alt.starts_with('<') && alt.ends_with('>') {
            let inner = &alt[1..alt.len() - 1];
            let mut parts = inner.split(':');
            let first = parts.next().unwrap_or("");
            let first_type = first.parse::<SvType>()?;
            let subtypes: Vec<String> = parts.map(|s| s.to_string()).collect();
            return Ok(Allele::Symbolic {
                first_type,
                subtypes,
            });
        }
        if alt.contains('[') || alt.contains(']') {
            return Allele::breakend_parse(alt);
        }
        if alt.len() > 1 && (alt.starts_with('.') || alt.ends_with('.')) {
            return Allele::breakend_parse(alt);
        }
        Allele::seq(alt)
    }
}

/// Single breakend: `.t` or `t.` where t is a non-empty sequence.
fn is_single_breakend(s: &str) -> bool {
    if let Some(rest) = s.strip_prefix('.') {
        return is_seq(rest);
    }
    if let Some(rest) = s.strip_suffix('.') {
        return is_seq(rest);
    }
    false
}

/// Paired breakend: `t[p[`, `t]p]`, `[p[t`, `]p]t` where both brackets are the
/// same char, t is a (possibly empty) sequence, and p is `chr:pos`.
fn is_paired_breakend(s: &str) -> bool {
    let bytes = s.as_bytes();
    let open = match bytes.iter().position(|&b| b == b'[' || b == b']') {
        Some(i) => i,
        None => return false,
    };
    let bracket = bytes[open];
    let close = match bytes.iter().rposition(|&b| b == b'[' || b == b']') {
        Some(i) if i != open => i,
        _ => return false,
    };
    if bytes[close] != bracket {
        return false;
    }
    let left = &s[..open];
    let mate = &s[open + 1..close];
    let right = &s[close + 1..];
    // Exactly one of left/right is the sequence side; the other is empty.
    if !is_seq_or_empty(left) || !is_seq_or_empty(right) {
        return false;
    }
    if left.is_empty() == right.is_empty() {
        return false; // need exactly one side with the replacement sequence
    }
    valid_mate(mate)
}

/// `chr:pos` with pos all-digits and a non-empty contig containing no brackets.
fn valid_mate(mate: &str) -> bool {
    match mate.rsplit_once(':') {
        Some((chrom, pos)) => {
            !chrom.is_empty()
                && !chrom.contains('[')
                && !chrom.contains(']')
                && !pos.is_empty()
                && pos.bytes().all(|b| b.is_ascii_digit())
        }
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq_validates() {
        assert_eq!(Allele::seq("GAT").unwrap().render(), "GAT");
        assert!(Allele::seq("GX").is_err());
    }

    #[test]
    fn symbolic_render_and_type_str() {
        let dup = Allele::duplication(["TANDEM"]);
        assert_eq!(dup.render(), "<DUP:TANDEM>");
        assert_eq!(dup.symbolic_type_str().as_deref(), Some("DUP:TANDEM"));
    }

    #[test]
    fn parse_dispatch() {
        assert!(matches!(Allele::parse("*").unwrap(), Allele::Star));
        assert!(matches!(Allele::parse("<*>").unwrap(), Allele::Unspecified));
        assert!(matches!(
            Allele::parse("<DEL>").unwrap(),
            Allele::Symbolic { .. }
        ));
        assert!(matches!(
            Allele::parse("T[chr2:5[").unwrap(),
            Allele::Breakend { single: false, .. }
        ));
        assert!(matches!(
            Allele::parse(".A").unwrap(),
            Allele::Breakend { single: true, .. }
        ));
        assert!(matches!(Allele::parse("ACGT").unwrap(), Allele::Seq(_)));
    }

    #[test]
    fn breakend_parse_rejects_junk() {
        assert!(Allele::breakend_parse("not-a-breakend").is_err());
        assert!(Allele::breakend_parse("G[chr2:321[").is_ok());
        assert!(Allele::breakend_parse("A.").is_ok());
    }

    #[test]
    fn parse_rejects_unsupported() {
        assert!(matches!(
            Allele::parse("<UNKNOWN>"),
            Err(crate::error::BuildError::BadSvType(_))
        ));
        assert!(matches!(
            Allele::parse("xyz"),
            Err(crate::error::BuildError::BadAlleleBases(_))
        ));
        assert!(Allele::parse("").is_err());
        // valid cases still pass:
        assert!(matches!(Allele::parse("*").unwrap(), Allele::Star));
        assert!(matches!(Allele::parse("<*>").unwrap(), Allele::Unspecified));
        assert!(matches!(
            Allele::parse("<DEL>").unwrap(),
            Allele::Symbolic { .. }
        ));
        assert!(matches!(
            Allele::parse("T[chr2:5[").unwrap(),
            Allele::Breakend { .. }
        ));
        assert!(matches!(Allele::parse("ACGT").unwrap(), Allele::Seq(_)));
    }
}
