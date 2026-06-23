use std::str::FromStr;

use crate::error::BuildError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SvType {
    Del,
    Ins,
    Dup,
    Inv,
    Cnv,
}

impl SvType {
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
    pub fn seq(bases: impl Into<String>) -> Result<Allele, BuildError> {
        let bases = bases.into();
        if !is_seq(&bases) {
            return Err(BuildError::BadAlleleBases(bases));
        }
        Ok(Allele::Seq(bases))
    }

    pub fn star() -> Allele {
        Allele::Star
    }

    pub fn unspecified() -> Allele {
        Allele::Unspecified
    }

    fn symbolic(first: SvType, subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::Symbolic {
            first_type: first,
            subtypes: subtypes.into_iter().map(Into::into).collect(),
        }
    }

    pub fn deletion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Del, subtypes)
    }

    pub fn insertion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Ins, subtypes)
    }

    pub fn duplication(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Dup, subtypes)
    }

    pub fn inversion(subtypes: impl IntoIterator<Item = impl Into<String>>) -> Allele {
        Allele::symbolic(SvType::Inv, subtypes)
    }

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

    pub fn render(&self) -> String {
        match self {
            Allele::Seq(b) => b.clone(),
            Allele::Star => "*".to_string(),
            Allele::Unspecified => "<*>".to_string(),
            Allele::Breakend { raw, .. } => raw.clone(),
            Allele::Symbolic { .. } => format!("<{}>", self.symbolic_type_str().unwrap()),
        }
    }

    /// Syntactic dispatch from a raw ALT string (never fails: junk falls back
    /// to a sequence allele, matching the Python `classify_allele`).
    pub fn parse(alt: &str) -> Allele {
        if alt == "*" {
            return Allele::Star;
        }
        if alt == "<*>" {
            return Allele::Unspecified;
        }
        if alt.starts_with('<') && alt.ends_with('>') {
            let inner = &alt[1..alt.len() - 1];
            let mut parts = inner.split(':');
            let first = parts.next().unwrap_or("");
            let first_type = first.parse::<SvType>().unwrap_or(SvType::Del);
            return Allele::Symbolic {
                first_type,
                subtypes: parts.map(|s| s.to_string()).collect(),
            };
        }
        if alt.contains('[') || alt.contains(']') {
            if let Ok(b) = Allele::breakend_parse(alt) {
                return b;
            }
        }
        if alt.len() > 1 && (alt.starts_with('.') || alt.ends_with('.')) {
            if let Ok(b) = Allele::breakend_parse(alt) {
                return b;
            }
        }
        Allele::Seq(alt.to_string())
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
        assert!(matches!(Allele::parse("*"), Allele::Star));
        assert!(matches!(Allele::parse("<*>"), Allele::Unspecified));
        assert!(matches!(Allele::parse("<DEL>"), Allele::Symbolic { .. }));
        assert!(matches!(
            Allele::parse("T[chr2:5["),
            Allele::Breakend { single: false, .. }
        ));
        assert!(matches!(
            Allele::parse(".A"),
            Allele::Breakend { single: true, .. }
        ));
        assert!(matches!(Allele::parse("ACGT"), Allele::Seq(_)));
    }

    #[test]
    fn breakend_parse_rejects_junk() {
        assert!(Allele::breakend_parse("not-a-breakend").is_err());
        assert!(Allele::breakend_parse("G[chr2:321[").is_ok());
        assert!(Allele::breakend_parse("A.").is_ok());
    }
}
