#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// Allele indices in call order; `None` is a missing allele (`.`).
    pub alleles: Vec<Option<u32>>,
    /// One bool per separator; `true` = phased (`|`). Length == ploidy - 1.
    pub phased: Vec<bool>,
}

impl Genotype {
    pub fn parse(s: &str) -> Genotype {
        let mut alleles = Vec::new();
        let mut phased = Vec::new();
        let mut token = String::new();
        for ch in s.chars() {
            if ch == '|' || ch == '/' {
                alleles.push(parse_allele(&token));
                token.clear();
                phased.push(ch == '|');
            } else {
                token.push(ch);
            }
        }
        alleles.push(parse_allele(&token));
        Genotype { alleles, phased }
    }

    pub fn ploidy(&self) -> usize {
        self.alleles.len()
    }

    pub fn is_phased(&self) -> bool {
        !self.phased.is_empty() && self.phased.iter().all(|&p| p)
    }

    pub fn render(&self) -> String {
        let mut out = String::new();
        for (i, a) in self.alleles.iter().enumerate() {
            if i > 0 {
                out.push(if self.phased[i - 1] { '|' } else { '/' });
            }
            match a {
                Some(v) => out.push_str(&v.to_string()),
                None => out.push('.'),
            }
        }
        out
    }
}

fn parse_allele(tok: &str) -> Option<u32> {
    if tok == "." {
        None
    } else {
        Some(
            tok.parse()
                .expect("genotype allele index must be an integer or '.'"),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_render_roundtrip() {
        for s in ["0|1", "1/1", "./.", "0", ".|1"] {
            assert_eq!(Genotype::parse(s).render(), s);
        }
    }

    #[test]
    fn phasing_and_ploidy() {
        assert!(Genotype::parse("0|1").is_phased());
        assert!(!Genotype::parse("0/1").is_phased());
        assert!(!Genotype::parse("0").is_phased());
        assert_eq!(Genotype::parse("0|1|1").ploidy(), 3);
        assert_eq!(Genotype::parse("./.").alleles, vec![None, None]);
    }
}
