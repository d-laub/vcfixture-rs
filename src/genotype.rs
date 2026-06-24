use crate::error::BuildError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// Allele indices in call order; `None` is a missing allele (`.`).
    pub alleles: Vec<Option<u32>>,
    /// One bool per separator; `true` = phased (`|`). Length == ploidy - 1.
    pub phased: Vec<bool>,
}

impl Genotype {
    pub fn parse(s: &str) -> Result<Genotype, BuildError> {
        let mut alleles = Vec::new();
        let mut phased = Vec::new();
        let mut token = String::new();
        for ch in s.chars() {
            if ch == '|' || ch == '/' {
                alleles.push(parse_allele(&token, s)?);
                token.clear();
                phased.push(ch == '|');
            } else {
                token.push(ch);
            }
        }
        alleles.push(parse_allele(&token, s)?);
        Ok(Genotype { alleles, phased })
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

fn parse_allele(tok: &str, full_gt: &str) -> Result<Option<u32>, BuildError> {
    if tok == "." {
        Ok(None)
    } else {
        tok.parse::<u32>()
            .map(Some)
            .map_err(|_| BuildError::BadGenotype(full_gt.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_render_roundtrip() {
        for s in ["0|1", "1/1", "./.", "0", ".|1"] {
            assert_eq!(Genotype::parse(s).unwrap().render(), s);
        }
    }

    #[test]
    fn phasing_and_ploidy() {
        assert!(Genotype::parse("0|1").unwrap().is_phased());
        assert!(!Genotype::parse("0/1").unwrap().is_phased());
        assert!(!Genotype::parse("0").unwrap().is_phased());
        assert_eq!(Genotype::parse("0|1|1").unwrap().ploidy(), 3);
        assert_eq!(Genotype::parse("./.").unwrap().alleles, vec![None, None]);
    }

    #[test]
    fn bad_genotype_errors() {
        assert!(matches!(
            Genotype::parse("0|x"),
            Err(crate::error::BuildError::BadGenotype(_))
        ));
        assert!(matches!(
            Genotype::parse("A/A"),
            Err(crate::error::BuildError::BadGenotype(_))
        ));
        // valid: "./." is Ok with [None, None]
        assert!(Genotype::parse("./.").is_ok());
    }
}
