use crate::error::BuildError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NumberKind {
    Fixed,
    A,
    R,
    G,
    Dot,
    Flag,
}

/// VCF `Number=` cardinality descriptor. `count` is set only for `Fixed`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Number {
    pub kind: NumberKind,
    pub count: Option<u32>,
}

impl Number {
    pub const ONE: Number = Number {
        kind: NumberKind::Fixed,
        count: Some(1),
    };
    pub const A: Number = Number {
        kind: NumberKind::A,
        count: None,
    };
    pub const R: Number = Number {
        kind: NumberKind::R,
        count: None,
    };
    pub const G: Number = Number {
        kind: NumberKind::G,
        count: None,
    };
    pub const DOT: Number = Number {
        kind: NumberKind::Dot,
        count: None,
    };
    pub const FLAG: Number = Number {
        kind: NumberKind::Flag,
        count: None,
    };

    pub fn fixed(n: u32) -> Result<Number, BuildError> {
        // u32 cannot be negative; kept fallible to mirror the spec API and to
        // allow a future signed source. Always Ok for now.
        Ok(Number {
            kind: NumberKind::Fixed,
            count: Some(n),
        })
    }

    pub fn header_str(&self) -> String {
        match self.kind {
            NumberKind::Fixed => self.count.unwrap_or(0).to_string(),
            NumberKind::Flag => "0".to_string(),
            NumberKind::A => "A".to_string(),
            NumberKind::R => "R".to_string(),
            NumberKind::G => "G".to_string(),
            NumberKind::Dot => ".".to_string(),
        }
    }

    /// Resolve to a concrete value count for one record, or `None` when
    /// unbounded (`Number=.`).
    pub fn cardinality(&self, n_alt: usize, ploidy: usize) -> Option<usize> {
        match self.kind {
            NumberKind::Fixed => Some(self.count.unwrap_or(0) as usize),
            NumberKind::Flag => Some(0),
            NumberKind::A => Some(n_alt),
            NumberKind::R => Some(n_alt + 1),
            NumberKind::G => {
                let n_alleles = n_alt + 1;
                Some(binom(n_alleles + ploidy - 1, ploidy))
            }
            NumberKind::Dot => None,
        }
    }
}

/// Binomial coefficient C(n, k), computed iteratively to avoid overflow for
/// the small values used here.
fn binom(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result: usize = 1;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cardinalities() {
        assert_eq!(Number::A.cardinality(3, 2), Some(3));
        assert_eq!(Number::R.cardinality(3, 2), Some(4));
        assert_eq!(Number::ONE.cardinality(3, 2), Some(1));
        assert_eq!(Number::DOT.cardinality(3, 2), None);
        // Number=G: diploid, n_alleles = n_alt+1 = 3 => C(3+2-1,2)=C(4,2)=6
        assert_eq!(Number::G.cardinality(2, 2), Some(6));
        // haploid G => n_alleles
        assert_eq!(Number::G.cardinality(2, 1), Some(3));
    }

    #[test]
    fn header_tokens() {
        assert_eq!(Number::A.header_str(), "A");
        assert_eq!(Number::fixed(2).unwrap().header_str(), "2");
        assert_eq!(Number::FLAG.header_str(), "0");
        assert_eq!(Number::DOT.header_str(), ".");
    }
}
