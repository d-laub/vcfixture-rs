/// Ordered genotypes per the VCF `Number=G` ordering.
pub fn genotype_ordering(ploidy: usize, n_alleles: usize) -> Vec<Vec<u32>> {
    assert!(ploidy >= 1, "ploidy must be >= 1");
    rec(ploidy, n_alleles)
}

fn rec(p: usize, n_alleles: usize) -> Vec<Vec<u32>> {
    if p == 1 {
        return (0..n_alleles as u32).map(|a| vec![a]).collect();
    }
    let mut out = Vec::new();
    for a in 0..n_alleles as u32 {
        for prefix in rec(p - 1, n_alleles) {
            if *prefix.last().unwrap() <= a {
                let mut g = prefix.clone();
                g.push(a);
                out.push(g);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diploid_biallelic_order() {
        // VCF Number=G order for ploidy 2, 2 alleles: 0/0, 0/1, 1/1
        assert_eq!(
            genotype_ordering(2, 2),
            vec![vec![0, 0], vec![0, 1], vec![1, 1]]
        );
    }

    #[test]
    fn count_matches_binomial() {
        // ploidy 2, 3 alleles => 6 genotypes
        assert_eq!(genotype_ordering(2, 3).len(), 6);
    }
}
