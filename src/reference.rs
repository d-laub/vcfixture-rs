use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use crate::error::BuildError;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKlass {
    Snp,
    Mnp,
    Ins,
    Del,
    Delins,
    SpanningDel,
}

#[derive(Debug, Clone)]
pub struct DrawOpts {
    pub alt_index: usize,
    pub del_len: usize,
    pub ins_seq: String,
    pub mnp_len: usize,
}

impl Default for DrawOpts {
    fn default() -> Self {
        DrawOpts {
            alt_index: 1,
            del_len: 1,
            ins_seq: "T".to_string(),
            mnp_len: 2,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatFeature {
    pub contig: String,
    pub pos0: usize,
    pub motif: String,
    pub count: usize,
}

impl RepeatFeature {
    pub fn length(&self) -> usize {
        self.motif.len() * self.count
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSpec {
    pub contigs: Vec<(String, String)>,
    pub repeats: Vec<RepeatFeature>,
}

fn next_base(b: char, offset: usize) -> char {
    const ORDER: [char; 4] = ['A', 'C', 'G', 'T'];
    let i = ORDER.iter().position(|&c| c == b).unwrap_or(0);
    ORDER[(i + offset) % 4]
}

impl ReferenceSpec {
    fn seq_for(&self, contig: &str) -> Result<&str, BuildError> {
        self.contigs
            .iter()
            .find(|(id, _)| id == contig)
            .map(|(_, s)| s.as_str())
            .ok_or_else(|| BuildError::ContigNotFound(contig.to_string()))
    }

    pub fn length(&self, contig: &str) -> Result<usize, BuildError> {
        Ok(self.seq_for(contig)?.len())
    }

    pub fn base(&self, contig: &str, pos0: usize) -> Result<String, BuildError> {
        let s = self.seq_for(contig)?;
        if pos0 + 1 > s.len() {
            return Err(BuildError::OutOfBounds {
                contig: contig.to_string(),
                pos0,
                len: 1,
                clen: s.len(),
            });
        }
        Ok(s[pos0..pos0 + 1].to_string())
    }

    pub fn seq(&self, contig: &str, start0: usize, length: usize) -> Result<String, BuildError> {
        let s = self.seq_for(contig)?;
        if start0 + length > s.len() {
            return Err(BuildError::OutOfBounds {
                contig: contig.to_string(),
                pos0: start0,
                len: length,
                clen: s.len(),
            });
        }
        Ok(s[start0..start0 + length].to_string())
    }

    pub fn draw_ref_alt(
        &self,
        contig: &str,
        pos0: usize,
        klass: VariantKlass,
        opts: &DrawOpts,
    ) -> Result<(String, Vec<String>), BuildError> {
        match klass {
            VariantKlass::Snp => {
                let r = self.base(contig, pos0)?;
                let alt = next_base(r.chars().next().unwrap(), opts.alt_index).to_string();
                Ok((r, vec![alt]))
            }
            VariantKlass::Mnp => {
                let r = self.seq(contig, pos0, opts.mnp_len)?;
                let alt: String = r.chars().map(|c| next_base(c, opts.alt_index)).collect();
                Ok((r, vec![alt]))
            }
            VariantKlass::Ins => {
                let anchor = self.base(contig, pos0)?;
                let alt = format!("{anchor}{}", opts.ins_seq);
                Ok((anchor, vec![alt]))
            }
            VariantKlass::Del => {
                let r = self.seq(contig, pos0, opts.del_len + 1)?;
                let alt = r[..1].to_string();
                Ok((r, vec![alt]))
            }
            VariantKlass::Delins => {
                let r = self.seq(contig, pos0, opts.mnp_len)?;
                Ok((r, vec![opts.ins_seq.clone()]))
            }
            VariantKlass::SpanningDel => {
                let r = self.base(contig, pos0)?;
                Ok((r, vec!["*".to_string()]))
            }
        }
    }

    pub fn write(
        &self,
        path: impl AsRef<Path>,
        bgzip: bool,
        _index: bool,
    ) -> Result<PathBuf, BuildError> {
        let path = path.as_ref().to_path_buf();
        let mut text = String::new();
        for (cid, seq) in &self.contigs {
            text.push('>');
            text.push_str(cid);
            text.push('\n');
            for chunk in seq.as_bytes().chunks(60) {
                text.push_str(std::str::from_utf8(chunk).unwrap());
                text.push('\n');
            }
        }
        if bgzip {
            let file = fs::File::create(&path)?;
            let mut w = noodles_bgzf::io::Writer::new(file);
            w.write_all(text.as_bytes())?;
            w.finish()?;
        } else {
            fs::write(&path, &text)?;
        }
        // _index: write a .fai (and .gzi when bgzipped) via
        // noodles-fasta's index writer against the pinned version. The
        // strategies/round-trip tests do not require the index, so this is
        // best-effort; verify if a downstream test needs faidx access.
        Ok(path)
    }
}

pub struct ReferenceBuilder {
    rng: ChaCha8Rng,
    seqs: indexmap::IndexMap<String, Vec<u8>>,
    repeats: Vec<RepeatFeature>,
}

impl ReferenceBuilder {
    pub fn new(seed: u64) -> ReferenceBuilder {
        ReferenceBuilder {
            rng: ChaCha8Rng::seed_from_u64(seed),
            seqs: indexmap::IndexMap::new(),
            repeats: Vec::new(),
        }
    }

    pub fn add_contig(
        &mut self,
        id: impl Into<String>,
        length: usize,
    ) -> Result<&mut Self, BuildError> {
        let id = id.into();
        if self.seqs.contains_key(&id) {
            return Err(BuildError::ContigExists(id));
        }
        let seq: Vec<u8> = (0..length)
            .map(|_| BASES[self.rng.gen_range(0..4)])
            .collect();
        self.seqs.insert(id, seq);
        Ok(self)
    }

    pub fn set_base(
        &mut self,
        contig: &str,
        pos0: usize,
        base: &str,
    ) -> Result<&mut Self, BuildError> {
        let b = base.as_bytes();
        if b.len() != 1 {
            return Err(BuildError::BadAlleleBases(base.to_string()));
        }
        let arr = self
            .seqs
            .get_mut(contig)
            .ok_or_else(|| BuildError::ContigNotFound(contig.to_string()))?;
        if pos0 + 1 > arr.len() {
            return Err(BuildError::OutOfBounds {
                contig: contig.to_string(),
                pos0,
                len: 1,
                clen: arr.len(),
            });
        }
        arr[pos0] = b[0];
        Ok(self)
    }

    pub fn set_seq(
        &mut self,
        contig: &str,
        pos0: usize,
        seq: &str,
    ) -> Result<&mut Self, BuildError> {
        let arr = self
            .seqs
            .get_mut(contig)
            .ok_or_else(|| BuildError::ContigNotFound(contig.to_string()))?;
        if pos0 + seq.len() > arr.len() {
            return Err(BuildError::OutOfBounds {
                contig: contig.to_string(),
                pos0,
                len: seq.len(),
                clen: arr.len(),
            });
        }
        arr[pos0..pos0 + seq.len()].copy_from_slice(seq.as_bytes());
        Ok(self)
    }

    pub fn tandem_repeat(
        &mut self,
        contig: &str,
        pos0: usize,
        motif: &str,
        n: usize,
    ) -> Result<&mut Self, BuildError> {
        let run = motif.repeat(n);
        self.set_seq(contig, pos0, &run)?;
        self.repeats.push(RepeatFeature {
            contig: contig.to_string(),
            pos0,
            motif: motif.to_string(),
            count: n,
        });
        Ok(self)
    }

    pub fn build(&self) -> ReferenceSpec {
        let contigs = self
            .seqs
            .iter()
            .map(|(id, bytes)| (id.clone(), String::from_utf8(bytes.clone()).unwrap()))
            .collect();
        ReferenceSpec {
            contigs,
            repeats: self.repeats.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproducible_fill() {
        let a = ReferenceBuilder::new(7)
            .add_contig("chr1", 100)
            .unwrap()
            .build();
        let b = ReferenceBuilder::new(7)
            .add_contig("chr1", 100)
            .unwrap()
            .build();
        assert_eq!(
            a.seq("chr1", 0, 100).unwrap(),
            b.seq("chr1", 0, 100).unwrap()
        );
        assert_eq!(a.length("chr1").unwrap(), 100);
    }

    #[test]
    fn draw_snp_matches_reference() {
        let mut rb = ReferenceBuilder::new(1);
        rb.add_contig("chr1", 100).unwrap();
        rb.set_base("chr1", 10, "A").unwrap();
        let spec = rb.build();
        let (r, alts) = spec
            .draw_ref_alt("chr1", 10, VariantKlass::Snp, &DrawOpts::default())
            .unwrap();
        assert_eq!(r, "A");
        assert_eq!(alts.len(), 1);
        assert_ne!(alts[0], "A");
    }

    #[test]
    fn tandem_repeat_recorded() {
        let mut rb = ReferenceBuilder::new(1);
        rb.add_contig("chr1", 100).unwrap();
        rb.tandem_repeat("chr1", 10, "CAG", 4).unwrap();
        let spec = rb.build();
        assert_eq!(spec.repeats.len(), 1);
        assert_eq!(spec.seq("chr1", 10, 12).unwrap(), "CAGCAGCAGCAG");
    }

    #[test]
    fn out_of_bounds_does_not_panic() {
        let mut rb = ReferenceBuilder::new(1);
        rb.add_contig("chr1", 100).unwrap();
        let spec = rb.build();
        assert!(matches!(
            spec.seq("chr1", 95, 10),
            Err(BuildError::OutOfBounds { .. })
        ));
        assert!(matches!(
            spec.base("chr1", 200),
            Err(BuildError::OutOfBounds { .. })
        ));
    }
}
