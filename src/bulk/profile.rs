//! The [`Profile`] schema: fitted statistics plus a dialed payload choice.
//!
//! A profile is deliberately split into two parts that must never be
//! conflated:
//!
//! - [`Fitted`] — statistics estimated from a real cohort (see
//!   `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`). Every
//!   value here must trace back to a fit; hand-picked numbers do not belong
//!   in this struct.
//! - [`Dialed`] — knobs a user picks explicitly, independent of any fit.

use crate::bulk::BulkError;

const GERMLINE_1KGP: &str = include_str!("../../profiles/germline-1kgp.json");

/// A named bundle of fitted statistics and dialed generation choices.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Profile {
    pub name: String,
    pub provenance: Provenance,
    pub fitted: Fitted,
    pub dialed: Dialed,
}

/// Where a profile's fitted statistics came from.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Provenance {
    pub source: String,
    pub n_samples_source: usize,
    pub n_variants_source: u64,
    pub fitted_on: String,
    pub fit_tool_version: String,
}

/// Statistics estimated from a real cohort. Never hand-pick values here.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Fitted {
    pub contigs: Vec<ContigStat>,
    pub gap_dist: Histogram,
    pub sfs: Histogram,
    pub variant_classes: ClassMix,
    pub indel_length: Histogram,
    pub titv: f64,
    pub multiallelic_rate: f64,
    pub missing_rate: f64,
    pub phased_rate: f64,
    pub ploidy: u8,
}

/// Per-contig variant count and density, as observed in the source cohort.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ContigStat {
    pub id: String,
    pub n_variants: u64,
    pub density_per_kb: f64,
}

/// A histogram with `weights.len() == edges.len() - 1`.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Histogram {
    pub edges: Vec<f64>,
    pub weights: Vec<f64>,
}

/// Relative frequency of each variant class; must sum to 1.0.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ClassMix {
    pub snp: f64,
    pub insertion: f64,
    pub deletion: f64,
    pub mnp: f64,
    pub complex: f64,
    pub symbolic: f64,
}

/// Generation choices a user dials in, independent of any fit.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Dialed {
    pub payload: Payload,
}

/// Which per-sample/per-record fields to synthesize.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Payload {
    GtOnly,
    GtVaf,
    Gatk,
    Mutect2,
}

impl Profile {
    /// Load a profile bundled with the crate by name.
    pub fn builtin(name: &str) -> Result<Profile, BulkError> {
        let src = match name {
            "germline-1kgp" => GERMLINE_1KGP,
            other => return Err(BulkError::UnknownProfile(other.to_string())),
        };
        let p = Profile::from_json(src)?;
        p.validate()?;
        Ok(p)
    }

    /// Parse a profile from JSON text.
    pub fn from_json(s: &str) -> Result<Profile, BulkError> {
        Ok(serde_json::from_str(s)?)
    }

    /// Check internal consistency of the fitted statistics.
    pub fn validate(&self) -> Result<(), BulkError> {
        self.fitted.gap_dist.validate()?;
        self.fitted.sfs.validate()?;
        self.fitted.indel_length.validate()?;
        self.fitted.variant_classes.validate()?;
        if self.fitted.ploidy == 0 {
            return Err(BulkError::Invalid("ploidy must be >= 1".into()));
        }
        for (label, v) in [
            ("multiallelic_rate", self.fitted.multiallelic_rate),
            ("missing_rate", self.fitted.missing_rate),
            ("phased_rate", self.fitted.phased_rate),
        ] {
            if !(0.0..=1.0).contains(&v) {
                return Err(BulkError::Invalid(format!("{label} must be in [0, 1]")));
            }
        }
        if self.fitted.contigs.is_empty() {
            return Err(BulkError::Invalid("need >= 1 contig".into()));
        }
        Ok(())
    }
}

impl Histogram {
    /// Check edge/weight shape and basic sanity (increasing edges,
    /// non-negative weights summing to a positive total).
    pub fn validate(&self) -> Result<(), BulkError> {
        if self.edges.len() < 2 {
            return Err(BulkError::Invalid("histogram needs >= 2 edges".into()));
        }
        if self.weights.len() + 1 != self.edges.len() {
            return Err(BulkError::Invalid(format!(
                "histogram weights ({}) must be edges ({}) - 1",
                self.weights.len(),
                self.edges.len()
            )));
        }
        // Reject NaN/Inf before any comparison-based check below: NaN
        // comparisons are always false, so `< 0.0`, `<= 0.0`, etc. would
        // otherwise let a NaN-poisoned histogram through silently.
        if self.weights.iter().any(|w| !w.is_finite()) {
            return Err(BulkError::Invalid(
                "histogram weights must be finite (no NaN or Inf)".into(),
            ));
        }
        if self.edges.iter().any(|e| !e.is_finite()) {
            return Err(BulkError::Invalid(
                "histogram edges must be finite (no NaN or Inf)".into(),
            ));
        }
        if self.weights.iter().any(|w| *w < 0.0) {
            return Err(BulkError::Invalid("histogram weights must be >= 0".into()));
        }
        if self.weights.iter().sum::<f64>() <= 0.0 {
            return Err(BulkError::Invalid("histogram weights must sum > 0".into()));
        }
        if self.edges.windows(2).any(|w| w[1] <= w[0]) {
            return Err(BulkError::Invalid(
                "histogram edges must be increasing".into(),
            ));
        }
        Ok(())
    }
}

impl ClassMix {
    /// Check that class frequencies sum to 1.0 (within floating-point tolerance).
    pub fn validate(&self) -> Result<(), BulkError> {
        // NaN poisons the sum below, but `(sum - 1.0).abs() > 1e-6` is false
        // for NaN, so check each component explicitly first.
        for (label, v) in [
            ("snp", self.snp),
            ("insertion", self.insertion),
            ("deletion", self.deletion),
            ("mnp", self.mnp),
            ("complex", self.complex),
            ("symbolic", self.symbolic),
        ] {
            if v.is_nan() {
                return Err(BulkError::Invalid(format!(
                    "variant_classes.{label} must not be NaN"
                )));
            }
        }
        let sum =
            self.snp + self.insertion + self.deletion + self.mnp + self.complex + self.symbolic;
        if (sum - 1.0).abs() > 1e-6 {
            return Err(BulkError::Invalid(format!(
                "variant_classes must sum to 1.0, got {sum}"
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_germline_loads_and_validates() {
        let p = Profile::builtin("germline-1kgp").unwrap();
        assert_eq!(p.name, "germline-1kgp");
        assert_eq!(p.dialed.payload, Payload::GtOnly);
        assert_eq!(p.fitted.ploidy, 2);
        p.validate().unwrap();
    }

    #[test]
    fn unknown_builtin_errors() {
        assert!(Profile::builtin("nope").is_err());
    }

    #[test]
    fn histogram_length_mismatch_is_rejected() {
        // weights must have exactly edges.len() - 1 entries
        let h = Histogram {
            edges: vec![0.0, 1.0, 2.0],
            weights: vec![1.0],
        };
        assert!(h.validate().is_err());
    }

    #[test]
    fn class_mix_must_sum_to_one() {
        let m = ClassMix {
            snp: 0.5,
            insertion: 0.1,
            deletion: 0.1,
            mnp: 0.0,
            complex: 0.0,
            symbolic: 0.0,
        };
        assert!(m.validate().is_err());
    }

    #[test]
    fn histogram_rejects_nan_weight() {
        // NaN weights must not slip past validation into sampling.
        let h = Histogram {
            edges: vec![0.0, 1.0, 2.0],
            weights: vec![1.0, f64::NAN],
        };
        assert!(h.validate().is_err());
    }

    #[test]
    fn histogram_rejects_nan_edge() {
        // NaN edges must not slip past validation into sampling.
        let h = Histogram {
            edges: vec![0.0, f64::NAN, 2.0],
            weights: vec![1.0, 1.0],
        };
        assert!(h.validate().is_err());
    }

    #[test]
    fn class_mix_rejects_nan_component() {
        // NaN in a component poisons the sum, but (sum - 1.0).abs() > 1e-6
        // is false for NaN, so the old check let this through silently.
        let m = ClassMix {
            snp: f64::NAN,
            insertion: 0.2,
            deletion: 0.2,
            mnp: 0.2,
            complex: 0.2,
            symbolic: 0.2,
        };
        assert!(m.validate().is_err());
    }

    #[test]
    fn payload_round_trips_through_serde() {
        let json = r#""mutect2""#;
        let p: Payload = serde_json::from_str(json).unwrap();
        assert_eq!(p, Payload::Mutect2);
        assert_eq!(serde_json::to_string(&p).unwrap(), json);
    }
}
