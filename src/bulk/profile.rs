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
const GERMLINE_1KGP_UNPHASED: &str = include_str!("../../profiles/germline-1kgp-unphased.json");
const SOMATIC_GDC: &str = include_str!("../../profiles/somatic-gdc.json");

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
    pub ploidy: u8,
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
            "germline-1kgp-unphased" => GERMLINE_1KGP_UNPHASED,
            "somatic-gdc" => SOMATIC_GDC,
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
        if self.dialed.ploidy == 0 {
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
        // Reject NaN/Inf before the sum check below: NaN poisons the sum,
        // and infinities (e.g. opposite-sign components) can also poison it
        // to NaN, but `(sum - 1.0).abs() > 1e-6` is false for NaN either
        // way, so check each component explicitly first.
        for (label, v) in [
            ("snp", self.snp),
            ("insertion", self.insertion),
            ("deletion", self.deletion),
            ("mnp", self.mnp),
            ("complex", self.complex),
            ("symbolic", self.symbolic),
        ] {
            if !v.is_finite() {
                return Err(BulkError::Invalid(format!(
                    "variant_classes.{label} must be finite (no NaN or Inf)"
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
        assert_eq!(p.dialed.ploidy, 2);
        p.validate().unwrap();
    }

    #[test]
    fn builtin_somatic_loads_and_validates() {
        let p = Profile::builtin("somatic-gdc").unwrap();
        assert_eq!(p.name, "somatic-gdc");
        assert_eq!(p.dialed.payload, Payload::GtVaf);
        p.validate().unwrap();
        assert_eq!(p.provenance.n_samples_source, 16007);
    }

    #[test]
    fn germline_profile_is_really_fitted_not_placeholder() {
        let p = Profile::builtin("germline-1kgp").unwrap();
        assert_eq!(p.provenance.n_samples_source, 3202);
        assert!(!p.provenance.source.contains("PLACEHOLDER"));
    }

    #[test]
    fn germline_sfs_is_empirical_not_neutral() {
        // A neutral 1/x SFS gives ~12% singletons; the real *unphased* 1kGP callset is
        // ~36%. This test is the guard that we fitted data rather than theory.
        //
        // It deliberately targets `germline-1kgp-unphased`, NOT `germline-1kgp`: the
        // phased panel has a singleton fraction of 0.0 because phasing is precisely what
        // removes unphaseable singletons. Pointing this guard at the phased profile would
        // assert something no phased file can satisfy.
        let p = Profile::builtin("germline-1kgp-unphased").unwrap();
        let total: f64 = p.fitted.sfs.weights.iter().sum();
        let singleton_frac = p.fitted.sfs.weights[0] / total;
        assert!(
            singleton_frac > 0.3,
            "singleton fraction {singleton_frac} looks neutral, not empirical"
        );
    }

    #[test]
    fn germline_variants_differ_in_phasing() {
        // The two germline profiles exist to capture a real trade-off that no single file
        // provides: phased data has no singletons, unphased data has no phase.
        let phased = Profile::builtin("germline-1kgp").unwrap();
        let unphased = Profile::builtin("germline-1kgp-unphased").unwrap();
        assert_eq!(phased.fitted.phased_rate, 1.0);
        assert_eq!(unphased.fitted.phased_rate, 0.0);
        assert_eq!(
            phased.provenance.n_samples_source,
            unphased.provenance.n_samples_source
        );
        let sf = |p: &Profile| p.fitted.sfs.weights[0] / p.fitted.sfs.weights.iter().sum::<f64>();
        assert!(
            sf(&phased) < 0.01,
            "phased panel should have ~no singletons"
        );
        assert!(
            sf(&unphased) > 0.3,
            "unphased callset should be singleton-rich"
        );
    }

    #[test]
    fn builtin_unphased_germline_loads_and_validates() {
        let p = Profile::builtin("germline-1kgp-unphased").unwrap();
        assert_eq!(p.name, "germline-1kgp-unphased");
        p.validate().unwrap();
        assert_eq!(p.provenance.n_samples_source, 3202);
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
    fn class_mix_rejects_infinite_components() {
        // Opposite-sign infinities poison the sum to NaN, but
        // (sum - 1.0).abs() > 1e-6 is false for NaN, so an is_nan()-only
        // per-component check lets this through silently.
        let m = ClassMix {
            snp: f64::INFINITY,
            insertion: f64::NEG_INFINITY,
            deletion: 0.0,
            mnp: 0.0,
            complex: 0.0,
            symbolic: 0.0,
        };
        assert!(m.validate().is_err());
    }

    #[test]
    fn ploidy_lives_in_dialed_not_fitted() {
        let p = Profile::builtin("germline-1kgp").unwrap();
        assert_eq!(p.dialed.ploidy, 2);
    }

    #[test]
    fn payload_round_trips_through_serde() {
        let json = r#""mutect2""#;
        let p: Payload = serde_json::from_str(json).unwrap();
        assert_eq!(p, Payload::Mutect2);
        assert_eq!(serde_json::to_string(&p).unwrap(), json);
    }
}
