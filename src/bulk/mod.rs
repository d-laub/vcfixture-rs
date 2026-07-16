//! Bulk generation of realistic-enough VCF/BCF at benchmark scale.
//!
//! Unlike the fixture path ([`crate::build::VcfBuilder`]), bulk generation
//! streams records and derives no per-genotype oracle — see
//! `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`.
//!
//! [`BulkSpec`] is the public entry point: a builder over a [`Profile`] that
//! ties the samplers ([`sample`]), the record generator ([`gen`]), the
//! streaming writer ([`writer`]), and the summary truth ([`summary`])
//! together into one `write(path)` call.

pub mod gen;
pub mod profile;
pub mod sample;
pub mod summary;
pub mod writer;

use std::num::NonZero;
use std::path::Path;

use rand::Rng;
use rayon::prelude::*;

use noodles_vcf::{
    self as vcf,
    header::record::value::{
        map::{
            format::{Number as FormatNumber, Type as FormatType},
            Contig as ContigMap, Format as HeaderFormatMap,
        },
        Map,
    },
};

use gen::{block_rng, gen_record, to_record_buf, GenRecord};
use profile::{ContigStat, Fitted};
use sample::Samplers;
use writer::BulkWriter;

pub use profile::{Payload, Profile};
pub use summary::Summary;
pub use writer::Format;

/// Errors from bulk generation.
#[derive(Debug, thiserror::Error)]
pub enum BulkError {
    #[error("unknown builtin profile: {0}")]
    UnknownProfile(String),
    #[error("invalid profile: {0}")]
    Invalid(String),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}

/// How many records to generate, and how that maps onto per-contig counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Size {
    /// Exactly `n` records total, split across contigs proportional to each
    /// contig's fitted density (see [`resolve_contig_stat`]).
    Records(u64),
    /// Exactly `n` records for *each* requested contig.
    RecordsPerContig(u64),
    /// Grow the output until its compressed size is `>= n` bytes, then stop.
    /// May overshoot (never undershoot) — see [`BulkSpec::write`].
    Target(u64),
}

/// One generated record plus this call's phasing draw.
///
/// [`GenRecord`] is Task 6's type ([`crate::bulk::gen`]) and out of scope to
/// modify here, and it has no `phased` field (phasing is a per-record
/// decision only [`to_record_buf`] needs, not part of the site/genotype
/// generation `gen_record` performs) — so it is tracked alongside, not
/// inside, the generated record.
struct Rec {
    g: GenRecord,
    phased: bool,
}

/// Records generated for one block (`BulkSpec::BLOCK_SIZE` records, or fewer
/// for a contig's final partial block), with positions still relative to
/// the block's own start (see [`BulkSpec::generate_contig`]).
type BlockOutput = (Vec<Rec>, u64);

/// A builder for one bulk-generation run: samples, contigs, size, payload,
/// seed, output format, and worker count, ending in [`BulkSpec::write`].
///
/// # Contig length rule
///
/// The output header's `##contig` `length` is never a real chromosome
/// length — it is the *populated span* (`pos_max`) of whatever was actually
/// generated for that contig. See the design doc's "Contigs are declared at
/// fake lengths equal to the populated span" for why: declaring a real
/// hg38 length over a sparse, prefix-only population would make region
/// queries outside that prefix return nothing, which is exactly the
/// pathological case a benchmark must not trip over.
///
/// # Contig name resolution (spec-critical)
///
/// The contig names passed to [`BulkSpec::contigs`] are the **output**
/// names and are authoritative for what gets written — they are *not*
/// required to match `profile.fitted.contigs[].id`. This matters because
/// committed profiles are fit from real pvar files using bare contig ids
/// (`"1"`, `"2"`, ...), while callers (and this crate's own defaults)
/// request `chr1`-style names. See [`resolve_contig_stat`] for the exact
/// resolution rule.
pub struct BulkSpec {
    profile: Profile,
    n_samples: usize,
    contig_ids: Vec<String>,
    size: Size,
    payload: Payload,
    seed: u64,
    format: Format,
    workers: NonZero<usize>,
    compression_level: u8,
}

impl BulkSpec {
    /// Records generated per parallel unit of work ("block"). A block's RNG
    /// stream is a pure function of `(seed, block_idx)` via [`block_rng`],
    /// so this is also the granularity at which thread-count independence
    /// is achieved: rayon may compute blocks on any thread in any order,
    /// but [`Vec::from_par_iter`]/`.collect()` always assembles the results
    /// back in index order.
    const BLOCK_SIZE: u64 = 500;

    /// `block_idx` is derived as `contig_idx * CONTIG_BLOCK_STRIDE +
    /// local_block`, so a contig's stream never depends on how many
    /// contigs precede it. At `BLOCK_SIZE` records per block this allows up
    /// to `CONTIG_BLOCK_STRIDE * BLOCK_SIZE` (500 billion) records per
    /// contig before colliding with the next contig's block-index space —
    /// far beyond any realistic run (a 100 MB benchmark BCF is ~265k
    /// records total).
    const CONTIG_BLOCK_STRIDE: u64 = 1_000_000;

    /// Builds a spec from a profile, with defaults matching a small smoke
    /// test rather than a benchmark-scale run: 1 sample, `chr1..chr3`, 1000
    /// records per contig, the profile's own dialed payload, seed 0, BCF
    /// output, all available cores, and bgzf compression level 6.
    pub fn new(profile: Profile) -> BulkSpec {
        let payload = profile.dialed.payload.clone();
        let workers = std::thread::available_parallelism().unwrap_or(NonZero::new(1).unwrap());
        BulkSpec {
            profile,
            n_samples: 1,
            contig_ids: vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
            size: Size::RecordsPerContig(1000),
            payload,
            seed: 0,
            format: Format::Bcf,
            workers,
            compression_level: 6,
        }
    }

    /// Sets the sample count. Sample names are generated as `s0..s{n-1}`.
    pub fn samples(mut self, n: usize) -> BulkSpec {
        self.n_samples = n;
        self
    }

    /// Sets the output contig names, in the order they will be written.
    /// See the [`BulkSpec`] doc comment for how these resolve against the
    /// profile's fitted per-contig stats.
    pub fn contigs<I, S>(mut self, ids: I) -> BulkSpec
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.contig_ids = ids.into_iter().map(Into::into).collect();
        self
    }

    /// Sets how many records to generate (see [`Size`]).
    pub fn size(mut self, size: Size) -> BulkSpec {
        self.size = size;
        self
    }

    /// Sets which FORMAT fields to synthesize.
    pub fn payload(mut self, p: Payload) -> BulkSpec {
        self.payload = p;
        self
    }

    /// Sets the PRNG seed. Same seed + profile + spec => byte-identical
    /// output, regardless of `workers`.
    pub fn seed(mut self, seed: u64) -> BulkSpec {
        self.seed = seed;
        self
    }

    /// Sets the output container format.
    pub fn format(mut self, f: Format) -> BulkSpec {
        self.format = f;
        self
    }

    /// Sets the bgzf compression worker count. Never affects output bytes
    /// (see [`BulkSpec::generate_contig`] and `writer::tests::
    /// output_is_byte_identical_regardless_of_worker_count`), only speed.
    pub fn workers(mut self, n: NonZero<usize>) -> BulkSpec {
        self.workers = n;
        self
    }

    /// Sets the bgzf compression level (0-9).
    pub fn compression_level(mut self, level: u8) -> BulkSpec {
        self.compression_level = level;
        self
    }

    /// Generates and writes the file, then a CSI index and a
    /// `<path>.summary.json` alongside it.
    ///
    /// Records are buffered one contig at a time (never the whole file) so
    /// that contig length can be finalized before the header — which must
    /// precede any record — is written. See the [`BulkSpec`] doc comment
    /// for the contig-length and contig-name-resolution rules this
    /// enforces.
    pub fn write(self, path: impl AsRef<Path>) -> Result<Summary, BulkError> {
        let path = path.as_ref();
        self.profile.validate()?;
        if self.contig_ids.is_empty() {
            return Err(BulkError::Invalid("need >= 1 output contig".into()));
        }
        if self.n_samples == 0 {
            return Err(BulkError::Invalid("need >= 1 sample".into()));
        }

        let fitted = &self.profile.fitted;
        let samplers = Samplers::new(fitted)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.workers.get())
            .build()
            .map_err(|e| BulkError::Invalid(format!("failed to build worker pool: {e}")))?;

        let per_contig: Vec<Vec<Rec>> = match self.size {
            Size::RecordsPerContig(n) => self
                .contig_ids
                .iter()
                .enumerate()
                .map(|(i, id)| self.generate_contig(&pool, &samplers, fitted, id, i as u64, n))
                .collect(),
            Size::Records(total) => {
                let counts = distribute_by_density(fitted, &self.contig_ids, total);
                self.contig_ids
                    .iter()
                    .zip(counts)
                    .enumerate()
                    .map(|(i, (id, n))| {
                        self.generate_contig(&pool, &samplers, fitted, id, i as u64, n)
                    })
                    .collect()
            }
            Size::Target(target_bytes) => {
                self.generate_for_target(&pool, &samplers, fitted, target_bytes)?
            }
        };

        let header = self.build_header(&per_contig);
        let mut writer = BulkWriter::create(
            path,
            self.format,
            &header,
            self.compression_level,
            self.workers,
        )?;
        let mut summary = Summary::new(self.n_samples);
        for (id, recs) in self.contig_ids.iter().zip(&per_contig) {
            for r in recs {
                let buf = to_record_buf(&r.g, self.payload.clone(), r.phased);
                writer.write(&header, &buf)?;
                summary.observe(id, r.g.pos, r.g.class, &r.g.gts);
            }
            // Flush once per contig (a "record block" in the writer's
            // terms), not per record: `MultithreadedWriter` only dispatches
            // a compressed block once its ~64 KiB staging buffer fills, so
            // without any flush `compressed_bytes()` would read 0 for a
            // long time, and flushing every record would fragment output
            // into one bgzf block per record and defeat both compression
            // and parallelism.
            writer.flush()?;
        }
        writer.finish_and_index(path)?;

        let json = summary.to_json()?;
        let mut summary_path = path.as_os_str().to_os_string();
        summary_path.push(".summary.json");
        std::fs::write(&summary_path, json)?;

        Ok(summary)
    }

    /// Builds the header: sample names, one `##FORMAT` line per key the
    /// spec's [`Payload`] preset uses, and one `##contig` per requested
    /// contig with `length` set to that contig's populated span
    /// (`pos_max`, or `1` if the contig ended up with zero records).
    fn build_header(&self, per_contig: &[Vec<Rec>]) -> vcf::Header {
        let mut hb = vcf::Header::builder();
        for i in 0..self.n_samples {
            hb = hb.add_sample_name(format!("s{i}"));
        }
        for &key in payload_keys(&self.payload) {
            hb = hb.add_format(key, format_map(key));
        }
        for (id, recs) in self.contig_ids.iter().zip(per_contig) {
            let span = recs.iter().map(|r| r.g.pos).max().unwrap_or(1);
            let mut contig = Map::<ContigMap>::new();
            *contig.length_mut() = Some(span as usize);
            hb = hb.add_contig(id.clone(), contig);
        }
        hb.build()
    }

    /// Generates one contig's records.
    ///
    /// Parallelizes across `Self::BLOCK_SIZE`-record blocks with rayon
    /// (`into_par_iter` over block indices, each seeded independently via
    /// [`block_rng`]), then reassembles them in index order — this, not
    /// anything about `pool`'s thread count, is what makes output
    /// thread-count independent: `.collect()` on a rayon parallel iterator
    /// always preserves index order regardless of which thread computed
    /// which item, and each block's content is a pure function of
    /// `(seed, block_idx)`, never of thread identity or a shared mutable
    /// RNG.
    ///
    /// Positions must be strictly increasing across the *whole* contig
    /// (VCF requires sorted records), but each block can only compute
    /// positions relative to its own start while running in parallel with
    /// no knowledge of the previous block's total span. So each block
    /// generates with a block-local position starting at 0 (first record's
    /// position is its own first gap draw, `>= 1`), returning both its
    /// records and its local span (the last record's local position); a
    /// cheap sequential prefix sum over blocks then turns those into
    /// absolute contig positions.
    fn generate_contig(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        chrom: &str,
        contig_idx: u64,
        n_records: u64,
    ) -> Vec<Rec> {
        if n_records == 0 {
            return Vec::new();
        }

        let ploidy = fitted.ploidy;
        let n_samples = self.n_samples;
        let seed = self.seed;
        let n_blocks = n_records.div_ceil(Self::BLOCK_SIZE);

        let blocks: Vec<BlockOutput> = pool.install(|| {
            (0..n_blocks)
                .into_par_iter()
                .map(|local_block| {
                    let block_idx = contig_idx * Self::CONTIG_BLOCK_STRIDE + local_block;
                    let mut rng = block_rng(seed, block_idx);
                    let start = local_block * Self::BLOCK_SIZE;
                    let count = Self::BLOCK_SIZE.min(n_records - start);

                    let mut local_pos = 0u64;
                    let mut recs = Vec::with_capacity(count as usize);
                    for _ in 0..count {
                        local_pos += samplers.gap(&mut rng);
                        let g = gen_record(
                            &mut rng, samplers, chrom, local_pos, n_samples, ploidy, fitted,
                        );
                        // Phasing is a per-record draw, not part of
                        // `gen_record` (see `Rec`'s doc comment) — drawn
                        // from the same block-local RNG, right after the
                        // record it applies to, so the block's stream stays
                        // a pure function of `(seed, block_idx)` alone.
                        let phased = rng.gen::<f64>() < fitted.phased_rate;
                        recs.push(Rec { g, phased });
                    }
                    (recs, local_pos)
                })
                .collect()
        });

        let mut out = Vec::with_capacity(n_records as usize);
        let mut offset = 0u64;
        for (mut recs, span) in blocks {
            for r in &mut recs {
                r.g.pos += offset;
            }
            out.extend(recs);
            offset += span;
        }
        out
    }

    /// Resolves [`Size::Target`] to per-contig record counts, generating as
    /// it goes.
    ///
    /// This is the plan's "two-pass" approach, generalized to as many
    /// rounds as needed: each round generates a candidate `per_contig` (the
    /// *actual* records that would be written — not a proxy), writes it to
    /// a temp file via the real [`BulkWriter`] (so measured bytes are
    /// exactly what the real write would produce, not an estimate — the
    /// same header, format, and compression level, so `finish_and_index`'d
    /// file size on disk is exact), and checks it against the target. If
    /// short, it extrapolates the additional records needed from the
    /// observed bytes/record ratio (with a 15% margin so successive rounds
    /// converge quickly rather than repeatedly undershooting) and retries.
    /// The winning round's `per_contig` is returned and reused directly for
    /// the real write in [`BulkSpec::write`] — no wasted final regeneration
    /// pass.
    fn generate_for_target(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        target_bytes: u64,
    ) -> Result<Vec<Vec<Rec>>, BulkError> {
        const INITIAL_PER_CONTIG: u64 = 500;
        const MAX_ROUNDS: usize = 25;

        let n_contigs = self.contig_ids.len() as u64;
        let mut per_contig_count = vec![INITIAL_PER_CONTIG; self.contig_ids.len()];

        for _round in 0..MAX_ROUNDS {
            let per_contig: Vec<Vec<Rec>> = self
                .contig_ids
                .iter()
                .zip(&per_contig_count)
                .enumerate()
                .map(|(i, (id, &n))| self.generate_contig(pool, samplers, fitted, id, i as u64, n))
                .collect();

            let total_records: u64 = per_contig.iter().map(|r| r.len() as u64).sum();
            let bytes = self.measure_compressed_bytes(&per_contig)?;

            if bytes >= target_bytes {
                return Ok(per_contig);
            }

            let bytes_per_record = (bytes as f64 / total_records.max(1) as f64).max(1.0);
            let shortfall = (target_bytes - bytes) as f64;
            let extra = ((shortfall / bytes_per_record) * 1.15).ceil() as u64 + 1;
            let extra_per_contig = extra.div_ceil(n_contigs).max(1);
            for c in per_contig_count.iter_mut() {
                *c += extra_per_contig;
            }
        }

        Err(BulkError::Invalid(format!(
            "could not reach target size {target_bytes} bytes within {MAX_ROUNDS} rounds"
        )))
    }

    /// Writes a candidate `per_contig` to a throwaway temp file through the
    /// real [`BulkWriter`] and returns its exact on-disk size after
    /// `finish_and_index` — not the live `compressed_bytes()` counter,
    /// which the writer documents as lagging until the writer is finished
    /// (dispatch to the compression thread pool is asynchronous).
    fn measure_compressed_bytes(&self, per_contig: &[Vec<Rec>]) -> Result<u64, BulkError> {
        let header = self.build_header(per_contig);
        let tmp = tempfile::NamedTempFile::new()?;
        let tmp_path = tmp.path().to_path_buf();

        let mut w = BulkWriter::create(
            &tmp_path,
            self.format,
            &header,
            self.compression_level,
            self.workers,
        )?;
        for recs in per_contig {
            for r in recs {
                let buf = to_record_buf(&r.g, self.payload.clone(), r.phased);
                w.write(&header, &buf)?;
            }
        }
        w.finish_and_index(&tmp_path)?;

        let bytes = std::fs::metadata(&tmp_path)?.len();

        // `finish_and_index` may have written a `<tmp_path>.csi` companion
        // (Bcf only) that `NamedTempFile`'s `Drop` does not know about;
        // best-effort clean it up so repeated measurement rounds don't
        // litter the temp dir.
        let mut csi_path = tmp_path.as_os_str().to_os_string();
        csi_path.push(".csi");
        let _ = std::fs::remove_file(csi_path);

        Ok(bytes)
    }
}

/// Resolves a per-contig fitted statistic for one requested *output*
/// contig.
///
/// Output contig names (via [`BulkSpec::contigs`]) are chosen by the caller
/// and are authoritative for what gets written — they are **not** required
/// to match `profile.fitted.contigs[].id`. This matters concretely: the
/// committed profiles are fit from real pvar files that use bare contig ids
/// (`"1"`, `"2"`, `"10"`, ...), while the placeholder profile and this
/// crate's own defaults use `chr1`-style names. A profile's contig ids must
/// never be required to agree with a caller's requested output names.
///
/// Resolution order for the `idx`-th requested contig `id`:
///
/// 1. **Exact id match** against `fitted.contigs[].id`. Covers a profile
///    whose ids already agree with the requested names (e.g. the
///    placeholder profile, whose ids are literally `chr1`/`chr2`/`chr3`).
/// 2. **Positional fallback**: `fitted.contigs[idx % fitted.contigs.len()]`.
///    Both the placeholder profile and a real 1000-Genomes/GDC fit list
///    contigs in chromosome order (`chr1..chrN` / `1..N`), so the `idx`-th
///    requested contig corresponds to the `idx`-th fitted contig even when
///    the ids themselves don't match textually. Wrapping (`%`) means
///    requesting more output contigs than the profile has fitted stats for
///    never panics and never silently zeroes out a stat — it cycles
///    through the fitted set, which is no worse a default than an
///    arbitrary unweighted one.
///
/// This never returns a synthetic zero/default stat: [`Profile::validate`]
/// (always run at the top of [`BulkSpec::write`]) rejects an empty
/// `fitted.contigs`, so the positional fallback always resolves to a real
/// fitted entry.
fn resolve_contig_stat<'a>(fitted: &'a Fitted, idx: usize, id: &str) -> &'a ContigStat {
    fitted
        .contigs
        .iter()
        .find(|c| c.id == id)
        .unwrap_or(&fitted.contigs[idx % fitted.contigs.len()])
}

/// Splits `total` records across `contig_ids` proportional to each
/// contig's fitted density (via [`resolve_contig_stat`]), using the
/// largest-remainder method so the per-contig counts sum to exactly
/// `total`. Falls back to an even split if every resolved weight is zero
/// (a degenerate profile), rather than dividing by zero.
fn distribute_by_density(fitted: &Fitted, contig_ids: &[String], total: u64) -> Vec<u64> {
    let weights: Vec<f64> = contig_ids
        .iter()
        .enumerate()
        .map(|(i, id)| resolve_contig_stat(fitted, i, id).density_per_kb.max(0.0))
        .collect();
    let weight_sum: f64 = weights.iter().sum();
    let n = contig_ids.len() as u64;

    if weight_sum <= 0.0 {
        let base = total / n;
        let mut rem = total % n;
        return (0..n)
            .map(|_| {
                if rem > 0 {
                    rem -= 1;
                    base + 1
                } else {
                    base
                }
            })
            .collect();
    }

    let raw: Vec<f64> = weights
        .iter()
        .map(|w| w / weight_sum * total as f64)
        .collect();
    let mut counts: Vec<u64> = raw.iter().map(|r| r.floor() as u64).collect();
    let assigned: u64 = counts.iter().sum();
    let mut remainder = total.saturating_sub(assigned);

    let mut order: Vec<usize> = (0..raw.len()).collect();
    order.sort_by(|&a, &b| {
        let fa = raw[a] - raw[a].floor();
        let fb = raw[b] - raw[b].floor();
        fb.partial_cmp(&fa).unwrap_or(std::cmp::Ordering::Equal)
    });
    for &i in &order {
        if remainder == 0 {
            break;
        }
        counts[i] += 1;
        remainder -= 1;
    }
    counts
}

/// The ordered FORMAT key list for one [`Payload`] preset.
///
/// Must stay in sync with `gen::to_record_buf`'s own (private) `key_names`
/// match — duplicated here rather than shared because `gen.rs` is out of
/// scope to modify for this task. `payload_presets_all_write_readable_files`
/// (`tests/bulk.rs`) is the safety net: any drift would surface as a
/// "missing FORMAT header record" write error, not a silent mismatch.
fn payload_keys(payload: &Payload) -> &'static [&'static str] {
    match payload {
        Payload::GtOnly => &["GT"],
        Payload::GtVaf => &["GT", "VAF"],
        Payload::Gatk => &["GT", "AD", "DP", "GQ", "PL"],
        Payload::Mutect2 => &["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
    }
}

/// The header `##FORMAT` definition for one FORMAT key.
///
/// `noodles-bcf`'s encoder looks up each key's declared `Number`/`Type` in
/// the header and dispatches its integer/float/string encoding on it (see
/// `noodles_bcf::record::codec::encoder::samples::values::write_values`),
/// so a wrong or default (`Count(1)`, `String`) declaration is not a
/// cosmetic mismatch — it is a write-time "type mismatch" error. `GT`,
/// `AD`, `DP`, `GQ`, `PL` are VCF-reserved keys and get a correct
/// definition for free from `Map::from(key)`. `VAF`/`AF` (as a per-sample
/// FORMAT field, not the INFO field), `F1R2`, `F2R1`, and `SB` are GATK/
/// Mutect2 conventions, not VCF-reserved, so they need an explicit
/// definition matching exactly what `gen::SampleStats::value_for` emits.
fn format_map(key: &str) -> Map<HeaderFormatMap> {
    match key {
        "VAF" | "AF" => Map::<HeaderFormatMap>::new(
            FormatNumber::Count(1),
            FormatType::Float,
            "Variant allele frequency",
        ),
        "F1R2" => Map::<HeaderFormatMap>::new(
            FormatNumber::Count(2),
            FormatType::Integer,
            "Count of reads in F1R2 pair orientation supporting each allele",
        ),
        "F2R1" => Map::<HeaderFormatMap>::new(
            FormatNumber::Count(2),
            FormatType::Integer,
            "Count of reads in F2R1 pair orientation supporting each allele",
        ),
        "SB" => Map::<HeaderFormatMap>::new(
            FormatNumber::Count(4),
            FormatType::Integer,
            "Per-sample component statistics for Fisher strand bias",
        ),
        other => Map::<HeaderFormatMap>::from(other),
    }
}
