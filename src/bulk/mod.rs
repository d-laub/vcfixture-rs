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
    /// May overshoot, but never undershoot: each candidate is measured by
    /// writing it to a temp file through the same writer, format, and
    /// compression settings — and the same (absent) mid-stream flush
    /// cadence — as the real write, so the measured size is exactly what
    /// the real write produces, not an estimate. Per-contig counts are
    /// split proportional to fitted density, like [`Size::Records`]. See
    /// [`BulkSpec::write`] and [`BulkSpec::resolve_target_counts`].
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
    /// Records are generated and buffered **one contig at a time** (never
    /// the whole file at once), via two passes over `self.contig_ids`:
    ///
    /// 1. **Span pass**: generate each contig just to learn its populated
    ///    span (`pos_max`, via [`contig_span`]), then drop its records
    ///    immediately — only the span (a `u64`) is retained per contig.
    /// 2. Build and write the header from those spans (it must declare
    ///    every contig's length before any record is written).
    /// 3. **Write pass**: regenerate each contig — byte-identical to the
    ///    span pass, since [`BulkSpec::generate_contig`] is a pure function
    ///    of `(seed, contig_idx, n_records)`, never of what a previous call
    ///    computed — and write it immediately, dropping its records before
    ///    moving to the next contig.
    ///
    /// Peak memory is therefore bounded by the largest single contig's
    /// records, not the sum across every contig. Regenerating trades extra
    /// CPU (each contig's records are computed twice) for that memory
    /// bound; per-contig generation is itself parallelized (see
    /// [`BulkSpec::generate_contig`]), so this is cheap relative to the
    /// memory it avoids. See the [`BulkSpec`] doc comment for the
    /// contig-length and contig-name-resolution rules this enforces.
    pub fn write(self, path: impl AsRef<Path>) -> Result<Summary, BulkError> {
        let path = path.as_ref();
        self.profile.validate()?;
        if self.contig_ids.is_empty() {
            return Err(BulkError::Invalid("need >= 1 output contig".into()));
        }
        if self.n_samples == 0 {
            return Err(BulkError::Invalid("need >= 1 sample".into()));
        }
        // A duplicate output contig name would give each occurrence its own
        // independent position stream starting back at 0, so positions run
        // backwards across the file even though noodles dedupes the
        // `##contig` header line to one entry. A CSI built over such
        // out-of-order records silently drops region-query hits rather than
        // erroring anywhere — exactly the malformed-file class this crate
        // exists to prevent — so reject it up front instead.
        {
            let mut seen = std::collections::HashSet::with_capacity(self.contig_ids.len());
            for id in &self.contig_ids {
                if !seen.insert(id.as_str()) {
                    return Err(BulkError::Invalid(format!(
                        "duplicate output contig name: {id:?} (each requested contig \
                         must be unique; duplicates produce backwards positions and a \
                         CSI that silently drops region-query hits)"
                    )));
                }
            }
        }

        let fitted = &self.profile.fitted;
        // `SampleStats::value_for` (`gen.rs`) hard-codes `AD`/`PL` values
        // sized for exactly 2 allele calls per sample (diploid): `AD` is a
        // 2-element `[n_ref, n_alt]`, `PL` a fixed 3-element diploid
        // likelihood triple. A profile with `ploidy != 2` would declare a
        // `Number=G` FORMAT field but emit values that don't actually match
        // that ploidy's genotype-likelihood cardinality, producing a
        // malformed file rather than an error. `Profile::validate` only
        // requires `ploidy >= 1`, so this must be checked here.
        let keys = payload_keys(&self.payload);
        if (keys.contains(&"PL") || keys.contains(&"AD")) && fitted.ploidy != 2 {
            return Err(BulkError::Invalid(format!(
                "payload {:?} declares PL and/or AD, which are hard-coded for \
                 diploid (ploidy 2) genotype calls, but the profile's ploidy is {}",
                self.payload, fitted.ploidy
            )));
        }

        let samplers = Samplers::new(fitted)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.workers.get())
            .build()
            .map_err(|e| BulkError::Invalid(format!("failed to build worker pool: {e}")))?;

        let counts: Vec<u64> = match self.size {
            Size::RecordsPerContig(n) => vec![n; self.contig_ids.len()],
            Size::Records(total) => distribute_by_density(fitted, &self.contig_ids, total),
            Size::Target(target_bytes) => {
                self.resolve_target_counts(&pool, &samplers, fitted, target_bytes)?
            }
        };

        let spans: Vec<u64> = self
            .contig_ids
            .iter()
            .zip(&counts)
            .enumerate()
            .map(|(i, (id, &n))| {
                let recs = self.generate_contig(&pool, &samplers, fitted, id, i as u64, n);
                contig_span(&recs)
            })
            .collect();

        let header = self.build_header(&spans);
        let mut writer = BulkWriter::create(
            path,
            self.format,
            &header,
            self.compression_level,
            self.workers,
        )?;
        let mut summary = Summary::new(self.n_samples);
        for (i, (id, &n)) in self.contig_ids.iter().zip(&counts).enumerate() {
            let recs = self.generate_contig(&pool, &samplers, fitted, id, i as u64, n);
            for r in &recs {
                let buf = to_record_buf(&r.g, self.payload.clone(), r.phased);
                writer.write(&header, &buf)?;
                summary.observe(id, r.g.pos, r.g.class, &r.g.gts);
            }
            // No `writer.flush()` here: `MultithreadedWriter` dispatches a
            // compressed bgzf block once its ~64 KiB staging buffer fills
            // regardless, and `write()` never polls `compressed_bytes()`
            // (only `Size::Target`'s `measure_compressed_bytes` needs an
            // exact count, via a temp-file `finish_and_index`, not this live
            // counter). A per-contig flush here would force a block boundary
            // at every contig, which fragments and hurts compression, and —
            // critically — makes this write's bgzf block layout differ from
            // `measure_compressed_bytes`'s (which does not flush per
            // contig), which was exactly the bug: `Size::Target` could
            // measure a byte count that the real write would not reproduce.
            // Removing this call makes the two structurally identical, so
            // the measurement is now exact, not just close.
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
    /// contig with `length` set to that contig's populated span (`spans`,
    /// parallel to `self.contig_ids`; see [`contig_span`] and
    /// [`BulkSpec::write`]'s two-pass structure for how spans are computed
    /// without holding every contig's records in memory at once).
    fn build_header(&self, spans: &[u64]) -> vcf::Header {
        let mut hb = vcf::Header::builder();
        for i in 0..self.n_samples {
            hb = hb.add_sample_name(format!("s{i}"));
        }
        for &key in payload_keys(&self.payload) {
            hb = hb.add_format(key, format_map(key));
        }
        for (id, &span) in self.contig_ids.iter().zip(spans) {
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

    /// Resolves [`Size::Target`] to per-contig record *counts* (not the
    /// records themselves), generating candidates as it goes to measure
    /// them.
    ///
    /// This is the plan's "two-pass" approach, generalized to as many
    /// rounds as needed: each round generates a candidate `per_contig` (the
    /// *actual* records that would be written — not a proxy), writes it to
    /// a temp file via the real [`BulkWriter`] (so measured bytes are
    /// exactly what the real write would produce — same header, format,
    /// compression level, and, since neither this nor
    /// [`BulkSpec::write`] calls [`BulkWriter::flush`] mid-stream, the same
    /// bgzf block layout too — so the `finish_and_index`'d temp file's size
    /// on disk is exact, not an estimate), and checks it against the
    /// target. If short, it extrapolates the additional records needed from
    /// the observed bytes/record ratio (with a 15% margin so successive
    /// rounds converge quickly rather than repeatedly undershooting) and
    /// retries.
    ///
    /// Only the winning round's per-contig *counts* are returned.
    /// [`BulkSpec::write`] regenerates from them (its own two-pass
    /// span/write structure), rather than this function returning the
    /// generated records directly, so that at most one round's
    /// `Vec<Vec<Rec>>` — not every round's data, and not a whole extra
    /// file's worth of records held alongside the real write — is ever
    /// live at once.
    ///
    /// Both the initial guess and each round's top-up are split
    /// proportional to each contig's fitted density via
    /// [`distribute_by_density`] — the same helper [`Size::Records`] uses —
    /// rather than an even split, so `Size::Target`'s per-contig realism
    /// matches `Size::Records`'s. This is pure arithmetic on already-fitted
    /// statistics (no new randomness), so it does not affect determinism.
    fn resolve_target_counts(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        target_bytes: u64,
    ) -> Result<Vec<u64>, BulkError> {
        const INITIAL_PER_CONTIG: u64 = 500;
        const MAX_ROUNDS: usize = 25;

        let n_contigs = self.contig_ids.len() as u64;
        let mut per_contig_count =
            distribute_by_density(fitted, &self.contig_ids, INITIAL_PER_CONTIG * n_contigs);

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
                return Ok(per_contig_count);
            }

            let bytes_per_record = (bytes as f64 / total_records.max(1) as f64).max(1.0);
            let shortfall = (target_bytes - bytes) as f64;
            let extra = ((shortfall / bytes_per_record) * 1.15).ceil() as u64 + 1;
            let extra_split = distribute_by_density(fitted, &self.contig_ids, extra);
            for (c, e) in per_contig_count.iter_mut().zip(&extra_split) {
                *c += e;
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
    /// (dispatch to the compression thread pool is asynchronous). Because
    /// [`BulkSpec::write`] does not call [`BulkWriter::flush`] mid-stream
    /// either, this temp-file write is structurally identical to the real
    /// write — same header, same records, same (absent) flush cadence — so
    /// the byte count returned here is exactly what the real write
    /// produces for the same `per_contig`, not merely close to it.
    fn measure_compressed_bytes(&self, per_contig: &[Vec<Rec>]) -> Result<u64, BulkError> {
        let spans: Vec<u64> = per_contig.iter().map(|recs| contig_span(recs)).collect();
        let header = self.build_header(&spans);
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

/// The populated span of one contig's generated records — the maximum
/// position among them, or `1` if the contig has zero records. Positions
/// are strictly increasing within a contig (see
/// [`BulkSpec::generate_contig`]'s doc comment), so this equals the last
/// record's position, but is computed via `.max()` rather than relying on
/// that ordering to hold forever.
fn contig_span(recs: &[Rec]) -> u64 {
    recs.iter().map(|r| r.g.pos).max().unwrap_or(1)
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
/// 2. **`chr`-normalized match** (see [`normalize_contig_id`]): compare
///    `id` and each fitted id case-insensitively after stripping a leading
///    `chr`/`Chr`/`CHR` prefix from both sides. This resolves the actual
///    motivating real-data case — requested `chr1` against a profile fit
///    with bare id `1` — *by name*, not by a coincidence of ordering, so
///    it stays correct even when the caller's requested contigs are
///    reordered or a subset relative to the profile's fitted contigs
///    (e.g. `.contigs(["chr22", "chr1"])` still pairs `chr22` with fitted
///    id `22`'s stats, not the first fitted entry's).
/// 3. **Positional fallback**: `fitted.contigs[idx % fitted.contigs.len()]`,
///    used only when neither name-based rule above finds a match (e.g.
///    scaffold or otherwise unrelated output names) — a genuine last
///    resort now, not the primary bare-id-profile path. Both the
///    placeholder profile and a real 1000-Genomes/GDC fit list contigs in
///    chromosome order (`chr1..chrN` / `1..N`), so the `idx`-th requested
///    contig often still corresponds to the `idx`-th fitted contig even
///    when the ids themselves don't match textually. Wrapping (`%`) means
///    requesting more output contigs than the profile has fitted stats for
///    never panics and never silently zeroes out a stat — it cycles
///    through the fitted set, which is no worse a default than an
///    arbitrary unweighted one.
///
/// This never returns a synthetic zero/default stat: [`Profile::validate`]
/// (always run at the top of [`BulkSpec::write`]) rejects an empty
/// `fitted.contigs`, so the fallback always resolves to a real fitted
/// entry.
fn resolve_contig_stat<'a>(fitted: &'a Fitted, idx: usize, id: &str) -> &'a ContigStat {
    fitted
        .contigs
        .iter()
        .find(|c| c.id == id)
        .or_else(|| {
            let norm = normalize_contig_id(id);
            fitted
                .contigs
                .iter()
                .find(|c| normalize_contig_id(&c.id) == norm)
        })
        .unwrap_or(&fitted.contigs[idx % fitted.contigs.len()])
}

/// Lowercases `id` and strips one leading `chr` prefix (case-insensitively,
/// via the lowercased form), so `"chr1"`, `"Chr1"`, `"CHR1"`, and `"1"` all
/// normalize to `"1"`. Used by [`resolve_contig_stat`] to match a requested
/// output contig name against a fitted contig id by name, ahead of the
/// positional fallback.
fn normalize_contig_id(id: &str) -> String {
    let lower = id.to_ascii_lowercase();
    lower.strip_prefix("chr").unwrap_or(&lower).to_string()
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
