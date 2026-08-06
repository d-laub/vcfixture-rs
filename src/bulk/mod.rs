//! Bulk generation of realistic-enough VCF/BCF at benchmark scale.
//!
//! Unlike the fixture path ([`crate::build::VcfBuilder`]), bulk generation
//! streams records and derives no per-genotype oracle — see
//! `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`.
//!
//! [`BulkSpec`] is the public entry point: a builder over a [`Profile`] that
//! ties the samplers ([`sample`]), the record generator ([`generate`]), the
//! streaming writer ([`writer`]), and the summary truth ([`summary`])
//! together into one `write(path)` call.

pub mod generate;
pub mod profile;
pub mod sample;
pub mod summary;
pub mod writer;

use std::collections::BTreeMap;
use std::num::NonZero;
use std::path::Path;

use rand::Rng;
use rayon::prelude::*;

use noodles_bgzf as bgzf;

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

use generate::{block_rng, gen_record, to_record_buf, GenRecord, Stream};
use profile::{ContigStat, Fitted};
use sample::Samplers;
use writer::BulkWriter;

pub use profile::{Payload, Profile};
pub use summary::Summary;
pub use writer::Format;

/// Errors from bulk generation.
///
/// Variants are grouped by *who* is at fault, because the message a user
/// sees depends on it: a malformed profile JSON, a spec the caller built
/// wrong, an unparseable argument, or a failure at generation time. The
/// `invalid profile:` prefix appears on exactly the first group -- those
/// messages name a profile field rather than describing themselves, so they
/// need the context; every other variant's message stands alone.
///
/// Profile-content failures share one `InvalidProfile(String)` rather than
/// getting a variant each: a caller cannot act differently on "histogram
/// edges must be increasing" than on "histogram weights must sum > 0", so
/// the extra variants would only ever reach `Display`.
///
/// # Non-exhaustive
///
/// This enum is `#[non_exhaustive]`: a downstream `match` must carry a
/// wildcard arm. The error classes here are still being discovered -- every
/// new validation check has so far meant a new variant, and each of those
/// would otherwise be a breaking change. Callers that want to branch should
/// do so on the variants they handle and fall through on the rest, which is
/// what a growing error type asks for anyway.
///
/// Branch on what you handle, and fall through on the rest:
///
/// ```
/// use vcfixture::bulk::BulkError;
///
/// fn exit_code(e: &BulkError) -> i32 {
///     match e {
///         BulkError::NoContigs => 2,
///         BulkError::NoSamples => 3,
///         _ => 1,
///     }
/// }
/// ```
///
/// Dropping that final arm does not compile, which
/// [`crate::compile_fail_guards`] pins.
///
/// **Adding a variant?** Add its arm to the guard block there too. Until you
/// do, that block fails to compile for the ordinary "forgot a variant"
/// reason and silently stops testing this attribute.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum BulkError {
    #[error("unknown builtin profile: {0}")]
    UnknownProfile(String),

    // --- profile content --------------------------------------------------
    /// A fitted or dialed statistic that fails [`Profile::validate`].
    #[error("invalid profile: {0}")]
    InvalidProfile(String),
    /// The one profile-content failure a caller can fix without editing the
    /// profile JSON -- by choosing a payload that doesn't emit AD/PL.
    #[error(
        "invalid profile: payload {payload:?} emits AD and/or PL, which are \
         hard-coded for diploid (ploidy 2) calls, but ploidy is {ploidy}"
    )]
    PayloadPloidy { payload: Payload, ploidy: u8 },

    // --- spec / caller validation -----------------------------------------
    #[error("need >= 1 output contig")]
    NoContigs,
    #[error("need >= 1 sample")]
    NoSamples,
    #[error(
        "duplicate output contig name: {0:?} (each requested contig must be \
         unique; duplicates produce backwards positions and a CSI that \
         silently drops region-query hits)"
    )]
    DuplicateContig(String),
    /// `Size::PerContig` omitted a contig that `.contigs()` requested. The
    /// names are the actionable payload, so they are structured rather than
    /// pre-formatted into a string.
    #[error(
        "Size::PerContig has no record count for requested contig(s) {0:?}; \
         every contig passed to .contigs() needs an entry (names are matched \
         exactly, with no chr-prefix normalization)"
    )]
    PerContigMissing(Vec<String>),
    /// `Size::PerContig` named a contig that `.contigs()` never requested.
    #[error(
        "Size::PerContig names contig(s) {0:?} that were not requested via \
         .contigs(); names are matched exactly, with no chr-prefix \
         normalization"
    )]
    PerContigUnknown(Vec<String>),
    /// A contig needs more blocks than `CONTIG_BLOCK_STRIDE` reserves,
    /// which would make its block indices collide with the next contig's
    /// and silently reuse that contig's PRNG streams. Smaller blocks (wide
    /// cohorts) bring this within reach, so it is checked rather than
    /// documented as unreachable.
    #[error(
        "contig {contig:?} needs {n_blocks} blocks, which exceeds the \
         {stride}-block stride reserved per contig; reduce its record count \
         or widen CONTIG_BLOCK_STRIDE"
    )]
    TooManyBlocks {
        contig: String,
        n_blocks: u64,
        stride: u64,
    },

    // --- argument parsing -------------------------------------------------
    #[error("bad size: {0:?} (expected a byte count, optionally suffixed KB/MB/GB)")]
    BadSize(String),
    /// A malformed `--records-for` token. One string variant rather than four:
    /// a caller cannot act differently on "missing `=`" than on "count is not
    /// a number" -- both mean *this token is malformed, here is how*.
    #[error("{0}")]
    BadRecordsFor(String),
    /// An out-of-range bgzf compression level, carrying the level the caller
    /// passed rather than a rendering of noodles' `TryFromU8Error` (which
    /// displays as `invalid input: 99`, giving the doubled and half-useful
    /// `invalid compression level: invalid input: 99`). The offending value
    /// is right there at the call site and a caller can branch on it, so it
    /// is kept structured.
    ///
    /// The accepted range is read off `noodles_bgzf`'s own public constants
    /// rather than written out as `0-9`: noodles gates its maximum on its
    /// `libdeflate` feature (9 without it, 12 with), and that feature can be
    /// switched on by *any* crate in a downstream build's graph through
    /// feature unification. Hard-coding the bound here would reintroduce, in
    /// a quieter form, the exact problem of asserting something about
    /// noodles that noodles can change.
    #[error("invalid compression level: {0} (expected {min}-{max})",
        min = bgzf::io::writer::CompressionLevel::NONE.get(),
        max = bgzf::io::writer::CompressionLevel::BEST.get())]
    CompressionLevel(u8),

    // --- profile loading --------------------------------------------------
    #[error("profile {path:?} is not a builtin name and could not be read as a file: {source}")]
    ProfileLoad {
        path: String,
        source: std::io::Error,
    },

    // --- runtime ----------------------------------------------------------
    #[error("failed to build worker pool: {0}")]
    WorkerPool(#[from] rayon::ThreadPoolBuildError),
    #[error(
        "could not reach target size {target_bytes} bytes within {corrections} corrective rounds"
    )]
    TargetNotReached {
        target_bytes: u64,
        corrections: usize,
    },

    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}

/// How many records to generate, and how that maps onto per-contig counts.
// Not `Copy`: `PerContig` carries a `BTreeMap`. `BulkSpec::size` takes its
// argument by value, so the builder API is unaffected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Size {
    /// Exactly `n` records total, split across contigs proportional to
    /// fitted per-contig variant count (`n_variants`) (see
    /// [`resolve_contig_stat`]).
    Records(u64),
    /// Exactly `n` records for *each* requested contig.
    RecordsPerContig(u64),
    /// Exactly the given number of records for each named contig.
    ///
    /// The profile's fitted per-contig statistics are not consulted at all:
    /// this is the escape hatch for reproducing a specific cohort's
    /// per-contig shape against a profile that was fit on a different one —
    /// the scheduling-benchmark case that [`Size::Records`]'s
    /// profile-derived split cannot express.
    ///
    /// Keys are matched **exactly** against the names passed to
    /// [`BulkSpec::contigs`] — unlike [`resolve_contig_stat`], there is no
    /// `chr`-prefix normalization here, because both lists come from the
    /// same caller in the same call. A requested contig with no entry, or
    /// an entry naming a contig that was not requested, is an error rather
    /// than a silent empty contig or a silently ignored key (see
    /// [`per_contig_counts`]).
    ///
    /// A count of `0` is legal and means "generate nothing for this
    /// contig"; it still gets a `##contig` header entry.
    PerContig(BTreeMap<String, u64>),
    /// Grow the output until its compressed size is `>= n` bytes, then stop.
    /// May overshoot, but never undershoot: the record count is calibrated
    /// from two cheap byte measurements (fitting `bytes ~= b0 + k*records`),
    /// corrected by a few slope-based rounds if needed, each of which writes
    /// a real candidate to a temp file through the same writer, format, and
    /// compression settings — and the same (absent) mid-stream flush
    /// cadence — as the real write, so the measured size is exactly what
    /// the real write produces, not an estimate. The winning temp file is
    /// then promoted (moved, not regenerated) onto the real destination.
    /// Per-contig counts are split proportional to fitted per-contig
    /// variant count (`n_variants`), like [`Size::Records`]. See
    /// [`BulkSpec::write`] and [`BulkSpec::resolve_target_counts`].
    Target(u64),
}

/// One generated record plus this call's phasing draw.
///
/// [`GenRecord`] is Task 6's type ([`crate::bulk::generate`]) and out of scope to
/// modify here, and it has no `phased` field (phasing is a per-record
/// decision only [`to_record_buf`] needs, not part of the site/genotype
/// generation `gen_record` performs) — so it is tracked alongside, not
/// inside, the generated record.
struct Rec {
    g: GenRecord,
    phased: bool,
}

/// Records generated for one block (`BulkSpec::MAX_BLOCK_RECORDS` records,
/// or fewer for a contig's final partial block), with positions still
/// relative to the block's own start (see [`BulkSpec::generate_contig`]).
// TODO(#22 task 6): deleted along with `generate_contig` and `contig_span`
// once the block pipeline replaces the regenerate-twice structure.
type BlockOutput = (Vec<Rec>, u64);

/// How one contig decomposes into blocks, and how long each block's slice
/// of the contig is.
///
/// `block_spans[i]` is the sum of block `i`'s gap draws — computed by a
/// pass that draws from [`Stream::Position`] only, never generating a
/// genotype (issue #22). The contig's populated span is their sum, and
/// block `i`'s absolute position offset is their exclusive prefix sum, so
/// one `Vec<u64>` carries everything the write pass needs.
// `#[allow(dead_code)]`: nothing calls `ContigLayout` or `compute_layouts`
// yet -- the block-pipeline task (#22 task 6) wires them into `write` and
// removes this allow.
#[allow(dead_code)]
pub(crate) struct ContigLayout {
    /// Records per block for this run (constant except for the final
    /// partial block).
    block_records: u64,
    n_records: u64,
    block_spans: Vec<u64>,
}

#[allow(dead_code)]
impl ContigLayout {
    fn n_blocks(&self) -> usize {
        self.block_spans.len()
    }

    /// Records in block `i` — `block_records`, or the remainder for the
    /// final partial block.
    fn block_len(&self, i: usize) -> u64 {
        debug_assert!(
            i < self.n_blocks(),
            "block index {i} out of range (n_blocks={})",
            self.n_blocks()
        );
        let start = (i as u64) * self.block_records;
        self.block_records.min(self.n_records - start)
    }

    /// The contig's populated span (its last record's position).
    fn span(&self) -> u64 {
        self.block_spans.iter().sum()
    }

    /// Absolute position offset for each block: the exclusive prefix sum of
    /// block spans.
    fn offsets(&self) -> Vec<u64> {
        let mut acc = 0u64;
        self.block_spans
            .iter()
            .map(|&s| {
                let o = acc;
                acc += s;
                o
            })
            .collect()
    }
}

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
    /// Upper bound on records per block. Below this, block size is set by
    /// cells (see [`BulkSpec::block_records`]).
    // `pub` (not just `pub(crate)`) so that `tests/bulk.rs` -- a separate
    // crate compiled against the public API -- can reference the real
    // constant instead of mirroring its value as a literal, which has
    // already silently regressed a determinism test's vacuity guard twice.
    // `#[doc(hidden)]` keeps it out of rendered docs.
    #[doc(hidden)]
    pub const MAX_BLOCK_RECORDS: u64 = 500;

    /// Target genotype cells per block. A block holds its records' full
    /// genotype vectors plus their encoded bytes, so its memory scales with
    /// `records * n_samples * ploidy`, not with records alone — a flat 500
    /// records is ~32 MB of genotypes at 32,000 samples, and the pipeline
    /// keeps `2 * workers` blocks in flight. 4M cells is ~4 MB of encoded
    /// gt-only output per block, which bounds in-flight memory at a few
    /// hundred MB on a 48-core run while leaving each block far more work
    /// than its per-block overhead (one buffer rewind, two RNG inits).
    const TARGET_CELLS_PER_BLOCK: u64 = 4_000_000;

    /// Records per parallel unit of work ("block"), sized so a block holds
    /// about [`BulkSpec::TARGET_CELLS_PER_BLOCK`] genotype cells, capped at
    /// [`BulkSpec::MAX_BLOCK_RECORDS`] and never below 1.
    ///
    /// A block's RNG streams are a pure function of
    /// `(seed, block_idx, stream)`, so this is also the granularity at
    /// which thread-count independence is achieved: rayon may compute
    /// blocks on any thread in any order, but `.collect()` assembles them
    /// back in index order. Block boundaries depend on `n_samples`,
    /// `ploidy`, and the record count — deliberately never on `workers`,
    /// which is what keeps output independent of thread count.
    #[doc(hidden)]
    pub fn block_records(n_samples: usize, ploidy: u8) -> u64 {
        let cells = (n_samples as u64)
            .saturating_mul(ploidy.max(1) as u64)
            .max(1);
        (Self::TARGET_CELLS_PER_BLOCK / cells).clamp(1, Self::MAX_BLOCK_RECORDS)
    }

    /// `block_idx` is derived as `contig_idx * CONTIG_BLOCK_STRIDE +
    /// local_block`, so a contig's stream never depends on how many
    /// contigs precede it. At `MAX_BLOCK_RECORDS` records per block this
    /// allows up to `CONTIG_BLOCK_STRIDE * MAX_BLOCK_RECORDS` (500 billion)
    /// records per contig before colliding with the next contig's
    /// block-index space — far beyond any realistic run (a 100 MB benchmark
    /// BCF is ~265k records total). [`BulkSpec::compute_layouts`] checks
    /// this bound explicitly rather than leaving it merely documented,
    /// since a wide cohort's smaller [`BulkSpec::block_records`] brings a
    /// pathologically large single-contig request within reach.
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
    ///
    /// # Coupling with the fitted site-frequency spectrum
    ///
    /// The profile's `fitted.sfs` histogram is fit against the *source*
    /// cohort's native size (`profile.provenance.n_samples_source`) — its
    /// edges are absolute allele counts observed in that cohort, not
    /// frequencies. Requesting a different sample count here does not
    /// clamp those absolute counts (which would silently saturate every
    /// high-AC bin to "every genotype is alt" whenever `n < n_samples_source`);
    /// instead each drawn allele count is rescaled to a frequency against
    /// the source cohort's `AN` and re-applied to this run's `AN`, so the
    /// realized alt-allele density matches the source cohort's regardless of
    /// how many samples are requested. See [`crate::bulk::sample::Samplers::
    /// allele_count`] for the exact rescaling formula.
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
            return Err(BulkError::NoContigs);
        }
        if self.n_samples == 0 {
            return Err(BulkError::NoSamples);
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
                    return Err(BulkError::DuplicateContig(id.clone()));
                }
            }
        }

        let fitted = &self.profile.fitted;
        // The sfs histogram's edges are absolute allele counts against the
        // *source* cohort's AN, not frequencies (see `BulkSpec::samples`'s
        // doc comment); `Samplers::allele_count` needs that source AN to
        // rescale a drawn count to whatever cohort size this run requests.
        let an_source = 2 * self.profile.provenance.n_samples_source as u64;
        let samplers = Samplers::new(fitted, an_source)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.workers.get())
            .build()?;

        // `&self.size` rather than `self.size`: `Size` is no longer `Copy`
        // (`PerContig` carries a map), and moving it out of `self` here
        // would forbid the `&self` method calls further down.
        let counts: Vec<u64> = match &self.size {
            Size::RecordsPerContig(n) => vec![*n; self.contig_ids.len()],
            Size::Records(total) => distribute_by_n_variants(fitted, &self.contig_ids, *total),
            Size::PerContig(map) => per_contig_counts(map, &self.contig_ids)?,
            Size::Target(target_bytes) => {
                let (_counts, tmp, _bytes, summary) =
                    self.resolve_target_counts(&pool, &samplers, fitted, *target_bytes)?;
                Self::promote_temp(tmp, path, self.format)?;
                let json = summary.to_json()?;
                let mut summary_path = path.as_os_str().to_os_string();
                summary_path.push(".summary.json");
                std::fs::write(&summary_path, json)?;
                return Ok(summary);
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
                let buf = to_record_buf(&r.g, &self.payload, r.phased);
                writer.write(&header, &buf)?;
                summary.observe(id, r.g.pos, r.g.class, &r.g.gts);
            }
            // No mid-stream flush here: `MultithreadedWriter` dispatches a
            // compressed bgzf block once its ~64 KiB staging buffer fills
            // regardless, and forcing one at every contig boundary would
            // fragment the output and hurt compression. It would also make
            // this write's bgzf block layout differ from
            // `measure_compressed_bytes`'s (which likewise never flushes
            // per contig) — exactly the bug `Size::Target` used to have:
            // measuring a byte count the real write would not reproduce.
            // Neither path flushing mid-stream keeps the two structurally
            // identical, so the measurement stays exact, not just close.
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

    /// Computes every contig's [`ContigLayout`] in one flat parallel pass
    /// that draws **only** gaps.
    ///
    /// This is the fix for issue #22's Finding 1. Contig lengths must be
    /// known before the header is written, and the header must precede any
    /// record; the old structure got them by generating every contig in
    /// full and keeping one `u64`, which meant every genotype in the file
    /// was drawn twice. Because positions now come from
    /// [`Stream::Position`] (see [`block_rng`]), the same numbers fall out
    /// of `O(total_records)` gap draws with no per-sample work at all.
    // `#[allow(dead_code)]`: not called yet -- the block-pipeline task (#22
    // task 6) wires this into `write` and removes this allow.
    #[allow(dead_code)]
    fn compute_layouts(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        _fitted: &Fitted,
        counts: &[u64],
    ) -> Result<Vec<ContigLayout>, BulkError> {
        // Ploidy lives on `self.profile.dialed`, not `Fitted` -- the brief's
        // sketch of this signature named `fitted.ploidy`, which does not
        // exist on `Fitted` (see `src/bulk/profile.rs`); `fitted` is kept as
        // a parameter to match the brief's interface and `generate_contig`'s
        // sibling signature, but ploidy is read from `self` instead.
        let block_records = Self::block_records(self.n_samples, self.profile.dialed.ploidy);
        let seed = self.seed;

        // Build each contig's `ContigLayout` skeleton (`block_spans` filled
        // in below, once the parallel gap-sum pass returns) up front, so
        // `ContigLayout::block_len` is the single place a block's record
        // count is computed — not re-derived a second time here to build
        // the job list. Two independent computations of the same partition
        // is exactly the silent-divergence mode `BulkError::TooManyBlocks`
        // exists to prevent.
        let mut layouts: Vec<ContigLayout> = Vec::with_capacity(counts.len());
        // (contig_idx, local_block, records_in_block), flattened so every
        // block in the run is one parallel work item regardless of which
        // contig it belongs to.
        let mut jobs: Vec<(u64, u64, u64)> = Vec::new();
        for (ci, (&n, id)) in counts.iter().zip(&self.contig_ids).enumerate() {
            let n_blocks = n.div_ceil(block_records);
            // `>`, not `>=`: the largest local block index is `n_blocks -
            // 1`, so `block_idx` collides with the next contig's stream
            // only once `n_blocks` exceeds `CONTIG_BLOCK_STRIDE` --
            // `n_blocks == CONTIG_BLOCK_STRIDE` fits exactly.
            if n_blocks > Self::CONTIG_BLOCK_STRIDE {
                return Err(BulkError::TooManyBlocks {
                    contig: id.clone(),
                    n_blocks,
                    stride: Self::CONTIG_BLOCK_STRIDE,
                });
            }
            let layout = ContigLayout {
                block_records,
                n_records: n,
                block_spans: vec![0; n_blocks as usize],
            };
            for lb in 0..layout.n_blocks() {
                jobs.push((ci as u64, lb as u64, layout.block_len(lb)));
            }
            layouts.push(layout);
        }

        let spans: Vec<u64> = pool.install(|| {
            jobs.par_iter()
                .map(|&(ci, lb, count)| {
                    let block_idx = ci * Self::CONTIG_BLOCK_STRIDE + lb;
                    let mut rng = block_rng(seed, block_idx, Stream::Position);
                    (0..count).map(|_| samplers.gap(&mut rng)).sum()
                })
                .collect()
        });

        let mut at = 0usize;
        for layout in &mut layouts {
            let n_blocks = layout.n_blocks();
            layout
                .block_spans
                .copy_from_slice(&spans[at..at + n_blocks]);
            at += n_blocks;
        }
        Ok(layouts)
    }

    /// Generates one contig's records.
    ///
    /// Parallelizes across `Self::MAX_BLOCK_RECORDS`-record blocks with
    /// rayon (`into_par_iter` over block indices, each seeded independently
    /// via [`block_rng`]), then reassembles them in index order — this, not
    /// anything about `pool`'s thread count, is what makes output
    /// thread-count independent: `.collect()` on a rayon parallel iterator
    /// always preserves index order regardless of which thread computed
    /// which item, and each block's position and content streams are each a
    /// pure function of `(seed, block_idx, stream)`, never of thread
    /// identity or a shared mutable RNG.
    ///
    /// Positions must be strictly increasing across the *whole* contig
    /// (VCF requires sorted records), but each block can only compute
    /// positions relative to its own start while running in parallel with
    /// no knowledge of the previous block's total span. So each block
    /// generates with a block-local position starting at 0 (first record's
    /// position is its own first gap draw, `>= 1`), returning both its
    /// records and its local span (the last record's local position); a
    /// cheap sequential prefix sum over blocks then turns those into
    /// absolute contig positions. Because gap draws come from their own
    /// stream, a block's local positions are a pure function of `(seed,
    /// block_idx, count)` alone — independent of `n_samples`, ploidy, and
    /// payload — which is the form Task 5's positions-only span pass
    /// consumes.
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

        let ploidy = self.profile.dialed.ploidy;
        let n_samples = self.n_samples;
        let seed = self.seed;
        // Flat `MAX_BLOCK_RECORDS`, not cell-sized `Self::block_records`:
        // `generate_contig` is deleted in #22 task 6 once the block
        // pipeline replaces this regenerate-twice structure, so it keeps
        // its pre-task-5 sizing (previously the `BLOCK_SIZE` constant, same
        // value) rather than adopting cell-based sizing here too.
        let n_blocks = n_records.div_ceil(Self::MAX_BLOCK_RECORDS);

        let blocks: Vec<BlockOutput> = pool.install(|| {
            (0..n_blocks)
                .into_par_iter()
                .map(|local_block| {
                    let block_idx = contig_idx * Self::CONTIG_BLOCK_STRIDE + local_block;
                    let mut pos_rng = block_rng(seed, block_idx, Stream::Position);
                    let mut rng = block_rng(seed, block_idx, Stream::Content);
                    let start = local_block * Self::MAX_BLOCK_RECORDS;
                    let count = Self::MAX_BLOCK_RECORDS.min(n_records - start);

                    let mut local_pos = 0u64;
                    let mut recs = Vec::with_capacity(count as usize);
                    for _ in 0..count {
                        local_pos += samplers.gap(&mut pos_rng);
                        let g = gen_record(
                            &mut rng, samplers, chrom, local_pos, n_samples, ploidy, fitted,
                        );
                        // Phasing is a per-record draw, not part of
                        // `gen_record` (see `Rec`'s doc comment) — drawn
                        // from the same block-local Content stream, right
                        // after the record it applies to, so that stream
                        // stays a pure function of `(seed, block_idx,
                        // Stream::Content)` alone.
                        let phased = rng.random::<f64>() < fitted.phased_rate;
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

    /// Resolves [`Size::Target`] to per-contig record *counts*, the
    /// already-written temp file those counts produced (byte length `>=
    /// target_bytes`), its byte length, and its `Summary` — so the caller
    /// ([`BulkSpec::write`]) can promote the temp straight to the real
    /// destination instead of regenerating a third time.
    ///
    /// Two cheap calibration points (`1_000` and `2_000` records per
    /// contig, measured bytes-only via [`BulkSpec::measure_compressed_bytes`])
    /// fit `bytes ~= b0 + k*records` (`k` bytes/record, `b0` the fixed
    /// header/index cost), which gives a direct count estimate for
    /// `target_bytes` in one step rather than the old scheme's up-to-25
    /// rounds of repeated doubling (each round generating every contig
    /// *twice* to measure, then discarding the result). That estimate is
    /// then corrected, at most [`MAX_CORRECTIONS`] times: each round writes
    /// the current guess to a real temp file via [`BulkSpec::write_to_temp`]
    /// (building the `Summary` for free, since this write is no longer a
    /// throwaway measurement -- it may be the one promoted), and if it's
    /// still short, tops up every contig's count proportionally (via
    /// [`distribute_by_n_variants`], same as the initial split), plus a 2%
    /// margin so rounds converge instead of oscillating just under the
    /// target.
    ///
    /// The top-up itself does **not** keep using the calibration slope `k`:
    /// `k` is fitted at the 1-2k-records/contig calibration scale and can
    /// mis-estimate marginal bytes/record at the (usually much larger)
    /// target scale, which made the top-up systematically under-buy and
    /// approach `target_bytes` asymptotically from below without ever
    /// crossing it within the round budget. From round 1 onward, each
    /// round instead refits a *local* slope from the two most recent real
    /// measurements (`(bytes - bytes_prev) / (records - records_prev)`),
    /// which tracks the true marginal cost at the actual scale in play;
    /// round 0 and any round where the local refit comes out non-positive
    /// or non-finite fall back to `k`. As a second guard against
    /// stalling — a top-up of only 1-2 records near the target can land
    /// inside an existing bgzf block with zero or even negative byte
    /// delta, since deflate's output isn't strictly monotonic in input
    /// size — any round that made no byte progress at all over the
    /// previous round floors its own top-up at twice the previous round's.
    ///
    /// Both the initial guess and each round's top-up are split
    /// proportional to each contig's fitted per-contig variant count
    /// (`n_variants`) via [`distribute_by_n_variants`] — the same helper
    /// [`Size::Records`] uses — rather than an even split, so
    /// `Size::Target`'s per-contig realism matches `Size::Records`'s. This
    /// is pure arithmetic on already-fitted statistics (no new randomness),
    /// so it does not affect determinism.
    fn resolve_target_counts(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        target_bytes: u64,
    ) -> Result<(Vec<u64>, tempfile::NamedTempFile, u64, Summary), BulkError> {
        let n_contigs = self.contig_ids.len() as u64;

        // Two calibration points; c2 = 2*c1 so the slope is well-conditioned.
        let split1 = distribute_by_n_variants(fitted, &self.contig_ids, 1_000 * n_contigs);
        let split2 = distribute_by_n_variants(fitted, &self.contig_ids, 2_000 * n_contigs);
        let bytes1 = self.measure_compressed_bytes(pool, samplers, fitted, &split1)?;
        let bytes2 = self.measure_compressed_bytes(pool, samplers, fitted, &split2)?;

        let r1 = split1.iter().sum::<u64>() as f64;
        let r2 = split2.iter().sum::<u64>() as f64;
        // bytes ~= b0 + k*records ; k bytes/record, b0 the fixed header cost.
        let k = ((bytes2 as f64 - bytes1 as f64) / (r2 - r1)).max(1e-9);
        let b0 = bytes1 as f64 - k * r1;

        // Direct count; never below the larger calibration (a known-good
        // measurement) and never below 1 record/contig.
        let want = (((target_bytes as f64 - b0) / k).ceil() as i64).max(r2 as i64) as u64;
        let mut counts = distribute_by_n_variants(fitted, &self.contig_ids, want);

        // Slope-based correction; usually converges in 1-2 rounds, but the
        // local-refit and no-progress guards below (see the doc comment)
        // can need a few more at the margin.
        const MAX_CORRECTIONS: usize = 6;
        let mut prev: Option<(u64, u64, u64)> = None; // (records, bytes, extra) from the last round
        for _ in 0..MAX_CORRECTIONS {
            let (tmp, bytes, summary) = self.write_to_temp(pool, samplers, fitted, &counts)?;
            if bytes >= target_bytes {
                return Ok((counts, tmp, bytes, summary));
            }
            let records: u64 = counts.iter().sum();

            let k_eff = prev
                .and_then(|(prev_records, prev_bytes, _)| {
                    let dr = records as f64 - prev_records as f64;
                    let db = bytes as f64 - prev_bytes as f64;
                    let cand = db / dr;
                    (cand.is_finite() && cand > 0.0).then_some(cand)
                })
                .unwrap_or(k);

            let shortfall = (target_bytes - bytes) as f64;
            let mut extra = ((shortfall / k_eff) * 1.02).ceil() as u64 + 1;

            // No-progress guard: a top-up this small can be absorbed into
            // an existing bgzf block with ~zero byte delta, stalling the
            // slope estimate forever; force forward progress.
            if let Some((_, prev_bytes, prev_extra)) = prev {
                if bytes <= prev_bytes {
                    extra = extra.max(prev_extra.saturating_mul(2)).max(1);
                }
            }

            let extra_split = distribute_by_n_variants(fitted, &self.contig_ids, extra);
            for (c, e) in counts.iter_mut().zip(&extra_split) {
                *c += e;
            }
            prev = Some((records, bytes, extra));

            // `write_to_temp` may have left a `<tmp_path>.csi` companion
            // (Bcf only) that `NamedTempFile`'s `Drop` does not know about;
            // best-effort clean it up, mirroring
            // `BulkSpec::measure_compressed_bytes`, so repeated corrective
            // rounds don't litter the temp dir.
            if matches!(self.format, Format::Bcf) {
                let mut csi_path = tmp.path().as_os_str().to_os_string();
                csi_path.push(".csi");
                let _ = std::fs::remove_file(csi_path);
            }
            drop(tmp); // discard the under-target temp before regenerating
        }

        Err(BulkError::TargetNotReached {
            target_bytes,
            corrections: MAX_CORRECTIONS,
        })
    }

    /// Like [`BulkSpec::measure_compressed_bytes`], but builds the
    /// [`Summary`] during the write pass and returns the live temp file
    /// instead of deleting it, so the caller can promote it to the real
    /// destination via [`BulkSpec::promote_temp`]. Byte-exact: identical
    /// header, records, and (absent) flush cadence as [`BulkSpec::write`].
    fn write_to_temp(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        per_contig_count: &[u64],
    ) -> Result<(tempfile::NamedTempFile, u64, Summary), BulkError> {
        let spans: Vec<u64> = self
            .contig_ids
            .iter()
            .zip(per_contig_count)
            .enumerate()
            .map(|(i, (id, &n))| {
                let recs = self.generate_contig(pool, samplers, fitted, id, i as u64, n);
                contig_span(&recs)
            })
            .collect();

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
        let mut summary = Summary::new(self.n_samples);
        for (i, (id, &n)) in self.contig_ids.iter().zip(per_contig_count).enumerate() {
            let recs = self.generate_contig(pool, samplers, fitted, id, i as u64, n);
            for r in &recs {
                let buf = to_record_buf(&r.g, &self.payload, r.phased);
                w.write(&header, &buf)?;
                summary.observe(id, r.g.pos, r.g.class, &r.g.gts);
            }
        }
        w.finish_and_index(&tmp_path)?;
        let bytes = std::fs::metadata(&tmp_path)?.len();
        Ok((tmp, bytes, summary))
    }

    /// Moves a written temp file (and, for BCF, its `.csi` companion) onto
    /// the real destination. Rename when possible; falls back to copy across
    /// filesystems (`TMPDIR` may differ from the output dir) via
    /// [`move_file`]/`NamedTempFile::persist`.
    fn promote_temp(
        tmp: tempfile::NamedTempFile,
        dest: &Path,
        format: Format,
    ) -> Result<(), BulkError> {
        let tmp_path = tmp.path().to_path_buf();
        if matches!(format, Format::Bcf) {
            let mut src_csi = tmp_path.as_os_str().to_os_string();
            src_csi.push(".csi");
            let mut dst_csi = dest.as_os_str().to_os_string();
            dst_csi.push(".csi");
            move_file(Path::new(&src_csi), Path::new(&dst_csi))?;
        }
        match tmp.persist(dest) {
            Ok(_) => Ok(()),
            Err(e) => {
                std::fs::copy(e.file.path(), dest)?;
                Ok(())
            }
        }
    }

    /// Measures the exact compressed byte size that `per_contig_count`
    /// would produce, by actually generating and writing it to a throwaway
    /// temp file through the real [`BulkWriter`] — **one contig at a time**,
    /// exactly as [`BulkSpec::write`] does: a span pass that generates each
    /// contig only long enough to learn its populated span before dropping
    /// its records, then a write pass that regenerates each contig (a pure
    /// function of `(seed, contig_idx, n_records)`, so byte-identical to the
    /// span pass) and writes it immediately, dropping its records before
    /// moving to the next contig. Peak memory here is therefore bounded by
    /// the largest single contig's records, not the sum across every
    /// contig — the same bound `write()` gives the real output, and the
    /// fix for this function previously holding every contig's full record
    /// set (`Vec<Vec<Rec>>`) live at once.
    ///
    /// The returned size is read back from `finish_and_index`'d file
    /// metadata, not a live byte counter (dispatch to the compression
    /// thread pool is asynchronous and would otherwise lag). Because
    /// `write()` also never calls a mid-stream flush, this temp-file write
    /// is structurally identical to the real write — same header, same
    /// records, same (absent) flush cadence — so the byte count returned
    /// here is exactly what the real write produces for this
    /// `per_contig_count`, not merely close to it.
    fn measure_compressed_bytes(
        &self,
        pool: &rayon::ThreadPool,
        samplers: &Samplers,
        fitted: &Fitted,
        per_contig_count: &[u64],
    ) -> Result<u64, BulkError> {
        let spans: Vec<u64> = self
            .contig_ids
            .iter()
            .zip(per_contig_count)
            .enumerate()
            .map(|(i, (id, &n))| {
                let recs = self.generate_contig(pool, samplers, fitted, id, i as u64, n);
                contig_span(&recs)
            })
            .collect();

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
        for (i, (id, &n)) in self.contig_ids.iter().zip(per_contig_count).enumerate() {
            let recs = self.generate_contig(pool, samplers, fitted, id, i as u64, n);
            for r in &recs {
                let buf = to_record_buf(&r.g, &self.payload, r.phased);
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

/// Parse a byte size like `100MB`, `512KB`, `1GB`, or a bare byte count.
pub fn parse_size(s: &str) -> Result<u64, BulkError> {
    let t = s.trim();
    let (num, mult) = if let Some(p) = t.strip_suffix("GB") {
        (p, 1024 * 1024 * 1024)
    } else if let Some(p) = t.strip_suffix("MB") {
        (p, 1024 * 1024)
    } else if let Some(p) = t.strip_suffix("KB") {
        (p, 1024)
    } else {
        (t, 1)
    };
    num.trim()
        .parse::<u64>()
        .map(|n| n * mult)
        .map_err(|_| BulkError::BadSize(s.to_string()))
}

/// Parse `--records-for` tokens (`NAME=COUNT`) into ordered `(name, count)`
/// pairs.
///
/// Returns a `Vec`, not the `BTreeMap` that [`Size::PerContig`] wants,
/// because command-line order is load-bearing: when `--contigs` is omitted
/// these names supply the output contig order, which a map would discard.
/// The caller collects into a map once it has taken the order it needs.
///
/// Lives here rather than in the binary for the same reason [`parse_size`]
/// does — `tests/cli.rs` cannot import from a `[[bin]]` target.
pub fn parse_records_for(tokens: &[String]) -> Result<Vec<(String, u64)>, BulkError> {
    let mut out: Vec<(String, u64)> = Vec::with_capacity(tokens.len());
    for tok in tokens {
        let (name, count) = tok.split_once('=').ok_or_else(|| {
            BulkError::BadRecordsFor(format!("expected NAME=COUNT in --records-for, got {tok:?}"))
        })?;
        let name = name.trim();
        if name.is_empty() {
            return Err(BulkError::BadRecordsFor(format!(
                "empty contig name in --records-for entry {tok:?}"
            )));
        }
        let count: u64 = count.trim().parse().map_err(|_| {
            BulkError::BadRecordsFor(format!(
                "expected a non-negative integer record count in --records-for \
                 entry {tok:?}"
            ))
        })?;
        // Linear scan: a run has tens of contigs, not thousands, and this
        // preserves the order a `BTreeMap`-based dedupe would lose.
        if out.iter().any(|(n, _)| n == name) {
            return Err(BulkError::BadRecordsFor(format!(
                "duplicate contig name {name:?} in --records-for"
            )));
        }
        out.push((name.to_string(), count));
    }
    Ok(out)
}

/// Moves `src` to `dst`: a rename when both are on the same filesystem
/// (the common case, and atomic), falling back to copy-then-remove when
/// they aren't (e.g. `TMPDIR` on a different filesystem than the output
/// directory, where `rename` returns `EXDEV`). Used by
/// [`BulkSpec::promote_temp`] for the `.csi` companion, which
/// `NamedTempFile::persist` doesn't know how to move itself.
fn move_file(src: &Path, dst: &Path) -> Result<(), BulkError> {
    if std::fs::rename(src, dst).is_ok() {
        return Ok(());
    }
    std::fs::copy(src, dst)?;
    let _ = std::fs::remove_file(src);
    Ok(())
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

/// Resolves [`Size::PerContig`] to per-contig counts parallel to
/// `contig_ids`, validating the two name sets against each other first.
///
/// Both directions are errors, and both name the offending contigs:
///
/// - A requested contig with no entry would otherwise generate a silently
///   empty contig — indistinguishable in the output from a deliberate
///   zero, but not what the caller asked for.
/// - An entry naming a contig that was not requested is almost always a
///   typo (`"1"` for `"chr1"`, `"chrX"` in an autosome-only run). Ignoring
///   it silently hands back a corpus that does not match the request.
///
/// Names are compared exactly. [`resolve_contig_stat`]'s `chr`-prefix
/// normalization exists to reconcile caller-chosen output names against the
/// bare ids committed profiles were fit from; that is a different problem.
/// Here both lists come from the same caller in the same call, so a loud
/// error beats a second, fuzzier way to spell a contig name.
fn per_contig_counts(
    map: &BTreeMap<String, u64>,
    contig_ids: &[String],
) -> Result<Vec<u64>, BulkError> {
    let missing: Vec<String> = contig_ids
        .iter()
        .filter(|id| !map.contains_key(id.as_str()))
        .map(|id| id.to_string())
        .collect();
    if !missing.is_empty() {
        return Err(BulkError::PerContigMissing(missing));
    }

    let unknown: Vec<String> = map
        .keys()
        .filter(|k| !contig_ids.iter().any(|id| id == *k))
        .map(|k| k.to_string())
        .collect();
    if !unknown.is_empty() {
        return Err(BulkError::PerContigUnknown(unknown));
    }

    // Indexing is guarded by the `missing` check above.
    Ok(contig_ids.iter().map(|id| map[id.as_str()]).collect())
}

/// Splits `total` records across `contig_ids` proportional to each
/// contig's fitted per-contig variant count (`n_variants`, via
/// [`resolve_contig_stat`]), using the largest-remainder method so the
/// per-contig counts sum to exactly `total`. Falls back to an even split if
/// every resolved weight is zero (a degenerate profile), rather than
/// dividing by zero.
///
/// Weighting by `n_variants` rather than `density_per_kb`: output density is
/// a *global* `1/mean(gap)` draw (`gap_dist` is not fit per contig), so a
/// per-contig fitted density is never actually reproduced by this split --
/// and an outlier contig (e.g. MT at ~350/kb, ~12x the rest) would skew the
/// split for a statistic the output doesn't even follow. `n_variants`, by
/// contrast, reproduces the source's real per-contig variant distribution.
fn distribute_by_n_variants(fitted: &Fitted, contig_ids: &[String], total: u64) -> Vec<u64> {
    let weights: Vec<f64> = contig_ids
        .iter()
        .enumerate()
        .map(|(i, id)| resolve_contig_stat(fitted, i, id).n_variants as f64)
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
/// Must stay in sync with `generate::to_record_buf`'s own (private) `key_names`
/// match — duplicated here rather than shared because `generate.rs` is out of
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
/// definition matching exactly what `generate::SampleStats::value_for` emits.
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

#[cfg(test)]
mod tests {
    use super::*;

    // A profile where n_variants order != density order, so the two split
    // strategies give different answers.
    const DISCRIMINATING_PROFILE: &str = r#"{
      "name": "disc", "provenance": {"source":"x","n_samples_source":10,
        "n_variants_source":11000,"fitted_on":"2026-01-01","fit_tool_version":"t",
        "supplied":["ploidy"]},
      "fitted": { "contigs": [
          { "id": "big", "n_variants": 10000, "density_per_kb": 10.0 },
          { "id": "small", "n_variants": 1000, "density_per_kb": 90.0 }
        ],
        "gap_dist": {"edges":[1.0,2.0],"weights":[1.0]},
        "sfs": {"edges":[1.0,2.0],"weights":[1.0]},
        "variant_classes": {"snp":1.0,"insertion":0.0,"deletion":0.0,"mnp":0.0,"complex":0.0,"symbolic":0.0},
        "indel_length": {"edges":[1.0,2.0],"weights":[1.0]},
        "titv": 2.0, "multiallelic_rate": 0.0, "missing_rate": 0.0, "phased_rate": 1.0
      },
      "dialed": { "payload": "gt-only", "ploidy": 2 }
    }"#;

    #[test]
    fn records_split_follows_n_variants_not_density() {
        let p = Profile::from_json(DISCRIMINATING_PROFILE).unwrap();
        let ids = vec!["big".to_string(), "small".to_string()];
        let counts = distribute_by_n_variants(&p.fitted, &ids, 11_000);
        assert_eq!(counts, vec![10_000, 1_000]);
    }

    /// The whole point of the split: only profile-content errors carry the
    /// `invalid profile:` prefix. Every other class names its own problem.
    #[test]
    fn error_messages_are_classified() {
        // Profile content keeps the prefix -- these messages name a profile
        // field, not themselves, so the prefix is load-bearing.
        assert_eq!(
            BulkError::InvalidProfile("ploidy must be >= 1".into()).to_string(),
            "invalid profile: ploidy must be >= 1"
        );
        assert!(BulkError::PayloadPloidy {
            payload: Payload::Gatk,
            ploidy: 3,
        }
        .to_string()
        .starts_with("invalid profile: payload Gatk emits AD and/or PL"));

        // Everything else must NOT claim the profile is at fault.
        for e in [
            BulkError::NoContigs,
            BulkError::NoSamples,
            BulkError::DuplicateContig("chr1".into()),
            BulkError::BadSize("banana".into()),
            BulkError::CompressionLevel(99),
            BulkError::TargetNotReached {
                target_bytes: 1024,
                corrections: 4,
            },
            BulkError::ProfileLoad {
                path: "x".into(),
                source: std::io::Error::from(std::io::ErrorKind::NotFound),
            },
            BulkError::BadRecordsFor("expected NAME=COUNT in --records-for, got \"x\"".into()),
            BulkError::PerContigMissing(vec!["chr2".into()]),
            BulkError::PerContigUnknown(vec!["1".into()]),
            BulkError::TooManyBlocks {
                contig: "chr1".into(),
                // Not `stride` itself: `n_blocks == stride` fits exactly
                // (see the `>` check in `compute_layouts`) and would render
                // the nonsensical "needs 1000000 blocks, which exceeds the
                // 1000000-block stride". One past it is the smallest value
                // that actually triggers the error, so it's what the
                // message under test should read.
                n_blocks: 1_000_001,
                stride: 1_000_000,
            },
        ] {
            let msg = e.to_string();
            assert!(
                !msg.starts_with("invalid profile:"),
                "{e:?} must not be described as an invalid profile, got: {msg}"
            );
        }

        assert_eq!(BulkError::NoContigs.to_string(), "need >= 1 output contig");
        assert_eq!(BulkError::NoSamples.to_string(), "need >= 1 sample");
        assert!(BulkError::DuplicateContig("chr1".into())
            .to_string()
            .starts_with(r#"duplicate output contig name: "chr1""#));
        assert!(BulkError::BadSize("banana".into())
            .to_string()
            .starts_with(r#"bad size: "banana""#));
        assert_eq!(
            BulkError::TargetNotReached {
                target_bytes: 1024,
                corrections: 4,
            }
            .to_string(),
            "could not reach target size 1024 bytes within 4 corrective rounds"
        );
    }

    #[test]
    fn block_records_is_bounded_by_cells_and_by_the_record_cap() {
        // Narrow cohorts saturate at the record cap.
        assert_eq!(BulkSpec::block_records(1, 2), BulkSpec::MAX_BLOCK_RECORDS);
        assert_eq!(
            BulkSpec::block_records(4_000, 2),
            BulkSpec::MAX_BLOCK_RECORDS
        );
        // Wide cohorts are bounded by cells: 4e6 / (32000*2) = 62.
        assert_eq!(BulkSpec::block_records(32_000, 2), 62);
        // Never zero, however wide.
        assert_eq!(BulkSpec::block_records(100_000_000, 2), 1);
    }

    #[test]
    fn contig_layout_arithmetic_matches_block_records_and_offsets() {
        // `n_records == 0`: no blocks at all, and nothing to sum or offset.
        let empty = ContigLayout {
            block_records: 100,
            n_records: 0,
            block_spans: vec![],
        };
        assert_eq!(empty.n_blocks(), 0);
        assert_eq!(empty.span(), 0);
        assert_eq!(empty.offsets(), Vec::<u64>::new());

        // Exact multiple of `block_records`: every block, including the
        // last, is full-width -- no partial tail to get wrong.
        let exact = ContigLayout {
            block_records: 100,
            n_records: 300,
            block_spans: vec![1_000, 2_000, 3_000],
        };
        assert_eq!(exact.n_blocks(), 3);
        assert_eq!(exact.block_len(0), 100);
        assert_eq!(exact.block_len(1), 100);
        assert_eq!(exact.block_len(2), 100);
        assert_eq!(exact.span(), 6_000);
        assert_eq!(exact.offsets(), vec![0, 1_000, 3_000]);

        // Partial tail: 250 records at 100/block is 100, 100, 50 -- the
        // case `block_len`'s `.min(n_records - start)` exists for.
        let partial = ContigLayout {
            block_records: 100,
            n_records: 250,
            block_spans: vec![500, 700, 300],
        };
        assert_eq!(partial.n_blocks(), 3);
        assert_eq!(partial.block_len(0), 100);
        assert_eq!(partial.block_len(1), 100);
        assert_eq!(partial.block_len(2), 50);
        assert_eq!(partial.span(), 1_500);
        assert_eq!(partial.offsets(), vec![0, 500, 1_200]);
    }
}
