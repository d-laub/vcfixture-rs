//! `vcfixture` CLI: currently just the `bulk` subcommand (benchmark-scale
//! BCF/VCF generation). See `vcfixture::bulk` for the library API this
//! wraps, and `docs/book/src/bulk-generation.md` for the user-facing guide.

use std::num::NonZero;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use clap::{Parser, Subcommand, ValueEnum};

use vcfixture::bulk::{parse_size, BulkError, BulkSpec, Format, Payload, Profile, Size};

#[derive(Parser)]
#[command(name = "vcfixture", about = "Generate VCF/BCF test data")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Generate a bulk BCF for benchmarking.
    Bulk {
        /// Builtin profile name (germline-1kgp, germline-1kgp-unphased,
        /// somatic-gdc), or a path to a profile JSON.
        #[arg(long, default_value = "germline-1kgp")]
        profile: String,
        /// Number of samples. Sample names are generated as `s0..s{n-1}`.
        /// The profile's fitted site-frequency spectrum is measured against
        /// its own source cohort size (`provenance.n_samples_source`, e.g.
        /// 3202 for the 1kGP profiles); requesting a different sample count
        /// here does not change that -- each drawn allele count is rescaled
        /// to a frequency against the source cohort and reapplied to this
        /// run's cohort size, so alt-allele density stays realistic
        /// regardless of how many samples you ask for.
        #[arg(long, default_value_t = 1)]
        samples: usize,
        /// Output contig names, in the order they will be written.
        #[arg(
            long,
            value_delimiter = ',',
            default_values_t = ["chr1".to_string(), "chr2".to_string(), "chr3".to_string()]
        )]
        contigs: Vec<String>,
        /// Stop once the compressed output reaches this size (e.g. 100MB).
        #[arg(long, conflicts_with_all = ["records", "records_per_contig"])]
        target_size: Option<String>,
        /// Generate exactly this many records total, split across contigs
        /// proportional to the profile's fitted density.
        #[arg(long, conflicts_with_all = ["target_size", "records_per_contig"])]
        records: Option<u64>,
        /// Generate exactly this many records for each requested contig.
        #[arg(long, conflicts_with_all = ["target_size", "records"])]
        records_per_contig: Option<u64>,
        /// Override the profile's payload preset (which FORMAT fields to
        /// synthesize).
        #[arg(long)]
        payload: Option<PayloadArg>,
        /// Output container format.
        #[arg(long, default_value = "bcf")]
        format: FormatArg,
        /// PRNG seed. Same seed + profile + spec produces byte-identical
        /// output, regardless of thread count.
        #[arg(long, default_value_t = 0)]
        seed: u64,
        /// bgzf compression level (0-9).
        #[arg(long, default_value_t = 6)]
        compression_level: u8,
        /// Worker thread count for compression/generation. Defaults to all
        /// available cores.
        #[arg(long, value_parser = parse_threads)]
        threads: Option<NonZero<usize>>,
        /// Output path (e.g. `out.bcf`). A `.csi` index and a
        /// `.summary.json` are written alongside it.
        #[arg(short, long)]
        output: PathBuf,
    },
}

#[derive(Clone, Copy, ValueEnum)]
enum PayloadArg {
    GtOnly,
    GtVaf,
    Gatk,
    Mutect2,
}

impl From<PayloadArg> for Payload {
    fn from(p: PayloadArg) -> Payload {
        match p {
            PayloadArg::GtOnly => Payload::GtOnly,
            PayloadArg::GtVaf => Payload::GtVaf,
            PayloadArg::Gatk => Payload::Gatk,
            PayloadArg::Mutect2 => Payload::Mutect2,
        }
    }
}

#[derive(Clone, Copy, ValueEnum)]
enum FormatArg {
    Bcf,
    VcfGz,
    Vcf,
}

impl From<FormatArg> for Format {
    fn from(f: FormatArg) -> Format {
        match f {
            FormatArg::Bcf => Format::Bcf,
            FormatArg::VcfGz => Format::VcfGz,
            FormatArg::Vcf => Format::Vcf,
        }
    }
}

/// Resolves `--profile` as a builtin name first, falling back to reading it
/// as a path to a profile JSON file.
fn resolve_profile(name_or_path: &str) -> Result<Profile, BulkError> {
    match Profile::builtin(name_or_path) {
        Ok(p) => Ok(p),
        Err(BulkError::UnknownProfile(_)) => {
            let text =
                std::fs::read_to_string(name_or_path).map_err(|e| BulkError::ProfileLoad {
                    path: name_or_path.to_string(),
                    source: e,
                })?;
            Profile::from_json(&text)
        }
        Err(e) => Err(e),
    }
}

/// Parses `--threads` straight into a `NonZero<usize>` so zero is rejected
/// by clap's own usage error and never becomes a library error variant.
fn parse_threads(s: &str) -> Result<NonZero<usize>, String> {
    let n: usize = s
        .parse()
        .map_err(|_| format!("{s:?} is not a thread count"))?;
    NonZero::new(n).ok_or_else(|| "must be >= 1".to_string())
}

fn run() -> Result<(), BulkError> {
    let cli = Cli::parse();
    let Cmd::Bulk {
        profile,
        samples,
        contigs,
        target_size,
        records,
        records_per_contig,
        payload,
        format,
        seed,
        compression_level,
        threads,
        output,
    } = cli.cmd;

    let profile = resolve_profile(&profile)?;

    let size = if let Some(s) = target_size {
        Size::Target(parse_size(&s)?)
    } else if let Some(n) = records {
        Size::Records(n)
    } else if let Some(n) = records_per_contig {
        Size::RecordsPerContig(n)
    } else {
        Size::RecordsPerContig(1000)
    };

    let mut spec = BulkSpec::new(profile)
        .samples(samples)
        .contigs(contigs)
        .size(size)
        .format(format.into())
        .seed(seed)
        .compression_level(compression_level);

    if let Some(p) = payload {
        spec = spec.payload(p.into());
    }
    if let Some(n) = threads {
        spec = spec.workers(n);
    }

    let start = Instant::now();
    let summary = spec.write(&output)?;
    let elapsed = start.elapsed();

    let bytes = std::fs::metadata(&output).map(|m| m.len()).unwrap_or(0);
    eprintln!(
        "wrote {} ({} bytes, {} records) in {:.2?}",
        output.display(),
        bytes,
        summary.n_records_total(),
        elapsed
    );

    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
