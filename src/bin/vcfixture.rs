//! `vcfixture` CLI: currently just the `bulk` subcommand (benchmark-scale
//! BCF/VCF generation). See `vcfixture::bulk` for the library API this
//! wraps, and `docs/book/src/bulk-generation.md` for the user-facing guide.

use std::num::NonZero;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use clap::{Parser, Subcommand, ValueEnum};

use vcfixture::bulk::{
    parse_records_for, parse_size, BulkError, BulkSpec, Format, Payload, Profile, Size,
};

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
        /// Defaults to the `--records-for` names in the order given, or to
        /// `chr1,chr2,chr3` when neither is supplied.
        #[arg(long, value_delimiter = ',')]
        contigs: Option<Vec<String>>,
        /// Stop once the compressed output reaches this size (e.g. 100MB).
        #[arg(long, conflicts_with_all = ["records", "records_per_contig", "records_for"])]
        target_size: Option<String>,
        /// Generate exactly this many records total, split across contigs
        /// proportional to the profile's fitted per-contig variant counts
        /// (`n_variants`) -- not to its fitted density, which is close to
        /// uniform across the human autosomes.
        #[arg(long, conflicts_with_all = ["target_size", "records_per_contig", "records_for"])]
        records: Option<u64>,
        /// Generate exactly this many records for each requested contig.
        #[arg(long, conflicts_with_all = ["target_size", "records", "records_for"])]
        records_per_contig: Option<u64>,
        /// Generate exactly this many records for the named contig, e.g.
        /// `--records-for chr1=5759060,chr2=6088598`. Repeatable and
        /// comma-separated. When `--contigs` is omitted these names also
        /// set the output contig list and its order.
        #[arg(
            long,
            value_delimiter = ',',
            conflicts_with_all = ["target_size", "records", "records_per_contig"]
        )]
        records_for: Vec<String>,
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
        #[arg(long)]
        threads: Option<usize>,
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
            let text = std::fs::read_to_string(name_or_path).map_err(|e| {
                BulkError::Invalid(format!(
                    "profile {name_or_path:?} is not a builtin name and could not be \
                     read as a file: {e}"
                ))
            })?;
            Profile::from_json(&text)
        }
        Err(e) => Err(e),
    }
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
        records_for,
        payload,
        format,
        seed,
        compression_level,
        threads,
        output,
    } = cli.cmd;

    let profile = resolve_profile(&profile)?;

    let records_for = if records_for.is_empty() {
        None
    } else {
        Some(parse_records_for(&records_for)?)
    };

    // `--contigs` wins; otherwise `--records-for`'s names supply the output
    // order (so a 22-autosome run spells them out once, not twice, in two
    // places that could drift); otherwise the historical default.
    let contigs = match (contigs, &records_for) {
        (Some(c), _) => c,
        (None, Some(pairs)) => pairs.iter().map(|(id, _)| id.clone()).collect(),
        (None, None) => vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
    };

    let size = if let Some(s) = target_size {
        Size::Target(parse_size(&s)?)
    } else if let Some(n) = records {
        Size::Records(n)
    } else if let Some(n) = records_per_contig {
        Size::RecordsPerContig(n)
    } else if let Some(pairs) = records_for {
        // If `--contigs` was also given and disagrees with these keys,
        // `BulkSpec::write`'s `Size::PerContig` validation reports it with
        // the offending names -- no separate check needed here.
        Size::PerContig(pairs.into_iter().collect())
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
        let n = NonZero::new(n)
            .ok_or_else(|| BulkError::Invalid("--threads must be >= 1".to_string()))?;
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
