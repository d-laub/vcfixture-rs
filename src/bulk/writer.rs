//! Streaming BCF/VCF(.gz) output: a multithreaded bgzf writer plus a
//! second-pass CSI index.
//!
//! See `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`,
//! "Writer and sizing": `MultithreadedWriter` exposes neither
//! `virtual_position()` nor a byte count, so the index is a second pass via
//! [`noodles_bcf::fs::index`]. [`Size::Target`](crate::bulk::Size::Target)
//! measures compressed output size via `fs::metadata` on a finished temp
//! file (see [`crate::bulk::BulkSpec::measure_compressed_bytes`]), not a
//! live byte counter — an earlier live `CountingWriter` counter was dead
//! code outside this module's own tests and has been removed.

use std::num::NonZero;
use std::path::Path;

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_csi as csi;
use noodles_vcf::{self as vcf, variant::io::Write as _, variant::RecordBuf};

use crate::bulk::BulkError;

/// Output container format for bulk generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Bcf,
    VcfGz,
    Vcf,
}

enum Sink {
    Bcf(bcf::io::Writer<bgzf::io::MultithreadedWriter<std::fs::File>>),
    VcfGz(vcf::io::Writer<bgzf::io::MultithreadedWriter<std::fs::File>>),
    Vcf(vcf::io::Writer<std::fs::File>),
}

/// A streaming BCF/VCF(.gz) writer.
///
/// Indexing is a second pass ([`BulkWriter::finish_and_index`]): the
/// multithreaded bgzf writer underneath exposes no virtual position, so a
/// CSI cannot be built incrementally while writing.
pub struct BulkWriter {
    sink: Sink,
    format: Format,
}

impl BulkWriter {
    /// Creates a new writer at `path`, writing `header` immediately.
    ///
    /// `workers` sizes the bgzf compression thread pool for `Bcf`/`VcfGz`;
    /// it is ignored for uncompressed `Vcf`. `compression_level` is
    /// validated unconditionally, even for `Vcf` (where it goes unused): an
    /// out-of-range level is rejected the same way regardless of format
    /// rather than only surfacing once a caller switches to a compressed
    /// one.
    pub fn create(
        path: &Path,
        format: Format,
        header: &vcf::Header,
        compression_level: u8,
        workers: NonZero<usize>,
    ) -> Result<BulkWriter, BulkError> {
        let file = std::fs::File::create(path)?;

        // The rejected level, not noodles' rendering of it: `BulkError` names
        // the range itself, so the foreign error adds nothing here.
        let level = bgzf::io::writer::CompressionLevel::try_from(compression_level)
            .map_err(|_| BulkError::CompressionLevel(compression_level))?;

        let mut w = match format {
            Format::Bcf => {
                let inner = bgzf::io::multithreaded_writer::Builder::default()
                    .set_worker_count(workers)
                    .set_compression_level(level)
                    .build_from_writer(file);
                BulkWriter {
                    sink: Sink::Bcf(bcf::io::Writer::from(inner)),
                    format,
                }
            }
            Format::VcfGz => {
                let inner = bgzf::io::multithreaded_writer::Builder::default()
                    .set_worker_count(workers)
                    .set_compression_level(level)
                    .build_from_writer(file);
                BulkWriter {
                    sink: Sink::VcfGz(vcf::io::Writer::new(inner)),
                    format,
                }
            }
            Format::Vcf => BulkWriter {
                sink: Sink::Vcf(vcf::io::Writer::new(file)),
                format,
            },
        };

        match &mut w.sink {
            Sink::Bcf(x) => x.write_header(header)?,
            Sink::VcfGz(x) => x.write_header(header)?,
            Sink::Vcf(x) => x.write_header(header)?,
        }

        Ok(w)
    }

    /// Writes one record.
    ///
    /// This does not force a flush: the bgzf writer batches uncompressed
    /// data into ~64 KiB blocks before dispatching them to the compression
    /// thread pool, which is what the multithreaded throughput depends on.
    pub fn write(&mut self, header: &vcf::Header, record: &RecordBuf) -> Result<(), BulkError> {
        match &mut self.sink {
            Sink::Bcf(x) => x.write_variant_record(header, record)?,
            Sink::VcfGz(x) => x.write_variant_record(header, record)?,
            Sink::Vcf(x) => x.write_variant_record(header, record)?,
        }
        Ok(())
    }

    /// Finishes the stream, then writes `<path>.csi` via a second read pass.
    ///
    /// Indexing only applies to `Bcf`; `VcfGz`/`Vcf` are not indexed here.
    pub fn finish_and_index(self, path: &Path) -> Result<(), BulkError> {
        let format = self.format;

        // Dropping the sink drops the underlying `MultithreadedWriter`, whose
        // `Drop` impl calls `finish()` (flushing remaining buffers, joining
        // the compression/writer threads, and appending the BGZF EOF block).
        // For the uncompressed `Vcf` sink, dropping the inner `File` simply
        // flushes it on close.
        drop(self.sink);

        if format == Format::Bcf {
            let index = bcf::fs::index(path)?;
            let mut csi_path = path.as_os_str().to_os_string();
            csi_path.push(".csi");
            csi::fs::write(std::path::PathBuf::from(csi_path), &index)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_vcf::{self as vcf, variant::RecordBuf};
    use std::num::NonZero;

    fn header() -> vcf::Header {
        vcf::Header::builder()
            .add_contig(
                "chr1",
                vcf::header::record::value::Map::<vcf::header::record::value::map::Contig>::new(),
            )
            .add_sample_name("s1")
            .build()
    }

    /// An out-of-range bgzf compression level is an argument error. It is
    /// the caller's number that is wrong, not the profile.
    ///
    /// The rendered message is asserted, not just the variant: the whole
    /// point of carrying the `u8` is that the user sees their own number and
    /// the accepted range, instead of noodles' `invalid input: 99` doubled
    /// behind our own prefix.
    #[test]
    fn out_of_range_compression_level_is_an_argument_error() {
        let dir = tempfile::tempdir().unwrap();
        // `header()` is the existing helper at the top of this test module.
        let result = BulkWriter::create(
            &dir.path().join("a.bcf"),
            Format::Bcf,
            &header(),
            99,
            NonZero::new(1).unwrap(),
        );
        assert!(
            matches!(result, Err(BulkError::CompressionLevel(99))),
            "compression level 99 is out of range and must be an argument \
             error carrying the offending level, not an invalid profile"
        );

        let msg = result.err().unwrap().to_string();
        assert!(!msg.starts_with("invalid profile:"));
        // The upper bound is whatever the linked noodles-bgzf accepts, so it
        // is read the same way `BulkError` renders it rather than written out
        // here -- noodles' `libdeflate` feature moves it from 9 to 12.
        let max = bgzf::io::writer::CompressionLevel::BEST.get();
        assert_eq!(
            msg,
            format!("invalid compression level: 99 (expected 0-{max})")
        );
        assert!(
            !msg.contains("invalid input"),
            "must not embed noodles' error wording, got: {msg}"
        );
    }

    /// The largest level noodles accepts must be accepted here too -- the
    /// boundary the error message advertises has to be real, or the message
    /// sends users at a level that then fails.
    #[test]
    fn the_advertised_maximum_compression_level_is_accepted() {
        let dir = tempfile::tempdir().unwrap();
        let max = bgzf::io::writer::CompressionLevel::BEST.get();
        let result = BulkWriter::create(
            &dir.path().join("max.bcf"),
            Format::Bcf,
            &header(),
            max,
            NonZero::new(1).unwrap(),
        );
        assert!(
            result.is_ok(),
            "level {max} is advertised as valid and must be accepted"
        );
    }

    #[test]
    fn writes_a_readable_indexed_bcf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bcf");
        let h = header();
        let mut w =
            BulkWriter::create(&path, Format::Bcf, &h, 6, NonZero::new(2).unwrap()).unwrap();
        for pos in [100usize, 200, 300] {
            let rec = RecordBuf::builder()
                .set_reference_sequence_name("chr1")
                .set_variant_start(noodles_core::Position::try_from(pos).unwrap())
                .set_reference_bases("A")
                .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![
                    String::from("T"),
                ]))
                .build();
            w.write(&h, &rec).unwrap();
        }

        // `finish_and_index` drops `self.sink`, whose `MultithreadedWriter`
        // `Drop` impl calls `finish()`, which synchronously joins the
        // deflater/writer threads before returning — so by the time this
        // returns, every byte is guaranteed to be on disk (no polling an
        // async background thread needed).
        w.finish_and_index(&path).unwrap();

        assert!(path.exists());
        assert!(
            std::fs::metadata(&path).unwrap().len() > 0,
            "output file must be non-empty"
        );
        assert!(
            path.with_extension("bcf.csi").exists(),
            "csi must be written"
        );

        // Read back through an independent path and confirm the records survive.
        let mut r = noodles_bcf::io::reader::Builder::default()
            .build_from_path(&path)
            .unwrap();
        let rh = r.read_header().unwrap();
        let n = r.records().count();
        assert_eq!(n, 3);
        assert_eq!(rh.sample_names().len(), 1);
    }

    #[test]
    fn output_is_byte_identical_regardless_of_worker_count() {
        // `MultithreadedWriter` only has something to *reorder* once a run
        // spans more than one bgzf block: its staging buffer dispatches a
        // block to the compression/writer thread pool once it holds
        // `MAX_BUF_SIZE` (~65,498) uncompressed bytes. A payload under that
        // threshold becomes a single block no matter how many workers are
        // configured, in which case this test would pass even if a real
        // reordering bug existed. The record count and padded ALT allele
        // below are sized to comfortably clear several block boundaries —
        // the payload-size assertion further down pins that down so the
        // test can't be silently shrunk back into vacuity.
        const RECORDS: usize = 1_000;
        const ALT_PAD_LEN: usize = 400;
        const MAX_BUF_SIZE: usize = 65_498;

        fn records(alt: &str) -> Vec<RecordBuf> {
            (1..=RECORDS)
                .map(|pos| {
                    RecordBuf::builder()
                        .set_reference_sequence_name("chr1")
                        .set_variant_start(noodles_core::Position::try_from(pos).unwrap())
                        .set_reference_bases("A")
                        .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![
                            alt.to_string(),
                        ]))
                        .build()
                })
                .collect()
        }

        fn run(workers: usize, h: &vcf::Header, recs: &[RecordBuf]) -> Vec<u8> {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("t.bcf");
            let mut w =
                BulkWriter::create(&path, Format::Bcf, h, 6, NonZero::new(workers).unwrap())
                    .unwrap();
            for rec in recs {
                w.write(h, rec).unwrap();
            }
            w.finish_and_index(&path).unwrap();
            std::fs::read(&path).unwrap()
        }

        let h = header();
        let alt = "T".repeat(ALT_PAD_LEN);
        let recs = records(&alt);

        // Measure the true *uncompressed* record payload directly (rather
        // than trusting arithmetic on record count and padding length),
        // independent of how well bgzf's DEFLATE happens to compress this
        // repetitive test data, by encoding the same records into an
        // in-memory, uncompressed BCF.
        let mut uncompressed = bcf::io::Writer::from(Vec::new());
        uncompressed.write_header(&h).unwrap();
        for rec in &recs {
            uncompressed.write_variant_record(&h, rec).unwrap();
        }
        let uncompressed_payload = uncompressed.get_ref().len();
        assert!(
            uncompressed_payload > 3 * MAX_BUF_SIZE,
            "test payload ({uncompressed_payload} bytes) must exceed several bgzf \
             blocks ({MAX_BUF_SIZE} bytes each) or this test cannot detect reordering bugs"
        );

        let a = run(1, &h, &recs);
        let b = run(4, &h, &recs);
        let c = run(16, &h, &recs);
        assert_eq!(
            a, b,
            "output must be byte-identical regardless of thread count (1 vs 4 workers)"
        );
        assert_eq!(
            a, c,
            "output must be byte-identical regardless of thread count (1 vs 16 workers)"
        );
        assert!(!a.is_empty());
    }
}
