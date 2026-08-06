//! Streaming BCF/VCF(.gz) output: a multithreaded bgzf writer plus a
//! second-pass CSI index.
//!
//! See `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`,
//! "Writer and sizing": `MultithreadedWriter` exposes neither
//! `virtual_position()` nor a byte count, so the index is a second pass via
//! [`noodles_bcf::fs::index`]. [`Size::Target`](crate::bulk::Size::Target)
//! measures compressed output size via `fs::metadata` on a finished temp
//! file (see [`crate::bulk::BulkSpec::measured_bytes`]), not a
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

/// Encodes a block's records into an in-memory buffer, off the writer
/// thread.
///
/// One encoder is built per rayon worker (not per block) and its buffer is
/// reused: `noodles_bcf::io::Writer` derives its `StringMaps` inside
/// `write_header` and keeps it privately, so a header-less writer fails at
/// runtime with "chromosome not in string map". The header is therefore
/// written once to populate that map, `header_len` is remembered, and each
/// block rewinds to it with `truncate`. At 32,000 samples the header text
/// is ~200 KB of sample names, which is exactly why this is per worker
/// rather than per block.
///
/// Driven by [`crate::bulk::BulkSpec`]'s block pipeline, whose serial
/// consumer hands the resulting bytes to [`BulkWriter::write_encoded`].
pub(crate) enum BlockEncoder {
    Bcf {
        w: bcf::io::Writer<Vec<u8>>,
        header_len: usize,
    },
    Text {
        w: vcf::io::Writer<Vec<u8>>,
        header_len: usize,
    },
}

impl BlockEncoder {
    /// Builds an encoder whose output matches what [`BulkWriter`] with the
    /// same `format` and `header` would emit. `VcfGz` and `Vcf` share the
    /// same text encoding; they differ only in whether [`BulkWriter`]'s
    /// sink compresses.
    pub(crate) fn new(format: Format, header: &vcf::Header) -> Result<BlockEncoder, BulkError> {
        match format {
            Format::Bcf => {
                let mut w = bcf::io::Writer::from(Vec::new());
                w.write_header(header)?;
                let header_len = w.get_ref().len();
                Ok(BlockEncoder::Bcf { w, header_len })
            }
            Format::VcfGz | Format::Vcf => {
                let mut w = vcf::io::Writer::new(Vec::new());
                w.write_header(header)?;
                let header_len = w.get_ref().len();
                Ok(BlockEncoder::Text { w, header_len })
            }
        }
    }

    /// Starts a new block, discarding the previous one's bytes while
    /// keeping the string map the header write populated.
    pub(crate) fn begin(&mut self) {
        match self {
            BlockEncoder::Bcf { w, header_len } => w.get_mut().truncate(*header_len),
            BlockEncoder::Text { w, header_len } => w.get_mut().truncate(*header_len),
        }
    }

    /// Appends one record to the current block.
    pub(crate) fn push(
        &mut self,
        header: &vcf::Header,
        record: &RecordBuf,
    ) -> Result<(), BulkError> {
        match self {
            BlockEncoder::Bcf { w, .. } => w.write_variant_record(header, record)?,
            BlockEncoder::Text { w, .. } => w.write_variant_record(header, record)?,
        }
        Ok(())
    }

    /// The current block's encoded records, without the header prefix.
    pub(crate) fn bytes(&self) -> &[u8] {
        match self {
            BlockEncoder::Bcf { w, header_len } => &w.get_ref()[*header_len..],
            BlockEncoder::Text { w, header_len } => &w.get_ref()[*header_len..],
        }
    }
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

    /// Writes an already-encoded block of records straight to the sink.
    ///
    /// The counterpart to [`BlockEncoder`]: record encoding happens in a
    /// rayon worker, and the writer thread only concatenates the resulting
    /// bytes in block order. Like [`BulkWriter::write`], this does not
    /// force a flush — the bgzf writer batches into ~64 KiB blocks itself.
    pub fn write_encoded(&mut self, bytes: &[u8]) -> Result<(), BulkError> {
        use std::io::Write as _;
        match &mut self.sink {
            Sink::Bcf(x) => x.get_mut().write_all(bytes)?,
            Sink::VcfGz(x) => x.get_mut().write_all(bytes)?,
            Sink::Vcf(x) => x.get_mut().write_all(bytes)?,
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
            // Reserved GT FORMAT: `sample_records` (below) sets a "GT" value
            // on each record, and the BCF string map (built from the header
            // by `write_header`) must know about "GT" up front or encoding
            // fails at runtime with "genotype key not in string map".
            .add_format(
                vcf::variant::record::samples::keys::key::GENOTYPE,
                vcf::header::record::value::Map::<vcf::header::record::value::map::Format>::from(
                    vcf::variant::record::samples::keys::key::GENOTYPE,
                ),
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

    /// The range the message advertises must be the range actually enforced,
    /// checked from both sides of the boundary.
    ///
    /// Asserting only that the maximum is accepted would be a test that
    /// cannot fail: `create` hands the level straight to noodles, and
    /// `BEST` is by construction within noodles' own bounds. The claim worth
    /// pinning is the one that *can* drift -- that the bound named in the
    /// message is the bound the code applies. Hard-coding `0-9` in the
    /// message while a downstream feature-unified build accepted up to 12
    /// would fail here, and so would going back to noodles' wording.
    #[test]
    fn the_advertised_range_is_the_range_enforced() {
        let dir = tempfile::tempdir().unwrap();
        let max = bgzf::io::writer::CompressionLevel::BEST.get();

        // The top of the advertised range is accepted ...
        assert!(
            BulkWriter::create(
                &dir.path().join("max.bcf"),
                Format::Bcf,
                &header(),
                max,
                NonZero::new(1).unwrap(),
            )
            .is_ok(),
            "level {max} is advertised as valid and must be accepted"
        );

        // ... and one past it is rejected, by a message naming that same
        // bound. This is what ties the advertised range to the real one.
        let over = max + 1;
        // Not `expect_err`: `BulkWriter` is deliberately not `Debug`.
        let err = match BulkWriter::create(
            &dir.path().join("over.bcf"),
            Format::Bcf,
            &header(),
            over,
            NonZero::new(1).unwrap(),
        ) {
            Ok(_) => panic!("level {over} is past the advertised maximum and must be rejected"),
            Err(e) => e,
        };
        assert_eq!(
            err.to_string(),
            format!("invalid compression level: {over} (expected 0-{max})")
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

    fn sample_records(n: usize) -> Vec<RecordBuf> {
        (1..=n)
            .map(|pos| {
                let keys: vcf::variant::record_buf::samples::Keys =
                    ["GT"].iter().map(|k| k.to_string()).collect();
                let values = vec![
                    vec![Some(
                        vcf::variant::record_buf::samples::sample::Value::from("0|1".to_string(),)
                    )];
                    1
                ];
                RecordBuf::builder()
                    .set_reference_sequence_name("chr1")
                    .set_variant_start(noodles_core::Position::try_from(pos * 10).unwrap())
                    .set_reference_bases("A")
                    .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![
                        String::from("T"),
                    ]))
                    .set_samples(vcf::variant::record_buf::Samples::new(keys, values))
                    .build()
            })
            .collect()
    }

    /// The load-bearing assumption of the block pipeline: bytes produced by
    /// a reused `BlockEncoder` are exactly what the real writer would have
    /// emitted for the same records in the same order.
    #[test]
    fn block_encoded_bytes_match_a_direct_write() {
        for format in [Format::Bcf, Format::VcfGz, Format::Vcf] {
            let h = header();
            let records = sample_records(10);

            // Direct: every record straight through BulkWriter.
            let dir = tempfile::tempdir().unwrap();
            let direct_path = dir.path().join("direct.out");
            let mut w =
                BulkWriter::create(&direct_path, format, &h, 6, NonZero::new(1).unwrap()).unwrap();
            for r in &records {
                w.write(&h, r).unwrap();
            }
            w.finish_and_index(&direct_path).unwrap();

            // Blocked: encode in blocks of 3 through one reused encoder,
            // hand the bytes to BulkWriter.
            let blocked_path = dir.path().join("blocked.out");
            let mut w =
                BulkWriter::create(&blocked_path, format, &h, 6, NonZero::new(1).unwrap()).unwrap();
            let mut enc = BlockEncoder::new(format, &h).unwrap();
            for chunk in records.chunks(3) {
                enc.begin();
                for r in chunk {
                    enc.push(&h, r).unwrap();
                }
                w.write_encoded(enc.bytes()).unwrap();
            }
            w.finish_and_index(&blocked_path).unwrap();

            assert_eq!(
                std::fs::read(&direct_path).unwrap(),
                std::fs::read(&blocked_path).unwrap(),
                "block-encoded output must be byte-identical for {format:?}"
            );
        }
    }

    /// `begin()` must discard the previous block, not append to it.
    #[test]
    fn begin_rewinds_the_block_buffer() {
        let h = header();
        let records = sample_records(4);
        let mut enc = BlockEncoder::new(Format::Bcf, &h).unwrap();

        enc.begin();
        enc.push(&h, &records[0]).unwrap();
        let first = enc.bytes().to_vec();

        enc.begin();
        enc.push(&h, &records[0]).unwrap();
        assert_eq!(enc.bytes(), first.as_slice(), "begin must reset the buffer");

        enc.begin();
        assert!(enc.bytes().is_empty(), "begin must leave an empty block");
    }
}
