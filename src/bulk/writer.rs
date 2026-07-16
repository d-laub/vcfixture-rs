//! Streaming BCF/VCF(.gz) output: a byte-counting writer underneath a
//! multithreaded bgzf writer, and a second-pass CSI index.
//!
//! See `docs/superpowers/specs/2026-07-16-bulk-generation-design.md`,
//! "Writer and sizing": `MultithreadedWriter` exposes neither
//! `virtual_position()` nor a byte count, so the index is a second pass via
//! [`noodles_bcf::fs::index`], and the compressed byte count comes from a
//! [`CountingWriter`] placed *underneath* the bgzf writer.

use std::io::Write;
use std::num::NonZero;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_csi as csi;
use noodles_vcf::{self as vcf, variant::io::Write as _, variant::RecordBuf};

use crate::bulk::BulkError;

/// Wraps a writer and counts bytes passing through it.
///
/// Placed *underneath* the bgzf writer, so the count is the compressed size —
/// `MultithreadedWriter` exposes no position of its own.
pub struct CountingWriter<W> {
    inner: W,
    count: Arc<AtomicU64>,
}

impl<W: Write> CountingWriter<W> {
    pub fn new(inner: W) -> (CountingWriter<W>, Arc<AtomicU64>) {
        let count = Arc::new(AtomicU64::new(0));
        (
            CountingWriter {
                inner,
                count: Arc::clone(&count),
            },
            count,
        )
    }
}

impl<W: Write> Write for CountingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.count.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// Output container format for bulk generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Bcf,
    VcfGz,
    Vcf,
}

enum Sink {
    Bcf(bcf::io::Writer<bgzf::io::MultithreadedWriter<CountingWriter<std::fs::File>>>),
    VcfGz(vcf::io::Writer<bgzf::io::MultithreadedWriter<CountingWriter<std::fs::File>>>),
    Vcf(vcf::io::Writer<CountingWriter<std::fs::File>>),
}

/// A streaming BCF/VCF(.gz) writer with a live compressed-byte counter.
///
/// Indexing is a second pass ([`BulkWriter::finish_and_index`]): the
/// multithreaded bgzf writer underneath exposes no virtual position, so a
/// CSI cannot be built incrementally while writing.
pub struct BulkWriter {
    sink: Sink,
    count: Arc<AtomicU64>,
    format: Format,
}

impl BulkWriter {
    /// Creates a new writer at `path`, writing `header` immediately.
    ///
    /// `workers` sizes the bgzf compression thread pool for `Bcf`/`VcfGz`;
    /// it is ignored for uncompressed `Vcf`.
    pub fn create(
        path: &Path,
        format: Format,
        header: &vcf::Header,
        compression_level: u8,
        workers: NonZero<usize>,
    ) -> Result<BulkWriter, BulkError> {
        let file = std::fs::File::create(path)?;
        let (counting, count) = CountingWriter::new(file);

        let level = bgzf::io::writer::CompressionLevel::try_from(compression_level)
            .map_err(|e| BulkError::Invalid(format!("invalid compression level: {e}")))?;

        let mut w = match format {
            Format::Bcf => {
                let inner = bgzf::io::multithreaded_writer::Builder::default()
                    .set_worker_count(workers)
                    .set_compression_level(level)
                    .build_from_writer(counting);
                BulkWriter {
                    sink: Sink::Bcf(bcf::io::Writer::from(inner)),
                    count,
                    format,
                }
            }
            Format::VcfGz => {
                let inner = bgzf::io::multithreaded_writer::Builder::default()
                    .set_worker_count(workers)
                    .set_compression_level(level)
                    .build_from_writer(counting);
                BulkWriter {
                    sink: Sink::VcfGz(vcf::io::Writer::new(inner)),
                    count,
                    format,
                }
            }
            Format::Vcf => BulkWriter {
                sink: Sink::Vcf(vcf::io::Writer::new(counting)),
                count,
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
    /// Call [`BulkWriter::flush`] between batches of records to make
    /// [`BulkWriter::compressed_bytes`] reflect reality promptly.
    pub fn write(&mut self, header: &vcf::Header, record: &RecordBuf) -> Result<(), BulkError> {
        match &mut self.sink {
            Sink::Bcf(x) => x.write_variant_record(header, record)?,
            Sink::VcfGz(x) => x.write_variant_record(header, record)?,
            Sink::Vcf(x) => x.write_variant_record(header, record)?,
        }
        Ok(())
    }

    /// Flushes the underlying bgzf stream, dispatching any buffered data to
    /// the compression thread pool.
    ///
    /// The multithreaded bgzf writer only ships a block once its ~64 KiB
    /// uncompressed staging buffer fills, so [`BulkWriter::compressed_bytes`]
    /// can otherwise lag well behind what has actually been written. A
    /// size-targeting loop should call this between batches of records
    /// (the design's "record blocks") before checking
    /// [`BulkWriter::compressed_bytes`] — not after every single record,
    /// which would defeat both compression and parallelism by fragmenting
    /// the output into many tiny bgzf blocks.
    ///
    /// Dispatch to the compression/writer threads is asynchronous even
    /// after this call returns, so [`BulkWriter::compressed_bytes`] may
    /// still lag slightly; it becomes exact only after
    /// [`BulkWriter::finish_and_index`].
    pub fn flush(&mut self) -> Result<(), BulkError> {
        match &mut self.sink {
            Sink::Bcf(x) => x.get_mut().flush()?,
            Sink::VcfGz(x) => x.get_mut().flush()?,
            Sink::Vcf(x) => x.get_mut().flush()?,
        }
        Ok(())
    }

    /// Compressed bytes written so far (polled by the size-targeting loop).
    ///
    /// For uncompressed `Vcf` output this is simply the byte count.
    pub fn compressed_bytes(&self) -> u64 {
        self.count.load(Ordering::Relaxed)
    }

    /// Finishes the stream, then writes `<path>.csi` via a second read pass.
    ///
    /// Indexing only applies to `Bcf`; `VcfGz`/`Vcf` are not indexed here.
    pub fn finish_and_index(self, path: &Path) -> Result<(), BulkError> {
        let format = self.format;

        // Dropping the sink drops the underlying `MultithreadedWriter`, whose
        // `Drop` impl calls `finish()` (flushing remaining buffers, joining
        // the compression/writer threads, and appending the BGZF EOF block).
        // For the uncompressed `Vcf` sink, dropping the `CountingWriter`
        // simply drops the `File`, which flushes on close.
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

    #[test]
    fn counting_writer_counts_bytes_through() {
        let (mut w, count) = CountingWriter::new(Vec::new());
        use std::io::Write;
        w.write_all(b"hello").unwrap();
        w.flush().unwrap();
        assert_eq!(count.load(std::sync::atomic::Ordering::Relaxed), 5);
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
        // The multithreaded bgzf writer dispatches blocks to background
        // compression/writer threads; `flush()` starts that dispatch, but
        // completion is asynchronous, so poll briefly rather than assuming
        // it has already landed by the time `flush()` returns.
        w.flush().unwrap();
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(2);
        while w.compressed_bytes() == 0 && std::time::Instant::now() < deadline {
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        assert!(
            w.compressed_bytes() > 0,
            "counter should see compressed bytes"
        );
        w.finish_and_index(&path).unwrap();

        assert!(path.exists());
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
        fn run(workers: usize) -> Vec<u8> {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("t.bcf");
            let h = header();
            let mut w =
                BulkWriter::create(&path, Format::Bcf, &h, 6, NonZero::new(workers).unwrap())
                    .unwrap();
            for pos in 1usize..=500 {
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
            w.finish_and_index(&path).unwrap();
            std::fs::read(&path).unwrap()
        }

        let a = run(1);
        let b = run(8);
        assert_eq!(
            a, b,
            "output must be byte-identical regardless of thread count"
        );
        assert!(!a.is_empty());
    }
}
