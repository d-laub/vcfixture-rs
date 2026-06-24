use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::index::BinnedIndex};

use crate::error::BuildError;
use crate::model::{Document, Record};
use crate::value::{FieldValue, Scalar};

/// Options controlling how a VCF document is written to disk.
#[derive(Debug, Clone, Copy, Default)]
pub struct WriteOpts {
    pub bgzip: bool,
    pub index: bool,
}

impl WriteOpts {
    pub fn text() -> WriteOpts {
        WriteOpts {
            bgzip: false,
            index: false,
        }
    }
    pub fn bgzipped() -> WriteOpts {
        WriteOpts {
            bgzip: true,
            index: false,
        }
    }
    pub fn bgzipped_indexed() -> WriteOpts {
        WriteOpts {
            bgzip: true,
            index: true,
        }
    }
}

// ---------------------------------------------------------------------------
// Percent-encoding for reserved characters in string values. '%' must be first.
// ---------------------------------------------------------------------------
const PERCENT: &[(char, &str)] = &[
    ('%', "%25"),
    (':', "%3A"),
    (';', "%3B"),
    ('=', "%3D"),
    (',', "%2C"),
    ('\r', "%0D"),
    ('\n', "%0A"),
    ('\t', "%09"),
];

fn encode(s: &str) -> String {
    let mut out = s.to_string();
    for (ch, rep) in PERCENT {
        out = out.replace(*ch, rep);
    }
    out
}

fn fmt_scalar(s: &Scalar) -> String {
    match s {
        Scalar::Int(n) => n.to_string(),
        Scalar::Float(f) => {
            if f.is_nan() || f.is_infinite() {
                ".".to_string()
            } else {
                format!("{f}")
            }
        }
        Scalar::Char(c) => encode(&c.to_string()),
        Scalar::Str(s) => encode(s),
    }
}

fn fmt_opt_scalar(s: &Option<Scalar>) -> String {
    match s {
        Some(v) => fmt_scalar(v),
        None => ".".to_string(),
    }
}

fn fmt_value(v: &FieldValue) -> String {
    match v {
        FieldValue::Flag => "1".to_string(), // never reached for INFO rendering
        FieldValue::Scalar(s) => fmt_scalar(s),
        FieldValue::List(xs) => {
            if xs.is_empty() {
                ".".to_string()
            } else {
                xs.iter().map(fmt_opt_scalar).collect::<Vec<_>>().join(",")
            }
        }
    }
}

fn render_info(rec: &Record) -> String {
    if rec.info.is_empty() {
        return ".".to_string();
    }
    let mut parts = Vec::new();
    for (key, val) in &rec.info {
        match val {
            FieldValue::Flag => parts.push(key.clone()),
            other => parts.push(format!("{key}={}", fmt_value(other))),
        }
    }
    if parts.is_empty() {
        ".".to_string()
    } else {
        parts.join(";")
    }
}

fn render_sample(rec: &Record, si: usize) -> String {
    let sample = &rec.samples[si];
    let mut vals = Vec::with_capacity(rec.fmt_keys.len());
    for key in &rec.fmt_keys {
        if key == "GT" {
            vals.push(
                sample
                    .gt
                    .as_ref()
                    .map(|g| g.render())
                    .unwrap_or_else(|| ".".to_string()),
            );
        } else {
            match sample.values.get(key) {
                Some(v) => vals.push(fmt_value(v)),
                None => vals.push(".".to_string()),
            }
        }
    }
    vals.join(":")
}

fn render_record(rec: &Record) -> String {
    let ids = match &rec.ids {
        Some(v) if !v.is_empty() => v.join(";"),
        _ => ".".to_string(),
    };
    let alt = if rec.alts.is_empty() {
        ".".to_string()
    } else {
        rec.alts
            .iter()
            .map(|a| a.render())
            .collect::<Vec<_>>()
            .join(",")
    };
    let qual = match rec.qual {
        Some(q) => fmt_scalar(&Scalar::Float(q)),
        None => ".".to_string(),
    };
    let filt = match &rec.filters {
        None => ".".to_string(),
        Some(v) if v.is_empty() => "PASS".to_string(),
        Some(v) => v.join(";"),
    };
    let mut cols = vec![
        rec.chrom.clone(),
        rec.pos.to_string(),
        ids,
        rec.ref_.clone(),
        alt,
        qual,
        filt,
        render_info(rec),
    ];
    if !rec.fmt_keys.is_empty() {
        cols.push(rec.fmt_keys.join(":"));
        for si in 0..rec.samples.len() {
            cols.push(render_sample(rec, si));
        }
    }
    cols.join("\t")
}

/// Render a `Document` to a VCF text string.
pub fn render(doc: &Document) -> String {
    let mut lines = vec![format!("##fileformat={}", doc.version.as_str())];
    for c in &doc.contigs {
        lines.push(c.header_line());
    }
    for f in &doc.info_defs {
        lines.push(f.header_line());
    }
    for (id, desc) in &doc.filter_defs {
        lines.push(format!("##FILTER=<ID={id},Description=\"{desc}\">"));
    }
    for f in &doc.format_defs {
        lines.push(f.header_line());
    }
    for ad in &doc.alt_defs {
        lines.push(ad.header_line());
    }
    let mut header = vec![
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let has_fmt = !doc.format_defs.is_empty() || doc.records.iter().any(|r| !r.fmt_keys.is_empty());
    if has_fmt {
        header.push("FORMAT".to_string());
        header.extend(doc.samples.iter().cloned());
    }
    lines.push(header.join("\t"));
    for rec in &doc.records {
        lines.push(render_record(rec));
    }
    let mut text = lines.join("\n");
    text.push('\n');
    text
}

/// Write a `Document` to `path`, applying `opts` (plain text, bgzipped, or bgzipped+indexed).
///
/// Returns the actual path written (may gain a `.gz` suffix when bgzip is requested).
pub fn write(
    doc: &Document,
    path: impl AsRef<Path>,
    opts: WriteOpts,
) -> Result<PathBuf, BuildError> {
    let text = render(doc);
    let mut path = path.as_ref().to_path_buf();
    if !opts.bgzip {
        fs::write(&path, text)?;
        return Ok(path);
    }
    // Ensure path ends with .gz
    if path.extension().and_then(|e| e.to_str()) != Some("gz") {
        let mut name = path.into_os_string();
        name.push(".gz");
        path = PathBuf::from(name);
    }
    // bgzip-compress the text
    let file = fs::File::create(&path)?;
    let mut writer = bgzf::io::Writer::new(file);
    writer.write_all(text.as_bytes())?;
    writer.finish()?;
    if opts.index {
        write_csi(doc, &path)?;
    }
    Ok(path)
}

/// Build and write a CSI index alongside a bgzipped VCF.
///
/// Uses `noodles_csi` (v0.49) with a `BinnedIndex` (the CSI-native index type).
/// Writes `<path>.csi` in the standard noodles CSI binary format.
///
/// API used:
///   - `csi::binning_index::Indexer::<BinnedIndex>::default()`
///   - `indexer.set_header(csi::binning_index::index::header::Builder::vcf().build())`
///   - `indexer.add_record(Some((ref_id, start, end, true)), chunk)` — one call per record
///   - `indexer.build(n_ref_seqs)` → `csi::Index`
///   - `csi::fs::write(csi_path, &index)` — writes bgzf-compressed CSI file
fn write_csi(doc: &Document, bgzf_path: &Path) -> Result<(), BuildError> {
    // Build a contig-name → id mapping from the document's contigs.
    let contig_ids: std::collections::HashMap<&str, usize> = doc
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.id.as_str(), i))
        .collect();

    let header = csi::binning_index::index::header::Builder::vcf().build();

    let mut indexer = csi::binning_index::Indexer::<BinnedIndex>::default().set_header(header);

    // Re-open the bgzf file to track virtual positions per record.
    // We reproduce the write in order to know the start/end virtual positions.
    // Strategy: re-render the text, then replay through a bgzf writer into a memory buffer
    // while recording positions.
    let text = render(doc);
    let lines: Vec<&str> = text.lines().collect();

    // Count header lines (those starting with '#') to find the first data line.
    let first_data_line = lines.iter().take_while(|l| l.starts_with('#')).count();

    // Replay through a bgzf writer into Vec<u8> to track virtual positions.
    let mut mem_writer = bgzf::io::Writer::new(Vec::<u8>::new());

    // Write header lines.
    for line in &lines[..first_data_line] {
        mem_writer.write_all(line.as_bytes())?;
        mem_writer.write_all(b"\n")?;
    }

    // Write data lines and record virtual positions.
    for (rec, line) in doc.records.iter().zip(lines[first_data_line..].iter()) {
        let start_vpos = mem_writer.virtual_position();
        mem_writer.write_all(line.as_bytes())?;
        mem_writer.write_all(b"\n")?;
        let end_vpos = mem_writer.virtual_position();

        let chunk =
            csi::binning_index::index::reference_sequence::bin::Chunk::new(start_vpos, end_vpos);

        // VCF POS is 1-based; noodles Position is also 1-based.
        let start = Position::try_from(rec.pos as usize)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
        // For SNPs/short variants: end = pos + len(REF) - 1
        let end_pos = rec.pos as usize + rec.ref_.len().saturating_sub(1);
        let end = Position::try_from(end_pos)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

        let ref_id = contig_ids.get(rec.chrom.as_str()).copied().unwrap_or(0); // fallback to 0 if contig not declared

        indexer.add_record(Some((ref_id, start, end, true)), chunk)?;
    }

    let n_ref = doc.contigs.len();
    let index = indexer.build(n_ref);

    // Write <path>.csi
    let csi_path = {
        let mut p = bgzf_path.as_os_str().to_os_string();
        p.push(".csi");
        PathBuf::from(p)
    };
    csi::fs::write(&csi_path, &index)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::allele::Allele;
    use crate::build::{RecordSpec, VcfBuilder};
    use crate::spec::version::LATEST;
    use crate::value::FieldValue;

    fn base() -> VcfBuilder {
        VcfBuilder::new(["s1"], [("chr1", Some(100_000u64))], LATEST)
    }

    /// Returns the single data line (the last non-empty line) of a rendered doc.
    fn data_line(text: &str) -> &str {
        text.lines()
            .rfind(|l| !l.starts_with('#') && !l.is_empty())
            .expect("expected at least one data line")
    }

    /// Returns the 0-based tab-separated column of the single data line.
    fn data_col(text: &str, idx: usize) -> &str {
        data_line(text).split('\t').nth(idx).expect("column exists")
    }

    // --- FILTER three-way -------------------------------------------------

    #[test]
    fn filter_empty_renders_pass() {
        let text = base()
            .format("GT")
            .record(
                RecordSpec::at("chr1", 100)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1"])
                    .filter(Vec::<&str>::new()),
            )
            .render()
            .unwrap();
        // FILTER is column index 6 (CHROM POS ID REF ALT QUAL FILTER ...)
        assert_eq!(data_col(&text, 6), "PASS");
    }

    #[test]
    fn filter_unset_renders_dot() {
        let text = base()
            .format("GT")
            .record(
                RecordSpec::at("chr1", 100)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1"]),
            )
            .render()
            .unwrap();
        assert_eq!(data_col(&text, 6), ".");
    }

    #[test]
    fn filter_named_renders_value() {
        let text = base()
            .filter("q10", "Quality below 10")
            .format("GT")
            .record(
                RecordSpec::at("chr1", 100)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1"])
                    .filter(["q10"]),
            )
            .render()
            .unwrap();
        assert_eq!(data_col(&text, 6), "q10");
        // And the FILTER header definition is emitted.
        assert!(text.contains("##FILTER=<ID=q10,Description=\"Quality below 10\">"));
    }

    // --- Flag INFO bare key ----------------------------------------------

    #[test]
    fn flag_info_renders_bare_key() {
        let text = base()
            .info("DB") // reserved DB = Flag
            .format("GT")
            .record(
                RecordSpec::at("chr1", 100)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1"])
                    .info("DB", FieldValue::Flag),
            )
            .render()
            .unwrap();
        // INFO is column index 7.
        assert_eq!(data_col(&text, 7), "DB");
        assert!(!text.contains("DB="));
    }

    // --- Percent-encoding -------------------------------------------------

    #[test]
    fn string_info_percent_encodes_semicolon() {
        use crate::build::Field;
        use crate::spec::number::Number;
        use crate::spec::types::Type;
        // A String INFO value containing ';' must render as '%3B'.
        let text = base()
            .info(Field::typed("NOTE", Number::ONE, Type::String))
            .format("GT")
            .record(
                RecordSpec::at("chr1", 100)
                    .ref_("A")
                    .alt([Allele::seq("T").unwrap()])
                    .gt(["0|1"])
                    .info("NOTE", FieldValue::strings(["a;b"])),
            )
            .render()
            .unwrap();
        let info = data_col(&text, 7);
        assert!(info.contains("NOTE=a%3Bb"), "INFO was: {info}");
        // The raw reserved char must not survive inside the value.
        assert_eq!(info, "NOTE=a%3Bb");
    }

    // --- Symbolic ALT rendering ------------------------------------------

    #[test]
    fn symbolic_deletion_renders_angle_del() {
        // DEL at LATEST (>= 4.4) needs single-base REF + SVLEN + SVCLAIM.
        let text = base()
            .info("SVLEN")
            .info("SVCLAIM")
            .format("GT")
            .record(
                RecordSpec::at("chr1", 100)
                    .ref_("A")
                    .alt([Allele::deletion(Vec::<&str>::new())])
                    .gt(["0|1"])
                    .info("SVLEN", FieldValue::ints([100]))
                    .info("SVCLAIM", FieldValue::strings(["D"])),
            )
            .render()
            .unwrap();
        // ALT is column index 4.
        assert_eq!(data_col(&text, 4), "<DEL>");
    }
}
