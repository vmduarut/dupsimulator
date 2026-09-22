use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

use flate2::{Compression, write::GzEncoder};

pub(crate) type Report = Box<dyn Write>;

pub(crate) fn open_report(path: Option<&Path>) -> io::Result<Option<Report>> {
    let mut report = path.map(output_writer).transpose()?;
    if let Some(report) = report.as_deref_mut() {
        // Write the header once before simulation appends one row per family member.
        writeln!(report, "family_id\trole\tcopy_number\tread_id")?;
    }
    Ok(report)
}

pub(crate) fn flush_report(report: &mut Option<Report>) -> io::Result<()> {
    if let Some(report) = report.as_deref_mut() {
        report.flush()?;
    }
    Ok(())
}

pub(crate) fn output_writer(path: &Path) -> io::Result<Box<dyn Write>> {
    let output = BufWriter::new(File::create(path)?);
    // Preserve the established behavior: any path containing "gz" is gzip-compressed.
    if path.to_string_lossy().contains("gz") {
        Ok(Box::new(GzEncoder::new(output, Compression::default())))
    } else {
        Ok(Box::new(output))
    }
}

pub(crate) fn resolve_output_path(output: Option<&Path>, input_r1: &Path, mate: &str) -> PathBuf {
    // Explicit paths override the predictable, sample-derived default names.
    output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| default_output_path(input_r1, mate))
}

fn sample_id(input_r1: &Path) -> String {
    const EXTENSIONS: [&str; 3] = [".fastq", ".fq", ".gz"];
    const READ_SUFFIXES: [&str; 14] = [
        "_1", "_2", ".1", ".2", "_R1", "_R2", "_r1", "_r2", ".R1", ".R2", ".r1", ".r2", "/1", "/2",
    ];
    let mut name = input_r1
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    // Repeatedly remove extensions so names such as "sample_R1.fastq.gz" work.
    while let Some(stripped) = EXTENSIONS
        .iter()
        .find_map(|extension| name.strip_suffix(extension))
    {
        name = stripped.to_string();
    }
    // Remove one mate suffix, but retain meaningful sample-name suffixes such as "_10".
    if let Some(stripped) = READ_SUFFIXES
        .iter()
        .find_map(|suffix| name.strip_suffix(suffix))
    {
        name = stripped.to_string();
    }
    name
}

fn default_output_path(input_r1: &Path, mate: &str) -> PathBuf {
    PathBuf::from(format!("{}_dup_{mate}.fastq.gz", sample_id(input_r1)))
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};
    use std::path::Path;

    use flate2::read::GzDecoder;

    use super::{default_output_path, output_writer, sample_id};

    #[test]
    fn derives_sample_id_from_read_one_filename() {
        assert_eq!(sample_id(Path::new("ERR2538964_1.fastq.gz")), "ERR2538964");
        assert_eq!(sample_id(Path::new("reads_R1.fq")), "reads");
        assert_eq!(sample_id(Path::new("reads.r1.gz")), "reads");
        assert_eq!(sample_id(Path::new("foo.fastq")), "foo");
        assert_eq!(sample_id(Path::new("sample_10.fastq")), "sample_10");
    }

    #[test]
    fn defaults_outputs_to_sample_id_dup_mate_fastq_gz() {
        let input = Path::new("ERR2538964_1.fastq.gz");
        assert_eq!(
            default_output_path(input, "R1"),
            Path::new("ERR2538964_dup_R1.fastq.gz")
        );
        assert_eq!(
            default_output_path(input, "R2"),
            Path::new("ERR2538964_dup_R2.fastq.gz")
        );
    }

    #[test]
    fn compresses_outputs_with_gz_in_their_name() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("output-gz.fastq");

        {
            let mut output = output_writer(&path).unwrap();
            output.write_all(b"FASTQ").unwrap();
        }

        let mut contents = String::new();
        GzDecoder::new(std::fs::File::open(path).unwrap())
            .read_to_string(&mut contents)
            .unwrap();
        assert_eq!(contents, "FASTQ");
    }
}
