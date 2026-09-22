use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

use clap::Parser;
use flate2::{Compression, write::GzEncoder};
use rand::{SeedableRng, rngs::StdRng};

mod fastq;
mod input;

#[derive(Parser)]
#[command(about = "Simulate duplicate FASTQ reads from plain or gzip-compressed input")]
struct Args {
    /// Path to the R1 FASTQ file, optionally gzip-compressed
    fastq_r1: PathBuf,

    /// Path to the R2 FASTQ file for paired-end simulation
    fastq_r2: Option<PathBuf>,

    /// Destination R1 FASTQ path; defaults to {sampleID}_dup_R1.fastq.gz
    #[arg(long)]
    output_r1: Option<PathBuf>,

    /// Destination R2 FASTQ path; defaults to {sampleID}_dup_R2.fastq.gz
    #[arg(long)]
    output_r2: Option<PathBuf>,

    /// TSV destination describing each simulated duplicate family
    #[arg(long)]
    report: Option<PathBuf>,

    /// Target expected fraction of output records that are duplicates
    #[arg(long, default_value_t = 0.0)]
    duplicate_rate: f64,

    /// Mean total records per duplicate group, including the original read
    #[arg(long, default_value_t = 2.0)]
    duplicate_group_size: f64,

    /// Random seed for reproducible duplicate simulation
    #[arg(long)]
    seed: Option<u64>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut report = args.report.as_deref().map(output_writer).transpose()?;
    if let Some(report) = report.as_deref_mut() {
        writeln!(report, "family_id\trole\tcopy_number\tread_id")?;
    }
    let source_selection_probability =
        source_selection_probability(args.duplicate_rate, args.duplicate_group_size)
            .map_err(|message| io::Error::new(io::ErrorKind::InvalidInput, message))?;
    let mut rng = match args.seed {
        Some(seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_os_rng(),
    };

    let (stats, unit_label) = match &args.fastq_r2 {
        None => {
            if args.output_r2.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--output-r2 requires an R2 input",
                )
                .into());
            }
            let output_r1 = args
                .output_r1
                .clone()
                .unwrap_or_else(|| default_output_path(&args.fastq_r1, "R1"));
            log_start(&args, "single-end", &output_r1, None);
            (
                fastq::simulate_records_with_report(
                    input::open_fastq(&args.fastq_r1)?,
                    output_writer(&output_r1)?,
                    source_selection_probability,
                    args.duplicate_group_size - 1.0,
                    &mut rng,
                    report.as_mut(),
                )?,
                "records",
            )
        }
        Some(fastq_r2) => {
            let output_r1 = args
                .output_r1
                .clone()
                .unwrap_or_else(|| default_output_path(&args.fastq_r1, "R1"));
            let output_r2 = args
                .output_r2
                .clone()
                .unwrap_or_else(|| default_output_path(&args.fastq_r1, "R2"));
            log_start(&args, "paired-end", &output_r1, Some((fastq_r2, &output_r2)));
            (
                fastq::simulate_paired_records_with_report(
                    input::open_fastq(&args.fastq_r1)?,
                    input::open_fastq(fastq_r2)?,
                    output_writer(&output_r1)?,
                    output_writer(&output_r2)?,
                    source_selection_probability,
                    args.duplicate_group_size - 1.0,
                    &mut rng,
                    report.as_mut(),
                )?,
                "read pairs",
            )
        }
    };
    if let Some(report) = report.as_deref_mut() {
        report.flush()?;
    }

    eprintln!(
        "Simulation complete: input {unit_label}={}, duplicate groups={}, duplicate {unit_label}={}, output {unit_label}={}",
        stats.source_units,
        stats.duplicate_groups,
        stats.duplicate_units,
        stats.output_units(),
    );

    Ok(())
}

fn output_writer(path: &Path) -> io::Result<Box<dyn Write>> {
    let output = BufWriter::new(File::create(path)?);
    if path.to_string_lossy().contains("gz") {
        Ok(Box::new(GzEncoder::new(output, Compression::default())))
    } else {
        Ok(Box::new(output))
    }
}

fn sample_id(input_r1: &Path) -> String {
    const EXTENSIONS: [&str; 3] = [".fastq", ".fq", ".gz"];
    const READ_SUFFIXES: [&str; 14] = [
        "_1", "_2", ".1", ".2", "_R1", "_R2", "_r1", "_r2", ".R1", ".R2", ".r1", ".r2", "/1",
        "/2",
    ];
    let mut name = input_r1
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    while let Some(stripped) = EXTENSIONS
        .iter()
        .find_map(|extension| name.strip_suffix(extension))
    {
        name = stripped.to_string();
    }
    if let Some(stripped) = READ_SUFFIXES.iter().find_map(|suffix| name.strip_suffix(suffix)) {
        name = stripped.to_string();
    }
    name
}

fn default_output_path(input_r1: &Path, mate: &str) -> PathBuf {
    PathBuf::from(format!("{}_dup_{mate}.fastq.gz", sample_id(input_r1)))
}

fn log_start(
    args: &Args,
    mode: &str,
    output_r1: &Path,
    r2_paths: Option<(&Path, &Path)>,
) {
    let seed = args
        .seed
        .map(|seed| seed.to_string())
        .unwrap_or_else(|| "random".into());
    eprintln!(
        "Starting {mode} simulation: input R1={}, output R1={}, duplicate rate={}, duplicate group size={}, seed={seed}",
        args.fastq_r1.display(),
        output_r1.display(),
        args.duplicate_rate,
        args.duplicate_group_size,
    );
    if let Some((input_r2, output_r2)) = r2_paths {
        eprintln!(
            "Paired-end paths: input R2={}, output R2={}",
            input_r2.display(),
            output_r2.display(),
        );
    }
}

fn source_selection_probability(
    duplicate_rate: f64,
    duplicate_group_size: f64,
) -> Result<f64, String> {
    if !duplicate_rate.is_finite() || !(0.0..1.0).contains(&duplicate_rate) {
        return Err(
            "--duplicate-rate must be a finite value from 0 (inclusive) to 1 (exclusive)".into(),
        );
    }
    if !duplicate_group_size.is_finite() || duplicate_group_size < 2.0 {
        return Err("--duplicate-group-size must be a finite value of at least 2".into());
    }

    let probability = duplicate_rate / ((1.0 - duplicate_rate) * (duplicate_group_size - 1.0));
    if probability > 1.0 {
        return Err("the requested duplicate rate requires a larger --duplicate-group-size".into());
    }

    Ok(probability)
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};
    use std::path::Path;

    use flate2::read::GzDecoder;

    use super::{default_output_path, output_writer, sample_id, source_selection_probability};

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
    fn calibrates_the_source_selection_probability() {
        assert!((source_selection_probability(0.2, 3.0).unwrap() - 0.125).abs() < f64::EPSILON);
    }

    #[test]
    fn rejects_unachievable_rate_and_group_size() {
        assert!(source_selection_probability(0.8, 2.0).is_err());
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
