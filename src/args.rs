use std::path::PathBuf;

use clap::Parser;

#[derive(Parser)]
#[command(about = "Simulate duplicate FASTQ reads from plain or gzip-compressed input")]
pub(crate) struct Args {
    /// Path to the R1 FASTQ file, optionally gzip-compressed
    pub(crate) fastq_r1: PathBuf,

    /// Path to the R2 FASTQ file for paired-end simulation
    pub(crate) fastq_r2: Option<PathBuf>,

    /// Destination R1 FASTQ path; defaults to {sampleID}_dup_R1.fastq.gz
    #[arg(long)]
    pub(crate) output_r1: Option<PathBuf>,

    /// Destination R2 FASTQ path; defaults to {sampleID}_dup_R2.fastq.gz
    #[arg(long)]
    pub(crate) output_r2: Option<PathBuf>,

    /// TSV destination describing each simulated duplicate family
    #[arg(long)]
    pub(crate) report: Option<PathBuf>,

    /// Target expected fraction of output records that are duplicates
    #[arg(long, default_value_t = 0.0)]
    pub(crate) duplicate_rate: f64,

    /// Mean total records per duplicate group, including the original read
    #[arg(long, default_value_t = 2.0)]
    pub(crate) duplicate_group_size: f64,

    /// Random seed for reproducible duplicate simulation
    #[arg(long)]
    pub(crate) seed: Option<u64>,
}

pub(crate) fn source_selection_probability(
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

    // Solve for the fraction of source reads to select so the generated output reaches the target
    // duplicate fraction in expectation.
    let probability = duplicate_rate / ((1.0 - duplicate_rate) * (duplicate_group_size - 1.0));
    if probability > 1.0 {
        return Err("the requested duplicate rate requires a larger --duplicate-group-size".into());
    }

    Ok(probability)
}

#[cfg(test)]
mod tests {
    use super::source_selection_probability;

    #[test]
    fn calibrates_the_source_selection_probability() {
        assert!((source_selection_probability(0.2, 3.0).unwrap() - 0.125).abs() < f64::EPSILON);
    }

    #[test]
    fn rejects_unachievable_rate_and_group_size() {
        assert!(source_selection_probability(0.8, 2.0).is_err());
    }
}
