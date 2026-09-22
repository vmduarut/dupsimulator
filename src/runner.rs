use std::{io, path::Path};

use rand::{SeedableRng, rngs::StdRng};

use crate::{
    args::Args,
    fastq, input,
    output::{self, Report},
};

pub(crate) fn random_generator(seed: Option<u64>) -> StdRng {
    match seed {
        Some(seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_os_rng(),
    }
}

pub(crate) fn run(
    args: &Args,
    source_selection_probability: f64,
    rng: &mut StdRng,
    report: Option<&mut Report>,
) -> io::Result<(fastq::SimulationStats, &'static str)> {
    // The R2 input determines whether source units are records or read pairs.
    match args.fastq_r2.as_deref() {
        None => Ok((
            run_single_end(args, source_selection_probability, rng, report)?,
            "records",
        )),
        Some(fastq_r2) => Ok((
            run_paired_end(args, fastq_r2, source_selection_probability, rng, report)?,
            "read pairs",
        )),
    }
}

pub(crate) fn log_completion(stats: &fastq::SimulationStats, unit_label: &str) {
    eprintln!(
        "Simulation complete: input {unit_label}={}, duplicate groups={}, duplicate {unit_label}={}, output {unit_label}={}",
        stats.source_units,
        stats.duplicate_groups,
        stats.duplicate_units,
        stats.output_units(),
    );
}

fn run_single_end(
    args: &Args,
    source_selection_probability: f64,
    rng: &mut StdRng,
    report: Option<&mut Report>,
) -> io::Result<fastq::SimulationStats> {
    // An R2 destination has no meaning unless paired input was supplied.
    if args.output_r2.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--output-r2 requires an R2 input",
        ));
    }
    let output_r1 = output::resolve_output_path(args.output_r1.as_deref(), &args.fastq_r1, "R1");
    log_start(args, "single-end", &output_r1, None);

    // The simulation samples extra copies; the CLI value includes the original record.
    fastq::simulate_records_with_report(
        input::open_fastq(&args.fastq_r1)?,
        output::output_writer(&output_r1)?,
        source_selection_probability,
        args.duplicate_group_size - 1.0,
        rng,
        report,
    )
}

fn run_paired_end(
    args: &Args,
    fastq_r2: &Path,
    source_selection_probability: f64,
    rng: &mut StdRng,
    report: Option<&mut Report>,
) -> io::Result<fastq::SimulationStats> {
    let output_r1 = output::resolve_output_path(args.output_r1.as_deref(), &args.fastq_r1, "R1");
    let output_r2 = output::resolve_output_path(args.output_r2.as_deref(), &args.fastq_r1, "R2");
    log_start(args, "paired-end", &output_r1, Some((fastq_r2, &output_r2)));

    // A paired family is emitted together so both mates receive matching duplicate numbers.
    fastq::simulate_paired_records_with_report(
        input::open_fastq(&args.fastq_r1)?,
        input::open_fastq(fastq_r2)?,
        output::output_writer(&output_r1)?,
        output::output_writer(&output_r2)?,
        source_selection_probability,
        args.duplicate_group_size - 1.0,
        rng,
        report,
    )
}

fn log_start(args: &Args, mode: &str, output_r1: &Path, r2_paths: Option<(&Path, &Path)>) {
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
