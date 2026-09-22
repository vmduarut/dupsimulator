use std::io;

use clap::Parser;

mod args;
mod fastq;
mod input;
mod output;
mod runner;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = args::Args::parse();
    let mut report = output::open_report(args.report.as_deref())?;
    // Validation stays independent of I/O; convert its user-facing message at the CLI boundary.
    let source_selection_probability =
        args::source_selection_probability(args.duplicate_rate, args.duplicate_group_size)
            .map_err(|message| io::Error::new(io::ErrorKind::InvalidInput, message))?;
    let mut rng = runner::random_generator(args.seed);
    let (stats, unit_label) = runner::run(
        &args,
        source_selection_probability,
        &mut rng,
        report.as_mut(),
    )?;
    output::flush_report(&mut report)?;
    runner::log_completion(&stats, unit_label);

    Ok(())
}
