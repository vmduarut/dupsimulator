use std::io::{self, BufRead, Write};

use rand::Rng;

#[derive(Debug, Default, Eq, PartialEq)]
pub struct SimulationStats {
    pub source_units: usize,
    pub duplicate_groups: usize,
    pub duplicate_units: usize,
}

impl SimulationStats {
    pub fn output_units(&self) -> usize {
        self.source_units + self.duplicate_units
    }
}

#[derive(Clone)]
struct FastqRecord {
    header: String,
    sequence: String,
    separator: String,
    quality: String,
}

#[cfg(test)]
pub fn simulate_records(
    input: impl BufRead,
    output: impl Write,
    source_selection_probability: f64,
    mean_extra_copies: f64,
    rng: &mut impl Rng,
) -> io::Result<SimulationStats> {
    simulate_records_with_report(
        input,
        output,
        source_selection_probability,
        mean_extra_copies,
        rng,
        None::<&mut Vec<u8>>,
    )
}

pub fn simulate_records_with_report<R: Write>(
    input: impl BufRead,
    mut output: impl Write,
    source_selection_probability: f64,
    mean_extra_copies: f64,
    rng: &mut impl Rng,
    mut report: Option<&mut R>,
) -> io::Result<SimulationStats> {
    let mut records = FastqReader::new(input);
    let mut stats = SimulationStats::default();
    while let Some(record) = records.next_record()? {
        stats.source_units += 1;

        if rng.random_bool(source_selection_probability) {
            let extra_copies = sample_geometric(mean_extra_copies, rng);
            stats.duplicate_groups += 1;
            stats.duplicate_units += extra_copies;
            write_family_member(&record, 0, &mut output, &mut report)?;
            for duplicate_number in 1..=extra_copies {
                write_family_member(&record, duplicate_number, &mut output, &mut report)?;
            }
        } else {
            write_record(&record, &mut output)?;
        }
    }

    output.flush()?;
    Ok(stats)
}

#[cfg(test)]
pub fn simulate_paired_records(
    input_r1: impl BufRead,
    input_r2: impl BufRead,
    output_r1: impl Write,
    output_r2: impl Write,
    source_selection_probability: f64,
    mean_extra_copies: f64,
    rng: &mut impl Rng,
) -> io::Result<SimulationStats> {
    simulate_paired_records_with_report(
        input_r1,
        input_r2,
        output_r1,
        output_r2,
        source_selection_probability,
        mean_extra_copies,
        rng,
        None::<&mut Vec<u8>>,
    )
}

pub fn simulate_paired_records_with_report<R: Write>(
    input_r1: impl BufRead,
    input_r2: impl BufRead,
    mut output_r1: impl Write,
    mut output_r2: impl Write,
    source_selection_probability: f64,
    mean_extra_copies: f64,
    rng: &mut impl Rng,
    mut report: Option<&mut R>,
) -> io::Result<SimulationStats> {
    let mut r1_records = FastqReader::new(input_r1);
    let mut r2_records = FastqReader::new(input_r2);
    let mut stats = SimulationStats::default();

    loop {
        let record_number = r1_records.record_number.max(r2_records.record_number) + 1;
        let r1 = r1_records.next_record()?;
        let r2 = r2_records.next_record()?;

        let (r1, r2) = match (r1, r2) {
            (None, None) => break,
            (Some(r1), Some(r2)) => (r1, r2),
            _ => {
                return Err(invalid_record(
                    record_number,
                    "paired FASTQ inputs have different record counts",
                ));
            }
        };

        if paired_identifier(&r1.header) != paired_identifier(&r2.header) {
            return Err(invalid_record(
                record_number,
                "paired FASTQ headers have different identifiers",
            ));
        }

        stats.source_units += 1;

        if rng.random_bool(source_selection_probability) {
            let extra_copies = sample_geometric(mean_extra_copies, rng);
            stats.duplicate_groups += 1;
            stats.duplicate_units += extra_copies;
            write_family_member(&r1, 0, &mut output_r1, &mut report)?;
            write_family_member(&r2, 0, &mut output_r2, &mut report)?;
            for duplicate_number in 1..=extra_copies {
                write_family_member(&r1, duplicate_number, &mut output_r1, &mut report)?;
                write_family_member(&r2, duplicate_number, &mut output_r2, &mut report)?;
            }
        } else {
            write_record(&r1, &mut output_r1)?;
            write_record(&r2, &mut output_r2)?;
        }
    }

    output_r1.flush()?;
    output_r2.flush()?;
    Ok(stats)
}

struct FastqReader<R> {
    input: R,
    record_number: usize,
}

impl<R: BufRead> FastqReader<R> {
    fn new(input: R) -> Self {
        Self {
            input,
            record_number: 0,
        }
    }

    fn next_record(&mut self) -> io::Result<Option<FastqRecord>> {
        let record_number = self.record_number + 1;
        let Some(header) = read_line(&mut self.input)? else {
            return Ok(None);
        };
        let sequence = required_line(&mut self.input, record_number)?;
        let separator = required_line(&mut self.input, record_number)?;
        let quality = required_line(&mut self.input, record_number)?;
        let record = FastqRecord {
            header,
            sequence,
            separator,
            quality,
        };

        validate_record(&record, record_number)?;
        self.record_number = record_number;
        Ok(Some(record))
    }
}

fn read_line(input: &mut impl BufRead) -> io::Result<Option<String>> {
    let mut line = String::new();
    if input.read_line(&mut line)? == 0 {
        return Ok(None);
    }

    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }

    Ok(Some(line))
}

fn required_line(input: &mut impl BufRead, record_number: usize) -> io::Result<String> {
    read_line(input)?.ok_or_else(|| invalid_record(record_number, "is incomplete"))
}

fn write_record(record: &FastqRecord, output: &mut impl Write) -> io::Result<()> {
    writeln!(output, "{}", record.header)?;
    writeln!(output, "{}", record.sequence)?;
    writeln!(output, "{}", record.separator)?;
    writeln!(output, "{}", record.quality)
}

fn duplicate_header(header: &str, duplicate_number: usize) -> String {
    let identifier_end = header.find(char::is_whitespace).unwrap_or(header.len());
    let identifier = &header[..identifier_end];
    let (identifier, mate_suffix) = match identifier
        .strip_suffix("/1")
        .or_else(|| identifier.strip_suffix("/2"))
    {
        Some(identifier) => (
            &header[..identifier.len()],
            &header[identifier.len()..identifier_end],
        ),
        None => (identifier, ""),
    };
    format!(
        "{}__dup={}{}{}",
        identifier,
        duplicate_number,
        mate_suffix,
        &header[identifier_end..]
    )
}

fn write_family_member<R: Write>(
    record: &FastqRecord,
    copy_number: usize,
    output: &mut impl Write,
    report: &mut Option<&mut R>,
) -> io::Result<()> {
    let mut member = record.clone();
    member.header = duplicate_header(&record.header, copy_number);
    write_record(&member, output)?;
    if let Some(report) = report.as_deref_mut() {
        writeln!(
            report,
            "{}\t{}\t{copy_number}\t{}",
            paired_identifier(&record.header),
            if copy_number == 0 {
                "original"
            } else {
                "duplicate"
            },
            read_identifier(&member.header),
        )?;
    }
    Ok(())
}

fn paired_identifier(header: &str) -> &str {
    let identifier = header[1..].split_whitespace().next().unwrap_or_default();
    identifier
        .strip_suffix("/1")
        .or_else(|| identifier.strip_suffix("/2"))
        .unwrap_or(identifier)
}

fn read_identifier(header: &str) -> &str {
    header[1..].split_whitespace().next().unwrap_or_default()
}

fn sample_geometric(mean: f64, rng: &mut impl Rng) -> usize {
    if mean == 1.0 {
        return 1;
    }

    let success_probability = 1.0 / mean;
    ((-rng.random::<f64>()).ln_1p() / (-success_probability).ln_1p()).floor() as usize + 1
}

fn validate_record(record: &FastqRecord, record_number: usize) -> io::Result<()> {
    if !record.header.starts_with('@') {
        return Err(invalid_record(record_number, "header must start with '@'"));
    }
    if !record.separator.starts_with('+') {
        return Err(invalid_record(
            record_number,
            "separator must start with '+'",
        ));
    }
    if record.sequence.len() != record.quality.len() {
        return Err(invalid_record(
            record_number,
            "sequence and quality lengths differ",
        ));
    }

    Ok(())
}

fn invalid_record(record_number: usize, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("FASTQ record {record_number} {message}"),
    )
}

#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor, Write};

    use rand::{SeedableRng, rngs::StdRng};

    use super::{
        SimulationStats, simulate_paired_records, simulate_records, simulate_records_with_report,
    };

    const FASTQ: &[u8] = b"@read-1\nACGT\n+\n!#$%\n@read-2\nTA\n+\n&'\n";
    const R1: &[u8] = b"@read-1/1\nAC\n+\n!!\n@read-2/1\nGT\n+\n##\n";
    const R2: &[u8] = b"@read-1/2\nTG\n+\n$$\n@read-2/2\nCA\n+\n%%\n";

    #[test]
    fn prints_plain_fastq_records() {
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(42);

        simulate_records(
            BufReader::new(Cursor::new(FASTQ)),
            &mut output,
            0.0,
            1.0,
            &mut rng,
        )
        .unwrap();

        assert_eq!(output, FASTQ);
    }

    #[test]
    fn rejects_incomplete_record() {
        let mut rng = StdRng::seed_from_u64(42);
        let error = simulate_records(
            BufReader::new(Cursor::new(b"@read\nACGT\n+\n")),
            Vec::new(),
            0.0,
            1.0,
            &mut rng,
        )
        .unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("record 1 is incomplete"));
    }

    #[test]
    fn rejects_mismatched_sequence_and_quality_lengths() {
        let mut rng = StdRng::seed_from_u64(42);
        let error = simulate_records(
            BufReader::new(Cursor::new(b"@read\nACGT\n+\n!!!\n")),
            Vec::new(),
            0.0,
            1.0,
            &mut rng,
        )
        .unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(
            error
                .to_string()
                .contains("sequence and quality lengths differ")
        );
    }

    #[test]
    fn adds_suffixed_duplicates_without_changing_read_data() {
        let input = b"@read-1 instrument metadata\nACGT\n+\n!#$%\n";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(42);

        let stats = simulate_records(
            BufReader::new(Cursor::new(input)),
            &mut output,
            1.0,
            1.0,
            &mut rng,
        )
        .unwrap();

        assert_eq!(
            stats,
            SimulationStats {
                source_units: 1,
                duplicate_groups: 1,
                duplicate_units: 1,
            }
        );
        assert_eq!(
            output,
            b"@read-1__dup=0 instrument metadata\nACGT\n+\n!#$%\n@read-1__dup=1 instrument metadata\nACGT\n+\n!#$%\n"
        );
    }

    #[test]
    fn reports_original_and_duplicate_family_members() {
        let input = b"@read-1 instrument metadata\nACGT\n+\n!#$%\n";
        let mut output = Vec::new();
        let mut report = Vec::new();
        let mut rng = StdRng::seed_from_u64(42);

        simulate_records_with_report(
            BufReader::new(Cursor::new(input)),
            &mut output,
            1.0,
            1.0,
            &mut rng,
            Some(&mut report),
        )
        .unwrap();

        assert_eq!(
            String::from_utf8(report).unwrap(),
            "read-1\toriginal\t0\tread-1__dup=0\nread-1\tduplicate\t1\tread-1__dup=1\n"
        );
    }

    #[test]
    fn produces_reproducible_duplicates_with_the_same_seed() {
        let mut first_output = Vec::new();
        let mut second_output = Vec::new();
        let mut first_rng = StdRng::seed_from_u64(42);
        let mut second_rng = StdRng::seed_from_u64(42);

        simulate_records(
            BufReader::new(Cursor::new(FASTQ)),
            &mut first_output,
            0.75,
            2.0,
            &mut first_rng,
        )
        .unwrap();
        simulate_records(
            BufReader::new(Cursor::new(FASTQ)),
            &mut second_output,
            0.75,
            2.0,
            &mut second_rng,
        )
        .unwrap();

        assert_eq!(first_output, second_output);
    }

    #[test]
    fn duplicates_paired_reads_together_with_matching_suffixes() {
        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        let mut rng = StdRng::seed_from_u64(42);

        let stats = simulate_paired_records(
            BufReader::new(Cursor::new(R1)),
            BufReader::new(Cursor::new(R2)),
            &mut r1_output,
            &mut r2_output,
            1.0,
            1.0,
            &mut rng,
        )
        .unwrap();

        assert_eq!(
            stats,
            SimulationStats {
                source_units: 2,
                duplicate_groups: 2,
                duplicate_units: 2,
            }
        );
        assert_eq!(
            r1_output,
            b"@read-1__dup=0/1\nAC\n+\n!!\n@read-1__dup=1/1\nAC\n+\n!!\n@read-2__dup=0/1\nGT\n+\n##\n@read-2__dup=1/1\nGT\n+\n##\n"
        );
        assert_eq!(
            r2_output,
            b"@read-1__dup=0/2\nTG\n+\n$$\n@read-1__dup=1/2\nTG\n+\n$$\n@read-2__dup=0/2\nCA\n+\n%%\n@read-2__dup=1/2\nCA\n+\n%%\n"
        );
    }

    #[test]
    fn accepts_illumina_style_paired_headers() {
        let r1 = b"@cluster-1 1:N:0:1\nAC\n+\n!!\n";
        let r2 = b"@cluster-1 2:N:0:1\nTG\n+\n$$\n";
        let mut rng = StdRng::seed_from_u64(42);

        simulate_paired_records(
            BufReader::new(Cursor::new(r1)),
            BufReader::new(Cursor::new(r2)),
            Vec::new(),
            Vec::new(),
            0.0,
            1.0,
            &mut rng,
        )
        .unwrap();
    }

    #[test]
    fn rejects_paired_files_with_different_record_counts() {
        let mut rng = StdRng::seed_from_u64(42);
        let error = simulate_paired_records(
            BufReader::new(Cursor::new(R1)),
            BufReader::new(Cursor::new(b"@read-1/2\nTG\n+\n$$\n")),
            Vec::new(),
            Vec::new(),
            0.0,
            1.0,
            &mut rng,
        )
        .unwrap_err();

        assert!(error.to_string().contains("different record counts"));
    }

    #[test]
    fn rejects_paired_files_with_different_identifiers() {
        let mut rng = StdRng::seed_from_u64(42);
        let error = simulate_paired_records(
            BufReader::new(Cursor::new(R1)),
            BufReader::new(Cursor::new(b"@other/2\nTG\n+\n$$\n")),
            Vec::new(),
            Vec::new(),
            0.0,
            1.0,
            &mut rng,
        )
        .unwrap_err();

        assert!(error.to_string().contains("different identifiers"));
    }

    #[test]
    fn reaches_the_requested_duplicate_fraction_in_expectation() {
        let original_records = 10_000;
        let mut input = Vec::new();
        for index in 0..original_records {
            writeln!(input, "@read-{index}\nAC\n+\n!!").unwrap();
        }

        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(7);
        simulate_records(
            BufReader::new(Cursor::new(input)),
            &mut output,
            0.125,
            2.0,
            &mut rng,
        )
        .unwrap();

        let output_records = output.iter().filter(|&&byte| byte == b'\n').count() / 4;
        let duplicate_fraction = (output_records - original_records) as f64 / output_records as f64;

        assert!((duplicate_fraction - 0.2).abs() < 0.02);
    }
}
