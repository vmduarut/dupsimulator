use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use flate2::read::GzDecoder;

pub fn open_fastq(path: &Path) -> io::Result<Box<dyn BufRead>> {
    reader_for_input(BufReader::new(File::open(path)?))
}

fn reader_for_input(mut input: impl BufRead + 'static) -> io::Result<Box<dyn BufRead>> {
    // Compression is determined from gzip's magic bytes, not the filename extension.
    if input.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::new(GzDecoder::new(input))))
    } else {
        Ok(Box::new(input))
    }
}

#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor, Write};

    use flate2::{Compression, write::GzEncoder};
    use rand::{SeedableRng, rngs::StdRng};

    use super::reader_for_input;
    use crate::fastq::simulate_records;

    const FASTQ: &[u8] = b"@read-1\nACGT\n+\n!#$%\n";

    #[test]
    fn detects_and_decompresses_gzip_fastq() {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(FASTQ).unwrap();
        let compressed = encoder.finish().unwrap();
        let input = reader_for_input(BufReader::new(Cursor::new(compressed))).unwrap();
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(42);

        simulate_records(input, &mut output, 0.0, 1.0, &mut rng).unwrap();

        assert_eq!(output, FASTQ);
    }
}
