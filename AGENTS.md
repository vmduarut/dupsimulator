# AGENTS.md

Single-binary Rust crate (edition 2024) that simulates PCR/optical duplicate FASTQ reads. CLI args live in `src/args.rs`; orchestration in `src/runner.rs`; output/report handling in `src/output.rs`; simulation core in `src/fastq.rs`; transparent gzip open in `src/input.rs`.

## Commands

- Unit tests: `cargo test` (plain `#[cfg(test)]` modules per file, no framework)
- Lint: `cargo clippy --all-targets`. `fastq.rs:96 simulate_paired_records_with_report` triggers a pre-existing `too_many_arguments` warning — leave the signature alone.
- No formatter/typecheck config beyond rustfmt defaults.

## Architecture & conventions

- `args.rs`: `clap` derive args; validation via `source_selection_probability`, a pure fn returning `Result<f64, String>` (never `io::Error`). `main` converts it with `.map_err(|m| io::Error::new(io::ErrorKind::InvalidInput, m))` — keep error concerns at the boundary. `main` returns `Result<(), Box<dyn std::error::Error>>` (no `anyhow`).
- `runner.rs`: owns RNG construction, single/paired dispatch, and lifecycle logs. `main.rs` only wires args, report setup, validation, execution, and report flush.
- `output.rs`: owns report creation/flushing, output writers, and default output-path naming.
- `fastq.rs`: single-end and paired simulation. The `#[cfg(test)]`-only fns `simulate_records`/`simulate_paired_records` are the public single/paired entrypoints used by tests; binaries call the `_with_report` variants.
- `rand` 0.9: build a `StdRng` via `SeedableRng`; `rng.random_bool(p)` picks records, `sample_geometric` draws copy counts. RNG is threaded as `&mut impl Rng`. Tests seed `StdRng::seed_from_u64(42)` for determinism; keep that reproducible pattern.
- `mean_extra_copies` passed to simulation is `duplicate_group_size - 1.0`.
- All output args (`--output-r1`, `--output-r2`) are optional. Omitted ones default to `{sampleID}_dup_R1.fastq.gz` / `_dup_R2.fastq.gz`, where `sampleID` is derived from the R1 filename in `sample_id` (strip `.fastq`/`.fq`/`.gz`, then one read-suffix like `_1`/`_R1`/`.r2`). The program always writes files; there is no stdout mode. `--output-r2` without an R2 input is an error.

## Gotchas

- Output compression is decided by substring match: `path.to_string_lossy().contains("gz")`, not file extension (`output.rs` `output_writer`). Default names end in `.gz`.
- Input compression is sniffed by magic bytes (`0x1f 0x8b`) via `fill_buf`, so the input filename extension is irrelevant.
- ~4 GB of `*.fastq.gz` test fixtures sit in the repo root and are gitignored. For manual runs, write a tiny FASTQ to `/tmp` — never create large files in the repo.
