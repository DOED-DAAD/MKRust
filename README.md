# MKRust

- [Introduction](#introduction)
  - [Citation](#citation)
  - [Contact](#contact)
- [Install](#install)
  - [Compatibility](#compatibility)
- [Getting Started](#getting-started)
  - [Usage](#usage)
  - [Configuration and Settings](#configuration-and-settings)
  - [Data Input](#data-input)
  - [Data Output](#data-output)
- [Troubleshooting and FAQs](#troubleshooting-and-faqs)
- [Other information](#other-information)
- [Legal and Compliance Information](#legal-and-compliance-information)
- [Updates and Release Notes](#updates-and-release-notes)

# Introduction

**MKRust** is a collection of high-performance Rust utilities that act as
companions to [mikrokondo](https://github.com/phac-nml/mikrokondo). The tools
are bundled into a single command-line binary (`mkrust`) exposing several
subcommands that handle common bioinformatics post-processing tasks in the
mikrokondo workflow, including:

- **FASTQ quality metrics** – summarise read-level and base-level quality
  statistics across one or more FASTQ files.
- **Consensus masking** – mask low-support positions in a consensus FASTA using
  a BAM pileup.
- **Mash hit parsing** – identify the top taxonomic hit from a Mash screen and
  detect potential metagenomic (mixed-taxa) samples.
- **Kraken2 report handling** – extract top taxonomic hits and bin contigs by
  taxonomic classification.

These utilities are written in Rust for speed and low memory overhead when
processing large sequencing datasets.

## Citation

If you use these tools as part of mikrokondo, please use the [mikrokonodo](https://github.com/phac-nml/mikrokondo) citation.

Otherwise you can cite this repository as:
Matthew Wells, Christy-Lynn Peterson, Github https://github.com/DOED-DAAD/MKRust


## Contact

[Matthew Wells] : <matthew.wells@phac-aspc.gc.ca>

# Install

MKRust is built with [Cargo](https://doc.rust-lang.org/cargo/)

```bash
# Clone the repository
git clone https://github.com/phac-nml/MKRust.git
cd MKRust

# Build an optimized release binary
cargo build --release

# The compiled binary will be located at:
./target/release/mkrust --help
```

Optionally install the binary onto your `PATH`:

```bash
cargo install --path .
```

## Compatibility

MKRust requires a recent stable Rust toolchain (edition 2024 features such as
`let`-chains and `with_added_extension` are used, so **Rust 1.80+** is
recommended).

> [!Note]
> Because `rust-htslib` links against **htslib**, a C compiler and the htslib
> system libraries must be available at build time.

# Getting Started

## Usage

MKRust is invoked through subcommands:

```bash
mkrust <SUBCOMMAND> [OPTIONS]
```

Available subcommands:

| Subcommand          | Description                                                           |
| ------------------- | -------------------------------------------------------------------- |
| `mask`              | Mask low-support positions in a FASTA using a BAM pileup.            |
| `fastq-quality`     | Compute per-file and combined FASTQ quality statistics.             |
| `mash-top-hit`      | Report the top taxonomic hit from a Mash screen.                    |
| `mash-metagenomic`  | Determine if a sample is metagenomic based on Mash hits.            |
| `kraken2-top-hit`   | Report the top Kraken2 hit at a given taxonomic level.             |
| `kraken2-bin`       | Bin contigs into separate FASTA files by taxonomic classification. |

**Examples**

Compute FASTQ quality statistics for reads and write to a file:

```bash
mkrust fastq-quality \
  --fastq sample_R1.fastq.gz sample_R2.fastq.gz sample_SE.fastq.gz \
  --names R1 R2 SE \
  --output sample_quality.json
```

Mask a consensus FASTA based on BAM support:

```bash
mkrust mask \
  --bam alignment.bam \
  --fasta consensus.fasta \
  --fwd-depth 1 \
  --rev-depth 1 \
  --pileup-support 0.75 \
  --min-depth 2
```

Get the top Mash hit above a 90% identity cutoff:

```bash
mkrust mash-top-hit --input mash_screen.tab --percent-identity-cutoff 0.90
```

Check whether a sample is metagenomic:

```bash
mkrust mash-metagenomic \
  --input mash_screen.tab \
  --equivalent-taxa equivalent_taxa.json \
  --taxonomic-level g \
  --percent-identity-cutoff 0.90
```

Extract the top Kraken2 hit at the species level:

```bash
mkrust kraken2-top-hit --report kraken2.report --taxa-level S
```

Bin contigs by species-level classification:

```bash
mkrust kraken2-bin \
  --report kraken2.report \
  --kraken-output kraken2.output \
  --contigs contigs.fasta \
  --prefix sample \
  --taxa-level S \
  --keep-parents
```

## Configuration and Settings

### `mask`

| Parameter           | Flag             | Default | Description                                                                 |
| ------------------- | ---------------- | ------- | --------------------------------------------------------------------------- |
| BAM file            | `-b, --bam`      | —       | Input BAM alignment path.                                                   |
| FASTA file          | `-f, --fasta`    | —       | Reference/consensus FASTA (bgzip supported, plain gzip is **not**).         |
| Forward depth       | `-w, --fwd-depth`| `1`     | Minimum read depth required from forward reads.                             |
| Reverse depth       | `-r, --rev-depth`| `1`     | Minimum read depth required from reverse reads.                            |
| Pileup support      | `-p, --pileup-support` | `0.75` | Minimum proportion (0–1) of pileup supporting the base to avoid masking. |
| Minimum depth       | `-m, --min-depth`| `2`     | Minimum total depth below which a position is masked.                       |

### `fastq-quality`

| Parameter    | Flag            | Default | Description                                              |
| ------------ | --------------- | ------- | ------------------------------------------------------- |
| FASTQ files  | `-f, --fastq`   | —       | One or more FASTQ files (plain or gzipped).            |
| Names        | `-n, --names`   | —       | Optional display names matching the FASTQ files order. |
| Output       | `-o, --output`  | stdout  | Optional output file; defaults to standard output.     |

### `mash-top-hit`

| Parameter              | Flag                          | Default | Description                       |
| ---------------------- | ----------------------------- | ------- | --------------------------------- |
| Input                  | `-i, --input`                 | —       | Mash screen output file.          |
| Percent identity cutoff| `-p, --percent-identity-cutoff` | `0.90` | Minimum identity to consider a hit. |

### `mash-metagenomic`

| Parameter              | Flag                          | Default | Description                                   |
| ---------------------- | ----------------------------- | ------- | --------------------------------------------- |
| Input                  | `-i, --input`                 | —       | Mash screen output file.                      |
| Equivalent taxa        | `-e, --equivalent-taxa`       | —       | JSON mapping of equivalent taxa.              |
| Taxonomic level        | `-t, --taxonomic-level`       | `g`     | Single-character taxonomic level to evaluate. |
| Percent identity cutoff| `-p, --percent-identity-cutoff` | `0.90` | Minimum identity to consider a hit.        |

### `kraken2-top-hit`

| Parameter    | Flag              | Default | Description                              |
| ------------ | ----------------- | ------- | ---------------------------------------- |
| Report       | `-r, --report`    | —       | Input Kraken2 report.                    |
| Taxa level   | `-t, --taxa-level`| `S`     | Taxonomic level to classify at.          |

### `kraken2-bin`

| Parameter      | Flag                | Default | Description                                         |
| -------------- | ------------------- | ------- | --------------------------------------------------- |
| Report         | `-r, --report`      | —       | Input Kraken2 report.                               |
| Kraken output  | `-k, --kraken-output` | —     | Input Kraken2 per-read/contig output.               |
| Contigs        | `-c, --contigs`     | —       | Input FASTA of contigs (indexed via faidx).         |
| Prefix         | `-p, --prefix`      | —       | Output filename prefix for binned FASTA files.      |
| Taxa level     | `-t, --taxa-level`  | `S`     | Taxonomic level to bin at.                          |
| Keep parents   | `--keep-parents`    | `false` | Also include contigs classified at a higher level.  |

## Data Input

- **FASTQ** (`fastq-quality`): standard FASTQ, optionally gzip-compressed
  (`.gz`). Quality scores are assumed to use Phred+33 encoding.
- **FASTA** (`mask`, `kraken2-bin`): indexed FASTA. An index (`.fai`) is built
  automatically if absent. For `mask`, only **bgzip** compression is supported
  (plain gzip is not).
- **BAM** (`mask`): sorted BAM alignments. A `.bai` index is built automatically
  if absent.
- **Mash screen output** (`mash-top-hit`, `mash-metagenomic`): tab-delimited
  Mash output with at least 6 fields (identity, shared-hashes, median
  multiplicity, p-value, query-id, query-note).
- **Equivalent taxa JSON** (`mash-metagenomic`): JSON object mapping a
  representative name to a list of semicolon-delimited taxonomy strings.
- **Kraken2 report / output** (`kraken2-top-hit`, `kraken2-bin`): standard
  Kraken2 report and read-classification output files.

## Data Output

### `fastq-quality`

Emits a JSON object with one entry per input file plus a `combined` entry.
Example fields:

| Field                   | Description                                      |
| ----------------------- | ------------------------------------------------ |
| `total_bp`              | Total number of bases.                           |
| `total_reads`           | Total number of reads.                           |
| `qual_min` / `qual_max` | Minimum / maximum observed Phred score.          |
| `qual_sum`              | Sum of all Phred scores.                         |
| `qual_mean` / `qual_std`| Mean / standard deviation of per-base quality.   |
| `read_qual_mean` / `read_qual_std` | Mean / std of per-read average quality. |
| `mean_sequence_length`  | Mean read length.                                |
| `min_sequence_length` / `max_sequence_length` | Min / max read length.     |
| `std_sequence_length`   | Standard deviation of read length.               |

### `mask`

Writes a masked FASTA (`<input>.masked.fasta`) where unsupported positions are
replaced with `N`. A per-position log of masked positions is written to
standard output.

### `mash-top-hit`

Prints the top taxonomic hit (or `No Species Identified`) to standard output.

### `mash-metagenomic`

Prints `true` if the sample appears metagenomic (more than one distinct taxon at
the selected level) or `false` otherwise.

### `kraken2-top-hit`

Prints the top taxon name at the selected level (or `No Species Identified`).

### `kraken2-bin`

Writes one gzip-compressed FASTA per taxon
(`<prefix>_<taxon>.binned.fasta.gz`) in the current working directory.

# Troubleshooting and FAQs

- **`Failed to build fasta index` / htslib errors** – ensure htslib system
  libraries are installed and that your FASTA uses bgzip (not plain gzip) for
  `mask`.
- **`Could not determine file extension`** – `fastq-quality` relies on the file
  extension (`.gz`) to decide whether to decompress; ensure files are named
  appropriately.
- **`Sequence in header is not present in bam file`** – the FASTA and BAM must
  share identical sequence names.

# Other information

# Legal and Compliance Information

Copyright Government of Canada 2026

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

# Updates and Release Notes

Please see the `CHANGELOG.md`.
