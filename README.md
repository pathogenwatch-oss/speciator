# Speciator

## About

Speciator determines the species of microbial genome FASTA (assembled) or FASTQ (raw sequence) files by
using [mash](https://github.com/marbl/Mash) to compare it to a hierarchical library of references.

Speciator is both fast and accurate, with a flexible library architecture that incorporates multiple sources of curated
representatives and automated extraction from GenBank.

Speciator consists of two core programs:

1. The library build tool. This provides automated builds of the Speciator library, along with builds of individual
   components.
2. The genome search tool. Speciate FASTA or FASTQ files, with output provided in JSON format.

Both can be run locally using python or Docker. If running outside Docker, it also requires `mash` to be installed.

Speciator has three primary commands:

1. `speciator fasta <fasta file>` - Search a FASTA file
2. `speciator fastq <fastq file 1> <fastq file 2>` - Search a FASTQ file
3. `speciator build` - Build the library

## Installing

Clone the repository to your local disk and then follow the instructions for one of:

1. [pixi](#using-pixi) - install into an anaconda workspace
2. [uv](#using-uv) - run directly or install as a python module
3. [docker](#using-docker) – create an image for distribution

The library must be built before running Speciator for the firstm time. This can be done directly, e.g. with `uv`, or by
using Docker to create a code image and running the code image.

### Using pixi

[Pixi]() will install all dependencies and create an Anaconda environment with Speciator installed.

#### Installing with pixi

```bash
pixi install
pixi shell
```

#### Building the library

```bash
speciator-build-lib -o library full v3.2.4
```

#### Using speciator

```bash
speciator fasta query.fasta
```

### Using uv

With [uv]() installed, it's possible to run speciator and the library builder scripts directly from within the
directory or install the module into a virtual environment or system environment.

#### Running directly with uv

```bash
# Build the library
uv run speciator-build-lib -o library full v3.2.4

# Use speciator
uv run speciator fasta query.fasta
```

#### Installing with uv

```bash
uv pip install .
```

#### Building the library

```bash
speciator-build-lib -o library full v3.2.4
```

#### Using speciator

```bash
speciator fasta query.fasta
```

### Using Docker

There are two Docker build targets available:

1. `speciator` - Standard interface with full file support (requires mounting)
2. `pathogenwatch` - Pathogenwatch-style interface (STDIN/STDOUT, plain text FASTA only)

#### Building the Docker images

```bash
# Build the code image
docker build --rm --pull --target speciator -t speciator-code:v5.0.0 .

# Build the Pathogenwatch image
docker build --rm --pull --target pathogenwatch -t speciator:v5.0.0 .
```

#### Building the library

```bash
# Mount paths for output and build cache (recommended for caching intermediate steps)
docker run --rm -v /path/to/final/:/final -v /path/to/build_space/:/build_space speciator-code:v5.0.0 -b /build_space -o /final full v3.2.4
```

#### Using the Docker images

```bash
# Using the standard interface
docker run --rm -v /path/to/data:/data speciator-code:v5.0.0 fasta /data/query.fasta

# Using the Pathogenwatch interface (STDIN/STDOUT, plain text FASTA only)
cat query.fasta | docker run --rm -i speciator:v5.0.0 > my_result.json

# For gzipped files, decompress first
gunzip -c query.fasta.gz | docker run --rm -i speciator:v5.0.0 > my_result.json
```

### Using Singularity

There are three Singularity definition files available:

1. `speciator.def` - Full image with a pre-built library included
2. `speciator-lite.def` - Lightweight image requiring external library mount
3. `speciator-library-builder.def` - Image for building the reference library

#### Building the Singularity images

```bash
# Build the library builder image
singularity build --fakeroot speciator-library-builder.sif speciator-library-builder.def

# Build the library (mount paths for output and build cache)
singularity run -B /path/to/final:/final -B /path/to/build_space:/build_space speciator-library-builder.sif -b /build_space -o /final full v3.2.4

# Build the lite image (requires external library)
singularity build --fakeroot speciator-lite.sif speciator-lite.def

# Build the full image (includes library - must be built after creating the library)
singularity build --fakeroot speciator.sif speciator.def
```

#### Using the Singularity images

```bash
# Using the lite image (requires library mount)
singularity run -B /path/to/library:/speciator/library speciator-lite.sif fasta query.fasta

# Using the full image (library included)
singularity run speciator.sif fasta query.fasta
```

## The library builder

To build the complete library, use the run_library_builder project script, and provide the `full` command along with
the version of `kleborate` to incorporate. The `speciator-build-lib` script also accepts `--log-level`. For complete
options see [below](#library-build-commands-and-options).

e.g. with Kleborate

```bash
> speciator-build-lib -L debug full v3.2.4 2>&1 | tee logs.out
```

To limit NCBI extraction by assembly level, use `--ncbi-extraction-level` (or `-x`) with one of `complete`,
`chromosome`, `scaffold`, or `contig`.

### Library build commands and options

`speciator-build-lib` is a Typer CLI with global options and subcommands.

Global options (apply to all subcommands):

- `-c, --config-file` (default: `config.toml`) TOML config with thresholds.
- `-m, --mash-path` (default: `mash`) path to the `mash` executable.
- `-t, --taxonkit-path` (default: `taxonkit`) path to the `taxonkit` executable.
- `-b, --build-dir` (default: `build_space`) working directory for intermediate files.
- `-o, --out-dir` (default: `final`) output directory for final libraries.
- `-L, --log-level` (default: `WARNING`) logging level.
- `-d, --download-threads` (default: `2`) download/sketch worker threads.

Subcommands:

- `kleborate [--version/-v VERSION] [DB_NAME]` (default version: `v3.2.4`, db name: `Kleborate`) downloads and builds
  the Kleborate library.
- `curated [DB_NAME] [--clean]` (default db name: `Curated`) imports curated libraries and optionally cleans added
  files.
- `flu [--batch-size N] [--scaling-factor F] [--max-reps N] [DB_NAME]` (defaults: `500`, `0.005`, `2000`, db name:
  `Influenza`) builds the influenza library.
-

`ncbi [--update-metadata] [--clean] [--qc-metrics FILE] [--batch-size N] [--cluster-threshold F] [--scaling-factor F] [--max-reps N] [--skip-unclassified] [--assembly-level/-A LEVEL]` (
defaults: `filtered_metrics.csv`, `500`, `0.0000001`, `0.005`, `2000`, `contig`) builds the NCBI-derived libraries for
bacteria, archaea, fungi, viruses.

-

`full KLEBORATE_VERSION [--clean/-c] [--qc-metrics/-q FILE] [--batch-size/-b N] [--ncbi-scaling-factor/-s F] [--ncbi-max-reps/-m N] [--ncbi-extraction-level/-x LEVEL] [--flu-batch-size/-B N] [--flu-scaling-factor/-S F] [--flu-max-reps/-M N] [--curated-name NAME]`
orchestrates `kleborate`, `curated`, `ncbi`, and `flu` in one run.

## The search tool

`speciator` is a Typer CLI with global options and subcommands.

Global options (apply to all subcommands):

- `-c, --config-file` (default: `config.toml`) configuration file.
- `-l, --library-location` path to the library directory; overrides the config file value.
- `-L, --log-level` (default: `WARNING`) logging level.
- `-t, --threads` (default: `1`) thread count for search parallelisation.
- `--install-completion` install shell completion for the current shell.
- `--show-completion` print shell completion for the current shell.

Subcommands:

- `fasta FASTA` search a single FASTA file (plain or gzipped; `-` allowed for STDIN).
- `fastq FASTQ_1 FASTQ_2` search a paired FASTQ (plain or gzipped; `-` allowed for STDIN).
- `version` print the installed speciator version as JSON.

### Search a FASTA file

Run the speciator script and provide the `fasta` command along with the path to the (optionally gzipped) FASTA file.

```
speciator fasta ESC_EB8749AA_AS.fasta.gz
```

### Search a FASTQ file

Run the speciator script and provide the `fastq` command along with the paths both of to the (optionally gzipped) FASTQ
files.

```
speciator fastq ESC_EB8749AA_AS_1.fastq.gz ESC_EB8749AA_AS_2.fastq.gz
```

### Search using an explicit library directory

Override the library location from the config file.

```
speciator --library-location /path/to/library fasta query.fasta
```

### Check the installed version

```
speciator version
```

## Configuring Speciator

!TODO: Describe the config file

## How it works

### Library Structure

The library contains a set of sub-libraries and a summary file containing all the species lineages found in the
sub-libraries: `lineages.pqt`

Each sub-library consists of two files:

1. the compiled mash library used to search the query.
2. a parquet-format file containing the accession and species code for each representative.

The sub-libraries are searched in order. If a match is obtained for one, then the script exits and prints the result.

#### Curated

#### Kleborate

#### Flu

#### RefSeq/GenBank

### The algorithm

1. The query sequence is compiled into a mash profile.
2. It is then compared to each library in the order specified in the config file. The default order is:
    1. Curated,
    2. Kleborate
    3. Bacteria
    4. Influenza
    5. Viruses
    6. Fungi
    7. Archaea
3. If a match is found to a library, then the search is stopped and the result reported.
4. If no match is found ...
5. TODO

### The match process

Instead of simply taking the best matching genome, Speciator uses
the [Bactinspector](https://gitlab.com/antunderwood/bactinspector) method to select the consensus best match amongst the
top close hits. For some genera, where species classifications are unclear, or they have been recently redefined, there
can be mislabelled RefSeq genomes (e.g. amongst the Klebsiella). The consensus approach has proven robust against such
errors.

## Helper scripts

### Installed with speciator

These tools are all included if the speciator module is installed in the environment, or can be run using `uv run`.

#### speciator-select-representatives

This tool extracts a set of representative sketches from a collection of FASTA/FASTQ items. It uses a two-phase
selection process: Phase A performs strict initial clustering, and Phase B iteratively relaxes the threshold until a
target maximum number of representatives is reached. Its primary intended use is for generating a representative
library for a species.

Example:

```bash
speciator-select-representatives [OPTIONS] ITEMS...
```

#### speciator-lineages

Generates a Parquet file containing full NCBI lineage information for a list of species taxon IDs using taxonkit.

Example:

```bash
speciator-lineages species_ids.txt --output lineages.pqt --taxonkit-path /usr/bin/taxonkit
```

#### speciator-cat-parquet

A utility to concatenate multiple lineage Parquet files into one, automatically removing duplicate rows. Useful when
combining lineages from different sources or parallel build runs.

Example:

```bash
speciator-cat-parquet lineages_bacteria.pqt lineages_viruses.pqt --output lineages_combined.pqt
```

### Generic

This is a standalone script that can be moved outside the project and run anywhere using `uv`.

#### format_converter.py

Helper CLI to convert between Parquet (`.pqt`) and CSV/TSV using Polars. It is an `uv` script, so you can run it with
`uv run --script utils/format_converter.py ...`.

Commands:

- `p2c PARQUET_FILE` writes `PARQUET_FILE.stem.csv` in the current directory.
- `c2p CSV_FILE [--separator/-s SEP] [--schema JSON]` writes `CSV_FILE.stem.pqt` in the current directory.

`c2p` schema format: a JSON object of `{"column": "PolarsType"}` (e.g. `Utf8`, `Int64`, `Float64`).

Examples:

```bash
uv run --script utils/format_converter.py p2c lineages.pqt
uv run --script utils/format_converter.py c2p lineages.csv --separator ',' \
  --schema '{"record_key": "Utf8", "organism_code": "Int64", "accession": "Utf8", "represents": "Int64"}'
```

### Acknowledgements

The current version of speciator was written by Corin Yeats and draws on the Bactinspector method developed by Anthony
Underwood. The species QC thresholds used for selecting references were provided by Nabil-Fareed Alikhan.

Khalil Abudahab wrote the original version of speciator.
