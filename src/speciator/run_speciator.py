import json
import logging
import sys
import tempfile
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as import_version
from pathlib import Path
from typing import Annotated, Tuple

import typer

from speciator.mash_functions import sketch_input
from speciator.search import SearchEnvironment, SearchResult, assign_species

app = typer.Typer(no_args_is_help=True)

searchConfig = SearchEnvironment.default_search_environment()

# Set up logging
logger = logging.getLogger(__name__)


@app.command()
def fastq(
    fastq_paths: Annotated[
        Tuple[Path, Path],
        typer.Argument(
            help="Path to the FASTQ files. Each file must be specified separately.",
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            allow_dash=True,
        ),
    ],
):
    run(searchConfig, [str(fastq_path) for fastq_path in fastq_paths])


@app.command()
def fasta(
    fasta_path: Annotated[
        Path,
        typer.Argument(
            help="Path to the FASTA file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            allow_dash=True,
        ),
    ],
):
    run(searchConfig, [str(fasta_path)])


def run(search_config: SearchEnvironment, files: list[str]):
    logger.info(f"Running with configuration: {search_config}")
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".msh", delete=False, prefix="speciator_"
    ) as tmp:
        mash_file_path = Path(tmp.name)
    try:
        sketch_input(search_config.mash, files, mash_file_path)
        logger.info(f"Created mash sketch file: {str(mash_file_path)}")
        result: SearchResult = assign_species(mash_file_path, searchConfig)
        logger.debug(f"Species assignment result: {result}")
    finally:
        if mash_file_path.exists():
            mash_file_path.unlink()
    output = {**result.__dict__, "speciator_version": get_version()}
    print(json.dumps(output, default=lambda o: o.__dict__), file=sys.stdout)
    logger.info("Finished running speciator")


@app.command()
def version():
    print(json.dumps({"version": get_version()}), file=sys.stdout)


@app.callback()
def main(
    config_file: Annotated[
        Path,
        typer.Option(
            "--config-file",
            "-c",
            help="Path to the configuration file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ] = searchConfig.config_file,
    library_location: Annotated[
        Path,
        typer.Option(
            "--library-location",
            "-l",
            help="Path to the library directory. By default it will be taken from the configuration file.",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
        ),
    ] = None,
    log_level: Annotated[
        str,
        typer.Option(
            "--log-level",
            "-L",
            help="Set the logging level",
            case_sensitive=False,
        ),
    ] = "WARNING",
    threads: Annotated[int, typer.Option(
        "--threads",
        "-t",
        help="Set the number of threads for parallelization. This may or may not help improve performance, depending.",
    )] = 1
):
    global searchConfig

    # Set up logging
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")
    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logger.info(f"Initializing with log level: {log_level}")

    searchConfig = SearchEnvironment.from_config_file(
        config_file, threads=threads, library=library_location
    )
    logger.info(f"SearchEnvironment initialized: {searchConfig}")


def get_version() -> str:
    """Get the package version from installed metadata."""
    try:
        return import_version("speciator")
    except PackageNotFoundError:
        return "unknown"


if __name__ == "__main__":
    app()
