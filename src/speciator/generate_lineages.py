import dataclasses
import logging
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from speciator.lineages import (
    _UNKNOWN_CODE,
    _UNKNOWN_SPECIES_CODE,
    TaxonKit,
)  # Import the necessary classes

app = typer.Typer(no_args_is_help=True)
logger = logging.getLogger(__name__)

@app.command()
def main(
        species_ids_file: Annotated[
            Path,
            typer.Argument(
                help="Path to a file containing species IDs, one per line.",
                exists=True,
                file_okay=True,
                dir_okay=False,
                readable=True,
                resolve_path=True,
            ),
        ],
        output_file: Annotated[
            Path,
            typer.Option(
                "--output", "-o",
                help="Path to the output Parquet file (e.g., lineages.pqt).",
                file_okay=True,
                dir_okay=False,
                writable=True,
                resolve_path=True,
            ),
        ] = Path("lineages.pqt"),
        taxonkit_path: Annotated[
            str,
            typer.Option(
                "--taxonkit-path",
                "-t",
                help="Path to the taxonkit executable (or 'taxonkit' if on PATH).",
            ),
        ] = "taxonkit",
        log_level: Annotated[
            str,
            typer.Option(
                "--log-level",
                "-L",
                help="Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).",
                case_sensitive=False,
            ),
        ] = "INFO",
) -> None:
    """
    Generates a Parquet file containing lineage information for a given list of species IDs.

    This script reads a list of species IDs from an input file, uses the 'taxonkit'
    tool to retrieve their full lineage information, and then saves this data
    into a highly efficient Parquet format.
    """
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise typer.BadParameter(f"Invalid --log-level: {log_level}")
    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger.info(f"Starting lineage generation for species IDs in {species_ids_file}")
    species_ids = set()
    try:
        with open(species_ids_file, "r") as f:
            for line in f:
                stripped_line = line.strip()
                if stripped_line: # Only add non-empty lines
                    species_ids.add(stripped_line)
    except FileNotFoundError:
        logger.error(f"Input file not found: {species_ids_file}")
        raise typer.Exit(code=1)

    if not species_ids:
        logger.warning("No species IDs found in the input file. Creating an empty output file.")
        df = pl.DataFrame()
        df.write_parquet(output_file)
        logger.info(f"Empty Parquet file created at {output_file}")
        return

    logger.info(f"Found {len(species_ids)} unique species IDs.")

    taxonkit_instance = TaxonKit(taxonkit_path)
    logger.info(f"Extracting lineages using TaxonKit. Ensure 'taxonkit' is installed and accessible at '{taxonkit_path}'...")

    try:
        # The extract_lineages method expects a set of strings as taxon_ids
        lineages_dict = taxonkit_instance.extract_lineages(species_ids)
    except Exception as e:
        logger.error(f"Error extracting lineages with TaxonKit: {e}")
        logger.error("Please ensure 'taxonkit' is correctly installed and configured.")
        raise typer.Exit(code=1)

    if not lineages_dict:
        logger.warning("No lineage information could be extracted for the provided IDs. Creating an empty output file.")
        df = pl.DataFrame()
        df.write_parquet(output_file)
        logger.info(f"Empty Parquet file created at {output_file}")
        return

    if len(lineages_dict) != len(species_ids):
        extra = {ext for ext in lineages_dict.keys() if ext not in species_ids}
        if extra and next(iter(extra)) not in [_UNKNOWN_CODE, _UNKNOWN_SPECIES_CODE]:
            logger.warning(f"Number of lineages does not match the number of species IDs. Expected: {len(species_ids)}, found: {len(lineages_dict)}")
            missing = {exp for exp in species_ids if exp not in lineages_dict}
            logger.warning(f"Missing: {missing}")
            logger.warning(f"Extra: {extra}")
            raise typer.Exit(code=1)

    # Convert the dictionary of Lineage objects to a list of dictionaries for Polars
    # The Lineage dataclass maps directly to column names.
    lineage_records = [dataclasses.asdict(lineage) for lineage in lineages_dict.values()]

    # Create a Polars DataFrame
    df = pl.DataFrame(lineage_records)

    logger.info(f"Writing {len(lineage_records)} lineage records to Parquet file: {output_file}")
    df.write_parquet(output_file)
    logger.info("Lineage data successfully written to Parquet.")

if __name__ == "__main__":
    app()