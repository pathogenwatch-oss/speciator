#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
#     "typer",
# ]
# ///

import csv
import json
from pathlib import Path
from typing import Annotated, Optional

import polars as pl
import typer

app = typer.Typer(no_args_is_help=True)


@app.command()
def p2c(
    parquet_file: Annotated[
        Path,
        typer.Argument(
            help="Path to the Parquet format file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
        ),
    ],
):
    df = pl.read_parquet(parquet_file)
    df.write_csv(f"{parquet_file.stem}.csv")


@app.command()
def c2p(
    csv_file: Annotated[
        Path,
        typer.Argument(
            help="Path to the TSV/CSV file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
        ),
    ],
    separator: Annotated[
        str,
        typer.Option(
            "--separator", "-s", help="Specify a separator for the CSV/TSV file."
        ),
    ] = ",",
    schema: Annotated[
        Optional[str],
        typer.Option(
            "--schema",
            help='JSON string defining the schema for the output Parquet file. For example: `python format_converter.py c2p your_file.csv --schema \'{"record_key": "Utf8", "organism_code": "Int64", "accession": "Utf8", "represents": "Int64"}',
        ),
    ] = None,
):
    # Read with Python's csv module to preserve quoted strings
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=separator)
        data = list(reader)
    
    # Convert to Polars DataFrame - this preserves string types for quoted fields
    df = pl.DataFrame(data)
    
    # # Now apply smart type inference that respects the original string nature
    # for col in df.columns:
    #     series = df[col]
    #     # Try to convert to numeric only if ALL values can be converted
    #     # and they weren't originally quoted strings from CSV
    #     try:
    #         # Check if all non-null values are numeric
    #         numeric_series = series.cast(pl.Int64, strict=False)
    #         if numeric_series.null_count() == series.null_count():
    #             # All values successfully converted, use numeric type
    #             df = df.with_columns(numeric_series.alias(col))
    #     except:
    #         # Keep as string if conversion fails
    #         pass

    if schema:
        try:
            schema_dict = json.loads(schema)
            typed_schema = {k: getattr(pl, v) for k, v in schema_dict.items()}
            df = df.select(
                [pl.col(col).cast(dtype) for col, dtype in typed_schema.items()]
            )
        except json.JSONDecodeError:
            typer.echo("Error: Invalid JSON schema provided", err=True)
            raise typer.Exit(code=1)
        except AttributeError as e:
            typer.echo(f"Error: Invalid data type in schema: {str(e)}", err=True)
            raise typer.Exit(code=1)

    df.write_parquet(f"{csv_file.stem}.pqt")
    typer.echo(f"Parquet file created: {csv_file.stem}.pqt")


if __name__ == "__main__":
    app()
