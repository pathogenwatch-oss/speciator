from pathlib import Path
from typing import List

import polars as pl
import typer

app = typer.Typer()

@app.command()
def main(
        input_files: List[Path] = typer.Argument(..., help="List of lineage.pqt files to concatenate"),
        output: Path = typer.Option("lineages.pqt", "--output", "-o", help="Path for the output Parquet file")
):
    """
    Concatenate multiple lineage Parquet files and remove duplicate rows.
    """
    existing_files = [f for f in input_files if f.exists()]

    if not existing_files:
        typer.echo("No existing input files found.", err=True)
        raise typer.Exit(code=1)

    # Load and concatenate
    df = pl.concat([pl.read_parquet(f) for f in existing_files])

    # Remove duplicates
    unique_df = df.unique()

    unique_df.write_parquet(output)
    typer.echo(f"Successfully concatenated {len(existing_files)} files into {output}")
    typer.echo(f"Total unique rows: {unique_df.height}")

if __name__ == "__main__":
    app()