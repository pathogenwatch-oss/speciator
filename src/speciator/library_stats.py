from __future__ import annotations

import dataclasses
from pathlib import Path

import polars as pl
import typer

app = typer.Typer()


@dataclasses.dataclass
class LibraryStats:
    name: str
    counts: dict[str, int] = dataclasses.field(default_factory=dict)

    @staticmethod
    def extract_stats(metadata_file: Path) -> LibraryStats:
        name = metadata_file.stem
        data = pl.read_parquet(metadata_file)
        counts = data.group_by("organism_code").len().to_dict(as_series=False)
        counts_dict = dict(zip(counts["organism_code"], counts["len"]))
        return LibraryStats(name=name, counts=counts_dict)


@app.command()
def main(stats_file: Path) -> None:
    """This is a utility tool for exploring the libraries."""
    stats = LibraryStats.extract_stats(stats_file)
    print("organism_code,count")
    for organism, count in stats.counts.items():
        print(f"{organism},{count}")


if __name__ == "__main__":
    app()
