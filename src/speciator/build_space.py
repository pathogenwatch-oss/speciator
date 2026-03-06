import dataclasses
import os
from collections.abc import Iterable
from pathlib import Path

import polars as pl

from speciator.lineages import Lineage


@dataclasses.dataclass
class RecordInfo:
    """Information for a record in the output mash databases."""
    record_key: str
    organism_code: str
    accession: str
    represents: int


@dataclasses.dataclass
class BuildSpace:
    location: Path
    output_location: Path

    def __post_init__(self):
        self.fasta_dir = self.location / "fasta"
        self.library_dir = self.location / "library"
        self.sketch_dir = self.location / "mash"
        self.build_dir = self.location / "build"
        self.lineages_path = self.library_dir / "lineages.pqt"
        os.makedirs(self.fasta_dir, exist_ok=True)
        os.makedirs(self.library_dir, exist_ok=True)
        os.makedirs(self.sketch_dir, exist_ok=True)
        os.makedirs(self.build_dir, exist_ok=True)
        os.makedirs(self.output_location, exist_ok=True)

    def __str__(self) -> str:
        return f"BuildSpace(location={self.location}, output_location={self.output_location}): {self.fasta_dir}, {self.library_dir}, {self.sketch_dir}, {self.build_dir}, {self.lineages_path}"

    @staticmethod
    def default():
        return BuildSpace(Path("build_space"), Path("."))

    def library_path(self, library_name: str) -> tuple[Path,Path]:
        return self.library_dir / f"{library_name}.msh", self.library_dir / f"{library_name}.pqt"

    def finalise_library(
        self,
        library_name: str,
        mash_archive: Path | str,
        records: list[RecordInfo],
        lineages: dict[str, Lineage],
    ) -> None:
        """Finalise a library."""
        library_mash_path, library_desc_path = self.library_path(library_name)
        os.rename(mash_archive, library_mash_path)
        pl.from_dicts([dataclasses.asdict(record) for record in records]).write_parquet(
            library_desc_path
        )
        self.update_lineage_file(lineages.values())

    def update_lineage_file(self, lineages: Iterable[Lineage]) -> None:

        if not os.path.exists(self.lineages_path):
            pl.from_dicts(
                [dataclasses.asdict(record) for record in lineages]
            ).write_parquet(self.lineages_path)
        else:
            pl.concat(
                [
                    pl.read_parquet(self.lineages_path),
                    (
                        pl.from_dicts(
                            [dataclasses.asdict(record) for record in lineages]
                        )
                    ),
                ]
            ).unique().write_parquet(self.lineages_path)
