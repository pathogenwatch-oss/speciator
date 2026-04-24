from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from speciator.mash_functions import paste_sketches

app = typer.Typer(no_args_is_help=True)
logger = logging.getLogger(__name__)


def _resolve_sketch_for_row(
    row: dict[str, str], sketches_dir: Path
) -> Path | None:
    candidates = [row.get("record_key"), row.get("accession")]
    for value in candidates:
        if not value:
            continue
        sketch = sketches_dir / f"{value}.msh"
        if sketch.exists():
            return sketch
    return None


@app.command()
def split(
    library_name: Annotated[
        str,
        typer.Argument(
            help="Logical library name (e.g. Bacteria for Bacteria.msh/Bacteria.pqt)."
        ),
    ],
    max_sub_library_size: Annotated[
        int,
        typer.Option(
            "--max-sub-library-size",
            "-s",
            min=1,
            help="Maximum records per output shard.",
        ),
    ] = 25000,
    library_dir: Annotated[
        Path,
        typer.Option(
            "--library-dir",
            "-l",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            help="Directory containing the unsplit library files.",
        ),
    ] = Path("library"),
    build_space_dir: Annotated[
        Path,
        typer.Option(
            "--build-space-dir",
            "-b",
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            help="Build-space directory containing mash sketches at <build_space>/mash.",
        ),
    ] = Path("build_space"),
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output-dir",
            "-o",
            file_okay=False,
            dir_okay=True,
            help="Directory to write sharded files. Defaults to --library-dir.",
        ),
    ] = None,
    mash_path: Annotated[
        str,
        typer.Option(
            "--mash-path",
            "-m",
            help="Path to mash executable if not on PATH.",
        ),
    ] = "mash",
    overwrite: Annotated[
        bool,
        typer.Option(
            "--overwrite",
            help="Overwrite existing shard outputs if they already exist.",
        ),
    ] = False,
    keep_original: Annotated[
        bool,
        typer.Option(
            "--keep-original",
            help="Keep original <library>.msh/.pqt when writing shards in-place.",
        ),
    ] = False,
) -> None:
    numeric_level = getattr(logging, "INFO", None)
    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    output_dir = output_dir if output_dir is not None else library_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    sketches_dir = build_space_dir / "mash"
    if not sketches_dir.exists():
        raise FileNotFoundError(f"Sketch directory not found: {sketches_dir}")

    source_msh = library_dir / f"{library_name}.msh"
    source_pqt = library_dir / f"{library_name}.pqt"
    if not source_msh.exists() or not source_pqt.exists():
        raise FileNotFoundError(
            f"Missing source library files. Expected {source_msh} and {source_pqt}"
        )

    records = pl.read_parquet(source_pqt)
    if records.height == 0:
        raise ValueError(f"No rows found in {source_pqt}")
    if records.height <= max_sub_library_size:
        logger.info(
            f"{library_name} has {records.height} records (<= {max_sub_library_size}); no split required."
        )
        return

    rows = records.to_dicts()
    sketch_paths: list[Path] = []
    missing: list[str] = []
    for row in rows:
        sketch = _resolve_sketch_for_row(row, sketches_dir)
        if sketch is None:
            missing.append(str(row.get("record_key") or row.get("accession") or row))
        else:
            sketch_paths.append(sketch)

    if missing:
        preview = ", ".join(missing[:10])
        raise FileNotFoundError(
            f"Missing {len(missing)} sketch files in {sketches_dir}. First missing: {preview}"
        )

    if len(sketch_paths) != records.height:
        raise ValueError(
            f"Resolved {len(sketch_paths)} sketches for {records.height} records"
        )

    n_shards = (records.height + max_sub_library_size - 1) // max_sub_library_size
    logger.info(
        f"Splitting {library_name} ({records.height} records) into {n_shards} shards of <= {max_sub_library_size}"
    )
    for shard_idx, start in enumerate(
        range(0, records.height, max_sub_library_size), start=1
    ):
        end = min(start + max_sub_library_size, records.height)
        shard_msh = output_dir / f"{library_name}.{shard_idx}.msh"
        shard_pqt = output_dir / f"{library_name}.{shard_idx}.pqt"
        if not overwrite and (shard_msh.exists() or shard_pqt.exists()):
            raise FileExistsError(
                f"Shard already exists: {shard_msh if shard_msh.exists() else shard_pqt}. Use --overwrite to replace."
            )
        paste_sketches(mash_path, sketch_paths[start:end], shard_msh)
        records.slice(start, end - start).write_parquet(shard_pqt)
        logger.info(f"Wrote {shard_msh.name} + {shard_pqt.name} ({end - start} rows)")

    in_place = output_dir.resolve() == library_dir.resolve()
    if in_place and not keep_original:
        archived_msh = library_dir / f"{library_name}_full.msh"
        archived_pqt = library_dir / f"{library_name}_full.pqt"
        if not overwrite and (archived_msh.exists() or archived_pqt.exists()):
            raise FileExistsError(
                f"Archive targets already exist: {archived_msh} / {archived_pqt}. Use --overwrite or remove them."
            )
        if archived_msh.exists():
            archived_msh.unlink()
        if archived_pqt.exists():
            archived_pqt.unlink()
        source_msh.rename(archived_msh)
        source_pqt.rename(archived_pqt)
        logger.info(
            f"Archived original library to {archived_msh.name} and {archived_pqt.name}"
        )
    elif in_place and keep_original:
        logger.warning(
            "Original unsplit library remains active and will be preferred by search over shards."
        )


if __name__ == "__main__":
    app()
