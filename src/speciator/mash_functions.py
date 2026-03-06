import csv
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import uuid
from collections.abc import Iterator
from io import StringIO
from pathlib import Path
from typing import TextIO

import polars as pl

logger = logging.getLogger(__name__)


def sketch_stream(
    mash: str, query_name, extra_args: list[str] = None, source: TextIO = sys.stdin
) -> None:
    if extra_args is None:
        extra_args: list[str] = []

    with subprocess.Popen(
        [mash, "sketch", "-o", query_name, "-"] + extra_args,
        stdin=subprocess.PIPE,
        universal_newlines=True,
    ) as shell:
        stdout, stderr = shell.communicate(source.read())
        if shell.returncode != 0:
            raise Exception(stderr)


def get_best_mash_matches(
    sample_sketch,
    ref_seq_sketch,
    mash_path="",
    number_of_best_matches=10,
    distance_threshold=0.5,
    threads=None,
):
    """
    run mash dist for the sample sketch file vs the specified library and return the best matches
    """
    mash_dists = mash_polars(
        mash_path, ref_seq_sketch, sample_sketch, distance_threshold, threads
    )
    # filter columns
    mash_dists = mash_dists.select(
        ["query", "subject", "distance", "p-value", "shared-hashes"]
    )
    # sort by distance and output the subjects (match in refseq)
    matches = (
        mash_dists.sort("distance", descending=False)
        .head(number_of_best_matches)
        .rename({"query": "filename"})
        .select(["filename", "distance", "p-value", "shared-hashes"])
    )
    return matches


def execute_mashing(
    mash_path,
    reference_sketch,
    sample_sketch,
    distance_threshold,
    threads: int  = 1,
) -> Iterator[str]:
    """
    Run `mash dist` and yield stdout lines.

    Sorting is intentionally not done here; consumers can sort results as needed.
    """
    threads = threads if threads is not None else 1
    logger.debug(f"Running mash with {threads} threads")
    mash_command = [
        "mash" if mash_path == "mash" else os.path.join(mash_path, "mash"),
        "dist",
        "-d",
        str(distance_threshold),
        "-p",
        str(threads),
        reference_sketch,
        sample_sketch,
    ]

    try:
        mash_process = subprocess.Popen(mash_command, stdout=subprocess.PIPE, text=True)

        for line in mash_process.stdout:
            yield line

        mash_process.wait()

        if mash_process.returncode != 0:
            raise subprocess.CalledProcessError(mash_process.returncode, mash_command)

    except subprocess.CalledProcessError as e:
        logger.error(e.output)
        logger.error(
            f"Error executing mash dist for {sample_sketch.replace('.msh', '')}",
        )
        raise e


def mash_polars(mash_path, ref_seq_sketch, sample_sketch, distance_threshold, threads):
    mash_output = execute_mashing(
        mash_path,
        ref_seq_sketch,
        sample_sketch,
        distance_threshold,
        threads
    )

    # Collect all output into a list
    mash_output_lines = list(mash_output)

    # Create an empty DataFrame with the correct columns if there's no output
    if not mash_output_lines:
        logger.info(f"No output from mash for sketch {ref_seq_sketch}")
        logger.info(f"mash_polars rows returned: {0}")
        return pl.DataFrame(
            schema={
                "query": pl.Utf8,
                "subject": pl.Utf8,
                "distance": pl.Float64,
                "p-value": pl.Float64,
                "shared-hashes": pl.Utf8,
            }
        )

    mash_output_str = "".join(mash_output_lines)

    try:
        mash_dists = pl.read_csv(
            StringIO(mash_output_str),
            separator="\t",
            has_header=False,
            new_columns=["query", "subject", "distance", "p-value", "shared-hashes"],
        ).sort("distance", descending=False)

        logger.info(f"mash_polars rows returned: {mash_dists.height}")
        return mash_dists
    except Exception as e:
        logger.error(f"Error parsing mash output: {str(e)}")
        logger.debug(f"Mash output: {mash_output_str}")
        raise


def sketch_file(
    mash_path: str, fasta_gz: Path | str, filename: Path, record_name: str = None
) -> None:
    if not record_name:
        record_name = filename.name.replace(".msh", "")
    subprocess.check_output(
        f'{mash_path} sketch -o "{str(filename)}" -I {record_name} "{fasta_gz}"',
        shell=True,
    )


def sketch_input(
    mash_path: str, files: list[str], temp_file: Path = Path(f"{str(uuid.uuid4())}.msh")
) -> Path:
    extra_args: list[str] = ["-m", "3"] if len(files) > 1 else []
    if "-" in files:
        for file in files:
            if file != "-":
                logger.error(
                    "FASTQs can either be read from file or streamed from stdin, not a mixture of the two."
                )
                exit(1)
        sketch_stream(mash_path, temp_file, extra_args)
    else:
        command = f'cat {" ".join(['"' + file + '"' for file in files])} | "{mash_path}" sketch -o "{temp_file}" {" ".join(extra_args)} -'
        subprocess.check_output(command, shell=True)
    return temp_file


def paste_sketches(
    mash_path: str, sketch_files: list[Path | str], output_file: Path
) -> None:
    if output_file.exists():
        logger.warning(f"{output_file} already exists, overwriting.")
        output_file.unlink()
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
        for sketch in sketch_files:
            temp_file.write(f"{str(sketch)}\n")
        temp_file_name = temp_file.name
    try:
        command = f'"{mash_path}" paste "{str(output_file)}" -l "{temp_file_name}"'
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing mash paste: {str(e)}")
        logger.error(f"Copying temporary file {temp_file_name} to current directory")
        shutil.copy(temp_file_name, "./")
        raise e
    finally:
        os.remove(temp_file_name)


def self_search(
    mash_path: str, sketches: str, threads, distance_threshold: float = 0.5
) -> Iterator[list[str]]:
    """Returns scores line-by-line sorted"""
    scores = execute_mashing(mash_path, sketches, sketches, distance_threshold, threads)
    reader = csv.reader(scores, delimiter="\t")
    for row in reader:
        yield row


def _extract_single_record_name_from_mash_info(output: str) -> str | None:
    """
    Best-effort parser for `mash info <sketch.msh>` output.

    For single-record sketches, Mash typically prints a record line starting with '1000'
    (possibly with leading spaces). We return the 3rd token (the record name).
    """
    for line in output.splitlines():
        s = line.lstrip()
        if not s.startswith("1000"):
            continue
        parts = s.split()
        if len(parts) >= 3:
            return parts[2]
    return None


def mash_record_id(mash_path: str, sketch_path: Path) -> str:
    """
    Return the sketch record identifier as Mash reports it (via `mash info`).

    If parsing fails, falls back to `sketch_path.stem`.
    """
    try:
        proc = subprocess.run(
            [mash_path, "info", str(sketch_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.warning(
            "Failed to run '%s info' for %s: %s", mash_path, sketch_path, str(e)
        )
        return sketch_path.stem

    record_name = _extract_single_record_name_from_mash_info(proc.stdout)
    if record_name is None:
        logger.warning(
            "Could not parse mash record id from '%s info' output for %s",
            mash_path,
            sketch_path,
        )
        return sketch_path.stem
    return record_name
