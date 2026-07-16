from __future__ import annotations

import dataclasses
import logging
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated, Any, Optional

import polars as pl
import toml
import typer
from tenacity import RetryError

from speciator.build_space import BuildSpace, RecordInfo
from speciator.kleborate import download_kleborate, extract_kleborate_lineages
from speciator.lineages import Lineage, TaxonKit
from speciator.linking import Linker
from speciator.mash_functions import paste_sketches, sketch_file
from speciator.metrics import Metric, read_metrics
from speciator.ncbi_utils import (
    Genome,
    get_exemplars,
    mash_distances,
    read_genbank_flu,
    read_kingdom_metadata,
)

app = typer.Typer(no_args_is_help=True)

# Set up logging
logger = logging.getLogger(__name__)

mash: str = "mash"
taxonkit: str = "taxonkit"
dwnld_threads: int = 2
mash_threads: int = os.cpu_count()
max_sub_library_size: int = 0
data_paths: BuildSpace = BuildSpace.default()
FLU_SPECIES = ["2955744", "2955935", "2955291", "2809131", "2955465"]
config: dict[str, Any] = {}


def move_files(source: Path, basename: str, dest: Path):
    for item in source.iterdir():
        if item.is_file() and item.name.startswith(basename):
            dest.mkdir(parents=True, exist_ok=True)
            shutil.move(item, dest / item.name)


def merge_parquet_files(parquet_files: list[Path], out_name: str) -> None:
    """Merges parquet files into a single parquet file."""
    if not parquet_files:
        logger.warning("No parquet files provided for merging")
        return

    # Filter out files that don't exist
    existing_files = [f for f in parquet_files if f.exists()]
    if not existing_files:
        logger.warning("No existing parquet files found for merging")
        return

    # Read and concatenate all parquet files
    dataframes = [pl.read_parquet(file) for file in existing_files]
    merged_df = pl.concat(dataframes, how="vertical")

    # Write merged the dataframe to the output location
    output_path = data_paths.library_dir / f"{out_name}.pqt"
    merged_df.write_parquet(output_path)

    logger.info(f"Merged {len(existing_files)} parquet files into {output_path}")


def finalise_from_representatives(
    library_name: str,
    representatives: list[Genome],
    records: list[RecordInfo],
    db_metadata: dict[str, Lineage],
    split_size: int = 0,
) -> None:
    if len(representatives) != len(records):
        raise ValueError(
            f"Representative count ({len(representatives)}) does not match record count ({len(records)})"
        )

    if split_size <= 0 or len(representatives) <= split_size:
        db_path = data_paths.build_dir / f"{library_name}.msh"
        paste_sketches(mash, [rep.sketch_path for rep in representatives], db_path)
        data_paths.finalise_library(library_name, db_path, records, db_metadata)
        move_files(data_paths.library_dir, library_name, data_paths.output_location)
        return

    n_shards = (len(representatives) + split_size - 1) // split_size
    logger.info(
        f"Splitting {library_name} into {n_shards} sub-libraries with max size {split_size}"
    )
    for shard_idx, start in enumerate(
        range(0, len(representatives), split_size), start=1
    ):
        end = min(start + split_size, len(representatives))
        shard_name = f"{library_name}.{shard_idx}"
        shard_path = data_paths.build_dir / f"{shard_name}.msh"
        reps_chunk = representatives[start:end]
        records_chunk = records[start:end]
        paste_sketches(mash, [rep.sketch_path for rep in reps_chunk], shard_path)
        data_paths.finalise_library(shard_name, shard_path, records_chunk, db_metadata)
        move_files(data_paths.library_dir, shard_name, data_paths.output_location)


@app.command()
def kleborate(
    version: Annotated[
        str, typer.Option("--version", "-v", help="Specify a kleborate version")
    ] = "v3.2.4",
    db_name: str = "Kleborate",
):
    # This prepares the standard Kleborate library
    download_kleborate(version, data_paths, mash=mash)
    logger.info(f"Downloaded Kleborate version {version}")

    # Now to make modifications
    modified_reps = [
        {
            "record_key": "Salmonella_enterica/GCF_003717755.fna.gz",
            "organism_name": "Salmonella enterica subsp. enterica serovar Typhi",
            "organism_code": "90370",
            "accession": "Salmonella_enterica/GCF_003717755",
        }
    ]

    lineages: dict[str, Lineage] = extract_kleborate_lineages(
        {rep["organism_name"]: rep["organism_code"] for rep in modified_reps},
        TaxonKit(taxonkit),
    )
    for rep in modified_reps:
        lineages[rep["organism_code"]].OrganismId = rep["organism_code"]
        lineages[rep["organism_code"]].OrganismName = rep["organism_name"]

    kleborate_df = pl.read_parquet(data_paths.library_dir / f"{db_name}.pqt")
    for rep in modified_reps:
        kleborate_df = kleborate_df.with_columns(
            organism_code=pl.when(pl.col("record_key") == rep["record_key"])
            .then(pl.lit(rep["organism_code"]))
            .otherwise(pl.col("organism_code"))
        )

    # Replace the kleborate parquet file
    if len(modified_reps):
        kleborate_df.write_parquet(data_paths.library_dir / f"{db_name}.pqt")

    # Update lineages file
    data_paths.update_lineage_file(lineages.values())

    # Copy to the final output directory.
    move_files(data_paths.library_dir, db_name, data_paths.output_location)
    shutil.copy(data_paths.lineages_path, data_paths.output_location / "lineages.pqt")


@app.command()
def flu(
    batch_size: int = 500,
    scaling_factor: float = 0.05,
    max_reps: int = 2000,
    db_name: str = "Influenza",
):
    qc_metrics = {}
    # Download/read genome metadata
    genomes, organism_ids = read_genbank_flu(data_paths, qc_metrics, FLU_SPECIES)
    species_mapping: dict[str, list[Genome]] = defaultdict(list)
    lineages: dict[str, Lineage] = TaxonKit(taxonkit).extract_lineages(organism_ids)

    logger.info(f"Building Flu library with species IDs {FLU_SPECIES}")

    for genome in genomes.values():
        taxid = str(genome.organism_taxid)
        if taxid not in lineages:
            logger.warning(f"No lineage for {taxid}")
            continue
        if taxid not in FLU_SPECIES:
            logger.debug(f"Skipping {taxid} due to not being flu virus")
            continue
        species_mapping[lineages[taxid].OrganismId].append(genome)

    # Download FASTA files
    logger.info("Downloading FASTA files.")
    sketched: dict[str, Genome] = download_and_sketch(
        genomes, qc_metrics, dwnld_threads
    )
    if len(sketched) != genomes:
        logger.error(
            f"Not all flu genomes were downloaded: {len(sketched)}/{len(genomes)}"
        )
        for species, sp_genomes in species_mapping.items():
            species_mapping[species] = [
                genome for genome in sp_genomes if genome.name in sketched
            ]
    reps: list[Genome] = []
    species_list: list[str] = []
    cluster_sizes: list[int] = []
    # db_to_reps: dict[str, list[Genome]] = defaultdict(list)  # Genomes in each DB
    # db_to_species: dict[str, list[str]] = defaultdict(list)  # Species in each DB
    # db_to_cluster_sizes: dict[str, list[int]] = defaultdict(list)  #
    for species, sp_genomes in species_mapping.items():
        logger.info(
            f"Selecting representatives for {species} from {len(sp_genomes)} genomes"
        )
        if species not in lineages:
            logger.warning(
                f"Error: Unexpectedly no lineage for {species}: {str([str(genome) for genome in sp_genomes])}"
            )
            continue
        if species not in FLU_SPECIES:
            logger.debug(
                f"Skipping {species} due to not being one of the flu viruses: {', '.join(FLU_SPECIES)}"
            )
            continue

        # Get all the relevant pairwise distances for creating the nr databases
        temp_lib: Path = data_paths.build_dir / f"{species}_full.msh"
        current_genomes: list[Genome] = sp_genomes
        threshold = config["thresholds"][db_name] * scaling_factor
        if len(sp_genomes) <= max_reps:
            nr_exemplars, nr_counts = (
                [genome.name for genome in sp_genomes],
                [1] * len(sp_genomes),
            )
        else:
            nr_exemplars, nr_counts = [], [1] * len(sp_genomes)
        while (
            len(current_genomes) > max_reps
            and threshold < config["thresholds"][db_name]
        ):
            logger.info(
                f"Building species library for {db_name} -> {species} with {len(current_genomes)} genomes and threshold {threshold}",
            )
            paste_sketches(
                mash, [genome.sketch_path for genome in current_genomes], temp_lib
            )
            nr_exemplars, nr_counts = extract_exemplars(
                batch_size, current_genomes, temp_lib, threshold, nr_counts
            )
            exemplar_order = {name: idx for idx, name in enumerate(nr_exemplars)}
            current_genomes = sorted(
                [genome for genome in current_genomes if genome.name in exemplar_order],
                key=lambda g: exemplar_order[g.name],
            )
            threshold = threshold * 2
            logger.debug(
                f"{len(current_genomes)} genomes remaining. Threshold adjusted to {threshold} - max threshold: {config['thresholds'][db_name]}"
            )
        else:
            # Create a mapping for efficient lookup
            reps.extend(
                [genome for genome in sp_genomes if genome.name in nr_exemplars]
            )
            species_list.append(species)
            cluster_sizes.extend(nr_counts)

    db_path = data_paths.build_dir / f"{db_name}.msh"
    logger.info(f"Building library for {db_name}, starting with {len(reps)} genomes")
    if len(reps) == 0:
        logger.exit(f"No genomes found for {db_name}")
    paste_sketches(mash, [rep.sketch_path for rep in reps], db_path)
    db_metadata = {str(species): lineages[str(species)] for species in species_list}
    data_paths.finalise_library(
        db_name,
        db_path,
        [
            RecordInfo(
                rep.name,
                rep.organism_taxid,
                rep.assembly_accession,
                cluster_sizes[idx],
            )
            for idx, rep in enumerate(reps)
        ],
        db_metadata,
    )
    move_files(data_paths.library_dir, db_name, data_paths.output_location)
    shutil.copy(data_paths.lineages_path, data_paths.output_location / "lineages.pqt")


def extract_exemplars(
    batch_size: int,
    current_genomes: list[Genome],
    mash_lib: Path,
    threshold: float,
    previous_counts: list[int] = None,
) -> tuple[list[str | int], list[int]]:
    if previous_counts is None:
        previous_counts = [1] * len(current_genomes)
    logger.debug(f"Calculating pairwise distances for {len(current_genomes)} genomes")

    linked_groups = build_linked_groups(batch_size, mash_lib, threshold)

    mash_lib.unlink()
    current_exemplars, current_counts = [], []
    # Select the nr representatives
    for group in linked_groups.groups:
        input_weights = [
            previous_counts[idx]
            for idx, genome in enumerate(current_genomes)
            if genome.assembly_accession in group.members
        ]
        group_exemplars, group_counts = get_exemplars(
            group.scores, threshold=threshold, input_weights=input_weights
        )
        current_exemplars.extend(group_exemplars)
        current_counts.extend(group_counts)
    return current_exemplars, current_counts


def build_linked_groups(batch_size: int, temp_lib: Path, threshold: float) -> Linker:
    linked_groups = Linker()
    scores_batch = []

    for genome1, genome2, distance in mash_distances(
        temp_lib, mash, threshold=threshold, threads=mash_threads
    ):
        scores_batch.append((genome1, genome2, distance))
        if len(scores_batch) >= batch_size:
            linked_groups.add_scores_batch(scores_batch)
            scores_batch = []

    if scores_batch:  # Add any remaining scores
        linked_groups.add_scores_batch(scores_batch)
    return linked_groups


@app.command()
def curated(db_name: str = "Curated", clean: bool = False) -> None:
    """Validates and imports the pre-computed curated libraries."""
    curated_libraries: dict[str, dict[str, str]] = {
        "Ecoli_Shigella": {
            "source_dir": "Curated_Ecoli_Shigella",
            "source_name": "Ecoli_Shigella",
        },
        "Vch_complex": {
            "source_dir": "Curated_Vch_complex",
            "source_name": "Vch_complex",
        },
        "Viridans": {
            "source_dir": "Curated_Viridans",
            "source_name": "Viridans",
        },
        # "Added": {
        #     "source_dir": data_paths.library_dir,
        #     "source_name": "Added",
        # },
    }
    added_reps = []

    for import_db, source_info in curated_libraries.items():
        if import_db == "Added":
            # Set up the "Added" library
            # Download added reps, sketch and add sketches to the curated library
            added_library: Path = data_paths.library_dir / f"{import_db}.msh"
            for rep in added_reps:
                subprocess.run(
                    f"wget -O {data_paths.fasta_dir / rep['file_name']} {rep['file_link']}",
                    shell=True,
                    check=True,
                )
                if (
                    not rep["file_name"].endswith(".gz")
                    and not (data_paths.fasta_dir / f"{rep['file_name']}.gz").exists()
                ):
                    subprocess.run(
                        f"gzip {data_paths.fasta_dir / rep['file_name']}",
                        shell=True,
                        check=True,
                    )
                rep_sketch = data_paths.sketch_dir / f"{rep['file_name']}.gz.msh"
                subprocess.run(
                    f"{mash} sketch -o {str(rep_sketch)} {data_paths.fasta_dir / rep['file_name']}.gz",
                    shell=True,
                    check=True,
                )
                if added_library.exists():
                    subprocess.run(
                        f"{mash} paste temp {str(added_library)} {str(rep_sketch)}",
                        shell=True,
                        check=True,
                    )
                    os.rename("temp.msh", f"{data_paths.library_dir}/{import_db}.msh")
                else:
                    shutil.copy(rep_sketch, added_library)

                # metadata
                pl.from_dicts(
                    [
                        dataclasses.asdict(
                            RecordInfo(
                                record["record_key"],
                                record["organism_code"],
                                record["accession"],
                                1,
                            )
                        )
                        for record in added_reps
                    ]
                ).write_parquet(data_paths.library_dir / f"{import_db}.pqt")
                # new lineages
                lineages = extract_kleborate_lineages(
                    {rep["organism_name"]: rep["organism_code"] for rep in added_reps},
                    TaxonKit(taxonkit),
                )
                data_paths.update_lineage_file(lineages.values())
                if clean:
                    os.unlink(f"{data_paths.sketch_dir / rep['file_name']}.gz.msh")
                    os.unlink(f"{data_paths.fasta_dir / rep['file_name']}.gz")

        else:
            import_library(
                source_dir=Path(source_info["source_dir"]),
                source_name=source_info["source_name"],
                db_name=import_db,
            )

    # Merge mash and parquet files.
    paste_sketches(
        mash,
        [
            data_paths.library_dir / f"{import_db}.msh"
            for import_db in curated_libraries.keys()
        ],
        data_paths.library_dir / f"{db_name}.msh",
    )
    merge_parquet_files(
        [
            data_paths.library_dir / f"{import_db}.pqt"
            for import_db in curated_libraries.keys()
        ],
        db_name,
    )

    move_files(data_paths.library_dir, db_name, data_paths.output_location)
    shutil.copy(data_paths.lineages_path, data_paths.output_location / "lineages.pqt")


def import_library(
    source_dir: Path = Path("Ecoli_Shigella_curated_lib"),
    source_name: str = "Ecoli_Shigella",
    db_name: str = "Ecoli_Shigella",
) -> None:
    """Validates and imports the pre-computed E. coli & Shigella library."""
    if not source_dir.exists():
        raise FileNotFoundError(f"Source directory {source_dir} does not exist.")

    @dataclasses.dataclass
    class DBFile:
        type: str
        source: Path
        destination: Path

    files: dict[str, DBFile] = {
        "mash": DBFile(
            type="mash",
            source=source_dir / f"{source_name}.msh.gz",
            destination=data_paths.library_dir / f"{db_name}.msh",
        ),
        "metadata": DBFile(
            type="pqt",
            source=source_dir / f"{source_name}.pqt.gz",
            destination=data_paths.library_dir / f"{db_name}.pqt",
        ),
    }

    for file in files.values():
        if file.source.exists():
            subprocess.run(["gunzip", file.source], check=True)
        unzipped_file = source_dir / file.source.name.replace(".gz", "")
        if not unzipped_file.exists():
            raise FileNotFoundError(
                f"{unzipped_file} does not exist. Either the gzipped file or uncompressed file must exist"
            )
        shutil.copy(unzipped_file, file.destination)
        subprocess.run(["gzip", unzipped_file], check=True)

    metadata = pl.read_parquet(files["metadata"].destination)
    organism_ids = set(metadata.select("organism_code").to_series().cast(str))
    lineages: dict[str, Lineage] = TaxonKit(taxonkit).extract_lineages(organism_ids)
    data_paths.update_lineage_file(lineages.values())

    logger.info(
        f"Imported {db_name} library from {source_dir} using name {source_name}"
    )


@app.command()
def ncbi(
    update_metadata: bool = False,
    clean: bool = False,
    qc_metrics_file: Path = Path("filtered_metrics.csv"),
    batch_size: int = 500,
    cluster_threshold: float = 0.0001,
    scaling_factor: float = 0.1,
    max_reps: int = 2000,
    skip_unclassified: bool = True,
) -> None:
    """Generates a library from NCBI genomes using the provided QC metrics. Note that currently it skips flu virus
    genomes, which are handled by a separate method."""
    qc_metrics: dict[str, dict[str, Metric]] = read_metrics(qc_metrics_file)

    logger.info("Starting NCBI data extraction")
    # Download/read genome metadata
    for kingdom in ["Bacteria", "Archaea", "Fungi", "Viruses"]:
        logger.info(f"Starting with {kingdom} genomes")
        genomes, organism_ids = read_kingdom_metadata(
            kingdom,
            data_paths,
            qc_metrics,
            update_metadata,
            skip_unclassified,
        )

        logger.info(f"Starting with {len(genomes)} genomes")

        # Update all the lineages
        logger.info("Updating lineages")
        lineages: dict[str, Lineage] = TaxonKit(taxonkit).extract_lineages(organism_ids)
        logger.info(f"Extracted {len(lineages)} lineages")
        data_paths.update_lineage_file(lineages.values())

        # Download FASTA files
        logger.info("Downloading FASTA files.")
        genomes_after_qc: dict[str, Genome] = download_and_sketch(
            genomes, qc_metrics, dwnld_threads
        )

        # Group genomes by species
        species_mapping: dict[str, list[Genome]] = defaultdict(list)
        for genome in genomes_after_qc.values():
            taxid = str(genome.organism_taxid)
            if taxid not in lineages:
                logger.warning(f"No lineage for {taxid}")
                continue
            if taxid in FLU_SPECIES:  # Alpha- & Betainfluenzavirus
                logger.debug(
                    f"Skipping {taxid} due to being one of the flu virus species: {', '.join(FLU_SPECIES)}"
                )
                continue
            species_mapping[lineages[taxid].OrganismId].append(genome)

        # Extract representative genomes per species
        kingdom_reps: list[Genome] = []
        kingdom_species: list[str] = []
        cluster_sizes: list[int] = []
        for species, sp_genomes in species_mapping.items():
            if species not in lineages:
                logger.warning(
                    f"Error: Unexpectedly no lineage for {species}: {str([str(genome) for genome in sp_genomes])}"
                )
                continue
            if species in FLU_SPECIES:
                logger.debug(
                    f"Skipping {species} due to being one of the flu virus species: {', '.join(FLU_SPECIES)}"
                )
                continue
            if (db_name := lineages[species].SuperkingdomName) is None or db_name == "":
                logger.warning(
                    f"Error: No superkingdom for {species} ({lineages[species].SpeciesName}): {str(lineages[species])}"
                )
                continue
            # if db_name != kingdom:
            #     raise ValueError(f"Unexpected superkingdom for {species} ({lineages[species].SpeciesName}): {str(lineages[species])}")
            temp_lib = data_paths.build_dir / f"{species}_full.msh"
            logger.info(
                f"Building species library for {db_name} -> {species} with {len(sp_genomes)} genomes"
            )
            if len(sp_genomes) == 1:
                # shutil.copy(sp_genomes[0].sketch_path, temp_lib)
                nr_exemplars, nr_counts = [sp_genomes[0].assembly_accession], [1]
                exemplar_genomes = sp_genomes
            else:
                # Run the initial clustering pass once (was duplicated across branches)
                paste_sketches(
                    mash, [genome.sketch_path for genome in sp_genomes], temp_lib
                )
                nr_exemplars, nr_counts = extract_exemplars(
                    batch_size,
                    sp_genomes,
                    temp_lib,
                    cluster_threshold,
                )
                exemplar_order = {
                    accession: idx for idx, accession in enumerate(nr_exemplars)
                }
                exemplar_genomes = sorted(
                    [
                        genome
                        for genome in sp_genomes
                        if genome.assembly_accession in exemplar_order
                    ],
                    key=lambda g: exemplar_order[g.assembly_accession],
                )

                # Only do iterative refinement when there are too many reps
                if len(sp_genomes) > max_reps:
                    threshold = min(
                        config["thresholds"][db_name] * scaling_factor,
                        config["thresholds"][db_name],
                    )
                    while (
                        len(exemplar_genomes) > max_reps
                        and threshold <= config["thresholds"][db_name]
                    ):
                        logger.info(
                            f"Building species library for {db_name} -> {species} with {len(exemplar_genomes)} genomes and threshold {threshold}",
                        )
                        paste_sketches(
                            mash,
                            [genome.sketch_path for genome in exemplar_genomes],
                            temp_lib,
                        )
                        nr_exemplars, nr_counts = extract_exemplars(
                            batch_size, exemplar_genomes, temp_lib, threshold, nr_counts
                        )
                        exemplar_order = {
                            accession: idx for idx, accession in enumerate(nr_exemplars)
                        }
                        exemplar_genomes = sorted(
                            [
                                genome
                                for genome in exemplar_genomes
                                if genome.assembly_accession in exemplar_order
                            ],
                            key=lambda g: exemplar_order[g.assembly_accession],
                        )
                        # Exit the loop if the maximum threshold has been reached
                        if threshold == config["thresholds"][db_name]:
                            break
                        threshold = min(threshold * 2, config["thresholds"][db_name])
                        logger.debug(
                            f"{len(exemplar_genomes)} genomes remaining. Threshold adjusted to {threshold} - max threshold: {config['thresholds'][db_name]}"
                        )
            # Create a mapping for efficient lookup
            kingdom_reps.extend(exemplar_genomes)
            kingdom_species.append(species)
            cluster_sizes.extend(nr_counts)
            if not nr_exemplars:
                logger.warning(f"No exemplars found for {db_name} -> {species}")
                continue
            logger.info(
                f"Selected {len(nr_exemplars)} exemplars for {db_name} -> {species}"
            )

        logger.info(f"Building kingdom library for {kingdom}")
        db_metadata = {
            str(species): lineages[str(species)] for species in set(kingdom_species)
        }
        records = [
            RecordInfo(
                rep.assembly_accession,
                rep.organism_taxid,
                rep.assembly_accession,
                cluster_sizes[idx],
            )
            for idx, rep in enumerate(kingdom_reps)
        ]
        finalise_from_representatives(
            kingdom,
            kingdom_reps,
            records,
            db_metadata,
            split_size=max_sub_library_size,
        )
    logger.info("Finalising lineages file")
    shutil.copy(data_paths.lineages_path, data_paths.output_location / "lineages.pqt")

    if clean:
        shutil.rmtree(data_paths.fasta_dir)
        shutil.rmtree(data_paths.sketch_dir)
        shutil.rmtree(data_paths.build_dir)
        shutil.rmtree(data_paths.library_dir)


def process_genome(
    genome: Genome,
    qc_metrics: dict[str, dict[str, Metric]],
    mash_path: str,
) -> Optional[tuple[str, Genome]]:
    """Process a single genome: check if the sketch exists, download FASTA if needed, and sketch it."""
    key = genome.name

    # Check if the sketch already exists
    if genome.sketch_path.exists() and genome.sketch_path.stat().st_size != 0:
        if genome.fasta_path.exists():
            genome.fasta_path.unlink()
        return key, genome

    # Check QC
    if str(genome.species_taxid) in qc_metrics and not genome.passes_qc(
        qc_metrics[str(genome.species_taxid)]
    ):
        logger.debug(f"Skipping {genome.name} due to QC failure")
        return None

    # Download FASTA
    try:
        genome.download_fasta()
    except RetryError as e:
        logger.warning(f"Error downloading FASTA for {genome.fasta_path}: {str(e)}.")
        try:
            genome.fasta_path.unlink()
        except FileNotFoundError:
            pass
        return None

    # Sketch the genome
    if not genome.sketch_path.exists() or genome.sketch_path.stat().st_size == 0:
        sketch_file(mash_path, genome.fasta_path, genome.sketch_path)

    # Clean up FASTA
    genome.fasta_path.unlink()

    return key, genome


def download_and_sketch(
    genomes: dict[str, Genome],
    qc_metrics: dict[str, dict[str, Metric]],
    n_threads: int = 4,
) -> dict[str, Genome]:
    """Download and sketch genomes in parallel using n threads."""
    sketched_genomes: dict[str, Genome] = {}

    print(f"Sketching genomes in parallel using {n_threads} threads...")

    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        # Submit all tasks
        future_to_key = {
            executor.submit(process_genome, genome, qc_metrics, mash): key
            for key, genome in genomes.items()
        }

        # Process completed tasks
        for future in as_completed(future_to_key):
            try:
                result = future.result()
                if result is not None:
                    key, genome = result
                    sketched_genomes[key] = genome
            except Exception as e:
                original_key = future_to_key[future]
                logger.error(f"Error processing genome {original_key}: {str(e)}")

    return sketched_genomes


@app.command()
def full(
    kleborate_version: Annotated[
        str,
        typer.Argument(help="Specify a kleborate version"),
    ],
    clean: Annotated[
        bool,
        typer.Option(
            "--clean", "-c", help="Clean the output directory before running."
        ),
    ] = False,
    qc_metrics_file: Annotated[
        Path,
        typer.Option("--qc-metrics", "-q", help="Location of the QC metrics file."),
    ] = Path("filtered_metrics.csv"),
    ncbi_batch_size: Annotated[
        int,
        typer.Option(
            "--batch-size",
            "-b",
            help="Batch size for adding scores. Increasing this value may improve speed but may also use more memory.",
        ),
    ] = 500,
    ncbi_cluster_threshold: Annotated[
        float,
        typer.Option(
            "--ncbi-cluster-threshold",
            "-c",
            help="Clustering threshold for NCBI genomes. Smaller values will result in fewer genomes being used.",
        ),
    ] = 0.0001,
    ncbi_scaling_factor: Annotated[
        float,
        typer.Option(
            "--ncbi-scaling-factor",
            "-s",
            help="Scaling factor for the NCBI genomes sketches.",
        ),
    ] = 0.1,
    ncbi_max_reps: Annotated[
        int,
        typer.Option(
            "--ncbi-max-reps",
            "-m",
            help="Maximum number of NCBI genomes to use for building the library",
        ),
    ] = 2000,
    ncbi_skip_unclassified: Annotated[
        bool,
        typer.Option(
            "--ncbi-skip-unclassified",
            "-u",
            help="Skip unclassified genomes from NCBI",
        ),
    ] = True,
    curated_name: Annotated[
        str, typer.Option(help="Name of the curated library.")
    ] = "Curated",
    flu_batch_size: Annotated[
        int,
        typer.Option("--flu-batch-size", "-B", help="Batch size for adding scores"),
    ] = 500,
    flu_scaling_factor: Annotated[
        float,
        typer.Option(
            "--flu-scaling-factor",
            "-S",
            help="Scaling factor for the NCBI genomes sketches. Smaller takes longer but uses less RAM.",
        ),
    ] = 0.05,
    flu_max_reps: Annotated[
        int,
        typer.Option(
            "--flu-max-reps",
            "-M",
            help="Maximum number of reps to use for building the flu library",
        ),
    ] = 2000,
):
    print(f"{str(data_paths)}", file=sys.stderr)
    if data_paths.lineages_path.exists():
        data_paths.lineages_path.unlink()

    if clean:
        shutil.rmtree(data_paths.library_dir)
        os.makedirs(data_paths.library_dir)
        shutil.rmtree(data_paths.output_location)
        os.makedirs(data_paths.output_location)

    # Individual builds.
    kleborate(version=kleborate_version)
    curated(curated_name)
    ncbi(
        update_metadata=False,
        clean=False,
        qc_metrics_file=qc_metrics_file,
        batch_size=ncbi_batch_size,
        scaling_factor=ncbi_scaling_factor,
        max_reps=ncbi_max_reps,
        skip_unclassified=ncbi_skip_unclassified,
        cluster_threshold=ncbi_cluster_threshold,
    )
    flu(
        batch_size=flu_batch_size,
        scaling_factor=flu_scaling_factor,
        max_reps=flu_max_reps,
    )


@app.callback()
def main(
    config_file: Annotated[
        Path,
        typer.Option(
            "--config-file",
            "-c",
            help="Path to the TOML configuration file",
            dir_okay=False,
            exists=True,
            file_okay=True,
        ),
    ] = Path("config.toml"),
    mash_path: Annotated[
        str,
        typer.Option(
            "--mash-path", "-m", help="Path to the mash executable (if not on $PATH)"
        ),
    ] = "mash",
    taxonkit_path: Annotated[
        str,
        typer.Option(
            "--taxonkit-path",
            "-t",
            help="Path to the taxonkit executable (if not on $PATH)",
        ),
    ] = "taxonkit",
    build_dir: Annotated[
        Path,
        typer.Option(
            "--build-space-dir",
            "-b",
            help="Directory for the build process",
            dir_okay=True,
            file_okay=False,
        ),
    ] = Path("build_space"),
    out_dir: Annotated[
        Path,
        typer.Option(
            "--output-dir",
            "-o",
            help="Output directory for the final libraries",
            dir_okay=True,
            file_okay=False,
        ),
    ] = "final",
    log_level: Annotated[
        str,
        typer.Option(
            "--log-level",
            "-L",
            help="Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)",
            case_sensitive=False,
        ),
    ] = "WARNING",
    download_threads: Annotated[
        int,
        typer.Option(
            "--download-threads",
            "-d",
            help="Number of threads for downloading genomes",
        ),
    ] = 2,
    search_threads: Annotated[
        int,
        typer.Option(
            "--search-threads",
            "-T",
            help="Number of threads for sketching genomes",
        ),
    ] = os.cpu_count(),
    max_sub_library_size_override: Annotated[
        Optional[int],
        typer.Option(
            "--max-sub-library-size",
            help="Maximum representatives per output mash sub-library. 0 disables splitting.",
            min=0,
        ),
    ] = 25000,
):
    global mash
    mash = mash_path
    global taxonkit
    taxonkit = taxonkit_path
    global data_paths
    data_paths = BuildSpace(build_dir, out_dir)
    global dwnld_threads
    dwnld_threads = download_threads
    global mash_threads
    mash_threads = search_threads
    global config
    config = toml.load(config_file)
    global max_sub_library_size
    configured_size = int(config.get("parameters", {}).get("max_sub_library_size", 25000))
    max_sub_library_size = (
        max_sub_library_size_override
        if max_sub_library_size_override is not None
        else configured_size
    )
    global logger
    # Set up logging
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")
    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger.info(f"Initializing with log level: {log_level}")
    logger.info(f"Max sub-library size: {max_sub_library_size}")
    logger.info(f"{str(data_paths)}")


if __name__ == "__main__":
    app()
