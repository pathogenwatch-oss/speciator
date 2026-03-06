from __future__ import annotations

import dataclasses
import logging
import os
import subprocess

from speciator.build_space import BuildSpace, RecordInfo
from speciator.lineages import _UNKNOWN_CODE, _UNKNOWN_SPECIES_CODE, Lineage, TaxonKit

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class KleborateRecordInfo:
    """Information for a record in the Kleborate mash database."""

    record_key: str
    species_name: str
    accession: str


def map_kleborate_records(
    records: list[KleborateRecordInfo], species_name2code: dict[str, str]
) -> list[RecordInfo]:
    """Map a set of Kleborate records to a set of KleborateRecordInfo objects."""
    return [
        RecordInfo(
            record_key=record.record_key,
            organism_code=species_name2code[record.species_name],
            accession=record.accession,
            represents=1,
        )
        for record in records
    ]


def download_kleborate(
    version: str,
    build_space: BuildSpace,
    mash: str = "mash",
    taxonkit: TaxonKit = TaxonKit(
        "taxonkit",
    ),
):
    # Download the kleborate file
    kleborate_filename = "species_mash_sketches.msh"
    kleborate_mash = build_space.build_dir / kleborate_filename
    if not os.path.exists(kleborate_mash):
        kleborate_url = f"https://github.com/klebgenomics/Kleborate/raw/{version}/kleborate/modules/enterobacterales__species/data/{kleborate_filename}"
        subprocess.run(
            f"wget -O {kleborate_mash} {kleborate_url}", shell=True, check=True
        )

    if not os.path.exists(kleborate_mash):
        logger.error(f"Failed to download Kleborate mash database to {kleborate_mash}")

    mash_command = f"{mash} info {kleborate_mash}"
    try:
        kleborate_records = subprocess.run(
            mash_command,
            shell=True,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing mash info for command '{mash_command}'")
        logger.error(str(e))
        raise e

    # Example record:
    #  1000 4987107 Arsenophonus_nasoniae/GCF_004768525.fna.gz [18 seqs] NZ_CP038613.1 Arsenophonus nasoniae strain FIN chromosome, complete genome [...]
    species_names: set[str] = set()
    kleborate_recs: list[KleborateRecordInfo] = []

    for record in kleborate_records.stdout.splitlines():
        if not record.startswith("  1000"):
            continue
        hashes, length, name = record.strip().split()[0:3]
        species_raw, accession = name.split("/")
        accession = accession.replace(".fna.gz", "")
        species = (
            species_raw.replace("_", " ").replace("subsp", "subsp.").replace("NEW", "")
        )
        kleborate_recs.append(KleborateRecordInfo(name, species, accession))
        species_names.add(species)

    species_name2code: dict[str, str] = taxonkit.convert_names_to_taxid(species_names)

    kleborate_lineages = extract_kleborate_lineages(species_name2code, taxonkit)
    recs: list[RecordInfo] = map_kleborate_records(kleborate_recs, species_name2code)
    build_space.finalise_library(
        "Kleborate",
        kleborate_mash,
        recs,
        kleborate_lineages,
    )


def extract_kleborate_lineages(
    species_name2code: dict[str, str], taxonkit: TaxonKit
) -> dict[str, Lineage]:
    """Extract lineages for a list of taxon IDs."""
    lineages: dict[str, Lineage] = dict()
    filtered_taxids: set[str] = set()

    for taxon_name, taxon_id in species_name2code.items():
        if taxon_id == _UNKNOWN_CODE or taxon_id == _UNKNOWN_SPECIES_CODE:
            # Species not in NCBI taxonomy database
            lineages[taxon_id] = Lineage.build_unknown_species()
        else:
            filtered_taxids.add(taxon_id)

    # Test for missing lineages and create them if needed
    lineages: dict[str, Lineage] = lineages | taxonkit.extract_lineages(
        filtered_taxids | {"573"}, add_unknown=True
    )
    for taxon_name, taxon_id in species_name2code.items():
        if taxon_id not in lineages:
            if taxon_name.startswith("Klebsiella"):
                # Kleborate novel Klebsiella
                logger.info(
                    f"Creating Klebsiella lineage for novel Klebsiella species {taxon_name} ({taxon_id})"
                )
                lineages[taxon_id] = lineages[taxon_id] = Lineage.replace_species(
                    lineages["573"],
                    Lineage.species_only(taxon_name, taxon_id),
                    taxon_id,
                )
            else:
                raise ValueError(f"{taxon_name}/{taxon_id} is missing a record")
    return lineages
