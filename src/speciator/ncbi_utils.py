from __future__ import annotations

import csv
import dataclasses
import gzip
import logging
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Callable, Iterable, Iterator

import numpy as np
from sklearn.cluster import AgglomerativeClustering
from tenacity import retry, stop_after_attempt, wait_exponential

from speciator.build_space import BuildSpace
from speciator.mash_functions import self_search
from speciator.metrics import Metric

logger = logging.getLogger(__name__)

REFSEQ_URLS = {
    "Bacteria": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
    "Archaea": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt",
    "Fungi": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
    "Viruses": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt",
}

GENBANK_URLS = {
    "Bacteria": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt",
    "Archaea": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt",
    "Fungi": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt",
    "Viruses": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt",
}


@dataclasses.dataclass
class Genome:
    name: str
    assembly_accession: str
    shared_accession: str
    url: str
    assembly_status: str
    gc_percent: float
    species_taxid: int  # Ensures it's a valid taxid
    size: int
    contigs: int
    data_repository: BuildSpace
    organism_name: dataclasses.InitVar[str]
    taxid: dataclasses.InitVar[str]

    # a bit ugly maybe ...
    def __post_init__(self, organism_name: str, taxid: str):
        self.fasta_path: Path = (
            self.data_repository.fasta_dir / f"{self.assembly_accession}.fna.gz"
        )

        self.sketch_path: Path = (
            self.data_repository.sketch_dir / f"{self.assembly_accession}.msh"
        )
        if (
            organism_name.startswith("Salmonella") and "Typhi" in organism_name
        ) or taxid == "90370":
            self.organism_taxid = "90370"
        elif (
            taxid == "2697049"
            or organism_name == "Severe acute respiratory syndrome coronavirus 2"
        ):
            self.organism_taxid = "2697049"
        else:
            self.organism_taxid = self.species_taxid

    @retry(wait=wait_exponential(min=2, max=10), stop=stop_after_attempt(2))
    def download_fasta(self) -> Path:
        """Download the genome from NCBI. Names the downloaded file using self.fasta_path."""
        if not self.fasta_path.exists():
            subprocess.run(
                f"wget -c --retry-connrefused --tries=3 --timeout=5 -O {self.fasta_path} {self.url}/{self.url.split('/')[-1]}_genomic.fna.gz",
                shell=True,
                check=True,
            )
        if self.fasta_path.stat().st_size == 0:
            self.fasta_path.unlink()
            raise FileNotFoundError(f"Download of FASTA file failed: {self.url} ")
        return self.fasta_path

    def passes_qc(self, metrics: dict[str, Metric]) -> bool:
        fasta_metrics: dict[str, float] = self.fasta_stats()
        return (
            metrics["GC_Content"].lower_bounds
            <= float(self.gc_percent)
            <= metrics["GC_Content"].upper_bounds
            and metrics["Genome_Size"].lower_bounds
            <= int(self.size)
            <= metrics["Genome_Size"].upper_bounds
            and metrics["N50"].lower_bounds
            <= fasta_metrics["N50"]
            <= metrics["N50"].upper_bounds
            and metrics["no_of_contigs"].lower_bounds
            <= fasta_metrics["Contig_Count"]
            <= metrics["no_of_contigs"].upper_bounds
        )

    def fasta_stats(self) -> dict[str, float]:
        if not self.fasta_path.exists():
            raise FileNotFoundError(
                f"No FASTA file found at {self.fasta_path}. Try downloading it first using the `download_fasta()` function."
            )
        with gzip.open(self.fasta_path, "rt") as f:
            first = f.readline().strip()
            if not first.startswith(">"):
                raise ValueError(f"Invalid FASTA format: {self.fasta_path}")
            count = 1
            sizes = []
            current = ""
            for line in f:
                if line.startswith(">"):
                    sizes.append(len(current))
                    current = ""
                    count += 1
                else:
                    current += line.strip()
            sizes.append(len(current))
            count += 1
        sizes.sort(reverse=True)
        # Calculate N50
        total_length = sum(sizes)
        cumulative_length = 0
        n50 = 0
        for size in sizes:
            cumulative_length += size
            if cumulative_length >= total_length / 2:
                n50 = size
                break

        return {
            "Contig_Count": count,
            "Total_Length": total_length,
            "N50": n50,
            "Largest_Contig": sizes[0],
            "Smallest_Contig": sizes[-1],
        }


def fetch_assembly_summary(url: str, summary_file: Path, force: bool = False) -> Path:
    if force:
        summary_file.unlink()
    else:
        if summary_file.exists():
            return summary_file
    logger.debug(f"Writing {url} to {str(summary_file)}")
    subprocess.run(
        f"wget -c --retry-connrefused --tries=4 --timeout=5 -O {str(summary_file)} {url}",
        shell=True,
        check=True,
    )
    return summary_file


def read_genomes(
    summary_file: Path,
    data_repository: BuildSpace,
    condition: Callable[[dict[str, str]], bool] = lambda genome: True,
) -> tuple[dict[str, Genome], set[str]]:
    selected_genomes: dict[str, Genome] = dict()
    organism_ids: set[str] = set()

    with open(summary_file, "r") as f:
        f.readline()  # Skip initial comment line
        reader: csv.DictReader = csv.DictReader(f, delimiter="\t")

        for genome in reader:
            if condition(genome):
                name = genome["#assembly_accession"]
                if (
                    name not in selected_genomes
                    and genome["gbrs_paired_asm"] not in selected_genomes
                    and genome["#assembly_accession"] not in selected_genomes
                ):
                    selected_genomes[name] = Genome(
                        name=name,
                        assembly_accession=genome["#assembly_accession"],
                        shared_accession=genome["gbrs_paired_asm"],
                        url=genome["ftp_path"],
                        assembly_status=genome["assembly_level"],
                        gc_percent=genome["gc_percent"],
                        species_taxid=genome["species_taxid"],
                        size=genome["genome_size"],
                        contigs=genome["contig_count"],
                        data_repository=data_repository,
                        organism_name=genome["organism_name"],
                        taxid=genome["taxid"],
                    )
                    organism_ids.add(selected_genomes[name].organism_taxid)
                # logger.debug(f"Currently keeping {len(genomes)} genomes")
    return selected_genomes, organism_ids


def read_genbank_flu(
    build_space: BuildSpace,
    qc_metrics: dict[str, dict[str, Metric]],
    species: list[str],
) -> tuple[dict[str, Genome], set[str]]:
    def select(query_genome: dict[str, str]) -> bool:
        if (
            query_genome["assembly_level"] != "Complete Genome"
            or query_genome["ftp_path"] == "na"
            or query_genome["species_taxid"] not in species
        ):
            return False
        metrics = qc_metrics.get(query_genome["species_taxid"])
        if metrics is None:
            # logger.debug(f"No QC metrics found for {query_genome['organism_name']}")
            return float(query_genome["genome_size"]) > 1000
        else:
            return (
                "GC_Content" in metrics
                and metrics["GC_Content"].lower_bounds
                <= float(query_genome["gc_percent"])
                <= metrics["GC_Content"].upper_bounds
            ) and (
                "Genome_Size" in metrics
                and metrics["Genome_Size"].lower_bounds
                <= float(query_genome["genome_size"])
                <= metrics["Genome_Size"].upper_bounds
            )

    genbank_summary_file = fetch_assembly_summary(
        "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt",
        build_space.build_dir / "genbank_viral.txt",
    )
    return read_genomes(
        genbank_summary_file,
        build_space,
        select,
    )


def read_kingdom_metadata(
    kingdom: str,
    data_repository: BuildSpace,
    qc_metrics: dict[str, dict[str, Metric]],
    force: bool = False,
    skip_unclassified: bool = False,
) -> tuple[dict[str, Genome], set[str]]:
    if kingdom not in REFSEQ_URLS:
        logger.error(
            f"Invalid kingdom: {kingdom}.\nShould be one of: {', '.join(REFSEQ_URLS.keys())}."
        )
    refseq_url = REFSEQ_URLS[kingdom]
    genbank_url = GENBANK_URLS[kingdom]

    def refseq_filter(query_genome: dict[str, str]) -> bool:
        if (
            query_genome["assembly_level"] == "Transcriptome"
            or query_genome["ftp_path"] == "na"
        ):
            return False
        elif skip_unclassified and (
            " sp." in query_genome["organism_name"]
            or " bacterium" in query_genome["organism_name"]
        ):
            return False

        metrics = qc_metrics.get(query_genome["species_taxid"])
        if metrics is None:
            return float(query_genome["genome_size"]) > 1000.0

        gc_ok = (
            "GC_Content" in metrics
            and metrics["GC_Content"].lower_bounds
            <= float(query_genome["gc_percent"])
            <= metrics["GC_Content"].upper_bounds
        )
        size_ok = (
            "Genome_Size" in metrics
            and metrics["Genome_Size"].lower_bounds
            <= float(query_genome["genome_size"])
            <= metrics["Genome_Size"].upper_bounds
        )
        return gc_ok and size_ok

    genomes, organisms = from_ncbi(
        "refseq", {kingdom: refseq_url}, refseq_filter, data_repository, force
    )

    logger.info(
        f"Read {len(genomes)} genomes and {len(organisms)} species from RefSeq {kingdom}"
    )

    refseq_genome_names = set(genomes.keys())

    def genbank_filter(query_genome: dict[str, str]) -> bool:
        return (
            # Only use complete genomes from GenBank
            query_genome["assembly_level"] == "Complete Genome"
            # query_genome["assembly_level"] != "Transcriptome"
            and refseq_filter(query_genome)
            and query_genome["gbrs_paired_asm"] not in refseq_genome_names
        )

    genbank_genomes, genbank_organisms = from_ncbi(
        "genbank", {kingdom: genbank_url}, genbank_filter, data_repository, force
    )

    logger.info(
        f"Read {len(genbank_genomes)} genomes and {len(genbank_organisms)} species from GenBank {kingdom}"
    )

    genomes.update(genbank_genomes)
    organisms.update(genbank_organisms)

    logger.info(
        f"Read {len(genomes)} genomes and {len(organisms)} species from RefSeq and GenBank for kingdom {kingdom}"
    )

    return genomes, organisms

def read_kingdom_incremental(
        kingdom: str,
        data_repository: BuildSpace,
        qc_metrics: dict[str, dict[str, Metric]],
        force: bool = False,
        skip_unclassified: bool = False,
        extraction_level: str = "contig",
        all_complete: bool = True
):
    extracted: dict[str, Genome] = dict()
    extracted_species: set[str] = set()
    if kingdom not in REFSEQ_URLS:
        raise ValueError(f"Kingdom {kingdom} not supported. Supported kingdoms are: {', '.join(REFSEQ_URLS.keys())}")

    logger.info(f"Beginning incremental species extraction for {kingdom}")
    extraction_level = extraction_level.lower()
    extraction_order = ["complete", "chromosome", "scaffold", "contig"]
    if extraction_level not in extraction_order:
        raise ValueError(
            "extraction_level must be one of: "
            + ", ".join(extraction_order)
        )
    max_level_index = extraction_order.index(extraction_level)

    refseq_url = REFSEQ_URLS[kingdom]
    genbank_url = GENBANK_URLS[kingdom]

    def base_filter(query_genome):
        if  skip_unclassified and (
            " sp." in query_genome["organism_name"]
            or " bacterium" in query_genome["organism_name"]
        ):
            return False
        if query_genome["ftp_path"] == "na":
            return False
        metrics = qc_metrics.get(query_genome["species_taxid"])
        if metrics is None:
            return float(query_genome["genome_size"]) > 1000.0
        gc_ok = (
                "GC_Content" in metrics
                and metrics["GC_Content"].lower_bounds
                <= float(query_genome["gc_percent"])
                <= metrics["GC_Content"].upper_bounds
        )
        size_ok = (
                "Genome_Size" in metrics
                and metrics["Genome_Size"].lower_bounds
                <= float(query_genome["genome_size"])
                <= metrics["Genome_Size"].upper_bounds
        )
        return gc_ok and size_ok

    def duplicate_filter(query_genome):
        return query_genome["gbrs_paired_asm"] not in extracted and base_filter(query_genome)

    def additive_filter(query_genome):
        return query_genome["taxid"] not in extracted_species and query_genome["gbrs_paired_asm"] not in extracted and base_filter(query_genome)

    def refseq_complete(query_genome):
        return query_genome["assembly_level"] == "Complete Genome" and base_filter(query_genome)

    complete_genome_filter = duplicate_filter if all_complete else additive_filter

    def genbank_complete(query_genome):
        return query_genome["assembly_level"] == "Complete Genome" and complete_genome_filter(query_genome)

    def refseq_chromosome(query_genome):
        return query_genome["assembly_level"] == "Chromosome" and additive_filter(query_genome)

    def refseq_scaffold(query_genome):
        return query_genome["assembly_level"] == "Scaffold" and additive_filter(query_genome)

    def refseq_contig(query_genome):
        return query_genome["assembly_level"] == "Contig" and additive_filter(query_genome)

    # for name, url in REFSEQ_URLS.items():
    metadata_file = data_repository.build_dir / f"refseq_{kingdom}.txt"
    logger.debug(f"Searching for {metadata_file} from {refseq_url}")

    refseq_summary_file = fetch_assembly_summary(refseq_url, metadata_file, force=force)

    metadata_file = data_repository.build_dir / f"genbank_{kingdom}.txt"
    logger.debug(f"Searching for {metadata_file} from {genbank_url}")

    genbank_summary_file = fetch_assembly_summary(genbank_url, metadata_file, force=force)

    # complete genomes
    if max_level_index >= extraction_order.index("complete"):
        refseq_complete_genomes, refseq_complete_orgs = read_genomes(
            refseq_summary_file,
            data_repository,
            refseq_complete,
        )

        logger.info(
            f"Read {len(refseq_complete_genomes)} genomes and {len(refseq_complete_orgs)} species for {kingdom} with complete genomes in Refseq"
        )

        extracted.update(refseq_complete_genomes)
        extracted_species.update(refseq_complete_orgs)

        extra_genomes, extra_orgs = read_genomes(
            genbank_summary_file,
            data_repository,
            genbank_complete,
        )

        logger.info(
            f"Read {len(extra_genomes)} genomes and {len(extra_orgs)} species for {kingdom} with complete genomes in Genbank"
        )

        extracted.update(extra_genomes)
        extracted_species.update(extra_orgs)

    # Refseq chromosomes
    if max_level_index >= extraction_order.index("chromosome"):
        extra_genomes, extra_orgs = read_genomes(
            refseq_summary_file,
            data_repository,
            refseq_chromosome
        )

        logger.info(
            f"Read {len(extra_genomes)} genomes and {len(extra_orgs)} species for {kingdom} ... chromosomes ... RefSeq"
        )

        extracted.update(extra_genomes)
        extracted_species.update(extra_orgs)

    # Refseq scaffolds
    if max_level_index >= extraction_order.index("scaffold"):
        extra_genomes, extra_orgs = read_genomes(
            refseq_summary_file,
            data_repository,
            refseq_scaffold
        )

        logger.info(
            f"Read {len(extra_genomes)} genomes and {len(extra_orgs)} species for {kingdom} ... scaffolds ... RefSeq"
        )

        extracted.update(extra_genomes)
        extracted_species.update(extra_orgs)

    # Refseq contigs
    if max_level_index >= extraction_order.index("contig"):
        extra_genomes, extra_orgs = read_genomes(
            refseq_summary_file,
            data_repository,
            refseq_contig
        )

        logger.info(
            f"Read {len(extra_genomes)} genomes and {len(extra_orgs)} species for {kingdom} ... contigs ... RefSeq"
        )

        extracted.update(extra_genomes)
        extracted_species.update(extra_orgs)

    return extracted, extracted_species


def from_ncbi(
    tag: str,
    urls: dict[str, str],
    select_genome: Callable[[dict[str, str]], bool],
    data_repository: BuildSpace,
    force: bool = False,
) -> tuple[dict[str, Genome], set[str]]:
    """Used for downloading RefSeq/GenBank genomes. Urls is a list of NCBI RefSeq/GenBank assembly_summary.txt URLs."""

    extracted: dict[str, Genome] = dict()
    extracted_species: set[str] = set()
    for name, url in urls.items():
        metadata_file = data_repository.build_dir / f"{tag}_{name}.txt"
        logger.debug(f"Searching for {metadata_file} from {url}")
        refseq_summary_file = fetch_assembly_summary(url, metadata_file, force=force)
        subset_refseq_genomes, subset_refseq_organisms = read_genomes(
            refseq_summary_file,
            data_repository,
            select_genome,
        )
        extracted.update(subset_refseq_genomes)
        extracted_species.update(subset_refseq_organisms)
        logger.info(
            f"Read {len(subset_refseq_genomes)} genomes and {len(subset_refseq_organisms)} species for {name}"
        )

    return extracted, extracted_species


class BiMap:
    def __init__(self):
        self.name_to_index: dict[str, int] = {}
        self.index_to_name: list[str] = []

    def add(self, name: str) -> int:
        if name not in self.name_to_index:
            index = len(self.index_to_name)
            self.name_to_index[name] = index
            self.index_to_name.append(name)
        return self.name_to_index[name]

    def get_index(self, name: str) -> int:
        return self.name_to_index[name]

    def get_name(self, index: int) -> str:
        return self.index_to_name[index]

    def __len__(self) -> int:
        return len(self.index_to_name)

    @staticmethod
    def map_score_dict(
        scores: dict[tuple[str | int, str | int], float],
    ) -> tuple[dict[tuple[str | int, str | int], float], BiMap]:
        bimap = BiMap()
        mapped_scores: dict[tuple[str | int, str | int], float] = {}
        for id_pair, score in scores.items():
            bimap.add(id_pair[0])
            bimap.add(id_pair[1])
            mapped_scores[
                (bimap.get_index(id_pair[0]), bimap.get_index(id_pair[1]))
            ] = score
        return mapped_scores, bimap

    def invert_list(self, mapped_list: list[int]) -> list[str | int]:
        return [self.get_name(index) for index in mapped_list]


def select_exemplars(
    data_matrix, cluster_assignments, input_weights: list[int]
) -> tuple[list[int], list[int]]:
    """Selects exemplars from a set of cluster assignments. Returns a list of exemplar indices and the corresponding
    count of members in that cluster."""
    clusters = defaultdict(list)
    for i, cluster in enumerate(cluster_assignments):
        clusters[cluster].append(i)

    exemplars = []
    exemplar_weights = []
    for members in clusters.values():
        exemplar_weights.append(sum(input_weights[i] for i in members))
        if len(members) == 1:
            exemplars.append(members[0])
        else:
            cluster_matrix = data_matrix[np.ix_(members, members)]
            avg_distances = np.mean(cluster_matrix, axis=1)
            exemplars.append(members[np.argmin(avg_distances)])

    return exemplars, exemplar_weights


def mash_fields(
    library: Path, mash_path: str, fields: Iterable[int], threshold: float = 0.005, threads: int = 1
):
    """Returns the selected fields from the mash result."""
    # Run self-self mash comparison
    # Get the number of sketches
    for row in self_search(mash_path, str(library), distance_threshold=threshold, threads=threads):
        yield (
            (
                g_pair := sorted(
                    [(row[0].replace(".fna.gz", "")), (row[1].replace(".fna.gz", ""))]
                )
            )[0],
            g_pair[1],
            *[row[field] for field in fields],
        )


def mash_distances(
    library: Path, mash_path: str, threshold: float = 0.005, threads: int = 1
) -> Iterator[tuple[str, str, float]]:
    """Returns the distances at a given threshold."""
    # Run self-self mash comparison
    # Get the number of sketches
    for result in mash_fields(library, mash_path, [2], threshold=threshold, threads=threads):
        yield result[0], result[1], float(result[2])


def get_exemplars(
    scores: dict[tuple[str | int, str | int], float],
    threshold: float,
    input_weights: list[int] = None,
) -> tuple[list[str | int], list[int]]:
    if len(scores) < 1:
        raise ValueError("No data provided to cluster.")

    mapped_scores, bimap = BiMap.map_score_dict(scores)
    matrix_size = len(bimap)

    if input_weights is None:
        input_weights = [1] * matrix_size

    if matrix_size == 1:
        return [bimap.get_name(0)], [input_weights[0]]

    input_matrix = np.ones((matrix_size, matrix_size), dtype=np.float32)
    np.fill_diagonal(input_matrix, 0)
    for id_pair, score in mapped_scores.items():
        input_matrix[id_pair[0]][id_pair[1]] = score
        input_matrix[id_pair[1]][id_pair[0]] = score

    # Convert scores to a distance matrix
    # cluster to the threshold (multi-linkage)
    clustering = AgglomerativeClustering(
        metric="precomputed",
        distance_threshold=threshold,
        linkage="complete",
        n_clusters=None,
    ).fit(input_matrix)

    # select exemplar for each cluster
    raw_exemplars, ex_weights = select_exemplars(
        input_matrix, clustering.labels_, input_weights
    )
    remapped_exemplars: list[str | int] = bimap.invert_list(raw_exemplars)
    return remapped_exemplars, ex_weights
