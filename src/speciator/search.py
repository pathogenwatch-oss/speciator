from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path

import polars as pl
import toml
from polars import DataFrame

from speciator.mash_functions import get_best_mash_matches

logger = logging.getLogger(__name__)

HEAVY_WEIGHT = 1000000
AMBIGUOUS_SP = [" sp.", " bacterium"]


@dataclass
class SearchEnvironment:
    """Search configuration"""

    config_file: Path
    library_dir: Path
    mash: str
    num_matches: dict[str, int]
    thresholds: dict[str, float]
    mash_processes: int
    sketch_extension: str = ".msh"
    info_extension: str = ".pqt"
    lineage_file: Path = field(init=False)

    def __post_init__(self):
        self.lineage_file = self.library_dir / "lineages.pqt"

    @staticmethod
    def default_search_environment() -> SearchEnvironment:
        return SearchEnvironment(
            config_file=Path("config.toml"),
            library_dir=Path("library"),
            mash="mash",
            num_matches={"superkingdom": 10},
            sketch_extension=".msh",
            info_extension=".pqt",
            thresholds={"superkingdom": 0.05},
            mash_processes=1,
        )

    @staticmethod
    def from_config_file(
        config_file: Path, mash_processes: int, library: Path = None
    ) -> SearchEnvironment:
        try:
            config = toml.load(config_file)
            logger.info(f"Successfully loaded configuration from {config_file}")
        except Exception as e:
            logger.error(f"Failed to load configuration from {config_file}: {str(e)}")
            raise e
        if library is not None:
            logger.info(
                f"Overriding library directory with provided path {str(library)}"
            )
        return SearchEnvironment(
            config_file=config_file,
            library_dir=Path(config["paths"]["library"])
            if library is None
            else library,
            mash=config["paths"]["mash"],
            num_matches=config["num_matches"],
            sketch_extension=config["parameters"]["sketch_extension"],
            info_extension=config["parameters"]["info_extension"],
            thresholds=config["thresholds"],
            mash_processes=mash_processes,
        )


@dataclass
class Confidence:
    totalMatches: int = 0
    speciesCounts: dict[str, int] = field(default_factory=dict)
    consistency: float = 0.0
    diversity: int = 0
    matches: list[dict[str, float | int | str]] = field(default_factory=list)

    @staticmethod
    def default_confidence() -> Confidence:
        return Confidence()

    @staticmethod
    def build(species_name: str, matches: pl.DataFrame) -> Confidence:
        total_matches: int = matches.shape[0]
        if total_matches == 0:
            return Confidence.default_confidence()
        species_counts_df = matches.group_by("OrganismName").agg(
            pl.col("adjusted_represents").sum().alias("count")
        )

        # Convert to dictionary
        species_counts: dict[str, int] = dict(
            zip(
                species_counts_df["OrganismName"].to_list(),
                species_counts_df["count"].to_list(),
            )
        )
        total_representatives = species_counts_df["count"].sum()
        consistency: float = round(species_counts[species_name] / total_representatives, 2)
        diversity: int = len(species_counts)
        reduced_matches = [
            {
                "distance": float(row["distance"]),
                "assembly_accession": str(row["accession"]),
                "represents": int(row["adjusted_represents"]),
                "organism": str(row["OrganismName"]),
            }
            for row in matches.sort("distance").iter_rows(named=True)
        ]
        return Confidence(
            totalMatches=total_matches,
            speciesCounts=species_counts,
            consistency=consistency,
            diversity=diversity,
            matches=reduced_matches,
        )


@dataclass
class SearchResult:
    taxId: str
    taxName: str
    speciesId: str
    speciesName: str
    genusId: str
    genusName: str
    superkingdomId: str
    superkingdomName: str
    referenceId: str
    mashDistance: float
    pValue: float
    matchingHashes: str
    confidence: str
    source: str
    confidence: Confidence

    @staticmethod
    def default_result() -> SearchResult:
        return SearchResult(
            taxId="32644",
            taxName="Unidentified",
            speciesId="32644",
            speciesName="Unidentified",
            genusId="9999999",
            genusName="Unclassified",
            superkingdomId="12908",
            superkingdomName="Unclassified Sequences",
            referenceId="None found",
            mashDistance=1.0,
            pValue=0,
            matchingHashes="0/1000",
            confidence=Confidence.default_confidence(),
            source="None",
        )


def determine_threshold(sorted_matches: pl.DataFrame, max_threshold: int):
    """Determine the dynamic threshold based on species representation"""
    if max_threshold == 1:
        return max_threshold  # No threshold adjustment needed if max_threshold is 1
    if max_threshold < 0:
        raise ValueError("max_threshold must be a positive integer")
    # Initialise the threshold to the maximum allowed
    threshold: int = max_threshold

    # Step through matches in distance order up to max_threshold
    for hit_index in range(min(len(sorted_matches), max_threshold)):
        # Get the species count for this hit (number of representatives for this species)
        species_count = sorted_matches[hit_index, "weighted_len"]

        # If species count is less than the current threshold, we may need to adjust
        if species_count < threshold:
            # If species count <= hit index (0-based), set {threshold} to hit index + 1
            if (
                species_count <= hit_index + 1
            ):  # +1 because hit_index is 0-based, but we want 1-based positioning
                threshold = hit_index + 1
                break
            # Otherwise, set {threshold} to the species count
            else:
                threshold = species_count
                break
        # If species count >= threshold, continue to next match
    return threshold


def remove_ambiguous_taxa(matches: pl.DataFrame) -> pl.DataFrame:
    """
    If there is a species with full assignment (i.e. not 'bacterium' or 'sp.'), then all ambiguous species are removed.
    """
    if "OrganismName" not in matches.columns:
        return matches

    # Define ambiguous patterns
    ambiguous_patterns = r" sp\.| bacterium"

    # Check if there are any unambiguous species names
    unambiguous_present = (
        matches.filter(~pl.col("OrganismName").str.contains(ambiguous_patterns)).height
        > 0
    )

    if unambiguous_present:
        # If unambiguous names are present, filter out the ambiguous ones
        return matches.filter(~pl.col("OrganismName").str.contains(ambiguous_patterns))
    else:
        # Otherwise, return the original dataframe
        return matches


def get_most_frequent_species_match(
    matches: pl.DataFrame, taxonomy_links: pl.DataFrame, max_threshold: int = 10
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Use polars to merge the best match file with ref species info and report the most frequent species
    return species and count, using the top matches until the threshold is reached
    """
    taxonomy_links = taxonomy_links.rename({"record_key": "filename"})

    if matches.is_empty():
        matches = matches.join(taxonomy_links, on="filename")
        return matches.with_columns(
            [
                pl.lit(None).alias("most_frequent_species_name"),
                pl.lit(None).alias("most_frequent_species_count"),
                pl.lit(None).alias("count"),
            ]
        ), matches
    # else:
    # We need to cover the case where there a few representative for a species that is closely related
    # to a more heavily sampled one (e.g. V. cholerae & V. tarriae). In this case the threshold is set to 2n-1 if
    # 2n-1 is less than the default threshold.
    # Process for determining the threshold.
    # Set {threshold} = max_threshold
    # Step through the matches in sorted order up to the max threshold.
    # If the hit has {species count} < {threshold}:
    #    if {species count} <= {hit index}: then the threshold is set to {hit index} + 1.
    #    else {threshold} is set to {species count}.
    # else continue to the next match.
    matches = matches.join(taxonomy_links, on="filename")
    if matches.is_empty():
        raise ValueError(f"No species links found for matches in {str(matches)}")
    # The values of "n" for each species.
    matches = matches.join(
        taxonomy_links.group_by("LinkCode").len(), on="LinkCode", how="left"
    )

    # Step 2: Get the list of matches up to the threshold, taking into account how many cluster members each
    # match represents.
    # Sort matches by distance and calculate the cumulative sum of 'represents'
    matches = add_weights(matches)
    threshold = determine_threshold(matches, max_threshold)

    matches = matches.drop("weighted_len").with_columns(
        pl.col("represents").cum_sum().alias("cumulative_represents")
    )

    # Find the first row that exceeds the threshold
    threshold_exceeded = matches.with_row_index().filter(
        pl.col("cumulative_represents") > threshold
    )

    if threshold_exceeded.is_empty():
        # If no row exceeds the threshold, use all matches
        filtered_matches = matches
    else:
        # Include all rows up to and including the first one that equals or exceeds the threshold
        threshold_index = threshold_exceeded["index"][0]
        if matches[threshold_index - 1, "cumulative_represents"] == threshold:
            # Adjust as the previous match equals the threshold, so the extra match is excluded.
            threshold_index = threshold_index - 1
        filtered_matches = matches.head(threshold_index + 1)

    # Adjust the 'represents' value of the last match to match the threshold. Otherwise, the last match could
    # be overweighted.
    last_match_adjustment = (
        threshold - filtered_matches[-2, "cumulative_represents"]
        if filtered_matches.height > 1
        else threshold
    )

    filtered_matches = filtered_matches.with_columns(
        pl.when(pl.col("cumulative_represents") > threshold)
        .then(last_match_adjustment)
        .otherwise(pl.col("represents"))
        .alias("adjusted_represents")
    )

    species_filtered_matches = remove_ambiguous_taxa(filtered_matches)

    # Calculate species counts using the filtered and adjusted matches
    # and sort by the most represented species first.
    species_counts: pl.DataFrame = (
        species_filtered_matches.group_by("LinkCode")
        .agg(pl.col("adjusted_represents").sum().alias("count"))
        .sort("count", descending=True)
    )

    if (
        species_counts.height == 1
        or species_counts["count"][0] != species_counts["count"][1]
    ):
        # Only one species, or one species with the highest representation
        most_frequent_species = species_counts["LinkCode"][0]
        most_frequent_species_count = species_counts["count"][0]
    else:
        equal_often = species_counts.filter(
            pl.col("count") == species_counts["count"][0]
        )["LinkCode"].to_list()
        best_match = (
            pl.concat(
                [
                    filtered_matches.filter(pl.col("LinkCode") == test_name)
                    .sort("distance")
                    .head(1)
                    for test_name in equal_often
                ]
            )
            .sort("distance")
            .head(1)
        )
        most_frequent_species = best_match["LinkCode"][0]
        most_frequent_species_count = species_counts["count"][0]

    # Get top hit of the most frequent species as measured by distance
    top_hit = (
        filtered_matches.filter(pl.col("LinkCode") == most_frequent_species)
        .sort("distance")
        .head(1)
    )

    top_hit = top_hit.with_columns(
        [
            pl.lit(most_frequent_species).alias("most_frequent_species_name"),
            pl.lit(most_frequent_species_count).alias("most_frequent_species_count"),
            pl.lit(filtered_matches.height).alias("count"),
        ]
    )

    return top_hit, filtered_matches


def add_weights(matches: DataFrame) -> DataFrame:
    matches = matches.sort("distance", descending=False).with_columns(
        pl.when(
            (pl.col("SpeciesName").str.contains(AMBIGUOUS_SP[0]))
            | (pl.col("SpeciesName").str.contains(AMBIGUOUS_SP[1]))
        )
        .then(HEAVY_WEIGHT)
        .otherwise(pl.col("len"))
        .alias("weighted_len")
    )
    return matches


def build_path(library: str, directory: Path | str, extension: str) -> Path:
    return Path(os.path.join(directory, f"{library}{extension}"))


def library_paths(
    library: str, search_environment: SearchEnvironment
) -> tuple[Path, Path]:
    return build_path(
        library, search_environment.library_dir, search_environment.info_extension
    ), build_path(
        library, search_environment.library_dir, search_environment.sketch_extension
    )


def _library_shard_sort_key(
    path: Path, library: str, extension: str
) -> tuple[int, int | str]:
    stem = path.name[: -len(extension)]
    shard_suffix = stem[len(library) + 1 :]
    if shard_suffix.isdigit():
        return 0, int(shard_suffix)
    return 1, shard_suffix


def library_shard_paths(
    library: str, search_environment: SearchEnvironment
) -> list[tuple[Path, Path]]:
    """
    Resolve one logical library to one or more physical sketch+metadata shard pairs.
    """
    single_info, single_sketch = library_paths(library, search_environment)
    if single_info.exists() and single_sketch.exists():
        return [(single_info, single_sketch)]
    if single_info.exists() != single_sketch.exists():
        raise FileNotFoundError(
            f"Incomplete library files for {library}: expected both {single_info.name} and {single_sketch.name}"
        )

    shard_pattern = f"{library}.*{search_environment.sketch_extension}"
    sketch_shards = sorted(
        search_environment.library_dir.glob(shard_pattern),
        key=lambda p: _library_shard_sort_key(
            p, library, search_environment.sketch_extension
        ),
    )
    if not sketch_shards:
        raise FileNotFoundError(
            f"No sketch files found for library {library} in {search_environment.library_dir}"
        )

    shards: list[tuple[Path, Path]] = []
    for sketch_path in sketch_shards:
        stem = sketch_path.name[: -len(search_environment.sketch_extension)]
        info_path = (
            search_environment.library_dir
            / f"{stem}{search_environment.info_extension}"
        )
        if not info_path.exists():
            raise FileNotFoundError(
                f"Missing metadata file for shard {sketch_path.name}: expected {info_path.name}"
            )
        shards.append((info_path, sketch_path))
    return shards


def _top_n_matches(matches: pl.DataFrame, n: int) -> pl.DataFrame:
    """
    Keep top-n rows by distance while preserving first-seen order for ties.
    """
    if matches.is_empty() or n <= 0:
        return matches.head(0)
    return (
        matches.with_row_index(name="match_order")
        .sort(by=["distance", "match_order"], descending=[False, False])
        .drop("match_order")
        .head(n)
    )


def run_search(
    sketch_path: Path,
    library: str,
    search_environment: SearchEnvironment,
) -> tuple[None | dict[str, str | float], None | pl.DataFrame]:
    """
    Run search on a single sketch
    :param sketch_path: path to the sketch
    :param library: library name
    :param search_environment: search configuration
    :return: The top reference from the matching species plus the filtered matches
    """
    shards = library_shard_paths(library, search_environment)
    lineage = pl.read_parquet(search_environment.lineage_file)
    taxon_links = pl.concat(
        [
            pl.read_parquet(
                info_path,
                schema={
                    "record_key": pl.Utf8,
                    "organism_code": pl.Utf8,
                    "accession": pl.Utf8,
                    "represents": pl.Int64,
                },
            )
            for info_path, _ in shards
        ],
        how="vertical",
    )
    taxon_links = taxon_links.rename({"organism_code": "LinkCode"}).join(
        lineage, on="LinkCode", how="left"
    )

    n_matches = search_environment.num_matches[library]
    max_workers = min(len(shards), max(1, search_environment.mash_processes))
    mash_threads_per_process = 1

    shard_matches: list[pl.DataFrame] = [pl.DataFrame()] * len(shards)

    def run_shard(shard_idx: int, shard_sketch: Path) -> tuple[int, pl.DataFrame]:
        return shard_idx, get_best_mash_matches(
            sketch_path,
            shard_sketch,
            search_environment.mash,
            n_matches,
            search_environment.thresholds[library],
            mash_threads_per_process,
        )

    if len(shards) == 1:
        _, only_sketch = shards[0]
        shard_matches[0] = run_shard(0, only_sketch)[1]
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(run_shard, idx, sketch_lib)
                for idx, (_, sketch_lib) in enumerate(shards)
            ]
            for future in futures:
                idx, matches = future.result()
                shard_matches[idx] = matches

    sample_matches = _top_n_matches(
        pl.concat(shard_matches, how="vertical"), n_matches
    )

    species_assignment, matches = get_most_frequent_species_match(
        sample_matches, taxon_links, n_matches
    )

    if not species_assignment.is_empty():
        species_assignment = species_assignment.with_columns(
            [
                pl.lit(os.path.basename(sketch_path).replace(".msh", "")).alias("file"),
                (pl.col("most_frequent_species_count") / pl.col("count") * 100).alias(
                    "percentage"
                ),
            ]
        )
        species_assignment = species_assignment.rename(
            {
                "percentage": f"%_of_{n_matches}_best_matches=species"
            }
        )

    # species_assignment = species_assignment.rename({"organism_code": "OrganismId"})
    # species_assignment = species_assignment.rename({"OrganismId": "LinkCode"})
    # species_assignment = species_assignment.join(pl.read_parquet(search_environment.lineage_file), on="LinkCode", how="left")

    return species_assignment.row(
        0, named=True
    ) if not species_assignment.is_empty() else None, matches


def assign_species(
    sketch_path: Path, search_environment: SearchEnvironment
) -> SearchResult:
    run = partial(
        run_search, sketch_path=sketch_path, search_environment=search_environment
    )
    # Libraries are run in the order provided in the configuration file.
    for library, thresholds in search_environment.thresholds.items():
        species_assignment, matches = run(
            library=library,
        )
        if (
            species_assignment is not None
            and species_assignment["OrganismId"] is not None
            and species_assignment["OrganismName"] != "Unknown"
            and not (
                library == "Kleborate"
                and not (
                    species_assignment["GenusName"].startswith("Klebsiella")
                    or species_assignment["GenusName"].startswith("Salmonella")
                    or species_assignment["GenusName"].startswith("Raoultella")
                    or species_assignment["SpeciesName"]
                    == "Enterobacter quasihormaechei"
                )
            )
        ):
            result = format_result(species_assignment, matches, library)
            break
    else:
        result = SearchResult.default_result()
    return result


def format_result(
    result: dict[str, str | float], matches: pl.DataFrame, library_name: str
) -> SearchResult:
    if result is None or matches is None:
        return SearchResult.default_result()
    else:
        return SearchResult(
            str(result["OrganismId"]) if "OrganismId" in result else "",
            result["OrganismName"] if "OrganismName" in result else "",
            str(result["SpeciesId"]) if "SpeciesId" in result else "",
            result["SpeciesName"] if "SpeciesName" in result else "",
            str(result["GenusId"])
            if "GenusId" in result and result["GenusId"] != ""
            else "9999999",
            result["GenusName"]
            if "GenusName" in result and result["GenusName"] != ""
            else "Unclassified",
            str(result["SuperkingdomId"]) if "SuperkingdomId" in result else "",
            result["SuperkingdomName"] if "SuperkingdomName" in result else "",
            str(result["accession"]) if "filename" in result else "",
            float(result["distance"]) if "distance" in result else 1.0,
            float(result["p-value"]) if "p-value" in result else 0.0,
            result["shared-hashes"] if "shared-hashes" in result else "",
            Confidence.build(result["OrganismName"], matches)
            if matches is not None
            else Confidence.default_confidence(),
            library_name,
        )
