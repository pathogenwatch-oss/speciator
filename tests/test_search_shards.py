from pathlib import Path

import polars as pl
import pytest

from speciator import search as search_mod
from speciator.search import (
    SearchEnvironment,
    _top_n_matches,
    library_shard_paths,
    run_search,
)


def make_search_environment(library_dir: Path) -> SearchEnvironment:
    return SearchEnvironment(
        config_file=library_dir / "config.toml",
        library_dir=library_dir,
        mash="mash",
        num_matches={"Bacteria": 10},
        thresholds={"Bacteria": 0.05},
        mash_processes=4,
        sketch_extension=".msh",
        info_extension=".pqt",
    )


def test_library_shard_paths_single_library(tmp_path: Path):
    (tmp_path / "Bacteria.msh").touch()
    (tmp_path / "Bacteria.pqt").touch()

    env = make_search_environment(tmp_path)
    shards = library_shard_paths("Bacteria", env)

    assert shards == [(tmp_path / "Bacteria.pqt", tmp_path / "Bacteria.msh")]


def test_library_shard_paths_discovers_and_orders_numeric_suffixes(tmp_path: Path):
    for idx in [10, 2, 1]:
        (tmp_path / f"Bacteria.{idx}.msh").touch()
        (tmp_path / f"Bacteria.{idx}.pqt").touch()

    env = make_search_environment(tmp_path)
    shards = library_shard_paths("Bacteria", env)

    assert [sketch.name for _, sketch in shards] == [
        "Bacteria.1.msh",
        "Bacteria.2.msh",
        "Bacteria.10.msh",
    ]


def test_library_shard_paths_raises_on_missing_metadata_pair(tmp_path: Path):
    (tmp_path / "Bacteria.1.msh").touch()

    env = make_search_environment(tmp_path)
    with pytest.raises(FileNotFoundError, match="Missing metadata file"):
        library_shard_paths("Bacteria", env)


def test_top_n_matches_legacy_preserves_first_seen_order_for_ties():
    merged_matches = pl.DataFrame(
        {
            "filename": ["A", "B", "C", "D", "E"],
            "distance": [0.01, 0.01, 0.01, 0.02, 0.02],
            "p-value": [0.0] * 5,
            "shared-hashes": ["1000/1000"] * 5,
        }
    )

    top = _top_n_matches(merged_matches, 3)

    assert top["filename"].to_list() == ["A", "B", "C"]


def _make_lineages() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "LinkCode": ["101", "202"],
            "OrganismId": ["101", "202"],
            "OrganismName": ["Species A", "Species C"],
            "SpeciesId": ["101", "202"],
            "SpeciesName": ["Species A", "Species C"],
            "GenusId": ["11", "22"],
            "GenusName": ["GenusA", "GenusC"],
            "SuperkingdomId": ["2", "2"],
            "SuperkingdomName": ["Bacteria", "Bacteria"],
        }
    )


def _make_metadata() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "record_key": ["A", "B", "C", "D"],
            "organism_code": ["101", "101", "202", "202"],
            "accession": ["A", "B", "C", "D"],
            "represents": [1, 1, 1, 1],
        }
    )


def test_run_search_split_and_unsplit_equivalent(monkeypatch, tmp_path: Path):
    _make_lineages().write_parquet(tmp_path / "lineages.pqt")

    unsplit_env = SearchEnvironment(
        config_file=tmp_path / "cfg.toml",
        library_dir=tmp_path,
        mash="mash",
        num_matches={"TestLib": 3},
        thresholds={"TestLib": 0.05},
        mash_processes=4,
    )
    split_env = SearchEnvironment(
        config_file=tmp_path / "cfg.toml",
        library_dir=tmp_path,
        mash="mash",
        num_matches={"TestLibSplit": 3},
        thresholds={"TestLibSplit": 0.05},
        mash_processes=4,
    )

    (tmp_path / "TestLib.msh").touch()
    _make_metadata().write_parquet(tmp_path / "TestLib.pqt")

    (tmp_path / "TestLibSplit.1.msh").touch()
    (tmp_path / "TestLibSplit.2.msh").touch()
    _make_metadata().slice(0, 2).write_parquet(tmp_path / "TestLibSplit.1.pqt")
    _make_metadata().slice(2, 2).write_parquet(tmp_path / "TestLibSplit.2.pqt")

    def fake_best_matches(
        sample_sketch,
        ref_seq_sketch,
        mash_path="",
        number_of_best_matches=10,
        distance_threshold=0.5,
        threads=None,
    ):
        lib = Path(ref_seq_sketch).name
        if lib == "TestLib.msh":
            return pl.DataFrame(
                {
                    "filename": ["A", "B", "C"],
                    "distance": [0.01, 0.02, 0.03],
                    "p-value": [0.0, 0.0, 0.0],
                    "shared-hashes": ["1000/1000", "980/1000", "960/1000"],
                }
            )
        if lib == "TestLibSplit.1.msh":
            return pl.DataFrame(
                {
                    "filename": ["A", "D"],
                    "distance": [0.01, 0.04],
                    "p-value": [0.0, 0.0],
                    "shared-hashes": ["1000/1000", "920/1000"],
                }
            )
        if lib == "TestLibSplit.2.msh":
            return pl.DataFrame(
                {
                    "filename": ["B", "C"],
                    "distance": [0.02, 0.03],
                    "p-value": [0.0, 0.0],
                    "shared-hashes": ["980/1000", "960/1000"],
                }
            )
        raise AssertionError(f"Unexpected library sketch: {lib}")

    monkeypatch.setattr(search_mod, "get_best_mash_matches", fake_best_matches)

    query = tmp_path / "query.msh"
    query.touch()
    unsplit_assignment, unsplit_matches = run_search(query, "TestLib", unsplit_env)
    split_assignment, split_matches = run_search(query, "TestLibSplit", split_env)

    assert unsplit_assignment["OrganismId"] == split_assignment["OrganismId"]
    assert unsplit_assignment["OrganismName"] == split_assignment["OrganismName"]
    assert unsplit_assignment["distance"] == split_assignment["distance"]
    assert unsplit_matches.select("filename", "distance").rows() == split_matches.select(
        "filename", "distance"
    ).rows()


def test_run_search_uses_single_mash_thread_per_shard(monkeypatch, tmp_path: Path):
    _make_lineages().write_parquet(tmp_path / "lineages.pqt")
    env = SearchEnvironment(
        config_file=tmp_path / "cfg.toml",
        library_dir=tmp_path,
        mash="mash",
        num_matches={"TestLib": 2},
        thresholds={"TestLib": 0.05},
        mash_processes=8,
    )
    (tmp_path / "TestLib.1.msh").touch()
    (tmp_path / "TestLib.2.msh").touch()
    _make_metadata().slice(0, 2).write_parquet(tmp_path / "TestLib.1.pqt")
    _make_metadata().slice(2, 2).write_parquet(tmp_path / "TestLib.2.pqt")

    seen_threads: list[int | None] = []

    def fake_best_matches(
        sample_sketch,
        ref_seq_sketch,
        mash_path="",
        number_of_best_matches=10,
        distance_threshold=0.5,
        threads=None,
    ):
        seen_threads.append(threads)
        return pl.DataFrame(
            {
                "filename": ["A"],
                "distance": [0.01],
                "p-value": [0.0],
                "shared-hashes": ["1000/1000"],
            }
        )

    monkeypatch.setattr(search_mod, "get_best_mash_matches", fake_best_matches)
    query = tmp_path / "query.msh"
    query.touch()
    run_search(query, "TestLib", env)

    assert seen_threads == [1, 1]
