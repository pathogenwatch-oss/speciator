from pathlib import Path

import polars as pl
import pytest

from speciator import run_split_library as split_mod


def _write_library_records(path: Path, record_keys: list[str]) -> None:
    pl.DataFrame(
        {
            "record_key": record_keys,
            "organism_code": ["1"] * len(record_keys),
            "accession": record_keys,
            "represents": [1] * len(record_keys),
        }
    ).write_parquet(path)


def _touch_sketches(sketch_dir: Path, names: list[str]) -> None:
    sketch_dir.mkdir(parents=True, exist_ok=True)
    for name in names:
        (sketch_dir / f"{name}.msh").touch()


def test_split_generates_shards_and_archives_in_place(monkeypatch, tmp_path: Path):
    library_dir = tmp_path / "library"
    build_space = tmp_path / "build_space"
    mash_dir = build_space / "mash"
    library_dir.mkdir()
    _touch_sketches(mash_dir, ["A", "B", "C", "D"])

    (library_dir / "Bacteria.msh").write_text("original")
    _write_library_records(library_dir / "Bacteria.pqt", ["A", "B", "C", "D"])

    pasted: list[tuple[str, list[str]]] = []

    def fake_paste(mash_path: str, sketch_files: list[Path], output_file: Path):
        pasted.append((output_file.name, [Path(s).name for s in sketch_files]))
        output_file.write_text("shard")

    monkeypatch.setattr(split_mod, "paste_sketches", fake_paste)

    split_mod.split(
        library_name="Bacteria",
        max_sub_library_size=2,
        library_dir=library_dir,
        build_space_dir=build_space,
        output_dir=library_dir,
        mash_path="mash",
        overwrite=False,
        keep_original=False,
    )

    assert (library_dir / "Bacteria.1.msh").exists()
    assert (library_dir / "Bacteria.2.msh").exists()
    assert (library_dir / "Bacteria.1.pqt").exists()
    assert (library_dir / "Bacteria.2.pqt").exists()
    assert (library_dir / "Bacteria_full.msh").exists()
    assert (library_dir / "Bacteria_full.pqt").exists()
    assert not (library_dir / "Bacteria.msh").exists()
    assert not (library_dir / "Bacteria.pqt").exists()

    assert pasted == [
        ("Bacteria.1.msh", ["A.msh", "B.msh"]),
        ("Bacteria.2.msh", ["C.msh", "D.msh"]),
    ]

    shard1 = pl.read_parquet(library_dir / "Bacteria.1.pqt")
    shard2 = pl.read_parquet(library_dir / "Bacteria.2.pqt")
    assert shard1["record_key"].to_list() == ["A", "B"]
    assert shard2["record_key"].to_list() == ["C", "D"]


def test_split_raises_if_sketch_missing(monkeypatch, tmp_path: Path):
    library_dir = tmp_path / "library"
    build_space = tmp_path / "build_space"
    mash_dir = build_space / "mash"
    library_dir.mkdir()
    _touch_sketches(mash_dir, ["A"])

    (library_dir / "Viruses.msh").write_text("original")
    _write_library_records(library_dir / "Viruses.pqt", ["A", "B"])

    monkeypatch.setattr(split_mod, "paste_sketches", lambda *args, **kwargs: None)

    with pytest.raises(FileNotFoundError, match="Missing 1 sketch files"):
        split_mod.split(
            library_name="Viruses",
            max_sub_library_size=1,
            library_dir=library_dir,
            build_space_dir=build_space,
            output_dir=library_dir,
            mash_path="mash",
            overwrite=False,
            keep_original=True,
        )


def test_split_noop_when_library_is_small(monkeypatch, tmp_path: Path):
    library_dir = tmp_path / "library"
    build_space = tmp_path / "build_space"
    mash_dir = build_space / "mash"
    library_dir.mkdir()
    _touch_sketches(mash_dir, ["A", "B"])

    (library_dir / "Archaea.msh").write_text("original")
    _write_library_records(library_dir / "Archaea.pqt", ["A", "B"])

    called = {"paste": 0}

    def fake_paste(*args, **kwargs):
        called["paste"] += 1

    monkeypatch.setattr(split_mod, "paste_sketches", fake_paste)

    split_mod.split(
        library_name="Archaea",
        max_sub_library_size=2,
        library_dir=library_dir,
        build_space_dir=build_space,
        output_dir=library_dir,
        mash_path="mash",
        overwrite=False,
        keep_original=True,
    )

    assert called["paste"] == 0
    assert (library_dir / "Archaea.msh").exists()
    assert (library_dir / "Archaea.pqt").exists()
    assert not (library_dir / "Archaea_full.msh").exists()
