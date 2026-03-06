from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import numpy as np
from sklearn.cluster import AgglomerativeClustering

from speciator.linking import Linker
from speciator.mash_functions import paste_sketches
from speciator.ncbi_utils import BiMap, mash_distances, select_exemplars

logger = logging.getLogger(__name__)

Item = tuple[str, Path]  # (id_used_in_clustering, sketch_path)


def build_linked_groups(
    *,
    batch_size: int,
    temp_lib: Path,
    threshold: float,
    mash_path: str,
    search_threads: int =  1
) -> Linker:
    """
    Build connected components ("linked groups") of items based on pairwise Mash distances
    <= threshold, batching edges into the Linker.

    `temp_lib` must be a compiled mash library file (e.g. created via `paste_sketches`).
    """
    linked_groups = Linker()
    scores_batch: list[tuple[str, str, float]] = []

    for id1, id2, distance in mash_distances(temp_lib, mash_path, threshold=threshold, threads=search_threads):
        scores_batch.append((str(id1), str(id2), float(distance)))
        if len(scores_batch) >= batch_size:
            linked_groups.add_scores_batch(scores_batch)
            scores_batch = []

    if scores_batch:
        linked_groups.add_scores_batch(scores_batch)

    return linked_groups


def get_exemplars_weighted(
    scores: dict[tuple[str, str], float],
    *,
    threshold: float,
    weight_by_id: dict[str, int] | None = None,
) -> tuple[list[str], list[int]]:
    """
    Weighted wrapper around the same core approach used elsewhere:
    - map ids to dense indices
    - build a precomputed distance matrix
    - agglomerative clustering (complete linkage) at the given threshold
    - choose an exemplar per cluster (min average intra-cluster distance)
    - compute exemplar weights by summing member weights

    Unlike passing `input_weights` positionally, this function aligns weights by id
    to avoid relying on incidental dict insertion order.
    """
    if len(scores) < 1:
        raise ValueError("No data provided to cluster.")

    mapped_scores, bimap = BiMap.map_score_dict(scores)
    matrix_size = len(bimap)

    if matrix_size == 1:
        only_id = str(bimap.get_name(0))
        only_weight = 1 if weight_by_id is None else int(weight_by_id.get(only_id, 1))
        return [only_id], [only_weight]

    if weight_by_id is None:
        input_weights = [1] * matrix_size
    else:
        input_weights = [
            int(weight_by_id.get(str(bimap.get_name(i)), 1)) for i in range(matrix_size)
        ]

    input_matrix = np.ones((matrix_size, matrix_size), dtype=np.float32)
    np.fill_diagonal(input_matrix, 0.0)
    for (i, j), dist in mapped_scores.items():
        input_matrix[i][j] = dist
        input_matrix[j][i] = dist

    clustering = AgglomerativeClustering(
        metric="precomputed",
        distance_threshold=threshold,
        linkage="complete",
        n_clusters=None,
    ).fit(input_matrix)

    raw_exemplars, exemplar_weights = select_exemplars(
        input_matrix, clustering.labels_, input_weights
    )
    exemplar_ids = [str(bimap.get_name(i)) for i in raw_exemplars]
    return exemplar_ids, [int(w) for w in exemplar_weights]


def _make_temp_lib_path(*, tmp_dir: Path | None, prefix: str = "reps_") -> Path:
    if tmp_dir is not None:
        tmp_dir.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            dir=str(tmp_dir), prefix=prefix, suffix=".msh", delete=False
        )
    else:
        handle = tempfile.NamedTemporaryFile(prefix=prefix, suffix=".msh", delete=False)

    # We only want the name; mash will create/overwrite the file.
    handle.close()
    return Path(handle.name)


def extract_exemplars_from_items(
    *,
    items: list[Item],
    threshold: float,
    batch_size: int,
    mash_path: str,
    previous_counts: list[int] | None = None,
    tmp_dir: Path | None = None,
    search_threads: int = 1
) -> tuple[list[str], list[int]]:
    """
    Given a set of sketch items, build a temporary combined mash library, compute linked groups
    at `threshold`, and select exemplars per group.

    Returns:
      - rep_ids: list[str]
      - rep_counts: list[int] where each entry is the number of originals represented by that rep
    """
    if not items:
        raise ValueError("No items provided.")

    if previous_counts is None:
        previous_counts = [1] * len(items)
    if len(previous_counts) != len(items):
        raise ValueError(
            f"previous_counts length ({len(previous_counts)}) does not match items ({len(items)})."
        )

    logger.info(f"Initial rep count before clustering: {len(items)}")
    logger.debug(f"threshold={threshold}, batch_size={batch_size}")

    id_to_path = {item_id: sketch_path for item_id, sketch_path in items}
    weight_by_id = {item_id: int(previous_counts[i]) for i, (item_id, _) in enumerate(items)}

    temp_lib = _make_temp_lib_path(tmp_dir=tmp_dir, prefix="cluster_")
    try:
        paste_sketches(mash_path, [path for _, path in items], temp_lib)
        linked_groups = build_linked_groups(
            batch_size=batch_size, temp_lib=temp_lib, threshold=threshold, mash_path=mash_path, search_threads=search_threads
        )
    finally:
        try:
            if temp_lib.exists():
                temp_lib.unlink()
        except OSError:
            logger.warning("Failed to delete temporary mash library %s", temp_lib)

    rep_ids: list[str] = []
    rep_counts: list[int] = []

    for group in linked_groups.groups:
        group_scores: dict[tuple[str, str], float] = {
            (str(a), str(b)): float(d) for (a, b), d in group.scores.items()
        }
        group_exemplars, group_weights = get_exemplars_weighted(
            group_scores, threshold=threshold, weight_by_id=weight_by_id
        )
        rep_ids.extend(group_exemplars)
        rep_counts.extend(group_weights)

    # Handle any isolated nodes that never appeared in scores (e.g. if mash emits no distances)
    # In normal operation, self-self dist should include enough edges, but we keep this safety net.
    members_in_any_group = group_members(linked_groups)
    seen = set(rep_ids)
    for item_id, _ in items:
        if item_id not in members_in_any_group:
            # Not in any group => singleton
            if item_id not in seen:
                rep_ids.append(item_id)
                rep_counts.append(int(weight_by_id.get(item_id, 1)))
                seen.add(item_id)

    # Ensure stable output order matches input order when possible.
    order = {item_id: i for i, (item_id, _) in enumerate(items)}
    rep_ids_counts = sorted(zip(rep_ids, rep_counts), key=lambda rc: order.get(rc[0], 10**12))
    rep_ids, rep_counts = [r for r, _ in rep_ids_counts], [c for _, c in rep_ids_counts]

    # Also ensure reps have a known path (sanity)
    unknown = [rid for rid in rep_ids if rid not in id_to_path]
    if unknown:
        raise RuntimeError(
            "Selected representative ids not present in the provided items: "
            + ", ".join(map(str, unknown[:10]))
        )

    logger.info(f"Selected {len(rep_ids)} representative ids")

    return rep_ids, rep_counts


def group_members(linked_groups: Linker) -> set[str]:
    members: set[str] = set()
    for g in linked_groups.groups:
        members.update(map(str, g.members))
    return members


def iterative_select_representatives(
    *,
    items: list[Item],
    cluster_threshold: float,
    distance_threshold: float,
    scaling_factor: float,
    max_reps: int,
    batch_size: int,
    mash_path: str,
    tmp_dir: Path | None = None,
    search_threads: int = 1,
) -> tuple[list[str], list[int]]:
    """
    A script for extracting a set of representatives from a collection of FASTAs, using the same two-phase selection
     used by the standard library builder:

    Phase A: strict initial clustering at `cluster_threshold`.
    Phase B: if still > max_reps, iteratively relax threshold starting at
             min(distance_threshold * scaling_factor, distance_threshold),
             doubling each time (capped at distance_threshold),
             stopping once "threshold" hits the cap (warn if still > max_reps).

    Returns:
      (rep_ids, rep_counts)
    """
    if not items:
        raise ValueError("No items provided.")
    if max_reps < 1:
        raise ValueError("max_reps must be >= 1.")
    if distance_threshold <= 0:
        raise ValueError("distance_threshold must be > 0.")
    if scaling_factor <= 0:
        raise ValueError("scaling_factor must be > 0.")

    # Phase A: strict pass (only meaningful if there is more than 1 item)
    if len(items) == 1:
        rep_id, _ = items[0]
        return [rep_id], [1]

    rep_ids, rep_counts = extract_exemplars_from_items(
        items=items,
        threshold=cluster_threshold,
        batch_size=batch_size,
        mash_path=mash_path,
        previous_counts=None,
        tmp_dir=tmp_dir,
        search_threads=search_threads,
    )


    current_items = [(rid, dict(items)[rid]) for rid in rep_ids]

    # Phase B: iterative relaxation if needed
    if len(current_items) > max_reps:
        threshold = min(distance_threshold * scaling_factor, distance_threshold)
        if threshold <= 0:
            threshold = distance_threshold

        while len(current_items) > max_reps and threshold <= distance_threshold:
            logger.info(
                "Selecting representatives: %d items, threshold=%g (max=%g, target max_reps=%d)",
                len(current_items),
                threshold,
                distance_threshold,
                max_reps,
            )
            rep_ids, rep_counts = extract_exemplars_from_items(
                items=current_items,
                threshold=threshold,
                batch_size=batch_size,
                mash_path=mash_path,
                previous_counts=rep_counts,
                tmp_dir=tmp_dir,
                search_threads=search_threads,
            )
            current_lookup = dict(current_items)
            current_items = [(rid, current_lookup[rid]) for rid in rep_ids]

            if threshold == distance_threshold:
                break
            threshold = min(threshold * 2, distance_threshold)

        if len(current_items) > max_reps and threshold == distance_threshold:
            logger.warning(
                "Reached distance_threshold=%g but still have %d representatives (max_reps=%d). "
                "Proceeding with current representatives.",
                distance_threshold,
                len(current_items),
                max_reps,
            )

    # Final reps
    final_ids = [rid for rid, _ in current_items]
    final_counts = rep_counts
    logger.info(f"Selected {len(final_ids)} representatives")
    return final_ids, final_counts
