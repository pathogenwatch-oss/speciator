from __future__ import annotations

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from speciator.search import HEAVY_WEIGHT, add_weights, determine_threshold


def test_determine_threshold_no_adjustment_needed():
    """Test when all species have sufficient representatives (>= threshold)"""
    # Create test data where all species have 10+ representatives
    sorted_matches = pl.DataFrame(
        {
            "distance": [0.01, 0.02, 0.03, 0.04, 0.05],
            "weighted_len": [15, 12, 20, 18, 25],  # All species have >= 10 representatives
        }
    )

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 10, "Should return max_threshold when no adjustment needed"


def test_determine_threshold_species_count_less_than_hit_index():
    """Test when species count <= hit index, threshold should be hit_index + 1"""
    # Species at index 2 (3rd hit) has only 2 representatives, which is <= 3
    sorted_matches = pl.DataFrame(
        {
            "distance": [0.01, 0.02, 0.03, 0.04, 0.05],
            "weighted_len": [15, 12, 2, 18, 25],  # 3rd species has only 2 representatives
        }
    )

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 3, (
        "Should return hit_index + 1 when species_count <= hit_index + 1"
    )


def test_determine_threshold_species_count_greater_than_hit_index():
    """Test when species count > hit index, threshold should be species count"""
    # Species at index 1 (2nd hit) has 5 representatives, which is > 2
    sorted_matches = pl.DataFrame(
        {
            "distance": [0.01, 0.02, 0.03, 0.04, 0.05],
            "weighted_len": [15, 5, 20, 18, 25],  # 2nd species has 5 representatives
        }
    )

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 5, "Should return species_count when species_count > hit_index + 1"


def test_determine_threshold_first_hit_triggers_adjustment():
    """Test when first hit (index 0) has low species count"""
    # First species has only 1 representative
    sorted_matches = pl.DataFrame(
        {
            "distance": [0.01, 0.02, 0.03, 0.04, 0.05],
            "weighted_len": [1, 12, 20, 18, 25],  # 1st species has only 1 representative
        }
    )

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 1, (
        "Should return hit_index + 1 = 1 when first hit has species_count = 1"
    )


def test_determine_threshold_empty_matches():
    """Test with empty DataFrame"""
    sorted_matches = pl.DataFrame({"distance": [], "len": []})

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 10, "Should return max_threshold when no matches"


def test_determine_threshold_fewer_matches_than_threshold():
    """Test when we have fewer matches than max_threshold"""
    sorted_matches = pl.DataFrame(
        {
            "distance": [0.01, 0.02, 0.03],
            "weighted_len": [15, 12, 20],  # Only 3 matches, all with sufficient representatives
        }
    )

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 10, (
        "Should return max_threshold when all matches have sufficient representatives"
    )


def test_determine_threshold_edge_case_species_count_equals_hit_index():
    """Test edge case where species_count exactly equals hit_index + 1"""
    # Species at index 2 has exactly 3 representatives (hit_index + 1 = 3)
    sorted_matches = pl.DataFrame(
        {
            "distance": [0.01, 0.02, 0.03, 0.04, 0.05],
            "weighted_len": [15, 12, 3, 18, 25],  # 3rd species has exactly 3 representatives
        }
    )

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 3, (
        "Should return hit_index + 1 when species_count == hit_index + 1"
    )


def test_determine_threshold_beyond_threshold():
    """If a match sets the threshold to above its position in the sorted matches."""
    sorted_matches = pl.DataFrame(
        {
            "distance": [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07],
            "weighted_len": [
                15,
                15,
                15,
                15,
                15,
                15,
                2,
            ],  # All species have >= 10 representatives
        }
    )
    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 7, (
        "Should return position of match that sets threshold to above its position in sorted matches"
    )


def test_determine_threshold_vibrio_cholerae_scenario():
    # Simulate scenario where A. rare (few representatives) appears early
    # followed by A. common (many representatives)
    sorted_matches = pl.DataFrame(
        {
            "distance": [
                0.01,
                0.015,
                0.02,
                0.025,
                0.03,
                0.035,
                0.04,
                0.045,
                0.05,
                0.055,
            ],
            "weighted_len": [
                2,
                2,
                50,
                50,
                50,
                50,
                50,
                50,
                50,
                50,
            ],  # First 2 hits are A. rare (2 reps)
        }
    )

    result = determine_threshold(sorted_matches, max_threshold=10)
    assert result == 2, "Should return 2"


def test_determine_threshold_max_threshold_one():
    """Test with max_threshold = 1"""
    sorted_matches = pl.DataFrame({"distance": [0.01, 0.02, 0.03], "len": [5, 12, 20]})

    result = determine_threshold(sorted_matches, max_threshold=1)
    assert result == 1, "Should handle max_threshold = 1 correctly"


@pytest.fixture
def sample_matches_df() -> pl.DataFrame:
    """Fixture for a sample DataFrame to test add_weights."""
    return pl.DataFrame(
        {
            "SpeciesName": [
                "Salmonella enterica",
                "Bacterium eubacterium",
                "Hyphomicrobiales bacterium",
                "Rhizobium sp. GHKF11",
                "Staphylococcus aureus",
            ],
            "len": [10, 20, 30, 40, 50],
            "distance": [0.01, 0.02, 0.03, 0.04, 0.05],
        }
    )


def test_add_weights(sample_matches_df: pl.DataFrame):
    """
    Test the add_weights function to ensure it correctly assigns weights.

    - Standard species names should have `weighted_len` equal to `len`.
    - Generic species names containing ' bacterium' or ' sp.' should be assigned HEAVY_WEIGHT.
    """
    # Expected output DataFrame
    expected_df = sample_matches_df.with_columns(
        pl.Series(
            "weighted_len",
            [10, 20, HEAVY_WEIGHT, HEAVY_WEIGHT, 50],
            dtype=pl.Int64,
        )
    )

    # Apply the function to the test DataFrame
    result_df = add_weights(sample_matches_df)

    # Assert that the resulting DataFrame is equal to the expected one
    assert_frame_equal(result_df, expected_df)


def test_add_weights_complex_real():
    """Test add_weights with a more complex, realistic dataset."""
    csv = """filename,distance,p-value,shared-hashes
"GCF_013267135.1",0.0,0,"1000/1000"
"GCF_013267155.1",0.000582,0,"976/1000"
"GCF_003031575.1",0.0110882,0,"656/1000"
"GCF_029632585.1",0.0164989,0,"547/1000"
"GCA_001644925.1",0.0177122,0,"526/1000"
"GCA_900466835.1",0.0177716,0,"525/1000"
"GCA_900466415.1",0.0178312,0,"524/1000"
"GCA_900470725.1",0.0178312,0,"524/1000"
"GCA_900467925.1",0.0179508,0,"522/1000"
"GCA_900465765.1",0.0180108,0,"521/1000"
"""
    matches_df = pl.read_csv(
        source=csv.encode("utf-8"), has_header=True, quote_char='"'
    ).with_columns(
        pl.Series(
            "SpeciesName",
            [
                "Agrobacterium pusense",
                "Agrobacterium pusense",
                "Agrobacterium pusense",
                "Agrobacterium pusense",
                "Rhizobium sp. GHKF11",
                "Hyphomicrobiales bacterium",
                "Hyphomicrobiales bacterium",
                "Hyphomicrobiales bacterium",
                "Hyphomicrobiales bacterium",
                "Hyphomicrobiales bacterium",
            ],
        ),
        pl.Series("len", [1, 1, 1, 1, 1, 1, 1, 1, 1, 2]),
    )

    expected_df = matches_df.sort("distance", descending=False).with_columns(
        pl.Series(
            "weighted_len",
            [
                1,
                1,
                1,
                1,
                HEAVY_WEIGHT,
                HEAVY_WEIGHT,
                HEAVY_WEIGHT,
                HEAVY_WEIGHT,
                HEAVY_WEIGHT,
                HEAVY_WEIGHT,
            ],
            dtype=pl.Int64,
        )
    )

    result_df = add_weights(matches_df)

    assert_frame_equal(result_df, expected_df)


def test_add_weights_empty_dataframe():
    """Test that add_weights handles an empty DataFrame correctly."""
    empty_df = pl.DataFrame(
        {"SpeciesName": [], "len": [], "distance": 0},
        schema={"SpeciesName": pl.Utf8, "len": pl.Int64, "distance": pl.Float64},
    )
    expected_df = empty_df.with_columns(pl.Series("weighted_len", [], dtype=pl.Int64))

    result_df = add_weights(empty_df)

    assert_frame_equal(result_df, expected_df)
