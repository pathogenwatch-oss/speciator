
import os
import tempfile
from pathlib import Path

import polars as pl

from speciator.library_stats import LibraryStats


def test_extract_stats():
    # Create test data
    test_data = pl.DataFrame({
        "organism_code": ["123", "456", "123", "789", "456", "456"],
        "other_column": ["a", "b", "c", "d", "e", "f"]
    })
    
    # Create temporary parquet file
    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp_file:
        tmp_path = Path(tmp_file.name)
        test_data.write_parquet(tmp_path)
    
    try:
        # Test the extract_stats method
        result = LibraryStats.extract_stats(tmp_path)
        
        # Assertions
        assert result.name == tmp_path.stem
        assert isinstance(result.counts, dict)
        assert result.counts["123"] == 2
        assert result.counts["456"] == 3
        assert result.counts["789"] == 1
        assert len(result.counts) == 3
        
    finally:
        # Clean up temporary file
        os.unlink(tmp_path)


def test_extract_stats_empty_file():
    # Create empty test data
    test_data = pl.DataFrame({"organism_code": []})
    
    with tempfile.NamedTemporaryFile(suffix=".pqt", delete=False) as tmp_file:
        tmp_path = Path(tmp_file.name)
        test_data.write_parquet(tmp_path)
    
    try:
        result = LibraryStats.extract_stats(tmp_path)
        
        assert result.name == tmp_path.stem
        assert result.counts == {}
        
    finally:
        os.unlink(tmp_path)
