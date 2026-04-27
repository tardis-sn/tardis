import gzip

import pandas as pd
import pytest

# Import your actual class here
from tardis.io.continuum_reader.s92_reader import OpacityReader


@pytest.fixture
def mock_data_path(tmp_path):
    """Programmatically creates a synthetic dataset for testing."""
    # We add the .gz extension so the reader knows it's a zip
    path = tmp_path / "test_data.gz"
    content = (
        "Header Line 1\nHeader Line 2\n"  # 2 lines of junk
        "10 1.1 2.2 0.03\n"  # Good Row
        "20 1.2 2.3 0.04\n"  # Good Row
        "150 9.9 9.9 0.99\n"  # The 'Spike' (Index 150)
        "30 1.3 2.4 0.05\n"  # Post-Spike Good Row
    )
    with gzip.open(path, "wt") as f:
        f.write(content)
    return str(path)  # Return string path for the reader


class TestS92Reader:
    """Test class to ensure the module works under stress and with datasets."""

    def test_file_not_found(self):
        """Ensures the reader raises FileNotFoundError for missing files."""
        reader = OpacityReader("non_existent_file.gz")
        with pytest.raises(FileNotFoundError):
            reader.read()

    def test_malformed_columns_fallback(self, tmp_path):
        """Verifies that malformed columns return a RAW table instead of crashing. It is to showcase the resilience and
        stability of the program to not crash in the face of unknown data formay."""
        path = tmp_path / "broken.gz"
        content = "Header\n10 1.1 2.2 \n20 1.2 2.3\n"
        with gzip.open(path, "wt") as f:
            f.write(content)
        reader = OpacityReader(str(path), skip_rows=1)
        df = reader.read()
        # PROVE the user gets a Raw Table instead of a crash
        assert isinstance(df, pd.DataFrame)
        # Verify it's the raw table (doesn't have the named 'index' column)
        assert "index" not in df.columns
    def test_reader_parsing_and_masking(self, mock_data_path):
        """File to ensure no warning and to showcase that the test works as expected."""
        reader = OpacityReader(mock_data_path, skip_rows=2, masking_value = 120)
        data = reader.read()
        assert isinstance(data, pd.DataFrame)
        assert not data.empty


    def test_skippers(self, mock_data_path):
        reader = OpacityReader(mock_data_path, skip_rows =2)
        df = reader.read()
        assert len(df) == 4
        assert list(df.columns) == ["index", "log_t", "log_r", "opacity"]
        assert df.iloc[0]["index"] == 10
        assert df['index'].dtype == 'int64'
    def test_maskordered(self, mock_data_path):
        reader = OpacityReader(mock_data_path, skip_rows=2, masking_value = 120)
        df = reader.read()
        assert len(df) == 3
        assert 150 not in df['index'].values
