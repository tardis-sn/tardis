import os
import numpy as np
import pytest
from pathlib import Path
from astropy import units as u

from tardis.io.model.readers.snec.base import xg_block_size, read_xg_timestamps_and_data

def test_xg_block_size(minimal_xg_file):
    num_blocks, block_size = xg_block_size(minimal_xg_file)
    
    assert num_blocks == 2
    assert block_size == 4


def test_read_xg_timestamps_and_data(minimal_xg_file):
    column_names = ["mass", "radius", "temperature"]
    timestamps, data_blocks = read_xg_timestamps_and_data(minimal_xg_file, column_names)
    
    assert len(timestamps) == 2
    assert timestamps.unit == u.s
    assert timestamps[0].value == 0.0
    assert timestamps[1].value == 1.0
    
    assert len(data_blocks) == 2
    assert data_blocks[0].shape == (3, 3)
    assert data_blocks[1].shape == (3, 3)
    
    # Basic check for data content
    assert np.isclose(data_blocks[0].iloc[0, 0], 0.1)
    assert np.isclose(data_blocks[1].iloc[0, 1], 1.1)


def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        read_xg_timestamps_and_data("non_existent_file.xg", ["column1"])
