import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from astropy import units as u

from tardis.io.model.readers.snec.xg_files import xg_block_size, read_xg_file, XGData

def test_xg_block_size(minimal_xg_file):
    num_blocks, block_size = xg_block_size(minimal_xg_file)
    
    assert num_blocks == 2
    assert block_size == 4

def test_read_xg_file(minimal_xg_file):
    column_names = ["mass", "radius", "temperature"]
    
    xg_data = read_xg_file(minimal_xg_file, column_names)
    
    assert isinstance(xg_data, XGData)
    assert len(xg_data.timestamps) == 2
    assert xg_data.timestamps.unit == u.s
    assert len(xg_data.data_blocks) == 2
    
    # Check one data block
    assert isinstance(xg_data.data_blocks[0], pd.DataFrame)
    assert list(xg_data.data_blocks[0].columns) == column_names
    
    # Check a value
    assert np.isclose(xg_data.data_blocks[1].iloc[1, 2], 210.0)


def test_xgdata_to_xarray(minimal_xg_file):
    column_names = ["mass", "radius", "temperature"]
    xg_data = read_xg_file(minimal_xg_file, column_names)
    
    xarray_data = xg_data.to_xarray()
    
    assert xarray_data.dims == ('time', 'cell_id', 'quantity')
    assert list(xarray_data.coords['quantity'].values) == column_names
    assert len(xarray_data.coords['time'].values) == 2
