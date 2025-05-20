import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from astropy import units as u

from tardis.io.model.readers.snec.xg_files import xg_block_size, read_xg_file, XGData

@pytest.fixture
def regression_test_xg_file(regression_data):
    return regression_data.regression_data_path / "MESA_STIR_MESA_SNEC" / "output" / "cs2.xg"


def test_xg_block_size(regression_test_xg_file):
    num_blocks, block_size = xg_block_size(regression_test_xg_file)
    
    assert num_blocks == 1001
    assert block_size == np.int64(1141)


def test_read_xg_file(regression_test_xg_file):
    xg_data = read_xg_file(regression_test_xg_file, column_names=["Time"])
    
    assert isinstance(xg_data, XGData)
    assert len(xg_data.timestamps) == 1001
    assert xg_data.timestamps.unit == u.s
    assert len(xg_data.data_blocks) == 1001
    
    assert isinstance(xg_data.data_blocks[0], pd.DataFrame)
    assert list(xg_data.data_blocks[0].columns) == ["Time"]
    
    assert np.isclose(xg_data.data_blocks[1].iloc[1, 0], np.float64(425927938160107.2))


def test_xgdata_to_xarray(regression_test_xg_file):
    xg_data = read_xg_file(regression_test_xg_file, column_names=["Time"])
    
    xarray_data = xg_data.to_xarray()
    
    assert xarray_data.dims == ('time', 'cell_id', 'quantity')
    assert list(xarray_data.coords['quantity'].values) == ["Time"]
    assert len(xarray_data.coords['time'].values) == 1001
