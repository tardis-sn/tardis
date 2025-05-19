import pytest
from pathlib import Path
import numpy as np
import shutil

from tardis.io.model.readers.snec.read_snec_output import read_snec_output_xg
from tardis.io.model.readers.snec.xg_files import read_xg_file


def test_read_xg_components(minimal_xg_file):
    xg_data = read_xg_file(minimal_xg_file, column_names=["mass", "radius", "temperature"])
    assert len(xg_data.timestamps) == 2
    assert len(xg_data.data_blocks) == 2    
    assert np.isclose(xg_data.data_blocks[0].iloc[0]['radius'], 1.0)
    assert np.isclose(xg_data.data_blocks[0].iloc[0]['temperature'], 100.0)
