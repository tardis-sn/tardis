import tardis
from tardis.io.model_reader import read_artis_density

from astropy import units as u
import numpy as np

import pytest

import os

data_path = os.path.join(tardis.__path__[0], 'io', 'tests', 'data')

@pytest.fixture
def artis_density_fname():
    return os.path.join(data_path, 'artis_model.dat')

def test_simple_read_artis_density(artis_density_fname):
    (time_of_model, index, v_inner, v_outer,
     mean_density) = read_artis_density(artis_density_fname)

    assert np.isclose(0.00114661 * u.day, time_of_model, atol=1e-7 * u.day)
    assert len(mean_density) == 69
    assert len(mean_density) == len(v_inner)


