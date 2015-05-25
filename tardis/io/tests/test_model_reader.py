import tardis
from tardis.io.model_reader import (read_artis_density,
read_simple_ascii_abundances)

from astropy import units as u
import numpy as np

import pytest

import os

data_path = os.path.join(tardis.__path__[0], 'io', 'tests', 'data')

@pytest.fixture
def artis_density_fname():
    return os.path.join(data_path, 'artis_model.dat')
@pytest.fixture
def artis_abundances_fname():
    return os.path.join(data_path, 'artis_abundances.dat')

def test_simple_read_artis_density(artis_density_fname):
    (time_of_model, index, v_inner, v_outer,
     mean_density) = read_artis_density(artis_density_fname)

    assert np.isclose(0.00114661 * u.day, time_of_model, atol=1e-7 * u.day)
    assert np.isclose(mean_density[23], 0.2250048 * u.g / u.cm**3, atol=1.e-6 * 
        u.g / u.cm**3)
    assert len(mean_density) == 69
    assert len(mean_density) == len(v_inner)

#Artis files are currently read with read ascii files function
def test_read_simple_ascii_abundances(artis_abundances_fname):
    index, abundances = read_simple_ascii_abundances(artis_abundances_fname)
    assert len(abundances.columns) == 69
    assert np.isclose(abundances[23].ix[2], 2.672351e-08 , atol=1.e-12)
