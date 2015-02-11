import pytest
import numpy as np

import pandas as pd

from astropy import units as u

from tardis.plasma.standard_plasmas import LTEPlasma


@pytest.fixture
def number_of_cells():
    return 20

@pytest.fixture
def abundance(number_of_cells):
    selected_atoms = [2]
    abundances = [1.0]
    return pd.DataFrame(data=1.0, index=selected_atoms,
                        columns=range(20), dtype=np.float64)

@pytest.fixture
def density(number_of_cells):
    return np.ones(number_of_cells) * 1e-14

@pytest.fixture
def w(number_of_cells):
    return np.ones(number_of_cells) * 0.5


@pytest.fixture
def time_explosion():
    return (19 * u.day).to(u.s).value

@pytest.fixture
def t_rad(number_of_cells):
    return np.ones(number_of_cells) * 10000


@pytest.fixture
def standard_lte_plasma_he_db(t_rad, abundance, density, time_explosion,
                              included_he_atomic_data):
    return LTEPlasma(t_rad, abundance, density, time_explosion,
                     included_he_atomic_data)