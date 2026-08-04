import numpy as np
import pytest
from astropy import constants as const

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def test_beta_rad(t_rad, beta_rad):
    assert np.allclose(beta_rad, 1 / (const.k_B.cgs.value * t_rad))


def test_g_electron(g_electron):
    assert np.allclose(g_electron, 2.4146828342691716e21)


def test_number_density(number_density):
    assert np.isclose(number_density[0].loc[2], 1504556808.6958313)


def test_selected_atoms(selected_atoms):
    assert selected_atoms == [2]


def test_electron_temperature(t_rad, link_t_rad_t_electron, t_electrons):
    assert np.allclose(t_electrons, t_rad * link_t_rad_t_electron)


def test_beta_electron(beta_electron, t_electrons):
    assert np.allclose(beta_electron, 1 / (const.k_B.cgs.value * t_electrons))
