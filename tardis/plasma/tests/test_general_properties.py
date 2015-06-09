import numpy as np
from astropy import constants as const

def test_beta_rad(t_rad, beta_rad):
    assert np.allclose(beta_rad, 1 / (const.k_B.cgs.value * t_rad))

def test_g_electron(g_electron):
    assert np.allclose(g_electron, 2.4146828342691716e+21)

def test_number_density(number_density):
    assert np.isclose(number_density[0].loc[2], 1504556808.6958313)