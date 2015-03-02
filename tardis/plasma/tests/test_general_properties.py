import numpy as np
from astropy import constants as const

from tardis.plasma.properties.general import (BetaRadiation, GElectron,
                                              NumberDensity, SelectedAtoms)



def test_beta_rad1(t_rad):
    beta_rad_module = BetaRadiation(None)
    assert beta_rad_module.calculate(5000) < beta_rad_module.calculate(2000)
    assert np.allclose(beta_rad_module.calculate(t_rad), 1 /
                       (const.k_B.cgs.value * t_rad))

def test_g_electron1(t_rad):
    g_electron_module = GElectron(None)
    g_electron = g_electron_module.calculate(1 / (const.k_B.cgs.value * 10000.))
    assert np.isclose(g_electron, 2.4146828342691716e+21)

def test_number_density(included_he_atomic_data, abundance, density):
    nb_density_module = NumberDensity(None)
    included_he_atomic_data.atom_data.mass

    he_nb_density = nb_density_module.calculate(
        included_he_atomic_data.atom_data.mass, abundance, density)

    assert np.isclose(he_nb_density[0].loc[2], 1504556808.6958313)

def test_number_density(abundance):
    selected_atoms_module = SelectedAtoms(None)
    selected_atoms = selected_atoms_module.calculate(abundance)
    assert np.isclose(selected_atoms, abundance.index)

