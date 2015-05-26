import pytest

import numpy as np

from tardis.plasma.properties.ion_population import (PhiSahaLTE,
IonNumberDensity)
from tardis.plasma.properties.general import (BetaRadiation, GElectron,
NumberDensity)
from tardis.plasma.properties.partition_function import (LevelBoltzmannFactor,
LTEPartitionFunction)
from tardis.plasma.properties.atomic import Levels, AtomicMass, IonizationData

def test_phi_saha_lte(t_rad, beta_rad, g_electron, ionization_data,
        phi_saha_lte):
    assert(phi_saha_lte.shape == (2,20))
    assert np.all(np.isclose(g_electron * 4 * np.exp(
        -ionization_data.ionization_energy.ix[2].ix[1] * beta_rad),
        phi_saha_lte.ix[2].ix[1]))

def test_ion_number_density(phi_saha_lte, ion_number_density):
    ion_numbers = ion_number_density.index.get_level_values(1).values
    ion_numbers = ion_numbers.reshape((ion_numbers.shape[0], 1))
    n_electron_from_ions = (ion_number_density.values * ion_numbers).sum(
        axis=0)
    n_electron_from_saha = (ion_number_density.ix[2].ix[0]/
        ion_number_density.ix[2].ix[1]) * phi_saha_lte.ix[2].ix[1]
    tolerance = 0.01 * n_electron_from_ions
    assert np.allclose(n_electron_from_ions, n_electron_from_saha,
        atol=tolerance)