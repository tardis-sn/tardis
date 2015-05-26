import numpy as np

from tardis.plasma.properties.ion_population import (PhiSahaLTE,
IonNumberDensity)
from tardis.plasma.properties.general import (BetaRadiation, GElectron,
NumberDensity)
from tardis.plasma.properties.partition_function import (LevelBoltzmannFactor,
LTEPartitionFunction)
from tardis.plasma.properties.atomic import Levels, AtomicMass, IonizationData

def test_phi_saha_lte(t_rad, included_he_atomic_data):
    beta_rad_module = BetaRadiation(None)
    beta_rad = beta_rad_module.calculate(t_rad)
    g_electron_module = GElectron(None)
    g_electron = g_electron_module.calculate(beta_rad)
    levels_module = Levels(None)
    selected_atoms = [2]
    levels = levels_module.calculate(included_he_atomic_data, selected_atoms)
    level_boltzmann_factor_module = LevelBoltzmannFactor(None)
    level_boltzmann_factor = level_boltzmann_factor_module.calculate(levels,
        beta_rad)
    partition_function_module = LTEPartitionFunction(None)
    partition_function = partition_function_module.calculate(levels,
        level_boltzmann_factor)
    ionization_data_module = IonizationData(None)
    ionization_data = ionization_data_module.calculate(included_he_atomic_data,
        selected_atoms)
    phi_module = PhiSahaLTE(None)
    phi = phi_module.calculate(g_electron, beta_rad, partition_function,
        ionization_data)
    assert(phi.shape == (2,20))
    assert np.all(np.isclose(g_electron * 4 * np.exp(
        -ionization_data.ionization_energy.ix[2].ix[1] * beta_rad),
        phi.ix[2].ix[1]))

def test_ion_number_density(t_rad, abundance, density,
        included_he_atomic_data):
    beta_rad_module = BetaRadiation(None)
    beta_rad = beta_rad_module.calculate(t_rad)
    g_electron_module = GElectron(None)
    g_electron = g_electron_module.calculate(beta_rad)
    selected_atoms = [2]
    atomic_mass_module = AtomicMass(None)
    atomic_mass = atomic_mass_module.calculate(included_he_atomic_data,
        selected_atoms)
    number_density_module = NumberDensity(None)
    number_density = number_density_module.calculate(atomic_mass,
        abundance, density)
    levels_module = Levels(None)
    selected_atoms = [2]
    levels = levels_module.calculate(included_he_atomic_data, selected_atoms)
    level_boltzmann_factor_module = LevelBoltzmannFactor(None)
    level_boltzmann_factor = level_boltzmann_factor_module.calculate(levels,
        beta_rad)
    partition_function_module = LTEPartitionFunction(None)
    partition_function = partition_function_module.calculate(levels,
        level_boltzmann_factor)
    ionization_data_module = IonizationData(None)
    ionization_data = ionization_data_module.calculate(included_he_atomic_data,
        selected_atoms)
    phi_module = PhiSahaLTE(None)
    phi = phi_module.calculate(g_electron, beta_rad, partition_function,
        ionization_data)
    ion_number_density_module = IonNumberDensity(None)
    ion_number_density = ion_number_density_module.calculate(phi,
        partition_function, number_density)
    ion_numbers = ion_number_density.index.get_level_values(1).values
    ion_numbers = ion_numbers.reshape((ion_numbers.shape[0], 1))
    n_electron_from_ions = (ion_number_density.values * ion_numbers).sum(
        axis=0)
    n_electron_from_saha = (ion_number_density.ix[2].ix[0]/
        ion_number_density.ix[2].ix[1]) * phi.ix[2].ix[1]
    tolerance = 0.01 * n_electron_from_ions
    assert np.allclose(n_electron_from_ions, n_electron_from_saha,
        atol=tolerance)