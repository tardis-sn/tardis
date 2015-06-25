import os

import numpy as np
import pandas as pd
from astropy import units as u
import pytest

import tardis
from tardis.atomic import AtomData

from tardis.plasma.properties.ion_population import (PhiGeneral,
PhiSahaLTE, IonNumberDensity, ElectronDensity, PhiSahaNebular,
RadiationFieldCorrection)
from tardis.plasma.properties.general import (BetaRadiation, GElectron,
NumberDensity, ElectronTemperature, BetaElectron)
from tardis.plasma.properties.partition_function import (
LevelBoltzmannFactorLTE, PartitionFunction, LevelBoltzmannFactorDiluteLTE)
from tardis.plasma.properties.atomic import (Levels, Lines, AtomicMass,
IonizationData, LinesUpperLevelIndex, LinesLowerLevelIndex, ZetaData)
from tardis.plasma.standard_plasmas import LTEPlasma
from tardis.plasma.properties.level_population import (LevelPopulation,
LevelNumberDensity)
from tardis.plasma.properties.radiative_properties import (TauSobolev,
StimulatedEmissionFactor, BetaSobolev, TransitionProbabilities)

@pytest.fixture
def atomic_data():
    atomic_db_fname = os.path.join(tardis.__path__[0], 'tests', 'data',
                                   'chianti_he_db.h5')
    return AtomData.from_hdf5(atomic_db_fname)

@pytest.fixture
def number_of_cells():
    return 20

@pytest.fixture
def selected_atoms():
    return [2]

@pytest.fixture
def abundance(number_of_cells, selected_atoms):
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
def beta_rad(t_rad):
    beta_rad_module = BetaRadiation(None)
    return beta_rad_module.calculate(t_rad)

@pytest.fixture
def g_electron(beta_rad):
    g_electron_module = GElectron(None)
    return g_electron_module.calculate(beta_rad)

@pytest.fixture
def levels(atomic_data, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_data, selected_atoms)

@pytest.fixture
def lines(atomic_data, selected_atoms):
    lines_module = Lines(None)
    return lines_module.calculate(atomic_data, selected_atoms)

@pytest.fixture
def j_blues(lines):
    return pd.DataFrame(1.e-5, index=lines.index, columns=range(20))

@pytest.fixture
def link_t_rad_t_electron():
    return 0.9

@pytest.fixture
def level_boltzmann_factor_lte(levels, beta_rad):
    level_boltzmann_factor_module = LevelBoltzmannFactorLTE(None)
    return level_boltzmann_factor_module.calculate(levels, beta_rad)

@pytest.fixture
def level_boltzmann_factor_dilute_lte(levels, beta_rad, w):
    level_boltzmann_factor_module = LevelBoltzmannFactorDiluteLTE(None)
    return level_boltzmann_factor_module.calculate(levels, beta_rad, w)

@pytest.fixture
def partition_function(levels, level_boltzmann_factor_lte):
    partition_function_module = PartitionFunction(None)
    return partition_function_module.calculate(levels,
        level_boltzmann_factor_lte)

@pytest.fixture
def ionization_data(atomic_data, selected_atoms):
    ionization_data_module = IonizationData(None)
    return ionization_data_module.calculate(atomic_data,
        selected_atoms)

@pytest.fixture
def general_phi(g_electron, beta_rad, partition_function, ionization_data):
    phi_general_module = PhiGeneral(None)
    return phi_general_module.calculate(g_electron, beta_rad,
        partition_function, ionization_data)

@pytest.fixture
def phi_saha_lte(general_phi):
    phi_module = PhiSahaLTE(None)
    return phi_module.calculate(general_phi)

@pytest.fixture
def atomic_mass(atomic_data, selected_atoms):
    atomic_mass_module = AtomicMass(None)
    return atomic_mass_module.calculate(atomic_data,
        selected_atoms)

@pytest.fixture
def zeta_data(atomic_data, selected_atoms):
    zeta_data_module = ZetaData(None)
    return zeta_data_module.calculate(atomic_data, selected_atoms)

@pytest.fixture
def number_density(atomic_mass, abundance, density):
    number_density_module = NumberDensity(None)
    return number_density_module.calculate(atomic_mass, abundance, density)

@pytest.fixture
def ion_number_density(phi_saha_lte, partition_function,
    number_density):
    ion_number_density_module = IonNumberDensity(None)
    return ion_number_density_module.calculate(phi_saha_lte,
        partition_function, number_density)

@pytest.fixture
def electron_densities(ion_number_density, phi_saha_lte):
    electron_density_module = ElectronDensity(None)
    return electron_density_module.calculate(ion_number_density,
        phi_saha_lte)

@pytest.fixture
def level_population(levels, partition_function,
    level_boltzmann_factor_lte):
    level_population_module = LevelPopulation(None)
    return level_population_module.calculate(levels,
        partition_function, level_boltzmann_factor_lte)

@pytest.fixture
def level_number_density(level_population, ion_number_density):
    level_number_density_module = LevelNumberDensity(None)
    return level_number_density_module.calculate(level_population,
    ion_number_density)

@pytest.fixture
def lines_upper_level_index(lines, levels):
    upper_level_index_module = LinesUpperLevelIndex(None)
    return upper_level_index_module.calculate(levels, lines)

@pytest.fixture
def lines_lower_level_index(lines, levels):
    lower_level_index_module = LinesLowerLevelIndex(None)
    return lower_level_index_module.calculate(levels, lines)

@pytest.fixture
def stimulated_emission_factor(levels, level_number_density,
    lines_lower_level_index, lines_upper_level_index):
    factor_module = StimulatedEmissionFactor(None)
    return factor_module.calculate(levels, level_number_density,
        lines_lower_level_index,lines_upper_level_index)

@pytest.fixture
def tau_sobolev(lines, level_number_density, lines_lower_level_index,
    time_explosion, stimulated_emission_factor, j_blues):
    tau_sobolev_module = TauSobolev(None)
    return tau_sobolev_module.calculate(lines, level_number_density,
        lines_lower_level_index, time_explosion, stimulated_emission_factor,
        j_blues)

@pytest.fixture
def t_electron(t_rad, link_t_rad_t_electron):
    electron_temperature_module = ElectronTemperature(None)
    return electron_temperature_module.calculate(t_rad, link_t_rad_t_electron)

@pytest.fixture
def beta_electron(t_electron):
    beta_electron_module = BetaElectron(None)
    return beta_electron_module.calculate(t_electron)

@pytest.fixture
def electron_densities(ion_number_density, phi_saha_lte):
    electron_densities_module = ElectronDensity(None)
    return electron_densities_module.calculate(ion_number_density,
        phi_saha_lte)

@pytest.fixture
def delta(w, ionization_data, beta_rad, t_electron, t_rad, beta_electron,
    levels):
    delta_input = None
    delta_module = RadiationFieldCorrection(None)
    delta_module.chi_0_species = (2,2)
    return delta_module.calculate(w, ionization_data, beta_rad, t_electron,
        t_rad, beta_electron, levels, delta_input)

@pytest.fixture
def phi_saha_nebular(general_phi, t_rad, w, zeta_data, t_electron, delta):
    phi_saha_nebular_module = PhiSahaNebular(None)
    return phi_saha_nebular_module.calculate(general_phi, t_rad, w, zeta_data,
                                   t_electron, delta)

@pytest.fixture
def beta_sobolev(tau_sobolev):
    beta_sobolev_module = BetaSobolev(None)
    return beta_sobolev_module.calculate(tau_sobolev)

@pytest.fixture
def transition_probabilities(atomic_data, beta_sobolev, j_blues,
    stimulated_emission_factor, tau_sobolev):
    transition_probabilities_module = TransitionProbabilities(None)
    return transition_probabilities_module.calculate(atomic_data, beta_sobolev,
        j_blues, stimulated_emission_factor, tau_sobolev)

@pytest.fixture
def standard_lte_plasma_he_db(t_rad, abundance, density, time_explosion,
                              atomic_data, j_blues,
                              link_t_rad_t_electron):
    return LTEPlasma(t_rad, abundance, density, time_explosion,
                     atomic_data, j_blues, link_t_rad_t_electron)