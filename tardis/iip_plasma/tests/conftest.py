import os

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

import tardis
from tardis.iip_plasma.properties import *
from tardis.io.atom_data import AtomData

# INPUTS


@pytest.fixture
def atomic_data():
    atomic_db_fname = os.path.join(
        tardis.__path__[0], "tests", "data", "chianti_he_db.h5"
    )
    return AtomData.from_hdf(atomic_db_fname)


@pytest.fixture
def number_of_cells():
    return 20


@pytest.fixture
def abundance(number_of_cells):
    return pd.DataFrame(
        data=1.0, index=[2], columns=range(number_of_cells), dtype=np.float64
    )


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
def j_blues(lines):
    return pd.DataFrame(1.0e-5, index=lines.index, columns=range(20))


@pytest.fixture
def link_t_rad_t_electron():
    return 0.9


# GENERAL PROPERTIES


@pytest.fixture
def selected_atoms(abundance):
    selected_atoms_module = SelectedAtoms(None)
    return selected_atoms_module.calculate(abundance)


@pytest.fixture
def beta_rad(t_rad):
    beta_rad_module = BetaRadiation(None)
    return beta_rad_module.calculate(t_rad)


@pytest.fixture
def g_electron(beta_rad):
    g_electron_module = GElectron(None)
    return g_electron_module.calculate(beta_rad)


@pytest.fixture
def number_density(atomic_mass, abundance, density):
    number_density_module = NumberDensity(None)
    return number_density_module.calculate(atomic_mass, abundance, density)


@pytest.fixture
def t_electrons(t_rad, link_t_rad_t_electron):
    electron_temperature_module = ElectronTemperature(None)
    return electron_temperature_module.calculate(t_rad, link_t_rad_t_electron)


@pytest.fixture
def beta_electron(t_electrons):
    beta_electron_module = BetaElectron(None)
    return beta_electron_module.calculate(t_electrons)


# ATOMIC PROPERTIES


@pytest.fixture
def levels(atomic_data, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_data, selected_atoms)[0]


@pytest.fixture
def excitation_energy(atomic_data, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_data, selected_atoms)[1]


@pytest.fixture
def metastability(atomic_data, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_data, selected_atoms)[2]


@pytest.fixture
def g(atomic_data, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_data, selected_atoms)[3]


@pytest.fixture
def lines(atomic_data, selected_atoms):
    lines_module = Lines(None)
    return lines_module.calculate(atomic_data, selected_atoms)[0]


@pytest.fixture
def nu(atomic_data, selected_atoms):
    lines_module = Lines(None)
    return lines_module.calculate(atomic_data, selected_atoms)[1]


@pytest.fixture
def f_lu(atomic_data, selected_atoms):
    lines_module = Lines(None)
    return lines_module.calculate(atomic_data, selected_atoms)[2]


@pytest.fixture
def wavelength_cm(atomic_data, selected_atoms):
    lines_module = Lines(None)
    return lines_module.calculate(atomic_data, selected_atoms)[3]


@pytest.fixture
def ionization_data(atomic_data, selected_atoms):
    ionization_data_module = IonizationData(None)
    return ionization_data_module.calculate(atomic_data, selected_atoms)


@pytest.fixture
def atomic_mass(atomic_data, selected_atoms):
    atomic_mass_module = AtomicMass(None)
    return atomic_mass_module.calculate(atomic_data, selected_atoms)


@pytest.fixture
def zeta_data(atomic_data, selected_atoms):
    zeta_data_module = ZetaData(None)
    return zeta_data_module.calculate(atomic_data, selected_atoms)


@pytest.fixture
def lines_upper_level_index(lines, levels):
    upper_level_index_module = LinesUpperLevelIndex(None)
    return upper_level_index_module.calculate(levels, lines)


@pytest.fixture
def lines_lower_level_index(lines, levels):
    lower_level_index_module = LinesLowerLevelIndex(None)
    return lower_level_index_module.calculate(levels, lines)


@pytest.fixture
def chi_0(atomic_data):
    chi_0_module = Chi0(None)
    return chi_0_module.calculate(atomic_data)


# PARTITION FUNCTION PROPERTIES


@pytest.fixture
def level_boltzmann_factor_lte(excitation_energy, g, beta_rad, levels):
    level_boltzmann_factor_module = LevelBoltzmannFactorLTE(None)
    return level_boltzmann_factor_module.calculate(
        excitation_energy, g, beta_rad, levels
    )


@pytest.fixture
def level_boltzmann_factor_dilute_lte(
    levels, g, excitation_energy, beta_rad, w, metastability
):
    level_boltzmann_factor_module = LevelBoltzmannFactorDiluteLTE(None)
    return level_boltzmann_factor_module.calculate(
        levels, g, excitation_energy, beta_rad, w, metastability
    )


@pytest.fixture
def partition_function(level_boltzmann_factor_lte):
    partition_function_module = PartitionFunction(None)
    return partition_function_module.calculate(level_boltzmann_factor_lte)


# ION POPULATION PROPERTIES


@pytest.fixture
def phi_saha_lte(g_electron, beta_rad, partition_function, ionization_data):
    phi_saha_lte_module = PhiSahaLTE(None)
    return phi_saha_lte_module.calculate(
        g_electron, beta_rad, partition_function, ionization_data
    )


@pytest.fixture
def phi_saha_nebular(
    t_rad,
    w,
    zeta_data,
    t_electrons,
    delta,
    g_electron,
    beta_rad,
    partition_function,
    ionization_data,
):
    phi_saha_nebular_module = PhiSahaNebular(None)
    return phi_saha_nebular_module.calculate(
        t_rad,
        w,
        zeta_data,
        t_electrons,
        delta,
        g_electron,
        beta_rad,
        partition_function,
        ionization_data,
    )


@pytest.fixture
def ion_number_density(phi_saha_lte, partition_function, number_density):
    ion_number_density_module = IonNumberDensity(None)
    ion_number_density, electron_densities = (
        ion_number_density_module.calculate(
            phi_saha_lte, partition_function, number_density
        )
    )
    return ion_number_density


@pytest.fixture
def electron_densities(phi_saha_lte, partition_function, number_density):
    electron_density_module = IonNumberDensity(None)
    ion_number_density, electron_densities = electron_density_module.calculate(
        phi_saha_lte, partition_function, number_density
    )
    return electron_densities


@pytest.fixture
def delta(
    w, ionization_data, beta_rad, t_electrons, t_rad, beta_electron, levels
):
    delta_module = RadiationFieldCorrection(chi_0_species=(2, 1))
    return delta_module.calculate(
        w, ionization_data, beta_rad, t_electrons, t_rad, beta_electron
    )


# LEVEL POPULATION PROPERTIES


@pytest.fixture
def level_number_density(
    level_boltzmann_factor_lte, ion_number_density, levels, partition_function
):
    level_number_density_module = LevelNumberDensity(None)
    return level_number_density_module.calculate(
        level_boltzmann_factor_lte,
        ion_number_density,
        levels,
        partition_function,
    )


# RADIATIVE PROPERTIES


@pytest.fixture
def stimulated_emission_factor(
    g,
    level_number_density,
    lines_lower_level_index,
    lines_upper_level_index,
    metastability,
    lines,
):
    factor_module = StimulatedEmissionFactor(nlte_species=None)
    return factor_module.calculate(
        g,
        level_number_density,
        lines_lower_level_index,
        lines_upper_level_index,
        metastability,
        lines,
    )


@pytest.fixture
def tau_sobolev(
    lines,
    level_number_density,
    lines_lower_level_index,
    time_explosion,
    stimulated_emission_factor,
    j_blues,
    f_lu,
    wavelength_cm,
):
    tau_sobolev_module = TauSobolev(None)
    return tau_sobolev_module.calculate(
        lines,
        level_number_density,
        lines_lower_level_index,
        time_explosion,
        stimulated_emission_factor,
        j_blues,
        f_lu,
        wavelength_cm,
    )


@pytest.fixture
def beta_sobolev(tau_sobolev):
    beta_sobolev_module = BetaSobolev(None)
    return beta_sobolev_module.calculate(tau_sobolev)


@pytest.fixture
def transition_probabilities(
    atomic_data, beta_sobolev, j_blues, stimulated_emission_factor, tau_sobolev
):
    transition_probabilities_module = TransitionProbabilities(None)
    return transition_probabilities_module.calculate(
        atomic_data,
        beta_sobolev,
        j_blues,
        stimulated_emission_factor,
        tau_sobolev,
    )
