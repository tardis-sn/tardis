import os

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

import tardis
from tardis.plasma.properties import *
from tardis.io.atom_data import AtomData
from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)

# INPUTS

@pytest.fixture
def number_of_cells():
    return 20

@pytest.fixture
def abundance(number_of_cells):
    return pd.DataFrame(
        data=[[0.8] * number_of_cells, [0.2] * number_of_cells],
        index=pd.Index([1, 2], name="atomic_number"),
        columns=range(number_of_cells),
        dtype=np.float64,
    )


@pytest.fixture
def density(number_of_cells):
    return np.ones(number_of_cells) * 1e-14


@pytest.fixture
def t_rad(number_of_cells):
    return np.ones(number_of_cells) * 10000

@pytest.fixture
def w(number_of_cells):
    return np.ones(number_of_cells) * 1.0

@pytest.fixture
def dilute_planckian_radiation_field(t_rad, w):
    return DilutePlanckianRadiationField(t_rad * u.K, w)

# GENERAL PROPERTIES

@pytest.fixture
def t_rad_calculated(dilute_planckian_radiation_field):
    t_rad_module = TRadiative(None)
    return t_rad_module.calculate(dilute_planckian_radiation_field)

@pytest.fixture
def beta_rad(t_rad):
    beta_rad_module = BetaRadiation(None)
    return beta_rad_module.calculate(t_rad)

@pytest.fixture
def g_electron(beta_rad):
    g_electron_module = GElectron(None)
    return g_electron_module.calculate(beta_rad)

@pytest.fixture
def link_t_rad_t_electron():
    return 0.9

@pytest.fixture
def t_electrons(t_rad, link_t_rad_t_electron):
    electron_temperature_module = ElectronTemperature(None)
    return electron_temperature_module.calculate(t_rad, link_t_rad_t_electron)

@pytest.fixture
def number_density(atomic_dataset, abundance, density):
    atomic_mass = atomic_dataset.atom_data.reindex(abundance.index)["mass"]
    return abundance.mul(density, axis=1).div(atomic_mass, axis=0)

@pytest.fixture
def selected_atoms(number_density):
    selected_atoms_module = SelectedAtoms(None)
    return selected_atoms_module.calculate(number_density)

# ATOMIC PROPERTIES

@pytest.fixture
def levels(atomic_dataset, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_dataset, selected_atoms)[0]


@pytest.fixture
def excitation_energy(atomic_dataset, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_dataset, selected_atoms)[1]

@pytest.fixture
def g(atomic_dataset, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_dataset, selected_atoms)[3]

# PARTITION FUNCTION PROPERTIES

@pytest.fixture
def level_boltzmann_factor_lte(excitation_energy, g, beta_rad, levels):
    level_boltzmann_factor_module = LevelBoltzmannFactorLTE(None)
    return level_boltzmann_factor_module.calculate(
        excitation_energy, g, beta_rad, levels
    )