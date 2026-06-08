import os

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

import tardis
from tardis.plasma.properties import *
from tardis.io.atom_data import AtomData

# INPUTS

@pytest.fixture
def number_of_cells():
    return 20

@pytest.fixture
def t_rad(number_of_cells):
    return np.ones(number_of_cells) * 10000

# GENERAL PROPERTIES

@pytest.fixture
def beta_rad(t_rad):
    beta_rad_module = BetaRadiation(None)
    return beta_rad_module.calculate(t_rad)

@pytest.fixture
def g_electron(beta_rad):
    g_electron_module = GElectron(None)
    return g_electron_module.calculate(beta_rad)

@pytest.fixture
def number_density(number_of_cells):
    return pd.DataFrame(
        data=1.0, index=[1, 2], columns=range(number_of_cells), dtype=np.float64
    )

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