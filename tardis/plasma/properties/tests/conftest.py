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
def g(atomic_data, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_data, selected_atoms)[3]

# PARTITION FUNCTION PROPERTIES

@pytest.fixture
def level_boltzmann_factor_lte(excitation_energy, g, beta_rad, levels):
    level_boltzmann_factor_module = LevelBoltzmannFactorLTE(None)
    return level_boltzmann_factor_module.calculate(
        excitation_energy, g, beta_rad, levels
    )