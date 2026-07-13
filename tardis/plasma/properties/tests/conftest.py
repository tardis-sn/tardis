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
    return np.linspace(0.7e-14, 1.3e-14, number_of_cells)


@pytest.fixture
def t_rad(number_of_cells):
    return np.linspace(9500.0, 10500.0, number_of_cells)

@pytest.fixture
def w(number_of_cells):
    return np.linspace(0.85, 1.15, number_of_cells)

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
def beta_electron(t_electrons):
    beta_electron_module = BetaElectron(None)
    return beta_electron_module.calculate(t_electrons)

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
def metastability(atomic_dataset, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_dataset, selected_atoms)[2]

@pytest.fixture
def g(atomic_dataset, selected_atoms):
    levels_module = Levels(None)
    return levels_module.calculate(atomic_dataset, selected_atoms)[3]

@pytest.fixture
def lines(atomic_dataset, selected_atoms):
    lines_module = Lines(None)
    lines, _, _, _ = lines_module.calculate(atomic_dataset, selected_atoms)
    return lines[lines.index.isin(selected_atoms, level="atomic_number")]

@pytest.fixture
def lines_lower_level_index(levels, lines):
    lines_lower_level_index_module = LinesLowerLevelIndex(None)
    return lines_lower_level_index_module.calculate(levels, lines)

@pytest.fixture
def lines_upper_level_index(levels, lines):
    lines_upper_level_index_module = LinesUpperLevelIndex(None)
    return lines_upper_level_index_module.calculate(levels, lines)

# PARTITION FUNCTION PROPERTIES

@pytest.fixture
def level_boltzmann_factor_lte(excitation_energy, g, beta_rad, levels):
    level_boltzmann_factor_module = LevelBoltzmannFactorLTE(None)
    return level_boltzmann_factor_module.calculate(
        excitation_energy, g, beta_rad, levels
    )

@pytest.fixture
def partition_function(level_boltzmann_factor_lte):
    partition_function_module = PartitionFunction(None)
    return partition_function_module.calculate(level_boltzmann_factor_lte)

# ION / LEVEL POPULATION PROPERTIES

@pytest.fixture
def ionization_data(atomic_dataset, selected_atoms):
    ionization_data_module = IonizationData(None)
    return ionization_data_module.calculate(atomic_dataset, selected_atoms)

@pytest.fixture
def phi(g_electron, beta_rad, partition_function, ionization_data):
    phi_module = PhiSahaLTE(None)
    return phi_module.calculate(
        g_electron, beta_rad, partition_function, ionization_data
    )

@pytest.fixture
def ion_number_density(phi, partition_function, number_density):
    ion_number_density_module = IonNumberDensity(None)
    ion_number_density, _ = ion_number_density_module.calculate(
        phi, partition_function, number_density
    )
    return ion_number_density

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
    stimulated_emission_factor_module = StimulatedEmissionFactor(nlte_species=None)
    return stimulated_emission_factor_module.calculate(
        g,
        level_number_density,
        lines_lower_level_index,
        lines_upper_level_index,
        metastability,
        lines,
    )

# STIMULATED EMISSION FACTORS TWO LEVEL INPUT
@pytest.fixture
def two_level_inputs():
    """Return a builder for minimal 2-level, 1-zone, 1-line edge-case inputs."""

    def _build(
        n_lower: float,
        n_upper: float,
        g_lower: float,
        g_upper: float,
        metastable_upper: bool = False,
    ) -> dict:
        levels_index = pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        )
        lines_index = pd.MultiIndex.from_tuples(
            [(1, 0, 0, 1)],
            names=[
                "atomic_number",
                "ion_number",
                "level_number_lower",
                "level_number_upper",
            ],
        )
        return dict(
            g=pd.Series([g_lower, g_upper], index=levels_index),
            level_number_density=pd.DataFrame(
                [[n_lower], [n_upper]], index=levels_index, columns=[0]
            ),
            metastability=pd.Series([False, metastable_upper], index=levels_index),
            lines=pd.DataFrame({"wavelength": [6562.8]}, index=lines_index),
            lines_lower_level_index=np.array([0], dtype=np.int64),
            lines_upper_level_index=np.array([1], dtype=np.int64),
        )

    return _build