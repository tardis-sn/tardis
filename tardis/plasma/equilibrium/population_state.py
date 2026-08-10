from dataclasses import dataclass

import pandas as pd

from tardis.plasma.equilibrium.rate_matrix import ElementalStateIndex


@dataclass(frozen=True)
class SingleElementPopulationState:
    """Element-normalized and absolute populations for one element.

    Attributes
    ----------
    normalized_state_populations : pandas.DataFrame
        Population for every explicit level and total-ion state in the
        matrix ordering described by ``state_index.states``. Columns are
        shells.
    normalized_level_populations : pandas.DataFrame
        Element-normalized populations for explicit levels, indexed by
        ``(atomic_number, ion_number, level_number)``.
    normalized_ion_populations : pandas.DataFrame
        Element-normalized total population for every retained ion stage,
        including the terminal bare-nucleus state.
    level_populations : pandas.DataFrame
        Absolute explicit-level populations reconstructed from the elemental
        number density.
    ion_populations : pandas.DataFrame
        Absolute total-ion populations reconstructed from the elemental
        number density.
    state_index : ElementalStateIndex
        Mapping between physical level/ion labels and matrix positions.
    """

    normalized_state_populations: pd.DataFrame
    normalized_level_populations: pd.DataFrame
    normalized_ion_populations: pd.DataFrame
    level_populations: pd.DataFrame
    ion_populations: pd.DataFrame
    state_index: ElementalStateIndex


@dataclass(frozen=True)
class PopulationState:
    """Complete population state at one electron-density solution.

    Attributes
    ----------
    electron_densities : pandas.Series
        Electron number density in each shell.
    elemental_populations : dict[int, SingleElementPopulationState]
        Normalized and absolute solutions for every element.
    ion_number_density : pandas.DataFrame
        Absolute populations for every ion stage of every element.
    level_number_density : pandas.DataFrame
        Absolute populations for every atomic level.
    level_boltzmann_factor : pandas.DataFrame
        Level factors consistent with ``level_number_density``.
    """

    electron_densities: pd.Series
    elemental_populations: dict[int, SingleElementPopulationState]
    ion_number_density: pd.DataFrame
    level_number_density: pd.DataFrame
    level_boltzmann_factor: pd.DataFrame
