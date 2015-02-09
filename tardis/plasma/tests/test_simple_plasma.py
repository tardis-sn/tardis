import pytest

import pandas as pd
import numpy as np

from tardis.plasma.base import BasePlasma
from tardis.plasma.base_properties import (BetaRadiation,
                                             LevelBoltzmannFactor,
                                             AtomicLevels, AtomicLines,
                                             SelectedAtoms, AtomicMass,
                                             DiluteLTEPartitionFunction)

from tardis.plasma.plasma_input import TRadiative, AtomicData, Abundance, DilutionFactor

from tardis.plasma.standard_plasmas import LTEPlasma


@pytest.fixture
def number_of_cells():
    return 20

@pytest.fixture
def abundance(number_of_cells):
    selected_atoms = [8, 12, 14, 16, 18, 20]
    return pd.DataFrame(data=np.random.random((6, number_of_cells)), index=selected_atoms,
                        columns=range(20), dtype=np.float64)

@pytest.fixture
def t_rad(number_of_cells):
    return np.ones(number_of_cells) * 10000.


@pytest.fixture
def w(number_of_cells):
    return np.ones(number_of_cells) * 0.5



def test_simple_networkx_test1(atomic_data, abundance, t_rad):
    modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
               AtomicLevels, AtomicLines, AtomicData, Abundance, SelectedAtoms,
               AtomicMass]
    bp = BasePlasma(modules, t_rad=t_rad, atomic_data=atomic_data,
                    abundance=abundance)
    assert bp.t_rad[0] == 5000


def test_simple_lte_plasma(atomic_data, abundance, t_rad):
    lte_plasma = LTEPlasma(t_rad, abundance, atomic_data)
    lte_plasma.write_to_dot('test.dot')
    1/0