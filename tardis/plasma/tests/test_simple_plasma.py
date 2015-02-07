import pytest

import pandas as pd
import numpy as np

from tardis.plasma.base import BasePlasma
from tardis.plasma.plasma_properties import (BetaRadiation,
                                             LevelBoltzmannFactor,
                                             AtomicLevels, AtomicLines,
                                             SelectedAtoms, AtomicMass,
                                             LTEPartitionFunction)

from tardis.plasma.plasma_input import TRadiative, AtomicData, Abundance


@pytest.fixture
def abundance():
    selected_atoms = [8, 12, 14, 16, 18, 20]
    return pd.DataFrame(data=np.random.random((6, 20)), index=selected_atoms,
                        columns=range(20), dtype=np.float64)

@pytest.fixture
def t_rad():
    return np.ones(20) * 10000.


def test_simple_networkx_test1(atomic_data, abundance, t_rad):
    modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
               AtomicLevels, AtomicLines, AtomicData, Abundance, SelectedAtoms,
               AtomicMass, LTEPartitionFunction]
    bp = BasePlasma(modules, t_rad=t_rad, atomic_data=atomic_data,
                    abundance=abundance)
    assert bp.t_rad[0] == 5000
