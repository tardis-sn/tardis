import pytest

import pandas as pd
import numpy as np

from tardis.plasma.base import BasePlasma
from tardis.plasma.plasma_properties import (BetaRadiation,
                                             LevelBoltzmannFactor,
                                             AtomicLevels, AtomicLines, SelectedAtoms)

from tardis.plasma.plasma_input import TRadiative, AtomicData, Abundance


@pytest.fixture
def abundance():
    selected_atoms = [8, 12, 14, 16, 18, 20]
    return pd.DataFrame(data=np.random.random((6, 20)), index=selected_atoms,
                        columns=range(20), dtype=np.float64)


def test_simple_networkx_test1(atomic_data, abundance):
    modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
               AtomicLevels, AtomicLines, AtomicData, Abundance, SelectedAtoms]
    bp = BasePlasma(modules, t_rad=5000, atomic_data=atomic_data,
                    abundance=abundance)
    assert bp.t_rad == 5000
    1/0