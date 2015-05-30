import numpy as np

from tardis.plasma.properties import (TRadiative, BetaRadiation,
LevelBoltzmannFactor, Levels, Lines, AtomicData, Abundance,
SelectedAtoms, AtomicMass)
from tardis.plasma import BasePlasma

def test_simple_networkx_test(included_he_atomic_data, abundance, t_rad):
    modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
               Levels, Lines, AtomicData, Abundance, SelectedAtoms,
               AtomicMass]
    bp = BasePlasma(modules, t_rad=t_rad, atomic_data=included_he_atomic_data,
                    abundance=abundance)
    assert bp.t_rad[0] == 10000

def test_simple_lte_plasma(standard_lte_plasma_he_db):
    assert np.allclose(standard_lte_plasma_he_db.tau_sobolev.ix[564954],
        0.123817)