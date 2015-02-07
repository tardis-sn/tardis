from tardis.plasma.base import BasePlasma
from tardis.plasma.plasma_properties import (BetaRadiation,
                                             LevelBoltzmannFactor,
                                             AtomicLevels, AtomicLines)

from tardis.plasma.plasma_input import TRadiative, AtomicData


def test_simple_networkx_test1(atomic_data):
    modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
               AtomicLevels, AtomicLines, AtomicData]
    bp = BasePlasma(modules, t_rad=5000, atomic_data=atomic_data)
    assert bp.t_rad == 5000
    1/0