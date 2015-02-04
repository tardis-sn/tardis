from tardis.plasma.base import BasePlasma
from tardis.plasma.plasma_properties import BetaRadiation, LevelBoltzmannFactor
def test_simple_networkx_test1():
    bp = BasePlasma([BetaRadiation, LevelBoltzmannFactor])
