import pytest

import pandas as pd
import numpy as np

from astropy import units as u

from tardis.plasma.plasma_input import TRadiative, AtomicData, Abundance, DilutionFactor

from tardis.plasma.standard_plasmas import LTEPlasma



#def test_simple_networkx_test1(atomic_data, abundance, t_rad):
#    modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
#               AtomicLevels, AtomicLines, AtomicData, Abundance, SelectedAtoms,
#               AtomicMass]
#    bp = BasePlasma(modules, t_rad=t_rad, atomic_data=atomic_data,
#                    abundance=abundance)
#    assert bp.t_rad[0] == 5000


def test_simple_lte_plasma(included_he_atomic_data, abundance, t_rad, density,
                           time_explosion):
    lte_plasma = LTEPlasma(t_rad, abundance, density, time_explosion,
                           included_he_atomic_data)

