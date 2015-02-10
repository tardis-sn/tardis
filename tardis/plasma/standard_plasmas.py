import logging

import numpy as np
import pandas as pd

from tardis.plasma import BasePlasma
from tardis.plasma.plasma_input import (TRadiative, AtomicData, Abundance,
                                        Density, TimeExplosion)
from tardis.plasma.base_properties import (
    BetaRadiation, LevelBoltzmannFactor, Levels, Lines,
    SelectedAtoms, AtomicMass, LTEPartitionFunction, LevelPopulationLTE,
    LevelNumberDensity, PhiSahaLTE, GElectron, IonizationData, NumberDensity,
    IonNumberDensity, LinesLowerLevelIndex, LinesUpperLevelIndex, TauSobolev)

logger = logging.getLogger(__name__)

class LTEPlasma(BasePlasma):

    def __init__(self, t_rad, abundance, density, time_explosion, atomic_data):
        plasma_modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
               Levels, Lines, AtomicData, Abundance, SelectedAtoms,
               AtomicMass, LTEPartitionFunction, LevelPopulationLTE, PhiSahaLTE,
               GElectron, IonizationData, Density, NumberDensity,
               IonNumberDensity, LevelNumberDensity, LinesLowerLevelIndex,
               LinesUpperLevelIndex, TauSobolev, TimeExplosion]
        super(LTEPlasma, self).__init__(plasma_modules=plasma_modules,
                                        t_rad=t_rad, abundance=abundance,
                                        atomic_data=atomic_data,
                                        density=density,
                                        time_explosion=time_explosion)

class LegacyPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte'):
        pass
