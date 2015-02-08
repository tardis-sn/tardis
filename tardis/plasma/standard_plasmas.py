import logging

import numpy as np
import pandas as pd

from tardis.plasma import BasePlasma
from tardis.plasma.plasma_input import TRadiative, AtomicData, Abundance
from tardis.plasma.plasma_properties import (
    BetaRadiation, LevelBoltzmannFactor, AtomicLevels, AtomicLines,
    SelectedAtoms, AtomicMass, LTEPartitionFunction)

logger = logging.getLogger(__name__)

class LTEPlasma(BasePlasma):

    def __init__(self, t_rad, abundance, atomic_data):
        plasma_modules = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
               AtomicLevels, AtomicLines, AtomicData, Abundance, SelectedAtoms,
               AtomicMass, LTEPartitionFunction]
        super(LTEPlasma, self).__init__(plasma_modules=plasma_modules,
                                        t_rad=t_rad, abundance=abundance,
                                        atomic_data=atomic_data)

