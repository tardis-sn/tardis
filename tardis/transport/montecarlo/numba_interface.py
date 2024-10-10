from enum import IntEnum

import numpy as np
from numba import float64, int64
from numba.experimental import jitclass

from tardis import constants as const
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.opacities.opacity_state import opacity_state_initialize

C_SPEED_OF_LIGHT = const.c.to("cm/s").value



class LineInteractionType(IntEnum):
    SCATTER = 0
    DOWNBRANCH = 1
    MACROATOM = 2
