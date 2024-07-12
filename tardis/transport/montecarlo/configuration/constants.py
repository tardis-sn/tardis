from enum import IntEnum
from tardis import constants as const

SIGMA_THOMSON = const.sigma_T.to("cm^2").value
CLOSE_LINE_THRESHOLD = 1e-14
C_SPEED_OF_LIGHT = const.c.to("cm/s").value
MISS_DISTANCE = 1e99
KB = const.k_B.cgs.value
H = const.h.cgs.value


class LineInteractionType(IntEnum):
    SCATTER = 0
    DOWNBRANCH = 1
    MACROATOM = 2
