from tardis import constants as const

SIGMA_THOMSON = const.sigma_T.to("cm^2").value


def sigma_thomson_for_transport(disable_electron_scattering: bool) -> float:
    """Effective Thomson cross section for MC and formal-integral transport."""
    if disable_electron_scattering:
        return 1e-200
    return SIGMA_THOMSON
CLOSE_LINE_THRESHOLD = 1e-14
C_SPEED_OF_LIGHT = const.c.to("cm/s").value
MISS_DISTANCE = 1e99
KB = const.k_B.cgs.value
H = const.h.cgs.value
