import numpy as np
from astropy import units as u

from tardis.transport.montecarlo.numba_interface import (
    LineInteractionType,
)

ENABLE_FULL_RELATIVITY = False
TEMPORARY_V_PACKET_BINS = 0
NUMBER_OF_VPACKETS = 0
MONTECARLO_SEED = 0
LINE_INTERACTION_TYPE = 0
PACKET_SEEDS = np.empty(1, dtype=np.int64)
DISABLE_ELECTRON_SCATTERING = False
DISABLE_LINE_SCATTERING = False
SURVIVAL_PROBABILITY = 0.0
VPACKET_TAU_RUSSIAN = 10.0

INITIAL_TRACKING_ARRAY_LENGTH = 0
LEGACY_MODE_ENABLED = False

ENABLE_RPACKET_TRACKING = False
CONTINUUM_PROCESSES_ENABLED = False

VPACKET_SPAWN_START_FREQUENCY = 0
VPACKET_SPAWN_END_FREQUENCY = 1e200
ENABLE_VPACKET_TRACKING = False


def configuration_initialize(config, transport, number_of_vpackets):
    if transport.line_interaction_type == "macroatom":
        config.LINE_INTERACTION_TYPE = LineInteractionType.MACROATOM
    elif transport.line_interaction_type == "downbranch":
        config.LINE_INTERACTION_TYPE = LineInteractionType.DOWNBRANCH
    elif transport.line_interaction_type == "scatter":
        config.LINE_INTERACTION_TYPE = LineInteractionType.SCATTER
    else:
        raise ValueError(
            f'Line interaction type must be one of "macroatom",'
            f'"downbranch", or "scatter" but is '
            f"{transport.line_interaction_type}"
        )
    NUMBER_OF_VPACKETS = number_of_vpackets
    TEMPORARY_V_PACKET_BINS = number_of_vpackets
    ENABLE_FULL_RELATIVITY = transport.enable_full_relativity
    MONTECARLO_SEED = transport.packet_source.base_seed
    VPACKET_SPAWN_START_FREQUENCY = (
        transport.virtual_spectrum_spawn_range.end.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    VPACKET_SPAWN_END_FREQUENCY = (
        transport.virtual_spectrum_spawn_range.start.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    ENABLE_VPACKET_TRACKING = transport.enable_vpacket_tracking
    ENABLE_RPACKET_TRACKING = transport.enable_rpacket_tracking
