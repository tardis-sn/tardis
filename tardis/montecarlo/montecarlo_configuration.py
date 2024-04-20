from astropy import units as u
from numba import float64, int64, boolean
from numba.experimental import jitclass
import numpy as np

from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType,
)

numba_config_spec = [
    ("ENABLE_FULL_RELATIVITY", boolean),
    ("TEMPORARY_V_PACKET_BINS", int64),
    ("NUMBER_OF_VPACKETS", int64),
    ("MONTECARLO_SEED", int64),
    ("LINE_INTERACTION_TYPE", int64),
    ("PACKET_SEEDS", int64[:]),
    ("DISABLE_ELECTRON_SCATTERING", boolean),
    ("DISABLE_LINE_SCATTERING", boolean),
    ("SURVIVAL_PROBABILITY", float64),
    ("VPACKET_TAU_RUSSIAN", float64),
    ("INITIAL_TRACKING_ARRAY_LENGTH", int64),
    ("LEGACY_MODE_ENABLED", boolean),
    ("ENABLE_RPACKET_TRACKING", boolean),
    ("CONTINUUM_PROCESSES_ENABLED", boolean),
    ("VPACKET_SPAWN_START_FREQUENCY", float64),
    ("VPACKET_SPAWN_END_FREQUENCY", float64),
    ("ENABLE_VPACKET_TRACKING", boolean),
]


@jitclass(numba_config_spec)
class MonteCarloConfiguration(object):
    def __init__(self):
        self.ENABLE_FULL_RELATIVITY = False
        self.TEMPORARY_V_PACKET_BINS = 0
        self.NUMBER_OF_VPACKETS = 0
        self.MONTECARLO_SEED = 0
        self.LINE_INTERACTION_TYPE = 0
        self.PACKET_SEEDS = np.empty(1, dtype=np.int64)
        self.DISABLE_ELECTRON_SCATTERING = False
        self.DISABLE_LINE_SCATTERING = False
        self.SURVIVAL_PROBABILITY = 0.0
        self.VPACKET_TAU_RUSSIAN = 10.0

        self.INITIAL_TRACKING_ARRAY_LENGTH = 0
        self.LEGACY_MODE_ENABLED = False

        self.ENABLE_RPACKET_TRACKING = False
        self.CONTINUUM_PROCESSES_ENABLED = False

        self.VPACKET_SPAWN_START_FREQUENCY = 0
        self.VPACKET_SPAWN_END_FREQUENCY = 1e200
        self.ENABLE_VPACKET_TRACKING = False


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
    config.NUMBER_OF_VPACKETS = number_of_vpackets
    config.TEMPORARY_V_PACKET_BINS = number_of_vpackets
    config.ENABLE_FULL_RELATIVITY = transport.enable_full_relativity
    config.MONTECARLO_SEED = transport.packet_source.base_seed
    config.VPACKET_SPAWN_START_FREQUENCY = (
        transport.virtual_spectrum_spawn_range.end.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    config.VPACKET_SPAWN_END_FREQUENCY = (
        transport.virtual_spectrum_spawn_range.start.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    config.ENABLE_VPACKET_TRACKING = transport.enable_vpacket_tracking
    config.ENABLE_RPACKET_TRACKING = transport.enable_rpacket_tracking
