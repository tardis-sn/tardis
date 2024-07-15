from numba import float64, int64, boolean
from numba.experimental import jitclass
import numpy as np
from astropy import units as u

from tardis.transport.montecarlo.configuration.constants import (
    LineInteractionType,
)

from tardis.transport.montecarlo.configuration import montecarlo_globals

montecarlo_config_spec = [
    ("TEMPORARY_V_PACKET_BINS", int64),
    ("NUMBER_OF_VPACKETS", int64),
    ("MONTECARLO_SEED", int64),
    ("PACKET_SEEDS", int64[:]),
    ("SURVIVAL_PROBABILITY", float64),
    ("VPACKET_TAU_RUSSIAN", float64),
    ("INITIAL_TRACKING_ARRAY_LENGTH", int64),
    ("VPACKET_SPAWN_START_FREQUENCY", float64),
    ("VPACKET_SPAWN_END_FREQUENCY", float64),
]


@jitclass(montecarlo_config_spec)
class MonteCarloConfiguration:
    def __init__(self):
        self.TEMPORARY_V_PACKET_BINS = 0
        self.NUMBER_OF_VPACKETS = 0
        self.MONTECARLO_SEED = 0

        self.PACKET_SEEDS = np.empty(1, dtype=np.int64)
        self.SURVIVAL_PROBABILITY = 0.0
        self.VPACKET_TAU_RUSSIAN = 10.0

        self.INITIAL_TRACKING_ARRAY_LENGTH = 1

        self.VPACKET_SPAWN_START_FREQUENCY = 0
        self.VPACKET_SPAWN_END_FREQUENCY = 1e200


def initialize_configuration(
    montecarlo_configuration, config, enable_virtual_packet_logging
):
    # set globals
    if config.plasma.line_interaction_type == "macroatom":
        montecarlo_globals.LINE_INTERACTION_TYPE = LineInteractionType.MACROATOM
    elif config.plasma.line_interaction_type == "downbranch":
        montecarlo_globals.LINE_INTERACTION_TYPE = (
            LineInteractionType.DOWNBRANCH
        )
    elif config.plasma.line_interaction_type == "scatter":
        montecarlo_globals.LINE_INTERACTION_TYPE = LineInteractionType.SCATTER
    else:
        raise ValueError(
            f'Line interaction type must be one of "macroatom",'
            f'"downbranch", or "scatter" but is '
            f"{config.plasma.line_interaction_type}"
        )
    montecarlo_globals.ENABLE_FULL_RELATIVITY = (
        config.montecarlo.enable_full_relativity
    )
    montecarlo_globals.ENABLE_VPACKET_TRACKING = (
        config.spectrum.virtual.virtual_packet_logging
        | enable_virtual_packet_logging
    )
    montecarlo_globals.ENABLE_RPACKET_TRACKING = (
        config.montecarlo.tracking.track_rpacket
    )

    montecarlo_globals.DISABLE_ELECTRON_SCATTERING = (
        config.plasma.disable_electron_scattering
    )
    montecarlo_globals.DISABLE_LINE_SCATTERING = (
        config.plasma.disable_line_scattering
    )

    # non-globals
    montecarlo_configuration.MONTECARLO_SEED = config.montecarlo.seed
    montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY = (
        config.montecarlo.virtual_spectrum_spawn_range.end.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY = (
        config.montecarlo.virtual_spectrum_spawn_range.start.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
