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

        self.INITIAL_TRACKING_ARRAY_LENGTH = 0

        self.VPACKET_SPAWN_START_FREQUENCY = 0
        self.VPACKET_SPAWN_END_FREQUENCY = 1e200


def configuration_initialize(
    montecarlo_configuration, transport_solver, number_of_vpackets
):
    if transport_solver.line_interaction_type == "macroatom":
        montecarlo_globals.LINE_INTERACTION_TYPE = LineInteractionType.MACROATOM
    elif transport_solver.line_interaction_type == "downbranch":
        montecarlo_globals.LINE_INTERACTION_TYPE = (
            LineInteractionType.DOWNBRANCH
        )
    elif transport_solver.line_interaction_type == "scatter":
        montecarlo_globals.LINE_INTERACTION_TYPE = LineInteractionType.SCATTER
    else:
        raise ValueError(
            f'Line interaction type must be one of "macroatom",'
            f'"downbranch", or "scatter" but is '
            f"{transport_solver.line_interaction_type}"
        )
    montecarlo_configuration.NUMBER_OF_VPACKETS = number_of_vpackets
    montecarlo_configuration.TEMPORARY_V_PACKET_BINS = number_of_vpackets
    montecarlo_globals.ENABLE_FULL_RELATIVITY = (
        transport_solver.enable_full_relativity
    )
    montecarlo_configuration.MONTECARLO_SEED = (
        transport_solver.packet_source.base_seed
    )
    montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY = (
        transport_solver.virtual_spectrum_spawn_range.end.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY = (
        transport_solver.virtual_spectrum_spawn_range.start.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_globals.ENABLE_VPACKET_TRACKING = (
        transport_solver.enable_vpacket_tracking
    )
    montecarlo_globals.ENABLE_RPACKET_TRACKING = (
        transport_solver.enable_rpacket_tracking
    )

    montecarlo_globals.LEGACY_MODE_ENABLED = transport_solver.enable_legacy_mode

    montecarlo_globals.DISABLE_ELECTRON_SCATTERING = (
        transport_solver.disable_electron_scattering
    )
    montecarlo_globals.DISABLE_LINE_SCATTERING = (
        transport_solver.disable_line_scattering
    )
