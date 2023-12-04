from astropy import units as u
from tardis import constants as const
from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType,
)

ENABLE_FULL_RELATIVITY = False
TEMPORARY_V_PACKET_BINS = 0
NUMBER_OF_VPACKETS = 0
MONTECARLO_SEED = 0
LINE_INTERACTION_TYPE = None
PACKET_SEEDS = []
DISABLE_ELECTRON_SCATTERING = False
DISABLE_LINE_SCATTERING = False
SURVIVAL_PROBABILITY = 0.0
TAU_RUSSIAN = 10.0


INITIAL_TRACKING_ARRAY_LENGTH = None
LEGACY_MODE_ENABLED = False

ENABLE_RPACKET_TRACKING = False
CONTINUUM_PROCESSES_ENABLED = False


VPACKET_SPAWN_START_FREQUENCY = 0
VPACKET_SPAWN_END_FREQUENCY = 1e200
ENABLE_VPACKET_TRACKING = False


def configuration_initialize(transport, number_of_vpackets):
    if transport.line_interaction_type == "macroatom":
        montecarlo_configuration.LINE_INTERACTION_TYPE = (
            LineInteractionType.MACROATOM
        )
    elif transport.line_interaction_type == "downbranch":
        montecarlo_configuration.LINE_INTERACTION_TYPE = (
            LineInteractionType.DOWNBRANCH
        )
    elif transport.line_interaction_type == "scatter":
        montecarlo_configuration.LINE_INTERACTION_TYPE = (
            LineInteractionType.SCATTER
        )
    else:
        raise ValueError(
            f'Line interaction type must be one of "macroatom",'
            f'"downbranch", or "scatter" but is '
            f"{transport.line_interaction_type}"
        )
    montecarlo_configuration.NUMBER_OF_VPACKETS = number_of_vpackets
    montecarlo_configuration.TEMPORARY_V_PACKET_BINS = number_of_vpackets
    montecarlo_configuration.ENABLE_FULL_RELATIVITY = (
        transport.enable_full_relativity
    )
    montecarlo_configuration.MONTECARLO_SEED = transport.packet_source.base_seed
    montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY = (
        transport.virtual_spectrum_spawn_range.end.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY = (
        transport.virtual_spectrum_spawn_range.start.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.ENABLE_VPACKET_TRACKING = (
        transport.enable_vpacket_tracking
    )
    montecarlo_configuration.ENABLE_RPACKET_TRACKING = (
        transport.enable_rpacket_tracking
    )
