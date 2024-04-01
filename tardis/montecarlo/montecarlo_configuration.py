from astropy import units as u
import functools
from numba import float64, int64, boolean, jit
from numba.experimental import jitclass
import numpy as np

from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType,
)


class MonteCarloConfiguration:
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


@functools.cache
def obj2strkeydict(obj, config_name):
    # unpack object to freevars and close over them
    tmp_ENABLE_FULL_RELATIVITY = obj.ENABLE_FULL_RELATIVITY
    tmp_TEMPORARY_V_PACKET_BINS = obj.TEMPORARY_V_PACKET_BINS
    tmp_NUMBER_OF_VPACKETS = obj.NUMBER_OF_VPACKETS
    tmp_MONTECARLO_SEED = obj.MONTECARLO_SEED
    tmp_LINE_INTERACTION_TYPE = obj.LINE_INTERACTION_TYPE
    tmp_PACKET_SEEDS = obj.PACKET_SEEDS
    tmp_DISABLE_ELECTRON_SCATTERING = obj.DISABLE_ELECTRON_SCATTERING
    tmp_DISABLE_LINE_SCATTERING = obj.DISABLE_LINE_SCATTERING
    tmp_SURVIVAL_PROBABILITY = obj.SURVIVAL_PROBABILITY
    tmp_VPACKET_TAU_RUSSIAN = obj.VPACKET_TAU_RUSSIAN
    tmp_INITIAL_TRACKING_ARRAY_LENGTH = obj.INITIAL_TRACKING_ARRAY_LENGTH
    tmp_LEGACY_MODE_ENABLED = obj.LEGACY_MODE_ENABLED
    tmp_ENABLE_RPACKET_TRACKING = obj.ENABLE_RPACKET_TRACKING
    tmp_CONTINUUM_PROCESSES_ENABLED = obj.CONTINUUM_PROCESSES_ENABLED
    tmp_VPACKET_SPAWN_START_FREQUENCY = obj.VPACKET_SPAWN_START_FREQUENCY
    tmp_VPACKET_SPAWN_END_FREQUENCY = obj.VPACKET_SPAWN_END_FREQUENCY
    tmp_ENABLE_VPACKET_TRACKING = obj.ENABLE_VPACKET_TRACKING

    assert isinstance(config_name, str)
    tmp_force_heterogeneous = config_name

    @jit
    def configurator():
        dict = {
            "ENABLE_FULL_RELATIVITY": tmp_ENABLE_FULL_RELATIVITY,
            "TEMPORARY_V_PACKET_BINS": tmp_TEMPORARY_V_PACKET_BINS,
            "NUMBER_OF_VPACKETS": tmp_NUMBER_OF_VPACKETS,
            "MONTECARLO_SEED": tmp_MONTECARLO_SEED,
            "LINE_INTERACTION_TYPE": tmp_LINE_INTERACTION_TYPE,
            "PACKET_SEEDS": tmp_PACKET_SEEDS,
            "DISABLE_ELECTRON_SCATTERING": tmp_DISABLE_ELECTRON_SCATTERING,
            "DISABLE_LINE_SCATTERING": tmp_DISABLE_LINE_SCATTERING,
            "SURVIVAL_PROBABILITY": tmp_SURVIVAL_PROBABILITY,
            "VPACKET_TAU_RUSSIAN": tmp_VPACKET_TAU_RUSSIAN,
            "INITIAL_TRACKING_ARRAY_LENGTH": tmp_INITIAL_TRACKING_ARRAY_LENGTH,
            "LEGACY_MODE_ENABLED": tmp_LEGACY_MODE_ENABLED,
            "ENABLE_RPACKET_TRACKING": tmp_ENABLE_RPACKET_TRACKING,
            "CONTINUUM_PROCESSES_ENABLED": tmp_CONTINUUM_PROCESSES_ENABLED,
            "VPACKET_SPAWN_START_FREQUENCY": tmp_VPACKET_SPAWN_START_FREQUENCY,
            "VPACKET_SPAWN_END_FREQUENCY": tmp_VPACKET_SPAWN_END_FREQUENCY,
            "ENABLE_VPACKET_TRACKING": tmp_ENABLE_VPACKET_TRACKING,
            "config_name": tmp_force_heterogeneous,
        }
        return dict

    # return a configuration function that returns a string-key-dict
    # representation of the configuration object.
    return configurator
