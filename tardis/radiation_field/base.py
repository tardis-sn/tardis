import numpy as np
from astropy import units as u

from tardis.montecarlo.packet_source import BasePacketSource
from tardis.montecarlo.montecarlo_numba.numba_interface import OpacityState


class RadiationField:
    """_summary_

    Parameters
    ----------
    t_rad : numpy.ndarray
        Radiative temperature in each shell
    w : numpy.ndarray
        Dilution Factors in each shell
    opacities : OpacityState
        Opacity container object
    source_function : SourceFunction
        Source function for radiative transfer, for example a packet_source
    """

    def __init__(self, t_rad, w, opacities, source_function):
        self.t_rad = t_rad
        self.w = w
        self.opacities = opacities
        self.source_function = source_function


class MonteCarloRadiationFieldState:
    """_summary_

    Parameters
    ----------
    t_rad : numpy.ndarray
        Radiative temperature in each shell
    w : numpy.ndarray
        Dilution Factors in each shell
    opacities : OpacityState
        Opacity container object
    packet_source : SourceFunction
        Source function for radiative transfer, for example a packet_source
    """

    def __init__(
        self,
        t_radiative: np.ndarray,
        dilution_factor: np.ndarray,
        opacities: OpacityState,
        packet_source: BasePacketSource,
    ):
        self.t_radiative = t_radiative
        self.dilution_factor = dilution_factor
        self.t_rad = self.t_radiative
        self.w = self.dilution_factor
        self.opacities = opacities
        self.source_function = packet_source


def convert_config_to_radiationfield_state(
    config, no_of_shells, temperature=None
):  ### make it a to_method in the end
    if temperature:
        t_radiative = temperature
    elif config.plasma.initial_t_rad > 0 * u.K:
        t_radiative = np.ones(no_of_shells + 1) * config.plasma.initial_t_rad
    else:
        t_radiative = None

    return MonteCarloRadiationFieldState(t_radiative, None, None, None)
