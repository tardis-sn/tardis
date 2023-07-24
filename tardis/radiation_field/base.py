import numpy as np
from astropy import units as u

from tardis.montecarlo.packet_source import BasePacketSource
from tardis.montecarlo.montecarlo_numba.numba_interface import OpacityState


class RadiationField:
    """_summary_

    Parameters
    ----------
    t_radiative : numpy.ndarray
        Radiative temperature in each shell
    dilution_factor : numpy.ndarray
        Dilution Factors in each shell
    opacities : OpacityState
        Opacity container object
    source_function : SourceFunction
        Source function for radiative transfer, for example a packet_source
    """

    def __init__(self, t_radiative, dilution_factor, opacities, source_function):
        self.t_radiative = t_radiative
        self.dilution_factor = dilution_factor
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
        self.opacities = opacities
        self.source_function = packet_source


def convert_config_to_radiation_field_state(
    config
):  ### make it a to_method in the end
    if t_radiative is None:
        lambda_wien_inner = (
            constants.b_wien / self.blackbody_packet_source.temperature
        )
        self._t_radiative = constants.b_wien / (
            lambda_wien_inner
            * (1 + (self.v_middle - self.v_boundary_inner) / constants.c)
        )
    else:
        # self._t_radiative = t_radiative[self.v_boundary_inner_index + 1:self.v_boundary_outer_index]
        self._t_radiative = t_radiative

    if dilution_factor is None:
        self._dilution_factor = 0.5 * (
            1
            - np.sqrt(
                1 - (self.r_inner[0] ** 2 / self.r_middle**2).to(1).value
            )
        )
    else:
        # self.dilution_factor = dilution_factor[self.v_boundary_inner_index + 1:self.v_boundary_outer_index]
        self._dilution_factor = dilution_factor
    return MonteCarloRadiationFieldState(t_radiative, None, None, None)
