import numpy as np
from astropy import units as u

from tardis.transport.montecarlo.packet_source import BasePacketSource
from tardis.opacities.opacity_state import OpacityState


class MonteCarloRadiationFieldState:
    """_summary_

    Parameters
    ----------
    t_radiative : u.Quantity
        Radiative temperature in each shell
    dilution_factor : numpy.ndarray
        Dilution Factors in each shell
    opacities : OpacityState
        Opacity container object
    packet_source : SourceFunction
        Source function for radiative transfer, for example a packet_source
    """

    def __init__(
        self,
        t_radiative: u.Quantity,
        dilution_factor: np.ndarray,
        opacities: OpacityState,
        packet_source: BasePacketSource,
    ):
        self.t_radiative = t_radiative
        self.dilution_factor = dilution_factor
        self.t_rad = self.t_radiative
        self.opacities = opacities
        self.packet_source = packet_source
