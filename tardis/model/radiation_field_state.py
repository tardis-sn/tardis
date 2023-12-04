from tardis.montecarlo.montecarlo_numba.numba_interface import OpacityState


import numpy as np
from astropy import units as u


class DiluteThermalRadiationFieldState:
    """
    Represents the state of a dilute thermal radiation field.


    Parameters
    ----------
    t_radiative : u.Quantity
        Radiative temperature in each shell
    dilution_factor : numpy.ndarray
        Dilution Factors in each shell
    """

    def __init__(
        self,
        t_radiative: u.Quantity,
        dilution_factor: np.ndarray,
    ):
        # ensuring that the radiation_field has both
        # dilution_factor and t_radiative equal length
        assert len(t_radiative) == len(dilution_factor)
        assert np.all(t_radiative > 0 * u.K)
        assert np.all(dilution_factor > 0)
        self.t_radiative = t_radiative
        self.dilution_factor = dilution_factor
