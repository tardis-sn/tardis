import numpy as np
from astropy import units as u

from tardis.util.base import intensity_black_body

from typing import Union


class DiluteBlackBodyRadiationFieldState:
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

    def calculate_mean_intensity(self, nu: Union[u.Quantity, np.ndarray]):
        """
        Calculate the intensity of the radiation field at a given frequency.

        Parameters
        ----------
        nu : u.Quantity
            Frequency at which the intensity is to be calculated

        Returns
        -------
        intensity : u.Quantity
            Intensity of the radiation field at the given frequency
        """
        return self.dilution_factor * intensity_black_body(
            nu[np.newaxis].T, self.t_radiative
        )
