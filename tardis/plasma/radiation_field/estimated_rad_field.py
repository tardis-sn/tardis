from typing import Union

import numpy as np
from astropy import units as u


class EstimatedRadiationField:
    """
    Represents the state of an estimated radiation field.


    Parameters
    ----------
    j_blues : numpy.ndarray
        J_blues in each shell, which is the estimator for the mean intensity of the radiation field
    """

    def __init__(self, j_blues):
        self.j_blues = j_blues

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
        return self.j_blues
