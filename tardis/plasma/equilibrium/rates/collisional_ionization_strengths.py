import astropy.units as u
import numpy as np
import pandas as pd

from tardis import constants as const

H = const.h.cgs
K_B = const.k_B.cgs


class CollisionalIonizationSeaton:
    """Solver for collisional ionization rate coefficients in the Seaton approximation."""

    def __init__(self, photoionization_cross_sections):
        self.photoionization_cross_sections = photoionization_cross_sections

    def solve(self, electron_temperature):
        """
        Parameters
        ----------
        electron_temperature : u.Quantity
            The electron temperature in K.

        Returns
        -------
        pandas.DataFrame, dtype float
            The rate coefficient for collisional ionization in the Seaton
            approximation. Multiply with the electron density and the
            level number density to obtain the total rate.

        Notes
        -----
        The rate coefficient for collisional ionization in the Seaton approximation
        is calculated according to Eq. 9.60 in [1].

        References
        ----------
        .. [1] Hubeny, I. and Mihalas, D., "Theory of Stellar Atmospheres". 2014.
        """
        photo_ion_cross_sections_threshold = (
            self.photoionization_cross_sections.groupby(level=[0, 1, 2]).first()
        )
        nu_i = photo_ion_cross_sections_threshold["nu"]
        u0s = (
            nu_i.values[np.newaxis].T * u.Hz / electron_temperature * (H / K_B)
        )
        factor = np.exp(-u0s) / u0s
        factor = pd.DataFrame(factor, index=nu_i.index)
        coll_ion_coeff = 1.55e13 * photo_ion_cross_sections_threshold["x_sect"]
        coll_ion_coeff = factor.multiply(coll_ion_coeff, axis=0)
        coll_ion_coeff = coll_ion_coeff.divide(
            np.sqrt(electron_temperature), axis=1
        )

        ion_number = coll_ion_coeff.index.get_level_values("ion_number").values
        coll_ion_coeff[ion_number == 0] *= 0.1
        coll_ion_coeff[ion_number == 1] *= 0.2
        coll_ion_coeff[ion_number >= 2] *= 0.3
        return coll_ion_coeff
