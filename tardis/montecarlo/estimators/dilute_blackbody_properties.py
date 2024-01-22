import numpy as np
from astropy import units as u
from scipy.special import zeta

from tardis import constants as const
from tardis.model.radiation_field_state import (
    DiluteBlackBodyRadiationFieldState,
)

DILUTION_FACTOR_ESTIMATOR_CONSTANT = (
    (const.c**2 / (2 * const.h))
    * (15 / np.pi**4)
    * (const.h / const.k_B) ** 4
    / (4 * np.pi)
).cgs.value
T_RADIATIVE_ESTIMATOR_CONSTANT = (
    (np.pi**4 / (15 * 24 * zeta(5, 1))) * (const.h / const.k_B)
).cgs.value


class MCDiluteBlackBodyRadFieldSolver:
    def __init__(self) -> None:
        pass

    def solve(self, radfield_mc_estimators, time_of_simulation, volume):
        """
        Calculate an updated radiation field from the :math:
        `\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation.
        The details of the calculation can be found in the documentation.

        Parameters
        ----------
        nubar_estimator : np.ndarray (float)
        j_estimator : np.ndarray (float)

        Returns
        -------
        t_radiative : astropy.units.Quantity (float)
        dilution_factor : numpy.ndarray (float)
        """
        temperature_radiative = (
            T_RADIATIVE_ESTIMATOR_CONSTANT
            * radfield_mc_estimators.nu_bar_estimator
            / radfield_mc_estimators.j_estimator
        ) * u.K
        dilution_factor = radfield_mc_estimators.j_estimator / (
            4
            * const.sigma_sb.cgs.value
            * temperature_radiative.value**4
            * time_of_simulation.value
            * volume
        )

        return DiluteBlackBodyRadiationFieldState(
            temperature_radiative, dilution_factor
        )
