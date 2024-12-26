import numpy as np
from astropy import units as u
from scipy.special import zeta

from tardis import constants as const
from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)
from tardis.transport.montecarlo.estimators.base import (
    EstimatedRadiationFieldProperties,
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


class MCRadiationFieldPropertiesSolver:
    w_epsilon = 1e-10

    def __init__(self, w_epsilon=1e-10) -> None:
        self.w_epsilon = w_epsilon

    def solve(
        self,
        radfield_mc_estimators,
        time_explosion,
        time_of_simulation,
        volume,
        line_list_nu,
    ):
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
        dilute_planck_rad_field = self.estimate_dilute_planck_radiation_field(
            radfield_mc_estimators, time_of_simulation, volume
        )
        j_blues = self.estimate_jblues(
            radfield_mc_estimators.j_blue_estimator,
            dilute_planck_rad_field,
            time_explosion,
            time_of_simulation,
            volume,
            line_list_nu,
        )

        return EstimatedRadiationFieldProperties(
            dilute_blackbody_radiationfield_state=dilute_planck_rad_field,
            j_blues=j_blues,
        )

    def estimate_dilute_planck_radiation_field(
        self, radfield_mc_estimators, time_of_simulation, volume
    ):
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
        return DilutePlanckianRadiationField(
            temperature_radiative, dilution_factor
        )

    def estimate_jblues(
        self,
        j_blue_estimator,
        estimated_radfield_state,
        time_explosion,
        time_of_simulation,
        volume,
        line_list_nu,
    ):
        print("[J_BLUE_DEBUG] Input j_blue_estimator:", j_blue_estimator)
        
        j_blues_norm_factor = (
            const.c.cgs
            * time_explosion
            / (4 * np.pi * time_of_simulation * volume)
        )
        print("[J_BLUE_DEBUG] Normalization factor:", j_blues_norm_factor.cgs.value)
        
        j_blues = j_blue_estimator * j_blues_norm_factor.cgs.value
        print("[J_BLUE_DEBUG] After normalization:", j_blues)
        print("[J_BLUE_DEBUG] line_list_nu:", line_list_nu)
        
        planck_j_blues = estimated_radfield_state.calculate_mean_intensity(
            line_list_nu
        )
        print("[J_BLUE_DEBUG] Planck j_blues:", planck_j_blues)
        
        zero_j_blues = j_blues == 0.0
        j_blues[zero_j_blues] = self.w_epsilon * planck_j_blues[zero_j_blues]
        print("[J_BLUE_DEBUG] Final j_blues after zero handling:", j_blues)
        
        return j_blues
