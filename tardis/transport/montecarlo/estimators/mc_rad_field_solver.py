import numpy as np
from astropy import units as u
from astropy.units import Quantity
from scipy.special import zeta

from tardis import constants as const
from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)
from tardis.transport.montecarlo.estimators.base import (
    EstimatedRadiationFieldProperties,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
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

    def __init__(self, w_epsilon: float = 1e-10) -> None:
        self.w_epsilon = w_epsilon

    def solve(
        self,
        estimators_bulk: EstimatorsBulk,
        estimators_line: EstimatorsLine,
        time_explosion: Quantity,
        time_of_simulation: Quantity,
        volume: np.ndarray,
        line_list_nu: np.ndarray,
    ) -> EstimatedRadiationFieldProperties:
        """
        Calculate an updated radiation field from the :math:
        `\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation.
        The details of the calculation can be found in the documentation.

        Parameters
        ----------
        estimators_bulk
            Bulk radiation field estimators
        estimators_line
            Line interaction estimators
        time_explosion
            Time since explosion
        time_of_simulation
            Time of simulation
        volume
            Volume of each cell
        line_list_nu
            Frequency list for lines

        Returns
        -------
        Radiation field properties including t_radiative and dilution_factor
        """
        dilute_planck_rad_field = self.estimate_dilute_planck_radiation_field(
            estimators_bulk, time_of_simulation, volume
        )
        j_blues = self.estimate_jblues(
            estimators_line.mean_intensity_blueward,
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
        self,
        estimators_bulk: EstimatorsBulk,
        time_of_simulation: Quantity,
        volume: np.ndarray,
    ) -> DilutePlanckianRadiationField:
        temperature_radiative = (
            T_RADIATIVE_ESTIMATOR_CONSTANT
            * estimators_bulk.mean_frequency
            / estimators_bulk.mean_intensity_total
        ) * u.K
        dilution_factor = estimators_bulk.mean_intensity_total / (
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
        j_blue_estimator: np.ndarray,
        estimated_radfield_state: DilutePlanckianRadiationField,
        time_explosion: Quantity,
        time_of_simulation: Quantity,
        volume: np.ndarray,
        line_list_nu: np.ndarray,
    ) -> np.ndarray:
        j_blues_norm_factor = (
            const.c.cgs
            * time_explosion
            / (4 * np.pi * time_of_simulation * volume)
        )
        j_blues = j_blue_estimator * j_blues_norm_factor.cgs.value
        planck_j_blues = estimated_radfield_state.calculate_mean_intensity(
            line_list_nu
        )
        zero_j_blues = j_blues == 0.0
        j_blues[zero_j_blues] = self.w_epsilon * planck_j_blues[zero_j_blues]

        return j_blues
