import numpy as np
import pandas as pd

from tardis.plasma.properties.base import Input
from tardis.plasma.properties.continuum_processes.rates import H
from tardis.transport.montecarlo.estimators.util import (
    ProcessingPlasmaProperty,
    bound_free_estimator_array2frame,
    integrate_array_by_blocks,
)


class PhotoIonRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    gamma_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the rate coefficient for radiative
        ionization.
    """

    outputs = ("gamma_estimator",)
    latex_name = (r"\gamma_\textrm{estim}",)


class PhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : pandas.DataFrame, dtype float
        The rate coefficient for radiative ionization.
    """

    outputs = ("gamma",)
    latex_name = (r"\gamma",)

    def calculate(
        self,
        photo_ion_cross_sections,
        gamma_estimator,
        photo_ion_norm_factor,
        photo_ion_block_references,
        photo_ion_index,
        dilute_planckian_radiation_field,
        level2continuum_idx,
    ):
        # Used for initialization
        if gamma_estimator is None:
            gamma = self.calculate_from_dilute_bb(
                photo_ion_cross_sections,
                photo_ion_block_references,
                photo_ion_index,
                dilute_planckian_radiation_field,
            )
        else:
            gamma_estimator = bound_free_estimator_array2frame(
                gamma_estimator, level2continuum_idx
            )
            gamma = gamma_estimator * photo_ion_norm_factor.value

        return gamma

    @staticmethod
    def calculate_from_dilute_bb(
        photo_ion_cross_sections,
        photo_ion_block_references,
        photo_ion_index,
        dilute_planckian_radiation_field,
    ):
        nu = photo_ion_cross_sections["nu"]
        x_sect = photo_ion_cross_sections["x_sect"]
        j_nus = dilute_planckian_radiation_field.calculate_mean_intensity(
            nu,
        )
        gamma = j_nus.multiply(4.0 * np.pi * x_sect / nu / H, axis=0)
        gamma = integrate_array_by_blocks(
            gamma.values, nu.values, photo_ion_block_references
        )
        gamma = pd.DataFrame(gamma, index=photo_ion_index)
        return gamma
