import numpy as np
import pandas as pd
from tardis.plasma.properties.base import Input
from tardis.plasma.properties.continuum_processes.rates import H
from tardis.transport.montecarlo.estimators.util import (
    bound_free_estimator_array2frame,
    integrate_array_by_blocks,
    ProcessingPlasmaProperty,
)


class StimRecombRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    alpha_stim_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the rate coefficient for stimulated
        recombination.
    """

    outputs = ("alpha_stim_estimator",)
    latex_name = (r"\alpha^{\textrm{stim}}_\textrm{estim}",)


class StimRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_stim : pandas.DataFrame, dtype float
        The rate coefficient for stimulated recombination.
    """

    outputs = ("alpha_stim",)
    latex_name = (r"\alpha^{\textrm{stim}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        alpha_stim_estimator,
        photo_ion_norm_factor,
        photo_ion_block_references,
        photo_ion_index,
        dilute_planckian_radiation_field,
        phi_ik,
        t_electrons,
        boltzmann_factor_photo_ion,
        level2continuum_idx,
    ):
        # Used for initialization
        if alpha_stim_estimator is None:
            alpha_stim = self.calculate_from_dilute_bb(
                photo_ion_cross_sections,
                photo_ion_block_references,
                photo_ion_index,
                dilute_planckian_radiation_field,
                t_electrons,
                boltzmann_factor_photo_ion,
            )
        else:
            alpha_stim_estimator = bound_free_estimator_array2frame(
                alpha_stim_estimator, level2continuum_idx
            )
            alpha_stim = alpha_stim_estimator * photo_ion_norm_factor
        alpha_stim *= phi_ik.loc[alpha_stim.index]
        return alpha_stim

    @staticmethod
    def calculate_from_dilute_bb(
        photo_ion_cross_sections,
        photo_ion_block_references,
        photo_ion_index,
        dilute_planckian_radiation_field,
        t_electrons,
        boltzmann_factor_photo_ion,
    ):
        nu = photo_ion_cross_sections["nu"]
        x_sect = photo_ion_cross_sections["x_sect"]
        j_nus = dilute_planckian_radiation_field.calculate_mean_intensity(nu)
        j_nus *= boltzmann_factor_photo_ion
        alpha_stim = j_nus.multiply(4.0 * np.pi * x_sect / nu / H, axis=0)
        alpha_stim = integrate_array_by_blocks(
            alpha_stim.values, nu.values, photo_ion_block_references
        )
        alpha_stim = pd.DataFrame(alpha_stim, index=photo_ion_index)
        return alpha_stim
