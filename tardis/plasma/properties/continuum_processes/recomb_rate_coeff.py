import numpy as np
import pandas as pd

from tardis.plasma.properties.base import Input, ProcessingPlasmaProperty
from tardis.plasma.properties.continuum_processes.rates import C, H
from tardis.transport.montecarlo.estimators.util import (
    bound_free_estimator_array2frame,
    integrate_array_by_blocks,
)


class StimRecombRateFactor(Input):
    """
    Attributes
    ----------
    alpha_stim_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the rate coefficient for stimulated
        recombination.
    """

    outputs = ("alpha_stim_factor",)
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
        alpha_stim_factor,
        phi_ik,
    ):
        return alpha_stim_factor * phi_ik.loc[alpha_stim_factor.index]


class SpontRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_sp : pandas.DataFrame, dtype float
        The rate coefficient for spontaneous recombination.
    """

    outputs = ("alpha_sp",)
    latex_name = (r"\alpha^{\textrm{sp}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        photo_ion_block_references,
        photo_ion_index,
        phi_ik,
        boltzmann_factor_photo_ion,
    ):
        cross_section = photo_ion_cross_sections["x_sect"].values
        nu = photo_ion_cross_sections["nu"].values

        alpha_sp = 8 * np.pi * cross_section * nu**2 / C**2
        alpha_sp = alpha_sp[:, np.newaxis]
        alpha_sp = alpha_sp * boltzmann_factor_photo_ion
        alpha_sp = integrate_array_by_blocks(
            alpha_sp, nu, photo_ion_block_references
        )
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik.loc[alpha_sp.index]
