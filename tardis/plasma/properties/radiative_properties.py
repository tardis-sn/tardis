import logging

import numpy as np
import pandas as pd
from astropy import units as u
from numba import jit, prange

from tardis import constants as const
from tardis.opacities.macro_atom.base import TransitionProbabilities
from tardis.plasma.properties.base import (
    ProcessingPlasmaProperty,
    TransitionProbabilitiesProperty,
)

logger = logging.getLogger(__name__)

__all__ = [
    "StimulatedEmissionFactor",
    "BetaSobolev",
    "RawRadBoundBoundTransProbs",
]

C_EINSTEIN = (
    4.0 * (np.pi * const.e.esu) ** 2 / (const.c.cgs * const.m_e.cgs)
).value  # See tardis/docs/physics/plasma/macroatom.rst


class StimulatedEmissionFactor(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    stimulated_emission_factor : Numpy Array, dtype float
         Indexed by lines, columns as zones.
    """

    outputs = ("stimulated_emission_factor",)
    latex_formula = (r"1-\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}",)

    def __init__(self, plasma_parent=None, nlte_species=None):
        super(StimulatedEmissionFactor, self).__init__(plasma_parent)
        self._g_upper = None
        self._g_lower = None
        self.nlte_species = nlte_species

    def get_g_lower(self, g, lines_lower_level_index):
        if self._g_lower is None:
            g_lower = np.array(
                g.iloc[lines_lower_level_index], dtype=np.float64
            )
            self._g_lower = g_lower[np.newaxis].T
        return self._g_lower

    def get_g_upper(self, g, lines_upper_level_index):
        if self._g_upper is None:
            g_upper = np.array(
                g.iloc[lines_upper_level_index], dtype=np.float64
            )
            self._g_upper = g_upper[np.newaxis].T
        return self._g_upper

    def get_metastable_upper(self, metastability, lines_upper_level_index):
        if getattr(self, "_meta_stable_upper", None) is None:
            self._meta_stable_upper = metastability.values[
                lines_upper_level_index
            ][np.newaxis].T
        return self._meta_stable_upper

    def calculate(
        self,
        g,
        level_number_density,
        lines_lower_level_index,
        lines_upper_level_index,
        metastability,
        lines,
    ):
        n_lower = level_number_density.values.take(
            lines_lower_level_index, axis=0, mode="raise"
        )
        n_upper = level_number_density.values.take(
            lines_upper_level_index, axis=0, mode="raise"
        )
        g_lower = self.get_g_lower(g, lines_lower_level_index)
        g_upper = self.get_g_upper(g, lines_upper_level_index)
        meta_stable_upper = self.get_metastable_upper(
            metastability, lines_upper_level_index
        )

        stimulated_emission_factor = 1 - (
            (g_lower * n_upper) / (g_upper * n_lower)
        )
        stimulated_emission_factor[n_lower == 0.0] = 0.0
        stimulated_emission_factor[
            np.isneginf(stimulated_emission_factor)
        ] = 0.0
        stimulated_emission_factor[
            meta_stable_upper & (stimulated_emission_factor < 0)
        ] = 0.0
        if self.nlte_species:
            nlte_lines_mask = (
                lines.reset_index()
                .apply(
                    lambda row: (row.atomic_number, row.ion_number)
                    in self.nlte_species,
                    axis=1,
                )
                .values
            )
            stimulated_emission_factor[
                (stimulated_emission_factor < 0) & nlte_lines_mask[np.newaxis].T
            ] = 0.0
        return stimulated_emission_factor


class BetaSobolev(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    beta_sobolev : Numpy Array, dtype float
    """

    outputs = ("beta_sobolev",)
    latex_name = (r"\beta_{\textrm{sobolev}}",)

    def calculate(self, tau_sobolevs):
        if getattr(self, "beta_sobolev", None) is None:
            initial = 0.0
        else:
            initial = self.beta_sobolev

        beta_sobolev = pd.DataFrame(
            initial, index=tau_sobolevs.index, columns=tau_sobolevs.columns
        )

        self.calculate_beta_sobolev(
            tau_sobolevs.values.ravel(), beta_sobolev.values.ravel()
        )
        return beta_sobolev

    @staticmethod
    @jit(nopython=True, parallel=True)
    def calculate_beta_sobolev(tau_sobolevs, beta_sobolevs):
        for i in prange(len(tau_sobolevs)):
            if tau_sobolevs[i] > 1e3:
                beta_sobolevs[i] = tau_sobolevs[i] ** -1
            elif tau_sobolevs[i] < 1e-4:
                beta_sobolevs[i] = 1 - 0.5 * tau_sobolevs[i]
            else:
                beta_sobolevs[i] = (1 - np.exp(-tau_sobolevs[i])) / (
                    tau_sobolevs[i]
                )
        return beta_sobolevs


class RawRadBoundBoundTransProbs(
    TransitionProbabilities, TransitionProbabilitiesProperty
):
    """
    Attributes
    ----------
    p_rad_bb : pandas.DataFrame, dtype float
        Unnormalized transition probabilities for radiative bound-bound
        transitions
    """

    outputs = ("p_rad_bb",)
    transition_probabilities_outputs = ("p_rad_bb",)

    def __init__(self, plasma_parent):
        super(RawRadBoundBoundTransProbs, self).__init__(plasma_parent)
        self.normalize = False

    def calculate(
        self,
        atomic_data,
        beta_sobolev,
        j_blues,
        stimulated_emission_factor,
        tau_sobolevs,
        continuum_interaction_species,
    ):
        p_rad_bb = super().calculate(
            atomic_data,
            beta_sobolev,
            j_blues,
            stimulated_emission_factor,
            tau_sobolevs,
        )
        transition_type = atomic_data.macro_atom_data.transition_type.replace(
            1, 0
        )
        index = pd.MultiIndex.from_arrays(
            [
                atomic_data.macro_atom_data.source_level_idx,
                atomic_data.macro_atom_data.destination_level_idx,
                transition_type,
            ]
        )
        mask_continuum_species = pd.MultiIndex.from_arrays(
            [
                atomic_data.macro_atom_data.atomic_number,
                atomic_data.macro_atom_data.ion_number,
            ]
        ).isin(continuum_interaction_species)
        p_rad_bb = p_rad_bb.set_index(index, drop=True)[mask_continuum_species]
        # To obtain energy-flow rates in cgs from the precomputed transition
        # probabilities in the atomic data, we have to multiply by the
        # constant C_EINSTEIN and convert from eV to erg.
        # See tardis/docs/physics/plasma/macroatom.rst
        p_rad_bb = p_rad_bb * C_EINSTEIN * u.eV.to(u.erg)
        return p_rad_bb
