import logging

import numpy as np
import pandas as pd
from astropy import units as u

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
        super().__init__(plasma_parent)
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

        # In theory the factor should be 1 for n_lower = 0, but in practice the opacity is reduced to 0 anyway
        stimulated_emission_factor = np.zeros(n_lower.shape, dtype=np.float64)

        n_lower_zero_mask = n_lower == 0.0
        stimulated_emission_factor[~n_lower_zero_mask] = 1 - (
            (g_lower * n_upper)[~n_lower_zero_mask]
            / (g_upper * n_lower)[~n_lower_zero_mask]
        )

        # the following line probably can be removed as well
        stimulated_emission_factor[np.isneginf(stimulated_emission_factor)] = (
            0.0
        )
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
        super().__init__(plasma_parent)
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
