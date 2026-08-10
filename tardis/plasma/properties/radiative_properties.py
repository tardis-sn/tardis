import logging

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.opacities.macro_atom.base import TransitionProbabilities
from tardis.opacities.sobolev import calculate_stimulated_emission_factor
from tardis.plasma.properties.base import (
    ProcessingPlasmaProperty,
    TransitionProbabilitiesProperty,
)

logger = logging.getLogger(__name__)

__all__ = [
    "RawRadBoundBoundTransProbs",
    "StimulatedEmissionFactor",
    "calculate_stimulated_emission_factor",
]

C_EINSTEIN = (
    4.0 * (np.pi * const.e.esu) ** 2 / (const.c.cgs * const.m_e.cgs)
).value  # See tardis/docs/physics/plasma/macroatom.rst


class StimulatedEmissionFactor(ProcessingPlasmaProperty):
    """Calculate stimulated-emission factors for atomic lines.

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
        """Cache statistical weights for lower line levels."""
        if self._g_lower is None:
            g_lower = np.array(
                g.iloc[lines_lower_level_index], dtype=np.float64
            )
            self._g_lower = g_lower[np.newaxis].T
        return self._g_lower

    def get_g_upper(self, g, lines_upper_level_index):
        """Cache statistical weights for upper line levels."""
        if self._g_upper is None:
            g_upper = np.array(
                g.iloc[lines_upper_level_index], dtype=np.float64
            )
            self._g_upper = g_upper[np.newaxis].T
        return self._g_upper

    def get_metastable_upper(self, metastability, lines_upper_level_index):
        """Cache metastability flags for upper line levels."""
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
        """Calculate stimulated-emission factors for the supplied state."""
        return calculate_stimulated_emission_factor(
            g,
            level_number_density,
            lines_lower_level_index,
            lines_upper_level_index,
            metastability,
            lines,
            self.nlte_species,
        )


class RawRadBoundBoundTransProbs(
    TransitionProbabilities, TransitionProbabilitiesProperty
):
    """Calculate unnormalized radiative bound-bound transition probabilities.

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
        """Calculate radiative transition probabilities for active species."""
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
