import logging

import numpy as np
import pandas as pd
import numexpr as ne
from astropy import units as u
from tardis import constants as const
from numba import jit, prange

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.properties.util import macro_atom

logger = logging.getLogger(__name__)

__all__ = [
    "StimulatedEmissionFactor",
    "TauSobolev",
    "BetaSobolev",
    "TransitionProbabilities",
]


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

        stimulated_emission_factor = ne.evaluate(
            "1 - ((g_lower * n_upper) / " "(g_upper * n_lower))"
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


class TauSobolev(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    tau_sobolev : Pandas DataFrame, dtype float
                  Sobolev optical depth for each line. Indexed by line.
                  Columns as zones.
    """

    outputs = ("tau_sobolevs",)
    latex_name = (r"\tau_{\textrm{sobolev}}",)
    latex_formula = (
        r"\dfrac{\pi e^{2}}{m_{e} c}f_{lu}\lambda t_{exp}\
        n_{lower} \Big(1-\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}\Big)",
    )

    def __init__(self, plasma_parent):
        super(TauSobolev, self).__init__(plasma_parent)
        self.sobolev_coefficient = (
            (
                ((np.pi * const.e.gauss ** 2) / (const.m_e.cgs * const.c.cgs))
                * u.cm
                * u.s
                / u.cm ** 3
            )
            .to(1)
            .value
        )

    def calculate(
        self,
        lines,
        level_number_density,
        lines_lower_level_index,
        time_explosion,
        stimulated_emission_factor,
        j_blues,
        f_lu,
        wavelength_cm,
    ):
        f_lu = f_lu.values[np.newaxis].T
        wavelength = wavelength_cm.values[np.newaxis].T
        n_lower = level_number_density.values.take(
            lines_lower_level_index, axis=0, mode="raise"
        )
        tau_sobolevs = (
            self.sobolev_coefficient
            * f_lu
            * wavelength
            * time_explosion
            * n_lower
            * stimulated_emission_factor
        )

        if np.any(np.isnan(tau_sobolevs)) or np.any(
            np.isinf(np.abs(tau_sobolevs))
        ):
            raise ValueError(
                "Some tau_sobolevs are nan, inf, -inf in tau_sobolevs."
                " Something went wrong!"
            )

        return pd.DataFrame(
            tau_sobolevs,
            index=lines.index,
            columns=np.array(level_number_density.columns),
        )


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


class TransitionProbabilities(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    transition_probabilities : Pandas DataFrame, dtype float
    """

    outputs = ("transition_probabilities",)

    def __init__(self, plasma_parent):
        super(TransitionProbabilities, self).__init__(plasma_parent)
        self.initialize = True

    def calculate(
        self,
        atomic_data,
        beta_sobolev,
        j_blues,
        stimulated_emission_factor,
        tau_sobolevs,
    ):
        # I wonder why?
        # Not sure who wrote this but the answer is that when the plasma is
        # first initialised (before the first iteration, without temperature
        # values etc.) there are no j_blues values so this just prevents
        # an error. Aoife.
        if len(j_blues) == 0:
            return None
        macro_atom_data = self._get_macro_atom_data(atomic_data)
        if self.initialize:
            self.initialize_macro_atom_transition_type_filters(
                atomic_data, macro_atom_data
            )
            self.transition_probability_coef = (
                self._get_transition_probability_coefs(macro_atom_data)
            )
            self.initialize = False
        transition_probabilities = self._calculate_transition_probability(
            macro_atom_data, beta_sobolev, j_blues, stimulated_emission_factor
        )
        transition_probabilities = pd.DataFrame(
            transition_probabilities,
            index=macro_atom_data.transition_line_id,
            columns=tau_sobolevs.columns,
        )
        return transition_probabilities

    def _calculate_transition_probability(
        self, macro_atom_data, beta_sobolev, j_blues, stimulated_emission_factor
    ):
        transition_probabilities = np.empty(
            (self.transition_probability_coef.shape[0], beta_sobolev.shape[1])
        )
        # trans_old = self.calculate_transition_probabilities(macro_atom_data, beta_sobolev, j_blues, stimulated_emission_factor)
        transition_type = macro_atom_data.transition_type.values
        lines_idx = macro_atom_data.lines_idx.values
        tpos = macro_atom_data.transition_probability.values
        macro_atom.calculate_transition_probabilities(
            tpos,
            beta_sobolev.values,
            j_blues.values,
            stimulated_emission_factor,
            transition_type,
            lines_idx,
            self.block_references,
            transition_probabilities,
        )
        return transition_probabilities

    def calculate_transition_probabilities(
        self, macro_atom_data, beta_sobolev, j_blues, stimulated_emission_factor
    ):
        transition_probabilities = self.prepare_transition_probabilities(
            macro_atom_data, beta_sobolev, j_blues, stimulated_emission_factor
        )
        return transition_probabilities

    def initialize_macro_atom_transition_type_filters(
        self, atomic_data, macro_atom_data
    ):
        self.transition_up_filter = macro_atom_data.transition_type.values == 1
        self.transition_up_line_filter = macro_atom_data.lines_idx.values[
            self.transition_up_filter
        ]
        self.block_references = np.hstack(
            (
                atomic_data.macro_atom_references.block_references,
                len(macro_atom_data),
            )
        )

    @staticmethod
    def _get_transition_probability_coefs(macro_atom_data):
        return macro_atom_data.transition_probability.values[np.newaxis].T

    def prepare_transition_probabilities(
        self, macro_atom_data, beta_sobolev, j_blues, stimulated_emission_factor
    ):
        current_beta_sobolev = beta_sobolev.values.take(
            macro_atom_data.lines_idx.values, axis=0, mode="raise"
        )
        transition_probabilities = (
            self.transition_probability_coef * current_beta_sobolev
        )
        j_blues = j_blues.take(
            self.transition_up_line_filter, axis=0, mode="raise"
        )
        macro_stimulated_emission = stimulated_emission_factor.take(
            self.transition_up_line_filter, axis=0, mode="raise"
        )
        transition_probabilities[self.transition_up_filter] *= (
            j_blues * macro_stimulated_emission
        )
        return transition_probabilities

    def _normalize_transition_probabilities(self, transition_probabilities):
        macro_atom.normalize_transition_probabilities(
            transition_probabilities, self.block_references
        )

    def _new_normalize_transition_probabilities(self, transition_probabilites):
        for i, start_id in enumerate(self.block_references[:-1]):
            end_id = self.block_references[i + 1]
            block = transition_probabilites[start_id:end_id]
            transition_probabilites[start_id:end_id] *= 1 / ne.evaluate(
                "sum(block, 0)"
            )

    @staticmethod
    def _get_macro_atom_data(atomic_data):
        try:
            return atomic_data.macro_atom_data
        except:
            return atomic_data.macro_atom_data_all
