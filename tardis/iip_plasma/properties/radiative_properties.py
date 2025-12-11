import logging

import numexpr as ne
import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u

from tardis.iip_plasma.properties.base import (
    ConvergedPlasmaProperty,
    ProcessingPlasmaProperty,
)
from tardis.opacities.tau_sobolev import calculate_beta_sobolev
from tardis.opacities.macro_atom.util import (
    fast_calculate_transition_probabilities,
)
from tardis.opacities.macro_atom.macroatom_solver import (
    BoundBoundMacroAtomSolver,
)

logger = logging.getLogger(__name__)

__all__ = [
    "BetaSobolev",
    "LTEJBlues",
    "StimulatedEmissionFactor",
    "TauSobolev",
    "TauSobolevDeriv",
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
    latex_formula = ("1-\\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}",)

    def __init__(self, plasma_parent=None, nlte_species=None):
        super(StimulatedEmissionFactor, self).__init__(plasma_parent)
        self._g_upper = None
        self._g_lower = None
        try:
            self.nlte_species = self.plasma_parent.nlte_species
        except:
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
            "1 - ((g_lower * n_upper) / (g_upper * n_lower))"
        )
        stimulated_emission_factor[n_lower == 0.0] = 0.0
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
                (stimulated_emission_factor < 0)
                & nlte_lines_mask[np.newaxis, :].T
            ] = 0.0
        return stimulated_emission_factor


class TauSobolev(ConvergedPlasmaProperty):
    """
    Attributes
    ----------
    tau_sobolev : Pandas DataFrame, dtype float
                  Sobolev optical depth for each line. Indexed by line.
                  Columns as zones.
    """

    outputs = ("tau_sobolevs",)
    latex_name = ("\\tau_{\\textrm{sobolev}}",)
    latex_formula = (
        "\\dfrac{\\pi e^{2}}{m_{e} c}f_{lu}\\lambda t_{exp}\
        n_{lower} \\Big(1-\\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}\\Big)",
    )

    def __init__(self, plasma_parent):
        super(TauSobolev, self).__init__(plasma_parent)
        self.sobolev_coefficient = (
            (
                ((np.pi * const.e.gauss**2) / (const.m_e.cgs * const.c.cgs))
                * u.cm
                * u.s
                / u.cm**3
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
            columns=level_number_density.columns,
        )


class BetaSobolev(ConvergedPlasmaProperty):
    """
    Attributes
    ----------
    beta_sobolev : Numpy Array, dtype float
    """

    outputs = ("beta_sobolev",)
    latex_name = ("\\beta_{\\textrm{sobolev}}",)

    def calculate(self, tau_sobolevs):
        if getattr(self, "beta_sobolev", None) is None:
            beta_sobolev = np.zeros_like(tau_sobolevs.values)
        else:
            beta_sobolev = self.beta_sobolev

        beta_sobolev = calculate_beta_sobolev(tau_sobolevs)

        # TODO: REvert
        # beta_sobolev = np.ones_like(beta_sobolev) * 0.05
        # self.plasma_parent.previous_beta_sobolev.update(beta_sobolev)
        # self.plasma_parent.plasma_properties_dict['PreviousBetaSobolev'].set_value(beta_sobolev)

        # conv = np.fabs(beta_sobolev - self.plasma_parent.previous_beta_sobolev)/beta_sobolev
        # print "Beta_sobolev conv:", conv.max(axis=0)
        # if (conv > 1e-2).sum().sum():
        #    self.plasma_parent.update(previous_beta_sobolev=beta_sobolev)
        return beta_sobolev


class TransitionProbabilities(ConvergedPlasmaProperty):
    """
    Attributes
    ----------
    transition_probabilities : Pandas DataFrame, dtype float
    """

    outputs = ("transition_probabilities",)

    def __init__(self, plasma_parent):
        super(TransitionProbabilities, self).__init__(plasma_parent)
        self.initialize = True
        try:
            self.continuum_treatment = plasma_parent.continuum_treatment
        except:
            self.continuum_treatment = False

    def calculate(
        self,
        atomic_data,
        beta_sobolev,
        j_blues,
        stimulated_emission_factor,
        tau_sobolevs,
    ):
        # print "Calculating Transition probabilities"
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
        # optimized_calculate_transition_probabilities(tpos, beta_sobolev, j_blues, stimulated_emission_factor, transition_type, lines_idx, self.block_references, transition_probabilities)
        fast_calculate_transition_probabilities(
            tpos,
            beta_sobolev.values,
            j_blues,
            stimulated_emission_factor,
            transition_type,
            lines_idx,
            self.block_references,
            transition_probabilities,
            self.continuum_treatment,
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
        current_beta_sobolev = beta_sobolev.take(
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
        self.block_references = (
            BoundBoundMacroAtomSolver.normalize_transition_probabilities(
                transition_probabilities
            )
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


class LTEJBlues(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    lte_j_blues : Pandas DataFrame, dtype float
                  J_blue values as calculated in LTE.
    """

    outputs = ("lte_j_blues",)
    latex_name = "J^{b}_{lu(LTE)}"

    @staticmethod
    def calculate(lines, nu, beta_rad):
        beta_rad = pd.Series(beta_rad)
        nu = pd.Series(nu)
        h = const.h.cgs.value
        c = const.c.cgs.value
        df = pd.DataFrame(1, index=nu.index, columns=beta_rad.index)
        df = df.mul(nu, axis="index") * beta_rad
        exponential = (np.exp(h * df) - 1) ** (-1)
        remainder = 2 * (h * nu.values**3) / (c**2)
        j_blues = exponential.mul(remainder, axis=0)
        return pd.DataFrame(j_blues, index=lines.index, columns=beta_rad.index)


class TauSobolevDeriv(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    tau_sobolev : Pandas DataFrame, dtype float
                  Sobolev optical depth for each line. Indexed by line.
                  Columns as zones.
    """

    outputs = ("dtau_dnl", "dtau_dnu")
    # TODO
    latex_name = ("\\tau_{\\textrm{sobolev}}",)
    latex_formula = (
        "\\dfrac{\\pi e^{2}}{m_{e} c}f_{lu}\\lambda t_{exp}\
        n_{lower} \\Big(1-\\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}\\Big)",
    )

    def __init__(self, plasma_parent):
        super(TauSobolevDeriv, self).__init__(plasma_parent)
        self.sobolev_coefficient = (
            (
                ((np.pi * const.e.gauss**2) / (const.m_e.cgs * const.c.cgs))
                * u.cm
                * u.s
                / u.cm**3
            )
            .to(1)
            .value
        )

    def calculate(
        self,
        lines,
        lines_lower_level_index,
        lines_upper_level_index,
        time_explosion,
        f_lu,
        wavelength_cm,
        g,
    ):
        g_lower = np.array(g.iloc[lines_lower_level_index], dtype=np.float64)
        g_upper = np.array(g.iloc[lines_upper_level_index], dtype=np.float64)

        f_lu = f_lu.values[np.newaxis].T
        wavelength = wavelength_cm.values[np.newaxis].T
        dtau_dnl = self.sobolev_coefficient * f_lu * wavelength * time_explosion

        dtau_dnl = pd.DataFrame(dtau_dnl, index=lines.index)

        dtau_dnu = 0
        return dtau_dnl, dtau_dnu
