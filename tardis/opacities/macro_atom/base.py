import logging

import numpy as np
import pandas as pd

from tardis.opacities.macro_atom import util
from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)


def initialize_transition_probabilities(
    atomic_data,
):
    """Convienience Function for initializing the transition probabilities

    Parameters
    ----------
    atomic_data : tardis.io.atom_data.AtomData
            Atomic Data

    Returns
    -------
    dict
        "transition_probability_coef" : np.ndarray
            Reshaped macro atom transition probabilities
        "block_references": np.ndarray
            macro atom block references
    """
    macro_atom_data = get_macro_atom_data(atomic_data)
    (
        transition_up_filter,
        transition_up_line_filter,
        block_references,
    ) = initialize_macro_atom_transition_type_filters(
        atomic_data, macro_atom_data
    )
    transition_probability_coef = get_transition_probability_coefs(
        macro_atom_data
    )

    return {
        "transition_probability_coef": transition_probability_coef,
        "block_references": block_references,
    }


def initialize_macro_atom_transition_type_filters(atomic_data, macro_atom_data):
    """Get the filters and block references from the macro atom

    Parameters
    ----------
    atomic_data : tardis.io.atom_data.AtomData
            Atomic Data
    macro_atom_data : pd.DataFrame
        Macro Atom Data

    Returns
    -------
    np.ndarray
        Mask where the transition type is 1
    np.ndarray
        index of lines at these locations
    pd.ndarray
        macro atom block references
    """
    transition_up_filter = macro_atom_data.transition_type.values == 1
    transition_up_line_filter = macro_atom_data.lines_idx.values[
        transition_up_filter
    ]
    block_references = np.hstack(
        (
            atomic_data.macro_atom_references.block_references,
            len(macro_atom_data),
        )
    )

    return transition_up_filter, transition_up_line_filter, block_references


def get_transition_probability_coefs(macro_atom_data):
    """Coefficients of the transition probabilities

    Parameters
    ----------
    macro_atom_data : pd.DataFrame
        Macro Atom Data

    Returns
    -------
    np.ndarray
        Reshaped macro atom transition probabilities
    """
    return macro_atom_data.transition_probability.values[np.newaxis].T


def get_macro_atom_data(atomic_data):
    """Get the macro atom data from the atomic data

    Parameters
    ----------
    atomic_data : tardis.io.atom_data.AtomData
            Atomic Data

    Returns
    -------
    pd.DataFrame
        The macro atom data in the plasma
    """
    try:
        return atomic_data.macro_atom_data
    except:
        logger.debug(
            "Macro Atom Data was not found. Instead returning All Macro Atom Data"
        )
        return atomic_data.macro_atom_data_all


class TransitionProbabilities(
    ProcessingPlasmaProperty
):  # Base MacroAtom Property
    """
    Attributes
    ----------
    transition_probabilities : Pandas DataFrame, dtype float
    """

    outputs = ("transition_probabilities",)

    def __init__(self, plasma_parent):
        super(TransitionProbabilities, self).__init__(plasma_parent)
        self.initialize = True
        self.normalize = True

    def calculate(
        self,
        atomic_data,
        beta_sobolev,
        j_blues,
        stimulated_emission_factor,
        tau_sobolevs,
    ):
        # Referenced in https://github.com/tardis-sn/tardis/issues/3009
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
        transition_type = macro_atom_data.transition_type.values
        lines_idx = macro_atom_data.lines_idx.values
        tpos = macro_atom_data.transition_probability.values
        util.fast_calculate_transition_probabilities(
            tpos,
            beta_sobolev.values,
            j_blues.values,
            stimulated_emission_factor,
            transition_type,
            lines_idx,
            self.block_references,
            transition_probabilities,
            self.normalize,
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

    @staticmethod
    def _get_macro_atom_data(atomic_data):
        try:
            return atomic_data.macro_atom_data
        except:
            logger.debug(
                "Macro Atom Data was not found. Instead returning All Macro Atom Data"
            )
            return atomic_data.macro_atom_data_all


class NonMarkovChainTransitionProbabilities(
    TransitionProbabilities
):  # Continuum Only
    outputs = ("non_markov_transition_probabilities",)
