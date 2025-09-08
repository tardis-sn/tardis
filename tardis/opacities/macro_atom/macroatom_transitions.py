from __future__ import annotations

from astropy import constants as const
import numpy as np
import pandas as pd

CONST_C_CGS: float = const.c.cgs.value
CONST_H_CGS: float = const.h.cgs.value

P_INTERNAL_UP: int = 1
P_INTERNAL_DOWN: int = 0
P_EMISSION_DOWN: int = -1


def line_transition_internal_up(
    line_f_lus: pd.Series,
    line_nus: pd.Series,
    energies_lower: pd.Series,
    mean_intensities_blue_wing: pd.Series,
    beta_sobolevs: pd.Series,
    stimulated_emission_factors: pd.Series,
    transition_a_i_l_u_array: np.ndarray,
    line_ids: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate internal upward transition probabilities for line transitions in macro atoms.

    This function computes the probability of internal upward transitions between
    energy levels in macro atoms due to line transitions. It calculates the
    transition probabilities and creates metadata for tracking the transitions.

    Parameters
    ----------
    line_f_lus : array_like
        Oscillator strengths for line transitions from lower to upper levels.
    line_nus : array_like
        Frequencies of the line transitions in Hz.
    energies_lower : array_like
        Energy values of the lower levels in the transitions.
    mean_intensities_blue_wing : array_like
        Mean radiation field intensities at the blue wing of the lines.
    beta_sobolevs : array_like
        Sobolev escape probabilities for the line transitions.
    stimulated_emission_factors : array_like
        Factors accounting for stimulated emission in the transitions.
    transition_a_i_l_u_array : array_like
        2D array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : array_like
        Unique identifiers for each line transition.

    Returns
    -------
    p_internal_up : pandas.DataFrame
        DataFrame containing the calculated internal upward transition probabilities
        with source information, indexed by transition.
    internal_up_metadata : pandas.DataFrame
        DataFrame containing metadata for the transitions including transition line IDs,
        source and destination level information, transition type, and line indices.

    Notes
    -----
    The transition probability is calculated using the formula:
    P = (f_lu / (h * nu)) * stimulated_emission_factor * mean_intensity * beta * E_lower

    The function uses the constant P_INTERNAL_UP to mark the transition type.
    """

    p_internal_up = (
        line_f_lus
        / (CONST_H_CGS * line_nus)
        * stimulated_emission_factors
        * mean_intensities_blue_wing
        * beta_sobolevs
        * energies_lower
    )

    transition_indices = np.arange(len(line_ids))
    sources = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 2]]))
    destinations = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 3]]))

    internal_up_metadata = pd.DataFrame(
        {
            "transition_line_id": line_ids,
            "source": sources,
            "destination": destinations,
            "transition_type": P_INTERNAL_UP,
            "transition_line_idx": transition_indices,
        },
        index=p_internal_up.index,
    )
    p_internal_up["source"] = sources

    return p_internal_up, internal_up_metadata


def line_transition_internal_down(
    line_f_uls: pd.Series,
    line_nus: pd.Series,
    energies_lower: pd.Series,
    beta_sobolevs: pd.Series,
    transition_a_i_l_u_array: np.ndarray,
    line_ids: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate internal downward transition probabilities for line transitions.

    This function computes the probability of internal downward transitions
    in macro atoms for line transitions, based on oscillator strengths,
    frequencies, and other atomic parameters.

    Parameters
    ----------
    line_f_uls : array_like
        Oscillator strengths for line transitions.
    line_nus : array_like
        Frequencies of line transitions.
    energies_lower : array_like
        Energies of the lower levels involved in transitions.
    beta_sobolevs : array_like
        Sobolev beta factors for the transitions.
    transition_a_i_l_u_array : array_like
        Array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : array_like
        Identifiers for each line transition.

    Returns
    -------
    p_internal_down : pandas.DataFrame
        DataFrame containing the calculated internal downward transition
        probabilities with source information.
    internal_down_metadata : pandas.DataFrame
        DataFrame containing metadata for the transitions including
        transition line IDs, source and destination level information,
        transition type, and transition line indices.

    Notes
    -----
    The internal downward transition probability is calculated using the formula:
    p = 2 * nu^2 * f_ul / c^2 * beta * E_lower

    The function uses the constant P_INTERNAL_DOWN to mark the transition type.
    """
    p_internal_down = (
        2
        * line_nus**2
        * line_f_uls
        / CONST_C_CGS**2
        * beta_sobolevs
        * energies_lower
    )

    transition_indices = np.arange(len(line_ids))
    sources = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 3]]))
    destinations = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 2]]))

    internal_down_metadata = pd.DataFrame(
        {
            "transition_line_id": line_ids,
            "source": sources,
            "destination": destinations,
            "transition_type": P_INTERNAL_DOWN,
            "transition_line_idx": transition_indices,
        },
        index=p_internal_down.index,
    )
    p_internal_down["source"] = sources

    return p_internal_down, internal_down_metadata


def line_transition_emission_down(
    line_f_uls: pd.Series,
    line_nus: pd.Series,
    energies_upper: pd.Series,
    energies_lower: pd.Series,
    beta_sobolevs: pd.Series,
    transition_a_i_l_u_array: np.ndarray,
    line_ids: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate emission down transition probabilities for line transitions.

    This function computes the probability of emission down transitions for
    atomic line transitions using oscillator strengths, frequencies, energy
    differences, and Sobolev optical depths.

    Parameters
    ----------
    line_f_uls : array_like
        Oscillator strengths for the line transitions.
    line_nus : array_like
        Frequencies of the line transitions in Hz.
    energies_upper : array_like
        Energy values of the upper levels in the transitions.
    energies_lower : array_like
        Energy values of the lower levels in the transitions.
    beta_sobolevs : array_like
        Sobolev escape probabilities for the transitions.
    transition_a_i_l_u_array : array_like
        Array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : array_like
        Unique identifiers for each line transition.

    Returns
    -------
    p_emission_down : pandas.DataFrame
        DataFrame containing the calculated emission down probabilities with
        source information. Contains columns for the probability values and
        'source' tuples.
    emission_down_metadata : pandas.DataFrame
        DataFrame containing metadata for the transitions including transition
        line IDs, source and destination level information, transition type,
        and transition line indices.

    Notes
    -----
    The emission down probability is calculated using the formula:
    P = 2 * nu^2 * f_ul / c^2 * beta * (E_upper - E_lower)

    The function uses the constant P_EMISSION_DOWN to mark the transition type.
    """
    p_emission_down = (
        2
        * line_nus**2
        * line_f_uls
        / CONST_C_CGS**2
        * beta_sobolevs
        * (energies_upper - energies_lower)
    )

    transition_indices = np.arange(len(line_ids))
    sources = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 3]]))
    destinations = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 2]]))

    emission_down_metadata = pd.DataFrame(
        {
            "transition_line_id": line_ids,
            "source": sources,
            "destination": destinations,
            "transition_type": P_EMISSION_DOWN,
            "transition_line_idx": transition_indices,
        },
        index=p_emission_down.index,
    )
    p_emission_down["source"] = sources

    return p_emission_down, emission_down_metadata
