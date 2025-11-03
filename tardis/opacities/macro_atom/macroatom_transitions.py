from tardis import constants as const
import numpy as np
import pandas as pd
from tardis.plasma.properties.continuum_processes.rates import (
    get_ground_state_multi_index,
)
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType

CONST_C_CGS: float = const.c.cgs.value
CONST_H_CGS: float = const.h.cgs.value


def line_transition_internal_up(
    line_f_lus: np.ndarray,
    line_nus: np.ndarray,
    energies_lower: np.ndarray,
    mean_intensities_blue_wing: pd.DataFrame,
    beta_sobolevs: pd.DataFrame,
    stimulated_emission_factors: np.ndarray,
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
    line_f_lus : np.ndarray
        Oscillator strengths for line transitions from lower to upper levels.
    line_nus : np.ndarray
        Frequencies of the line transitions in Hz.
    energies_lower : np.ndarray
        Energy values of the lower levels in the transitions.
    mean_intensities_blue_wing : pd.DataFrame
        Mean radiation field intensities at the blue wing of the lines.
    beta_sobolevs : pd.DataFrame
        Sobolev escape probabilities for the line transitions.
    stimulated_emission_factors : pd.DataFrame
        Factors accounting for stimulated emission in the transitions.
    transition_a_i_l_u_array : np.ndarray
        2D array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : np.ndarray
        Unique identifiers for each line transition.

    Returns
    -------
    p_internal_up : pd.DataFrame
        DataFrame containing the calculated internal upward transition probabilities
        with source information, indexed by transition.
    internal_up_metadata : pd.DataFrame
        DataFrame containing metadata for the transitions including transition line IDs,
        source and destination level information, transition type, and line indices.

    Notes
    -----
    The transition probability is calculated using the formula:
    P = (f_lu / (h * nu)) * stimulated_emission_factor * mean_intensity * beta * E_lower

    The function uses MacroAtomTransitionType to mark the transition type.
    """
    p_internal_up = probability_internal_up(
        beta_sobolevs,
        line_nus,
        line_f_lus,
        stimulated_emission_factors,
        mean_intensities_blue_wing,
        energies_lower,
    )

    sources = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 2]]))
    transition_indices = np.arange(len(line_ids))
    destinations = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 3]]))

    internal_up_metadata = pd.DataFrame(
        {
            "transition_line_id": line_ids,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.INTERNAL_UP,
            "transition_line_idx": transition_indices,
        },
        index=p_internal_up.index,
    )
    p_internal_up["source"] = sources

    return p_internal_up, internal_up_metadata


def probability_internal_up(
    beta_sobolevs: pd.DataFrame,
    line_nus: np.ndarray,
    line_f_lus: np.ndarray,
    stimulated_emission_factors: np.ndarray,
    mean_intensities_blue_wing: pd.DataFrame,
    energies_lower: np.ndarray,
) -> pd.DataFrame:
    """
    Calculate internal upward transition probabilities.

    This function computes the probability of internal upward transitions between
    energy levels in macro atoms. The calculation involves oscillator strengths,
    stimulated emission factors, mean radiation field intensities, and energy values.

    Parameters
    ----------
    beta_sobolevs : pd.DataFrame
        Sobolev escape probabilities for the line transitions.
    line_nus : np.ndarray
        Frequencies of the line transitions in Hz.
    line_f_lus : np.ndarray
        Oscillator strengths for line transitions from lower to upper levels.
    stimulated_emission_factors : np.ndarray
        Factors accounting for stimulated emission in the transitions.
    mean_intensities_blue_wing : pd.DataFrame
        Mean radiation field intensities at the blue wing of the lines.
    energies_lower : np.ndarray
        Energy values of the lower levels in the transitions.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated internal upward transition probabilities.
    """
    p_internal_up = beta_sobolevs * (
        line_f_lus
        / (CONST_H_CGS * line_nus)
        * stimulated_emission_factors
        * mean_intensities_blue_wing
        * energies_lower
    )
    return p_internal_up


def line_transition_internal_down(
    line_f_uls: np.ndarray,
    line_nus: np.ndarray,
    energies_lower: np.ndarray,
    beta_sobolevs: pd.DataFrame,
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
    line_f_uls : np.ndarray
        Oscillator strengths for line transitions.
    line_nus : np.ndarray
        Frequencies of line transitions.
    energies_lower : np.ndarray
        Energies of the lower levels involved in transitions.
    beta_sobolevs : pd.DataFrame
        Sobolev beta factors for the transitions.
    transition_a_i_l_u_array : np.ndarray
        Array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : np.ndarray
        Identifiers for each line transition.

    Returns
    -------
    p_internal_down : pd.DataFrame
        DataFrame containing the calculated internal downward transition
        probabilities with source information.
    internal_down_metadata : pd.DataFrame
        DataFrame containing metadata for the transitions including
        transition line IDs, source and destination level information,
        transition type, and transition line indices.

    Notes
    -----
    The internal downward transition probability is calculated using the formula:
    p = 2 * nu^2 * f_ul / c^2 * beta * E_lower

    The function uses MacroAtomTransitionType to mark the transition type.
    """
    p_internal_down = probability_internal_down(
        beta_sobolevs, line_nus, line_f_uls, energies_lower
    )

    sources = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 3]]))
    transition_indices = np.arange(len(line_ids))
    destinations = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 2]]))

    internal_down_metadata = pd.DataFrame(
        {
            "transition_line_id": line_ids,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.INTERNAL_DOWN,
            "transition_line_idx": transition_indices,
        },
        index=p_internal_down.index,
    )
    p_internal_down["source"] = sources

    return p_internal_down, internal_down_metadata


def probability_internal_down(
    beta_sobolevs: pd.DataFrame,
    line_nus: np.ndarray,
    line_f_uls: np.ndarray,
    energies_lower: np.ndarray,
) -> pd.DataFrame:
    """
    Calculate internal downward transition probabilities.

    This function computes the probability of internal downward transitions between
    energy levels in macro atoms. The calculation is based on oscillator strengths,
    line frequencies, and lower level energies.

    Parameters
    ----------
    beta_sobolevs : pd.DataFrame
        Sobolev escape probabilities for the line transitions.
    line_nus : np.ndarray
        Frequencies of the line transitions in Hz.
    line_f_uls : np.ndarray
        Oscillator strengths for line transitions from upper to lower levels.
    energies_lower : np.ndarray
        Energy values of the lower levels in the transitions.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated internal downward transition probabilities.
    """
    p_internal_down = beta_sobolevs * (
        2 * line_nus**2 * line_f_uls / CONST_C_CGS**2 * energies_lower
    )
    return p_internal_down


def line_transition_emission_down(
    line_f_uls: np.ndarray,
    line_nus: np.ndarray,
    energies_upper: np.ndarray,
    energies_lower: np.ndarray,
    beta_sobolevs: pd.DataFrame,
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
    line_f_uls : np.ndarray
        Oscillator strengths for the line transitions.
    line_nus : np.ndarray
        Frequencies of the line transitions in Hz.
    energies_upper : np.ndarray
        Energy values of the upper levels in the transitions.
    energies_lower : np.ndarray
        Energy values of the lower levels in the transitions.
    beta_sobolevs : pd.DataFrame
        Sobolev escape probabilities for the transitions.
    transition_a_i_l_u_array : np.ndarray
        Array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : np.ndarray
        Unique identifiers for each line transition.

    Returns
    -------
    p_emission_down : pd.DataFrame
        DataFrame containing the calculated emission down probabilities with
        source information. Contains columns for the probability values and
        'source' tuples.
    emission_down_metadata : pd.DataFrame
        DataFrame containing metadata for the transitions including transition
        line IDs, source and destination level information, transition type,
        and transition line indices.

    Notes
    -----
    The emission down probability is calculated using the formula:
    P = 2 * nu^2 * f_ul / c^2 * beta * (E_upper - E_lower)

    The function uses MacroAtomTransitionType to mark the transition type.
    """
    p_emission_down = probability_emission_down(
        beta_sobolevs, line_nus, line_f_uls, energies_upper, energies_lower
    )

    sources = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 3]]))
    transition_indices = np.arange(len(line_ids))
    destinations = list(map(tuple, transition_a_i_l_u_array[:, [0, 1, 2]]))

    emission_down_metadata = pd.DataFrame(
        {
            "transition_line_id": line_ids,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.BB_EMISSION,
            "transition_line_idx": transition_indices,
        },
        index=p_emission_down.index,
    )
    p_emission_down["source"] = sources

    return p_emission_down, emission_down_metadata


def probability_emission_down(
    beta_sobolevs: pd.DataFrame,
    line_nus: np.ndarray,
    line_f_uls: np.ndarray,
    energies_upper: np.ndarray,
    energies_lower: np.ndarray,
) -> pd.DataFrame:
    """
    Calculate emission down transition probabilities.

    This function computes the probability of emission down transitions for
    atomic line transitions. The calculation considers oscillator strengths,
    line frequencies, and the energy difference between upper and lower levels.

    Parameters
    ----------
    beta_sobolevs : pd.DataFrame
        Sobolev escape probabilities for the line transitions.
    line_nus : np.ndarray
        Frequencies of the line transitions in Hz.
    line_f_uls : np.ndarray
        Oscillator strengths for line transitions from upper to lower levels.
    energies_upper : np.ndarray
        Energy values of the upper levels in the transitions.
    energies_lower : np.ndarray
        Energy values of the lower levels in the transitions.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated emission down transition probabilities.
    """
    p_emission_down = beta_sobolevs * (
        2
        * line_nus**2
        * line_f_uls
        / CONST_C_CGS**2
        * (energies_upper - energies_lower)
    )
    return p_emission_down


### Below are continuum transition specific functions


def continuum_transition_recombination_internal(
    spontaneous_recombination_coeff: pd.DataFrame,  # These will all be changes to spontaneous recombination coefficient
    photoionization_data_level_energies: pd.Series,  # This is just energy of the excitation state in the atomic data
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate unnormalized probabilities of radiative recombination.

    Parameters
    ----------
    spontaneous_recombination_coeff : pd.Series
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    photoionization_data_level_energies : pd.Series
        Energies of levels with bound-free transitions.

    Returns
    -------
    p_recombination : pd.DataFrame
        DataFrame containing recombination probabilities.
    recombination_metadata : pd.DataFrame
        DataFrame containing metadata for the recombination transitions.
    """
    p_recomb_internal = probability_recombination_internal(
        spontaneous_recombination_coeff, photoionization_data_level_energies
    )

    destinations = p_recomb_internal.index.values
    sources = get_ground_state_multi_index(p_recomb_internal.index).values

    recombination_internal_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.RECOMB_INTERNAL,
            "transition_line_idx": -99,
            "continuum_transition_idx": range(
                len(photoionization_data_level_energies)
            ),
        },
        index=p_recomb_internal.index,
    )

    p_recomb_internal["source"] = sources

    return p_recomb_internal, recombination_internal_metadata


def probability_recombination_internal(
    spontaneous_recombination_coeff: pd.DataFrame,
    photoionization_data_level_energies: pd.Series,
) -> pd.DataFrame:
    """
    Calculate unnormalized probabilities of radiative recombination (internal transitions).

    Parameters
    ----------
    spontaneous_recombination_coeff : pd.DataFrame
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    photoionization_data_level_energies : pd.Series
        Energies of levels with bound-free transitions.

    Returns
    -------
    pd.DataFrame
        DataFrame containing unnormalized recombination probabilities for internal transitions.
    """
    p_recomb_internal = spontaneous_recombination_coeff.multiply(
        photoionization_data_level_energies, axis=0
    )
    return p_recomb_internal


def continuum_transition_recombination_emission(
    spontaneous_recombination_coeff: pd.DataFrame,
    photoionization_data_frequencies: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate unnormalized probabilities and metadata for radiative recombination emission transitions (continuum).

    Parameters
    ----------
    spontaneous_recombination_coeff : pd.DataFrame
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    photoionization_data_frequencies : pd.Series
        Ionization threshold frequencies for the levels.

    Returns
    -------
    p_recomb_emission : pd.DataFrame
        DataFrame containing unnormalized recombination emission probabilities.
    recombination_emission_metadata : pd.DataFrame
        DataFrame containing metadata for the recombination emission transitions.
    """
    p_recomb_emission = probability_recombination_emission(
        spontaneous_recombination_coeff, photoionization_data_frequencies
    )

    destinations = p_recomb_emission.index.values
    sources = get_ground_state_multi_index(p_recomb_emission.index).values

    recombination_emission_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.RECOMB_EMISSION,
            "transition_line_idx": -99,
            "continuum_transition_idx": range(
                len(photoionization_data_frequencies)
            ),
        },
        index=p_recomb_emission.index,
    )

    p_recomb_emission["source"] = sources

    return p_recomb_emission, recombination_emission_metadata


def probability_recombination_emission(
    spontaneous_recombination_coeff: pd.DataFrame,
    photoionization_data_frequencies: pd.Series,
) -> pd.DataFrame:
    """
    Calculate unnormalized probabilities of radiative recombination emission transitions.

    Parameters
    ----------
    spontaneous_recombination_coeff : pd.DataFrame
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    photoionization_data_frequencies : pd.Series
        Ionization threshold frequencies for the levels.

    Returns
    -------
    pd.DataFrame
        DataFrame containing unnormalized recombination emission probabilities.
    """
    p_recomb_emission = (
        spontaneous_recombination_coeff.multiply(
            photoionization_data_frequencies, axis=0
        )
        * CONST_H_CGS
    )  # photoionization_data_frequencies is ionization threshold, so I guess delta e is just photoionization_data_frequencies * h?
    return p_recomb_emission


def continuum_transition_photoionization(
    stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
    photoionization_data_level_energies: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate photoionization probability unnormalized.

    Parameters
    ----------
    stim_recomb_corrected_photoionization_rate_coeff : pd.Series
        Corrected photoionization rate coefficient from level `i` to `k`.
    photoionization_data_level_energies : pd.Series
        Energies of the levels involved in photoionization.

    Returns
    -------
    p_photoionization : pd.DataFrame
        DataFrame containing photoionization probabilities.
    photoionization_metadata : pd.DataFrame
        DataFrame containing metadata for the photoionization transitions.
    """
    p_photoionization = probability_photoionization(
        stim_recomb_corrected_photoionization_rate_coeff,
        photoionization_data_level_energies,
    )

    sources = p_photoionization.index.values
    destinations = get_ground_state_multi_index(p_photoionization.index).values

    photoionization_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.PHOTOIONIZATION,
            "transition_line_idx": -99,
            "continuum_transition_idx": range(
                len(photoionization_data_level_energies)
            ),
        },
        index=p_photoionization.index,
    )

    p_photoionization["source"] = sources

    return p_photoionization, photoionization_metadata


def probability_photoionization(
    stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
    photoionization_data_level_energies: pd.Series,
) -> pd.DataFrame:
    """
    Calculate photoionization probability unnormalized.

    Parameters
    ----------
    stim_recomb_corrected_photoionization_rate_coeff : pd.Series
        Corrected photoionization rate coefficient from level `i` to `k`.
    photoionization_data_level_energies : pd.Series
        Energies of the levels involved in photoionization.

    Returns
    -------
    pd.DataFrame
        DataFrame containing photoionization probabilities.
    """
    p_photoionization = (
        stim_recomb_corrected_photoionization_rate_coeff.multiply(
            photoionization_data_level_energies, axis=0
        )
    )
    return p_photoionization


# def continuum_transition_collisional(
#     self,
#     coll_exc_coeff: pd.DataFrame,
#     coll_deexc_coeff: pd.DataFrame,
#     yg_idx: pd.DataFrame,
#     electron_densities: pd.DataFrame,
#     delta_E_yg: pd.Series,
#     atomic_data,
#     level_number_density: pd.Series,
# ) -> pd.DataFrame:
#     """
#     Calculate collisional transition probabilities.

#     Parameters
#     ----------
#     self
#         The instance of the class containing this method.
#     coll_exc_coeff : pd.DataFrame
#         Collisional excitation coefficients.
#     coll_deexc_coeff : pd.DataFrame
#         Collisional de-excitation coefficients.
#     yg_idx : pd.DataFrame
#         Index DataFrame containing level information.
#     electron_densities : pd.DataFrame
#         Electron densities.
#     delta_E_yg : pd.Series
#         Energy differences.
#     atomic_data
#         Atomic data containing level information.
#     level_number_density : pd.Series
#         Number densities of levels.

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame containing collisional transition probabilities.
#     """
#     p_deexc_deactivation = (coll_deexc_coeff * electron_densities).multiply(
#         delta_E_yg.values, axis=0
#     )
#     p_deexc_deactivation = self.set_index(
#         p_deexc_deactivation, yg_idx, reverse=True
#     )
#     p_deexc_deactivation = p_deexc_deactivation.groupby(level=[0]).sum()
#     index_dd = pd.MultiIndex.from_product(
#         [p_deexc_deactivation.index.values, ["k"], [0]],
#         names=list(yg_idx.columns) + ["transition_type"],
#     )
#     p_deexc_deactivation = p_deexc_deactivation.set_index(index_dd)

#     level_lower_index = coll_deexc_coeff.index.droplevel("level_number_upper")
#     energy_lower = atomic_data.levels.energy.loc[level_lower_index]
#     p_deexc_internal = (coll_deexc_coeff * electron_densities).multiply(
#         energy_lower.values, axis=0
#     )
#     p_deexc_internal = self.set_index(
#         p_deexc_internal, yg_idx, transition_type=0, reverse=True
#     )

#     p_exc_internal = (coll_exc_coeff * electron_densities).multiply(
#         energy_lower.values, axis=0
#     )
#     p_exc_internal = self.set_index(
#         p_exc_internal, yg_idx, transition_type=0, reverse=False
#     )
#     p_exc_cool = (coll_exc_coeff * electron_densities).multiply(
#         delta_E_yg.values, axis=0
#     )
#     p_exc_cool = p_exc_cool * level_number_density.loc[level_lower_index].values
#     p_exc_cool = self.set_index(p_exc_cool, yg_idx, reverse=False)
#     p_exc_cool = p_exc_cool.groupby(level="destination_level_idx").sum()
#     exc_cool_index = pd.MultiIndex.from_product(
#         [["k"], p_exc_cool.index.values, [0]],
#         names=list(yg_idx.columns) + ["transition_type"],
#     )
#     p_exc_cool = p_exc_cool.set_index(exc_cool_index)
#     p_coll = pd.concat(
#         [p_deexc_deactivation, p_deexc_internal, p_exc_internal, p_exc_cool]
#     )
#     return p_coll

#     # def calculate(self, electron_densities, t_electrons, time_explosion):
#     #     cool_rate_adiabatic = (
#     #         3.0 * electron_densities * K_B * t_electrons
#     #     ) / time_explosion

#     #     cool_rate_adiabatic = cooling_rate_series2dataframe(
#     #         cool_rate_adiabatic, destination_level_idx="adiabatic"
#     #     )
#     #     return cool_rate_adiabatic


# def collisional_transition_deexcite_deactivate(
#     coll_deexc_coeff, yg_idx, electron_densities, delta_E_yg
# ):
#     p_deexc_deactivation = probability_collision_deexc_deactivate(
#         coll_deexc_coeff, electron_densities, delta_E_yg
#     )
#     # yg idx is just source_level_idx, destination_level_idx
#     # p_deexc_deactivation = p_deexc_deactivation.groupby(level=[0]).sum()
#     # index_dd = pd.MultiIndex.from_product(
#     #     [p_deexc_deactivation.index.values, ["k"], [0]],
#     #     names=list(yg_idx.columns) + ["transition_type"],
#     # )
#     # p_deexc_deactivation = p_deexc_deactivation.set_index(index_dd)

#     coll_deexc_deactivate_metadata = pd.DataFrame(
#         {
#             "transition_line_id": -99,
#             "source": p_deexc_deactivation.index.values,
#             "destination": get_ground_state_multi_index(
#                 p_deexc_deactivation.index
#             ).values,
#             "transition_type": MacroAtomTransitionType.COLL_DEEXCITATION_DEACTIVATION,
#             "transition_line_idx": -99,
#             "continuum_transition_idx": range(len(delta_E_yg)),
#         },
#         index=p_deexc_deactivation.index,
#     )

#     return p_deexc_deactivation, coll_deexc_deactivate_metadata


# def probability_collision_deexc_deactivate(
#     coll_deexc_coeff, electron_densities, delta_E_yg
# ) -> pd.DataFrame:
#     """
#     Calculate collisional de-excitation and deactivation probabilities.

#     Parameters
#     ----------
#     coll_deexc_coeff : pd.DataFrame
#         Collisional de-excitation coefficients.
#     electron_densities : pd.DataFrame
#         Electron densities.
#     delta_E_yg : pd.Series
#         Energy differences.

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame containing collisional de-excitation and deactivation probabilities.
#     """
#     p_deexc_deactivation = (coll_deexc_coeff * electron_densities).multiply(
#         delta_E_yg.values, axis=0
#     )
#     p_deexc_deactivation = p_deexc_deactivation.groupby(level=[0]).sum()

#     return p_deexc_deactivation
