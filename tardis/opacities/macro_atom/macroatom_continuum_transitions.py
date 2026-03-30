import numpy as np
import pandas as pd

from tardis.plasma.properties.continuum_processes.rates import (
    F_K,
    K_B,
    get_ground_state_multi_index,
)
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType


def probability_recombination_internal(
    spontaneous_recombination_coeff: pd.DataFrame,
    photoionization_data_level_energies: pd.Series,
) -> pd.DataFrame:
    """
    Calculate unnormalized probabilities of radiative recombination (internal transitions).

    Parameters
    ----------
    spontaneous_recombination_coeff
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    photoionization_data_level_energies
        Energies of levels with bound-free transitions.

    Returns
    -------
    p_recomb_internal
        DataFrame containing unnormalized recombination probabilities for internal transitions.
    """
    p_recomb_internal = spontaneous_recombination_coeff.multiply(
        photoionization_data_level_energies, axis=0
    )
    return p_recomb_internal


def continuum_transition_recombination_internal(
    spontaneous_recombination_coeff: pd.DataFrame,  # These will all be changes to spontaneous recombination coefficient
    photoionization_data_level_energies: pd.Series,  # This is just energy of the excitation state in the atomic data
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate unnormalized probabilities of radiative recombination.

    Parameters
    ----------
    spontaneous_recombination_coeff
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    photoionization_data_level_energies
        Energies of levels with bound-free transitions.

    Returns
    -------
    p_recombination
        DataFrame containing recombination probabilities.
    recombination_metadata
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
            "photoionization_key_idx": range(
                len(photoionization_data_level_energies)
            ),
            "collision_key_idx": -99,
        },
        index=p_recomb_internal.index,
    )

    return p_recomb_internal, recombination_internal_metadata


def probability_recombination_emission(
    spontaneous_recombination_coeff: pd.DataFrame,
    energies_diff_bound_free: pd.Series,
) -> pd.DataFrame:
    """
    Calculate unnormalized probabilities of radiative recombination emission transitions.

    Parameters
    ----------
    spontaneous_recombination_coeff
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    energies_diff_bound_free
        Ionization threshold frequencies for the levels.

    Returns
    -------
    p_recomb_emission
        DataFrame containing unnormalized recombination emission probabilities.
    """
    p_recomb_emission = spontaneous_recombination_coeff.multiply(
        energies_diff_bound_free.values, axis=0
    )

    return p_recomb_emission


def continuum_transition_recombination_emission(
    spontaneous_recombination_coeff: pd.DataFrame,
    energies_diff_bound_free: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate unnormalized probabilities and metadata for radiative recombination emission transitions (continuum).

    Parameters
    ----------
    spontaneous_recombination_coeff
        Rate coefficient for spontaneous recombination from `k` to level `i`.
    energies_diff_bound_free
        Ionization threshold frequencies for the levels.

    Returns
    -------
    p_recomb_emission
        DataFrame containing unnormalized recombination emission probabilities.
    recombination_emission_metadata
        DataFrame containing metadata for the recombination emission transitions.
    """
    p_recomb_emission = probability_recombination_emission(
        spontaneous_recombination_coeff, energies_diff_bound_free
    )

    destinations = p_recomb_emission.index.values
    sources = get_ground_state_multi_index(p_recomb_emission.index).values

    recomb_emission_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.RECOMB_EMISSION,
            "transition_line_idx": -99,
            "photoionization_key_idx": range(len(energies_diff_bound_free)),
            "collision_key_idx": -99,
        },
        index=p_recomb_emission.index,
    )

    return p_recomb_emission, recomb_emission_metadata


def probability_photoionization_internal(
    stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
    photoionization_data_level_energies: pd.Series,
) -> pd.DataFrame:
    """
    Calculate photoionization probability unnormalized.

    Parameters
    ----------
    stim_recomb_corrected_photoionization_rate_coeff
        Corrected photoionization rate coefficient from level `i` to `k`.
    photoionization_data_level_energies
        Energies of the levels involved in photoionization.

    Returns
    -------
    p_photoion_int
        DataFrame containing photoionization probabilities.
    """
    p_photoion_int = stim_recomb_corrected_photoionization_rate_coeff.multiply(
        photoionization_data_level_energies, axis=0
    )
    return p_photoion_int


def continuum_transition_photoionization_internal(
    stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
    photoionization_data_level_energies: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate photoionization probability unnormalized.

    Parameters
    ----------
    stim_recomb_corrected_photoionization_rate_coeff
        Corrected photoionization rate coefficient from level `i` to `k`.
    photoionization_data_level_energies
        Energies of the levels involved in photoionization.

    Returns
    -------
    p_photoion_int
        DataFrame containing photoionization probabilities.
    photoion_internal_metadata
        DataFrame containing metadata for the photoionization transitions.
    """
    p_photoion_int = probability_photoionization_internal(
        stim_recomb_corrected_photoionization_rate_coeff,
        photoionization_data_level_energies,
    )

    sources = p_photoion_int.index.values
    destinations = get_ground_state_multi_index(p_photoion_int.index).values

    photoion_internal_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.PHOTOIONIZATION_INTERNAL,
            "transition_line_idx": -99,
            "photoionization_key_idx": range(
                len(photoionization_data_level_energies)
            ),
            "collision_key_idx": -99,
        },
        index=p_photoion_int.index,
    )

    return p_photoion_int, photoion_internal_metadata


def probability_photoionization_to_k_packet(
    stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
    energies_diff_bound_free: pd.Series,
) -> pd.DataFrame:
    """
    Calculate photoionization probability unnormalized.

    Parameters
    ----------
    stim_recomb_corrected_photoionization_rate_coeff
        Corrected photoionization rate coefficient from level `i` to `k`.
    energies_diff_bound_free
        Energies of the levels involved in photoionization.

    Returns
    -------
    p_photoion_to_k_packet
        DataFrame containing photoionization probabilities.
    """
    p_photoion_to_k_packet = (
        stim_recomb_corrected_photoionization_rate_coeff.multiply(
            energies_diff_bound_free.values, axis=0
        )
    )
    return p_photoion_to_k_packet


def continuum_transition_photoionization_to_k_packet(
    stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
    energies_bound_free_lower: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate photoionization probability unnormalized.

    Parameters
    ----------
    stim_recomb_corrected_photoionization_rate_coeff
        Corrected photoionization rate coefficient from level `i` to `k`.
    energies_bound_free_lower
        Energies of the levels involved in photoionization.

    Returns
    -------
    p_photoion_to_k_packet
        DataFrame containing photoionization probabilities.
    photoion_to_k_packet_metadata
        DataFrame containing metadata for the photoionization transitions.
    """
    p_photoion_to_k_packet = probability_photoionization_to_k_packet(
        stim_recomb_corrected_photoionization_rate_coeff,
        energies_bound_free_lower,
    )

    sources = p_photoion_to_k_packet.index.values
    destinations = get_ground_state_multi_index(
        p_photoion_to_k_packet.index
    ).values

    photoion_to_k_packet_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.PHOTOIONIZATION_TO_K_PACKET,
            "transition_line_idx": -99,
            "photoionization_key_idx": range(len(energies_bound_free_lower)),
            "collision_key_idx": -99,
        },
        index=p_photoion_to_k_packet.index,
    )

    return p_photoion_to_k_packet, photoion_to_k_packet_metadata


def probability_adiabatic_cooling(
    electron_densities: pd.DataFrame,
    t_electrons: pd.Series,
    time_explosion: float,
) -> pd.DataFrame:
    """
    Calculate the adiabatic cooling rate (unnormalized).

    This function computes a simple adiabatic cooling term used as an
    unnormalized probability for macro-atom transitions that represent
    adiabatic losses.

    Parameters
    ----------
    electron_densities
        Electron number densities (per shell / spatial zone).
    t_electrons
        Electron temperatures (same index/shape as electron_densities).
    time_explosion
        Time since explosion (seconds).

    Returns
    -------
    p_adiabatic_cool
        Unnormalized adiabatic cooling rates with the same index as
        `electron_densities` / `t_electrons`.
    """
    # Calculate the probability (unnormalized cooling rate)
    adiabatic_cooling_rate = (
        3.0 * electron_densities * K_B * t_electrons
    ) / time_explosion
    raise NotImplementedError  # Adiabatic cooling rate needs to be transformed to probability
    # return p_adiabatic_cooling


def continuum_adiabatic_cooling(
    electron_densities: pd.Series,
    t_electrons: pd.Series,
    time_explosion: float,
) -> tuple[pd.Series, pd.DataFrame]:
    """
    Calculate the adiabatic cooling rate.

    Parameters
    ----------
    electron_densities
        Electron densities.
    t_electrons
        Electron temperatures.
    time_explosion
        Time since explosion.

    Returns
    -------
    p_adiabatic_cool
        Unnormalized adiabatic cooling rates.
    adiabatic_cool_metadata
        Metadata DataFrame describing the cooling transitions.
    """
    p_adiabatic_cool = probability_adiabatic_cooling(
        electron_densities, t_electrons, time_explosion
    )
    sources = [("k", -99, -99)] * len(
        p_adiabatic_cool.index
    )  # k packets could use a better convention
    destinations = [("adiabatic", -99, -99)] * len(p_adiabatic_cool.index)
    adiabatic_cool_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,  # there are adiabatic, ff, bf
            "transition_type": MacroAtomTransitionType.ADIABATIC_COOLING,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": -99,
        },
        index=p_adiabatic_cool.index,
    )

    return p_adiabatic_cool, adiabatic_cool_metadata


def probability_free_free_cooling(
    ion_number_density: pd.Series | pd.DataFrame,
    electron_densities: pd.Series,
    t_electrons: pd.Series,
) -> pd.Series:
    """
    Calculate an unnormalized free-free (bremsstrahlung) cooling probability.

    This routine estimates a free-free cooling rate used as an unnormalized
    probability for macro-atom cooling transitions. The implementation uses
    an approximate scaling with electron temperature and ion charge.

    Parameters
    ----------
    ion_number_density
        Index or Series containing ion number densities (levelled by ion).
    electron_densities
        Electron number densities per zone.
    t_electrons
        Electron temperatures per zone.

    Returns
    -------
    p_free_free_cool
        Unnormalized free-free cooling rates with the same index as
        `electron_densities` / `t_electrons`.
    """
    ion_charge = ion_number_density.index.get_level_values(1).values
    cooling_factor = (
        electron_densities
        * ion_number_density.multiply(ion_charge**2, axis=0).sum()
    )
    ff_cool_rate = F_K * np.sqrt(t_electrons) * cooling_factor
    # NOT DONE, return probability, not the rate
    raise NotImplementedError
    # return probability_free_free_cooling


def continuum_free_free_cooling(
    ion_number_density: pd.Series | pd.DataFrame,
    electron_densities: pd.Series,
    t_electrons: pd.Series,
) -> tuple[pd.Series, pd.DataFrame]:
    """
    Wrap free-free cooling rates with metadata for macro-atom use.
    NOTE: This I believe is a deactivation transition. It is not implemented or hooked up to the solver yet.

    Parameters
    ----------
    ion_number_density
        Index or Series containing ion number densities (levelled by ion).
    electron_densities
        Electron number densities per zone.
    t_electrons
        Electron temperatures per zone.

    Returns
    -------
    p_free_free_cool
        Unnormalized free-free cooling rates.
    free_free_cool_metadata
        Metadata DataFrame describing the cooling transitions. The index
        matches `p_free_free_cool` and includes fields such as
        `transition_type` and `photoionization_key_idx`.
    """
    p_free_free_cool = probability_free_free_cooling(
        ion_number_density, electron_densities, t_electrons
    )
    sources = [("k", -99, -99)] * len(p_free_free_cool.index)
    destinations = [("free_free", -99, -99)] * len(p_free_free_cool.index)
    free_free_cool_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.FF_EMISSION,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": -99,
        },
        index=p_free_free_cool.index,
    )
    return p_free_free_cool, free_free_cool_metadata


# Collisional transitions below


def probability_collision_deexc_to_k_packet(
    coll_deexc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    delta_E_yg: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional de-excitation to k packet probabilities.

    Parameters
    ----------
    coll_deexc_coeff
        Collisional de-excitation coefficients.
    electron_densities
        Electron densities.
    delta_E_yg
        Energy differences.

    Returns
    -------
    p_coll_down_to_k_packet
        DataFrame containing collisional de-excitation to k packet probabilities.
    """
    p_coll_down_to_k_packet = (coll_deexc_coeff * electron_densities).multiply(
        delta_E_yg, axis=0
    )

    return p_coll_down_to_k_packet


def collisional_transition_deexc_to_k_packet(
    coll_deexc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    delta_E_yg: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create collisional de-excitation and deactivation transitions to a 'k' packet.

    This function wraps the `probability_collision_deexc_to_k_packet` probability
    calculation and produces a metadata DataFrame describing the resulting
    de-excitation transitions that deposit energy into the k-packet
    channel.

    Parameters
    ----------
    coll_deexc_coeff
        Collisional de-excitation coefficients indexed by transition.
    electron_densities
        Electron number densities per zone.
    delta_E_yg
        Energy differences for the transitions.

    Returns
    -------
    p_coll_down_to_k_packet
        Unnormalized probabilities (rates) for de-excitation/deactivation.
    coll_down_to_k_packet_metadata
        Metadata for the de-excitation transitions; index matches
        `p_coll_down_to_k_packet`.
    """
    p_coll_down_to_k_packet = probability_collision_deexc_to_k_packet(
        coll_deexc_coeff, electron_densities, delta_E_yg
    )

    sources = p_coll_down_to_k_packet.index.droplevel(
        "level_number_lower"
    ).values
    destinations = [("k", -99, -99)] * len(p_coll_down_to_k_packet.index)

    coll_down_to_k_packet_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_DOWN_TO_K_PACKET,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": range(len(coll_deexc_coeff)),
        },
        index=p_coll_down_to_k_packet.index,
    )

    return p_coll_down_to_k_packet, coll_down_to_k_packet_metadata


def probability_collision_internal_down(
    coll_deexc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_lower: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional internal de-excitation probabilities.

    Parameters
    ----------
    coll_deexc_coeff
        Collisional de-excitation coefficients.
    electron_densities
        Electron densities.
    energies_lower
        Energy values of the lower levels.

    Returns
    -------
    p_coll_internal_down
        DataFrame containing collisional de-excitation probabilities.
    """
    p_coll_internal_down = (coll_deexc_coeff * electron_densities).multiply(
        energies_lower.values, axis=0
    )
    return p_coll_internal_down


# Might not be a used transition
def collisional_transition_internal_down(
    coll_deexc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_lower: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build collisional internal-down transition probabilities and metadata.

    Parameters
    ----------
    coll_deexc_coeff
        Collisional de-excitation coefficients indexed by transitions.
    electron_densities
        Electron number densities per zone.
    energies_lower
        Energy values of the lower levels.

    Returns
    -------
    p_coll_internal_down
        Unnormalized collisional internal down transition probabilities.
    coll_internal_down_metadata
        Metadata for the collisional internal down transitions.
    """
    sources = coll_deexc_coeff.index.droplevel("level_number_lower").values
    destinations = coll_deexc_coeff.index.droplevel("level_number_upper").values
    p_coll_internal_down = probability_collision_internal_down(
        coll_deexc_coeff, electron_densities, energies_lower
    )

    coll_internal_down_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_DOWN_INTERNAL,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": range(len(coll_deexc_coeff)),
        },
        index=p_coll_internal_down.index,
    )
    return p_coll_internal_down, coll_internal_down_metadata


def probability_collision_exc_internal(
    coll_exc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_lower: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional internal up probabilities.

    Parameters
    ----------
    coll_exc_coeff
        Collisional excitation coefficients.
    electron_densities
        Electron densities.
    energies_lower
        Energy values of the lower levels.

    Returns
    -------
    p_coll_exc_internal
        DataFrame containing collisional internal up probabilities.
    """
    p_coll_exc_internal = (coll_exc_coeff * electron_densities).multiply(
        energies_lower.values, axis=0
    )

    return p_coll_exc_internal


def collisional_transition_excitation_internal(
    coll_exc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_lower: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build collisional internal-up transition probabilities and metadata.

    Parameters
    ----------
    coll_exc_coeff
        Collisional excitation coefficients indexed by transitions.
    electron_densities
        Electron number densities per zone.
    energies_lower
        Energy values of the lower levels.

    Returns
    -------
    p_coll_up_internal
        Unnormalized collisional excitation to macro transition probabilities.
    coll_up_internal_metadata
        Metadata for the collisional internal up transitions.
    """
    sources = coll_exc_coeff.index.droplevel("level_number_upper").values
    destinations = coll_exc_coeff.index.droplevel("level_number_lower").values

    p_coll_up_internal = probability_collision_exc_internal(
        coll_exc_coeff, electron_densities, energies_lower
    )
    coll_up_internal_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_UP_INTERNAL,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": range(len(coll_exc_coeff)),
        },
        index=p_coll_up_internal.index,
    )
    return p_coll_up_internal, coll_up_internal_metadata


def probability_collision_excitation_cool(
    coll_exc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    delta_E_yg: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional excitation cooling rates (unnormalized) and aggregate by destination level.

    Parameters
    ----------
    coll_exc_coeff
        Collisional excitation coefficients indexed by transitions.
    electron_densities
        Electron number densities per zone.
    delta_E_yg
        Energy differences for the excitation transitions (level -> g).

    Returns
    -------
    p_coll_excitation_cool
        Aggregated collisional excitation cooling rates grouped by
        `destination_level_idx`.
    """
    p_coll_excitation_cool = (coll_exc_coeff * electron_densities).multiply(
        delta_E_yg, axis=0
    )

    p_coll_excitation_cool = p_coll_excitation_cool.groupby(
        level=["atomic_number", "ion_number", "level_number_upper"]
    ).sum()
    return p_coll_excitation_cool


def collisional_transition_excitation_cool(
    coll_exc_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    delta_E_yg: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build collisional excitation cooling transitions and metadata.

    This function computes the collisional excitation cooling rates (an
    unnormalized probability) and packages them together with metadata
    describing the cooling transitions.

    Parameters
    ----------
    coll_exc_coeff
        Collisional excitation coefficients indexed by transitions.
    electron_densities
        Electron number densities per zone.
    delta_E_yg
        Energy differences for the excitations.

    Returns
    -------
    p_coll_excitation_cool
        Unnormalized collisional excitation cooling rates.
    coll_excitation_cool_metadata
        Metadata for the excitation cooling transitions.
    """
    p_coll_excitation_cool = probability_collision_excitation_cool(
        coll_exc_coeff,
        electron_densities,
        delta_E_yg,
    )
    sources = [("k", -99, -99)] * len(p_coll_excitation_cool)
    destinations = p_coll_excitation_cool.index.values
    coll_excitation_cool_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_EXC_COOL_TO_MACRO,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": -99,  # Collapsed along level_number_lower, so can't be mapped to collision keys
        },
        index=p_coll_excitation_cool.index,
    )

    return p_coll_excitation_cool, coll_excitation_cool_metadata


def probability_collision_ionization_internal(
    coll_ion_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_coll_ion: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional ionization internal probabilities.

    Parameters
    ----------
    coll_ion_coeff
        Collisional ionization coefficients.
    electron_densities
        Electron number densities.
    energies_coll_ion
        Ionization energy terms.

    Returns
    -------
    p_coll_ionization_internal
        Unnormalized collisional ionization internal probabilities.
    """
    p_coll_ionization_internal = (coll_ion_coeff * electron_densities).multiply(
        energies_coll_ion, axis=0
    )

    return p_coll_ionization_internal


def collisional_transition_ionization_internal(
    coll_ion_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_coll_ion: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build collisional ionization internal probabilities and metadata.

    Parameters
    ----------
    coll_ion_coeff
        Collisional ionization coefficients.
    electron_densities
        Electron number densities.
    energies_coll_ion
        Ionization energy terms.

    Returns
    -------
    p_coll_ionization_internal
        Unnormalized collisional ionization internal probabilities.
    coll_ionzation_internal_metadata
        Metadata for collisional ionization internal transitions.
    """
    p_coll_ionization_internal = probability_collision_ionization_internal(
        coll_ion_coeff,
        electron_densities,
        energies_coll_ion,
    )
    sources = coll_ion_coeff.index.values
    destinations = [("i", -99, -99)] * len(p_coll_ionization_internal)
    coll_ionzation_internal_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_ION_INTERNAL,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": -99,  # Collapsed along level_number_lower, so can't be mapped to collision keys
        },
        index=p_coll_ionization_internal.index,
    )

    return p_coll_ionization_internal, coll_ionzation_internal_metadata


def probability_collision_ionization_emission(
    coll_ion_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_diff_bound_free: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional ionization emission probabilities.

    Parameters
    ----------
    coll_ion_coeff
        Collisional ionization coefficients.
    electron_densities
        Electron number densities.
    energies_diff_bound_free
        Bound-free energy differences.

    Returns
    -------
    p_coll_ionization_internal
        Unnormalized collisional ionization emission probabilities.
    """
    p_coll_ionization_internal = (coll_ion_coeff * electron_densities).multiply(
        energies_diff_bound_free.values, axis=0
    )

    return p_coll_ionization_internal


def collisional_transition_ionization_emission(
    coll_ion_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_diff_bound_free: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build collisional ionization emission probabilities and metadata.

    Parameters
    ----------
    coll_ion_coeff
        Collisional ionization coefficients.
    electron_densities
        Electron number densities.
    energies_diff_bound_free
        Bound-free energy differences.

    Returns
    -------
    p_coll_ionization_emission
        Unnormalized collisional ionization emission probabilities.
    coll_ionization_emission_metadata
        Metadata for collisional ionization emission transitions.
    """
    p_coll_ionization_emission = probability_collision_ionization_emission(
        coll_ion_coeff,
        electron_densities,
        energies_diff_bound_free,
    )
    sources = coll_ion_coeff.index.values
    destinations = [("i", -99, -99)] * len(
        p_coll_ionization_emission
    )  # Maybe supposed to go to K block?
    coll_ionization_emission_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_ION_EMISSION,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": -99,
        },
        index=p_coll_ionization_emission.index,
    )

    return p_coll_ionization_emission, coll_ionization_emission_metadata


def probability_collision_recombination_internal(
    coll_recomb_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_coll_ion: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional recombination internal probabilities.

    Parameters
    ----------
    coll_recomb_coeff
        Collisional recombination coefficients.
    electron_densities
        Electron number densities.
    energies_coll_ion
        Energy terms.

    Returns
    -------
    p_coll_recomb_internal
        Unnormalized collisional recombination internal probabilities.
    """
    p_coll_ionization_internal = (
        coll_recomb_coeff * electron_densities
    ).multiply(energies_coll_ion.values, axis=0)

    return p_coll_ionization_internal


def collisional_transition_recombination_internal(
    coll_recomb_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_coll_ion: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build collisional recombination internal probabilities and metadata.

    Parameters
    ----------
    coll_recomb_coeff
        Collisional recombination coefficients.
    electron_densities
        Electron number densities.
    energies_coll_ion
        Energy terms.

    Returns
    -------
    p_coll_recomb_internal
        Unnormalized collisional recombination internal probabilities.
    coll_recomb_internal
        Metadata for collisional recombination internal transitions.
    """
    p_coll_recomb_internal = probability_collision_recombination_internal(
        coll_recomb_coeff,
        electron_densities,
        energies_coll_ion,
    )
    sources = [("i", -99, -99)] * len(
        p_coll_recomb_internal
    )  # Might also be from K?
    destinations = coll_recomb_coeff.index.values
    coll_recomb_internal = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_RECOMB_INTERNAL,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": -99,
        },
        index=p_coll_recomb_internal.index,
    )

    return p_coll_recomb_internal, coll_recomb_internal


def probability_collision_recombination_emission(
    coll_recomb_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_diff_bound_free: pd.Series,
) -> pd.DataFrame:
    """
    Calculate collisional recombination emission probabilities.

    Parameters
    ----------
    coll_recomb_coeff
        Collisional recombination coefficients.
    electron_densities
        Electron number densities.
    energies_diff_bound_free
        Bound-free energy differences.

    Returns
    -------
    p_coll_recomb_emission
        Unnormalized collisional recombination emission probabilities.
    """
    p_coll_ionization_internal = (
        coll_recomb_coeff * electron_densities
    ).multiply(energies_diff_bound_free.values, axis=0)

    return p_coll_ionization_internal


def collisional_transition_recombination_emission(
    coll_recomb_coeff: pd.DataFrame,
    electron_densities: pd.Series,
    energies_diff_bound_free: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build collisional recombination emission probabilities and metadata.

    Parameters
    ----------
    coll_recomb_coeff
        Collisional recombination coefficients.
    electron_densities
        Electron number densities.
    energies_diff_bound_free
        Bound-free energy differences.

    Returns
    -------
    p_coll_recomb_emission
        Unnormalized collisional recombination emission probabilities.
    coll_recomb_emission_metadata
        Metadata for collisional recombination emission transitions.
    """
    p_coll_recomb_emission = probability_collision_recombination_emission(
        coll_recomb_coeff,
        electron_densities,
        energies_diff_bound_free,
    )
    sources = [("i", -99, -99)] * len(p_coll_recomb_emission)
    destinations = coll_recomb_coeff.index.values  # Not sure this is right
    coll_recomb_emission_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_RECOMB_EMISSION,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": -99,
        },
        index=p_coll_recomb_emission.index,
    )

    return p_coll_recomb_emission, coll_recomb_emission_metadata
