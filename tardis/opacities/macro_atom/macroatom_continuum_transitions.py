import numpy as np
import pandas as pd

from tardis.opacities.macro_atom.macroatom_line_transitions import (
    CONST_H_CGS,
)
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
            "photoionization_key_idx": range(
                len(photoionization_data_level_energies)
            ),
            "collision_key_idx": -99,
        },
        index=p_recomb_internal.index,
    )

    p_recomb_internal["source"] = sources

    return p_recomb_internal, recombination_internal_metadata


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
    # The above comment comes from the source code, which was a plasma property
    return p_recomb_emission


def probability_recombination_cooling(
    spontaneous_recombination_cooling_coeff,
    electron_densities,
    ion_number_density,
):
    # Check BoundFreeThermalRates in plasma.equilibrium.rates.heating_cooling_rates
    # Implement from FreeBoundCoolingRate in plasma.properties.continuum_processes.rates
    # This will use MacroAtomTransitionType.BF_COOLING
    raise NotImplementedError


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
    ## TODO: FIX
    sources = get_ground_state_multi_index(
        p_recomb_emission.index
    ).values  # Should be k packet?

    unique_idx = photoionization_data_frequencies.index.unique()
    mapping = {idx: i for i, idx in enumerate(unique_idx)}
    ids = pd.Series(
        photoionization_data_frequencies.index.map(mapping),
        index=photoionization_data_frequencies.index,
    ).astype(int)

    recombination_emission_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.BF_EMISSION,
            "transition_line_idx": ids,  # this is used to indicate which block of the photoionization data the transition relates to - make it a different column later
            "photoionization_key_idx": range(
                len(photoionization_data_frequencies)
            ),
            "collision_key_idx": -99,
        },
        index=p_recomb_emission.index,
    )

    p_recomb_emission["source"] = sources

    return p_recomb_emission, recombination_emission_metadata


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
            "photoionization_key_idx": range(
                len(photoionization_data_level_energies)
            ),
            "collision_key_idx": -99,
        },
        index=p_photoionization.index,
    )

    p_photoionization["source"] = sources

    return p_photoionization, photoionization_metadata


def probability_adiabatic_cooling(
    electron_densities, t_electrons, time_explosion
):
    """
    Calculate the adiabatic cooling rate (unnormalized).

    This function computes a simple adiabatic cooling term used as an
    unnormalized probability for macro-atom transitions that represent
    adiabatic losses.

    Parameters
    ----------
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities (per shell / spatial zone).
    t_electrons : pd.DataFrame or pd.Series
        Electron temperatures (same index/shape as electron_densities).
    time_explosion : float
        Time since explosion (seconds).

    Returns
    -------
    pd.DataFrame or pd.Series
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
    electron_densities, t_electrons, time_explosion
):
    """
    Calculate the adiabatic cooling rate.

    Parameters
    ----------
    electron_densities : pd.DataFrame
        Electron densities.
    t_electrons : pd.DataFrame
        Electron temperatures.
    time_explosion : float
        Time since explosion.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the adiabatic cooling rates.
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
    ion_number_density, electron_densities, t_electrons
):
    """
    Calculate an unnormalized free-free (bremsstrahlung) cooling probability.

    This routine estimates a free-free cooling rate used as an unnormalized
    probability for macro-atom cooling transitions. The implementation uses
    an approximate scaling with electron temperature and ion charge.

    Parameters
    ----------
    ion_number_density : pd.MultiIndex or pd.DataFrame index
        Index or Series containing ion number densities (levelled by ion).
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities per zone.
    t_electrons : pd.DataFrame or pd.Series
        Electron temperatures per zone.

    Returns
    -------
    pd.DataFrame or pd.Series
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
    ion_number_density, electron_densities, t_electrons
):
    """
    Wrap free-free cooling rates with metadata for macro-atom use.
    NOTE: This I believe is a deactivation transition. It is not implemented or hooked up to the solver yet.

    Parameters
    ----------
    ion_number_density : pd.MultiIndex or pd.DataFrame index
        Index or Series containing ion number densities (levelled by ion).
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities per zone.
    t_electrons : pd.DataFrame or pd.Series
        Electron temperatures per zone.

    Returns
    -------
    p_free_free_cool : pd.DataFrame or pd.Series
        Unnormalized free-free cooling rates.
    free_free_cool_metadata : pd.DataFrame
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
    coll_deexc_coeff, electron_densities, delta_E_yg
):
    """
    Calculate collisional de-excitation to k packet probabilities.

    Parameters
    ----------
    coll_deexc_coeff : pd.DataFrame
        Collisional de-excitation coefficients.
    electron_densities : pd.DataFrame
        Electron densities.
    delta_E_yg : pd.Series
        Energy differences.

    Returns
    -------
    pd.DataFrame
        DataFrame containing collisional de-excitation to k packet probabilities.
    """
    p_coll_down_to_k_packet = (coll_deexc_coeff * electron_densities).multiply(
        delta_E_yg.values, axis=0
    )

    return p_coll_down_to_k_packet


def collisional_transition_deexc_to_k_packet(
    coll_deexc_coeff, electron_densities, delta_E_yg
):
    """
    Create collisional de-excitation and deactivation transitions to a 'k' packet.

    This function wraps the `probability_collision_deexc_to_k_packet` probability
    calculation and produces a metadata DataFrame describing the resulting
    de-excitation transitions that deposit energy into the k-packet
    channel.

    Parameters
    ----------
    coll_deexc_coeff : pd.DataFrame
        Collisional de-excitation coefficients indexed by transition.
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities per zone.
    delta_E_yg : pd.Series
        Energy differences for the transitions.

    Returns
    -------
    p_coll_down_to_k_packet : pd.DataFrame
        Unnormalized probabilities (rates) for de-excitation/deactivation.
    coll_down_to_k_packet_metadata : pd.DataFrame
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
    coll_deexc_coeff, electron_densities, energies_lowers
):
    """
    Calculate collisional internal de-excitation probabilities.

    Parameters
    ----------
    coll_deexc_coeff : pd.DataFrame
        Collisional de-excitation coefficients.
    electron_densities : pd.DataFrame
        Electron densities.
    energies_lowers : pd.Series
        Energy values of the lower levels.

    Returns
    -------
    pd.DataFrame
        DataFrame containing collisional de-excitation probabilities.
    """
    p_coll_internal_down = (coll_deexc_coeff * electron_densities).multiply(
        energies_lowers.values, axis=0
    )
    return p_coll_internal_down


# Might not be a used transition
def collisional_transition_internal_down(
    coll_deexc_coeff, electron_densities, energies_lowers
):
    """
    Build collisional internal-down transition probabilities and metadata.

    Parameters
    ----------
    coll_deexc_coeff : pd.DataFrame
        Collisional de-excitation coefficients indexed by transitions.
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities per zone.
    energies_lowers : pd.Series
        Energy values of the lower levels.

    Returns
    -------
    p_coll_internal_down : pd.DataFrame
        Unnormalized collisional internal down transition probabilities.
    coll_internal_down_metadata : pd.DataFrame
        Metadata for the collisional internal down transitions.
    """
    sources = coll_deexc_coeff.index.droplevel("level_number_lower").values
    destinations = coll_deexc_coeff.index.droplevel("level_number_upper").values
    p_coll_internal_down = probability_collision_internal_down(
        coll_deexc_coeff, electron_densities, energies_lowers
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


def probability_collision_exc_to_macro(
    coll_exc_coeff, electron_densities, energy_lowers
):
    """
    Calculate collisional internal up probabilities.

    Parameters
    ----------
    coll_exc_coeff : pd.DataFrame
        Collisional excitation coefficients.
    electron_densities : pd.DataFrame
        Electron densities.
    energy_lowers : pd.Series
        Energy values of the lower levels.

    Returns
    -------
    pd.DataFrame
        DataFrame containing collisional internal up probabilities.
    """
    p_coll_exc_to_macro = (coll_exc_coeff * electron_densities).multiply(
        energy_lowers.values, axis=0
    )

    return p_coll_exc_to_macro


def collisional_transition_excitation_to_macro_atom(
    coll_exc_coeff, electron_densities, energies_lowers
):
    """
    Build collisional internal-up transition probabilities and metadata.

    Parameters
    ----------
    coll_exc_coeff : pd.DataFrame
        Collisional excitation coefficients indexed by transitions.
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities per zone.
    energies_lowers : object
        Atomic data object providing `levels` with energy values.

    Returns
    -------
    p_coll_exc_to_macro : pd.DataFrame
        Unnormalized collisional excitation to macro transition probabilities.
    coll_excite_to_macro_metadata : pd.DataFrame
        Metadata for the collisional internal up transitions.
    """
    sources = coll_exc_coeff.index.droplevel("level_number_upper").values
    destinations = coll_exc_coeff.index.droplevel("level_number_lower").values

    p_coll_exc_to_macro = probability_collision_exc_to_macro(
        coll_exc_coeff, electron_densities, energies_lowers
    )
    coll_excite_to_macro_metadata = pd.DataFrame(
        {
            "transition_line_id": -99,
            "source": sources,
            "destination": destinations,
            "transition_type": MacroAtomTransitionType.COLL_EXC_COOL_TO_MACRO,
            "transition_line_idx": -99,
            "photoionization_key_idx": -99,
            "collision_key_idx": range(len(coll_exc_coeff)),
        },
        index=p_coll_exc_to_macro.index,
    )
    return p_coll_exc_to_macro, coll_excite_to_macro_metadata


def probability_collision_excitation_cool(
    coll_exc_coeff,
    electron_densities,
    delta_E_yg,
    level_number_density,
    lower_indices,
):
    """
    Calculate collisional excitation cooling rates (unnormalized) and aggregate by destination level.

    Parameters
    ----------
    coll_exc_coeff : pd.DataFrame
        Collisional excitation coefficients indexed by transitions.
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities per zone.
    delta_E_yg : pd.Series
        Energy differences for the excitation transitions (level -> g).
    level_number_density : pd.Series
        Number density for the lower levels involved in the transitions.
    lower_indices : pd.Index
        Index of lower level identifiers that aligns `level_number_density` with
        the rows of `coll_exc_coeff`.

    Returns
    -------
    pd.DataFrame or pd.Series
        Aggregated collisional excitation cooling rates grouped by
        `destination_level_idx`.
    """
    p_coll_excitation_cool = (coll_exc_coeff * electron_densities).multiply(
        delta_E_yg.values, axis=0
    ) * level_number_density.loc[lower_indices].values

    p_coll_excitation_cool = p_coll_excitation_cool.groupby(
        level=["atomic_number", "ion_number", "level_number_upper"]
    ).sum()
    return p_coll_excitation_cool


def collisional_transition_excitation_cool(
    coll_exc_coeff,
    electron_densities,
    delta_E_yg,
    level_number_density,
    lower_indices,
):
    """
    Build collisional excitation cooling transitions and metadata.

    This function computes the collisional excitation cooling rates (an
    unnormalized probability) and packages them together with metadata
    describing the cooling transitions.

    Parameters
    ----------
    coll_exc_coeff : pd.DataFrame
        Collisional excitation coefficients indexed by transitions.
    electron_densities : pd.DataFrame or pd.Series
        Electron number densities per zone.
    delta_E_yg : pd.Series
        Energy differences for the excitations.
    level_number_density : pd.Series
        Number density for the lower levels involved in the transitions.
    lower_indices : pd.Index
        Index of lower level identifiers that aligns `level_number_density` with
        the rows of `coll_exc_coeff`.

    Returns
    -------
    p_coll_excitation_cool : pd.DataFrame
        Unnormalized collisional excitation cooling rates.
    coll_excitation_cool_metadata : pd.DataFrame
        Metadata for the excitation cooling transitions.
    """
    p_coll_excitation_cool = probability_collision_excitation_cool(
        coll_exc_coeff,
        electron_densities,
        delta_E_yg,
        level_number_density,
        lower_indices,
    )
    sources = [("k", -99, -99)] * len(p_coll_excitation_cool)
    destinations = (
        coll_exc_coeff.index.droplevel(["level_number_lower"]).unique().values
    )
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
