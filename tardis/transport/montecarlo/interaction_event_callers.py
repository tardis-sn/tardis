import numpy as np
from numba import njit

import tardis.transport.montecarlo.configuration.montecarlo_globals as montecarlo_globals
from tardis.transport.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.macro_atom import (
    MacroAtomTransitionType,
    macro_atom_interaction,
)
from tardis.transport.montecarlo.interaction_event_handlers import (
    LineInteractionType,
    adiabatic_cooling,
    bf_cooling,
    bound_free_emission,
    determine_bf_macro_activation_idx,
    free_free_emission,
    line_emission,
    recombination_emission,
)
from tardis.transport.montecarlo.utils import get_random_mu


@njit(**njit_dict_no_parallel)
def macro_atom_event(
    destination_level_idx,
    r_packet,
    time_explosion,
    opacity_state,
    enable_full_relativity,
):
    """
    Macroatom event handler - run the macroatom and handle the result

    Parameters
    ----------
    destination_level_idx : int
    r_packet : tardis.transport.montecarlo.r_packet.RPacket
    time_explosion : float
    opacity_state : tardis.transport.montecarlo.numba_interface.OpacityState
    """
    transition_id, transition_type = macro_atom_interaction(
        destination_level_idx, r_packet.current_shell_id, opacity_state
    )

    if (
        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.FF_EMISSION
    ):
        free_free_emission(
            r_packet, time_explosion, opacity_state, enable_full_relativity
        )

    elif (
        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.BF_EMISSION
    ):
        bound_free_emission(
            r_packet,
            time_explosion,
            opacity_state,
            transition_id,
            enable_full_relativity,
        )
    elif (
        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.BF_COOLING
    ):
        bf_cooling(
            r_packet, time_explosion, opacity_state, enable_full_relativity
        )

    elif (
        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.ADIABATIC_COOLING
    ):
        adiabatic_cooling(r_packet)

    elif transition_type == MacroAtomTransitionType.BB_EMISSION:
        line_emission(
            r_packet,
            transition_id,
            time_explosion,
            opacity_state,
            enable_full_relativity,
        )
    elif transition_type == MacroAtomTransitionType.RECOMB_EMISSION:
        recombination_emission()
    else:
        raise Exception(f"Interaction {transition_type} not Found!")


@njit(**njit_dict_no_parallel)
def determine_continuum_macro_activation_idx(
    opacity_state, nu, chi_bf, chi_ff, chi_bf_contributions, active_continua
):
    """
    Determine the macro atom activation level after a continuum absorption.

    Parameters
    ----------
    nu : float
        Comoving frequency of the r-packet.
    chi_bf : numpy.ndarray, dtype float
        Bound-free opacity.
    chi_bf : numpy.ndarray, dtype float
        Free-free opacity.
    chi_bf_contributions : numpy.ndarray, dtype float
        Cumulative distribution of bound-free opacities at frequency
        `nu`.
    active_continua : numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.

    Returns
    -------
    float
        Macro atom activation idx.
    """
    fraction_bf = chi_bf / (chi_bf + chi_ff)
    # TODO: In principle, we can also decide here whether a Thomson
    # scattering event happens and need one less RNG call.
    if np.random.random() < fraction_bf:  # Bound-free absorption
        destination_level_idx = determine_bf_macro_activation_idx(
            opacity_state, nu, chi_bf_contributions, active_continua
        )
    else:  # Free-free absorption (i.e. k-packet creation)
        destination_level_idx = opacity_state.k_packet_idx
    return destination_level_idx


@njit(**njit_dict_no_parallel)
def continuum_event(
    r_packet,
    time_explosion,
    opacity_state,
    chi_bf_tot,
    chi_ff,
    chi_bf_contributions,
    current_continua,
    enable_full_relativity,
):
    """
    continuum event handler - activate the macroatom and run the handler

    Parameters
    ----------
    r_packet : tardis.transport.montecarlo.r_packet.RPacket
    time_explosion : float
    opacity_state : tardis.transport.montecarlo.numba_interface.OpacityState
    continuum : tardis.transport.montecarlo.numba_interface.Continuum
    """
    old_doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion, enable_full_relativity
    )

    r_packet.mu = get_random_mu()
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion, enable_full_relativity
    )
    comov_energy = r_packet.energy * old_doppler_factor
    comov_nu = (
        r_packet.nu * old_doppler_factor
    )  # make sure frequency should be updated
    r_packet.energy = comov_energy * inverse_doppler_factor
    r_packet.nu = comov_nu * inverse_doppler_factor

    destination_level_idx = determine_continuum_macro_activation_idx(
        opacity_state,
        comov_nu,
        chi_bf_tot,
        chi_ff,
        chi_bf_contributions,
        current_continua,
    )

    macro_atom_event(
        destination_level_idx,
        r_packet,
        time_explosion,
        opacity_state,
        enable_full_relativity,
    )


@njit(**njit_dict_no_parallel)
def line_scatter_event(
    r_packet,
    time_explosion,
    line_interaction_type,
    opacity_state,
    enable_full_relativity,
):
    """
    Line scatter function that handles the scattering itself, including new angle drawn, and calculating nu out using macro atom

    Parameters
    ----------
    r_packet : tardis.transport.montecarlo.r_packet.RPacket
    time_explosion : float
    line_interaction_type : enum
    opacity_state : tardis.transport.montecarlo.numba_interface.OpacityState
    """
    old_doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion, enable_full_relativity
    )
    r_packet.mu = get_random_mu()

    inverse_new_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion, enable_full_relativity
    )

    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.energy = comov_energy * inverse_new_doppler_factor

    if line_interaction_type == LineInteractionType.SCATTER:
        line_emission(
            r_packet,
            r_packet.next_line_id,
            time_explosion,
            opacity_state,
            enable_full_relativity,
        )
    else:  # includes both macro atom and downbranch - encoded in the transition probabilities
        comov_nu = r_packet.nu * old_doppler_factor  # Is this necessary?
        r_packet.nu = comov_nu * inverse_new_doppler_factor
        activation_level_id = opacity_state.line2macro_level_upper[
            r_packet.next_line_id
        ]
        macro_atom_event(
            activation_level_id,
            r_packet,
            time_explosion,
            opacity_state,
            enable_full_relativity,
        )
