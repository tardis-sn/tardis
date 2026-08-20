import numpy as np
from numba import njit

from tardis.opacities.continuum.continuum_state_numba import (
    ContinuumOpacityStateNumba,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.interaction_events import (
    LineInteractionType,
    adiabatic_cooling,
    bound_free_emission,
    determine_bf_macro_activation_idx,
    free_free_emission,
    line_emission,
)
from tardis.transport.montecarlo.macro_atom import (
    MacroAtomTransitionType,
    macro_atom_interaction_iip,
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket
from tardis.transport.montecarlo.utils import get_random_mu


@njit(**njit_dict_no_parallel)
def macro_atom_event(
    destination_level_idx: int,
    r_packet: RPacket,
    time_explosion: float,
    opacity_state: OpacityStateNumba,
    enable_full_relativity: bool,
    continuum_state: ContinuumOpacityStateNumba,
) -> None:
    """
    Macroatom event handler - run the macroatom and handle the result

    Parameters
    ----------
    destination_level_idx
    r_packet
    time_explosion
    opacity_state
    """
    transition_id, transition_type = macro_atom_interaction_iip(
        destination_level_idx,
        r_packet.current_shell_id,
        opacity_state,
        continuum_state,
    )

    if transition_type in (
        MacroAtomTransitionType.FF_EMISSION,
        MacroAtomTransitionType.FF_COOLING,
    ):
        free_free_emission(
            r_packet, time_explosion, opacity_state, enable_full_relativity
        )
    elif transition_type in (
        MacroAtomTransitionType.BF_EMISSION,
        MacroAtomTransitionType.FB_COOLING,
        MacroAtomTransitionType.PHOTO_RECOMB_EMISSION,
    ):
        bound_free_emission(
            r_packet,
            time_explosion,
            opacity_state,
            continuum_state,
            transition_id,
            enable_full_relativity,
        )
    elif transition_type == MacroAtomTransitionType.ADIABATIC_COOLING:
        adiabatic_cooling(r_packet)
    elif transition_type == MacroAtomTransitionType.BB_EMISSION:
        line_emission(
            r_packet,
            transition_id,
            time_explosion,
            opacity_state,
            enable_full_relativity,
        )
    else:
        raise Exception(f"Interaction {transition_type} not known or implemented!")


@njit(**njit_dict_no_parallel)
def determine_continuum_macro_activation_idx(
    continuum_state,
    nu,
    chi_bf,
    chi_ff,
    chi_bf_contributions,
    active_continua,
):
    """Determine the macro-atom activation level after continuum absorption."""
    fraction_bf = chi_bf / (chi_bf + chi_ff)
    if np.random.random() < fraction_bf:
        return determine_bf_macro_activation_idx(
            continuum_state, nu, chi_bf_contributions, active_continua
        )
    return continuum_state.k_packet_idx


@njit(**njit_dict_no_parallel)
def continuum_event(
    r_packet,
    time_explosion,
    opacity_state,
    continuum_state,
    chi_bf_tot,
    chi_ff,
    chi_bf_contributions,
    current_continua,
    enable_full_relativity,
):
    """Handle a continuum absorption event in IIP transport."""
    velocity = r_packet.r / time_explosion
    old_doppler_factor = get_doppler_factor(
        velocity, r_packet.mu, enable_full_relativity
    )
    r_packet.mu = get_random_mu()
    inverse_doppler_factor = get_inverse_doppler_factor(
        velocity, r_packet.mu, enable_full_relativity
    )
    r_packet.energy *= old_doppler_factor * inverse_doppler_factor
    comov_nu = r_packet.nu * old_doppler_factor

    destination_level_idx = determine_continuum_macro_activation_idx(
        continuum_state,
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
        continuum_state,
    )


@njit(**njit_dict_no_parallel)
def line_scatter_event(
    r_packet: RPacket,
    time_explosion: float,
    line_interaction_type: LineInteractionType,
    opacity_state: OpacityStateNumba,
    enable_full_relativity: bool,
    continuum_state: ContinuumOpacityStateNumba,
) -> None:
    """Handle a line interaction in IIP transport."""
    velocity = r_packet.r / time_explosion
    old_doppler_factor = get_doppler_factor(
        velocity, r_packet.mu, enable_full_relativity
    )
    r_packet.mu = get_random_mu()
    inverse_new_doppler_factor = get_inverse_doppler_factor(
        velocity, r_packet.mu, enable_full_relativity
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
    else:
        comov_nu = r_packet.nu * old_doppler_factor
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
            continuum_state,
        )
