from numba import njit

from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.interaction_events import (
    LineInteractionType,
    adiabatic_cooling,
    line_emission,
)
from tardis.transport.montecarlo.macro_atom import (
    MacroAtomTransitionType,
    macro_atom_interaction,
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
    transition_id, transition_type = macro_atom_interaction(
        destination_level_idx, r_packet.current_shell_id, opacity_state
    )
    if transition_type == MacroAtomTransitionType.BB_EMISSION:
        line_emission(
            r_packet,
            transition_id,
            time_explosion,
            opacity_state,
            enable_full_relativity,
        )
    elif transition_type == MacroAtomTransitionType.ADIABATIC_COOLING:
        adiabatic_cooling(r_packet)
    else:
        raise Exception(
            f"Interaction {transition_type} not known or implemented!"
        )


@njit(**njit_dict_no_parallel)
def line_scatter_event(
    r_packet,
    time_explosion,
    line_interaction_type,
    opacity_state,
    enable_full_relativity,
):
    """Handle a classic-mode line interaction."""
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
        )
