from numba import njit

from tardis.transport.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.macro_atom import (
    MacroAtomTransitionType,
    macro_atom_interaction,
)
from tardis.transport.montecarlo.modes.nonhomologous.interaction_events import (
    LineInteractionType,
    line_emission,
)
from tardis.transport.montecarlo.utils import get_random_mu


@njit(**njit_dict_no_parallel)
def macro_atom_event(
    destination_level_idx,
    r_packet,
    geometry,
    opacity_state,
    enable_full_relativity,
):
    """
    Macroatom event handler - run the macroatom and handle the result

    Parameters
    ----------
    destination_level_idx : int
    r_packet : tardis.transport.montecarlo.r_packet.RPacket
    geometry : NumbaRadial1DGeometry
    opacity_state : tardis.transport.montecarlo.numba_interface.OpacityState
    """
    transition_id, transition_type = macro_atom_interaction(
        destination_level_idx, r_packet.current_shell_id, opacity_state
    )

    if transition_type == MacroAtomTransitionType.BB_EMISSION:
        line_emission(
            r_packet,
            transition_id,
            geometry,
            opacity_state,
            enable_full_relativity,
        )
    else:
        # A future continuum path needs a separate continuum opacity state
        # and an event handler adapted from modes.iip.interaction_event_callers
        # to this mode's radial geometry.
        raise Exception(
            f"Interaction {transition_type} not known or implemented!"
        )


@njit(**njit_dict_no_parallel)
def line_scatter_event(
    r_packet,
    geometry,
    line_interaction_type,
    opacity_state,
    enable_full_relativity,
):
    """
    Line scatter function that handles the scattering itself, including new angle drawn, and calculating nu out using macro atom

    Parameters
    ----------
    r_packet : tardis.transport.montecarlo.r_packet.RPacket
    geometry : NumbaRadial1DGeometry
    line_interaction_type : enum
    opacity_state : tardis.transport.montecarlo.numba_interface.OpacityState
    """
    v = geometry.get_velocity(r_packet.r, r_packet.current_shell_id)
    old_doppler_factor = get_doppler_factor(
        v, r_packet.mu, enable_full_relativity
    )
    r_packet.mu = get_random_mu()

    inverse_new_doppler_factor = get_inverse_doppler_factor(
        v, r_packet.mu, enable_full_relativity
    )

    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.energy = comov_energy * inverse_new_doppler_factor

    if line_interaction_type == LineInteractionType.SCATTER:
        line_emission(
            r_packet,
            r_packet.next_line_id,
            geometry,
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
            geometry,
            opacity_state,
            enable_full_relativity,
        )
