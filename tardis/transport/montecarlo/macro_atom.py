from enum import IntEnum

import numpy as np
from numba import njit

from tardis.opacities.continuum.continuum_state_numba import (
    ContinuumOpacityStateNumba,
)
from tardis.opacities.opacity_state_numba import (
    OpacityStateNumba,
)
from tardis.transport.montecarlo import njit_dict_no_parallel


class MacroAtomError(ValueError):
    """Report invalid macro-atom probability or activation state."""


class MacroAtomTransitionType(IntEnum):
    """Identify internal jumps and negative-code deactivation processes."""

    # Bound bound transition block
    INTERNAL_UP = 1
    INTERNAL_DOWN = 0
    BB_EMISSION = -1

    # Currently unused transitions - bf and ff emission are implemented below by channel
    BF_EMISSION = -2
    FF_EMISSION = -3
    ADIABATIC_COOLING = -4
    BF_COOLING = -5  # TODO: Maybe merge this with BF_EMISSION - Yes this is taken care of by BF_EMISSION
    TWO_PHOTON = -6

    # Photoionization block
    PHOTOIONIZATION_INTERNAL = 2  # macro to i-packet
    PHOTO_RECOMB_INTERNAL = 3  # i-packet to macro
    PHOTO_RECOMB_EMISSION = -7  # i-packet to BF emission
    # Collisions block
    COLL_DOWN_TO_K_PACKET = 11  # macro to k-packet
    COLL_UP_INTERNAL = 4
    COLL_DOWN_INTERNAL = 5
    COLL_ION_INTERNAL = 6  # Macro to i-packet

    COLL_RECOMB_INTERNAL = 7  # i-packet to macro
    COLL_RECOMB_TO_K_PACKET = 14  # i-packet to k_packet creation
    # Cooling block
    FB_COOLING = -20  # k to bf emission
    FF_COOLING = -21  # k to ff emission
    COLL_EXC_COOL = 21  # k to macro
    COLL_ION_COOL = 22  # k to i-packet


@njit(**njit_dict_no_parallel)
def macro_atom_interaction(
    activation_level_id: int,
    current_shell_id: int,
    opacity_state: OpacityStateNumba,
) -> tuple[int, int]:
    """Sample original macro-atom transitions until a deactivation occurs.

    Parameters
    ----------
    activation_level_id
        Activation level idx of the macro atom.
    current_shell_id
    opacity_state

    Returns
    -------
    tuple[int, int]
        Emission channel identifier and negative transition type.
    """
    current_transition_type = 0
    while current_transition_type >= 0:
        probability = 0.0
        probability_event = np.random.random()

        block_start = opacity_state.macro_block_edge_index[activation_level_id]
        block_end = opacity_state.macro_block_edge_index[
            activation_level_id + 1
        ]

        # looping through the transition probabilities
        for transition_id in range(block_start, block_end):
            transition_probability = opacity_state.transition_probabilities[
                transition_id, current_shell_id
            ]

            probability += transition_probability

            if probability > probability_event:
                activation_level_id = opacity_state.destination_level_id[
                    transition_id
                ]
                current_transition_type = opacity_state.transition_type[
                    transition_id
                ]
                break
        else:
            raise MacroAtomError(
                "MacroAtom ran out of the block. This should not happen as "
                "the sum of probabilities is normalized to 1 and "
                "the probability_event should be less than 1"
            )

    return (
        opacity_state.transition_line_id[transition_id],
        current_transition_type,
    )


@njit(**njit_dict_no_parallel)
def macro_atom_interaction_iip(
    activation_level_idx: int,
    current_shell_id: int,
    opacity_state: OpacityStateNumba,
    continuum_state: ContinuumOpacityStateNumba,
) -> tuple[int, int]:
    """Alternate independent reduced visits and explicit retained transitions.

    Parameters
    ----------
    activation_level_idx : int
        Global activation level, including valid level zero and i/k pools.
    current_shell_id : int
        Active shell index.
    opacity_state : OpacityStateNumba
        Residual probabilities and original channel metadata.
    continuum_state : ContinuumOpacityStateNumba
        First-exit matrices and forward/inverse level maps.

    Returns
    -------
    emission_line_id : int
        Line or continuum ID for emitting process
    emission_process : int,
        Type of process emission defined by MacroAtomTransitionType in this file.
    """
    while True:
        if activation_level_idx < 0 or activation_level_idx >= len(
            continuum_state.macro_level_to_reduced_idx
        ):
            raise MacroAtomError(
                "MacroAtom activation lies outside the global level map"
            )
        reduced_idx = continuum_state.macro_level_to_reduced_idx[
            activation_level_idx
        ]
        if reduced_idx >= 0:
            probability_event = np.random.random()
            cumulative_probability = 0.0
            selected_idx = -1
            last_positive_idx = -1
            for exit_idx, state_probability in enumerate(
                continuum_state.absorbing_markov_probabilities[
                    current_shell_id, reduced_idx
                ]
            ):
                if state_probability > 0:
                    last_positive_idx = exit_idx
                cumulative_probability += state_probability
                if cumulative_probability > probability_event:
                    selected_idx = exit_idx
                    break
            if selected_idx < 0:
                if (
                    last_positive_idx >= 0
                    and abs(cumulative_probability - 1.0) <= 1e-10
                ):
                    selected_idx = last_positive_idx
                else:
                    raise MacroAtomError(
                        "MacroAtom first-exit probabilities are not normalized or activation has no outgoing probability"
                    )
            activation_level_idx = continuum_state.reduced_to_macro_level_idx[
                selected_idx
            ]

        block_start_idx = opacity_state.macro_block_edge_index[
            activation_level_idx
        ]
        block_end_idx = opacity_state.macro_block_edge_index[
            activation_level_idx + 1
        ]
        probability_event = np.random.random()
        cumulative_probability = 0.0
        selected_transition_idx = -1
        last_positive_idx = -1
        for transition_idx in range(block_start_idx, block_end_idx):
            transition_probability = opacity_state.transition_probabilities[
                transition_idx, current_shell_id
            ]
            if transition_probability > 0:
                last_positive_idx = transition_idx
            cumulative_probability += transition_probability
            if cumulative_probability > probability_event:
                selected_transition_idx = transition_idx
                break
        if selected_transition_idx < 0:
            if (
                last_positive_idx >= 0
                and abs(cumulative_probability - 1.0) <= 1e-10
            ):
                selected_transition_idx = last_positive_idx
            else:
                raise MacroAtomError(
                    "MacroAtom residual block is not normalized or activation has no outgoing probability"
                )
        transition_type = opacity_state.transition_type[selected_transition_idx]
        if transition_type < 0:
            return opacity_state.transition_line_id[
                selected_transition_idx
            ], transition_type
        activation_level_idx = opacity_state.destination_level_id[
            selected_transition_idx
        ]
