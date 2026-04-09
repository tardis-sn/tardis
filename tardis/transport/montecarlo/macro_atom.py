from enum import IntEnum

import numpy as np
from numba import njit

from tardis.opacities.opacity_state_numba import (
    OpacityStateNumba,
)
from tardis.opacities.opacity_state_numba_iip import (
    OpacityStateNumbaIIP,
)
from tardis.transport.montecarlo import njit_dict_no_parallel


class MacroAtomError(ValueError):
    pass


class MacroAtomTransitionType(IntEnum):
    # Bound bound transition block
    INTERNAL_UP = 1
    INTERNAL_DOWN = 0
    BB_EMISSION = -1

    # Currently unused transitions - bf and ff emission are implemented below
    BF_EMISSION = -2
    FF_EMISSION = -3
    ADIABATIC_COOLING = -4
    BF_COOLING = -5  # TODO: Maybe merge this with BF_EMISSION - Yes this is taken care of by BF_EMISSION
    TWO_PHOTON = -6

    # Photoionization block
    PHOTOIONIZATION_INTERNAL = 2  # macro to i-packet
    # PHOTOIONIZATION_TO_K_PACKET = (
    #     10  # macro to k, ctardis photoion deactivation
    # )
    PHOTO_RECOMB_INTERNAL = 3  # i-packet to macro
    PHOTO_RECOMB_EMISSION = -7  # i-packet to BF emission
    # Collisions block
    COLL_DOWN_TO_K_PACKET = 11  # macro to k-packet
    COLL_UP_INTERNAL = 4
    COLL_DOWN_INTERNAL = 5
    # COLL_EXC_TO_K_PACKET = 12  # k-packet to macro
    COLL_ION_INTERNAL = 6  # Macro to i-packet
    # COLL_ION_TO_K_PACKET = (
    #     13  # Macro to k-packet, ctardis coll ion deactivation
    # )
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
):
    """
    Parameters
    ----------
    activation_level_id
        Activation level idx of the macro atom.
    current_shell_id
    opacity_state

    Returns
    -------
    """
    current_transition_type = 0
    while current_transition_type >= 0:
        probability = 0.0
        probability_event = np.random.random()

        block_start = opacity_state.macro_block_references[activation_level_id]
        block_end = opacity_state.macro_block_references[
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
    opacity_state: OpacityStateNumbaIIP,
):
    """
    Parameters
    ----------
    activation_level_idx
        Activation level idx of the macro atom.
    current_shell_id
    opacity_state

    Returns
    -------
    emission_line_id : int
        Line or continuum ID for emitting process
    emission_process : int,
        Type of process emission defined by MacroAtomTransitionType in this file.
    """
    # step to absorbing state level
    current_transition_type = 0
    while current_transition_type >= 0:
        absorbing_state_probability = 0.0
        probability_event = np.random.random()

        if activation_level_idx == opacity_state.k_packet_idx:
            absorbing_activation_level_idx = activation_level_idx
        else:
            for to_state_index, state_probability in enumerate(
                opacity_state.absorbing_markov_probabilities[
                    current_shell_id, activation_level_idx
                ]
            ):
                absorbing_state_probability += state_probability

                if absorbing_state_probability > probability_event:
                    absorbing_activation_level_idx = to_state_index
                    break
            else:
                raise MacroAtomError(
                    "MacroAtom failed to select an absorbing state. The absorbing Markov"
                    "chain probabilities matrix may not be normalized or may contain only zeros. \n"
                    f"Attempted to activate from {activation_level_idx}."
                )

        # Handle second prob call for emission process from that state
        block_start_index = opacity_state.macro_block_references[
            absorbing_activation_level_idx
        ]
        block_end_index = opacity_state.macro_block_references[
            absorbing_activation_level_idx + 1
        ]
        emission_transition_probability = 0.0
        probability_emission_event = np.random.random()

        for deactivation_channel_index in range(
            block_start_index, block_end_index
        ):
            deactivation_probability = opacity_state.transition_probabilities[
                deactivation_channel_index, current_shell_id
            ]
            emission_transition_probability += deactivation_probability

            if emission_transition_probability > probability_emission_event:
                emission_process = opacity_state.transition_type[
                    deactivation_channel_index
                ]
                emission_line_id = opacity_state.transition_line_id[
                    deactivation_channel_index
                ]
                current_transition_type = opacity_state.transition_type[
                    deactivation_channel_index
                ]
                activation_level_idx = opacity_state.destination_level_id[
                    deactivation_channel_index
                ]  # Used only if you go back into the while loop
                break

        else:
            raise MacroAtomError(
                "MacroAtom ran out of the block. This should not happen as "
                "the sum of probabilities is normalized to 1 and "
                "the probability_event should be less than 1. \n"
                f"block indices are {block_start_index}, {block_end_index}"
            )
    return (
        emission_line_id,
        emission_process,
    )
