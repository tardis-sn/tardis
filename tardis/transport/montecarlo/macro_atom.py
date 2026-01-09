from enum import IntEnum

import numpy as np
from numba import njit

from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.configuration import montecarlo_globals


class MacroAtomError(ValueError):
    pass


class MacroAtomTransitionType(IntEnum):
    PHOTOIONIZATION = 3
    RECOMB_INTERNAL = 2
    INTERNAL_UP = 1
    INTERNAL_DOWN = 0
    BB_EMISSION = -1
    BF_EMISSION = -2  # This is recombination emission
    FF_EMISSION = -3
    ADIABATIC_COOLING = -4
    BF_COOLING = -5  # TODO: Maybe merge this with BF_EMISSION - Don't do that - need separate for cooling estimators
    TWO_PHOTON = -6
    COLL_DOWN_TO_K_PACKET = 9
    COLL_DOWN_INTERNAL = 10
    COLL_UP_INTERNAL = 11
    COLL_UP_COOLING = 12  # I think this should be a deactivation event


@njit(**njit_dict_no_parallel)
def macro_atom_interaction(
    activation_level_id, current_shell_id, opacity_state
):
    """
    Parameters
    ----------
    activation_level_id : int
        Activation level idx of the macro atom.
    current_shell_id : int
    opacity_state : tardis.transport.montecarlo.numba_interface.opacity_state.OpacityState

    Returns
    -------
    """
    current_transition_type = 0
    coll_down_to_k_packet_count_buffer = 0
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
                # Need to count the number of coll_deexc_down_to_k_packet interactions for the heating estimator
                # For the markov-chain macro atom, this won't work because this is an internal transition
                if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
                    if (
                        current_transition_type
                        == MacroAtomTransitionType.COLL_DOWN_TO_K_PACKET
                    ):
                        # This has to be multiplied by comoving energy, and then increment
                        # coll_deexc_heating_estimator for the shell
                        coll_down_to_k_packet_count_buffer += 1
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
        coll_down_to_k_packet_count_buffer,
    )
