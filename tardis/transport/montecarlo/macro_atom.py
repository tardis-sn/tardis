import numpy as np
from enum import IntEnum

from numba import njit
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel


class MacroAtomError(ValueError):
    pass


class MacroAtomTransitionType(IntEnum):
    INTERNAL_UP = 1
    INTERNAL_DOWN = 0
    BB_EMISSION = -1
    BF_EMISSION = -2
    FF_EMISSION = -3
    ADIABATIC_COOLING = -4
    BF_COOLING = -5  # TODO: Maybe merge this with BF_EMISSION
    TWO_PHOTON = -6


@njit(**njit_dict_no_parallel)
def macro_atom(activation_level_id, current_shell_id, opacity_state):
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

    # current_transition_type = MacroAtomTransitionType(current_transition_type)
    return (
        opacity_state.transition_line_id[transition_id],
        current_transition_type,
    )
