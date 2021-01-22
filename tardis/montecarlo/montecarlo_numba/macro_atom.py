import numpy as np
from enum import IntEnum

from numba import njit
from tardis.montecarlo.montecarlo_numba import njit_dict


class MacroAtomError(ValueError):
    pass


class MacroAtomTransitionType(IntEnum):
    INTERNAL_UP = 1
    INTERNAL_DOWN = 0
    BB_EMISSION = -1
    BF_EMISSION = -2
    FF_EMISSION = -3
    ADIABATIC_COOLING = -4


@njit(**njit_dict)
def macro_atom(r_packet, numba_plasma):
    """

    Parameters
    ----------
    r_packet: tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    numba_plasma: tardis.montecarlo.numba_interface.numba_plasma

    Returns
    -------

    """
    activation_level_id = numba_plasma.line2macro_level_upper[
        r_packet.next_line_id
    ]
    current_transition_type = 0

    while current_transition_type >= 0:
        probability = 0.0
        probability_event = np.random.random()

        block_start = numba_plasma.macro_block_references[activation_level_id]
        block_end = numba_plasma.macro_block_references[activation_level_id + 1]

        # looping through the transition probabilities
        for transition_id in range(block_start, block_end):

            transition_probability = numba_plasma.transition_probabilities[
                transition_id, r_packet.current_shell_id
            ]

            probability += transition_probability

            if probability > probability_event:
                activation_level_id = numba_plasma.destination_level_id[
                    transition_id
                ]
                current_transition_type = numba_plasma.transition_type[
                    transition_id
                ]
                break

        else:
            raise MacroAtomError(
                "MacroAtom ran out of the block. This should not happen as "
                "the sum of probabilities is normalized to 1 and "
                "the probability_event should be less than 1"
            )

    if current_transition_type == MacroAtomTransitionType.BB_EMISSION:
        return numba_plasma.transition_line_id[transition_id]
    else:
        raise MacroAtomError("MacroAtom currently only allows BB transitions")
