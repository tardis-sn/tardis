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
    r_packet: tardis.montecarlo.montecarlo_numba.rpacket.RPacket
    numba_plasma: tardis.montecarlo.numba_interface.numba_plasma

    Returns
    -------

    """
    activation_level_id = numba_plasma.line2macro_level_upper[
        r_packet.next_line_id]
    current_transition_type = 0

    while current_transition_type >= 0:
        probability = 0.0
        probability_event = np.random.random()

        block_start = numba_plasma.macro_block_references[activation_level_id]
        block_end = numba_plasma.macro_block_references[activation_level_id + 1]

        for transition_id in range(block_start, block_end):

            transition_probability = numba_plasma.transition_probabilities[
                transition_id, r_packet.current_shell_id]

            probability += transition_probability
            if probability > probability_event:
                current_transition_type = numba_plasma.transition_type[
                    transition_id]
                break

        else:
            raise MacroAtomError(
                'MacroAtom ran out of the block. This should not happen as the sum '
                'of probabilities is normalized to 1 and the probability_event '
                'should be less than 1')

    if current_transition_type == MacroAtomTransitionType.BB_EMISSION:
        return numba_plasma.transition_line_id[transition_id]
    else:
        raise MacroAtomError('MacroAtom currently only allows BB transitions')
"""
#void
#macro_atom (rpacket_t * packet, const storage_model_t * storage, rk_state *mt_state)
#{
  int emit = 0, i = 0, offset = -1;
  uint64_t activate_level = rpacket_get_macro_atom_activation_level (packet);
  while (emit >= 0)
    {
      double event_random = rk_double (mt_state);
      i = storage->macro_block_references[activate_level] - 1;
      double p = 0.0;
      offset = storage->transition_probabilities_nd *
                             rpacket_get_current_shell_id (packet);
      do
        {
          ++i;
          p += storage->transition_probabilities[offset + i];
        }
      while (p <= event_random);
      emit = storage->transition_type[i];
      activate_level = storage->destination_level_id[i];
    }
  switch (emit)
    {
      case BB_EMISSION:
        line_emission (packet, storage, storage->transition_line_id[i], mt_state);
        break;

      case BF_EMISSION:
        rpacket_set_current_continuum_id (packet, storage->transition_line_id[i]);
        storage->last_line_interaction_out_id[rpacket_get_id (packet)] =
          rpacket_get_current_continuum_id (packet);

        continuum_emission (packet, storage, mt_state, sample_nu_free_bound, 3);
        break;

      case FF_EMISSION:
        continuum_emission (packet, storage, mt_state, sample_nu_free_free, 4);
        break;

      case ADIABATIC_COOLING:
        storage->last_interaction_type[rpacket_get_id (packet)] = 5;
        rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
        break;

      default:
        fprintf (stderr, "This process for macro-atom deactivation should not exist! (emit = %d)\n", emit);
        exit(1);
    }
}
"""