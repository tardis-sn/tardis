from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.io.util import HDFWriterMixin
from tardis.opacities.macro_atom.transition_probabilities import (
    calculate_macro_block_references,
    calculate_combined_transition_probabilities,
    calculate_level_absorption_probs,
    calculate_deactivation_channel_probs,
    calculate_non_continuum_transitions_probs,
    calculate_markov_chain_transition_probs,
    calculate_fb_cooling_probs,
    calculate_non_continuum_trans_probs_mask,
    calculate_markov_chain_index,
    calculate_p_combined,
    calculate_macro_atom_info,
)
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState


class MacroAtomSolver:
    def __init__(self, *args, **kwargs):

        pass

    def solve(self, legacy_plasma, atomic_data, continuum_interaction_species):

        # TODO: Figure out how to calculate p_combined, Check TransitionProbabilitiesProperty in assemble_plasma, properties/base.py
        # Make the combined transition probabilities something that is configurable in the class
        non_markov_transition_probabilities = (
            legacy_plasma.non_markov_transition_probabilities
        )
        level_idxs2transition_idx = legacy_plasma.level_idxs2transition_idx
        cool_rate_fb = legacy_plasma.cool_rate_fb
        cool_rate_fb_tot = legacy_plasma.cool_rate_fb_tot
        level2continuum_idx = legacy_plasma.level2continuum_idx

        p_combined_args = [
            getattr(legacy_plasma, item)
            for item in ("p_rad_bb",)
            if hasattr(legacy_plasma, item)
        ]  # Maybe do this in the init
        p_combined_args = (legacy_plasma.p_rad_bb,)  # Do this for now
        p_combined = calculate_p_combined(p_combined_args)

        markov_chain_indices = calculate_markov_chain_index(
            atomic_data, continuum_interaction_species
        )
        idx2deactivation_idx = markov_chain_indices["idx2deactivation_idx"]
        k_packet_idx = markov_chain_indices["k_packet_idx"]
        idx2mkv_idx = markov_chain_indices["idx2mkv_idx"]

        markov_chain_transition_probs = calculate_markov_chain_transition_probs(
            p_combined, idx2mkv_idx
        )
        p_deactivation = markov_chain_transition_probs["p_deactivation"]
        B = markov_chain_transition_probs["B"]

        level_absorption_probs = calculate_level_absorption_probs(
            B, k_packet_idx, idx2deactivation_idx
        )
        fb_cooling_probs = calculate_fb_cooling_probs(
            cool_rate_fb,
            cool_rate_fb_tot,
            p_deactivation,
            level2continuum_idx,
            level_idxs2transition_idx,
        )
        deactivation_channel_probs = calculate_deactivation_channel_probs(
            level_idxs2transition_idx,
            p_deactivation,
            fb_cooling_probs,
            idx2deactivation_idx,
        )
        non_continuum_trans_probs_mask = (
            calculate_non_continuum_trans_probs_mask(
                atomic_data, continuum_interaction_species
            )
        )
        non_continuum_trans_probs = calculate_non_continuum_transitions_probs(
            atomic_data,
            non_markov_transition_probabilities,
            non_continuum_trans_probs_mask,
        )
        combined_transition_probs = calculate_combined_transition_probabilities(
            level_absorption_probs,
            deactivation_channel_probs,
            non_continuum_trans_probs,
        )
        if (
            montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        ):  # TODO: Unify this in the plasma solver

            macro_block_references = calculate_macro_block_references(
                combined_transition_probs
            )
        else:
            macro_block_references = atomic_data.macro_atom_references[
                "block_references"
            ]

        macro_atom_info = calculate_macro_atom_info(combined_transition_probs)
        transition_type = macro_atom_info["transition_type"]
        destination_level_id = macro_atom_info["destination_level_id"]
        transition_line_id = macro_atom_info["lines_idx"]

        return MacroAtomState(
            combined_transition_probs,
            transition_type,
            destination_level_id,
            macro_block_references,
            legacy_plasma.atomic_data.lines_upper2macro_reference_idx,
        )
