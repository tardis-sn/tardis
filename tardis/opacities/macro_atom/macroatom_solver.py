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
from tardis.opacities.macro_atom.base import (
    calculate_non_markov_transition_probabilities,
    initialize_non_markov_transition_probabilities,
)
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState


class MacroAtomSolver:  # Possibly make two classes, one for normal and one for continuum
    def __init__(self, initialize=True, normalize=True):

        self.initialize = initialize
        self.normalize = normalize

    def initialize_non_markov_transition_probabilities(self, atomic_data):

        coef_and_block_ref = initialize_non_markov_transition_probabilities(
            atomic_data
        )
        self.transition_probability_coef = coef_and_block_ref[
            "transition_probability_coef"
        ]
        self.block_references = coef_and_block_ref["block_references"]
        self.initialize = False

    def solve_non_markov_transition_probabilities(
        self,
        atomic_data,
        legacy_plasma,
        tau_sobolev,
        stimulated_emission_factor,
    ):
        if self.initialize:
            self.initialize_non_markov_transition_probabilities(atomic_data)

        non_markov_transition_probabilities = (
            calculate_non_markov_transition_probabilities(
                atomic_data,
                legacy_plasma.beta_sobolev,
                legacy_plasma.j_blues,
                stimulated_emission_factor,
                tau_sobolev,
                self.transition_probability_coef,
                self.block_references,
                normalize=self.normalize,
            )
        )

        return non_markov_transition_probabilities

    def solve(
        self,
        legacy_plasma,
        atomic_data,
        tau_sobolev,
        stimulated_emission_factor,
    ):

        # TODO: Figure out how to calculate p_combined, Check TransitionProbabilitiesProperty in assemble_plasma, properties/base.py
        # Make the combined transition probabilities something that is configurable in the class
        transition_probabilities = (
            self.solve_non_markov_transition_probabilities(
                atomic_data,
                legacy_plasma,
                tau_sobolev,
                stimulated_emission_factor,
            )
        )

        macro_block_references = atomic_data.macro_atom_references[
            "block_references"
        ]
        macro_atom_info = legacy_plasma.atomic_data.macro_atom_data

        return MacroAtomState(
            transition_probabilities,
            macro_atom_info["transition_type"],
            macro_atom_info["destination_level_idx"],
            macro_atom_info["lines_idx"],
            macro_block_references,
            legacy_plasma.atomic_data.lines_upper2macro_reference_idx,
        )


class MacroAtomContinuumSolver(MacroAtomSolver):
    def solve(
        self,
        legacy_plasma,
        atomic_data,
        tau_sobolev,
        stimulated_emission_factor,
        continuum_interaction_species=None,
    ):

        if continuum_interaction_species is None:  # TODO: Fix this
            continuum_interaction_species = (
                legacy_plasma.continuum_interaction_species
            )

        non_markov_transition_probabilities = (
            self.solve_non_markov_transition_probabilities(
                atomic_data,
                legacy_plasma,
                tau_sobolev,
                stimulated_emission_factor,
            )
        )

        level_idxs2transition_idx = legacy_plasma.level_idxs2transition_idx
        cool_rate_fb = legacy_plasma.cool_rate_fb
        cool_rate_fb_tot = legacy_plasma.cool_rate_fb_tot
        level2continuum_idx = legacy_plasma.level2continuum_idx

        p_combined_args = [
            getattr(legacy_plasma, item)
            for item in (
                "p_rad_bb",
                "p_recomb",
                "p_coll",
                "p_two_photon",
                "cool_rate_adiabatic",
                "cool_rate_ff",
                "cool_rate_fb_tot",
                "p_coll_ion",
                "p_coll_recomb",
                "cool_rate_coll_ion",
            )
            if hasattr(legacy_plasma, item)
        ]  # Maybe do this in the init
        p_combined = calculate_p_combined(
            non_markov_transition_probabilities, *p_combined_args
        )

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
            B,
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
        macro_block_references = calculate_macro_block_references(
            combined_transition_probs
        )
        macro_atom_info = calculate_macro_atom_info(combined_transition_probs)
        transition_probabilities = combined_transition_probs

        return MacroAtomState(
            transition_probabilities,
            macro_atom_info["transition_type"],
            macro_atom_info["destination_level_idx"],
            macro_atom_info["lines_idx"],
            macro_block_references,
            legacy_plasma.atomic_data.lines_upper2macro_reference_idx,
        )
