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
