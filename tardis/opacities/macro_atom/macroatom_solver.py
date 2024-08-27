from tardis.opacities.macro_atom.base import (
    calculate_transition_probabilities,
    initialize_transition_probabilities,
)
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState


class MacroAtomSolver(object):

    initialize: bool = True
    normalize: bool = True

    def __init__(self, initialize=True, normalize=True):
        """Solver class for Macro Atom related opacities

        Parameters
        ----------
        initialize: bool
            Whether or not to initialize the transition probabilitiy coefficients and block references when solving the first time (default True)
        normalize: bool
            Whether or not to normalize the transition probabilities to unity. Default True
        """

        self.initialize = initialize
        self.normalize = normalize

    def initialize_transition_probabilities(self, atomic_data):
        """initialize the transition probability coefficients and block references when solving the first time

        Parameters
        ----------
        atomic_data : tardis.io.atom_data.AtomData
            Atomic Data
        """

        coef_and_block_ref = initialize_transition_probabilities(atomic_data)
        self.transition_probability_coef = coef_and_block_ref[
            "transition_probability_coef"
        ]
        self.block_references = coef_and_block_ref["block_references"]
        self.initialize = False

    def solve_transition_probabilities(
        self,
        atomic_data,
        legacy_plasma,
        tau_sobolev,
        stimulated_emission_factor,
    ):
        """Solve the basic transition probabilities for the macroatom

        Parameters
        ----------
        atomic_data : tardis.io.atom_data.AtomData
            Atomic Data
        legacy_plasma : tarids.plasma.BasePlasma
            legacy base plasma
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths
        stimulated_emission_factor : np.ndarray

        Returns
        -------
        pd.DataFrame
            Transition Probabilities
        """
        if self.initialize:
            self.initialize_transition_probabilities(atomic_data)

        transition_probabilities = calculate_transition_probabilities(
            atomic_data,
            legacy_plasma.beta_sobolev,
            legacy_plasma.j_blues,
            stimulated_emission_factor,
            tau_sobolev,
            self.transition_probability_coef,
            self.block_references,
            normalize=self.normalize,
        )

        return transition_probabilities

    def solve(
        self,
        legacy_plasma,
        atomic_data,
        tau_sobolev,
        stimulated_emission_factor,
    ):
        """Solved the Macro Atom State

        Parameters
        ----------
        legacy_plasma : tarids.plasma.BasePlasma
            legacy base plasma
        atomic_data : tardis.io.atom_data.AtomData
            Atomic Data
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths
        stimulated_emission_factor : pd.DataFrame

        Returns
        -------
        tardis.opacities.macroatom_state.MacroAtomState
            State of the macro atom ready to be placed into the OpacityState
        """

        transition_probabilities = self.solve_transition_probabilities(
            atomic_data,
            legacy_plasma,
            tau_sobolev,
            stimulated_emission_factor,
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
