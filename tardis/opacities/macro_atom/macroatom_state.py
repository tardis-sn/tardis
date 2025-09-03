from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.transport.montecarlo.configuration import montecarlo_globals
import pandas as pd


class LegacyMacroAtomState(HDFWriterMixin):
    hdf_name = "macro_atom_state"

    hdf_properties = [
        "transition_probabilities",
        "transition_type",
        "destination_level_id",
        "transition_line_id",
        "macro_block_references",
        "line2macro_level_upper",
    ]

    def __init__(
        self,
        transition_probabilities,
        transition_type,
        destination_level_id,
        transition_line_id,
        macro_block_references,
        line2macro_level_upper,
    ):
        """
        Current State of the MacroAtom

        Parameters
        ----------
        transition_probabilities : pd.DataFrame
            Macro Atom Transition probabilities between levels
        transition_type : pd.DataFrame)
            Macro Atom transition types
        destination_level_id : pd.DataFrame
            ID of destination levels of the Macro Atom
        transition_line_id : pd.DataFrame
            ID of lines corresponding to Macro Atom transitions
        macro_block_references : pd.DataFrame or np.ndarray
            Index references to the Macro Atom blocks
        line2macro_level_upper : pd.DataFrame
            Mapping from lines to Macro Atom upper levels
        """
        self.transition_probabilities = transition_probabilities
        self.transition_type = transition_type
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id  # THESE ARE NOT TRANSITION LINE IDS. In the macro atom they were "lines_idx" and are the index locations of the lines in the lines dataframe.
        self.macro_block_references = macro_block_references
        self.line2macro_level_upper = line2macro_level_upper

    @classmethod
    def from_legacy_plasma(cls, plasma):
        """
        Generates a MacroAtomState object from a tardis BasePlasma

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            legacy base plasma

        Returns
        -------
        MacroAtomState
        """
        transition_probabilities = plasma.transition_probabilities
        transition_type = plasma.macro_atom_data["transition_type"]
        destination_level_id = plasma.macro_atom_data["destination_level_idx"]
        transition_line_id = plasma.macro_atom_data["lines_idx"]
        line2macro_level_upper = (
            plasma.atomic_data.lines_upper2macro_reference_idx
        )

        if (
            montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        ):  # TODO: Unify this in the plasma solver
            macro_block_references = plasma.macro_block_references
        else:
            macro_block_references = plasma.atomic_data.macro_atom_references[
                "block_references"
            ]

        return cls(
            transition_probabilities,
            transition_type,
            destination_level_id,
            transition_line_id,
            macro_block_references,
            line2macro_level_upper,
        )


class MacroAtomState:
    hdf_name = "macro_atom_state"

    hdf_properties = [
        "transition_probabilities",
        "transition_metadata",
        "line2macro_level_upper",
    ]

    def __init__(
        self,
        transition_probabilities: pd.DataFrame,
        transition_metadata: pd.DataFrame,
        line2macro_level_upper: pd.Series,
    ):
        """
        Current State of the MacroAtom

        Parameters
        ----------
        transition_probabilities : pd.DataFrame
            Transition probabilities for the macro atom, indexed by source and destination levels.
        transition_metadata : pd.DataFrame
            Metadata for the macro atom, including atomic number, ion number, level numbers for the transition, desination, and source.
        line2macro_level_upper : pd.Series
            Mapping from lines to the upper levels of the macro atom transitions.
        """
        self.transition_probabilities = transition_probabilities
        self.transition_metadata = transition_metadata
        self.line2macro_level_upper = line2macro_level_upper

    def to_legacy_format(self):
        """
        Convert the current state of the MacroAtom to legacy format.
        Returns
        -------


        """
        transition_probabilities = self.transition_probabilities
        transition_type = self.transition_metadata.transition_type
        destination_level_id = self.transition_metadata.destination_level_idx
        transition_line_id = self.transition_metadata.transition_line_idx
        unique_source_multi_index = pd.MultiIndex.from_tuples(
            self.transition_metadata.source.unique(),
            names=["atomic_number", "ion_number", "level_number"],
        )
        macro_data = (
            self.transition_metadata.reset_index()
            .groupby("source")
            .apply(lambda x: x.index[0])
        )
        macro_block_references = pd.Series(
            data=macro_data.values,
            index=unique_source_multi_index,
            name="macro_block_references",
        )
        line2macro_level_upper = self.line2macro_level_upper.values

        return LegacyMacroAtomState(
            transition_probabilities,
            transition_type,
            destination_level_id,
            transition_line_id,
            macro_block_references,
            line2macro_level_upper,
        )

    def sort_to_legacy(
        self, legacy_state: LegacyMacroAtomState, lines: pd.DataFrame
    ) -> "MacroAtomState":
        """
        Sort the current MacroAtomState to match the legacy MacroAtomState.

        Parameters
        ----------
        legacy_state : LegacyMacroAtomState
            The legacy state to sort to.
        lines : pd.DataFrame
            DataFrame containing line information.

        Returns
        -------
        MacroAtomState
            A new MacroAtomState sorted to match the legacy state.
        """

        legacy_sorting_frame = pd.DataFrame(legacy_state.transition_type)
        legacy_sorting_frame["transition_line_id"] = lines.iloc[
            legacy_state.transition_line_id
        ].line_id.values
        legacy_sorting_frame["match_key"] = legacy_sorting_frame.apply(
            lambda x: (x["transition_line_id"], x["transition_type"]), axis=1
        )  # match key uniquely identifies the transition

        resorting_frame = self.transition_metadata.copy()
        resorting_frame["match_key"] = resorting_frame.apply(
            lambda x: (x["transition_line_id"], x["transition_type"]), axis=1
        )

        resorted = (
            resorting_frame.reset_index()
            .set_index("match_key")
            .loc[legacy_sorting_frame["match_key"]]
            .set_index("macro_atom_transition_id")
        )

        resorted_transition_probabilities = self.transition_probabilities.loc[
            resorted.index
        ]
        resorted_metadata = self.transition_metadata.loc[resorted.index]

        return MacroAtomState(
            transition_probabilities=resorted_transition_probabilities,
            transition_metadata=resorted_metadata,
            line2macro_level_upper=self.line2macro_level_upper,
        )

    def recreate_legacy_macro_atom_state(
        self, legacy_state: LegacyMacroAtomState, lines: pd.DataFrame
    ) -> LegacyMacroAtomState:
        """
        Recreate the legacy MacroAtomState with new transition probabilities and new unique transition ids.

        Parameters
        ----------
        lines : pd.DataFrame
            DataFrame containing line information.

        Returns
        -------
        LegacyMacroAtomState
            The recreated legacy MacroAtomState.
        """
        legacy_sorted = self.sort_to_legacy(legacy_state, lines)
        legacy_macro_atom = legacy_sorted.to_legacy_format()
        legacy_macro_atom.macro_block_references = legacy_state.macro_block_references  # The old block references contained empty blocks. I can't recreate them from the new state so we just copy them over.
        legacy_macro_atom.line2macro_level_upper = (
            legacy_state.line2macro_level_upper
        )
        legacy_macro_atom.destination_level_id = (
            legacy_state.destination_level_id
        )

        return legacy_macro_atom
