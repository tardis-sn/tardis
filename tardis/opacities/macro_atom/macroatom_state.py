from typing import Self

import pandas as pd

from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.transport.montecarlo.configuration import montecarlo_globals


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
        transition_probabilities: pd.DataFrame,
        transition_type: pd.Series,
        destination_level_id: pd.Series,
        transition_line_id: pd.Series,
        macro_block_references: pd.DataFrame | pd.Series,
        line2macro_level_upper,
    ) -> None:
        """
        Initialize a LegacyMacroAtomState object.

        Parameters
        ----------
        transition_probabilities : pd.DataFrame
            Macro Atom transition probabilities between levels.
        transition_type : pd.Series
            Macro Atom transition types.
        destination_level_id : pd.Series
            ID of destination levels of the Macro Atom.
        transition_line_id : pd.Series
            ID of lines corresponding to Macro Atom transitions.
        macro_block_references : pd.DataFrame | pd.Series
            Index references to the Macro Atom blocks.
        line2macro_level_upper
            Mapping from lines to Macro Atom upper levels.
        """
        self.transition_probabilities = transition_probabilities
        self.transition_type = transition_type
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id  # THESE ARE NOT TRANSITION LINE IDS. In the macro atom they were "lines_idx" and are the index locations of the lines in the lines dataframe.
        self.macro_block_references = macro_block_references
        self.line2macro_level_upper = line2macro_level_upper

    @classmethod
    def from_legacy_plasma(cls, plasma) -> Self:
        """
        Generate a LegacyMacroAtomState object from a tardis BasePlasma.

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            Legacy base plasma object.

        Returns
        -------
        LegacyMacroAtomState
            A LegacyMacroAtomState object created from the plasma data.
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
        macro_block_references: pd.Series | None = None,
    ) -> None:
        """
        Initialize a MacroAtomState object.

        Parameters
        ----------
        transition_probabilities : pd.DataFrame
            Transition probabilities for the macro atom, indexed by source and destination levels.
        transition_metadata : pd.DataFrame
            Metadata for the macro atom, including atomic number, ion number, level numbers for the transition, destination, and source.
        line2macro_level_upper : pd.Series
            Mapping from lines to the upper levels of the macro atom transitions.
        macro_block_references : pd.Series, optional
            Index references to the Macro Atom blocks. Default is None.
        """
        self.transition_probabilities = transition_probabilities
        self.transition_metadata = transition_metadata
        self.line2macro_level_upper = line2macro_level_upper
        self.macro_block_references = macro_block_references

    def to_legacy_format(self) -> LegacyMacroAtomState:
        """
        Convert the current state of the MacroAtom to legacy format.

        This method transforms the modern MacroAtomState structure into the
        legacy LegacyMacroAtomState format for backward compatibility. It extracts
        individual components from the consolidated transition_metadata DataFrame
        and converts them to the separate arrays/DataFrames expected by the
        legacy format.

        Returns
        -------
        LegacyMacroAtomState
            The MacroAtomState converted to legacy format.
        """
        transition_probabilities = self.transition_probabilities
        transition_type = self.transition_metadata.transition_type
        destination_level_id = self.transition_metadata.destination_level_idx
        transition_line_id = self.transition_metadata.transition_line_idx
        macro_block_references = self.macro_block_references
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
    ):
        """
        Sort the current MacroAtomState to match the legacy MacroAtomState order.

        Parameters
        ----------
        legacy_state : LegacyMacroAtomState
            The legacy state to match the ordering of.
        lines : pd.DataFrame
            DataFrame containing line information.

        Returns
        -------
        MacroAtomState
            A new MacroAtomState with transitions sorted to match the legacy state order.
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
            macro_block_references=self.macro_block_references,
        )

    def recreate_legacy_macro_atom_state(
        self, legacy_state: LegacyMacroAtomState, lines: pd.DataFrame
    ) -> LegacyMacroAtomState:
        """
        Recreate the legacy MacroAtomState with new transition probabilities and unique transition IDs.

        Parameters
        ----------
        legacy_state : LegacyMacroAtomState
            The original legacy state to recreate from.
        lines : pd.DataFrame
            DataFrame containing line information.

        Returns
        -------
        LegacyMacroAtomState
            The recreated legacy MacroAtomState with updated transition data.
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
