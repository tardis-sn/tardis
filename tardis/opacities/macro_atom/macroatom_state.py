from tardis.io.util import HDFWriterMixin
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
        self.transition_line_id = transition_line_id  # THESE ARE NOT TRANSITION LINE IDS, THEY ARE INTEGER LOCATIONS OF THE LINES. IN THE MACRO ATOM DATA THEY WERE LINES_IDX
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
        transition_probabilities,
        transition_metadata,
        line2macro_level_upper,
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

    def to_legacy(self):
        """
        Convert the current state of the MacroAtom to legacy format.
        Returns
        -------


        """
        transition_probabilities = self.transition_probabilities
        transition_type = self.transition_metadata.transition_type
        destination_level_id = pd.Series(
            data=[level[2] for level in self.transition_metadata.destination],
            index=self.transition_metadata.index,
            name="destination_level_idx",
        )
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
