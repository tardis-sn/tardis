import numpy as np
import pandas as pd


class MacroAtomData:
    def __init__(
        self,
        transition_probabilities=None,
        block_references=None,
        lines_lower2macro_reference_idx=None,
        lines_upper2macro_reference_idx=None,
    ):
        self.transition_probability_data = transition_probabilities
        self.block_reference_data = block_references
        self.lines_lower2macro_reference_idx = lines_lower2macro_reference_idx
        self.lines_upper2macro_reference_idx = lines_upper2macro_reference_idx

    def active_data(
        self,
        selected_atomic_numbers,
        line_interaction_type,
        lines_index,
        lines_upper2level_idx,
        lines_lower2level_idx,
    ):
        new_transition_probabilities = self.transition_probability_data.loc[
            self.transition_probability_data["atomic_number"].isin(
                selected_atomic_numbers
            )
        ].copy()

        new_block_references = self.block_reference_data[
            self.block_reference_data.index.isin(
                selected_atomic_numbers, level="atomic_number"
            )
        ].copy()

        if line_interaction_type == "downbranch":
            new_transition_probabilities = new_transition_probabilities.loc[
                new_transition_probabilities["transition_type"] == -1
            ]
            new_block_references = new_block_references.loc[
                new_block_references["count_down"] > 0
            ]
            new_block_references.loc[:, "count_total"] = new_block_references[
                "count_down"
            ]
            new_block_references.loc[:, "block_references"] = np.hstack(
                (
                    0,
                    np.cumsum(new_block_references["count_down"].values[:-1]),
                )
            )

        elif line_interaction_type == "macroatom":
            new_block_references.loc[:, "block_references"] = np.hstack(
                (
                    0,
                    np.cumsum(new_block_references["count_total"].values[:-1]),
                )
            )

        new_block_references.loc[:, "references_idx"] = np.arange(
            len(new_block_references)
        )

        new_transition_probabilities.loc[:, "lines_idx"] = lines_index.loc[
            new_transition_probabilities["transition_line_id"]
        ].values

        lines_upper2macro_reference_idx = (
            new_block_references.loc[lines_upper2level_idx, "references_idx"]
            .astype(np.int64)
            .values
        )

        if line_interaction_type == "macroatom":
            lines_lower2macro_reference_idx = (
                new_block_references.loc[
                    lines_lower2level_idx, "references_idx"
                ]
                .astype(np.int64)
                .values
            )
            # Sets all
            tmp_macro_destination_level_idx = pd.MultiIndex.from_arrays(
                [
                    new_transition_probabilities["atomic_number"],
                    new_transition_probabilities["ion_number"],
                    new_transition_probabilities["destination_level_number"],
                ]
            )

            tmp_macro_source_level_idx = pd.MultiIndex.from_arrays(
                [
                    new_transition_probabilities["atomic_number"],
                    new_transition_probabilities["ion_number"],
                    new_transition_probabilities["source_level_number"],
                ]
            )

            new_transition_probabilities.loc[:, "destination_level_idx"] = (
                new_block_references.loc[
                    tmp_macro_destination_level_idx, "references_idx"
                ]
                .astype(np.int64)
                .values
            )

            new_transition_probabilities.loc[:, "source_level_idx"] = (
                new_block_references.loc[
                    tmp_macro_source_level_idx, "references_idx"
                ]
                .astype(np.int64)
                .values
            )

        elif line_interaction_type == "downbranch":
            # Sets all the destination levels to -1 to indicate that they
            # are not used in downbranch calculations
            new_transition_probabilities.loc[:, "destination_level_idx"] = -1

        return MacroAtomData(
            new_transition_probabilities,
            new_block_references,
            lines_lower2macro_reference_idx,
            lines_upper2macro_reference_idx,
        )
