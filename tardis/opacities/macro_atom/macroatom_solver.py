import numpy as np
import pandas as pd

from tardis.opacities.macro_atom import util
from tardis.opacities.macro_atom.base import (
    get_macro_atom_data,
    initialize_transition_probabilities,
)
from tardis.opacities.macro_atom.macroatom_transitions import (
    line_transition_emission_down,
    line_transition_internal_down,
    line_transition_internal_up,
)
from tardis.opacities.macro_atom.macroatom_state import (
    MacroAtomState,
    LegacyMacroAtomState,
)
from astropy import constants as const


class LegacyMacroAtomSolver:
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
        mean_intensities_lines_blue_wing,
        tau_sobolev,
        beta_sobolev,
        stimulated_emission_factor,
    ):
        """Solve the basic transition probabilities for the macroatom

        Parameters
        ----------
        atomic_data : tardis.io.atom_data.AtomData
            Atomic Data
        mean_intensities_lines_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell
            For more detail see Lucy 2003, https://doi.org/10.1051/0004-6361:20030357
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths
        beta_sobolev : pd.DataFrame
            Modified expansion Optical Depths
        stimulated_emission_factor : np.ndarray

        Returns
        -------
        pd.DataFrame
            Transition Probabilities
        """
        if self.initialize:
            self.initialize_transition_probabilities(atomic_data)
        # Referenced in https://github.com/tardis-sn/tardis/issues/3009
        if len(mean_intensities_lines_blue_wing) == 0:
            return None
        macro_atom_data = get_macro_atom_data(atomic_data)

        transition_probabilities = np.empty(
            (self.transition_probability_coef.shape[0], beta_sobolev.shape[1])
        )
        transition_type = macro_atom_data.transition_type.values
        lines_idx = macro_atom_data.lines_idx.values
        tpos = macro_atom_data.transition_probability.values
        # This function modifies transition_probabilities inplace
        util.fast_calculate_transition_probabilities(
            tpos,
            beta_sobolev.values,
            mean_intensities_lines_blue_wing.values,
            stimulated_emission_factor,
            transition_type,
            lines_idx,
            self.block_references,
            transition_probabilities,
            self.normalize,
        )
        transition_probabilities_df = pd.DataFrame(
            transition_probabilities,
            index=macro_atom_data.transition_line_id,
            columns=tau_sobolev.columns,
        )

        return transition_probabilities_df

    def solve(
        self,
        mean_intensities_lines_blue_wing,
        atomic_data,
        tau_sobolev,
        stimulated_emission_factor,
        beta_sobolev=None,
    ):
        """Solved the Macro Atom State

        Parameters
        ----------
        mean_intensities_lines_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell
        atomic_data : tardis.io.atom_data.AtomData
            Atomic Data
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths
        stimulated_emission_factor : pd.DataFrame
        beta_sobolev : pd.DataFrame

        Returns
        -------
        tardis.opacities.macroatom_state.MacroAtomState
            State of the macro atom ready to be placed into the OpacityState
        """
        transition_probabilities = self.solve_transition_probabilities(
            atomic_data,
            mean_intensities_lines_blue_wing,
            tau_sobolev,
            beta_sobolev,
            stimulated_emission_factor,
        )

        macro_block_references = atomic_data.macro_atom_references[
            "block_references"
        ]
        macro_atom_info = atomic_data.macro_atom_data

        return LegacyMacroAtomState(
            transition_probabilities,
            macro_atom_info["transition_type"],
            macro_atom_info["destination_level_idx"],
            macro_atom_info["lines_idx"],
            macro_block_references,
            atomic_data.lines_upper2macro_reference_idx,
        )


class BoundBoundMacroAtomSolver:
    levels: pd.DataFrame
    lines: pd.DataFrame

    def __init__(self, levels, lines):
        self.levels = levels
        self.lines = lines

    def solve(
        self,
        mean_intensities_blue_wing,
        beta_sobolevs,
        stimulated_emission_factors,
        normalize=True,
    ):
        """
        Solves the transition probabilities and returns a DataFrame with the probabilities and a DataFrame with the macro atom transition metadata.
        Referenced as $p_i$ in Lucy 2003, https://doi.org/10.1051/0004-6361:20030357
        Parameters
        ----------
        mean_intensities_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
            For more detail see Lucy 2003, https://doi.org/10.1051/0004-6361:20030357
            Referenced as 'J^b_{lu}' internally, or 'J^b_{ji}' in the original paper.
        beta_sobolevs : pd.DataFrame
            Escape probabilites for the Sobolev approximation.
        stimulated_emission_factors : pd.DataFrame
            Stimulated emission factors for the lines.
        normalize : bool, optional
            Whether to normalize the transition probabilities to unity. Default is True.
        Returns
        -------
        MacroAtomState
            A MacroAtomState object containing the transition probabilities, transition metadata, and a mapping from line IDs to macro atom level upper indices.
        """

        f_ul = self.lines.f_ul.values.reshape(-1, 1)
        f_lu = self.lines.f_lu.values.reshape(-1, 1)
        nus = self.lines.nu.values.reshape(-1, 1)
        line_ids = self.lines.line_id.values

        energies_upper = (
            self.levels[["energy"]]
            .rename(columns={"energy": "level_number_upper"})
            .reindex(self.lines.index.droplevel("level_number_lower"))
            .values
        )
        energies_lower = (
            self.levels[["energy"]]
            .rename(columns={"energy": "level_number_lower"})
            .reindex(self.lines.index.droplevel("level_number_upper"))
            .values
        )
        transition_a_i_l_u_array = (
            self.lines.reset_index()[
                [
                    "atomic_number",
                    "ion_number",
                    "level_number_lower",
                    "level_number_upper",
                ]
            ].values
        )  # This is a helper array to make the source and destination columns.

        lines_level_upper = self.lines.index.droplevel("level_number_lower")

        p_emission_down, emission_down_metadata = line_transition_emission_down(
            f_ul,
            nus,
            energies_upper,
            energies_lower,
            beta_sobolevs,
            transition_a_i_l_u_array,
            line_ids,
        )
        p_internal_down, internal_down_metadata = line_transition_internal_down(
            f_ul,
            nus,
            energies_lower,
            beta_sobolevs,
            transition_a_i_l_u_array,
            line_ids,
        )
        p_internal_up, internal_up_metadata = line_transition_internal_up(
            f_lu,
            nus,
            energies_lower,
            mean_intensities_blue_wing,
            beta_sobolevs,
            stimulated_emission_factors,
            transition_a_i_l_u_array,
            line_ids,
        )

        probabilities_df = pd.concat(
            [p_emission_down, p_internal_down, p_internal_up]
        )

        macro_atom_transition_metadata = pd.concat(
            [
                emission_down_metadata,
                internal_down_metadata,
                internal_up_metadata,
            ]
        )

        if normalize:
            # Normalize the probabilities by source.
            probabilities_df = probabilities_df.div(
                probabilities_df.groupby("source").transform("sum"),
                fill_value=0,
            )  # fill value for nans where the transition probabilites are all 0, which happens for bound levels that should never be accessed in the active macroatom.

        probabilities_df.drop(columns=["source"], inplace=True)
        probabilities_df = probabilities_df.reset_index(
            drop=True
        )  # Reset to create a unique macro_atom_transition_id.
        probabilities_df.index.rename("macro_atom_transition_id", inplace=True)

        macro_atom_transition_metadata = (
            macro_atom_transition_metadata.reset_index()
        )
        macro_atom_transition_metadata.index.rename(
            "macro_atom_transition_id", inplace=True
        )
        macro_atom_transition_metadata["source_level"] = (
            macro_atom_transition_metadata.source.apply(lambda x: x[2])
        )
        macro_atom_transition_metadata = (
            macro_atom_transition_metadata.sort_values(
                ["atomic_number", "ion_number", "source_level"]
            )
        )  # This is how carsus sorted the macro atom transitions.

        probabilities_df = probabilities_df.loc[
            macro_atom_transition_metadata.index
        ]  # Reorder to match the metadata, which was sorted to match carsus.

        # We have to create the line2macro object after sorting.
        unique_source_index = pd.MultiIndex.from_tuples(
            macro_atom_transition_metadata.source.unique(),
            names=["atomic_number", "ion_number", "level_number"],
        )
        unique_source_series = pd.Series(
            index=unique_source_index,
            data=range(len(macro_atom_transition_metadata.source.unique())),
        )
        line2macro_level_upper = unique_source_series.loc[lines_level_upper]

        macro_atom_transition_metadata.drop(
            columns=[
                "atomic_number",
                "ion_number",
                "level_number_lower",
                "level_number_upper",
                "source_level",
            ],
            inplace=True,
        )

        return MacroAtomState(
            probabilities_df,
            macro_atom_transition_metadata,
            line2macro_level_upper,
        )
