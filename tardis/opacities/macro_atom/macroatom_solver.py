import numpy as np
import pandas as pd

from tardis.opacities.macro_atom import util
from tardis.opacities.macro_atom.base import (
    get_macro_atom_data,
    initialize_transition_probabilities,
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


P_INTERNAL_UP = 1
P_INTERNAL_DOWN = 0
P_EMISSION_DOWN = -1


class MacroAtomSolver:
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
        """Solves the transition probabilities for the macro atom.
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
            A MacroAtomState object containing the transition probabilities and metadata.
        """

        # Solves the transition probabilities and returns a DataFrame with the probabilities and a DataFrame with the macro atom transition metadata.

        f_ul = self.lines.f_ul.values.reshape(-1, 1)
        f_lu = self.lines.f_lu.values.reshape(-1, 1)
        nus = self.lines.nu.values.reshape(-1, 1)

        # We rename the energies column so that it matches with the beta_sobolevs index.
        energies_upper = (
            self.levels[["energy"]]
            .rename(columns={"energy": "level_number_upper"})
            .reindex(beta_sobolevs.index.droplevel("level_number_lower"))
            .values
        )
        energies_lower = (
            self.levels[["energy"]]
            .rename(columns={"energy": "level_number_lower"})
            .reindex(beta_sobolevs.index.droplevel("level_number_upper"))
            .values
        )

        lines_level_upper = self.lines.index.droplevel("level_number_lower")

        p_emission_down = (
            2
            * nus**2
            * f_ul
            / const.c.cgs.value**2
            * beta_sobolevs
            * (energies_upper - energies_lower)
        )  # * gs_lower / gs_upper are considered by using f_ul instead of f_lu.
        p_internal_down = (
            2
            * nus**2
            * f_ul
            / const.c.cgs.value**2
            * beta_sobolevs
            * (energies_lower)
        )  # * gs_lower / gs_upper are considered by using f_ul instead of f_lu.
        p_internal_up = (
            f_lu
            / (const.h.cgs.value * nus)
            * stimulated_emission_factors
            * mean_intensities_blue_wing
            * beta_sobolevs
            * (energies_lower)
        )

        transition_a_i_l_u_array = (
            p_emission_down.reset_index()[
                [
                    "atomic_number",
                    "ion_number",
                    "level_number_lower",
                    "level_number_upper",
                ]
            ].values
        )  # This is a helper array to make the source and destination columns.

        p_internal_up["source"] = [
            tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 2]]
        ]
        p_internal_down["source"] = [
            tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 3]]
        ]
        p_emission_down["source"] = [
            tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 3]]
        ]

        probabilities_df = pd.concat(
            [p_internal_up, p_internal_down, p_emission_down]
        )
        probabilities_index = probabilities_df.index

        probabilities_df.index = probabilities_index

        probabilities_df = probabilities_df.reset_index(drop=True)

        if normalize:
            # Normalize the probabilities
            probabilities_df = probabilities_df.div(
                probabilities_df.groupby("source").transform("sum")
            )
            probabilities_df.replace(
                np.nan, 0, inplace=True
            )  # Some blocks have no transitions, so we replace NaN with 0.

        probabilities_df.drop(columns=["source"], inplace=True)
        probabilities_df = probabilities_df.reset_index(drop=True)

        p_internal_up["destination"] = [
            tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 3]]
        ]
        p_internal_down["destination"] = [
            tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 2]]
        ]
        p_emission_down["destination"] = [
            tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 2]]
        ]

        p_internal_up["transition_type"] = P_INTERNAL_UP
        p_internal_down["transition_type"] = P_INTERNAL_DOWN
        p_emission_down["transition_type"] = P_EMISSION_DOWN

        p_internal_up["transition_line_id"] = p_internal_down[
            "transition_line_id"
        ] = p_emission_down["transition_line_id"] = self.lines.line_id.values
        p_internal_up["transition_line_idx"] = p_internal_down[
            "transition_line_idx"
        ] = p_emission_down["transition_line_idx"] = range(len(self.lines))

        # Grab all the columns for the metadata.
        macro_atom_transition_metadata = pd.DataFrame(index=probabilities_index)
        for col in [
            "transition_line_id",
            "source",
            "destination",
            "transition_type",
            "transition_line_idx",
        ]:
            macro_atom_transition_metadata[col] = np.concat(
                [p_internal_up[col], p_internal_down[col], p_emission_down[col]]
            )

        macro_atom_transition_metadata = (
            macro_atom_transition_metadata.reset_index()
        )

        macro_atom_transition_metadata = (
            macro_atom_transition_metadata.sort_values(
                ["atomic_number", "ion_number", "source"]
            )
        )  # This is how carsus sorted the macro atom transitions.

        probabilities_df = probabilities_df.loc[
            macro_atom_transition_metadata.index
        ]

        probabilities_df.index.names = ["macro_atom_transition_id"]
        macro_atom_transition_metadata.index.names = [
            "macro_atom_transition_id"
        ]

        # We have to do this at the end so that it's sorted the same way.
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
            ],
            inplace=True,
        )

        return MacroAtomState(
            probabilities_df,
            macro_atom_transition_metadata,
            line2macro_level_upper,
        )
