import numpy as np
import pandas as pd

from tardis.configuration.sorting_globals import SORTING_ALGORITHM
from tardis.io.atom_data import AtomData
from tardis.opacities.macro_atom import util
from tardis.opacities.macro_atom.base import (
    get_macro_atom_data,
    initialize_transition_probabilities,
)
from tardis.opacities.macro_atom.macroatom_state import (
    LegacyMacroAtomState,
    MacroAtomState,
)
from tardis.opacities.macro_atom.macroatom_transitions import (
    line_transition_emission_down,
    line_transition_internal_down,
    line_transition_internal_up,
)


class LegacyMacroAtomSolver:
    initialize: bool = True
    normalize: bool = True

    def __init__(self, initialize: bool = True, normalize: bool = True) -> None:
        """
        Initialize the LegacyMacroAtomSolver.

        Parameters
        ----------
        initialize : bool, optional
            Whether or not to initialize the transition probability coefficients and block references when solving the first time. Default is True.
        normalize : bool, optional
            Whether or not to normalize the transition probabilities to unity. Default is True.
        """
        self.initialize = initialize
        self.normalize = normalize

    def initialize_transition_probabilities(
        self, atomic_data: AtomData
    ) -> None:
        """
        Initialize the transition probability coefficients and block references.

        This method should be called when solving for the first time to set up
        the necessary coefficients and block references.

        Parameters
        ----------
        atomic_data : AtomData
            Atomic data containing the necessary information for initialization.
        """
        coef_and_block_ref = initialize_transition_probabilities(atomic_data)
        self.transition_probability_coef = coef_and_block_ref[
            "transition_probability_coef"
        ]
        self.block_references = coef_and_block_ref["block_references"]
        self.initialize = False

    def solve_transition_probabilities(
        self,
        atomic_data: AtomData,
        mean_intensities_lines_blue_wing: pd.DataFrame,
        tau_sobolev: pd.DataFrame,
        beta_sobolev: pd.DataFrame | None,
        stimulated_emission_factor: pd.DataFrame | np.ndarray,
    ) -> pd.DataFrame | None:
        """
        Solve the basic transition probabilities for the macroatom.

        Parameters
        ----------
        atomic_data : AtomData
            Atomic data containing macro atom information.
        mean_intensities_lines_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
            For more detail see Lucy 2003, https://doi.org/10.1051/0004-6361:20030357
        tau_sobolev : pd.DataFrame
            Expansion optical depths.
        beta_sobolev : pd.DataFrame | None
            Modified expansion optical depths.
        stimulated_emission_factor : pd.DataFrame | np.ndarray
            Stimulated emission factors.

        Returns
        -------
        pd.DataFrame | None
            Transition probabilities. Returns None if mean_intensities_lines_blue_wing is empty.
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
        mean_intensities_lines_blue_wing: pd.DataFrame,
        atomic_data: AtomData,
        tau_sobolev: pd.DataFrame,
        stimulated_emission_factor: pd.DataFrame,
        beta_sobolev: pd.DataFrame | None = None,
    ) -> LegacyMacroAtomState:
        """
        Solve the Macro Atom State.

        Parameters
        ----------
        mean_intensities_lines_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
        atomic_data : AtomData
            Atomic data containing macro atom information.
        tau_sobolev : pd.DataFrame
            Expansion optical depths.
        stimulated_emission_factor : pd.DataFrame
            Stimulated emission factors.
        beta_sobolev : pd.DataFrame | None, optional
            Modified expansion optical depths. Default is None.

        Returns
        -------
        LegacyMacroAtomState
            State of the macro atom ready to be placed into the OpacityState.
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

    def __init__(
        self,
        levels: pd.DataFrame,
        lines: pd.DataFrame,
        line_interaction_type: str = "macroatom",
    ) -> None:
        """
        Initialize the BoundBoundMacroAtomSolver.

        Parameters
        ----------
        levels : pd.DataFrame
            DataFrame containing atomic level information.
        lines : pd.DataFrame
            DataFrame containing spectral line information.
        line_interaction_type : str, optional
            Type of line interaction to use. Default is "macroatom".
        """
        self.levels = levels
        self.lines = lines
        self.line_interaction_type = line_interaction_type

    def solve(
        self,
        mean_intensities_blue_wing: pd.DataFrame,
        beta_sobolevs: pd.DataFrame,
        stimulated_emission_factors: pd.DataFrame,
    ) -> MacroAtomState:
        """
        Solve the transition probabilities for the macroatom.

        This method calculates transition probabilities and returns a MacroAtomState object
        with the probabilities and macro atom transition metadata.
        Referenced as $p_i$ in Lucy 2003, https://doi.org/10.1051/0004-6361:20030357

        Parameters
        ----------
        mean_intensities_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
            For more detail see Lucy 2003, https://doi.org/10.1051/0004-6361:20030357
            Referenced as 'J^b_{lu}' internally, or 'J^b_{ji}' in the original paper.
        beta_sobolevs : pd.DataFrame
            Escape probabilities for the Sobolev approximation.
        stimulated_emission_factors : pd.DataFrame
            Stimulated emission factors for the lines.

        Returns
        -------
        MacroAtomState
            A MacroAtomState object containing the transition probabilities, transition metadata,
            and a mapping from line IDs to macro atom level upper indices.
        """
        oscillator_strength_ul = self.lines.f_ul.values.reshape(-1, 1)
        oscillator_strength_lu = self.lines.f_lu.values.reshape(-1, 1)
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
        transition_a_i_l_u_array = self.lines.reset_index()[
            [
                "atomic_number",
                "ion_number",
                "level_number_lower",
                "level_number_upper",
            ]
        ].values  # This is a helper array to make the source and destination columns. The letters stand for atomic_number, ion_number, lower level, upper level.

        lines_level_upper = self.lines.index.droplevel("level_number_lower")

        if self.line_interaction_type in ["downbranch", "macroatom"]:
            p_emission_down, emission_down_metadata = (
                line_transition_emission_down(
                    oscillator_strength_ul,
                    nus,
                    energies_upper,
                    energies_lower,
                    beta_sobolevs,
                    transition_a_i_l_u_array,
                    line_ids,
                )
            )
        else:
            raise ValueError(
                f"Unknown line interaction type: {self.line_interaction_type}"
            )
        probabilities_df = p_emission_down
        macro_atom_transition_metadata = emission_down_metadata

        if self.line_interaction_type == "macroatom":
            p_internal_down, internal_down_metadata = (
                line_transition_internal_down(
                    oscillator_strength_ul,
                    nus,
                    energies_lower,
                    beta_sobolevs,
                    transition_a_i_l_u_array,
                    line_ids,
                )
            )
            p_internal_up, internal_up_metadata = line_transition_internal_up(
                oscillator_strength_lu,
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

        # Normalize the probabilities by source. This used to be optional but is never not done in TARDIS. This also removes the source column from the probabilities DataFrame.
        normalized_probabilities = normalize_transition_probabilities(
            probabilities_df
        )

        normalized_probabilities, macro_atom_transition_metadata = (
            reindex_sort_and_clean_probabilities_and_metadata(
                normalized_probabilities, macro_atom_transition_metadata
            )
        )

        # We have to create the line2macro object after sorting.
        line2macro_level_upper = create_line2macro_level_upper(
            macro_atom_transition_metadata, lines_level_upper
        )

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
        source_to_index = {
            source: idx
            for idx, source in enumerate(
                macro_atom_transition_metadata.source.unique()
            )
        }
        # -99 should never be used downstream. The presence of it means the destination is not a source
        # , which means that the destination is only referenced from emission
        # (or macroatom deactivation) for the given macroatom configuration.
        macro_atom_transition_metadata["destination_level_idx"] = (
            (macro_atom_transition_metadata.destination.map(source_to_index))
            .fillna(-99)
            .astype(np.int64)
        )

        return MacroAtomState(
            normalized_probabilities,
            macro_atom_transition_metadata,
            line2macro_level_upper,
        )


def create_line2macro_level_upper(
    macro_atom_transition_metadata: pd.DataFrame,
    lines_level_upper: pd.MultiIndex,
) -> pd.Series:
    """
    Create a mapping from line transitions to macro atom level indices for upper levels.
    This method creates a mapping that connects line transition upper levels to their
    corresponding macro atom level indices. It first extracts unique source levels
    from the macro atom transition metadata and assigns sequential indices to them,
    then maps the line upper levels to these indices.

    Parameters
    ----------
    macro_atom_transition_metadata : pd.DataFrame
        DataFrame containing macro atom transition metadata
    lines_level_upper : pd.MultiIndex
        MultiIndex containing line upper level information

    Returns
    -------
    pd.Series
        Series mapping line transitions to macro atom level indices
    """
    unique_source_index = pd.MultiIndex.from_tuples(
        macro_atom_transition_metadata.source.unique(),
        names=["atomic_number", "ion_number", "level_number"],
    )
    unique_source_series = pd.Series(
        index=unique_source_index,
        data=range(len(macro_atom_transition_metadata.source.unique())),
    )
    line2macro_level_upper = unique_source_series.loc[lines_level_upper]

    return line2macro_level_upper


def normalize_transition_probabilities(
    probabilities_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Normalize transition probabilities by their source levels.

    Parameters
    ----------
    probabilities_df : pd.DataFrame
        DataFrame containing transition probabilities with a 'source' column
        for grouping.

    Returns
    -------
    pd.DataFrame
        Normalized probabilities where each source group sums to 1.0.
        NaN values are replaced with 0.0 for cases where all transition
        probabilities are zero (typically ground levels in macroatom).
    """
    # Normalize the probabilities by source. This used to be optional but is never not done in TARDIS.
    normalized_probabilities = probabilities_df.div(
        probabilities_df.groupby("source").transform("sum"),
    )
    normalized_probabilities.replace(np.nan, 0, inplace=True)

    return normalized_probabilities.drop(columns=["source"])


def reindex_sort_and_clean_probabilities_and_metadata(
    normalized_probabilities: pd.DataFrame,
    macro_atom_transition_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reindex and sort macro atom transition probabilities and metadata.

    Parameters
    ----------
    normalized_probabilities : pd.DataFrame
        DataFrame containing normalized transition probabilities.
    macro_atom_transition_metadata : pd.DataFrame
        DataFrame containing metadata for macro atom transitions.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Reindexed normalized probabilities and cleaned metadata sorted by
        atomic number, ion number, and source level.
    """
    normalized_probabilities = normalized_probabilities.reset_index(
        drop=True
    )  # Reset to create a unique macro_atom_transition_id.
    normalized_probabilities.index.rename(
        "macro_atom_transition_id", inplace=True
    )

    macro_atom_transition_metadata = (
        macro_atom_transition_metadata.reset_index()
    )
    macro_atom_transition_metadata.index.rename(
        "macro_atom_transition_id", inplace=True
    )
    macro_atom_transition_metadata["source_level"] = (
        macro_atom_transition_metadata.source.apply(lambda x: x[2])
    )
    macro_atom_transition_metadata = macro_atom_transition_metadata.sort_values(
        [
            "atomic_number",
            "ion_number",
            "source_level",
            "macro_atom_transition_id",
        ],
        kind=SORTING_ALGORITHM,
    )  # This is how carsus sorted the macro atom transitions, then also using macro_atom_transition_id to break ties.

    normalized_probabilities = normalized_probabilities.loc[
        macro_atom_transition_metadata.index
    ]  # Reorder to match the metadata, which was sorted to match carsus.

    return normalized_probabilities, macro_atom_transition_metadata
