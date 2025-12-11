import numpy as np
import pandas as pd

from tardis.configuration.sorting_globals import SORTING_ALGORITHM
from tardis.io.atom_data import AtomData
from tardis.opacities.macro_atom import util
from tardis.opacities.macro_atom.base import (
    get_macro_atom_data,
    initialize_transition_probabilities,
)
from tardis.opacities.macro_atom.macroatom_continuum_transitions import (
    collisional_transition_deexc_to_k_packet,
    collisional_transition_excitation_cool,
    collisional_transition_internal_down,
    collisional_transition_internal_up,
    continuum_transition_photoionization,
    continuum_transition_recombination_emission,
    continuum_transition_recombination_internal,
    probability_collision_deexc_to_k_packet,
    probability_collision_excitation_cool,
    probability_collision_internal_down,
    probability_collision_internal_up,
    probability_photoionization,
    probability_recombination_emission,
    probability_recombination_internal,
)
from tardis.opacities.macro_atom.macroatom_line_transitions import (
    line_transition_emission_down,
    line_transition_internal_down,
    line_transition_internal_up,
    probability_emission_down,
    probability_internal_down,
    probability_internal_up,
)
from tardis.opacities.macro_atom.macroatom_state import (
    LegacyMacroAtomState,
    MacroAtomState,
)
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType


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

        self._precompute_static_data()

    def _precompute_static_data(self):
        self._oscillator_strength_ul = self.lines.f_ul.to_numpy().reshape(-1, 1)
        self._oscillator_strength_lu = self.lines.f_lu.to_numpy().reshape(-1, 1)
        self._nus = self.lines.nu.to_numpy().reshape(-1, 1)
        self._energies_upper = (
            self.levels[["energy"]]
            .reindex(self.lines.index.droplevel("level_number_lower"))
            .to_numpy()
        )
        self._energies_lower = (
            self.levels[["energy"]]
            .reindex(self.lines.index.droplevel("level_number_upper"))
            .to_numpy()
        )
        self._transition_a_i_l_u_array = self.lines.reset_index()[
            [
                "atomic_number",
                "ion_number",
                "level_number_lower",
                "level_number_upper",
            ]
        ].to_numpy()  # This is a helper array to make the source and destination columns. The letters stand for atomic_number, ion_number, lower level, upper level.

        self._lines_level_upper = self.lines.index.droplevel(
            "level_number_lower"
        )

    def solve(
        self,
        mean_intensities_blue_wing: pd.DataFrame,
        beta_sobolevs: pd.DataFrame,
        stimulated_emission_factors: np.ndarray,
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
        stimulated_emission_factors : np.ndarray
            Stimulated emission factors for the lines.

        Returns
        -------
        MacroAtomState
            A MacroAtomState object containing the transition probabilities, transition metadata,
            and a mapping from line IDs to macro atom level upper indices.
        """
        is_first_iteration = not hasattr(self, "computed_metadata")

        if is_first_iteration:
            (
                normalized_probabilities,
                macro_atom_transition_metadata,
                line2macro_level_upper,
                macro_block_references,
                references_index,
            ) = self._solve_first_macroatom_iteration(
                mean_intensities_blue_wing,
                beta_sobolevs,
                stimulated_emission_factors,
            )
        else:
            normalized_probabilities = self._solve_next_macroatom_iteration(
                mean_intensities_blue_wing,
                beta_sobolevs,
                stimulated_emission_factors,
            )
            (
                macro_atom_transition_metadata,
                line2macro_level_upper,
                macro_block_references,
                references_index,
            ) = self.computed_metadata

        return MacroAtomState(
            normalized_probabilities,
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            references_index,
        )

    def _solve_first_macroatom_iteration(
        self,
        mean_intensities_blue_wing: pd.DataFrame,
        beta_sobolevs: pd.DataFrame,
        stimulated_emission_factors: np.ndarray,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, pd.Series]:
        """
        Handle the first iteration of the solve method.

        Fully computes all metadata for the macroatom and adds it to the class with
        the computed_metadata attribute. This method performs the complete calculation
        including transition probability computation, normalization, sorting, and
        metadata preparation.

        Parameters
        ----------
        mean_intensities_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
            For more detail see Lucy 2003, https://doi.org/10.1051/0004-6361:20030357.
            Referenced as 'J^b_{lu}' internally, or 'J^b_{ji}' in the original paper.
        beta_sobolevs : pd.DataFrame
            Escape probabilities for the Sobolev approximation. These probabilities
            represent the fraction of photons that escape the line formation region
            without being reabsorbed.
        stimulated_emission_factors : np.ndarray
            Factors accounting for stimulated emission in the transitions. These
            modify the transition probabilities based on the radiation field strength.

        Returns
        -------
        normalized_probabilities : pd.DataFrame
            DataFrame containing normalized transition probabilities where each source
            group sums to 1.0.
        macro_atom_transition_metadata : pd.DataFrame
            DataFrame containing metadata for transitions including source and
            destination levels, transition types, and line indices.
        line2macro_level_upper : pd.Series
            Series mapping line transitions to macro atom level indices for upper levels.
        macro_block_references : pd.Series
            Series with unique source levels as index and their first occurrence
            index in the metadata as values.
        """
        if self.line_interaction_type in ["downbranch", "macroatom"]:
            p_emission_down, emission_down_metadata = (
                line_transition_emission_down(
                    self._oscillator_strength_ul,
                    self._nus,
                    self._energies_upper,
                    self._energies_lower,
                    beta_sobolevs,
                    self._transition_a_i_l_u_array,
                    self.lines.line_id.to_numpy(),
                )
            )
        else:
            raise ValueError(
                f"Unknown line interaction type: {self.line_interaction_type}"
            )
        if self.line_interaction_type == "downbranch":
            probabilities_df = p_emission_down
            macro_atom_transition_metadata = emission_down_metadata

        elif self.line_interaction_type == "macroatom":
            p_internal_down, internal_down_metadata = (
                line_transition_internal_down(
                    self._oscillator_strength_ul,
                    self._nus,
                    self._energies_lower,
                    beta_sobolevs,
                    self._transition_a_i_l_u_array,
                    self.lines.line_id.to_numpy(),
                )
            )

            p_internal_up, internal_up_metadata = line_transition_internal_up(
                self._oscillator_strength_lu,
                self._nus,
                self._energies_lower,
                mean_intensities_blue_wing,
                beta_sobolevs,
                stimulated_emission_factors,
                self._transition_a_i_l_u_array,
                self.lines.line_id.to_numpy(),
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
        normalized_probabilities = self.normalize_transition_probabilities(
            probabilities_df
        )

        normalized_probabilities, macro_atom_transition_metadata = (
            self.reindex_sort_and_clean_probabilities_and_metadata(
                normalized_probabilities, macro_atom_transition_metadata
            )
        )

        # We have to create the line2macro object after sorting.
        line2macro_level_upper, reference_index = (
            self.create_line2macro_level_upper_and_reference_idx(
                macro_atom_transition_metadata, self._lines_level_upper
            )
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

        self.create_source_and_destination_idx_columns(
            macro_atom_transition_metadata
        )

        macro_block_references = self.create_macro_block_references(
            macro_atom_transition_metadata
        )

        self.computed_metadata = (
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            reference_index,
        )

        return (
            normalized_probabilities,
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            reference_index,
        )

    def _solve_next_macroatom_iteration(
        self,
        mean_intensities_blue_wing: pd.DataFrame,
        beta_sobolevs: pd.DataFrame,
        stimulated_emission_factors: np.ndarray,
    ) -> pd.DataFrame:
        """
        Handle subsequent iterations of the solve method.

        Uses precomputed metadata and only recalculates the probabilities. This method
        is optimized for speed by reusing the transition metadata, block references,
        and line mappings computed in the first iteration.

        Parameters
        ----------
        mean_intensities_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
            For more detail see Lucy 2003, https://doi.org/10.1051/0004-6361:20030357.
            Referenced as 'J^b_{lu}' internally, or 'J^b_{ji}' in the original paper.
            This parameter may have updated values compared to the first iteration.
        beta_sobolevs : pd.DataFrame
            Escape probabilities for the Sobolev approximation. These probabilities
            represent the fraction of photons that escape the line formation region
            without being reabsorbed. Values may be updated from the first iteration.
        stimulated_emission_factors : np.ndarray
            Factors accounting for stimulated emission in the transitions. These
            modify the transition probabilities based on the radiation field strength.
            May contain updated values from the radiation field calculation.

        Returns
        -------
        pd.DataFrame
            DataFrame containing normalized transition probabilities where each source
            group sums to 1.0. The structure matches the first iteration output but
            with updated probability values.
        """
        (
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            reference_index,
        ) = self.computed_metadata
        line_trans_internal_up_ids = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_UP
        ].transition_line_idx.to_numpy()
        line_trans_internal_down_ids = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_DOWN
        ].transition_line_idx.to_numpy()
        line_trans_emission_down_ids = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.BB_EMISSION
        ].transition_line_idx.to_numpy()

        probabilities_df = pd.DataFrame(
            np.zeros(
                (
                    macro_atom_transition_metadata.shape[0],
                    beta_sobolevs.shape[1],
                )
            ),
            index=macro_atom_transition_metadata.index,
            columns=beta_sobolevs.columns,
        )
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.BB_EMISSION
        ] = probability_emission_down(
            beta_sobolevs.iloc[line_trans_emission_down_ids],
            self._nus[line_trans_emission_down_ids],
            self._oscillator_strength_ul[line_trans_emission_down_ids],
            self._energies_upper[line_trans_emission_down_ids],
            self._energies_lower[line_trans_emission_down_ids],
        ).to_numpy()
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_DOWN
        ] = probability_internal_down(
            beta_sobolevs.iloc[line_trans_internal_down_ids],
            self._nus[line_trans_internal_down_ids],
            self._oscillator_strength_ul[line_trans_internal_down_ids],
            self._energies_lower[line_trans_internal_down_ids],
        ).to_numpy()

        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_UP
        ] = probability_internal_up(
            beta_sobolevs.iloc[line_trans_internal_up_ids],
            self._nus[line_trans_internal_up_ids],
            self._oscillator_strength_lu[line_trans_internal_up_ids],
            stimulated_emission_factors[line_trans_internal_up_ids],
            mean_intensities_blue_wing.iloc[line_trans_internal_up_ids],
            self._energies_lower[line_trans_internal_up_ids],
        ).to_numpy()

        probabilities_df["source"] = (
            macro_atom_transition_metadata.source.values
        )
        normalized_probabilities = self.normalize_transition_probabilities(
            probabilities_df
        )

        return normalized_probabilities

    def create_source_and_destination_idx_columns(
        self,
        macro_atom_transition_metadata,
    ):
        """
        This function creates numerical indices for source and destination levels
        by mapping unique source levels to sequential integers. The destination
        indices use -99 for destinations that are not sources (emission-only levels).

        Parameters
        ----------
        macro_atom_transition_metadata : pd.DataFrame
            DataFrame containing macro atom transition metadata with 'source' and
            'destination' columns.
        """
        source_to_index = {
            source: idx
            for idx, source in enumerate(
                macro_atom_transition_metadata.source.unique()
            )
        }
        # -99 should never be used downstream. The presence of it means the destination is not a source,
        # which means that the destination is only referenced from emission
        # (or macroatom deactivation) for the given macroatom configuration.
        macro_atom_transition_metadata["destination_level_idx"] = (
            (macro_atom_transition_metadata.destination.map(source_to_index))
            .fillna(-99)
            .astype(np.int64)
        )

        macro_atom_transition_metadata["source_level_idx"] = (
            macro_atom_transition_metadata.source.map(source_to_index)
        ).astype(np.int64)

    def create_macro_block_references(self, macro_atom_transition_metadata):
        """
        Create macro block references from the macro atom transition metadata.
        This method creates a mapping from unique source levels to their first occurrence index in the metadata.

        Parameters
        ----------
        macro_atom_transition_metadata : pandas.DataFrame
            DataFrame containing metadata for macro atom transitions.

        Returns
        -------
        pandas.Series
            Series with unique source levels as index and their first occurrence index in the metadata as values.
        """
        unique_source_multi_index = pd.MultiIndex.from_tuples(
            macro_atom_transition_metadata.source.unique(),
            names=["atomic_number", "ion_number", "level_number"],
        )
        macro_data = (
            macro_atom_transition_metadata.reset_index()
            .groupby("source")
            .apply(lambda x: x.index[0])
        )

        # Append a dummy index so that the interactions can access a "block end" if a packet activates the macroatom highest level of the heaviest element in the montecarlo.
        # Without this the kernel will crash trying to access an index that doesn't exist.
        macro_data = np.append(
            macro_data.values, len(macro_atom_transition_metadata)
        )
        unique_source_multi_index = unique_source_multi_index.append(
            pd.MultiIndex.from_tuples(
                [(-99, -99, -99)],
                names=["atomic_number", "ion_number", "level_number"],
            )
        )

        macro_block_references = pd.Series(
            data=macro_data,
            index=unique_source_multi_index,
            name="macro_block_references",
        )

        return macro_block_references

    def create_line2macro_level_upper_and_reference_idx(
        self,
        macro_atom_transition_metadata: pd.DataFrame,
        lines_level_upper: pd.MultiIndex,
    ) -> tuple[pd.Series, pd.Series]:
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
        pd.Series
            Series with unique source levels as index and their assigned indices as values
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
        return line2macro_level_upper, unique_source_series

    @staticmethod
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
        normalized_probabilities.replace(np.nan, 0.0, inplace=True)

        return normalized_probabilities.drop(columns=["source"])

    def reindex_sort_and_clean_probabilities_and_metadata(
        self,
        probabilities: pd.DataFrame,
        macro_atom_transition_metadata: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Reindex and sort macro atom transition probabilities and metadata. Also creates the unique metadata ID.

        Parameters
        ----------
        probabilities : pd.DataFrame
            DataFrame containing normalized transition probabilities.
        macro_atom_transition_metadata : pd.DataFrame
            DataFrame containing metadata for macro atom transitions.

        Returns
        -------
        tuple[pd.DataFrame, pd.DataFrame]
            Reindexed normalized probabilities and cleaned metadata sorted by
            atomic number, ion number, and source level.
        """
        probabilities = probabilities.reset_index(
            drop=True
        )  # Reset to create a unique macro_atom_transition_id.
        probabilities.index.rename("macro_atom_transition_id", inplace=True)

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

        probabilities = probabilities.loc[
            macro_atom_transition_metadata.index
        ]  # Reorder to match the metadata, which was sorted to match carsus.

        return probabilities, macro_atom_transition_metadata


class ContinuumMacroAtomSolver(BoundBoundMacroAtomSolver):
    levels: pd.DataFrame
    lines: pd.DataFrame
    line_interaction_type: str
    photoionization_data: pd.DataFrame
    selected_continuum_transitions: np.ndarray

    def __init__(
        self,
        levels: pd.DataFrame,
        lines: pd.DataFrame,
        photoionization_data: pd.DataFrame,
        selected_continuum_transitions: np.ndarray = np.array([]),
        line_interaction_type: str = "macroatom",
    ) -> None:
        """
        Initialize the ContinuumMacroAtomSolver.

        Parameters
        ----------
        levels : pd.DataFrame
            DataFrame containing atomic level information.
        lines : pd.DataFrame
            DataFrame containing spectral line information.
        photoionization_data : pd.DataFrame
            DataFrame containing photoionization cross-section information.
        line_interaction_type : str, optional
            Type of line interaction to use. Default is "macroatom".
        """
        super().__init__(
            lines=lines,
            levels=levels,
            line_interaction_type=line_interaction_type,
        )
        if line_interaction_type != "macroatom":
            raise NotImplementedError(
                "ContinuumMacroAtomSolver only supports line_interaction_type='macroatom' currently."
            )

        if selected_continuum_transitions.size > 0:
            included_species = photoionization_data.index.droplevel(
                "level_number"
            ).isin(selected_continuum_transitions)
            self.photoionization_data = photoionization_data[included_species]
        else:
            self.photoionization_data = photoionization_data

        # selected_continuum_transitions = [
        #     (1, 0),
        #     (1, 1),
        # ]  # Temporary hack to test the continuum macro atom implementation.
        # included_species = photoionization_data.index.droplevel(
        #     "level_number"
        # ).isin(selected_continuum_transitions)
        # self.photoionization_data = photoionization_data[included_species]

        # Here we probably want to check and throw an error if the photoionization data contains atoms not in the lines and levels dataframes.
        self.photoionization_data_level_energies = levels.loc[
            self.photoionization_data.index.unique()
        ].energy
        self.ionization_frequency_thresholds = (
            self.photoionization_data.groupby(level=[0, 1, 2]).first().nu
        )

    def solve(
        self,
        mean_intensities_blue_wing: pd.DataFrame,
        beta_sobolevs: pd.DataFrame,
        stimulated_emission_factors: np.ndarray,
        stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
        spontaneous_recombination_coeff: pd.DataFrame,
        coll_deexc_coeff: pd.DataFrame,
        coll_exc_coeff: pd.DataFrame,
        electron_densities: pd.DataFrame,
        level_number_density: pd.DataFrame,
        delta_E_yg: pd.DataFrame,
    ) -> MacroAtomState:
        """
        Solve the Macro Atom State including continuum transitions.

        This method calculates transition probabilities for both bound-bound (line) and continuum
        transitions and returns a MacroAtomState object with the probabilities and macro atom transition metadata.
        Referenced as $p_i$ in Lucy 2003, https://doi.org/10.1051/0004-6361:20030357

        Parameters
        ----------
        mean_intensities_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
            Referenced as 'J^b_{lu}' internally, or 'J^b_{ji}' in the original paper.
        beta_sobolevs : pd.DataFrame
            Escape probabilities for the Sobolev approximation.
        stimulated_emission_factors : np.ndarray
            Stimulated emission factors for the lines.
        stim_recomb_corrected_photoionization_rate_coeff : pd.DataFrame
            Corrected photoionization rate coefficients for continuum transitions.
        spontaneous_recombination_coeff : pd.DataFrame
            Spontaneous recombination coefficients for continuum transitions.

        Returns
        -------
        MacroAtomState
            State of the macro atom including continuum transitions, ready to be placed into the OpacityState.
        """
        is_first_iteration = not hasattr(self, "computed_metadata")

        if is_first_iteration:
            self._delta_E_yg = delta_E_yg  # This should be moved to the init, but can't yet because delta_E_yg comes from the plasma which isn't given to macroatom init
            (
                normalized_probabilities,
                macro_atom_transition_metadata,
                line2macro_level_upper,
                macro_block_references,
                references_index,
            ) = self._solve_first_macroatom_iteration(
                mean_intensities_blue_wing,
                beta_sobolevs,
                stimulated_emission_factors,
                stim_recomb_corrected_photoionization_rate_coeff,
                spontaneous_recombination_coeff,
                coll_deexc_coeff,
                coll_exc_coeff,
                electron_densities,
                level_number_density,
            )
        else:
            normalized_probabilities = self._solve_next_macroatom_iteration(
                mean_intensities_blue_wing,
                beta_sobolevs,
                stimulated_emission_factors,
                stim_recomb_corrected_photoionization_rate_coeff,
                spontaneous_recombination_coeff,
                coll_deexc_coeff,
                coll_exc_coeff,
                electron_densities,
                level_number_density,
            )
            (
                macro_atom_transition_metadata,
                line2macro_level_upper,
                macro_block_references,
                references_index,
            ) = self.computed_metadata

        return MacroAtomState(
            normalized_probabilities,
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            references_index,
        )

    def _solve_first_macroatom_iteration(
        self,
        mean_intensities_blue_wing: pd.DataFrame,
        beta_sobolevs: pd.DataFrame,
        stimulated_emission_factors: np.ndarray,
        stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
        spontaneous_recombination_coeff: pd.DataFrame,
        coll_deexc_coeff: pd.DataFrame,
        coll_exc_coeff: pd.DataFrame,
        electron_densities: pd.DataFrame,
        level_number_density: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, pd.Series]:
        """
        Handle the first iteration of the solve method for continuum macro atom.

        Fully computes all metadata for the macroatom including continuum transitions and adds it to the class
        with the computed_metadata attribute. This method performs the complete calculation including transition
        probability computation, normalization, sorting, and metadata preparation for both bound-bound and bound-free transitions.

        Parameters
        ----------
        mean_intensities_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
        beta_sobolevs : pd.DataFrame
            Escape probabilities for the Sobolev approximation.
        stimulated_emission_factors : np.ndarray
            Stimulated emission factors for the lines.
        stim_recomb_corrected_photoionization_rate_coeff : pd.DataFrame
            Corrected photoionization rate coefficients for continuum transitions.

        spontaneous_recombination_coeff : pd.DataFrame
            Spontaneous recombination coefficients for continuum transitions.

        Returns
        -------
        normalized_probabilities : pd.DataFrame
            DataFrame containing normalized transition probabilities where each source group sums to 1.0.
        macro_atom_transition_metadata : pd.DataFrame
            DataFrame containing metadata for transitions including source and destination levels, transition types, and line indices.
        line2macro_level_upper : pd.Series
            Series mapping line transitions to macro atom level indices for upper levels.
        macro_block_references : pd.Series
            Series with unique source levels as index and their first occurrence index in the metadata as values.
        references_index : pd.Series
            Series with unique source levels as index and their assigned indices as values.
        """
        # Assemble bound-bound transitions first.
        p_internal_up, internal_up_metadata = line_transition_internal_up(
            self._oscillator_strength_lu,
            self._nus,
            self._energies_lower,
            mean_intensities_blue_wing,
            beta_sobolevs,
            stimulated_emission_factors,
            self._transition_a_i_l_u_array,
            self.lines.line_id.to_numpy(),
        )
        internal_up_metadata[
            "photoionization_key_idx"
        ] = -99  # Bound-bound transitions don't have continuum ids
        internal_up_metadata[
            "collision_key_idx"
        ] = -99  # Bound-bound transitions don't have collision ids
        p_internal_down, internal_down_metadata = line_transition_internal_down(
            self._oscillator_strength_ul,
            self._nus,
            self._energies_lower,
            beta_sobolevs,
            self._transition_a_i_l_u_array,
            self.lines.line_id.to_numpy(),
        )
        internal_down_metadata[
            "photoionization_key_idx"
        ] = -99  # Bound-bound transitions don't have continuum ids
        internal_down_metadata[
            "collision_key_idx"
        ] = -99  # Bound-bound transitions don't have collision ids
        p_emission_down, emission_down_metadata = line_transition_emission_down(
            self._oscillator_strength_ul,
            self._nus,
            self._energies_upper,
            self._energies_lower,
            beta_sobolevs,
            self._transition_a_i_l_u_array,
            self.lines.line_id.to_numpy(),
        )
        emission_down_metadata[
            "photoionization_key_idx"
        ] = -99  # Bound-bound transitions don't have continuum ids
        emission_down_metadata[
            "collision_key_idx"
        ] = -99  # Bound-bound transitions don't have collision ids

        # Then assemble photoionization transitions
        p_photoionization, photoionization_metadata = (
            continuum_transition_photoionization(
                stim_recomb_corrected_photoionization_rate_coeff,
                self.photoionization_data_level_energies,
            )
        )
        p_recombination_emission, recombination_emission_metadata = (
            continuum_transition_recombination_emission(
                spontaneous_recombination_coeff, self.photoionization_data.nu
            )
        )
        p_recombination_internal, recombination_internal_metadata = (
            continuum_transition_recombination_internal(
                spontaneous_recombination_coeff,
                self.photoionization_data_level_energies,
            )
        )
        # Then assemble the collisional transitions
        self._coll_energies_lower = self.levels.energy.loc[
            coll_deexc_coeff.index.droplevel("level_number_upper")
        ]  # This should probably be moved to init - static data,
        self._coll_indices_lower = coll_exc_coeff.index.droplevel(
            "level_number_upper"
        )
        # but needs coll_deexc_coeff which does change per iteration

        p_coll_down_to_k_packet, coll_down_to_packet_metadata = (
            collisional_transition_deexc_to_k_packet(
                coll_deexc_coeff,
                electron_densities,
                self._delta_E_yg,
            )
        )
        p_coll_internal_down, coll_internal_down_metadata = (
            collisional_transition_internal_down(
                coll_deexc_coeff, electron_densities, self._coll_energies_lower
            )
        )
        p_coll_internal_up, coll_internal_up_metadata = (
            collisional_transition_internal_up(
                coll_exc_coeff, electron_densities, self._coll_energies_lower
            )
        )

        p_coll_excitation_cool, coll_excitation_cool_metadata = (
            collisional_transition_excitation_cool(
                coll_exc_coeff,
                electron_densities,
                self._delta_E_yg,
                level_number_density,
                self._coll_indices_lower,
            )
        )

        probabilities_df = pd.concat(
            [
                p_emission_down,
                p_internal_down,
                p_internal_up,
                p_photoionization,
                p_recombination_emission,
                p_recombination_internal,
                p_coll_down_to_k_packet,
                p_coll_internal_down,
                p_coll_internal_up,
                p_coll_excitation_cool,
            ],
            ignore_index=True,
        )

        macro_atom_transition_metadata = pd.concat(
            [
                emission_down_metadata,
                internal_down_metadata,
                internal_up_metadata,
                photoionization_metadata,
                recombination_emission_metadata,
                recombination_internal_metadata,
                coll_down_to_packet_metadata,
                coll_internal_down_metadata,
                coll_internal_up_metadata,
                coll_excitation_cool_metadata,
            ],
            ignore_index=True,
        )

        probabilities_df, macro_atom_transition_metadata = (
            self.reindex_sort_and_clean_probabilities_and_metadata(
                probabilities_df, macro_atom_transition_metadata
            )
        )
        probabilities_df["source"] = macro_atom_transition_metadata["source"]
        normalized_probabilities = self.normalize_transition_probabilities(
            probabilities_df
        )
        # normalized_probabilities = probabilities_df

        line2macro_level_upper, reference_index = (
            self.create_line2macro_level_upper_and_reference_idx(
                macro_atom_transition_metadata, self._lines_level_upper
            )
        )
        self.create_source_and_destination_idx_columns(
            macro_atom_transition_metadata
        )
        macro_block_references = self.create_macro_block_references(
            macro_atom_transition_metadata
        )

        self.computed_metadata = (
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            reference_index,
        )

        return (
            normalized_probabilities,
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            reference_index,
        )

    def _solve_next_macroatom_iteration(
        self,
        mean_intensities_blue_wing: pd.DataFrame,
        beta_sobolevs: pd.DataFrame,
        stimulated_emission_factors: np.ndarray,
        stim_recomb_corrected_photoionization_rate_coeff: pd.DataFrame,
        spontaneous_recombination_coeff: pd.DataFrame,
        coll_deexc_coeff: pd.DataFrame,
        coll_exc_coeff: pd.DataFrame,
        electron_densities: pd.DataFrame,
        level_number_density: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Handle subsequent iterations of the solve method.

        Uses precomputed metadata and only recalculates the probabilities. This method
        is optimized for speed by reusing the transition metadata, block references,
        and line mappings computed in the first iteration.

        Parameters
        ----------
        mean_intensities_blue_wing : pd.DataFrame
            Mean intensity of the radiation field of each line in the blue wing for each shell.
            For more detail see Lucy 2003, https://doi.org/10.1051/0004-6361:20030357.
            Referenced as 'J^b_{lu}' internally, or 'J^b_{ji}' in the original paper.
            This parameter may have updated values compared to the first iteration.
        beta_sobolevs : pd.DataFrame
            Escape probabilities for the Sobolev approximation. These probabilities
            represent the fraction of photons that escape the line formation region
            without being reabsorbed. Values may be updated from the first iteration.
        stimulated_emission_factors : np.ndarray
            Factors accounting for stimulated emission in the transitions. These
            modify the transition probabilities based on the radiation field strength.
            May contain updated values from the radiation field calculation.

        Returns
        -------
        pd.DataFrame
            DataFrame containing normalized transition probabilities where each source
            group sums to 1.0. The structure matches the first iteration output but
            with updated probability values.
        """
        (
            macro_atom_transition_metadata,
            line2macro_level_upper,
            macro_block_references,
            reference_index,
        ) = self.computed_metadata
        line_trans_internal_up_ids = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_UP
        ].transition_line_idx.to_numpy()
        line_trans_internal_down_ids = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_DOWN
        ].transition_line_idx.to_numpy()
        line_trans_emission_down_ids = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.BB_EMISSION
        ].transition_line_idx.to_numpy()

        # Photoionization and recombination indices
        continuum_photoionization_idxs = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.PHOTOIONIZATION
        ].photoionization_key_idx.to_numpy()

        # collisional indices
        collisional_down_k_idxs = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_DOWN_TO_K_PACKET
        ].collision_key_idx.to_numpy()
        collisional_up_cool_idxs = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_UP_COOLING  # TODO: Something with this - Currently not used
        ].collision_key_idx.to_numpy()
        collisional_down_internal_idxs = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_DOWN_INTERNAL
        ].collision_key_idx.to_numpy()
        collisional_up_internal_idxs = macro_atom_transition_metadata[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_UP_INTERNAL
        ].collision_key_idx.to_numpy()

        probabilities_df = pd.DataFrame(
            np.zeros(
                (
                    macro_atom_transition_metadata.shape[0],
                    beta_sobolevs.shape[1],
                )
            ),
            index=macro_atom_transition_metadata.index,
            columns=beta_sobolevs.columns,
        )
        # Update the line transitions first.
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.BB_EMISSION
        ] = probability_emission_down(
            beta_sobolevs.iloc[line_trans_emission_down_ids],
            self._nus[line_trans_emission_down_ids],
            self._oscillator_strength_ul[line_trans_emission_down_ids],
            self._energies_upper[line_trans_emission_down_ids],
            self._energies_lower[line_trans_emission_down_ids],
        ).to_numpy()
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_DOWN
        ] = probability_internal_down(
            beta_sobolevs.iloc[line_trans_internal_down_ids],
            self._nus[line_trans_internal_down_ids],
            self._oscillator_strength_ul[line_trans_internal_down_ids],
            self._energies_lower[line_trans_internal_down_ids],
        ).to_numpy()

        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.INTERNAL_UP
        ] = probability_internal_up(
            beta_sobolevs.iloc[line_trans_internal_up_ids],
            self._nus[line_trans_internal_up_ids],
            self._oscillator_strength_lu[line_trans_internal_up_ids],
            stimulated_emission_factors[line_trans_internal_up_ids],
            mean_intensities_blue_wing.iloc[line_trans_internal_up_ids],
            self._energies_lower[line_trans_internal_up_ids],
        ).to_numpy()

        # Then update the continuum transitions.
        photoionization_sources = pd.MultiIndex.from_tuples(
            macro_atom_transition_metadata[
                macro_atom_transition_metadata.transition_type
                == MacroAtomTransitionType.PHOTOIONIZATION
            ].source
        )
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.PHOTOIONIZATION
        ] = probability_photoionization(
            stim_recomb_corrected_photoionization_rate_coeff.loc[
                photoionization_sources
            ],
            self.photoionization_data_level_energies.iloc[
                continuum_photoionization_idxs
            ],
        ).to_numpy()

        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.BF_EMISSION
        ] = probability_recombination_emission(
            spontaneous_recombination_coeff,
            self.photoionization_data.nu,
        ).to_numpy()  # WARNING - This assumes that the recombination probabilities do not get reordered within a recomb_emission block,
        # because they are matched again against the photoionization nus.
        # There's no reason for them to get reordered, but they are not explicitly tracked.
        # Also, you'd have to add more metadata to track them explicitly, which could be done but isn't free.

        recomb_internal_destinations = pd.MultiIndex.from_tuples(
            macro_atom_transition_metadata[
                macro_atom_transition_metadata.transition_type
                == MacroAtomTransitionType.RECOMB_INTERNAL
            ].destination
        )
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.RECOMB_INTERNAL
        ] = probability_recombination_internal(
            spontaneous_recombination_coeff.loc[recomb_internal_destinations],
            self.photoionization_data_level_energies.loc[
                recomb_internal_destinations
            ],
        ).to_numpy()

        # verify these below are right
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_DOWN_TO_K_PACKET
        ] = probability_collision_deexc_to_k_packet(
            coll_deexc_coeff.iloc[collisional_down_k_idxs],
            electron_densities,
            self._delta_E_yg.iloc[collisional_down_k_idxs],
        ).to_numpy()

        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_DOWN_INTERNAL
        ] = probability_collision_internal_down(
            coll_deexc_coeff.iloc[collisional_down_internal_idxs],
            electron_densities,
            self._coll_energies_lower.iloc[collisional_down_internal_idxs],
        ).to_numpy()

        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_UP_INTERNAL
        ] = probability_collision_internal_up(
            coll_exc_coeff.iloc[collisional_up_internal_idxs],
            electron_densities,
            self._coll_energies_lower.iloc[collisional_up_internal_idxs],
        ).to_numpy()

        # collisional indices are hard here because we sum along an axis
        # if they don't get reordered this should be fine
        probabilities_df[
            macro_atom_transition_metadata.transition_type
            == MacroAtomTransitionType.COLL_UP_COOLING
        ] = probability_collision_excitation_cool(
            coll_exc_coeff,
            electron_densities,
            self._delta_E_yg,
            level_number_density,
            self._coll_indices_lower,
        ).to_numpy()

        probabilities_df["source"] = (
            macro_atom_transition_metadata.source.values
        )  # Normalize by source in the next line, so need source column.
        normalized_probabilities = self.normalize_transition_probabilities(
            probabilities_df
        )
        return normalized_probabilities

    def reindex_sort_and_clean_probabilities_and_metadata(
        self,
        probabilities: pd.DataFrame,
        macro_atom_transition_metadata: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Adapted for continuum macroatom, where continuum transitions do not have the same dataframe indices.
        Reindex and sort macro atom transition probabilities and metadata. Also creates the unique metadata ID.

        Parameters
        ----------
        probabilities : pd.DataFrame
            DataFrame containing normalized transition probabilities.
        macro_atom_transition_metadata : pd.DataFrame
            DataFrame containing metadata for macro atom transitions.

        Returns
        -------
        tuple[pd.DataFrame, pd.DataFrame]
            Reindexed normalized probabilities and cleaned metadata sorted by
            atomic number, ion number, and source level.
        """
        probabilities.index.rename("macro_atom_transition_id", inplace=True)
        macro_atom_transition_metadata.index.rename(
            "macro_atom_transition_id", inplace=True
        )

        macro_atom_transition_metadata = (
            macro_atom_transition_metadata.sort_values(
                [
                    "source",
                    "macro_atom_transition_id",
                ],
                kind=SORTING_ALGORITHM,
            )
        )

        probabilities = probabilities.loc[
            macro_atom_transition_metadata.index
        ]  # Reorder to match the metadata

        return probabilities, macro_atom_transition_metadata
