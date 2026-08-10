from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.sparse import coo_matrix


class LevelRateMatrix:
    """Construct bound-bound rate matrices from precomputed rates."""

    def __init__(
        self,
        levels: pd.DataFrame,
    ):
        """Initialize a level rate matrix builder from precomputed rate data.

        Parameters
        ----------
        levels : pd.DataFrame
            DataFrame of energy levels.
        """
        self.levels = levels

    @staticmethod
    def _validate_rates(rates: pd.DataFrame) -> None:
        """Check that a transition-rate table contains usable rates.

        Rate matrices assume finite, non-negative off-diagonal transition
        rates. The negative diagonal entries are created later from outgoing
        rates, so they are not expected in this input table.

        Parameters
        ----------
        rates : pandas.DataFrame
            Transition-rate coefficients or rates, with transitions in rows
            and cells in columns.
        """
        values = rates.to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("Transition rates must be finite and non-negative")

    def _combine_rates(
        self,
        radiative_rates_df: pd.DataFrame,
        collisional_rates_df: pd.DataFrame,
        electron_number_density: npt.ArrayLike,
    ) -> pd.DataFrame:
        """Align and combine radiative and electron-impact rates.

        Missing transitions are treated as zero in either input. The
        collisional coefficients are multiplied cell-by-cell by the electron
        number density before being added to the radiative rates, producing
        the total rate used in ``dn/dt = R @ n``.

        Parameters
        ----------
        radiative_rates_df, collisional_rates_df : pandas.DataFrame
            Transition tables with the same cell columns and compatible
            transition indexes.
        electron_number_density : array-like
            Electron number density for each cell.

        Returns
        -------
        pandas.DataFrame
            The union of the input transition indexes with total rates in
            each cell column.
        """
        electron_number_density = np.asarray(electron_number_density)
        if (
            not np.isfinite(electron_number_density).all()
            or (electron_number_density < 0).any()
        ):
            raise ValueError(
                "Electron number density must be finite and non-negative"
            )
        self._validate_rates(radiative_rates_df)
        self._validate_rates(collisional_rates_df)
        all_indexes = sorted(
            set(radiative_rates_df.index) | set(collisional_rates_df.index)
        )
        radiative_rates = radiative_rates_df.reindex(all_indexes, fill_value=0)
        collisional_rates = collisional_rates_df.reindex(
            all_indexes, fill_value=0
        ).mul(electron_number_density, axis="columns")
        rates = radiative_rates.add(collisional_rates, fill_value=0)
        self._validate_rates(rates)
        return rates

    def build_raw_rate_matrices(
        self,
        radiative_rates_df: pd.DataFrame,
        collisional_rates_df: pd.DataFrame,
        electron_number_density: npt.ArrayLike,
    ) -> pd.DataFrame:
        """Build unnormalized level transition rate matrices.

        The input rows describe transitions and the input columns describe
        cells. For each cell, this method puts the transition rate from level
        ``s`` to level ``d`` at ``R[d, s]``. The diagonal is the negative
        total outgoing rate, so each raw rate-matrix column sums to zero. For
        a level-population column vector ``n``, the statistical-equilibrium
        equation is therefore ``dn/dt = R @ n`` or, at steady state,
        ``R @ n = 0``. The total transition rate is assembled as
        ``radiative_rate + electron_number_density * collisional_coefficient``.

        The returned DataFrame is a container of matrices rather than a
        two-dimensional rate matrix itself. Its index identifies the ion and
        its columns identify the cell, for example::

            rate_matrices.loc[(6, 1), 0]
            # array([[-r_0_to_1, r_1_to_0],
            #        [ r_0_to_1, -r_1_to_0]])

        Here ``(6, 1)`` is ``(atomic_number, ion_number)`` and ``0`` is the
        cell label. Within each stored array, rows are destinations and
        columns are sources.

        Parameters
        ----------
        radiative_rates_df : pandas.DataFrame
            Bound-bound radiative rates.
        collisional_rates_df : pandas.DataFrame
            Bound-bound collisional rate coefficients.
        electron_number_density : numpy.ndarray
            Electron number density for each cell.

        Returns
        -------
        pandas.DataFrame
            A DataFrame indexed by ``(atomic_number, ion_number)`` and with
            one column per cell. Each value is an unnormalized square NumPy
            array whose rows are destination levels and columns are source
            levels.
        """
        rates = self._combine_rates(
            radiative_rates_df,
            collisional_rates_df,
            np.asarray(electron_number_density),
        )
        grouped_rates = rates.groupby(level=("atomic_number", "ion_number"))
        species_index = pd.MultiIndex.from_tuples(
            list(grouped_rates.groups.keys()),
            names=["atomic_number", "ion_number"],
        )
        raw_rate_matrices = pd.DataFrame(
            index=species_index, columns=rates.columns
        )

        for species_id, species_rates in grouped_rates:
            number_of_levels = len(self.levels.energy.loc[species_id])
            source = species_rates.index.get_level_values("level_number_source")
            destination = species_rates.index.get_level_values(
                "level_number_destination"
            )
            for cell in rates.columns:
                raw_rate_matrix = coo_matrix(
                    (
                        species_rates[cell].to_numpy(),
                        # Matrix rows are destinations; columns are sources.
                        (destination, source),
                    ),
                    shape=(number_of_levels, number_of_levels),
                ).toarray()
                np.fill_diagonal(
                    raw_rate_matrix,
                    -np.sum(raw_rate_matrix, axis=0),
                )
                raw_rate_matrices.loc[species_id, cell] = raw_rate_matrix

        raw_rate_matrices.index.names = ["atomic_number", "ion_number"]
        return raw_rate_matrices

    def build_raw_level_matrices_for_ion(
        self,
        radiative_rates_df: pd.DataFrame,
        collisional_rates_df: pd.DataFrame,
        electron_number_density: npt.ArrayLike,
        species_id: tuple[int, int],
    ) -> pd.DataFrame:
        """Build the raw level-rate block for one explicitly selected ion.

        This is the single-ion form of
        :meth:`build_raw_rate_matrices`. It preserves the same source-to-
        destination convention and returns a one-row DataFrame, so the block
        can be embedded into a larger elemental matrix without reconstructing
        its local level numbering. For example, for ``species_id == (6, 1)``
        and one cell, the result has the shape::

            result.index       result.columns
            MultiIndex([(6, 1)])   Index([0])

            result.loc[(6, 1), 0]
            # array([[-r_0_to_1, r_1_to_0],
            #        [ r_0_to_1, -r_1_to_0]])

        The block satisfies ``dn_levels/dt = R_level @ n_levels``. It is
        unnormalized: its columns sum to zero and no population-normalization
        row has been inserted.

        Parameters
        ----------
        radiative_rates_df : pandas.DataFrame
            Bound-bound radiative rates indexed by transition identifiers.
        collisional_rates_df : pandas.DataFrame
            Bound-bound collisional rate coefficients indexed by transition
            identifiers.
        electron_number_density : numpy.ndarray
            Electron number density for each cell.
        species_id : tuple[int, int]
            ``(atomic_number, ion_number)`` identifying the ion to select.

        Returns
        -------
        pandas.DataFrame
            A one-row DataFrame containing one raw square level-rate matrix
            per cell.
        """
        species_mask = (
            radiative_rates_df.index.get_level_values("atomic_number")
            == species_id[0]
        ) & (
            radiative_rates_df.index.get_level_values("ion_number")
            == species_id[1]
        )
        radiative_rates_df = radiative_rates_df.loc[species_mask]
        species_mask = (
            collisional_rates_df.index.get_level_values("atomic_number")
            == species_id[0]
        ) & (
            collisional_rates_df.index.get_level_values("ion_number")
            == species_id[1]
        )
        collisional_rates_df = collisional_rates_df.loc[species_mask]
        raw_rate_matrices = self.build_raw_rate_matrices(
            radiative_rates_df,
            collisional_rates_df,
            electron_number_density,
        )
        return raw_rate_matrices.loc[[species_id]]

    def solve(
        self,
        radiative_rates_df: pd.DataFrame,
        collisional_rates_df: pd.DataFrame,
        electron_number_density: npt.ArrayLike,
    ) -> pd.DataFrame:
        """Construct the compiled rate matrix dataframe.

        Parameters
        ----------
        radiative_rates_df : pd.DataFrame
            Precomputed radiative rates indexed by transition identifiers,
            with each column being a cell.
        collisional_rates_df : pd.DataFrame
            Precomputed thermal collisional rates indexed by transition
            identifiers, with each column being a cell.
        electron_number_density : pandas.Series | numpy.ndarray
            Electron number density indexed by cell.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        """
        raw_rate_matrices = self.build_raw_rate_matrices(
            radiative_rates_df,
            collisional_rates_df,
            electron_number_density,
        )
        normalized_rate_matrices = raw_rate_matrices.copy(deep=True)
        for species_id in raw_rate_matrices.index:
            for cell in raw_rate_matrices.columns:
                raw_rate_matrix = raw_rate_matrices.loc[species_id, cell].copy()
                raw_rate_matrix[0, :] = 1
                normalized_rate_matrices.loc[species_id, cell] = raw_rate_matrix

        normalized_rate_matrices.index.names = [
            "atomic_number",
            "ion_number",
        ]

        return normalized_rate_matrices


@dataclass(frozen=True)
class ElementalStateIndex:
    """Positions of explicit level and total-ion states in an element matrix."""

    atomic_number: int
    states: pd.MultiIndex
    level_positions: dict[tuple[int, int], int]
    ion_positions: dict[int, int]

    @property
    def size(self) -> int:
        """Return the number of states in the elemental matrix."""
        return len(self.states)


@dataclass(frozen=True)
class IonLevelRateMatrixSet:
    """Elemental rate matrices containing explicit levels and ion states."""

    normalized_rate_matrices: pd.DataFrame
    raw_elemental_rate_matrices: pd.DataFrame
    state_index: ElementalStateIndex


class IonRateMatrix:
    """Construct ionization rate matrices from precomputed rates."""

    def __init__(self):
        """Initialize the rate matrix builder."""

    @staticmethod
    def _transition_position(
        state_index: ElementalStateIndex,
        ion_number: int,
        level_number: int,
    ) -> int:
        """Map an ion/level coordinate to a global matrix position.

        Level-bearing ions use their explicit level position. An ion without
        an explicit level block, including the bare nucleus, is represented
        by one total-ion position, regardless of the level number supplied by
        a bound-free transition.
        """
        # Explicit level states have positions; remaining ions use one total state.
        level_position = state_index.level_positions.get(
            (ion_number, level_number)
        )
        if level_position is not None:
            return level_position
        return state_index.ion_positions[ion_number]

    @staticmethod
    def _remove_empty_rate_frames(
        rate_frames: tuple[pd.DataFrame | None, ...],
    ) -> tuple[pd.DataFrame, ...]:
        """Discard absent optional rate tables while preserving their order."""
        return tuple(frame for frame in rate_frames if frame is not None)

    def _build_state_index(
        self,
        atomic_number: int,
        raw_level_rate_matrices: pd.DataFrame | None,
        ion_stage_count: int,
    ) -> ElementalStateIndex:
        """Assign stable global positions to level and total-ion states.

        Explicit level blocks are placed first, grouped by ion and retaining
        each block's local level order. Every remaining charge state receives
        one total-population state, including the terminal bare nucleus. The
        resulting positions are used for both rows and columns, so a matrix
        entry ``A[d, s]`` always maps from the source state ``s`` to the
        destination state ``d``.

        Parameters
        ----------
        atomic_number : int
            Element whose states are being indexed.
        raw_level_rate_matrices : pandas.DataFrame, optional
            Per-ion raw level matrices. Only rows for ``atomic_number`` are
            considered.
        ion_stage_count : int
            Number of charge states, including the bare nucleus.

        Returns
        -------
        ElementalStateIndex
            Labels and lookup dictionaries for every global state.
        """
        level_positions: dict[tuple[int, int], int] = {}
        ion_positions: dict[int, int] = {}
        state_labels: list[tuple[str, int, int]] = []

        level_species = (
            []
            if raw_level_rate_matrices is None
            else sorted(
                species_id
                for species_id in raw_level_rate_matrices.index
                if species_id[0] == atomic_number
            )
        )
        position = 0
        for _, ion_number in level_species:
            # A level matrix is stored in one DataFrame cell; its array length
            # defines the local levels that become this ion's global block.
            block = raw_level_rate_matrices.loc[
                (atomic_number, ion_number)
            ].iloc[0]
            for level_number in range(len(block)):
                # Rows and columns use the same position for this state.
                level_positions[(ion_number, level_number)] = position
                state_labels.append(("level", ion_number, level_number))
                position += 1

        for ion_number in range(ion_stage_count):
            if any(ion_number == level_ion for _, level_ion in level_species):
                continue
            # Ions without explicit levels are represented by one total state.
            ion_positions[ion_number] = position
            state_labels.append(("ion", ion_number, -1))
            position += 1

        states = pd.MultiIndex.from_tuples(
            state_labels,
            names=["state_type", "ion_number", "level_number"],
        )
        return ElementalStateIndex(
            atomic_number=atomic_number,
            states=states,
            level_positions=level_positions,
            ion_positions=ion_positions,
        )

    @staticmethod
    def _validate_transition_values(rate_frame: pd.DataFrame) -> None:
        """Validate the non-negative rates about to be added to a matrix."""
        values = rate_frame.to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("Transition rates must be finite and non-negative")

    def _add_transition_rates(
        self,
        raw_elemental_rate_matrix: np.ndarray,
        state_index: ElementalStateIndex,
        rate_frames: tuple[pd.DataFrame, ...],
        cell_idx: int,
    ) -> None:
        """Add bound-free transitions to one elemental matrix in place.

        Each input row identifies a source and destination ion, optionally
        with level numbers. The state index resolves those coordinates to
        global positions. A transition rate ``r`` contributes ``+r`` at the
        destination row/source column and ``-r`` to the source diagonal,
        implementing the population equation ``dp/dt = A @ p`` while
        preserving zero column sums.

        Parameters
        ----------
        raw_elemental_rate_matrix : numpy.ndarray
            Matrix for one cell, modified in place.
        state_index : ElementalStateIndex
            Global level and total-ion position mappings.
        rate_frames : tuple[pandas.DataFrame, ...]
            Bound-free rate tables in source-to-destination index format.
        cell_idx : int
            Positional cell index to read from each rate table.
        """
        for rate_frame in rate_frames:
            self._validate_transition_values(rate_frame)
            for transition_index, values in rate_frame.iterrows():
                transition = dict(
                    zip(
                        rate_frame.index.names,
                        transition_index,
                        strict=True,
                    )
                )
                source_ion = int(transition["ion_number_source"])
                destination_ion = int(transition["ion_number_destination"])
                source_level = int(transition.get("level_number_source", 0))
                destination_level = int(
                    transition.get("level_number_destination", 0)
                )
                source_position = self._transition_position(
                    state_index, source_ion, source_level
                )
                destination_position = self._transition_position(
                    state_index, destination_ion, destination_level
                )
                rate = float(values.iloc[cell_idx])
                # Rows are destinations and columns are sources: +r transfers
                # population into the destination and -r removes it from the source.
                raw_elemental_rate_matrix[
                    destination_position, source_position
                ] += rate
                raw_elemental_rate_matrix[source_position, source_position] -= (
                    rate
                )

    def solve_ion_and_level(
        self,
        atomic_number: int,
        raw_level_rate_matrices: pd.DataFrame | None = None,
        photoion_rates_df: pd.DataFrame | None = None,
        recomb_rates_df: pd.DataFrame | None = None,
        collisional_ionization_rates_df: pd.DataFrame | None = None,
        collision_recombination_rates_df: pd.DataFrame | None = None,
        nebular_rates_df: pd.DataFrame | None = None,
        ion_stage_count: int | None = None,
    ) -> IonLevelRateMatrixSet:
        """Build an element-normalized matrix with explicit level blocks.

        ``raw_level_rate_matrices`` must be produced by
        :meth:`LevelRateMatrix.build_raw_rate_matrices`. All ionization-rate
        frames use the existing source-to-destination transition index.
        Missing higher ionization edges can be supplied as
        ``nebular_rates_df``.

        Every stored elemental matrix uses the same column-vector rate
        equation as the level matrices. If ``p`` contains the populations of
        the explicit levels and total-ion states, then each raw matrix ``A``
        represents ``dp/dt = A @ p``. For a source state ``s`` and destination
        state ``d``, a rate ``r_s_to_d`` contributes ``+r_s_to_d`` to
        ``A[d, s]`` and ``-r_s_to_d`` to ``A[s, s]``. Consequently, the
        stationary equations are ``A @ p = 0`` and every raw column sums to
        zero. The normalized matrix replaces one dependent balance row with
        ones, so the linear solve uses ``A_normalized @ p = b`` where ``b`` is
        zero except for a one in row 0.

        The returned object contains two DataFrames. Each has one row for the
        element and one column per cell; the value at ``.loc[6, 0]`` is the
        NumPy array for carbon in cell 0, for example::

            normalized_rate_matrices.loc[6, 0]  # shape (N_states, N_states)
            raw_elemental_rate_matrices.loc[6, 0]  # same state ordering

        ``state_index.states`` records the corresponding row and column
        ordering. For a carbon calculation with explicit C I levels and no
        other level blocks, it can look like::

            [('level', 0, 0), ('level', 0, 1),
             ('ion', 1, -1), ('ion', 2, -1), ..., ('ion', 6, -1)]

        Thus a bound-bound C I transition is placed within the leading level
        block, while a transition to C II uses the C II total-state position.

        Parameters
        ----------
        atomic_number : int
            Atomic number of the element to assemble.
        raw_level_rate_matrices : pandas.DataFrame, optional
            Raw level rate matrices indexed by (atomic_number, ion_number).
        photoion_rates_df, recomb_rates_df : pandas.DataFrame, optional
            Radiative bound-free rates.
        collisional_ionization_rates_df, collision_recombination_rates_df : pandas.DataFrame, optional
            Collisional bound-free rates.
        nebular_rates_df : pandas.DataFrame, optional
            Rates for uncovered ionization edges.
        ion_stage_count : int, optional
            Number of charge states, including the bare nucleus. Defaults to
            atomic_number + 1.

        Returns
        -------
        ElementalLevelIonRateMatrixSet
            Normalized and raw elemental rate-matrix DataFrames plus the
            explicit state-index metadata needed to interpret their arrays.
        """
        if ion_stage_count is None:
            ion_stage_count = atomic_number + 1
        if ion_stage_count < 1:
            raise ValueError("ion_stage_count must be positive")
        if raw_level_rate_matrices is not None:
            selected_raw_level_rate_matrices = raw_level_rate_matrices.loc[
                raw_level_rate_matrices.index.get_level_values("atomic_number")
                == atomic_number
            ]
        else:
            selected_raw_level_rate_matrices = None
        if selected_raw_level_rate_matrices is not None and any(
            ion_number < 0 or ion_number >= ion_stage_count - 1
            for _, ion_number in selected_raw_level_rate_matrices.index
        ):
            raise ValueError(
                "A level block must represent a non-terminal ion stage"
            )

        frames = self._remove_empty_rate_frames(
            (
                raw_level_rate_matrices,
                photoion_rates_df,
                recomb_rates_df,
                collisional_ionization_rates_df,
                collision_recombination_rates_df,
                nebular_rates_df,
            )
        )
        columns = frames[0].columns
        for frame in frames[1:]:
            if not frame.columns.equals(columns):
                raise ValueError(
                    "All rate frames must use the same cell columns"
                )

        state_index = self._build_state_index(
            atomic_number,
            selected_raw_level_rate_matrices,
            ion_stage_count,
        )
        raw_elemental_rate_matrices = pd.DataFrame(
            index=pd.Index([atomic_number], name="atomic_number"),
            columns=columns,
        )
        normalized_rate_matrices = raw_elemental_rate_matrices.copy(deep=True)
        rate_frames = self._remove_empty_rate_frames(
            (
                photoion_rates_df,
                recomb_rates_df,
                collisional_ionization_rates_df,
                collision_recombination_rates_df,
                nebular_rates_df,
            )
        )
        rate_frames = tuple(
            frame.loc[
                frame.index.get_level_values("atomic_number") == atomic_number
            ]
            for frame in rate_frames
        )

        for cell_idx, cell in enumerate(columns):
            raw_elemental_rate_matrix = np.zeros(
                (state_index.size, state_index.size)
            )
            if selected_raw_level_rate_matrices is not None:
                for species_id in selected_raw_level_rate_matrices.index:
                    block = np.asarray(
                        selected_raw_level_rate_matrices.loc[species_id, cell]
                    )
                    column_residual = block.sum(axis=0)
                    roundoff_tolerance = (
                        np.finfo(float).eps
                        * block.shape[0]
                        * np.abs(block).sum(axis=0)
                    )
                    if (
                        not np.isfinite(block).all()
                        or (block - np.diag(np.diag(block)) < 0).any()
                        or np.any(np.abs(column_residual) > roundoff_tolerance)
                    ):
                        raise ValueError(
                            f"Invalid raw level rate matrix for {species_id}"
                        )
                    level_ion = int(species_id[1])
                    level_positions = [
                        state_index.level_positions[(level_ion, level_number)]
                        for level_number in range(len(block))
                    ]
                    # Embed local level rows/columns into the global state matrix.
                    raw_elemental_rate_matrix[
                        np.ix_(level_positions, level_positions)
                    ] += block
            self._add_transition_rates(
                raw_elemental_rate_matrix,
                state_index,
                rate_frames,
                cell_idx,
            )
            if not np.isfinite(raw_elemental_rate_matrix).all():
                raise ValueError("The assembled rate matrix is not finite")
            raw_elemental_rate_matrices.loc[atomic_number, cell] = (
                raw_elemental_rate_matrix
            )
            normalized_rate_matrix = raw_elemental_rate_matrix.copy()
            # Replace the first balance equation with elemental normalization.
            normalized_rate_matrix[0, :] = 1.0
            normalized_rate_matrices.loc[atomic_number, cell] = (
                normalized_rate_matrix
            )

        raw_elemental_rate_matrices.attrs["state_index"] = state_index
        normalized_rate_matrices.attrs["state_index"] = state_index
        return IonLevelRateMatrixSet(
            normalized_rate_matrices,
            raw_elemental_rate_matrices,
            state_index,
        )

    def __calculate_total_grouped_rates(self, rates_df):
        """Calculate total rates from photoionization and recombination rates.

        Parameters
        ----------
        rates_df : pd.DataFrame
            DataFrame of rates indexed by atomic number and ion number,
            with each column being a cell.

        Returns
        -------
        pd.DataFrame
            A DataFrame of grouped total rates indexed by atomic number and ion number,
            with each column being a cell.
        """
        return (
            rates_df.groupby(
                level=(
                    "atomic_number",
                    "ion_number",
                    "ion_number_source",
                    "ion_number_destination",
                )
            )
            .sum()
            .groupby(level=("atomic_number"))
        )

    def __construct_rate_matrix(self, rate, cell, ion_states):
        """Construct a sparse rate matrix from the rates.

        Parameters
        ----------
        rate : pd.DataFrame
            Rate DataFrame indexed by atomic number and ion number
        cell : int
            Cell index
        ion_states : int
            Number of ion states for the atomic number

        Returns
        -------
        coo_matrix
            A sparse matrix representing the ionization rate for the given cell.
        """
        return coo_matrix(
            (
                rate[cell],
                (
                    rate.index.get_level_values("ion_number_destination"),
                    rate.index.get_level_values("ion_number_source"),
                ),
            ),
            shape=(ion_states, ion_states),
        )

    def solve(
        self,
        photoion_rates_df: pd.DataFrame,
        recomb_rates_df: pd.DataFrame,
        collisional_ionization_rates_df: pd.DataFrame,
        collision_recombination_rates_df: pd.DataFrame,
        charge_conservation: bool = False,
    ) -> pd.DataFrame:
        """Compute the ionization rate matrix.

        Parameters
        ----------
        photoion_rates_df : pd.DataFrame
            Precomputed photoionization rates.
        recomb_rates_df : pd.DataFrame
            Precomputed radiative recombination rates.
        collisional_ionization_rates_df : pd.DataFrame
            Precomputed collisional ionization rates.
        collision_recombination_rates_df : pd.DataFrame
            Precomputed collisional recombination rates.
        charge_conservation : bool, optional
            Whether to include a charge conservation row in the rate matrix.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell. Each entry is a numpy array.
        """
        grouped_photoion_rates_df = self.__calculate_total_grouped_rates(
            photoion_rates_df
        )
        grouped_recomb_rates_df = self.__calculate_total_grouped_rates(
            recomb_rates_df
        )

        grouped_collisional_ionization_rates_df = (
            self.__calculate_total_grouped_rates(
                collisional_ionization_rates_df
            )
        )
        grouped_collisional_recombination_rates_df = (
            self.__calculate_total_grouped_rates(
                collision_recombination_rates_df
            )
        )

        rate_matrices = pd.DataFrame(
            index=list(grouped_photoion_rates_df.groups.keys()),
            columns=photoion_rates_df.columns,
        )

        for atomic_number in grouped_photoion_rates_df.groups.keys():
            photoion_rates = grouped_photoion_rates_df.get_group(atomic_number)
            recomb_rates = grouped_recomb_rates_df.get_group(atomic_number)
            coll_ion_rates = grouped_collisional_ionization_rates_df.get_group(
                atomic_number
            )
            recomb_ion_rates = (
                grouped_collisional_recombination_rates_df.get_group(
                    atomic_number
                )
            )
            max_ion_number = max(
                photoion_rates.index.get_level_values(
                    "ion_number_destination"
                ).max(),
                photoion_rates.index.get_level_values(
                    "ion_number_source"
                ).max(),
            )
            ion_states = int(max_ion_number) + 1
            for shell in range(len(photoion_rates.columns)):
                photoion_matrix = self.__construct_rate_matrix(
                    photoion_rates, shell, ion_states
                )
                recomb_matrix = self.__construct_rate_matrix(
                    recomb_rates, shell, ion_states
                )
                coll_ion_matrix = self.__construct_rate_matrix(
                    coll_ion_rates, shell, ion_states
                )
                coll_recomb_matrix = self.__construct_rate_matrix(
                    recomb_ion_rates, shell, ion_states
                )

                matrix_array = (
                    photoion_matrix
                    + recomb_matrix
                    + coll_ion_matrix
                    + coll_recomb_matrix
                ).toarray()
                np.fill_diagonal(matrix_array, -np.sum(matrix_array, axis=0))
                matrix_array[1, :] = 1
                if charge_conservation:
                    charge_conservation_row = np.hstack(
                        (np.arange(0, ion_states), -1)
                    )
                    matrix_array = np.pad(matrix_array, ((0, 0), (0, 1)))
                    matrix_array = np.vstack(
                        (charge_conservation_row, matrix_array)
                    )
                rate_matrices.loc[atomic_number, shell] = matrix_array

        rate_matrices.index.names = ["atomic_number"]

        return rate_matrices


class LTEIonRateMatrix:
    """Constructs ionization rate matrices based on LTE ratios."""

    def __init__(self):
        """Initialize the solver."""

    @staticmethod
    def _prepare_phi(phi, ion_index):
        """Prepare phi (Saha ratios) by reindexing to full ionization structure.

        Parameters
        ----------
        phi : pd.DataFrame
            Saha ratios indexed by (atomic_number, ion_number), columns are cells.
        ion_index : pd.MultiIndex
            Full ionization state index.

        Returns
        -------
        pd.DataFrame
            Phi reindexed to ion_index, with NaNs filled.
        """
        # Check for NaNs
        no_nans = pd.isnull(phi).sum().sum()
        if no_nans:
            phi = phi.fillna(phi.min().min())

        # Zero phi values cause numerical issues. Replace with small values.
        phi_min = phi[phi > 0.0].min().min()
        phi = phi.copy()
        phi[phi == 0.0] = 1.0e-10 * phi_min

        # Shift ion number by 1 and reindex to match ion population structure
        atomic_number = phi.index.get_level_values(0).values
        ion_number = phi.index.get_level_values(1).values
        new_index = pd.MultiIndex.from_arrays([atomic_number, ion_number - 1])
        phi_prep = phi.set_index(new_index).reindex(ion_index).fillna(0.0)
        return phi_prep

    @staticmethod
    def _get_number_conservation_index(ion_index):
        """Get indices for number conservation constraint row.

        Parameters
        ----------
        ion_index : pd.MultiIndex
            Index with (atomic_number, ion_number).

        Returns
        -------
        tuple
            (row_indices, col_indices) for setting number conservation row.
        """
        atomic_number = np.unique(
            ion_index.get_level_values(0).values.astype(int)
        )
        sum1 = (atomic_number + 1).cumsum() - 1
        index1 = np.concatenate(
            [
                np.ones(j + 1, dtype=int) * i
                for i, j in zip(sum1, atomic_number, strict=True)
            ]
        )
        index2 = np.arange(len(ion_index), dtype=int)
        return (index1, index2)

    def solve(
        self,
        phi,
        partition_function,
        electron_number_density,
        charge_conservation=False,
    ):
        """Compute ionization rate matrices from LTE Saha ratios.

        This produces a single rate matrix for each cell that governs ionization equilibrium
        for all species. The rate matrix includes electron density as a semi-variable,
        with number and charge conservation constraints.

        Parameters
        ----------
        phi : pd.DataFrame
            Saha ratios indexed by (atomic_number, ion_number), columns are cells.
        partition_function : pd.DataFrame
            Partition functions indexed by (atomic_number, ion_number), columns are cells.
        electron_number_density : np.ndarray
            Electron number density indexed by cell.
        charge_conservation : bool, optional
            Whether to include a charge conservation row in the rate matrix.

        Returns
        -------
        list of np.ndarray
            List of rate matrices, one per cell. Each matrix has shape
            (num_ions + 1, num_ions + 1) where num_ions is the total number of ion states
            and +1 accounts for electron density.
        """
        ion_index = partition_function.index

        # Prepare phi for use in rate matrix construction
        phi_prep = self._prepare_phi(phi, ion_index)

        # Get constraint information
        number_conservation_index = self._get_number_conservation_index(
            ion_index
        )

        # Construct rate matrix for each cell
        rate_matrices = pd.DataFrame(
            index=phi.index.get_level_values(0),
            columns=phi_prep.columns,
        )

        for atomic_number in phi.index.get_level_values(0):
            ion_states = atomic_number + 1
            for shell in range(len(phi_prep.columns)):
                # Get Saha ratios for this cell
                lte_diag = (
                    -phi_prep[shell].values / electron_number_density[shell]
                )
                lte_offdiag = (lte_diag != 0).astype(float)[:-1]

                matrix_array = np.diag(lte_diag) + np.diag(lte_offdiag, k=1)

                matrix_array[number_conservation_index] = 1.0
                if charge_conservation:
                    charge_conservation_row = np.hstack(
                        (np.arange(0, ion_states), -1)
                    )
                    matrix_array = np.pad(matrix_array, ((0, 0), (0, 1)))
                    matrix_array = np.vstack(
                        (charge_conservation_row, matrix_array)
                    )

                rate_matrices.loc[atomic_number, shell] = matrix_array

        return rate_matrices
