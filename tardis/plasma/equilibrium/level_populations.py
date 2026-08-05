import numpy as np
import pandas as pd


class LevelPopulationSolver:
    """Solve normalized populations from level rate matrices."""

    def __init__(self, rates_matrices: pd.DataFrame, levels: pd.DataFrame):
        """Solve the normalized level population values from the rate matrices.

        Parameters
        ----------
        rates_matrices : pd.DataFrame
            DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        levels : pd.DataFrame
            DataFrame of energy levels.
        """
        self.rates_matrices = rates_matrices
        self.levels = levels

    def solve(self) -> pd.DataFrame:
        """Solve every level block independently for every cell.

        Each matrix is solved with a unit population normalization. The
        species and cell labels are retained in the returned DataFrame; no
        population is mixed between ions or cells.

        Returns
        -------
        pd.DataFrame
            Normalized level population values indexed by atomic number, ion
            number and level number. Columns are cells.
        """
        normalized_level_populations = pd.DataFrame(
            index=self.levels.index,
            columns=self.rates_matrices.columns,
            dtype=np.float64,
        )

        for species_id in self.rates_matrices.index:
            solved_matrices = []
            for cell in self.rates_matrices.columns:
                normalized_rate_matrix = self.rates_matrices.loc[species_id, cell]
                normalization = np.zeros(normalized_rate_matrix.shape[0])
                normalization[0] = 1.0
                population = np.linalg.solve(
                    normalized_rate_matrix,
                    normalization,
                )
                solved_matrices.append(population)
            normalized_level_populations.loc[species_id, :] = np.vstack(
                solved_matrices
            ).T
        return normalized_level_populations

    def solve_raw(
        self, raw_rate_matrices: pd.DataFrame
    ) -> pd.DataFrame:
        """Normalize and solve raw level blocks from ``LevelRateMatrix``.

        ``LevelRateMatrix.build_raw_rate_matrices`` returns generators whose
        columns sum to zero. This method copies those matrices, replaces the
        first equilibrium row with ones, and delegates the linear solves to
        :meth:`solve`. The input DataFrame and its matrix objects are not
        modified.

        Parameters
        ----------
        raw_rate_matrices : pandas.DataFrame
            Unnormalized level rate matrices indexed by species and with one
            matrix per cell.

        Returns
        -------
        pandas.DataFrame
            Element-normalized level fractions, indexed like ``levels``.
        """
        normalized_rate_matrices = raw_rate_matrices.copy(deep=True)
        for species_id in raw_rate_matrices.index:
            for cell in raw_rate_matrices.columns:
                normalized_rate_matrix = np.asarray(
                    raw_rate_matrices.loc[species_id, cell], dtype=float
                ).copy()
                # Replace one dependent equilibrium equation with normalization
                # without mutating the caller's raw generator.
                normalized_rate_matrix[0, :] = 1.0
                normalized_rate_matrices.loc[
                    species_id, cell
                ] = normalized_rate_matrix
        solver = LevelPopulationSolver(
            normalized_rate_matrices,
            self.levels,
        )
        return solver.solve()
