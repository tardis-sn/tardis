import numpy as np
import pandas as pd
import pytest

from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver


class TestLevelPopulationSolver:
    @pytest.fixture(autouse=True)
    def setup(self):
        # TODO: improve with regression data
        rates_matrices = pd.DataFrame(
            {
                0: [np.array([[1, 1], [2, -2]])],
            }
        )
        levels = pd.DataFrame({"energy": [0, 1]})
        self.solver = LevelPopulationSolver(rates_matrices, levels)

    def test_calculate_level_population_simple(self):
        """Test solving a 2-level ion."""
        rates_matrix = np.array([[1, 1], [2, -2]])
        expected_population = np.array([0.5, 0.5])
        result = self.solver._LevelPopulationSolver__calculate_level_population(
            rates_matrix
        )
        np.testing.assert_array_almost_equal(result, expected_population)

    def test_calculate_level_population_empty(self):
        """Test empty rate matrix."""
        rates_matrix = np.array([[]])
        with pytest.raises(np.linalg.LinAlgError):
            self.solver._LevelPopulationSolver__calculate_level_population(
                rates_matrix
            )

    def test_calculate_level_population_zeros(self):
        """Test zero rate matrix."""
        rates_matrix = np.array([[0, 0], [0, 0]])
        with pytest.raises(np.linalg.LinAlgError):
            self.solver._LevelPopulationSolver__calculate_level_population(
                rates_matrix
            )

    @pytest.mark.skip(reason="Need regression data")
    def test_solve(self):
        """Test the solve method."""
        expected_populations = False
        result = self.solver.solve()
        pd.testing.assert_frame_equal(result, expected_populations)
