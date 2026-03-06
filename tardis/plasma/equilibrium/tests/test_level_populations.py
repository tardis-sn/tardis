import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import RateMatrix
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


class TestLevelPopulationSolver:
    @pytest.fixture(autouse=True)
    def setup(
        self,
        rate_solver_list,
        new_chianti_atomic_dataset_si,
        collisional_simulation_state,
    ):
        rate_matrix_solver = RateMatrix(
            rate_solver_list, new_chianti_atomic_dataset_si.levels
        )

        rad_field = DilutePlanckianRadiationField(
            collisional_simulation_state.t_radiative,
            dilution_factor=np.zeros_like(
                collisional_simulation_state.t_radiative
            ),
        )
        electron_dist = ThermalElectronEnergyDistribution(
            0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
        )

        rates_matrices = rate_matrix_solver.solve(rad_field, electron_dist)
        self.solver = LevelPopulationSolver(
            rates_matrices, new_chianti_atomic_dataset_si.levels
        )

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

    def test_solve(self, regression_data):
        """Test the solve method."""
        result = self.solver.solve()

        # The order of parameters in the auto-generated filename has changed.
        # We manually specify the path to the existing regression file.
        expected_fname = "test_solve__collisional_rate_solver0-radiative_transitions0__.h5"
        expected_fpath = regression_data.fpath.parent / expected_fname
        expected_populations = pd.read_hdf(expected_fpath)
        pdt.assert_frame_equal(result, expected_populations, atol=0, rtol=1e-15)
