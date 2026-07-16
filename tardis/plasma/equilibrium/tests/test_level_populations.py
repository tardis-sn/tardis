import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import LevelRateMatrix
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


class TestLevelPopulationSolver:
    @pytest.fixture(autouse=True)
    def setup(
        self,
        radiative_rate_solver,
        collisional_rate_solver,
        new_chianti_atomic_dataset_si,
        collisional_simulation_state,
    ):
        rad_field = DilutePlanckianRadiationField(
            collisional_simulation_state.t_radiative,
            dilution_factor=np.zeros_like(
                collisional_simulation_state.t_radiative
            ),
        )
        electron_dist = ThermalElectronEnergyDistribution(
            0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
        )

        radiative_rates = radiative_rate_solver.solve(rad_field)
        collisional_rates = collisional_rate_solver.solve(electron_dist.temperature)
        rates_matrices = LevelRateMatrix(
            new_chianti_atomic_dataset_si.levels
        ).solve(
            radiative_rates,
            collisional_rates,
            electron_dist.number_density.value,
        )
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
        expected_populations = False
        result = self.solver.solve()
        expected_populations = regression_data.sync_dataframe(result)
        pdt.assert_frame_equal(result, expected_populations, atol=0, rtol=1e-15)
