import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix, RateMatrix
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def _assert_physical_rate_matrix(matrix, conservation_row):
    assert matrix.ndim == 2
    assert np.all(matrix[conservation_row] == 1.0)

    physical_matrix = matrix.copy()
    physical_matrix[conservation_row, :] = 0.0
    np.fill_diagonal(physical_matrix, 0.0)
    assert np.all(physical_matrix >= 0.0)

    diagonal = np.diag(matrix)
    diagonal = np.delete(diagonal, conservation_row)
    assert np.all(diagonal <= 0.0)


def test_bound_bound_rate_matrix_has_conservation_rows_and_physical_rates(
    new_chianti_atomic_dataset_si,
    rate_solver_list,
    collisional_simulation_state,
):
    rate_matrix_solver = RateMatrix(
        rate_solver_list, new_chianti_atomic_dataset_si.levels
    )
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0,
        collisional_simulation_state.t_radiative,
        1e6 * u.g / u.cm**3,
    )

    matrices = rate_matrix_solver.solve(rad_field, electron_dist)

    assert matrices.index.names == ["atomic_number", "ion_number"]
    assert matrices.columns.equals(
        pd.Index(range(len(collisional_simulation_state.t_radiative)))
    )
    for matrix in matrices.to_numpy().flat:
        _assert_physical_rate_matrix(matrix, conservation_row=0)


def test_bound_bound_rate_matrix_solves_normalized_balance_equations(
    new_chianti_atomic_dataset_si,
    rate_solver_list,
    collisional_simulation_state,
):
    rate_matrix_solver = RateMatrix(
        rate_solver_list, new_chianti_atomic_dataset_si.levels
    )
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0,
        collisional_simulation_state.t_radiative,
        1e6 * u.g / u.cm**3,
    )
    matrices = rate_matrix_solver.solve(rad_field, electron_dist)

    for matrix in matrices.to_numpy().flat:
        right_hand_side = np.zeros(matrix.shape[0])
        right_hand_side[0] = 1.0
        population = np.linalg.solve(matrix, right_hand_side)
        # The null-space residual is roundoff-sized for the normalized solve.
        npt.assert_allclose(
            matrix @ population, right_hand_side, rtol=1e-12, atol=1e-15
        )
        npt.assert_allclose(population.sum(), 1.0, rtol=1e-12)
        assert np.all(population >= 0.0)


def test_ion_rate_matrix_has_charge_and_normalization_rows(
    photoionization_rate_solver,
    collisional_ionization_rate_solver,
    collisional_simulation_state,
    mock_boltzmann_factor,
):
    rate_matrix_solver = IonRateMatrix(
        photoionization_rate_solver, collisional_ionization_rate_solver
    )
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0,
        collisional_simulation_state.t_radiative,
        1e6 * u.g / u.cm**3,
    )
    lte_level_population = pd.DataFrame(
        data=np.ones((2, 20)) * 1e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )
    lte_ion_population = pd.DataFrame(
        data=np.ones((2, 20)) * 1e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
    )

    matrices = rate_matrix_solver.solve(
        rad_field,
        electron_dist,
        lte_level_population,
        lte_level_population,
        lte_ion_population,
        lte_ion_population,
        1.0,
        mock_boltzmann_factor,
        charge_conservation=True,
    )

    assert matrices.index.names == ["atomic_number"]
    assert matrices.columns.equals(
        pd.Index(range(len(collisional_simulation_state.t_radiative)))
    )
    for matrix in matrices.to_numpy().flat:
        ion_states = matrix.shape[1] - 1
        expected_charge_row = np.hstack(
            (np.arange(ion_states), -1.0)
        )
        npt.assert_array_equal(matrix[0], expected_charge_row)
        # The extra column is the electron-density unknown and is zero in the
        # elemental normalization row.
        npt.assert_array_equal(matrix[-1], np.hstack((np.ones(ion_states), 0)))
        assert matrix.shape == (ion_states + 1, ion_states + 1)


def test_rate_matrix_solver(
    new_chianti_atomic_dataset_si,
    rate_solver_list,
    collisional_simulation_state,
    regression_data,
):
    rate_matrix_solver = RateMatrix(
        rate_solver_list, new_chianti_atomic_dataset_si.levels
    )

    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
    )

    actual = rate_matrix_solver.solve(rad_field, electron_dist)

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-15)


@pytest.mark.parametrize("charge_conservation", [True, False])
def test_ion_rate_matrix_solver(
    photoionization_rate_solver,
    collisional_ionization_rate_solver,
    collisional_simulation_state,
    mock_boltzmann_factor,
    charge_conservation,
    regression_data,
):
    rate_matrix_solver = IonRateMatrix(
        photoionization_rate_solver, collisional_ionization_rate_solver
    )

    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
    )

    lte_level_population = pd.DataFrame(
        data=np.ones((2, 20)) * 1e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )

    lte_ion_population = pd.DataFrame(
        data=np.ones((2, 20)) * 1e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
    )

    level_population = lte_level_population.copy() * 1.4
    ion_population = lte_ion_population.copy() * 3.0

    actual = rate_matrix_solver.solve(
        rad_field,
        electron_dist,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        1.0,
        mock_boltzmann_factor,
        charge_conservation,
    )

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-15)
