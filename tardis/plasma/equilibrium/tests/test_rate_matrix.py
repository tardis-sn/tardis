import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.io.atom_data import AtomData
from tardis.model.base import SimulationState
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix, RateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_bound_bound_rate_matrix_has_conservation_rows_and_physical_rates(
    new_chianti_atomic_dataset_si: AtomData,
    rate_solver_list: list[tuple[object, str]],
    collisional_simulation_state: SimulationState,
) -> None:
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
        assert matrix.ndim == 2
        assert np.all(matrix[0] == 1.0)

        physical_matrix = matrix.copy()
        physical_matrix[0, :] = 0.0
        np.fill_diagonal(physical_matrix, 0.0)
        assert np.all(physical_matrix >= 0.0)

        diagonal = np.delete(np.diag(matrix), 0)
        assert np.all(diagonal <= 0.0)


def test_bound_bound_rate_matrix_solves_normalized_balance_equations(
    new_chianti_atomic_dataset_si: AtomData,
    rate_solver_list: list[tuple[object, str]],
    collisional_simulation_state: SimulationState,
) -> None:
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


def test_ion_rate_matrix_preserves_electron_density_rate_powers(
    mock_photoionization_cross_sections: pd.DataFrame,
    mock_boltzmann_factor: pd.DataFrame,
) -> None:
    """Ion matrix rates retain their physical electron-density powers."""
    photoionization_cross_sections = pd.concat(
        [
            mock_photoionization_cross_sections,
            mock_photoionization_cross_sections.assign(
                nu=1.1 * mock_photoionization_cross_sections["nu"]
            ),
        ]
    ).sort_values(["atomic_number", "ion_number", "level_number", "nu"])
    rate_matrix_solver = IonRateMatrix(
        AnalyticPhotoionizationRateSolver(photoionization_cross_sections),
        CollisionalIonizationRateSolver(photoionization_cross_sections),
    )

    electron_densities = np.array([1.0e9, 1.0e9, 1.0e9, 1.0e9])
    columns = pd.RangeIndex(len(electron_densities))
    temperatures = np.full(len(columns), 1.0e4) * u.K
    radiation_field = DilutePlanckianRadiationField(
        temperatures, dilution_factor=np.full(len(columns), 0.5)
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        temperatures,
        electron_densities * u.cm**-3,
    )

    level_index = mock_photoionization_cross_sections.index
    lte_level_population = pd.DataFrame(
        1.0e5, index=level_index, columns=columns
    )
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)],
        names=["atomic_number", "ion_number"],
    )
    lte_ion_population = pd.DataFrame(1.0e5, index=ion_index, columns=columns)
    boltzmann_factor = pd.DataFrame(
        np.repeat(
            mock_boltzmann_factor.iloc[:, :1].to_numpy(),
            len(columns),
            axis=1,
        ),
        index=level_index,
        columns=columns,
    )
    level_to_continuum_saha_factor = pd.DataFrame(
        1.0e-15, index=level_index, columns=columns
    )
    partition_function = pd.DataFrame(
        1.0, index=ion_index, columns=columns
    )

    rate_matrices = rate_matrix_solver.solve(
        radiation_field,
        electron_distribution,
        lte_level_population,
        1.4 * lte_level_population,
        lte_ion_population,
        3.0 * lte_ion_population,
        partition_function,
        boltzmann_factor,
        level_to_continuum_saha_factor=level_to_continuum_saha_factor,
    )
    matrices = np.stack(rate_matrices.loc[1].to_numpy())

    assert np.all(np.isfinite(matrices))

    ionization_rates = -matrices[:, 0, 0]
    ionization_steps = np.diff(ionization_rates)
    assert ionization_rates[0] > 0.0
    assert ionization_steps[0] > 0.0
    npt.assert_allclose(
        ionization_steps, ionization_steps[0], rtol=1.0e-9, atol=0.0
    )

    recombination_rates = matrices[:, 0, 1]
    npt.assert_allclose(recombination_rates[0], 0.0, atol=0.0)
    recombination_second_differences = np.diff(recombination_rates, n=2)
    npt.assert_allclose(
        recombination_second_differences,
        recombination_second_differences[0],
        rtol=1.0e-9,
        atol=0.0,
    )
    assert recombination_second_differences[0] > 0.0
    radiative_recombination_at_one_step = (
        recombination_rates[1] - 0.5 * recombination_second_differences[0]
    )
    assert radiative_recombination_at_one_step > 0.0


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
    partition_function = pd.DataFrame(
        1.0,
        index=mock_boltzmann_factor.index,
        columns=mock_boltzmann_factor.columns,
    )

    actual = rate_matrix_solver.solve(
        rad_field,
        electron_dist,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        partition_function,
        mock_boltzmann_factor,
    )
    if charge_conservation:
        for matrix in actual.to_numpy().flat:
            right_hand_side = np.zeros(matrix.shape[0])
            right_hand_side[1] = 1.0
            population = np.linalg.solve(matrix, right_hand_side)
            npt.assert_allclose(
                matrix @ population,
                right_hand_side,
                rtol=1e-12,
                atol=1e-15,
            )
            npt.assert_allclose(population.sum(), 1.0, rtol=1e-12)
            assert np.all(population >= 0.0)
        pdt.assert_index_equal(
            rate_matrix_solver.ion_population_index,
            pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)],
                names=["atomic_number", "ion_number"],
            ),
        )
        return

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-15)
