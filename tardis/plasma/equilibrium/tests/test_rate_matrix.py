import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.io.atom_data import AtomData
from tardis.model.base import SimulationState
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rate_matrix import (
    IonRateMatrix,
    LevelRateMatrix,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_level_resolved_recombination_rate_dataframe,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_recombination_returns_population_to_its_bound_level() -> None:
    """Level-resolved recombination must populate its originating bound level."""
    rate = pd.DataFrame(
        [[2.0]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 7)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
        columns=[0],
    )

    indexed_rate = reindex_level_resolved_recombination_rate_dataframe(rate)

    assert indexed_rate.index[0] == (1, 0, 1, 0, 0, 7)


def test_level_rate_matrix_exposes_raw_rate_matrix_for_selected_ion() -> None:
    levels = pd.DataFrame(
        {"energy": [0.0, 1.0]},
        index=pd.MultiIndex.from_tuples(
            [(6, 1, 0), (6, 1, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )
    transition_index = pd.MultiIndex.from_tuples(
        [
            (6, 1, 1, 0, 0, 1),
            (6, 1, 1, 0, 1, 0),
        ],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ],
    )
    radiative_rates = pd.DataFrame(
        [2.0, 3.0], index=transition_index, columns=[0]
    )
    collisional_rates = pd.DataFrame(0.0, index=transition_index, columns=[0])

    raw_rate_matrices = LevelRateMatrix(
        levels
    ).build_raw_level_matrices_for_ion(
        radiative_rates,
        collisional_rates,
        np.array([1.0]),
        (6, 1),
    )
    raw_rate_matrix = raw_rate_matrices.loc[(6, 1), 0]

    # The distinct synthetic rates 2 and 3 expose the source/destination
    # orientation: columns sum to zero in this unnormalized generator.
    np.testing.assert_array_equal(raw_rate_matrix, [[-2.0, 3.0], [2.0, -3.0]])
    np.testing.assert_allclose(raw_rate_matrix.sum(axis=0), 0.0)


def test_elemental_level_ion_rate_matrix_set_retains_carbon_ion_states() -> (
    None
):
    levels = pd.DataFrame(
        {"energy": [0.0, 1.0]},
        index=pd.MultiIndex.from_tuples(
            [(6, 0, 0), (6, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )
    bound_bound_index = pd.MultiIndex.from_tuples(
        [
            (6, 0, 0, 0, 0, 1),
            (6, 0, 0, 0, 1, 0),
        ],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ],
    )
    raw_level_rate_matrices = LevelRateMatrix(levels).build_raw_rate_matrices(
        pd.DataFrame([2.0, 3.0], index=bound_bound_index, columns=[0]),
        pd.DataFrame(0.0, index=bound_bound_index, columns=[0]),
        np.array([1.0]),
    )
    # The unequal 5 and 7 bound-free rates make both off-diagonal directions
    # and their state-index positions independently observable.
    bound_free_index = pd.MultiIndex.from_tuples(
        [(6, 0, 0, 1, 0, 0)],
        names=bound_bound_index.names,
    )
    photoion_rates = pd.DataFrame([5.0], index=bound_free_index, columns=[0])
    recombination_index = pd.MultiIndex.from_tuples(
        [(6, 0, 1, 0, 0, 0)],
        names=bound_bound_index.names,
    )
    recombination_rates = pd.DataFrame(
        [7.0], index=recombination_index, columns=[0]
    )
    zero_rates = pd.DataFrame(0.0, index=bound_free_index, columns=[0])

    matrix_set = IonRateMatrix().solve_ion_and_level(
        atomic_number=6,
        raw_level_rate_matrices=raw_level_rate_matrices,
        photoion_rates_df=photoion_rates,
        recomb_rates_df=recombination_rates,
        collisional_ionization_rates_df=zero_rates,
        collision_recombination_rates_df=zero_rates,
        ion_stage_count=7,
    )

    # The state vector has two explicit C I levels, then six ion-stage totals,
    # ending with the bare C 6+ charge state that has no level block.
    assert matrix_set.state_index.level_positions == {(0, 0): 0, (0, 1): 1}
    assert matrix_set.state_index.ion_positions == {
        1: 2,
        2: 3,
        3: 4,
        4: 5,
        5: 6,
        6: 7,
    }
    normalized_rate_matrix = matrix_set.normalized_rate_matrices.loc[6, 0]
    assert normalized_rate_matrix.shape == (8, 8)
    np.testing.assert_array_equal(normalized_rate_matrix[0], np.ones(8))
    raw_rate_matrix = matrix_set.raw_elemental_rate_matrices.loc[6, 0]
    np.testing.assert_allclose(raw_rate_matrix.sum(axis=0), 0.0)
    assert raw_rate_matrix[2, 0] == 5.0
    assert raw_rate_matrix[0, 2] == 7.0


def test_elemental_matrix_accepts_scale_accurate_level_conservation() -> None:
    """Large physical rates may conserve population within roundoff, not absolute zero."""
    large_rate = 1e9
    small_rate = 1e-6
    raw_level_matrix = np.array(
        [
            [-(large_rate + small_rate), 1.0, 1.0],
            [large_rate, -1.0, 0.0],
            [small_rate, 0.0, -1.0],
        ]
    )
    raw_level_rate_matrices = pd.DataFrame(
        [[raw_level_matrix]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0)], names=["atomic_number", "ion_number"]
        ),
        columns=[0],
    )
    transition_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1, 0, 0)],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ],
    )
    zero_bound_free_rate = pd.DataFrame(
        [0.0], index=transition_index, columns=[0]
    )

    matrix_set = IonRateMatrix().solve_ion_and_level(
        atomic_number=1,
        raw_level_rate_matrices=raw_level_rate_matrices,
        photoion_rates_df=zero_bound_free_rate,
        ion_stage_count=2,
    )

    np.testing.assert_array_equal(
        matrix_set.raw_elemental_rate_matrices.loc[1, 0][:3, :3],
        raw_level_matrix,
    )


def test_bound_bound_rate_matrix_has_conservation_rows_and_physical_rates(
    new_chianti_atomic_dataset_si: AtomData,
    rate_solver_list: list[tuple[object, str]],
    collisional_simulation_state: SimulationState,
) -> None:
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0,
        collisional_simulation_state.t_radiative,
        1e6 * u.g / u.cm**3,
    )

    radiative_rates = rate_solver_list[0][0].solve(rad_field)
    collisional_rates = rate_solver_list[1][0].solve(electron_dist.temperature)
    matrices = LevelRateMatrix(new_chianti_atomic_dataset_si.levels).solve(
        radiative_rates,
        collisional_rates,
        electron_dist.number_density.value,
    )

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
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0,
        collisional_simulation_state.t_radiative,
        1e6 * u.g / u.cm**3,
    )

    radiative_rates = rate_solver_list[0][0].solve(rad_field)
    collisional_rates = rate_solver_list[1][0].solve(electron_dist.temperature)
    matrices = LevelRateMatrix(new_chianti_atomic_dataset_si.levels).solve(
        radiative_rates,
        collisional_rates,
        electron_dist.number_density.value,
    )

    for matrix in matrices.to_numpy().flat:
        right_hand_side = np.zeros(matrix.shape[0])
        right_hand_side[0] = 1.0
        population = np.linalg.solve(matrix, right_hand_side)
        np.testing.assert_allclose(
            matrix @ population,
            right_hand_side,
            rtol=1e-12,
            atol=1e-15,
        )
        np.testing.assert_allclose(population.sum(), 1.0, rtol=1e-12)
        assert np.all(population >= 0.0)


def test_ion_rate_matrix_has_charge_and_normalization_rows(
    photoionization_rate_solver: AnalyticPhotoionizationRateSolver,
    collisional_ionization_rate_solver: CollisionalIonizationRateSolver,
    collisional_simulation_state: SimulationState,
    mock_boltzmann_factor: pd.DataFrame,
) -> None:
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
    level_population = lte_level_population * 1.4
    ion_population = lte_ion_population * 3.0

    photoion_rates, recombination_rates = photoionization_rate_solver.solve(
        rad_field,
        electron_dist,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        1.0,
        mock_boltzmann_factor,
    )
    level_to_ion_population_factor = lte_level_population / (
        lte_ion_population.values * electron_dist.number_density.value
    )
    collisional_ionization_rates, collisional_recombination_rates = (
        collisional_ionization_rate_solver.solve(
            electron_dist,
            level_to_ion_population_factor,
            1.0,
            mock_boltzmann_factor,
        )
    )
    matrices = IonRateMatrix().solve(
        photoion_rates,
        recombination_rates,
        collisional_ionization_rates,
        collisional_recombination_rates,
        charge_conservation=True,
    )

    assert matrices.index.names == ["atomic_number"]
    assert matrices.columns.equals(
        pd.Index(range(len(collisional_simulation_state.t_radiative)))
    )
    for matrix in matrices.to_numpy().flat:
        ion_states = matrix.shape[1] - 1
        expected_charge_row = np.hstack((np.arange(ion_states), -1.0))
        np.testing.assert_array_equal(matrix[0], expected_charge_row)
        np.testing.assert_array_equal(
            matrix[-1], np.hstack((np.ones(ion_states), 0))
        )
        assert matrix.shape == (ion_states + 1, ion_states + 1)


def test_rate_matrix_solver(
    new_chianti_atomic_dataset_si,
    rate_solver_list,
    collisional_simulation_state,
    regression_data,
):
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
    )

    radiative_rates = rate_solver_list[0][0].solve(rad_field)
    collisional_rates = rate_solver_list[1][0].solve(electron_dist.temperature)
    actual = LevelRateMatrix(new_chianti_atomic_dataset_si.levels).solve(
        radiative_rates,
        collisional_rates,
        electron_dist.number_density.value,
    )

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

    photoion_rates, recombination_rates = photoionization_rate_solver.solve(
        rad_field,
        electron_dist,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        1.0,
        mock_boltzmann_factor,
    )
    level_to_ion_population_factor = lte_level_population / (
        lte_ion_population.values * electron_dist.number_density.value
    )
    collisional_ionization_rates, collisional_recombination_rates = (
        collisional_ionization_rate_solver.solve(
            electron_dist,
            level_to_ion_population_factor,
            1.0,
            mock_boltzmann_factor,
        )
    )
    actual = IonRateMatrix().solve(
        photoion_rates,
        recombination_rates,
        collisional_ionization_rates,
        collisional_recombination_rates,
        charge_conservation,
    )

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-15)
