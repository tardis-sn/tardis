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
from tardis.plasma.equilibrium.rate_matrix import (
    AnalyticIonRateMatrix,
    RateMatrix,
    assemble_ion_rate_matrices,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_bound_bound_rate_matrix_has_conservation_rows_and_physical_rates(
    new_chianti_atomic_dataset_si: AtomData,
    rate_solvers: tuple[RadiativeRatesSolver, ThermalCollisionalRateSolver],
    collisional_simulation_state: SimulationState,
) -> None:
    rate_matrix_solver = RateMatrix(
        *rate_solvers, new_chianti_atomic_dataset_si.levels
    )
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0,
        collisional_simulation_state.t_radiative,
        1e6 / u.cm**3,
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
    rate_solvers: tuple[RadiativeRatesSolver, ThermalCollisionalRateSolver],
    collisional_simulation_state: SimulationState,
) -> None:
    rate_matrix_solver = RateMatrix(
        *rate_solvers, new_chianti_atomic_dataset_si.levels
    )
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0,
        collisional_simulation_state.t_radiative,
        1e6 / u.cm**3,
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


def test_prepare_species_rates_preserves_fixed_rate_density_independence(
    new_chianti_atomic_dataset_si: AtomData,
    rate_solvers: tuple[RadiativeRatesSolver, ThermalCollisionalRateSolver],
    collisional_simulation_state: SimulationState,
) -> None:
    """Prepared fixed rates are independent of electron density."""
    rate_matrix_solver = RateMatrix(
        *rate_solvers, new_chianti_atomic_dataset_si.levels
    )
    species = tuple(rate_solvers[0].einstein_coefficients.index[0][:2])
    line_index = (
        rate_solvers[0]
        .einstein_coefficients.xs(
            species,
            level=("atomic_number", "ion_number"),
        )
        .index
    )
    line_index = pd.MultiIndex.from_tuples(
        [(species[0], species[1], *transition) for transition in line_index],
        names=rate_solvers[0].einstein_coefficients.index.names,
    )
    j_blues = pd.DataFrame(
        0.0,
        index=rate_solvers[0].einstein_coefficients.index,
        columns=[0, 1],
    )
    transition = (
        rate_solvers[0]
        .solve(j_blues)
        .xs(
            species,
            level=("atomic_number", "ion_number"),
        )
        .index[0]
    )
    fixed_index = pd.MultiIndex.from_tuples(
        [(species[0], species[1], *transition)],
        names=rate_solvers[0].solve(j_blues).index.names,
    )
    fixed_rates = pd.DataFrame([[3.0, 3.0]], index=fixed_index, columns=[0, 1])
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        collisional_simulation_state.t_radiative[:2],
        np.array([1.0, 1.0e12]) / u.cm**3,
    )

    prepared_rates = rate_matrix_solver.prepare_species_rates(
        j_blues,
        electron_distribution,
        species,
        line_index,
        fixed_rates,
    )

    assert prepared_rates[0].fixed_rate.sum() == 3.0
    assert prepared_rates[1].fixed_rate.sum() == 3.0


def test_assemble_ion_rate_matrices_adds_fixed_ionization() -> None:
    """Fixed ionization contributes to the global ion-stage matrix."""
    transition_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0, 1), (2, 1, 1, 2)],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
        ],
    )
    zero_rates = pd.DataFrame(0.0, index=transition_index, columns=[0])
    fixed_ionization = pd.DataFrame(
        [7.0, 11.0], index=transition_index, columns=[0]
    )

    matrices, _ = assemble_ion_rate_matrices(
        zero_rates,
        zero_rates,
        zero_rates,
        zero_rates,
        fixed_ionization_rates_df=fixed_ionization,
    )

    matrix = matrices.loc[2, 0]
    assert matrix[0, 0] == -7.0
    assert matrix[2, 1] == 11.0
    npt.assert_array_equal(matrix[1], np.ones(3))


def test_assemble_ion_rate_matrices_selects_active_fixed_rate_shell() -> None:
    """Fixed ionization remains local during a single-shell charge solve.

    Claim
    -----
    A rate matrix assembled for one active shell uses only that shell's frozen
    ionization rates.

    Regime
    ------
    Analytic rates contain one shell while the frozen-rate state retains two
    shells, as occurs inside the shell-by-shell charge-conservation root.

    Verification
    ------------
    The selected matrix has exactly the active shell label and contains its
    distinct fixed ionization coefficients, with no contribution from the
    inactive shell.
    """
    transition_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0, 1), (2, 1, 1, 2)],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
        ],
    )
    active_shell = "shell_1"
    zero_rates = pd.DataFrame(
        0.0, index=transition_index, columns=[active_shell]
    )
    fixed_ionization = pd.DataFrame(
        [[7.0, 13.0], [11.0, 17.0]],
        index=transition_index,
        columns=["shell_0", active_shell],
    )

    matrices, _ = assemble_ion_rate_matrices(
        zero_rates,
        zero_rates,
        zero_rates,
        zero_rates,
        fixed_ionization_rates_df=fixed_ionization,
    )

    assert matrices.columns.tolist() == [active_shell]
    matrix = matrices.at[2, active_shell]
    assert matrix[0, 0] == -13.0
    assert matrix[2, 1] == 17.0


def test_rate_matrix_solver(
    new_chianti_atomic_dataset_si,
    rate_solvers,
    collisional_simulation_state,
    regression_data,
):
    rate_matrix_solver = RateMatrix(
        *rate_solvers, new_chianti_atomic_dataset_si.levels
    )

    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 / u.cm**3
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
    rate_matrix_solver = AnalyticIonRateMatrix(
        AnalyticPhotoionizationRateSolver(photoionization_cross_sections),
        CollisionalIonizationRateSolver(photoionization_cross_sections),
    )

    electron_densities = np.array([0, 1e9, 2e9, 3e9])  # test a range of density
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
    helium_ionization_factor = pd.DataFrame(
        [3.0e9, 5.0e9],
        index=pd.MultiIndex.from_tuples(
            [(2, 0), (2, 1)],
            names=["atomic_number", "ion_number"],
        ),
        columns=[columns[0]],
    ).reindex(columns=columns, method="ffill")
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
    partition_function = pd.DataFrame(1.0, index=ion_index, columns=columns)

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
        lte_ionization_factor=helium_ionization_factor,
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

    right_hand_side = np.array([0.0, 1.0, 0.0])
    for shell_idx in range(1, len(columns)):
        helium_population = np.linalg.solve(
            rate_matrices.loc[2, columns[shell_idx]], right_hand_side
        )
        npt.assert_allclose(
            helium_population[1] / helium_population[0],
            helium_ionization_factor.loc[(2, 0), columns[shell_idx]]
            / electron_densities[shell_idx],
        )
        npt.assert_allclose(
            helium_population[2] / helium_population[1],
            helium_ionization_factor.loc[(2, 1), columns[shell_idx]]
            / electron_densities[shell_idx],
        )


def test_analytic_ion_rate_matrix_is_independent_of_shell_batching(
    mock_photoionization_cross_sections: pd.DataFrame,
    mock_boltzmann_factor: pd.DataFrame,
) -> None:
    """An isolated shell retains the ion rates from its batched solution.

    Claim
    -----
    Analytic ion-rate matrices depend on a shell's physical state, not on
    whether that shell is evaluated alone or alongside other shells.

    Regime
    ------
    Two hydrogen shells have distinct temperatures, radiation fields,
    electron densities, and LTE and estimated populations.

    Verification
    ------------
    The shell-1 matrix evaluated in isolation equals shell 1 selected from the
    two-shell solve. This batching invariance directly exercises dataframe
    shell-label alignment in the corrected photoionization rate.
    """
    photoionization_cross_sections = pd.concat(
        [
            mock_photoionization_cross_sections,
            mock_photoionization_cross_sections.assign(
                nu=1.1 * mock_photoionization_cross_sections["nu"]
            ),
        ]
    ).sort_values(["atomic_number", "ion_number", "level_number", "nu"])
    rate_matrix_solver = AnalyticIonRateMatrix(
        AnalyticPhotoionizationRateSolver(photoionization_cross_sections),
        CollisionalIonizationRateSolver(photoionization_cross_sections),
    )

    columns = pd.RangeIndex(2, name="shell")
    temperatures = np.array([9000.0, 12000.0]) * u.K
    dilution_factors = np.array([0.2, 0.7])
    electron_densities = np.array([1.0e7, 3.0e7])
    radiation_field = DilutePlanckianRadiationField(
        temperatures, dilution_factor=dilution_factors
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        temperatures,
        electron_densities * u.cm**-3,
    )
    level_index = mock_photoionization_cross_sections.index
    lte_level_population = pd.DataFrame(
        [[1.0e5, 2.0e5], [3.0e4, 4.0e4]],
        index=level_index,
        columns=columns,
    )
    estimated_level_population = pd.DataFrame(
        [[1.3e5, 2.8e5], [4.5e4, 5.2e4]],
        index=level_index,
        columns=columns,
    )
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)],
        names=["atomic_number", "ion_number"],
    )
    lte_ion_population = pd.DataFrame(
        [[2.0e5, 3.0e5], [8.0e5, 7.0e5]],
        index=ion_index,
        columns=columns,
    )
    estimated_ion_population = pd.DataFrame(
        [[3.0e5, 4.0e5], [7.0e5, 6.0e5]],
        index=ion_index,
        columns=columns,
    )
    partition_function = pd.DataFrame(1.0, index=ion_index, columns=columns)
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

    batched_matrices = rate_matrix_solver.solve(
        radiation_field,
        electron_distribution,
        lte_level_population,
        estimated_level_population,
        lte_ion_population,
        estimated_ion_population,
        partition_function,
        boltzmann_factor,
        level_to_continuum_saha_factor=level_to_continuum_saha_factor,
    )

    isolated_shell = columns[1]
    isolated_matrices = rate_matrix_solver.solve(
        DilutePlanckianRadiationField(
            temperatures[[1]], dilution_factor=dilution_factors[[1]]
        ),
        ThermalElectronEnergyDistribution(
            0.0 * u.erg,
            temperatures[[1]],
            electron_densities[[1]] * u.cm**-3,
        ),
        lte_level_population[[isolated_shell]],
        estimated_level_population[[isolated_shell]],
        lte_ion_population[[isolated_shell]],
        estimated_ion_population[[isolated_shell]],
        partition_function[[isolated_shell]],
        boltzmann_factor[[isolated_shell]],
        level_to_continuum_saha_factor=level_to_continuum_saha_factor[
            [isolated_shell]
        ],
    )

    npt.assert_allclose(
        isolated_matrices.at[1, isolated_shell],
        batched_matrices.at[1, isolated_shell],
        rtol=1.0e-14,
        atol=0.0,
    )


@pytest.mark.parametrize("charge_conservation", [True, False])
def test_ion_rate_matrix_solver(
    photoionization_rate_solver,
    collisional_ionization_rate_solver,
    collisional_simulation_state,
    mock_boltzmann_factor,
    charge_conservation,
    regression_data,
):
    rate_matrix_solver = AnalyticIonRateMatrix(
        photoionization_rate_solver, collisional_ionization_rate_solver
    )

    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 / u.cm**3
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
