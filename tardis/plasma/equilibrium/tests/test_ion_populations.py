import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.charge_conservation import (
    ChargeConservationSolver,
)
from tardis.plasma.equilibrium.ion_populations import (
    AnalyticEquilibriumIonPopulationSolver,
)
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix
from tardis.plasma.exceptions import PlasmaIonizationError
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_solve(
    photoionization_rate_solver,
    collisional_ionization_rate_solver,
    regression_data,
):

    radiation_field = DilutePlanckianRadiationField(
        np.ones(20) * 10000 * u.K, dilution_factor=np.ones(20) * 0.5
    )
    thermal_electron_energy_distribution = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 2e9 * u.cm**-3
    )
    lte_level_population = pd.DataFrame(
        data=np.vstack(
            [np.ones(20) * 1e5, np.ones(20) * 1e-1, np.ones(20) * 1e10]
        ),
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )

    lte_ion_population = pd.DataFrame(
        data=np.vstack([np.ones(20) * 1e5, np.ones(20) * 1e10]),
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
    )

    boltzmann_factor = pd.DataFrame(
        data=np.vstack(
            [np.ones(20) * 2.0, np.ones(20) * 0.000011, np.ones(20)]
        ),
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )

    elemental_number_density = pd.DataFrame(
        data=np.vstack([np.ones(20) * 1e5]),
        index=pd.Index([1], name="atomic_number"),
    )

    level_population = lte_level_population.copy() * 1.4
    ion_population = lte_ion_population.copy() * 1.1
    charge_conservation = False

    ion_population_solver = AnalyticEquilibriumIonPopulationSolver(
        IonRateMatrix(),
        photoionization_rate_solver,
        collisional_ionization_rate_solver,
        elemental_number_density,
    )

    actual_ion_population, actual_electron_density = (
        ion_population_solver.solve(
            radiation_field,
            thermal_electron_energy_distribution,
            lte_level_population,
            level_population,
            lte_ion_population,
            ion_population,
            1.0,
            boltzmann_factor,
            charge_conservation,
        )
    )

    expected_ion_population = regression_data.sync_dataframe(
        actual_ion_population, key="ion_population"
    )
    expected_electron_density = regression_data.sync_dataframe(
        actual_electron_density, key="electron_density"
    )

    pdt.assert_frame_equal(
        actual_ion_population, expected_ion_population, atol=0, rtol=1e-15
    )
    pdt.assert_series_equal(
        actual_electron_density, expected_electron_density, atol=0, rtol=1e-15
    )

    assert np.all(actual_ion_population.to_numpy() >= 0.0)
    number_density_from_ions = actual_ion_population.groupby(
        level="atomic_number"
    ).sum()
    pdt.assert_index_equal(
        number_density_from_ions.index,
        elemental_number_density.index,
        check_names=False,
    )
    npt.assert_allclose(
        number_density_from_ions.to_numpy(),
        elemental_number_density.to_numpy(),
        rtol=1e-12,
    )

    electron_density_from_ions = (
        actual_ion_population
        * actual_ion_population.index.get_level_values("ion_number").to_numpy()[
            :, None
        ]
    ).sum()
    npt.assert_allclose(
        actual_electron_density.to_numpy(),
        electron_density_from_ions.to_numpy(),
        rtol=1e-12,
    )

    for shell in actual_ion_population.columns:
        matrix = ion_population_solver.rates_matrices.loc[1, shell]
        population = actual_ion_population[shell].to_numpy()
        balance = np.array([0.0, elemental_number_density.loc[1, shell]])
        npt.assert_allclose(
            matrix @ population, balance, rtol=1e-12, atol=1e-12
        )
        assert np.isfinite(np.linalg.cond(matrix))


def test_build_elemental_population_solution_returns_populations() -> None:
    transition_names = [
        "atomic_number",
        "ion_number",
        "ion_number_source",
        "ion_number_destination",
        "level_number_source",
        "level_number_destination",
    ]
    transition_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1, 0, 0)], names=transition_names
    )
    raw_level_rate_matrices = pd.DataFrame(
        [[np.zeros((1, 1))]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0)], names=["atomic_number", "ion_number"]
        ),
        columns=[0],
    )
    photoionization_rates = pd.DataFrame(
        [4.0], index=transition_index, columns=[0]
    )
    recombination_index = pd.MultiIndex.from_tuples(
        [(1, 0, 1, 0, 0, 0)], names=transition_names
    )
    recombination_rates = pd.DataFrame(
        [2.0], index=recombination_index, columns=[0]
    )
    # The synthetic rates 4 and 2 give a 1:2 H I/H II stationary population
    # ratio; abundance 9 therefore produces the exact 3 and 6 populations.
    matrix_set = IonRateMatrix().solve_ion_and_level(
        atomic_number=1,
        raw_level_rate_matrices=raw_level_rate_matrices,
        photoion_rates_df=photoionization_rates,
        recomb_rates_df=recombination_rates,
        ion_stage_count=2,
    )

    solution = (
        AnalyticEquilibriumIonPopulationSolver._build_elemental_population_solution(
            matrix_set, pd.Series([9.0], index=[0])
        )
    )

    npt.assert_allclose(
        solution.normalized_level_populations.loc[(1, 0, 0)].to_numpy(),
        [1.0 / 3.0],
    )
    npt.assert_allclose(
        solution.normalized_ion_populations[0].to_numpy(),
        [1.0 / 3.0, 2.0 / 3.0],
    )
    npt.assert_allclose(
        solution.level_populations.loc[(1, 0, 0)].to_numpy(), [3.0]
    )
    npt.assert_allclose(solution.ion_populations[0].to_numpy(), [3.0, 6.0])
    npt.assert_allclose(solution.normalized_ion_populations[0].sum(), 1.0)


def test_build_elemental_population_solution_reports_singular_matrix() -> None:
    transition_names = [
        "atomic_number",
        "ion_number",
        "ion_number_source",
        "ion_number_destination",
        "level_number_source",
        "level_number_destination",
    ]
    transition_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1, 0, 0)], names=transition_names
    )
    raw_level_rate_matrices = pd.DataFrame(
        [[np.zeros((1, 1))]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0)], names=["atomic_number", "ion_number"]
        ),
        columns=[0],
    )
    # Zero photoionization and recombination rates leave the state disconnected.
    zero_rates = pd.DataFrame([0.0], index=transition_index, columns=[0])
    matrix_set = IonRateMatrix().solve_ion_and_level(
        atomic_number=1,
        raw_level_rate_matrices=raw_level_rate_matrices,
        photoion_rates_df=zero_rates,
        ion_stage_count=2,
    )

    with pytest.raises(np.linalg.LinAlgError):
        AnalyticEquilibriumIonPopulationSolver._build_elemental_population_solution(
            matrix_set, pd.Series([1.0], index=[0])
        )


def test_analytic_charge_conservation_uses_elemental_rate_matrices() -> None:
    transition_names = [
        "atomic_number",
        "ion_number",
        "ion_number_source",
        "ion_number_destination",
        "level_number_source",
        "level_number_destination",
    ]
    photoionization_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1, 0, 0)], names=transition_names
    )
    recombination_index = pd.MultiIndex.from_tuples(
        [(1, 0, 1, 0, 0, 0)], names=transition_names
    )

    # The synthetic photoionization rate 3 and recombination factor 2 give the
    # charge-conservation solve a closed-form quadratic reference solution.
    class PhotoionizationRates:
        def solve(
            self,
            _radiation_field: object,
            electron_distribution: ThermalElectronEnergyDistribution,
            *rate_args: object,
        ) -> tuple[pd.DataFrame, pd.DataFrame]:
            del rate_args
            recombination_rate = 2.0 * electron_distribution.number_density.value[
                0
            ]
            return (
                pd.DataFrame(
                    [3.0], index=photoionization_index, columns=[0]
                ),
                pd.DataFrame(
                    [recombination_rate],
                    index=recombination_index,
                    columns=[0],
                ),
            )

    class CollisionalRates:
        def solve(
            self, *rate_args: object
        ) -> tuple[pd.DataFrame, pd.DataFrame]:
            del rate_args
            return (
                pd.DataFrame([0.0], index=photoionization_index, columns=[0]),
                pd.DataFrame([0.0], index=recombination_index, columns=[0]),
            )

    lte_level_population = pd.DataFrame(
        [[1.0]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
        columns=[0],
    )
    lte_ion_population = pd.DataFrame(
        [[1.0], [1.0]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
        columns=[0],
    )
    estimated_ion_population = pd.DataFrame(
        [[1.0], [1.0], [0.5], [0.5]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1), (6, 0), (6, 6)],
            names=["atomic_number", "ion_number"],
        ),
        columns=[0],
    )
    solver = AnalyticEquilibriumIonPopulationSolver(
        IonRateMatrix(),
        PhotoionizationRates(),
        CollisionalRates(),
        pd.DataFrame(
            [[10.0], [1.0]],
            index=pd.Index([1, 6], name="atomic_number"),
            columns=[0],
        ),
    )

    ion_populations, electron_number_density = solver.solve(
        None,
        ThermalElectronEnergyDistribution(
            0 * u.erg,
            np.array([10000.0]) * u.K,
            np.array([1.0]) / u.cm**3,
        ),
        lte_level_population,
        lte_level_population,
        lte_ion_population,
        estimated_ion_population,
        pd.DataFrame([[1.0]], index=lte_level_population.index, columns=[0]),
        pd.DataFrame([[1.0]], index=lte_level_population.index, columns=[0]),
        charge_conservation=True,
        tolerance=1e-8,  # Looser than production to isolate the analytic root.
    )
    expected = (3.0 + np.sqrt(3.0**2 + 4 * 2.0 * 39.0)) / (2 * 2.0)

    npt.assert_allclose(electron_number_density, [expected])
    npt.assert_allclose(
        ion_populations.loc[(1, 1)].to_numpy(), [expected - 3.0]
    )


def test_charge_conservation_solver_solves_pure_hydrogen() -> None:
    # The abundance 10 and synthetic alpha/gamma values yield a simple
    # quadratic charge-closure root independent of atomic-data fixtures.
    elemental_number_density = pd.DataFrame(
        [[10.0]], index=pd.Index([1], name="atomic_number"), columns=[0]
    )
    alpha = 2.0
    gamma = 3.0

    def solve_trial(electron_density: pd.Series) -> pd.DataFrame:
        denominator = gamma + alpha * electron_density.iloc[0]
        ionized = gamma * 10.0 / denominator
        return pd.DataFrame(
            [[10.0 - ionized], [ionized]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)],
                names=["atomic_number", "ion_number"],
            ),
            columns=[0],
        )

    electron_number_density, ion_populations = ChargeConservationSolver(
        elemental_number_density, solve_trial
    ).solve()
    expected = (-gamma + np.sqrt(gamma**2 + 4 * alpha * gamma * 10.0)) / (
        2 * alpha
    )

    npt.assert_allclose(electron_number_density, [expected])
    npt.assert_allclose(
        ion_populations.loc[(1, 1)].to_numpy(), [expected]
    )


def test_charge_conservation_solver_includes_all_carbon_stages() -> None:
    elemental_number_density = pd.DataFrame(
        [[1.0], [1.0]],
        index=pd.Index([1, 6], name="atomic_number"),
        columns=[0],
    )
    carbon_populations = np.array([0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5])

    def solve_trial(electron_density: pd.Series) -> pd.DataFrame:
        # Carbon contributes 4.5 electrons; the 0.1 hydrogen response then
        # closes at the analytic electron density 5.0.
        hydrogen_ionized = 0.1 * electron_density.iloc[0]
        hydrogen = np.array([1.0 - hydrogen_ionized, hydrogen_ionized])
        rows = [(1, 0), (1, 1)] + [(6, ion_number) for ion_number in range(7)]
        return pd.DataFrame(
            np.concatenate((hydrogen, carbon_populations))[:, None],
            index=pd.MultiIndex.from_tuples(
                rows, names=["atomic_number", "ion_number"]
            ),
            columns=[0],
        )

    electron_number_density, _ = ChargeConservationSolver(
        elemental_number_density, solve_trial
    ).solve()

    npt.assert_allclose(electron_number_density, [5.0])


@pytest.mark.parametrize(
    ("ionized_at_zero", "ionized_at_max", "expected_electron_density"),
    [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)],
)
def test_charge_conservation_solver_accepts_physical_endpoints(
    ionized_at_zero: float,
    ionized_at_max: float,
    expected_electron_density: float,
) -> None:
    elemental_number_density = pd.DataFrame(
        [[1.0]], index=pd.Index([1], name="atomic_number"), columns=[0]
    )

    def solve_trial(electron_density: pd.Series) -> pd.DataFrame:
        # The 0 and 1 ion fractions place the exact root at each end of the
        # physical bracket, without calling the interior root finder.
        ionized = (
            ionized_at_zero
            if electron_density.iloc[0] == 0
            else ionized_at_max
        )
        return pd.DataFrame(
            [[1.0 - ionized], [ionized]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)],
                names=["atomic_number", "ion_number"],
            ),
            columns=[0],
        )

    electron_number_density, _ = ChargeConservationSolver(
        elemental_number_density, solve_trial
    ).solve()

    npt.assert_allclose(electron_number_density, [expected_electron_density])


def test_charge_conservation_solver_reports_nonfinite_residual() -> None:
    elemental_number_density = pd.DataFrame(
        [[1.0]], index=pd.Index([1], name="atomic_number"), columns=[0]
    )

    def solve_trial(_: pd.Series) -> pd.DataFrame:
        return pd.DataFrame(
            [[0.0], [1.0]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, np.inf)],
                names=["atomic_number", "ion_number"],
            ),
            columns=[0],
        )

    with pytest.raises(PlasmaIonizationError, match="nonfinite"):
        ChargeConservationSolver(
            elemental_number_density, solve_trial
        ).solve()


def test_charge_conservation_solver_reports_missing_bracket() -> None:
    elemental_number_density = pd.DataFrame(
        [[1.0]], index=pd.Index([1], name="atomic_number"), columns=[0]
    )

    def solve_trial(_: pd.Series) -> pd.DataFrame:
        return pd.DataFrame(
            [[0.0], [1.0]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 2)],
                names=["atomic_number", "ion_number"],
            ),
            columns=[0],
        )

    with pytest.raises(PlasmaIonizationError, match="does not bracket"):
        ChargeConservationSolver(
            elemental_number_density, solve_trial
        ).solve()
