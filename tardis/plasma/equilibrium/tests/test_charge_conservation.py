from collections.abc import Callable

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
from tardis.plasma.equilibrium.continuum import ContinuumPopulationSolver
from tardis.plasma.equilibrium.ion_populations import (
    ChargeConservingAnalyticIonPopulationSolver,
)
from tardis.plasma.equilibrium.population_state import PopulationState
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix
from tardis.plasma.exceptions import PlasmaIonizationError


class CallbackPopulationSolver:
    """Adapt a deterministic ion calculation to the charge-state contract."""

    def __init__(self, solve_ions: Callable[[pd.Series], pd.DataFrame]) -> None:
        self.solve_ions = solve_ions

    def solve_charge_state_at_electron_density(
        self,
        electron_number_density: pd.Series,
    ) -> PopulationState:
        electron_number_density = electron_number_density.copy()
        ion_number_density = self.solve_ions(electron_number_density)
        level_index = pd.MultiIndex.from_tuples(
            [], names=["atomic_number", "ion_number", "level_number"]
        )
        empty_levels = pd.DataFrame(
            index=level_index,
            columns=electron_number_density.index,
            dtype=float,
        )
        return PopulationState(
            electron_densities=electron_number_density,
            elemental_populations={},
            ion_number_density=ion_number_density,
            level_number_density=empty_levels,
            level_boltzmann_factor=empty_levels,
        )


def make_nebular_population_solver(
    number_density: pd.DataFrame,
    phi: pd.DataFrame,
) -> ContinuumPopulationSolver:
    """Construct a fixture-free nebular charge-state population solver."""
    columns = number_density.columns
    empty_level_index = pd.MultiIndex.from_tuples(
        [], names=["atomic_number", "ion_number", "level_number"]
    )
    empty_level_factor = pd.DataFrame(
        index=empty_level_index, columns=columns, dtype=float
    )
    empty_partition_function = pd.DataFrame(columns=columns, dtype=float)
    return ContinuumPopulationSolver(
        atomic_data=None,
        lines=pd.DataFrame(columns=columns),
        continuum_species=(),
        radiation_field=None,
        electron_temperatures=np.full(len(columns), 10_000.0),
        elemental_number_density=number_density,
        general_level_boltzmann_factor=empty_level_factor,
        general_partition_function=empty_partition_function,
        thermal_g_electron=np.full(len(columns), 1.0),
        beta_electron=np.full(len(columns), 1.0),
        thermal_lte_level_boltzmann_factor=empty_level_factor,
        thermal_lte_partition_function=empty_partition_function,
        ionization_data=pd.DataFrame(columns=columns),
        nebular_phi=phi,
        photoionization_rate_estimator=None,
        stimulated_recombination_rate_estimator=None,
    )


def test_charge_state_recomputes_nebular_carbon() -> None:
    """Verify carbon ion fractions respond to the trial electron density."""
    columns = pd.Index([0], name="shell")
    number_density = pd.DataFrame(
        [[127.0]],
        index=pd.Index([6], name="atomic_number"),
        columns=columns,
    )
    phi = pd.DataFrame(
        2.0,
        index=pd.MultiIndex.from_tuples(
            [(6, ion_number) for ion_number in range(1, 7)],
            names=["atomic_number", "ion_number"],
        ),
        columns=columns,
    )
    solver = make_nebular_population_solver(number_density, phi)

    low_density_state = solver.solve_charge_state_at_electron_density(
        pd.Series([1.0], index=columns)
    )
    high_density_state = solver.solve_charge_state_at_electron_density(
        pd.Series([2.0], index=columns)
    )

    npt.assert_allclose(
        low_density_state.ion_number_density[0],
        [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0],
    )
    npt.assert_allclose(
        high_density_state.ion_number_density[0], np.full(7, 127.0 / 7.0)
    )
    assert set(low_density_state.elemental_populations) == {6}
    carbon_solution = low_density_state.elemental_populations[6]
    pdt.assert_frame_equal(
        carbon_solution.ion_populations,
        low_density_state.ion_number_density,
        check_column_type=False,
    )


def test_continuum_charge_state_rejects_zero_abundance() -> None:
    """Require every charge-solved element to have positive abundance."""
    columns = pd.Index([0], name="shell")
    number_density = pd.DataFrame(
        [[0.0], [1.0]],
        index=pd.Index([1, 6], name="atomic_number"),
        columns=columns,
    )
    phi = pd.DataFrame(
        1.0,
        index=pd.MultiIndex.from_tuples(
            [(1, 1), *[(6, ion_number) for ion_number in range(1, 7)]],
            names=["atomic_number", "ion_number"],
        ),
        columns=columns,
    )
    population_solver = make_nebular_population_solver(number_density, phi)

    with pytest.raises(ValueError, match="positive abundance"):
        population_solver.solve_charge_state_at_electron_density(
            pd.Series([1.0], index=columns)
        )


def test_charge_conservation_returns_final_nebular_carbon_state() -> None:
    """Recover the charge-conserving uniform carbon ion distribution."""
    columns = pd.Index([0], name="shell")
    number_density = pd.DataFrame(
        [[1.0]],
        index=pd.Index([6], name="atomic_number"),
        columns=columns,
    )
    phi = pd.DataFrame(
        3.0,
        index=pd.MultiIndex.from_tuples(
            [(6, ion_number) for ion_number in range(1, 7)],
            names=["atomic_number", "ion_number"],
        ),
        columns=columns,
    )
    population_solver = make_nebular_population_solver(number_density, phi)

    state = ChargeConservationSolver(number_density, population_solver).solve()

    npt.assert_allclose(state.electron_densities, [3.0])
    npt.assert_allclose(state.ion_number_density[0], np.full(7, 1.0 / 7.0))


def test_nebular_zero_phi_uses_iip_input_regularization() -> None:
    """Preserve finite charge closure when an ionization factor vanishes."""
    columns = pd.Index([0, 1], name="shell")
    number_density = pd.DataFrame(
        [[1.0, 1.0]],
        index=pd.Index([6], name="atomic_number"),
        columns=columns,
    )
    phi = pd.DataFrame(
        [[1.0, 0.25], [0.0, 1.0], *[[1.0, 1.0]] * 4],
        index=pd.MultiIndex.from_tuples(
            [(6, ion_number) for ion_number in range(1, 7)],
            names=["atomic_number", "ion_number"],
        ),
        columns=columns,
    )
    population_solver = make_nebular_population_solver(number_density, phi)

    rates = population_solver._build_nebular_rates(
        pd.Series([1.0], index=columns[:1])
    )
    assert rates.loc[(6, 1, 1, 2, 0, 0), 0] == pytest.approx(2.5e-11)

    state = ChargeConservationSolver(number_density, population_solver).solve(
        pd.DataFrame(columns=columns)
    )
    npt.assert_allclose(state.electron_densities[0], 0.6180339938, rtol=1e-8)
    charge_density = (
        state.ion_number_density[0]
        * state.ion_number_density.index.get_level_values("ion_number")
    ).sum()
    npt.assert_allclose(charge_density, state.electron_densities[0], rtol=1e-10)


def test_analytic_charge_conservation_uses_charge_state_solver() -> None:
    """Recover the analytic hydrogen ionization-equilibrium charge root."""
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
        """Return electron-density-dependent recombination rates."""

        def solve(
            self,
            _radiation_field: object,
            electron_distribution: ThermalElectronEnergyDistribution,
            *rate_args: object,
        ) -> tuple[pd.DataFrame, pd.DataFrame]:
            del rate_args
            recombination_rate = (
                2.0 * electron_distribution.number_density.value[0]
            )
            return (
                pd.DataFrame([3.0], index=photoionization_index, columns=[0]),
                pd.DataFrame(
                    [recombination_rate],
                    index=recombination_index,
                    columns=[0],
                ),
            )

    class CollisionalRates:
        """Provide the zero collisional contribution for the analytic root."""

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
        [[1.0], [1.0]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
        columns=[0],
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        np.array([10000.0]) * u.K,
        np.array([1.0]) / u.cm**3,
    )
    partition_function = pd.DataFrame(
        [[1.0]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0)], names=["atomic_number", "ion_number"]
        ),
        columns=[0],
    )
    boltzmann_factor = pd.DataFrame(
        [[1.0]], index=lte_level_population.index, columns=[0]
    )
    solver = ChargeConservingAnalyticIonPopulationSolver(
        IonRateMatrix(),
        PhotoionizationRates(),
        CollisionalRates(),
        pd.DataFrame(
            [[10.0]],
            index=pd.Index([1], name="atomic_number"),
            columns=[0],
        ),
        radiation_field=None,
        thermal_electron_energy_distribution=electron_distribution,
        lte_level_population=lte_level_population,
        estimated_level_population=lte_level_population,
        lte_ion_population=lte_ion_population,
        estimated_ion_population=estimated_ion_population,
        partition_function=partition_function,
        boltzmann_factor=boltzmann_factor,
        tolerance=1e-8,  # Looser than production to isolate the analytic root.
    )

    ion_populations, electron_number_density = solver.solve()
    expected = (-3.0 + np.sqrt(3.0**2 + 4 * 2.0 * 3.0 * 10.0)) / (2 * 2.0)

    npt.assert_allclose(electron_number_density, [expected])
    npt.assert_allclose(ion_populations.loc[(1, 1)].to_numpy(), [expected])


def test_charge_conservation_solver_solves_pure_hydrogen() -> None:
    """Solve a one-element charge root with a closed-form expected density."""
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

    solver = ChargeConservationSolver(
        elemental_number_density, CallbackPopulationSolver(solve_trial)
    )
    state = solver.solve(pd.DataFrame(columns=[0]))
    seeded_state = solver.solve(
        pd.DataFrame(columns=[0]),
        pd.Series([9.0], index=[0]),
    )
    expected = (-gamma + np.sqrt(gamma**2 + 4 * alpha * gamma * 10.0)) / (
        2 * alpha
    )

    npt.assert_allclose(state.electron_densities, [expected])
    npt.assert_allclose(
        state.ion_number_density.loc[(1, 1)].to_numpy(), [expected]
    )
    pdt.assert_series_equal(
        seeded_state.electron_densities,
        state.electron_densities,
    )
    pdt.assert_frame_equal(
        seeded_state.ion_number_density,
        state.ion_number_density,
    )


def test_charge_conservation_solver_is_density_scale_independent() -> None:
    """Keep the normalized root solve stable for extremely small densities."""
    elemental_density = 1e-30
    elemental_number_density = pd.DataFrame(
        [[elemental_density]],
        index=pd.Index([1], name="atomic_number"),
        columns=[0],
    )

    def solve_trial(electron_density: pd.Series) -> pd.DataFrame:
        ionized = elemental_density - electron_density.iloc[0]
        return pd.DataFrame(
            [[elemental_density - ionized], [ionized]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)],
                names=["atomic_number", "ion_number"],
            ),
            columns=[0],
        )

    state = ChargeConservationSolver(
        elemental_number_density, CallbackPopulationSolver(solve_trial)
    ).solve(pd.DataFrame(columns=[0]))

    npt.assert_allclose(
        state.electron_densities.to_numpy(),
        [elemental_density / 2.0],
        rtol=1e-12,
        atol=0.0,
    )


def test_charge_conservation_solver_includes_all_carbon_stages() -> None:
    """Include every ion stage, including bare carbon, in charge density."""
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

    state = ChargeConservationSolver(
        elemental_number_density, CallbackPopulationSolver(solve_trial)
    ).solve(pd.DataFrame(columns=[0]))

    npt.assert_allclose(state.electron_densities, [5.0])


@pytest.mark.parametrize(
    ("ionized_at_zero", "ionized_at_max", "expected_electron_density"),
    [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)],
)
def test_charge_conservation_solver_accepts_physical_endpoints(
    ionized_at_zero: float,
    ionized_at_max: float,
    expected_electron_density: float,
) -> None:
    """Accept exact roots at either end of the physical density bracket."""
    elemental_number_density = pd.DataFrame(
        [[1.0]], index=pd.Index([1], name="atomic_number"), columns=[0]
    )

    def solve_trial(electron_density: pd.Series) -> pd.DataFrame:
        # The 0 and 1 ion fractions place the exact root at each end of the
        # physical bracket, without calling the interior root finder.
        ionized = (
            ionized_at_zero if electron_density.iloc[0] == 0 else ionized_at_max
        )
        return pd.DataFrame(
            [[1.0 - ionized], [ionized]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)],
                names=["atomic_number", "ion_number"],
            ),
            columns=[0],
        )

    state = ChargeConservationSolver(
        elemental_number_density, CallbackPopulationSolver(solve_trial)
    ).solve(pd.DataFrame(columns=[0]))

    npt.assert_allclose(state.electron_densities, [expected_electron_density])


def test_charge_conservation_solver_accepts_zero_maximum_density() -> None:
    """Return zero electron density when a shell contains no matter."""
    elemental_number_density = pd.DataFrame(
        [[0.0]], index=pd.Index([1], name="atomic_number"), columns=[0]
    )

    def solve_trial(electron_density: pd.Series) -> pd.DataFrame:
        return pd.DataFrame(
            [[0.0], [0.0]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)],
                names=["atomic_number", "ion_number"],
            ),
            columns=electron_density.index,
        )

    state = ChargeConservationSolver(
        elemental_number_density, CallbackPopulationSolver(solve_trial)
    ).solve(pd.DataFrame(columns=[0]))

    npt.assert_allclose(state.electron_densities, [0.0])


def test_charge_conservation_solver_validates_final_state() -> None:
    """Reject a final back-substitution that no longer conserves charge."""
    columns = pd.Index([0, 1], name="shell")
    elemental_number_density = pd.DataFrame(
        [[1.0, 1.0]],
        index=pd.Index([1], name="atomic_number"),
        columns=columns,
    )

    def solve_trial(electron_density: pd.Series) -> pd.DataFrame:
        ionized = np.zeros(len(electron_density))
        if len(electron_density) > 1:
            ionized.fill(1.0)
        return pd.DataFrame(
            [1.0 - ionized, ionized],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)],
                names=["atomic_number", "ion_number"],
            ),
            columns=electron_density.index,
        )

    with pytest.raises(PlasmaIonizationError, match="exceeds tolerance"):
        ChargeConservationSolver(
            elemental_number_density, CallbackPopulationSolver(solve_trial)
        ).solve(pd.DataFrame(columns=columns))


def test_charge_conservation_solver_reports_nonfinite_residual() -> None:
    """Reject trial ion populations that make the charge residual nonfinite."""
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
            elemental_number_density, CallbackPopulationSolver(solve_trial)
        ).solve(pd.DataFrame(columns=[0]))


def test_charge_conservation_solver_reports_missing_bracket() -> None:
    """Report a physical interval that cannot bracket charge conservation."""
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
            elemental_number_density, CallbackPopulationSolver(solve_trial)
        ).solve(pd.DataFrame(columns=[0]))
