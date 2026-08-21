from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluator,
)
from tardis.plasma.equilibrium.inputs import (
    NumberDensityPerShell,
    SobolevInputs,
)
from tardis.plasma.equilibrium.rate_matrix import RateMatrix
from tardis.plasma.equilibrium.rates import RadiativeRatesSolver
from tardis.transport.montecarlo.estimators import init_estimators_continuum


class ZeroElectronRateSolver:
    """Return no bound-bound electron rates for focused evaluator tests."""

    all_collisional_strengths_index = pd.MultiIndex.from_tuples(
        [],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ],
    )

    def solve(self, temperatures_electron: u.Quantity) -> pd.DataFrame:
        """Return an empty transition-rate frame."""
        return pd.DataFrame(
            index=self.all_collisional_strengths_index,
            columns=pd.RangeIndex(len(temperatures_electron)),
            dtype=np.float64,
        )


@pytest.fixture
def toy_evaluator() -> PlasmaEquilibriumEvaluator:
    """Build the shared one-shell evaluator used by focused tests."""
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    line_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )
    photoionization_cross_sections = pd.DataFrame(
        {"nu": [1.0e15, 1.1e15], "x_sect": [1.0e-18, 2.0e-18]},
        index=level_index,
    )
    levels = pd.DataFrame(
        {"energy": [0.0, 1.0e-12, 0.0], "g": [2.0, 4.0, 1.0]},
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
            names=level_index.names,
        ),
    )
    ionization_data = pd.Series(
        [2.0e-11],
        index=pd.MultiIndex.from_tuples(
            [(1, 1)], names=["atomic_number", "ion_number"]
        ),
    )
    estimators = init_estimators_continuum((2, 1), 1)
    estimators.photo_ion_estimator[:, 0] = [2.0, 3.0]
    estimators.stim_recomb_estimator[:, 0] = [5.0, 7.0]
    rate_matrix_solver = RateMatrix(
        RadiativeRatesSolver(
            pd.DataFrame(
                {"A_ul": [2.0], "B_ul": [0.0], "B_lu": [1.0], "nu": [1.0]},
                index=line_index,
            )
        ),
        ZeroElectronRateSolver(),
        levels,
    )
    population_geometry = NumberDensityPerShell(
        1.0e10, np.array([1.0e10, 0.0]), np.array([0, 1])
    )
    sobolev_inputs = SobolevInputs(
        np.array([0]),
        np.array([1]),
        np.array([2.0]),
        np.array([4.0]),
        np.array([False]),
        np.array([True]),
        np.array([1.0e-20]),
        np.array([0]),
        line_index,
    )
    return PlasmaEquilibriumEvaluator(
        photoionization_cross_sections,
        pd.Series([0, 1], index=level_index),
        estimators,
        2.0e5 * u.s,
        3.0e30 * u.cm**3,
        levels,
        ionization_data,
        rate_matrix_solver,
        pd.DataFrame([1.0], index=line_index),
        (population_geometry,),
        (sobolev_inputs,),
        level_index,
        (1, 0),
        pd.DataFrame([[1.0e10]], index=[1], columns=[0]),
        np.array([1.0e10]),
    )


def test_evaluator_uses_temperature_dependent_continuum_coefficients(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Keep fixed estimator coefficients while rebuilding thermal factors."""
    evaluator = toy_evaluator
    (rates,), _, _, _, *_ = evaluator.calculate_continuum_coefficients(
        np.array([1.0e4])
    )
    (hot_rates,), _, _, _, *_ = (
        evaluator.calculate_continuum_coefficients(np.array([2.0e4]))
    )
    npt.assert_allclose(hot_rates.photoionization, rates.photoionization)
    assert not np.array_equal(
        hot_rates.stimulated_recombination, rates.stimulated_recombination
    )
    assert not np.array_equal(
        hot_rates.collisional_ionization, rates.collisional_ionization
    )


def test_evaluator_solves_reduced_levels_and_rebuilds_sobolev_state(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Evaluate a small fixed candidate without a plasma object."""
    evaluator = toy_evaluator
    estimators = evaluator.estimators_continuum
    level_index = evaluator.level_population_index

    result = evaluator.evaluate(
        np.array([1.0e9]),  # Trial electron density (cm⁻³).
        np.array([1.0e4]),  # Shell electron temperature (K).
        pd.DataFrame([[0.5], [0.5]], index=level_index, columns=[0]),
    )

    npt.assert_allclose(result.normalized_population[0].sum(), 1.0)
    assert np.all(result.normalized_population[0] >= 0.0)
    assert result.diagnostic_ion_ratio.iloc[0] >= 0.0
    npt.assert_allclose(result.beta_sobolev.to_numpy(), [[1.0]])
    npt.assert_allclose(result.trial_level_residual.to_numpy()[0], [0.0])
    npt.assert_allclose(result.level_residual.to_numpy()[0], [0.0])

    original_photoionization_estimator = estimators.photo_ion_estimator.copy()
    candidate_seed = pd.DataFrame(
        [[0.5], [0.5]], index=level_index, columns=[0]
    )
    high_temperature_result = evaluator.evaluate(
        np.array([1.0e9]),  # Trial electron density (cm⁻³).
        np.array([2.0e4]),  # Hotter shell temperature (K).
        candidate_seed,
    )
    assert not np.array_equal(
        result.diagnostic_ion_ratio,
        high_temperature_result.diagnostic_ion_ratio,
    )
    np.testing.assert_array_equal(
        estimators.photo_ion_estimator,
        original_photoionization_estimator,
    )
    pdt.assert_frame_equal(
        candidate_seed,
        pd.DataFrame([[0.5], [0.5]], index=level_index, columns=[0]),
    )


def test_evaluator_accepts_physical_level_iterate_when_optimizer_stalls(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Match legacy acceptance of a physical, non-closing root iterate."""
    evaluator = toy_evaluator
    continuum_rates = evaluator.calculate_continuum_coefficients(
        np.array([1.0e4])
    )[0][0]
    arguments = (
        0,
        1.0e9,
        1.0e4,
        np.array([0.5, 0.5]),
        continuum_rates,
        evaluator.population_geometries[0],
        evaluator.sobolev_inputs[0],
    )
    stalled_fractions = np.array([0.4, 0.6])
    expected_residual, expected_beta, expected_ion_ratio = (
        evaluator._calculate_level_state(
            0,
            1.0e9,
            1.0e4,
            stalled_fractions,
            continuum_rates,
            evaluator.population_geometries[0],
            evaluator.sobolev_inputs[0],
        )
    )
    assert np.max(np.abs(expected_residual)) > 1e-10

    def stalled_root(*args: object, **kwargs: object) -> SimpleNamespace:
        del args, kwargs
        return SimpleNamespace(success=False, x=stalled_fractions)

    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.root", stalled_root
    )
    actual = evaluator._calculate_level_solution(*arguments)

    npt.assert_allclose(actual[0], stalled_fractions)
    npt.assert_allclose(actual[1], expected_ion_ratio)
    npt.assert_allclose(actual[2], expected_beta)
    npt.assert_allclose(actual[3], expected_residual)


def test_evaluator_rebuilds_final_residual_and_is_deterministic(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Use final charge density for closure without mutating inputs."""
    evaluator = toy_evaluator
    level_index = evaluator.level_population_index
    estimators = evaluator.estimators_continuum

    class ChargeSolver:
        def __init__(self) -> None:
            self.calls = 0

        def solve(self, **kwargs: object) -> tuple[pd.DataFrame, pd.Series]:
            del kwargs
            self.calls += 1
            ion_index = pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
            )
            return (
                # Ion densities and charge-solved electron density (cm⁻³).
                pd.DataFrame([[8.0e9], [2.0e9]], index=ion_index, columns=[0]),
                pd.Series([2.0e9], index=[0]),
            )

    class ThermalSolver:
        def solve(
            self,
            thermal_electron_distribution: ThermalElectronEnergyDistribution,
            **kwargs: object,
        ) -> tuple[pd.Series, pd.Series]:
            del kwargs
            density = pd.Series(
                thermal_electron_distribution.number_density.to_value("cm^-3"),
                index=[0],
            )
            return density, density

    charge_solver = ChargeSolver()
    evaluator.ion_population_solver = charge_solver
    evaluator.thermal_balance_solver = ThermalSolver()
    level_seed = pd.DataFrame([[0.5], [0.5]], index=level_index, columns=[0])
    original_seed = level_seed.copy(deep=True)

    first_result = evaluator.evaluate(
        np.array([1.0e9]),  # Trial electron density (cm⁻³).
        np.array([1.0e4]),  # Shell electron temperature (K).
        level_seed,
    )
    second_result = evaluator.evaluate(
        np.array([1.0e9]),  # Trial electron density (cm⁻³).
        np.array([1.0e4]),  # Shell electron temperature (K).
        level_seed,
    )

    expected_final_residual = evaluator._calculate_level_state(
        0,
        2.0e9,  # Final electron density (cm⁻³).
        1.0e4,  # Shell electron temperature (K).
        first_result.normalized_population[0].to_numpy(),
        evaluator.calculate_continuum_coefficients(np.array([1.0e4]))[0][
            0
        ],
        evaluator.population_geometries[0],
        evaluator.sobolev_inputs[0],
        first_result.absolute_level_population.iloc[:, 0].to_numpy(),
    )[0]
    npt.assert_allclose(
        first_result.level_residual[0].to_numpy(), expected_final_residual
    )
    assert not np.allclose(
        first_result.trial_level_residual.to_numpy(),
        first_result.level_residual.to_numpy(),
    )
    assert charge_solver.calls == 2
    npt.assert_allclose(first_result.charge_residual.to_numpy(), [0.0])
    npt.assert_allclose(first_result.electron_residual.to_numpy(), [1.0])
    npt.assert_allclose(first_result.fractional_heating.to_numpy(), [2.0e9])
    pdt.assert_frame_equal(
        first_result.absolute_level_population,
        second_result.absolute_level_population,
    )
    pdt.assert_frame_equal(
        first_result.level_residual, second_result.level_residual
    )
    pdt.assert_frame_equal(level_seed, original_seed)


def test_evaluator_closes_synthetic_one_shell_thermal_root(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Close electron and thermal residuals at a known one-shell root."""
    target_electron_density = 2.0e9
    target_electron_temperature = 8.0e3
    maximum_electron_density = 1.0e10
    radiation_temperature = 1.0e4

    class ChargeSolver:
        def solve(self, **kwargs: object) -> tuple[pd.DataFrame, pd.Series]:
            del kwargs
            ion_index = pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
            )
            return (
                pd.DataFrame(
                    [[8.0e9], [target_electron_density]],
                    index=ion_index,
                    columns=[0],
                ),
                pd.Series([target_electron_density], index=[0]),
            )

    class ThermalSolver:
        def solve(
            self,
            thermal_electron_distribution: ThermalElectronEnergyDistribution,
            **kwargs: object,
        ) -> tuple[pd.Series, pd.Series]:
            del kwargs
            fractional_heating = pd.Series(
                thermal_electron_distribution.temperature.to_value(u.K)
                / target_electron_temperature
                - 1.0,
                index=[0],
            )
            return fractional_heating, fractional_heating

    toy_evaluator.ion_population_solver = ChargeSolver()
    toy_evaluator.thermal_balance_solver = ThermalSolver()
    level_seed = pd.DataFrame(
        [[0.5], [0.5]],
        index=toy_evaluator.level_population_index,
        columns=[0],
    )

    final_evaluation = toy_evaluator.evaluate(
        np.array([0.2 * maximum_electron_density]),
        np.array([0.8 * radiation_temperature]),
        level_seed,
    )

    npt.assert_allclose(final_evaluation.electron_residual, [0.0], atol=1e-10)
    npt.assert_allclose(final_evaluation.fractional_heating, [0.0], atol=1e-10)
    npt.assert_allclose(final_evaluation.charge_residual, [0.0], atol=1e-10)


def test_evaluator_aligns_collisional_temperature_scaling_by_transition(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Apply candidate/reference ratios by transition label, not row order."""
    transition_index = pd.MultiIndex.from_tuples(
        [
            (1, 0, 0, 0, 0, 1),
            (1, 0, 0, 0, 1, 2),
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

    class ElectronRateSolver:
        all_collisional_strengths_index = transition_index

        def solve(self, temperature: u.Quantity) -> pd.DataFrame:
            if temperature.to_value(u.K)[0] == 1.0e4:
                values = [1.0, 1.0, 1.0, 1.0]
            else:
                values = [2.0, 3.0, 4.0, 5.0]
            return pd.DataFrame(values, columns=[0])

    class ThermalSolver:
        def solve(
            self,
            collisional_excitation_rate_coefficient: pd.DataFrame,
            collisional_deexcitation_rate_coefficient: pd.DataFrame,
            **kwargs: object,
        ) -> tuple[pd.Series, pd.Series]:
            del kwargs
            return (
                collisional_excitation_rate_coefficient.iloc[:, 0],
                collisional_deexcitation_rate_coefficient.iloc[:, 0],
            )

    reverse_index = transition_index[::-1]
    toy_evaluator.rate_matrix_solver.electron_rate_solver = ElectronRateSolver()
    toy_evaluator.thermal_balance_solver = ThermalSolver()
    toy_evaluator.reference_electron_temperature = np.array([1.0e4]) * u.K
    toy_evaluator.thermal_balance_arguments = {
        "collisional_excitation_rate_coefficient": pd.DataFrame(
            [10.0, 20.0], index=reverse_index, columns=[0]
        ),
        "collisional_deexcitation_rate_coefficient": pd.DataFrame(
            [100.0, 200.0], index=reverse_index, columns=[0]
        ),
    }
    state = SimpleNamespace(
        collisional_ionization_rate_coefficient=pd.DataFrame(columns=[0]),
        level_to_continuum_saha_factor=pd.DataFrame(columns=[0]),
    )

    excitation, deexcitation = toy_evaluator._calculate_heating(
        ThermalElectronEnergyDistribution(
            0.0 * u.erg,
            np.array([2.0e4]) * u.K,
            np.array([1.0e9]) / u.cm**3,
        ),
        pd.DataFrame([[1.0]], columns=[0]),
        None,
        state,
    )

    pdt.assert_series_equal(
        excitation,
        pd.Series([30.0, 40.0], index=reverse_index, name=0),
    )
    pdt.assert_series_equal(
        deexcitation,
        pd.Series([500.0, 800.0], index=reverse_index, name=0),
    )
