from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import root as scipy_root

from tardis import constants as const
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluator,
    calculate_lte_populations,
)
from tardis.plasma.equilibrium.inputs import (
    DiluteRadiationInputs,
    FixedRateState,
    ShellNumberDensity,
    SobolevInputs,
)
from tardis.plasma.equilibrium.rate_matrix import (
    RateMatrix,
    assemble_ion_rate_matrices,
)
from tardis.plasma.equilibrium.rates import RadiativeRatesSolver
from tardis.transport.montecarlo.estimators import init_estimators_continuum


@pytest.mark.parametrize("temperature_ratio", [0.5, 1.0, 1.5])
def test_nebular_factor_has_ground_recombination_limit(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    temperature_ratio: float,
) -> None:
    """For zeta=delta=1, nebular Phi is W sqrt(Te/Trad) times Saha Phi.

    Regime: one two-level hydrogen atom, fixed Trad=10000 K and W=0.2.
    The analytic Saha equation and the nebular limiting expression supply an
    independent oracle; rtol=1e-13 permits only floating-point arithmetic.
    Varying Te tests candidate dependence rather than a frozen LTE fallback.
    """
    temperature = 1e4
    dilution = 0.2
    partition = toy_evaluator.calculate_continuum_coefficients(
        np.array([temperature])
    )[5]
    toy_evaluator.dilute_radiation = DiluteRadiationInputs(
        np.array([temperature]),
        np.array([dilution]),
        pd.DataFrame(
            [[1.0, 1.0]],
            index=toy_evaluator.ionization_data.index,
            columns=[5000.0, 20000.0],
        ),
        pd.DataFrame(),
        partition,
        1.0,
    )
    actual = toy_evaluator.calculate_nebular_ionization_factor(
        np.array([temperature * temperature_ratio])
    )
    neutral_partition = 2.0 + 4.0 * np.exp(
        -1e-12 / (const.k_B.cgs.value * temperature)
    )
    saha = (
        2.0
        / neutral_partition
        * (
            2.0
            * np.pi
            * const.m_e.cgs.value
            * const.k_B.cgs.value
            * temperature
            / const.h.cgs.value**2
        )
        ** 1.5
        * np.exp(-2e-11 / (const.k_B.cgs.value * temperature))
    )
    npt.assert_allclose(
        actual.to_numpy(),
        dilution * np.sqrt(temperature_ratio) * saha,
        rtol=1e-13,
        atol=0.0,
    )


def test_background_excitation_uses_dilute_weights_not_input_fractions(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Rebuild background ions with prescribed dilute-LTE weights and conserve N.

    Regime: Li I has seed fractions [0.9,0.1], dilute weights [1,4], and new
    ion density 8 cm^-3. The exact normalized distribution is [1.6,6.4], not
    a rescaling of the seed. rtol=1e-14 covers only normalization roundoff.
    This verifies the excitation prescription wiring, not the atomic weights.
    """
    levels = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1), (3, 0, 0), (3, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    toy_evaluator.level_population_index = levels
    toy_evaluator.dilute_radiation = SimpleNamespace(
        level_boltzmann_factor=pd.DataFrame(
            [[1.0], [1.0], [1.0], [4.0]], index=levels, columns=[0]
        ),
    )
    state = SimpleNamespace(
        full_population=pd.DataFrame(
            [[1.0], [1.0], [9.0], [1.0]], index=levels, columns=[0]
        ),
        species_index=levels[:2],
        normalized_population=pd.DataFrame(
            [[0.25], [0.75]], index=levels[:2], columns=[0]
        ),
    )
    ions = pd.DataFrame(
        [[2.0], [8.0]],
        index=pd.MultiIndex.from_tuples([(1, 0), (3, 0)]),
        columns=[0],
    )
    actual = toy_evaluator._build_absolute_levels(state, ions)
    npt.assert_allclose(actual.loc[(3, 0), 0], [1.6, 6.4], rtol=1e-14, atol=0.0)
    npt.assert_allclose(actual.loc[(3, 0), 0].sum(), 8.0, rtol=1e-14, atol=0.0)


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
    population_geometry = ShellNumberDensity(
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
    (rates,), *_ = evaluator.calculate_continuum_coefficients(np.array([1.0e4]))
    (hot_rates,), *_ = evaluator.calculate_continuum_coefficients(
        np.array([2.0e4])
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


def test_evaluator_finds_same_unique_root_from_distinct_initial_guesses(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Find the same two-level equilibrium from distinct initial populations.

    Claim: A unique level equilibrium does not depend on its initial population.
    Regime: One shell, two neutral-hydrogen levels, fixed electron density and
    temperature, and two normalized initial guesses on opposite sides of the
    root.
    Verification: Initial-guess independence is a metamorphic relation; both
    solutions must also conserve population and close the level equations.
    """
    level_index = toy_evaluator.level_population_index
    low_excitation_initial_guess = pd.DataFrame(
        [[0.9], [0.1]], index=level_index, columns=[0]
    )
    high_excitation_initial_guess = pd.DataFrame(
        [[0.1], [0.9]], index=level_index, columns=[0]
    )

    low_excitation_result = toy_evaluator.evaluate(
        np.array([1.0e9]), np.array([1.0e4]), low_excitation_initial_guess
    )
    high_excitation_result = toy_evaluator.evaluate(
        np.array([1.0e9]), np.array([1.0e4]), high_excitation_initial_guess
    )

    npt.assert_allclose(
        low_excitation_result.normalized_population,
        high_excitation_result.normalized_population,
        rtol=1e-12,
        atol=1e-12,
        err_msg="the equilibrium depends on the initial level population",
    )
    npt.assert_allclose(
        low_excitation_result.normalized_population.sum(axis=0),
        1.0,
        rtol=0.0,
        atol=1e-15,
        err_msg="the equilibrium populations do not sum to one",
    )
    npt.assert_allclose(
        low_excitation_result.level_residual,
        0.0,
        rtol=0.0,
        atol=1e-10,
        err_msg="the equilibrium does not close the level equations",
    )


def test_evaluator_refines_physical_level_iterate_when_optimizer_stalls(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Refine a physical population when the level solver stops early.

    Claim: This stalled physical iterate can be refined to level equilibrium.
    Regime: One shell and two neutral-hydrogen levels at fixed density and
    temperature.
    Verification: Replace HYBR's result with a physical but unsolved iterate;
    the returned population must differ and close the level equations to the
    workflow's 1e-10 tolerance.
    """
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
        evaluator.shell_number_densities[0],
        evaluator.sobolev_inputs[0],
    )
    stalled_fractions = np.array([0.4, 0.6])

    # Replace SciPy's root finder with a controlled stopped result. The
    # SimpleNamespace supplies only the status and population returned by
    # SciPy, keeping this check focused on the evaluator's acceptance rule.
    def stalled_root(*args: object, **kwargs: object) -> SimpleNamespace:
        del args, kwargs
        return SimpleNamespace(success=False, x=stalled_fractions)

    # monkeypatch limits the replacement to this test and restores the real
    # root finder afterward.
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.root", stalled_root
    )
    actual = evaluator._calculate_level_solution(*arguments)

    assert actual[0].sum() == 1.0
    assert np.all(actual[0] >= 0.0)
    npt.assert_allclose(actual[3], 0.0, rtol=0.0, atol=1e-10)


def test_evaluator_retries_open_level_equations_despite_optimizer_success(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Select a retry by its physical residual, ignoring optimizer status.

    Claim: A near-ground seed that gives an open level solution must be retried
    from a less concentrated seed when another physical solution closes.
    Regime: One shell, two neutral-hydrogen levels, fixed 1e9 cm^-3 electron
    density and 10000 K temperature; the first two optimizers falsely report
    success and a later retry reports failure despite closing the equations.
    Verification: This is a control-flow check using the evaluator's physical
    equations. A separate SciPy root call establishes a closed candidate; the
    returned population is checked against those equations. The 1e-10
    tolerance is the workflow's level-closure limit.
    """
    evaluator = toy_evaluator
    electron_density = 1.0e9  # cm^-3
    electron_temperature = 1.0e4  # K
    continuum_rates = evaluator.calculate_continuum_coefficients(
        np.array([electron_temperature])
    )[0][0]
    population_geometry = evaluator.shell_number_densities[0]
    sobolev_inputs = evaluator.sobolev_inputs[0]

    def physical_residual(
        level_fractions: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        return evaluator._calculate_level_state(
            0,
            electron_density,
            electron_temperature,
            level_fractions,
            continuum_rates,
            population_geometry,
            sobolev_inputs,
        )[0]

    closed_fractions = scipy_root(
        physical_residual, np.array([0.5, 0.5]), options={"xtol": 1e-12}
    ).x
    stalled_fractions = np.array([0.4, 0.6])
    assert np.all(closed_fractions >= 0.0)
    npt.assert_allclose(closed_fractions.sum(), 1.0, rtol=0.0, atol=1e-10)
    npt.assert_allclose(
        physical_residual(closed_fractions), 0.0, rtol=0.0, atol=1e-10
    )
    assert np.max(np.abs(physical_residual(stalled_fractions))) > 1e-10

    # The synthetic 1e-10 split distinguishes the original 1e-12 excited
    # fraction from a materially less concentrated retry seed.
    def controlled_root(
        _residual: object, seed: ArrayLike, **kwargs: object
    ) -> SimpleNamespace:
        del kwargs
        if np.asarray(seed)[1] > 1e-10:
            return SimpleNamespace(success=False, x=closed_fractions)
        return SimpleNamespace(success=True, x=stalled_fractions)

    def controlled_log_solver(
        _residual: object, log_ratios: ArrayLike, **kwargs: object
    ) -> SimpleNamespace:
        del kwargs
        if np.exp(np.asarray(log_ratios)[0]) > 1e-10:
            fractions = closed_fractions
            success = False
        else:
            fractions = stalled_fractions
            success = True
        return SimpleNamespace(
            success=success,
            x=np.log(fractions[1:] / fractions[0]),
        )

    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.root", controlled_root
    )
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.least_squares",
        controlled_log_solver,
    )
    actual = evaluator._calculate_level_solution(
        0,
        electron_density,
        electron_temperature,
        np.array([1.0 - 1e-12, 1e-12]),
        continuum_rates,
        population_geometry,
        sobolev_inputs,
    )

    assert np.all(actual[0] >= 0.0)
    npt.assert_allclose(actual[0].sum(), 1.0, rtol=0.0, atol=1e-15)
    npt.assert_allclose(
        physical_residual(actual[0]),
        0.0,
        rtol=0.0,
        atol=1e-10,
        err_msg="optimizer success hid an open physical level solution",
    )
    npt.assert_allclose(actual[3], 0.0, rtol=0.0, atol=1e-10)


def test_evaluator_reports_open_level_solution_after_all_retries(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Expose the residual when no retry closes a physical trial iterate.

    Claim: An outer thermal trial can use a physical best-effort population,
    but its unsolved level balance must remain visible to final validation.
    Regime: One shell with two neutral-hydrogen levels at fixed density and
    temperature.
    Verification: Both optimizers return the same unsolved physical iterate;
    the evaluator returns its nonzero equation residual.
    """
    evaluator = toy_evaluator
    continuum_rates = evaluator.calculate_continuum_coefficients(
        np.array([1.0e4])
    )[0][0]
    stalled_fractions = np.array([0.4, 0.6])

    def stalled_root(*args: object, **kwargs: object) -> SimpleNamespace:
        del args, kwargs
        return SimpleNamespace(success=True, x=stalled_fractions)

    def stalled_log_solver(*args: object, **kwargs: object) -> SimpleNamespace:
        del args, kwargs
        return SimpleNamespace(success=True, x=np.array([np.log(1.5)]))

    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.root", stalled_root
    )
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.least_squares",
        stalled_log_solver,
    )
    actual = evaluator._calculate_level_solution(
        0,
        1.0e9,
        1.0e4,
        stalled_fractions.copy(),
        continuum_rates,
        evaluator.shell_number_densities[0],
        evaluator.sobolev_inputs[0],
    )
    npt.assert_array_equal(actual[0], stalled_fractions)
    assert np.max(np.abs(actual[3])) > 1e-10


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
    level_initial_guess = pd.DataFrame(
        [[0.5], [0.5]], index=level_index, columns=[0]
    )
    original_initial_guess = level_initial_guess.copy(deep=True)

    first_result = evaluator.evaluate(
        np.array([1.0e9]),  # Trial electron density (cm⁻³).
        np.array([1.0e4]),  # Shell electron temperature (K).
        level_initial_guess,
    )
    second_result = evaluator.evaluate(
        np.array([1.0e9]),  # Trial electron density (cm⁻³).
        np.array([1.0e4]),  # Shell electron temperature (K).
        level_initial_guess,
    )

    continuum_rate_coeff = evaluator.calculate_continuum_coefficients(
        np.array([1.0e4])
    )[0][0]

    expected_final_residual = evaluator._calculate_level_state(
        0,
        2.0e9,  # Final electron density (cm⁻³).
        1.0e4,  # Shell electron temperature (K).
        first_result.normalized_population[0].to_numpy(),
        continuum_rate_coeff,
        evaluator.shell_number_densities[0],
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
    pdt.assert_frame_equal(level_initial_guess, original_initial_guess)


def test_evaluator_closes_at_known_one_shell_thermal_root(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Close electron and thermal balance at an exact one-shell root.

    Claim: The electron and heating balances vanish at the prescribed density
    and temperature.
    Regime: One shell with a charge solution of 2e9 cm^-3 and a thermal root
    at 8000 K.
    Verification: Both values are fixed analytically by independent linear
    balance relations supplied to the evaluator.
    """
    target_electron_density = np.array([2.0e9])
    target_electron_temperature = np.array([8.0e3])
    maximum_electron_density = 1.0e10
    radiation_temperature = 1.0e4

    # Prescribe the ion and electron densities at the analytic charge root,
    # independent of the full ion-population calculation.
    class ChargeSolver:
        def solve(self, **kwargs: object) -> tuple[pd.DataFrame, pd.Series]:
            del kwargs
            ion_index = pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
            )
            return (
                pd.DataFrame(
                    [[8.0e9], target_electron_density],
                    index=ion_index,
                    columns=[0],
                ),
                pd.Series(target_electron_density, index=[0]),
            )

    # Use a linear heating relation with an exact zero at 8000 K so the
    # expected thermal root is known before the evaluator is called.
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
    level_initial_guess = pd.DataFrame(
        [[0.5], [0.5]],
        index=toy_evaluator.level_population_index,
        columns=[0],
    )

    final_evaluation = toy_evaluator.evaluate(
        target_electron_density,
        target_electron_temperature,
        level_initial_guess,
    )

    assert target_electron_density / maximum_electron_density == 0.2
    assert target_electron_temperature / radiation_temperature == 0.8
    npt.assert_allclose(
        final_evaluation.electron_residual, [0.0], rtol=0.0, atol=1e-15
    )
    npt.assert_allclose(
        final_evaluation.fractional_heating, [0.0], rtol=0.0, atol=1e-15
    )
    npt.assert_allclose(
        final_evaluation.charge_residual, [0.0], rtol=0.0, atol=1e-15
    )


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

    excitation, deexcitation, _ = toy_evaluator._calculate_heating(
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


def test_prepare_species_rates_uses_solved_line_amid_duplicate_lines() -> None:
    """Pair the solved species' line with its own J and beta row.

    Claim: duplicate level pairs of an unsolved species do not change which
    line-list row supplies the mean intensity and escape probability.
    Regime: one H I line between two helium lines that share a level pair.
    Verification: direct Einstein relations with distinct helium intensities,
    so reading a wrong row changes the rates; exact arithmetic is expected.
    """
    line_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0, 1), (1, 0, 0, 1), (2, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )
    levels = pd.DataFrame(
        {"energy": [0.0, 1.0e-12, 0.0], "g": [2.0, 4.0, 1.0]},
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )
    hydrogen_line = pd.DataFrame(
        {"A_ul": [2.0], "B_ul": [0.5], "B_lu": [1.0], "nu": [1.0]},
        index=line_index[[1]],
    )
    j_blues = pd.DataFrame([[7.0], [3.0], [11.0]], index=line_index)

    (rates,) = RateMatrix(
        RadiativeRatesSolver(hydrogen_line), ZeroElectronRateSolver(), levels
    ).prepare_species_rates(
        j_blues,
        ThermalElectronEnergyDistribution(
            0.0 * u.erg, np.array([1.0e4]) * u.K, np.array([1.0e9]) / u.cm**3
        ),
        (1, 0),
        line_index,
    )

    upward = (rates.source_level_idx == 0) & (rates.destination_level_idx == 1)
    downward = (rates.source_level_idx == 1) & (
        rates.destination_level_idx == 0
    )
    # B_lu * J and A_ul + B_ul * J with the H I row's J = 3.
    npt.assert_array_equal(rates.radiative_rate_coefficient[upward], [3.0])
    npt.assert_array_equal(rates.radiative_rate_coefficient[downward], [3.5])
    # Both directions read beta from the H I line's row in the full list.
    npt.assert_array_equal(rates.beta_line_idx, [1, 1])


@pytest.mark.parametrize("supply_background", [True, False])
def test_evaluator_repopulates_zeroed_background_ion(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    supply_background: bool,
) -> None:
    """Repopulate a zeroed trace ion from explicit background fractions.

    Regime: Li I has zero levels in shell 25, a populated [3, 1] distribution
    in shell 53, and zero ion density in shell 7. New ion totals are 2e-23,
    8, and 0 cm^-3. Given fractions [0.8, 0.2], the revived levels must be
    [1.6e-23, 4e-24]; the existing distribution must instead rescale to [6, 2].
    Verification uses these analytic values and per-ion conservation with
    rtol=1e-14 (roundoff only), no absolute tolerance that could hide the trace
    population. Without explicit fractions, legacy zero-background behavior
    is preserved. This verifies rebuilding, not the excitation prescription.
    """
    evaluator = toy_evaluator
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1), (3, 0, 0), (3, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    columns = pd.Index([25, 53, 7])
    evaluator.level_population_index = level_index
    evaluator.background_level_population_fraction = (
        pd.DataFrame(
            [[0.5] * 3, [0.5] * 3, [0.8] * 3, [0.2] * 3],
            index=level_index,
            columns=columns,
        )
        if supply_background
        else None
    )
    state = SimpleNamespace(
        species_index=level_index[:2],
        full_population=pd.DataFrame(
            [[1.0] * 3, [1.0] * 3, [0.0, 3.0, 0.0], [0.0, 1.0, 0.0]],
            index=level_index,
            columns=columns,
        ),
        normalized_population=pd.DataFrame(
            [[0.25] * 3, [0.75] * 3],
            index=level_index[:2],
            columns=columns,
        ),
    )
    ions = pd.DataFrame(
        [[2.0] * 3, [2.0e-23, 8.0, 0.0]],
        index=level_index.droplevel("level_number").unique(),
        columns=columns,
    )
    levels = evaluator._build_absolute_levels(state, ions)
    npt.assert_allclose(
        levels.loc[(3, 0), 25],
        [1.6e-23, 4e-24] if supply_background else [0.0, 0.0],
        rtol=1e-14,
        atol=0.0,
    )
    npt.assert_allclose(levels.loc[(3, 0), 53], [6.0, 2.0], rtol=1e-14)
    npt.assert_array_equal(levels.loc[(3, 0), 7], 0.0)
    npt.assert_allclose(levels.loc[(1, 0), 25], [0.5, 1.5], rtol=1e-14)
    if supply_background:
        pdt.assert_frame_equal(
            levels.groupby(level=["atomic_number", "ion_number"]).sum(),
            ions,
            rtol=1e-14,
            atol=0.0,
        )


def test_evaluator_reports_stage_timing(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    caplog: pytest.LogCaptureFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Expose progress through one complete equilibrium evaluation.

    Claim
    -----
    Each evaluator stage remains visible at DEBUG without flooding INFO.

    Regime
    ------
    The minimal one-shell evaluator exercises continuum, level, population,
    opacity, thermal, and final-residual orchestration.

    Verification
    ------------
    This is an orchestration assertion over stage labels only; it deliberately
    makes no machine-dependent assertion about elapsed time.
    """
    caplog.set_level("DEBUG", logger="tardis.plasma.equilibrium.evaluator")
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.perf_counter", lambda: 0.0
    )

    toy_evaluator.evaluate(
        np.array([1.0e9]),
        np.array([1.0e4]),
        pd.DataFrame([[0.8], [0.2]], columns=[0]),
    )

    for stage in (
        "continuum coefficients",
        "level populations",
        "ion populations",
        "absolute and LTE populations",
        "Sobolev opacities",
        "thermal balance",
        "final level residual",
    ):
        assert f"finished {stage}" in caplog.text
    assert "Finished plasma equilibrium evaluator call 1" in caplog.text
    assert not any(
        record.levelname == "INFO"
        and "plasma equilibrium evaluator call" in record.message.lower()
        for record in caplog.records
    )
    assert not any(
        record.levelname == "INFO"
        and "equilibrium level solve" in record.message
        for record in caplog.records
    )


def test_continuum_coefficients_preserve_nonconsecutive_shell_labels(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Keep physical shell labels when evaluating a restricted shell set."""
    shell_columns = pd.Index([53])
    toy_evaluator.elemental_number_density.columns = shell_columns

    (
        _,
        level_to_continuum_saha_factor,
        collisional_ionization_rate_coefficient,
        lte_ionization_factor,
        continuum_coefficients,
        partition_function,
        level_boltzmann_factor,
    ) = toy_evaluator.calculate_continuum_coefficients(np.array([1.0e4]))

    frames = (
        level_to_continuum_saha_factor,
        collisional_ionization_rate_coefficient,
        lte_ionization_factor,
        continuum_coefficients.photoionization,
        continuum_coefficients.stimulated_recombination,
        continuum_coefficients.spontaneous_recombination,
        continuum_coefficients.collisional_ionization,
        partition_function,
        level_boltzmann_factor,
    )
    for frame in frames:
        assert frame.columns.equals(shell_columns)

    lte_ion_population, lte_level_population = calculate_lte_populations(
        lte_ionization_factor,
        partition_function,
        toy_evaluator.elemental_number_density,
        pd.Series([1.0e8], index=shell_columns),
        level_boltzmann_factor,
        toy_evaluator.levels,
    )
    assert lte_ion_population.columns.equals(shell_columns)
    assert lte_level_population.columns.equals(shell_columns)


def test_fixed_ionization_enters_level_and_ion_stage_balances(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Use one fixed ionization coefficient in both coupled balances."""
    transition_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
        ],
    )
    fixed_ionization = pd.DataFrame(
        [[4.0]], index=transition_index, columns=[0]
    )
    toy_evaluator.fixed_rates = FixedRateState(
        ionization_rates=fixed_ionization
    )
    continuum_rates = toy_evaluator.calculate_continuum_coefficients(
        np.array([1.0e4])
    )[0][0]
    level_fractions = np.array([0.5, 0.5])

    _, _, fixed_ratio = toy_evaluator._calculate_level_state(
        0,
        1.0e9,
        1.0e4,
        level_fractions,
        continuum_rates,
        toy_evaluator.shell_number_densities[0],
        toy_evaluator.sobolev_inputs[0],
    )
    toy_evaluator.fixed_rates = None
    _, _, baseline_ratio = toy_evaluator._calculate_level_state(
        0,
        1.0e9,
        1.0e4,
        level_fractions,
        continuum_rates,
        toy_evaluator.shell_number_densities[0],
        toy_evaluator.sobolev_inputs[0],
    )

    recombination_rate = (
        continuum_rates.spontaneous_recombination[1]
        + continuum_rates.stimulated_recombination[1]
        + continuum_rates.collisional_recombination[1] * 1.0e9
    ) * 1.0e9
    npt.assert_allclose(
        fixed_ratio - baseline_ratio,
        4.0 / recombination_rate,
        rtol=1e-14,
    )

    zero_rates = pd.DataFrame(0.0, index=transition_index, columns=[0])
    ion_matrices, _ = assemble_ion_rate_matrices(
        zero_rates,
        zero_rates,
        zero_rates,
        zero_rates,
        fixed_ionization_rates_df=fixed_ionization,
    )
    assert ion_matrices.at[1, 0][0, 0] == -4.0


def test_case_b_suppresses_only_ground_state_recombination(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Switch between Case A and Case B without changing other rates."""
    continuum_rates = toy_evaluator.calculate_continuum_coefficients(
        np.array([1.0e4])
    )[0][0]
    arguments = (
        0,
        1.0e9,
        1.0e4,
        np.array([0.5, 0.5]),
        continuum_rates,
        toy_evaluator.shell_number_densities[0],
        toy_evaluator.sobolev_inputs[0],
    )

    toy_evaluator.case_b_species = frozenset()
    case_a_ratio = toy_evaluator._calculate_level_state(*arguments)[2]
    toy_evaluator.case_b_species = frozenset({(1, 0)})
    case_b_ratio = toy_evaluator._calculate_level_state(*arguments)[2]

    assert case_b_ratio > case_a_ratio


def test_evaluator_reports_adiabatic_estimate_without_changing_balance(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Compare electron expansion work with the unchanged heating deficit.

    Claim: diagnostic cooling equals electron pressure times 3/t and does
    not enter the thermal residual. Regime: one homologously expanding cell,
    with an age in days distinct from the Monte Carlo normalization time.
    Verification: the ideal-electron-gas energy density loses 2u/t; synthetic
    heating and cooling retain their prescribed net and fractional rates.
    """

    class ThermalSolver:
        def solve(self, **kwargs: object) -> tuple[pd.Series, pd.Series]:
            self.last_budget = pd.DataFrame(
                {"total_heating": [3.0e-11], "total_cooling": [1.0e-11]},
                index=kwargs["level_population"].columns,
            )
            return pd.Series([2.0e-11]), pd.Series([2.0])

    toy_evaluator.thermal_balance_solver = ThermalSolver()
    toy_evaluator.time_explosion = 2.0 * u.day
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        np.array([1.0e4]) * u.K,
        np.array([2.0e9]) / u.cm**3,
    )
    net_heating, fractional_heating, budget = toy_evaluator._calculate_heating(
        electron_distribution,
        pd.DataFrame([[1.0]], columns=[17]),
        None,
        SimpleNamespace(
            collisional_ionization_rate_coefficient=pd.DataFrame(columns=[17]),
            level_to_continuum_saha_factor=pd.DataFrame(columns=[17]),
        ),
    )
    electron_energy_density = (
        1.5 * const.k_B.cgs * 1.0e4 * u.K * 2.0e9 / u.cm**3
    )
    expected_cooling = (2 * electron_energy_density / (2 * u.day)).to_value(
        u.erg / (u.cm**3 * u.s)
    )
    npt.assert_allclose(budget["adiabatic_cooling_estimate"], expected_cooling)
    npt.assert_allclose(budget["electron_density"], [2.0e9])
    npt.assert_allclose(budget["electron_temperature"], [1.0e4])
    npt.assert_allclose(budget["time_explosion"], [172800.0])
    npt.assert_allclose(net_heating, [2.0e-11], rtol=1e-14)
    npt.assert_allclose(fractional_heating, [2.0], rtol=1e-14)
