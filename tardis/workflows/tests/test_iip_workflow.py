from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import NamedTuple

import astropy.units as u
import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis import constants as const
from tardis.conftest import assert_regression_dataframe
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.equilibrium.evaluator import calculate_lte_populations
from tardis.plasma.equilibrium.thermal_balance import (
    ElectronThermalBalanceSolver,
)
from tardis.plasma.exceptions import PlasmaIonizationError
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow

PLASMA_SOLVER_REGRESSION_OUTPUTS = (
    "electron_densities",
    "t_electrons",
    "link_t_rad_t_electron",
    "p_fb_deactivation",
    "chi_bf",
    "stimulated_emission_factor",
    "j_blues",
)

INITIAL_PLASMA_SOLVER_REGRESSION_OUTPUTS = (
    "ion_number_density",
    "tau_sobolevs",
    "beta_sobolev",
    "level_number_density",
    *PLASMA_SOLVER_REGRESSION_OUTPUTS,
)


@pytest.fixture
def ctardis_compare_config(
    tardis_regression_path: Path,
) -> Configuration:
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    config.plasma.nlte.species = [(1, 0)]
    return config


@pytest.fixture
def type_iip_workflow(
    ctardis_compare_config: Configuration,
) -> TypeIIPWorkflow:
    return TypeIIPWorkflow(ctardis_compare_config)


class FirstTransportThermalProblem(NamedTuple):
    """Comparison workflow ready for its first thermal-balance solve."""

    workflow: TypeIIPWorkflow
    initial_guess: npt.NDArray[np.float64]
    lower_bound: npt.NDArray[np.float64]
    upper_bound: npt.NDArray[np.float64]
    maximum_electron_density: npt.NDArray[np.float64]


@pytest.fixture
def first_transport_thermal_problem(
    type_iip_workflow: TypeIIPWorkflow,
) -> FirstTransportThermalProblem:
    """Return the comparison workflow after its first transport.

    The plasma holds the first Monte Carlo estimators, the thermal-balance
    evaluator is frozen, and no thermal solve has run.
    """
    workflow = type_iip_workflow
    workflow.show_progress_bars = False
    workflow.solve_montecarlo(
        workflow.solve_opacity(), workflow.real_packet_count
    )
    estimated_values, _ = workflow.get_convergence_estimates()
    workflow.solve_simulation_state(estimated_values)
    continuum_estimators, j_blues = workflow.update_estimators()
    workflow.solve_plasma(continuum_estimators, j_blues)
    return FirstTransportThermalProblem(
        workflow, *workflow.prepare_thermal_balance_problem()
    )


@pytest.fixture
def ctardis_reference_path(tardis_regression_path: Path) -> Path:
    """Return the directory containing C-TARDIS workflow references."""
    return tardis_regression_path / "tardis" / "workflows" / "tests"


def test_workflow_initial_populations_match_ctardis(
    type_iip_workflow: TypeIIPWorkflow,
    ctardis_reference_path: Path,
) -> None:
    """Compare standard initial plasma populations with C-TARDIS data."""
    # The implementations are independent; tolerances capture initialization
    # parity rather than require identical nonlinear iterates.
    plasma = type_iip_workflow.plasma_solver
    ion_number_density = pd.read_hdf(
        ctardis_reference_path / "ctardis_ion_density_init_nlte.h5",
        key="data",
    )
    level_number_density = pd.read_hdf(
        ctardis_reference_path / "ctardis_level_number_density_init_nlte.h5",
        key="data",
    )
    electron_densities = pd.read_hdf(
        ctardis_reference_path / "ctardis_electron_densities_init_nlte.h5",
        key="data",
    )
    electron_temperatures = pd.read_hdf(
        ctardis_reference_path / "ctardis_t_electrons_init_nlte.h5",
        key="data",
    )

    pd.testing.assert_frame_equal(
        plasma.ion_number_density,
        ion_number_density,
        rtol=1.5e-2,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        plasma.level_number_density,
        level_number_density,
        rtol=1.5e-2,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_series_equal(
        plasma.electron_densities,
        electron_densities,
        rtol=1e-4,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    np.testing.assert_allclose(
        plasma.t_electrons,
        electron_temperatures.to_numpy().ravel(),
        rtol=1e-9,
        atol=0.0,
    )


def test_workflow_initial_opacity_matches_ctardis(
    type_iip_workflow: TypeIIPWorkflow,
    ctardis_reference_path: Path,
) -> None:
    """Compare standard initial continuum and Sobolev opacity with C-TARDIS."""
    # The implementations are independent; tolerances capture initialization
    # parity rather than require identical nonlinear iterates.
    tau_sobolev = pd.read_hdf(
        ctardis_reference_path / "ctardis_tau_sobolevs_init_nlte.h5",
        key="data",
    )
    beta_sobolev = pd.read_hdf(
        ctardis_reference_path / "ctardis_beta_sobolevs_init_nlte.h5",
        key="data",
    )
    p_fb_deactivation = pd.read_hdf(
        ctardis_reference_path / "ctardis_p_fb_deactivation_init_nlte.h5",
        key="data",
    )
    chi_bf = pd.read_hdf(
        ctardis_reference_path / "ctardis_chi_bf_init_nlte.h5",
        key="data",
    )

    # Sobolev values are stored differently between codes, so comparing raw data instead
    np.testing.assert_allclose(
        type_iip_workflow._tau_sobolev.to_numpy(),
        tau_sobolev.to_numpy(),
        rtol=1.5e-2,
        atol=0.0,
    )
    np.testing.assert_allclose(
        type_iip_workflow._beta_sobolev.to_numpy(),
        beta_sobolev.to_numpy(),
        rtol=1e-2,
        atol=0.0,
    )
    np.testing.assert_allclose(
        type_iip_workflow.continuum_opacity_state.p_fb_deactivation.to_numpy(),
        p_fb_deactivation.to_numpy(),
        rtol=1e-5,
        atol=0.0,
    )
    np.testing.assert_allclose(
        type_iip_workflow.continuum_opacity_state.chi_bf.to_numpy(),
        chi_bf.to_numpy(),
        rtol=5e-3,
        atol=0.0,
    )


def test_type_iip_workflow_initial_plasma_regression(
    type_iip_workflow: TypeIIPWorkflow,
    regression_data: RegressionData,
) -> None:
    """Compare the standard dilute-LTE bootstrap with legacy IIP outputs.

    Claim: Initial populations, continuum opacity, and Sobolev quantities
    retain observable legacy parity before the first Monte Carlo estimators.
    Regime: The five-shell Type IIP comparison configuration.
    Verification: Stored IIP outputs are independent of the standard plasma
    graph; ``1e-4`` permits their distinct nonlinear initialization paths.
    """
    plasma = type_iip_workflow.plasma_solver
    outputs = {
        "ion_number_density": plasma.ion_number_density,
        "tau_sobolevs": type_iip_workflow._tau_sobolev,
        "beta_sobolev": type_iip_workflow._beta_sobolev,
        "level_number_density": plasma.level_number_density,
        "electron_densities": plasma.electron_densities,
        "t_electrons": plasma.t_electrons,
        "link_t_rad_t_electron": plasma.link_t_rad_t_electron,
        "p_fb_deactivation": (
            type_iip_workflow.continuum_opacity_state.p_fb_deactivation
        ),
        "chi_bf": type_iip_workflow.continuum_opacity_state.chi_bf,
        "stimulated_emission_factor": plasma.stimulated_emission_factor,
        # The standard graph labels lines by their physical transition index;
        # the legacy regression stored the same values positionally.
        "j_blues": plasma.j_blues.reset_index(drop=True),
    }
    for output_name in INITIAL_PLASMA_SOLVER_REGRESSION_OUTPUTS:
        assert_regression_dataframe(
            regression_data,
            f"workflow_init_{output_name}",
            outputs[output_name],
            # The standard dilute-LTE bootstrap and legacy IIP initialization
            # use different nonlinear owners before the first MC estimator
            # snapshot. Preserve observable parity without requiring identical
            # intermediate iterates.
            rtol=1e-4,
        )


def test_thermal_balance_residual_delegates_to_evaluator() -> None:
    """Map a two-shell thermal-balance candidate to evaluator quantities.

    Claim: Each shell's density fraction and temperature ratio are converted
    to physical values, and the two balance errors retain shell order.
    Regime: Two shells with unequal density limits and radiation temperatures.
    Verification: The expected values follow by direct multiplication and
    inserting the recorded evaluator result.
    """
    level_initial_guess = pd.DataFrame([[0.6, 0.7], [0.4, 0.3]])

    # Record the physical values received and return fixed balance errors so
    # this test checks scaling and shell order without running a plasma solve.
    class RecordingEvaluator:
        def __init__(self) -> None:
            self.call_count = 0
            self.arguments: tuple[object, ...] | None = None
            # SimpleNamespace contains only the two evaluator results consumed
            # by calculate_residual.
            self.result = SimpleNamespace(
                electron_residual=pd.Series([0.1, -0.2]),
                fractional_heating=pd.Series([0.3, -0.4]),
            )

        def evaluate(
            self,
            trial_electron_density: npt.ArrayLike,
            electron_temperature: npt.ArrayLike,
            candidate_level_initial_guess: pd.DataFrame,
        ) -> SimpleNamespace:
            self.call_count += 1
            self.arguments = (
                np.asarray(trial_electron_density),
                np.asarray(electron_temperature),
                candidate_level_initial_guess,
            )
            return self.result

    evaluator = RecordingEvaluator()
    candidate = np.array([0.25, 0.8, 0.5, 1.1])
    maximum_electron_density = np.array([4.0e9, 6.0e9])
    # The noise measurement is not exercised, so the batches are unused.
    solver = ElectronThermalBalanceSolver(
        evaluator,
        np.array([1.0e4, 2.0e4]),
        level_initial_guess,
        maximum_electron_density,
        batch_evaluators=lambda: (),
        batch_packet_counts=np.ones(2, dtype=np.int64),
    )

    residual = solver.calculate_residual(candidate)

    assert evaluator.call_count == 1
    assert evaluator.arguments is not None
    trial_density, temperature, actual_initial_guess = evaluator.arguments
    np.testing.assert_array_equal(trial_density, [1.0e9, 3.0e9])
    np.testing.assert_allclose(
        temperature,
        [8.0e3, 2.2e4],
        rtol=1e-15,  # Floating-point multiplication only.
    )
    pd.testing.assert_frame_equal(actual_initial_guess, level_initial_guess)
    np.testing.assert_array_equal(residual, [0.1, 0.3, -0.2, -0.4])


def test_packet_batches_sum_to_single_transport(
    type_iip_workflow: TypeIIPWorkflow,
) -> None:
    """Transport packets in consecutive batches without changing the result.

    Claim: The estimators and escaping packets of four packet batches equal
    those of one transport of the same packets.
    Regime: Five-shell pure-hydrogen comparison model, first transport.
    Verification: Every packet has its own random-number seed, so batching
    only changes the order in which packet contributions are summed. The
    1e-12 tolerance covers that floating-point reordering.
    """
    workflow = type_iip_workflow
    workflow.show_progress_bars = False
    opacity_states = workflow.solve_opacity()
    results = {}
    for estimator_batch_count in (1, 4):
        workflow.solve_montecarlo(
            opacity_states,
            workflow.real_packet_count,
            estimator_batch_count=estimator_batch_count,
        )
        state = workflow.transport_state
        results[estimator_batch_count] = {
            "output_nus": state.packet_collection.output_nus.copy(),
            "output_energies": state.packet_collection.output_energies.copy(),
            "mean_intensity_total": state.estimators_bulk.mean_intensity_total,
            "mean_intensity_blueward": (
                state.estimators_line.mean_intensity_blueward
            ),
            "photo_ion_estimator": (
                state.estimators_continuum.photo_ion_estimator
            ),
            "bf_heating_estimator": (
                state.estimators_continuum.bf_heating_estimator
            ),
            "ff_heating_estimator": (
                state.estimators_continuum.ff_heating_estimator
            ),
        }

    assert workflow.transport_state.estimator_batches.packet_counts.sum() == (
        workflow.real_packet_count
    )
    for name, single_transport in results[1].items():
        np.testing.assert_allclose(
            results[4][name],
            single_transport,
            rtol=1e-12,
            atol=0.0,
            err_msg=f"batched transport changed {name}",
        )


def calculate_candidate_uncertainty(
    solver: ElectronThermalBalanceSolver,
    candidate: npt.NDArray[np.float64],
    residual_limits: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Map residual limits to a bound on each candidate value.

    Each shell's 2x2 block of the residual Jacobian at the candidate is
    estimated by forward differences of 1e-6 in the electron fraction and the
    temperature ratio of every shell at once, since shells are independent.
    The absolute inverse block applied to the limits bounds the change of
    candidate values over any residual inside the limits.
    """
    relative_step = 1e-6
    residual = solver.calculate_residual(candidate)
    shell_jacobian = np.zeros((len(candidate) // 2, 2, 2))
    for variable_idx in (0, 1):
        stepped_candidate = candidate.copy()
        stepped_candidate[variable_idx::2] *= 1.0 + relative_step
        residual_change = (
            solver.calculate_residual(stepped_candidate) - residual
        ).reshape(-1, 2)
        shell_jacobian[:, :, variable_idx] = residual_change / (
            relative_step * candidate[variable_idx::2, np.newaxis]
        )
    shell_limits = np.asarray(residual_limits).reshape(-1, 2)
    return np.concatenate(
        [
            np.abs(np.linalg.inv(jacobian)).dot(limits)
            for jacobian, limits in zip(
                shell_jacobian, shell_limits, strict=True
            )
        ]
    )


def test_noise_limited_thermal_solution_is_consistent_with_tight_root(
    first_transport_thermal_problem: FirstTransportThermalProblem,
) -> None:
    """Accept a noise-limited state within its Monte Carlo uncertainty.

    Claim: In every shell, the electron density and temperature of the
    accepted noise-limited state lie within the uncertainty that the packet
    batches imply around the tight thermal root.
    Regime: Pure-hydrogen comparison model after its first Monte Carlo
    transport, using 16 independent packet batches.
    Verification: The tight root is the same joint solve with zero limits,
    run to the optimizer's own tolerance. A finite-difference Jacobian at the
    tight root maps the 2-sigma limits, plus the residual left at the tight
    root, into candidate space; the absolute inverse Jacobian bounds the
    change for any residual inside that box.
    """
    (
        workflow,
        initial_guess,
        lower_bound,
        upper_bound,
        maximum_electron_density,
    ) = first_transport_thermal_problem
    solver = ElectronThermalBalanceSolver(
        workflow._thermal_balance_evaluator,
        workflow._thermal_balance_radiation_temperature,
        workflow._thermal_balance_level_initial_guess,
        maximum_electron_density,
        batch_evaluators=lambda: (
            workflow.iterate_batch_thermal_balance_evaluators(
                maximum_electron_density
            )
        ),
        batch_packet_counts=workflow.estimator_batch_packet_counts,
    )

    noise_candidate, _ = solver.solve(initial_guess, lower_bound, upper_bound)
    tight_candidate = solver.solve_multiple_shells(
        initial_guess,
        lower_bound,
        upper_bound,
        np.zeros_like(initial_guess),
    )

    candidate_limits = calculate_candidate_uncertainty(
        solver,
        tight_candidate,
        solver.outer_residual_limits
        + np.abs(solver.calculate_residual(tight_candidate)),
    )
    np.testing.assert_array_less(
        np.abs(noise_candidate - tight_candidate),
        candidate_limits,
        err_msg=(
            "a shell moved beyond its propagated 2-sigma Monte Carlo "
            "uncertainty"
        ),
    )


def test_thermal_solve_rejects_failed_trial_and_keeps_root(
    first_transport_thermal_problem: FirstTransportThermalProblem,
) -> None:
    """Reject a trial point where the plasma solve fails.

    Claim: When the plasma evaluator raises PlasmaIonizationError at an
    optimizer trial point, the joint solve shrinks its step and still
    reaches the same thermal root.
    Regime: Pure-hydrogen comparison model after its first Monte Carlo
    transport, solved to the optimizer's own tolerance (zero limits).
    Verification: Metamorphic. The failure is armed only after the row scale
    is measured, so it fires inside the optimizer, at the first evaluation
    that moves the temperature by more than a finite-difference step. The
    undisturbed solve gives the expected root. The Jacobian maps the
    residuals left at both roots, each plus the ~1e-11 round-off scatter of
    the residual, into a bound on their difference.
    """
    (
        workflow,
        initial_guess,
        lower_bound,
        upper_bound,
        maximum_electron_density,
    ) = first_transport_thermal_problem
    evaluator = workflow._thermal_balance_evaluator
    initial_temperature = (
        workflow._thermal_balance_radiation_temperature * initial_guess[1::2]
    )

    class FailOnFirstTrialEvaluator:
        """Once armed, raise at the first trial far from the initial point."""

        def __init__(self) -> None:
            self.armed = False
            self.failed = False

        def evaluate(
            self,
            trial_electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            level_seed: pd.DataFrame,
        ) -> object:
            # Finite-difference columns move T_e by ~1e-8; a trial step by more.
            moved = np.max(
                np.abs(electron_temperature / initial_temperature - 1.0)
            )
            if self.armed and not self.failed and moved > 1e-6:
                self.failed = True
                raise PlasmaIonizationError("injected trial failure")
            return evaluator.evaluate(
                trial_electron_density, electron_temperature, level_seed
            )

    failing_evaluator = FailOnFirstTrialEvaluator()
    solvers = {
        label: ElectronThermalBalanceSolver(
            trial_evaluator,
            workflow._thermal_balance_radiation_temperature,
            workflow._thermal_balance_level_initial_guess,
            maximum_electron_density,
            batch_evaluators=lambda: (),
            batch_packet_counts=workflow.estimator_batch_packet_counts,
        )
        for label, trial_evaluator in (
            ("undisturbed", evaluator),
            ("rejected trial", failing_evaluator),
        )
    }
    measure_row_scale = solvers["rejected trial"].calculate_row_scale

    def measure_row_scale_then_arm(
        *args: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        row_scale = measure_row_scale(*args)
        failing_evaluator.armed = True
        return row_scale

    solvers["rejected trial"].calculate_row_scale = measure_row_scale_then_arm
    candidates = {
        label: solver.solve_multiple_shells(
            initial_guess,
            lower_bound,
            upper_bound,
            np.zeros_like(initial_guess),
        )
        for label, solver in solvers.items()
    }

    assert failing_evaluator.failed, "no trial point was rejected"
    np.testing.assert_array_less(
        np.abs(candidates["rejected trial"] - candidates["undisturbed"]),
        calculate_candidate_uncertainty(
            solvers["undisturbed"],
            candidates["undisturbed"],
            sum(
                np.abs(solvers["undisturbed"].calculate_residual(candidate))
                + 1e-11
                for candidate in candidates.values()
            ),
        ),
        err_msg="a rejected trial point changed the thermal root",
    )


def test_iip_level_solution_is_seed_independent_after_first_transport(
    first_transport_thermal_problem: FirstTransportThermalProblem,
) -> None:
    """Close the same IIP level equations from two initial populations.

    Claim: A fixed shell equilibrium is independent of the initial level
    population and closes the level equations.
    Regime: Pure-hydrogen comparison model, first Monte Carlo estimators, and
    the shell-3 thermal candidate that exposed seed sensitivity after #3704.
    Verification: Seed independence is a metamorphic relation; the direct
    level-equation residual prevents two matching but unsolved populations from
    passing. The 1e-10 residual tolerance matches the workflow's declared
    level-closure limit; 1e-8 relative and 1e-10 absolute population
    tolerances allow the 1e-12-step inner solve to amplify small changes in
    weakly populated excited levels.
    """
    workflow, _, _, _, maximum_electron_density = (
        first_transport_thermal_problem
    )
    plasma = workflow.plasma_solver
    evaluator = workflow._thermal_balance_evaluator
    shell_idx = 3
    # First thermal solve candidate: n_e in cm^-3 and T_e in K. Keep full
    # precision because rounding these ratios changes the level-root basin.
    electron_density = maximum_electron_density[shell_idx] * 0.998778561236385
    electron_temperature = (
        workflow._thermal_balance_radiation_temperature[shell_idx]
        * 0.8503012782179825
    )
    temperatures = np.asarray(plasma.t_electrons, dtype=np.float64).copy()
    temperatures[shell_idx] = electron_temperature
    electron_densities = plasma.electron_densities.to_numpy(
        dtype=np.float64
    ).copy()
    electron_densities[shell_idx] = electron_density
    initial_result = evaluator.evaluate(
        electron_densities,
        temperatures,
        workflow._thermal_balance_level_initial_guess,
    )
    accepted_result = evaluator.evaluate(
        electron_densities,
        temperatures,
        initial_result.normalized_population,
    )

    np.testing.assert_allclose(
        initial_result.normalized_population[shell_idx],
        accepted_result.normalized_population[shell_idx],
        rtol=1e-8,
        atol=1e-10,
        err_msg="the IIP level equilibrium depends on its initial population",
    )
    np.testing.assert_allclose(
        initial_result.relative_trial_level_residual[shell_idx],
        0.0,
        rtol=0.0,
        atol=1e-10,
        err_msg="the initial-seed IIP level equations do not close",
    )
    np.testing.assert_allclose(
        accepted_result.relative_trial_level_residual[shell_idx],
        0.0,
        rtol=0.0,
        atol=1e-10,
        err_msg="the accepted-seed IIP level equations do not close",
    )


def test_iip_level_equations_close_from_distinct_seeds_after_second_transport(
    type_iip_workflow: TypeIIPWorkflow,
) -> None:
    """Close second-transport hydrogen level equations from distinct seeds.

    Claim: Both initial populations lead to a physical solution of the fixed
    shell-23 level equations.
    Regime: Pure-hydrogen comparison model, two 100000-packet transport steps,
    and the second thermal candidate that exposed the near-ground seed trap.
    Verification: Recompute the level-equation residual independently of the
    solver's returned residual. The 1e-10 limit is the workflow's level-closure
    criterion. Distinct closed roots are allowed in this nonlinear model.
    """
    workflow = type_iip_workflow
    workflow.show_progress_bars = False
    for iteration_idx in range(2):
        opacity_states = workflow.solve_opacity()
        workflow.solve_montecarlo(opacity_states, workflow.real_packet_count)
        estimated_values, _ = workflow.get_convergence_estimates()
        workflow.solve_simulation_state(estimated_values)
        continuum_estimators, j_blues = workflow.update_estimators()
        workflow.damped_normalized_j_blues = j_blues
        workflow.solve_plasma(continuum_estimators, j_blues)
        if iteration_idx == 0:
            workflow.solve_thermal_balance()
            workflow.solve_continuum_state(continuum_estimators)
            workflow.converged = workflow.check_convergence(estimated_values)
            workflow.completed_iterations += 1

    plasma = workflow.plasma_solver
    maximum_electron_density = (
        plasma.number_density.multiply(
            plasma.number_density.index.to_numpy(), axis=0
        )
        .sum()
        .to_numpy()
    )
    workflow._initialize_thermal_balance_evaluator(maximum_electron_density)
    evaluator = workflow._thermal_balance_evaluator
    shell_idx = 23
    # These dimensionless candidate ratios reproduce the second thermal
    # iteration's open shell-23 level solution.
    electron_density = maximum_electron_density[shell_idx] * 0.69998
    electron_temperature = (
        workflow._thermal_balance_radiation_temperature[shell_idx] * 0.79428
    )
    temperatures = np.asarray(plasma.t_electrons, dtype=np.float64).copy()
    temperatures[shell_idx] = electron_temperature
    continuum_rates = evaluator.calculate_continuum_coefficients(temperatures)[
        0
    ][shell_idx]
    original_seed = workflow._thermal_balance_level_initial_guess.iloc[
        :, shell_idx
    ].to_numpy(dtype=np.float64)
    uniform_seed = np.full_like(original_seed, 1.0 / len(original_seed))

    original_solution = evaluator._calculate_level_solution(
        shell_idx,
        electron_density,
        electron_temperature,
        original_seed,
        continuum_rates,
        evaluator.shell_number_densities[shell_idx],
        evaluator.sobolev_inputs[shell_idx],
    )
    uniform_solution = evaluator._calculate_level_solution(
        shell_idx,
        electron_density,
        electron_temperature,
        uniform_seed,
        continuum_rates,
        evaluator.shell_number_densities[shell_idx],
        evaluator.sobolev_inputs[shell_idx],
    )

    for seed_name, solution in (
        ("uniform", uniform_solution),
        ("original", original_solution),
    ):
        fractions = solution[0]
        assert np.all(fractions >= 0.0)
        np.testing.assert_allclose(fractions.sum(), 1.0, rtol=0.0, atol=1e-12)
        physical_residual = evaluator._calculate_level_state(
            shell_idx,
            electron_density,
            electron_temperature,
            fractions,
            continuum_rates,
            evaluator.shell_number_densities[shell_idx],
            evaluator.sobolev_inputs[shell_idx],
        )[0]
        np.testing.assert_allclose(
            physical_residual,
            0.0,
            rtol=0.0,
            atol=1e-10,
            err_msg=f"the {seed_name} seed does not close the hydrogen equations",
        )


def test_iip_documentation_model_closes_every_shell_at_noise_limit(
    tardis_regression_path: Path,
) -> None:
    """Close every shell of the documented model at its Monte Carlo limits.

    Claim: The first thermal solve accepts a state in which every shell
    balances charge and heating to within twice its Monte Carlo standard
    errors.
    Regime: The ten-shell pure-hydrogen documentation model after its first
    200000-packet Monte Carlo transport step.
    Verification: Orchestration check of the solver against its own stopping
    rule: the electron and heating residuals of the accepted state are
    recomputed and compared with the limits measured from the transport's
    independent packet batches. The level and charge closure are checked by
    the solve itself, which raises if they fail.
    """
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/simple_iip_config.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    config.plasma.nlte.species = [(1, 0)]
    workflow = TypeIIPWorkflow(config)
    workflow.show_progress_bars = False
    workflow.solve_montecarlo(
        workflow.solve_opacity(), workflow.real_packet_count
    )
    estimated_values, _ = workflow.get_convergence_estimates()
    workflow.solve_simulation_state(estimated_values)
    continuum_estimators, j_blues = workflow.update_estimators()
    workflow.damped_normalized_j_blues = j_blues
    workflow.solve_plasma(continuum_estimators, j_blues)
    workflow.solve_thermal_balance()

    solver = workflow.thermal_balance_solver
    closed_shells = solver.find_converged_shells(
        solver.calculate_electron_and_thermal_residuals(
            workflow._thermal_balance_evaluation
        ),
        solver.outer_residual_limits,
    )
    assert closed_shells.all(), "a shell did not reach its Monte Carlo limits"


@pytest.fixture(scope="module")
def hydrogen_helium_config(tardis_regression_path: Path) -> Configuration:
    """Return the three-shell hydrogen-helium Type IIP configuration."""
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/iip_hydrogen_helium.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    return config


@pytest.fixture(scope="module")
def hydrogen_helium_workflow(
    hydrogen_helium_config: Configuration,
) -> TypeIIPWorkflow:
    """Run the three-shell hydrogen-helium Type IIP model once."""
    workflow = TypeIIPWorkflow(deepcopy(hydrogen_helium_config))
    workflow.show_progress_bars = False
    workflow.run()
    return workflow


@pytest.fixture(scope="module")
def hydrogen_helium_dilute_lte_level_fraction(
    hydrogen_helium_workflow: TypeIIPWorkflow,
) -> pd.DataFrame:
    """Return dilute-LTE level fractions at the accepted radiation field.

    The fraction of ion (Z, i) in level l is
    g_l W^(1 - m_l) exp(-E_l / k T_rad) / sum_l [g_l W^(1 - m_l)
    exp(-E_l / k T_rad)], where m_l is 1 for metastable levels.
    """
    plasma = hydrogen_helium_workflow.plasma_solver
    level_index = plasma.level_number_density.index
    levels = plasma.atomic_data.levels.loc[level_index]
    radiation_temperature = np.asarray(plasma.t_rad, dtype=np.float64)
    dilution_factor = np.asarray(plasma.w, dtype=np.float64)
    boltzmann_factor = levels.g.to_numpy()[:, None] * np.exp(
        -levels.energy.to_numpy()[:, None]
        / (const.k_B.cgs.value * radiation_temperature)
    )
    metastable = plasma.metastability.loc[level_index].to_numpy()[:, None]
    boltzmann_factor = np.where(
        metastable, boltzmann_factor, boltzmann_factor * dilution_factor
    )
    boltzmann_factor = pd.DataFrame(
        boltzmann_factor,
        index=level_index,
        columns=plasma.level_number_density.columns,
    )
    return boltzmann_factor / boltzmann_factor.groupby(
        level=["atomic_number", "ion_number"]
    ).transform("sum")


def test_iip_hydrogen_helium_workflow_closes_thermal_balance(
    hydrogen_helium_workflow: TypeIIPWorkflow,
) -> None:
    """Run a Type IIP model containing helium to its final iteration.

    Claim: A hydrogen-helium model completes the Type IIP iteration with only
    the atomic data in the configured atom file, and its accepted
    thermal-balance state closes every shell at its Monte Carlo limits.
    Regime: Three shells, 70 % hydrogen and 30 % helium by mass, with H I as
    the only continuum and NLTE species.
    Verification: Orchestration check on the run; the closed shells are
    recomputed from the accepted residuals and the measured limits.
    """
    workflow = hydrogen_helium_workflow
    solver = workflow.thermal_balance_solver
    closed_shells = solver.find_converged_shells(
        solver.calculate_electron_and_thermal_residuals(
            workflow._thermal_balance_evaluation
        ),
        solver.outer_residual_limits,
    )
    assert closed_shells.all(), "a shell did not reach its Monte Carlo limits"
    assert workflow.completed_iterations == 1
    output_energies = workflow.transport_state.packet_collection.output_energies
    assert len(output_energies) == workflow.final_iteration_packet_count


def test_iip_hydrogen_helium_hydrogen_levels_follow_nlte_rate_equations(
    hydrogen_helium_workflow: TypeIIPWorkflow,
    hydrogen_helium_dilute_lte_level_fraction: pd.DataFrame,
) -> None:
    """Solve H I level populations from the NLTE rate equations.

    Claim: H I is the only species whose level populations are solved from
    the NLTE rate equations, and those populations close the equations
    instead of following dilute-LTE excitation.
    Regime: Accepted final thermal-balance state of the hydrogen-helium model.
    Verification: The plasma H I level fractions must equal the solved NLTE
    fractions, the level equations must close to the workflow's 1e-10 of
    each equation's gross rate,
    and at least one H I level must depart from its directly evaluated
    dilute-LTE fraction by more than 10 %. The 1e-12 relative tolerance allows
    only floating-point differences in the population rebuild.
    """
    plasma = hydrogen_helium_workflow.plasma_solver
    evaluation = hydrogen_helium_workflow._thermal_balance_evaluation
    nlte_species = evaluation.normalized_population.index.droplevel(
        "level_number"
    ).unique()
    assert nlte_species.tolist() == [(1, 0)]

    hydrogen_level_fraction = (
        plasma.level_number_density.loc[(1, 0)]
        / plasma.ion_number_density.loc[(1, 0)]
    )
    pd.testing.assert_frame_equal(
        hydrogen_level_fraction,
        evaluation.normalized_population.loc[(1, 0)],
        rtol=1e-12,
        atol=0.0,
        check_names=False,
    )
    np.testing.assert_allclose(
        evaluation.relative_trial_level_residual,
        0.0,
        rtol=0.0,
        atol=1e-10,
        err_msg="the H I level populations do not close the NLTE equations",
    )
    departure_from_dilute_lte = (
        hydrogen_level_fraction
        / hydrogen_helium_dilute_lte_level_fraction.loc[(1, 0)]
    )
    assert (np.abs(departure_from_dilute_lte - 1.0) > 0.1).to_numpy().any(), (
        "H I level populations follow dilute-LTE excitation"
    )


def test_iip_hydrogen_helium_helium_ionization_is_nebular(
    hydrogen_helium_workflow: TypeIIPWorkflow,
) -> None:
    """Ionize helium with the nebular approximation of the reference code.

    Claim: Helium, which has no continuum rates, follows the C-TARDIS
    treatment of non-NLTE species: n(He II) n_e / n(He I) equals the nebular
    ionization factor of Mazzali & Lucy (1993) at the accepted T_rad, W, and
    T_e, with n_e set by charge conservation.
    Regime: Accepted final thermal-balance state of the hydrogen-helium model.
    Verification: The plasma graph evaluates the nebular factor independently
    of the equilibrium evaluator that solves the ion populations. The 2e-8
    relative tolerance is the workflow's electron-density closure limit.
    """
    plasma = hydrogen_helium_workflow.plasma_solver
    helium_ion_density = plasma.ion_number_density.loc[2]
    ionization_ratio = (
        helium_ion_density.loc[1]
        * plasma.electron_densities
        / helium_ion_density.loc[0]
    )
    np.testing.assert_allclose(
        ionization_ratio,
        plasma.phi.loc[(2, 1)],
        rtol=2e-8,
        atol=0.0,
        err_msg="He I/He II does not follow nebular ionization",
    )


def test_iip_hydrogen_helium_helium_excitation_is_dilute_lte(
    hydrogen_helium_workflow: TypeIIPWorkflow,
    hydrogen_helium_dilute_lte_level_fraction: pd.DataFrame,
) -> None:
    """Excite helium in dilute LTE at the accepted radiation field.

    Claim: He I and He II level populations follow the C-TARDIS dilute-LTE
    excitation at the radiation temperature and dilution factor used for
    the final thermal-balance solve.
    Regime: Accepted final thermal-balance state of the hydrogen-helium model.
    Verification: Level fractions are evaluated directly from level energies,
    statistical weights, metastability, T_rad, and W. The 1e-12 relative
    tolerance allows only floating-point differences.
    """
    plasma = hydrogen_helium_workflow.plasma_solver
    for species in ((2, 0), (2, 1)):
        helium_level_fraction = (
            plasma.level_number_density.loc[species]
            / plasma.ion_number_density.loc[species]
        )
        pd.testing.assert_frame_equal(
            helium_level_fraction,
            hydrogen_helium_dilute_lte_level_fraction.loc[species],
            rtol=1e-12,
            atol=0.0,
            check_names=False,
            obj=f"{species} level fractions",
        )


def test_iip_hydrogen_helium_helium_populations_match_ctardis(
    hydrogen_helium_config: Configuration,
    ctardis_reference_path: Path,
) -> None:
    """Reproduce C-TARDIS helium populations in its radiation field.

    Claim: At the radiation field, electron temperature, and electron density
    of a C-TARDIS run of the same hydrogen-helium model, the equilibrium
    evaluator gives the C-TARDIS helium ionization ratio and level fractions.
    Regime: Three-shell model after one C-TARDIS thermal-balance solve, with
    hydrogen continuum rates from one TARDIS Monte Carlo step.
    Verification: C-TARDIS is an independent implementation. Fixing its
    radiation field removes the Monte Carlo differences between the codes;
    the helium nebular ionization factor and dilute-LTE excitation depend on
    T_rad, W, and T_e but not on the hydrogen continuum estimators.
    The 5e-5 relative tolerance is above the 1.6e-5 offset between the
    C-TARDIS level and ion totals, and far below the factor-of-100 change
    from LTE helium ionization.
    """
    radiation_temperature, dilution_factor, electron_temperature = (
        pd.read_hdf(
            ctardis_reference_path / f"ctardis_hydrogen_helium_{name}.h5",
            key="data",
        ).to_numpy()
        for name in ("t_rad", "w", "t_electrons")
    )
    ctardis_electron_density = pd.read_hdf(
        ctardis_reference_path
        / "ctardis_hydrogen_helium_electron_densities.h5",
        key="data",
    ).to_numpy()
    ctardis_ion_density = pd.read_hdf(
        ctardis_reference_path / "ctardis_hydrogen_helium_ion_density.h5",
        key="data",
    )
    ctardis_level_density = pd.read_hdf(
        ctardis_reference_path
        / "ctardis_hydrogen_helium_level_number_density.h5",
        key="data",
    )

    workflow = TypeIIPWorkflow(deepcopy(hydrogen_helium_config))
    workflow.show_progress_bars = False
    workflow.solve_montecarlo(
        workflow.solve_opacity(), workflow.real_packet_count
    )
    continuum_estimators, j_blues = workflow.update_estimators()
    workflow.solve_plasma(continuum_estimators, j_blues)
    plasma = workflow.plasma_solver
    plasma.update(
        dilute_planckian_radiation_field=DilutePlanckianRadiationField(
            radiation_temperature * u.K, dilution_factor
        ),
        link_t_rad_t_electron=electron_temperature / radiation_temperature,
    )
    maximum_electron_density = (
        plasma.number_density.multiply(
            plasma.number_density.index.to_numpy(), axis=0
        )
        .sum()
        .to_numpy()
    )
    workflow._initialize_thermal_balance_evaluator(maximum_electron_density)
    evaluation = workflow._thermal_balance_evaluator.evaluate(
        ctardis_electron_density,
        electron_temperature,
        workflow._thermal_balance_level_initial_guess,
    )

    helium_ion_density = evaluation.ion_population.loc[2]
    ionization_ratio = (
        helium_ion_density.loc[1]
        * evaluation.charge_solved_electron_density
        / helium_ion_density.loc[0]
    )
    ctardis_ionization_ratio = (
        ctardis_ion_density.loc[(2, 1)]
        * ctardis_electron_density
        / ctardis_ion_density.loc[(2, 0)]
    )
    np.testing.assert_allclose(
        ionization_ratio,
        ctardis_ionization_ratio,
        rtol=5e-5,
        atol=0.0,
        err_msg="He I/He II differs from C-TARDIS in the same field",
    )
    for species in ((2, 0), (2, 1)):
        level_fraction = (
            evaluation.absolute_level_population.loc[species]
            / evaluation.ion_population.loc[species]
        )
        ctardis_level_fraction = (
            ctardis_level_density.loc[species]
            / ctardis_ion_density.loc[species]
        )
        pd.testing.assert_frame_equal(
            level_fraction,
            ctardis_level_fraction,
            rtol=5e-5,
            atol=0.0,
            check_names=False,
            obj=f"{species} level fractions",
        )


def test_solve_montecarlo(
    type_iip_workflow: TypeIIPWorkflow,
    regression_data: RegressionData,
) -> None:
    """Preserve the Type IIP emergent luminosity after standard initialization.

    The ``1e-5`` relative tolerance is the plan-wide legacy-parity allowance
    for the non-bitwise-identical standard plasma bootstrap.
    """
    opacity_states = type_iip_workflow.solve_opacity()
    type_iip_workflow.solve_montecarlo(opacity_states, 1000)
    type_iip_workflow.initialize_spectrum_solver()
    luminosity_density = type_iip_workflow.spectrum_solver.spectrum_real_packets.luminosity_density_lambda.value
    expected_luminosity_density = regression_data.sync_ndarray(
        luminosity_density
    )
    np.testing.assert_allclose(
        luminosity_density,
        expected_luminosity_density,
        atol=0,
        # The standard-plasma bootstrap is not bitwise identical to the
        # removed legacy IIP owner. Use the plan's general legacy-parity
        # contract rather than tuning this threshold to one packet sample.
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    "mode,include_thermal_pool,include_ionization_pool",
    [
        ("full", False, False),
        ("partial", False, False),
        ("partial", True, False),
        ("partial", False, True),
        ("partial", True, True),
        ("explicit", False, False),
    ],
)
def test_iip_outer_shell_population_cutoff_second_iteration_opacity(
    tardis_regression_path: Path,
    mode: str,
    include_thermal_pool: bool,
    include_ionization_pool: bool,
) -> None:
    """Evaluate finite outer-shell opacity at the 1500 K thermal floor.

    Claim: A zero LTE hydrogen-ion population produces the stimulated-
    recombination correction in ``chi_bf`` without invalid opacity values.
    Regime: Second opacity iteration in shells forced to the 1500 K floor.
    Verification: The expected opacity is evaluated directly from the
    bound-free population equation, independently of continuum-state assembly.
    All reduction modes retain legacy signed-coefficient compatibility and
    complete ten transport packets; completion is not probability validation.
    """
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/iip_population_cutoff.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    config.plasma.continuum_interaction.markov_chain.mode = mode
    config.plasma.continuum_interaction.markov_chain.elements = (
        ["H"] if mode == "partial" else []
    )
    config.plasma.continuum_interaction.markov_chain.include_thermal_pool = (
        include_thermal_pool
    )
    config.plasma.continuum_interaction.markov_chain.include_ionization_pool = (
        include_ionization_pool
    )
    workflow = TypeIIPWorkflow(config)
    workflow.show_progress_bars = False

    plasma_solver = workflow.plasma_solver
    # The reported outer shells hit the 1500 K thermal-balance floor before
    # their second opacity solve. Force that state directly so this regression
    # does not need the slow least-squares thermal solve.
    forced_t_electrons = np.full(
        workflow.simulation_state.geometry.no_of_shells_active, 1500.0
    )
    plasma_solver.update(
        link_t_rad_t_electron=(
            forced_t_electrons / np.asarray(plasma_solver.t_rad)
        ),
    )
    maximum_electron_density = (
        plasma_solver.number_density.multiply(
            plasma_solver.number_density.index.to_numpy(), axis=0
        )
        .sum()
        .to_numpy()
    )
    evaluator = workflow._build_thermal_balance_evaluator(
        maximum_electron_density, analytic=True
    )
    (
        _,
        level_to_continuum_saha_factor,
        _,
        _,
        coefficient_state,
        partition_function,
        level_boltzmann_factor,
    ) = evaluator.calculate_continuum_coefficients(forced_t_electrons)

    lte_ion_population, _ = calculate_lte_populations(
        plasma_solver.thermal_phi_lte,
        partition_function,
        plasma_solver.number_density,
        plasma_solver.electron_densities,
        level_boltzmann_factor,
        plasma_solver.atomic_data.levels.loc[
            plasma_solver.level_number_density.index
        ],
    )
    workflow._build_continuum_states(
        coefficient_state,
        level_to_continuum_saha_factor,
    )
    assert lte_ion_population.loc[(1, 1)].iloc[-1] == 0.0

    workflow.completed_iterations = 1
    opacity_states = workflow.solve_opacity()
    continuum_state = workflow.continuum_opacity_state
    assert np.isfinite(continuum_state.chi_bf.to_numpy()).all()
    cross_sections = plasma_solver.photo_ion_cross_sections
    upper_ion_index = pd.MultiIndex.from_arrays(
        [
            cross_sections.index.get_level_values("atomic_number"),
            cross_sections.index.get_level_values("ion_number") + 1,
        ],
        names=["atomic_number", "ion_number"],
    )
    stimulated_recombination_population = (
        level_to_continuum_saha_factor.loc[cross_sections.index].to_numpy()
        * plasma_solver.ion_number_density.loc[upper_ion_index].to_numpy()
        * plasma_solver.electron_densities.to_numpy()
    )
    boltzmann_factor = np.exp(
        -cross_sections.nu.to_numpy()[:, None]
        / forced_t_electrons
        * (const.h.cgs.value / const.k_B.cgs.value)
    )
    # chi_bf = [n_l - n_e n_(ion+1) Phi_lu exp(-h nu / k_B T_e)] sigma_bf.
    expected_chi_bf = (
        plasma_solver.level_number_density.loc[cross_sections.index]
        - stimulated_recombination_population * boltzmann_factor
    ).multiply(cross_sections.x_sect.to_numpy(), axis=0)
    pd.testing.assert_frame_equal(
        continuum_state.chi_bf,
        expected_chi_bf.loc[continuum_state.level2continuum_idx.index],
    )
    assert np.isfinite(continuum_state.p_fb_deactivation.to_numpy()).all()
    assert np.isfinite(continuum_state.emissivities.to_numpy()).all()
    workflow.solve_montecarlo(opacity_states, 10)
    assert len(workflow.transport_state.packet_collection.output_energies) == 10
