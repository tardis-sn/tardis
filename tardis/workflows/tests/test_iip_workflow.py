from pathlib import Path
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from astropy import units as u
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis import constants as const
from tardis.conftest import assert_regression_dataframe
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.equilibrium.evaluator import calculate_lte_populations
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


def test_thermal_balance_iteration_delegates_to_evaluator() -> None:
    """Map a two-shell thermal-balance candidate to evaluator quantities.

    Claim: Each shell's density fraction and temperature ratio are converted
    to physical values, and the two balance errors retain shell order.
    Regime: Two shells with unequal density limits and radiation temperatures.
    Verification: The expected values follow by direct multiplication and
    inserting the recorded evaluator result.
    """
    # Skip full workflow construction because this check needs only the state
    # read by thermal_balance_iteration, not a configured simulation.
    workflow = TypeIIPWorkflow.__new__(TypeIIPWorkflow)
    # SimpleNamespace provides the sole plasma value used here: each shell's
    # radiation temperature.
    workflow.plasma_solver = SimpleNamespace(t_rad=np.array([1.0e4, 2.0e4]))
    level_initial_guess = pd.DataFrame([[0.6, 0.7], [0.4, 0.3]])
    workflow._thermal_balance_level_initial_guess = level_initial_guess

    # Record the physical values received and return fixed balance errors so
    # this test checks scaling and shell order without running a plasma solve.
    class RecordingEvaluator:
        def __init__(self) -> None:
            self.call_count = 0
            self.arguments: tuple[object, ...] | None = None
            # SimpleNamespace contains only the two evaluator results consumed
            # by thermal_balance_iteration.
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
    workflow._thermal_balance_evaluator = evaluator
    workflow._thermal_balance_radiation_temperature = np.array([1.0e4, 2.0e4])
    candidate = np.array([0.25, 0.8, 0.5, 1.1])
    maximum_electron_density = np.array([4.0e9, 6.0e9])

    residual = workflow.thermal_balance_iteration(
        candidate, maximum_electron_density
    )

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
    assert not hasattr(workflow, "_thermal_balance_evaluation")


def test_thermal_balance_uses_joint_solution_without_local_retry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep the joint thermal solution when it passes final validation.

    Claim: Shell-local optimization is a fallback and does not add optimizer
    calls to a closed joint solution.
    Regime: One shell with a joint candidate inside the physical bounds.
    Verification: Record calls at the SciPy optimizer boundary and require the
    accepted plasma candidate to equal the joint result.
    """
    workflow = TypeIIPWorkflow.__new__(TypeIIPWorkflow)
    workflow.plasma_solver = SimpleNamespace(
        link_t_rad_t_electron=np.array([0.8]),
        number_density=pd.DataFrame([[2.0]], index=[1]),
        electron_densities=pd.Series([1.0]),
    )
    workflow.simulation_state = SimpleNamespace(
        geometry=SimpleNamespace(no_of_shells_active=1),
        t_radiative=np.array([1.0e4]) * u.K,
    )
    workflow._thermal_balance_radiation_temperature = np.array([1.0e4])
    workflow._thermal_balance_level_initial_guess = pd.DataFrame([[1.0]])
    evaluation = SimpleNamespace(normalized_population=pd.DataFrame([[1.0]]))
    workflow._thermal_balance_evaluator = SimpleNamespace(
        evaluate=lambda *_: evaluation
    )
    workflow._initialize_thermal_balance_evaluator = lambda _: None
    workflow._validate_thermal_balance_evaluation = lambda _: None
    accepted_candidates = []
    workflow._update_plasma_with_thermal_balance_state = lambda candidate: (
        accepted_candidates.append(candidate.copy())
    )
    optimizer_starts = []
    joint_candidate = np.array([0.6, 0.9])

    def recording_lsq(
        _objective: object, initial_guess: npt.ArrayLike, **_kwargs: object
    ) -> SimpleNamespace:
        optimizer_starts.append(np.asarray(initial_guess).copy())
        return SimpleNamespace(x=joint_candidate, fun=np.zeros(2), status=1)

    monkeypatch.setattr("tardis.workflows.type_iip_workflow.lsq", recording_lsq)

    workflow.solve_thermal_balance()

    assert len(optimizer_starts) == 1
    np.testing.assert_array_equal(accepted_candidates, [joint_candidate])


def test_thermal_balance_retries_failed_shell_from_alternate_start(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Recover a failed shell from an alternate bounded thermal state.

    Claim: A joint candidate that fails closure is retried locally, and a
    closing alternate is preferred over a lower-norm candidate that still
    violates an individual residual gate.
    Regime: Two shells, with the first joint electron fraction below the
    positive local bound and its first local result narrowly failing electron
    closure while the second shell is already closed.
    Verification: Record optimizer starts and supply residuals at the declared
    electron and heating limits independently of the workflow implementation.
    """
    workflow = TypeIIPWorkflow.__new__(TypeIIPWorkflow)
    workflow.plasma_solver = SimpleNamespace(
        link_t_rad_t_electron=np.array([0.8, 0.9]),
        number_density=pd.DataFrame([[2.0, 4.0]], index=[1]),
        electron_densities=pd.Series([1.0, 2.0]),
    )
    workflow.simulation_state = SimpleNamespace(
        geometry=SimpleNamespace(no_of_shells_active=2),
        t_radiative=np.array([1.0e4, 1.1e4]) * u.K,
    )
    workflow._thermal_balance_radiation_temperature = np.array([1.0e4, 1.1e4])
    workflow._thermal_balance_level_initial_guess = pd.DataFrame([[1.0, 1.0]])
    failed_evaluation = SimpleNamespace(
        normalized_population=pd.DataFrame([[1.0, 1.0]]),
        electron_residual=pd.Series([0.7, 0.0]),
        fractional_heating=pd.Series([0.06, 0.0]),
        total_heating=pd.Series([1.0e-8, 0.0]),
    )
    closed_evaluation = SimpleNamespace(
        normalized_population=pd.DataFrame([[1.0, 1.0]]),
        electron_residual=pd.Series([0.0, 0.0]),
        fractional_heating=pd.Series([0.0, 0.0]),
        total_heating=pd.Series([0.0, 0.0]),
    )
    evaluator_calls = 0

    def evaluate(*_args: object) -> SimpleNamespace:
        nonlocal evaluator_calls
        evaluator_calls += 1
        return failed_evaluation if evaluator_calls <= 2 else closed_evaluation

    workflow._thermal_balance_evaluator = SimpleNamespace(evaluate=evaluate)
    workflow._initialize_thermal_balance_evaluator = lambda _: None
    validation_calls = 0

    def validate(_evaluation: object) -> None:
        nonlocal validation_calls
        validation_calls += 1
        if validation_calls == 1:
            raise RuntimeError("joint closure failed")

    workflow._validate_thermal_balance_evaluation = validate
    accepted_candidates = []
    workflow._update_plasma_with_thermal_balance_state = lambda candidate: (
        accepted_candidates.append(candidate.copy())
    )
    optimizer_starts = []
    optimizer_results = (
        SimpleNamespace(
            x=np.array([0.0, 0.7, 0.7, 1.0]),
            fun=np.array([0.7, 0.06, 0.0, 0.0]),
            status=3,
        ),
        SimpleNamespace(
            x=np.array([0.4, 0.75]), fun=np.array([3.0e-8, 0.0]), status=3
        ),
        SimpleNamespace(
            x=np.array([0.6, 0.85]), fun=np.array([0.0, 1.0e-7]), status=1
        ),
    )

    def recording_lsq(
        _objective: object, initial_guess: npt.ArrayLike, **_kwargs: object
    ) -> SimpleNamespace:
        optimizer_starts.append(np.asarray(initial_guess).copy())
        return optimizer_results[len(optimizer_starts) - 1]

    monkeypatch.setattr("tardis.workflows.type_iip_workflow.lsq", recording_lsq)

    workflow.solve_thermal_balance()

    assert len(optimizer_starts) == 3
    np.testing.assert_array_equal(
        optimizer_starts[1], [np.finfo(float).eps, 0.7]
    )
    np.testing.assert_array_equal(optimizer_starts[2], [0.5, 0.8])
    np.testing.assert_array_equal(accepted_candidates, [[0.6, 0.85, 0.7, 1.0]])


def test_iip_level_solution_is_seed_independent_after_first_transport(
    type_iip_workflow: TypeIIPWorkflow,
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
    workflow = type_iip_workflow
    workflow.show_progress_bars = False
    workflow.solve_montecarlo(
        workflow.solve_opacity(), workflow.real_packet_count
    )
    estimated_values, _ = workflow.get_convergence_estimates()
    workflow.solve_simulation_state(estimated_values)
    continuum_estimators, j_blues = workflow.update_estimators()
    workflow.solve_plasma(continuum_estimators, j_blues)

    plasma = workflow.plasma_solver
    maximum_electron_density = (
        plasma.number_density.multiply(
            plasma.number_density.index.values, axis=0
        )
        .sum()
        .to_numpy()
    )
    workflow._initialize_thermal_balance_evaluator(maximum_electron_density)
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
        initial_result.trial_level_residual[shell_idx],
        0.0,
        rtol=0.0,
        atol=1e-10,
        err_msg="the initial-seed IIP level equations do not close",
    )
    np.testing.assert_allclose(
        accepted_result.trial_level_residual[shell_idx],
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
            plasma.number_density.index.values, axis=0
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


def test_iip_documentation_model_closes_first_thermal_balance(
    tardis_regression_path: Path,
) -> None:
    """Require charge and energy closure for the documented Type IIP model.

    Claim: The first thermal solve satisfies electron density and heating
    balance in every shell.
    Regime: The ten-shell pure-hydrogen documentation model after its first
    200000-packet Monte Carlo transport step.
    Verification: Compare the independently charge-solved electron density
    with the trial density and inspect the physical heating residuals after
    the workflow solve. The limits are the declared final closure thresholds;
    the 5e-13 heating floor allows cancellation of large thermal rates.
    """
    config = Configuration.from_yaml("../data/simple_iip_config.yml")
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

    evaluation = workflow._thermal_balance_evaluation
    np.testing.assert_allclose(
        evaluation.trial_electron_density,
        evaluation.charge_solved_electron_density,
        rtol=2e-8,
        atol=0.0,
        err_msg="the first thermal solve does not close electron density",
    )
    np.testing.assert_allclose(
        evaluation.total_heating,
        0.0,
        rtol=0.0,
        atol=5e-13,
        err_msg="the first thermal solve does not close total heating",
    )
    np.testing.assert_allclose(
        evaluation.fractional_heating,
        0.0,
        rtol=0.0,
        atol=2e-7,
        err_msg="the first thermal solve does not close fractional heating",
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


def test_iip_outer_shell_population_cutoff_second_iteration_opacity(
    tardis_regression_path: Path,
) -> None:
    """Evaluate finite outer-shell opacity at the 1500 K thermal floor.

    Claim: A zero LTE hydrogen-ion population produces the stimulated-
    recombination correction in ``chi_bf`` without invalid opacity values.
    Regime: Second opacity iteration in shells forced to the 1500 K floor.
    Verification: The expected opacity is evaluated directly from the
    bound-free population equation, independently of continuum-state assembly.
    """
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/iip_population_cutoff.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
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
            plasma_solver.number_density.index.values, axis=0
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
    assert np.isfinite(continuum_state.chi_bf.values).all()
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
    assert np.isfinite(continuum_state.p_fb_deactivation.values).all()
    assert np.isfinite(continuum_state.emissivities.values).all()
    workflow.solve_montecarlo(opacity_states, 10)
    assert len(workflow.transport_state.packet_collection.output_energies) == 10
