from pathlib import Path
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
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
    luminosity_density = (
        type_iip_workflow.spectrum_solver.spectrum_real_packets.luminosity_density_lambda.value
    )
    expected_luminosity_density = regression_data.sync_ndarray(luminosity_density)
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
    continuum_coefficients = evaluator.calculate_continuum_coefficients(
        forced_t_electrons
    )
    level_to_continuum_saha_factor = continuum_coefficients[1]
    partition_function = continuum_coefficients[5]
    level_boltzmann_factor = continuum_coefficients[6]
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
        continuum_coefficients,
        level_to_continuum_saha_factor,
    )
    assert lte_ion_population.loc[(1, 1)].iloc[-1] == 0.0

    workflow.completed_iterations = 1
    opacity_states = workflow.solve_opacity()
    continuum_state = opacity_states["opacity_state"].continuum_state
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
