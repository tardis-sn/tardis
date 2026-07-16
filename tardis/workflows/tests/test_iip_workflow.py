from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

from tardis.conftest import assert_regression_dataframe
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow

from collections.abc import Callable
from typing import Any


def _max_rel_diff(actual, expected):
    """Helper function to print relative diffs for checking broken tests"""
    actual_vals = actual.values
    expected_vals = expected.values

    relative_difference = np.abs(actual_vals - expected_vals) / np.abs(
        expected_vals
    )

    return float(np.nanmax(relative_difference))


def _as_flat_values(value: Any) -> np.ndarray:
    """Return plasma output values independent of pandas container shape."""
    values = value.values if hasattr(value, "values") else np.asarray(value)
    return np.asarray(values).reshape(-1)


def _assert_values_allclose(
    actual: Any,
    expected: Any,
    *,
    rtol: float,
    atol: float = 0.0,
) -> None:
    """Compare output values while allowing scalar/container shape differences."""
    actual_values = _as_flat_values(actual)
    expected_values = _as_flat_values(expected)
    if actual_values.size == 1:
        actual_values = np.full_like(expected_values, actual_values.item())
    elif expected_values.size == 1:
        expected_values = np.full_like(actual_values, expected_values.item())
    assert actual_values.size == expected_values.size
    np.testing.assert_allclose(
        actual_values,
        expected_values,
        rtol=rtol,
        atol=atol,
    )


def _as_regression_dataframe(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, pd.Series):
        return value.to_frame("value")

    values = np.asarray(value)
    if values.ndim == 1:
        return pd.DataFrame({"value": values})
    return pd.DataFrame(values)


def _assert_regression_dataframe(
    regression_data: Any,
    key: str,
    actual: Any,
    *,
    rtol: float = 1e-12,
    atol: float = 0.0,
) -> None:
    actual_frame = _as_regression_dataframe(actual)
    expected_frame = regression_data.sync_dataframe(actual_frame, key=key)
    pd.testing.assert_frame_equal(
        actual_frame,
        expected_frame,
        rtol=rtol,
        atol=atol,
        check_dtype=False,
        check_names=False,
    )


PLASMA_SOLVER_REGRESSION_OUTPUTS = (
    "electron_densities",
    "t_electrons",
    "link_t_rad_t_electron",
    "p_fb_deactivation",
    "chi_bf",
    "sp_fb_cooling_rates",
    "stimulated_emission_factor",
    "b",
    "ion_ratio",
    "j_blues",
)


INITIAL_PLASMA_SOLVER_REGRESSION_OUTPUTS = (
    "transition_probabilities",
    "ion_number_density",
    "tau_sobolevs",
    "beta_sobolev",
    "level_number_density",
    *PLASMA_SOLVER_REGRESSION_OUTPUTS,
)


@pytest.fixture
def iip_regression_path(tardis_regression_path):
    return tardis_regression_path / "tardis" / "workflows" / "tests"


@pytest.fixture
def thermal_balance_guess() -> Callable[[Any], tuple[np.ndarray, np.ndarray]]:
    def build_guess(plasma_solver: Any) -> tuple[np.ndarray, np.ndarray]:
        max_electron_number_density = (
            plasma_solver.number_density.multiply(
                plasma_solver.number_density.index.values,
                axis=0,
            )
            .sum()
            .values
        )
        electron_fraction = (
            plasma_solver.electron_densities / max_electron_number_density
        ).values

        link = np.asarray(plasma_solver.link_t_rad_t_electron)
        if link.ndim == 0:
            link = np.full_like(plasma_solver.electron_densities, link)
        guess = np.zeros(2 * len(link))
        guess[::2] = electron_fraction
        guess[1::2] = link

        return guess, max_electron_number_density

    return build_guess


@pytest.fixture
def ctardis_compare_config(tardis_regression_path):
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )

    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )

    config.plasma.nlte.species = [
        (1, 0)
    ]  # Hack to force config necessary for ctardis plasma
    return config


@pytest.fixture
def type_iip_workflow(ctardis_compare_config):
    workflow = TypeIIPWorkflow(ctardis_compare_config)
    return workflow


# identical to ctardis values
@pytest.fixture
def elemental_number_density(iip_regression_path):
    elemental_number_density = pd.read_hdf(
        iip_regression_path / "ctardis_elemental_number_density.h5",
        key="data",
    )
    elemental_number_density.columns = elemental_number_density.columns.astype(
        int
    )
    return elemental_number_density


# initial plasma setup matching ctardis
@pytest.fixture
def iip_plasma(type_iip_workflow):
    """Return the workflow's standard plasma for contract comparisons."""
    return type_iip_workflow.plasma_solver


# "NLTE init" is the first call to update_radiationfield to set up the plasma
@pytest.fixture
def iip_plasma_nlte_init(
    iip_regression_path, iip_plasma, ctardis_compare_config
):
    j_blues_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_j_blues_ctardis_init_nlte.h5",
        key="data",
    )

    # ctardis starts with a constant rad temperature in all cells
    radiation_temp = 9984.96131287 * np.ones(24)
    dilution_factor = np.array(
        [
            0.18635244,
            0.15938095,
            0.11736085,
            0.34665656,
            0.32265696,
            0.30224056,
            0.28436446,
            0.26841929,
            0.2540108,
            0.24086562,
            0.22878441,
            0.21761613,
            0.20724285,
            0.1975702,
            0.18852112,
            0.18003167,
            0.17204798,
            0.16452412,
            0.15742053,
            0.15070279,
            0.14434073,
            0.13830767,
            0.13257993,
            0.12856901,
        ]
    )

    iip_plasma.update_radiationfield(
        radiation_temp,
        dilution_factor,
        j_blues_ctardis,
        ctardis_compare_config.plasma.nlte,
        initialize_nlte=True,
        n_e_convergence_threshold=0.05,
    )
    return iip_plasma


@pytest.fixture
def iip_plasma_after_mc(
    iip_regression_path, iip_plasma_nlte_init, ctardis_compare_config
):
    j_blues_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_j_blues_ctardis_after_mc.h5",
        key="data",
    )

    radiation_temp = np.array(
        [
            9992.27229695,
            9992.59224105,
            9983.78000964,
            9980.58614386,
            9979.83477025,
            9968.05132981,
            9957.88724805,
            9949.36369847,
            9946.8743961,
            9937.71425418,
            9934.85610192,
            9928.23880919,
            9926.40535242,
            9916.93223133,
            9912.22246589,
            9911.051763,
            9910.26097021,
            9901.72775668,
            9895.9432972,
            9891.58754489,
            9886.70685954,
            9880.93185734,
            9876.00858684,
            9872.59842944,
        ]
    )

    dilution_factor = np.array(
        [
            0.3571996,
            0.31756545,
            0.27019532,
            0.36604569,
            0.33787167,
            0.31579601,
            0.29590609,
            0.27936991,
            0.2634541,
            0.24940025,
            0.23579985,
            0.22373621,
            0.21241799,
            0.20254584,
            0.19309261,
            0.18394483,
            0.1755579,
            0.16798016,
            0.16076174,
            0.15381029,
            0.14730572,
            0.14119434,
            0.13532174,
            0.13124624,
        ]
    )

    photo_ion_estimator = pd.read_hdf(
        iip_regression_path / "ctardis_photo_ion_estimator_after_mc.h5",
        key="data",
    )

    stim_recomb_estimator = pd.read_hdf(
        iip_regression_path / "ctardis_stim_recomb_estimator_after_mc.h5",
        key="data",
    )

    bf_heating_estimator = pd.read_hdf(
        iip_regression_path / "ctardis_bf_heating_estimator_after_mc.h5",
        key="data",
    )

    stim_recomb_cooling_estimator = pd.read_hdf(
        iip_regression_path
        / "ctardis_stim_recomb_cooling_estimator_after_mc.h5",
        key="data",
    )

    ff_heating_estimator = [
        4.89135279e-24,
        4.37696370e-24,
        3.75869301e-24,
        4.97847160e-24,
        4.52158002e-24,
        4.21024499e-24,
        3.94991540e-24,
        3.72915649e-24,
        3.58902110e-24,
        3.40170224e-24,
        3.20848519e-24,
        3.03540032e-24,
        2.87314722e-24,
        2.74328938e-24,
        2.61063140e-24,
        2.50640248e-24,
        2.38164559e-24,
        2.26967531e-24,
        2.24509826e-24,
        2.12378192e-24,
        2.02063266e-24,
        1.92509873e-24,
        1.83070678e-24,
        1.77346374e-24,
    ]

    continuum_estimators = {}

    continuum_estimators["photo_ion_estimator"] = photo_ion_estimator
    continuum_estimators["stim_recomb_estimator"] = stim_recomb_estimator
    continuum_estimators["bf_heating_estimator"] = bf_heating_estimator
    continuum_estimators["stim_recomb_cooling_estimator"] = (
        stim_recomb_cooling_estimator
    )
    continuum_estimators["ff_heating_estimator"] = ff_heating_estimator

    iip_plasma_nlte_init.update_radiationfield(
        radiation_temp,
        dilution_factor,
        j_blues_ctardis,
        photoionization_rate_estimator=photo_ion_estimator,
        stimulated_recombination_rate_estimator=stim_recomb_estimator,
        bound_free_heating_estimator=bf_heating_estimator,
        stimulated_recombination_cooling_estimator=stim_recomb_cooling_estimator,
        free_free_heating_estimator=ff_heating_estimator,
    )

    return iip_plasma_nlte_init


def test_type_iip_workflow_initial_plasma_regression(
    type_iip_workflow,
    regression_data,
):
    for attr in INITIAL_PLASMA_SOLVER_REGRESSION_OUTPUTS:
        assert_regression_dataframe(
            regression_data,
            f"workflow_init_{attr}",
            getattr(type_iip_workflow.plasma_solver, attr),
            rtol=1e-10,  # Mac ARM64 tolerance
        )


def test_iip_plasma_initialization(iip_plasma_nlte_init, iip_regression_path):
    tau_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_tau_sobolevs_init_nlte.h5",
        key="data",
    )
    beta_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_beta_sobolevs_init_nlte.h5",
        key="data",
    )
    ion_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_ion_density_init_nlte.h5",
        key="data",
    )
    level_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_level_number_density_init_nlte.h5",
        key="data",
    )
    electron_densities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_electron_densities_init_nlte.h5",
        key="data",
    )
    p_fb_deactivation_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_p_fb_deactivation_init_nlte.h5",
        key="data",
    )
    chi_bf_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_chi_bf_init_nlte.h5",
        key="data",
    )
    t_electrons_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_t_electrons_init_nlte.h5", key="data"
    )

    assert not iip_plasma_nlte_init.transition_probabilities.empty

    print(
        "init ion_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.ion_number_density,
                ion_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_nlte_init.ion_number_density,
        ion_number_density_ctardis,
        rtol=4e-8,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "init tau_sobolevs max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.tau_sobolevs,
                tau_sobolevs_ctardis,
            )
        )
    )
    # Sobolev values are stored differently between codes, so comparing raw data instead
    np.testing.assert_allclose(
        iip_plasma_nlte_init.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=4e-8,
        atol=0,
    )

    print(
        "init beta_sobolev max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.beta_sobolev,
                beta_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=3e-8,
        atol=0,
    )

    print(
        "init level_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.level_number_density,
                level_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_nlte_init.level_number_density,
        level_number_density_ctardis,
        rtol=4e-8,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "init electron_densities max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.electron_densities,
                electron_densities_ctardis,
            )
        )
    )
    pd.testing.assert_series_equal(
        iip_plasma_nlte_init.electron_densities,
        electron_densities_ctardis,
        rtol=2e-12,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "init t_electrons max rel diff: {:.3e}".format(
            _max_rel_diff(
                pd.DataFrame(iip_plasma_nlte_init.t_electrons),
                t_electrons_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.t_electrons,
        t_electrons_ctardis.values,
        rtol=2e-13,
        atol=0,
    )

    print(
        "init p_fb_deactivation max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.p_fb_deactivation,
                p_fb_deactivation_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.p_fb_deactivation.values,
        p_fb_deactivation_ctardis.values,
        rtol=2e-13,
        atol=0,
    )

    print(
        "init chi_bf max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.chi_bf,
                chi_bf_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.chi_bf.values,
        chi_bf_ctardis.values,
        rtol=4e-8,
        atol=0,
    )

    for attr, rtol in INITIAL_CTARDIS_PORT_RTOLS.items():
        _assert_values_allclose(
            getattr(equilibrium_iip_plasma_nlte_init, attr),
            getattr(iip_plasma_nlte_init, attr),
            rtol=rtol,
        )
    _assert_values_allclose(
        equilibrium_iip_plasma_nlte_init.t_rad,
        iip_plasma_nlte_init.t_rad,
        rtol=0.0,
    )
    _assert_values_allclose(
        equilibrium_iip_plasma_nlte_init.w,
        iip_plasma_nlte_init.w,
        rtol=0.0,
    )


# comparison of plasma after the Monte Carlo calculations have been performed
def test_iip_plasma_after_mc(
    iip_regression_path,
    iip_plasma_after_mc,
    regression_data,
):
    tau_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_tau_sobolevs_after_mc.h5",
        key="data",
    )

    beta_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_beta_sobolevs_after_mc.h5",
        key="data",
    )
    ion_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_ion_density_after_mc.h5",
        key="data",
    )
    level_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_level_number_density_after_mc.h5",
        key="data",
    )
    # tolerances are much worse than after init

    assert not iip_plasma_after_mc.transition_probabilities.empty

    print(
        "after MC ion_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.ion_number_density,
                ion_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_after_mc.ion_number_density,
        ion_number_density_ctardis,
        rtol=7e-6,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    # Sobolev values are stored differently between codes, so comparing raw data instead
    print(
        "after MC tau_sobolevs max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.tau_sobolevs,
                tau_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_after_mc.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=7e-6,
        atol=0,
    )

    print(
        "after MC beta_sobolev max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.beta_sobolev,
                beta_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_after_mc.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=7e-6,
        atol=0,
    )

    print(
        "after MC level_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.level_number_density,
                level_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_after_mc.level_number_density,
        level_number_density_ctardis,
        rtol=7e-6,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    assert_regression_dataframe(
        regression_data,
        "after_mc_fractional_heating",
        iip_plasma_after_mc.fractional_heating,
        rtol=2e-10,
    )

    for attr in PLASMA_SOLVER_REGRESSION_OUTPUTS:
        assert_regression_dataframe(
            regression_data,
            f"after_mc_{attr}",
            getattr(iip_plasma_after_mc, attr),
            rtol=4e-8,
        )


def thermal_balance_guess(
    plasma_solver: LegacyPlasmaArray,
) -> tuple[np.ndarray, np.ndarray]:
    """Used for test below to calculate a thermal balance guess from a plasma"""
    max_electron_number_density = (
        plasma_solver.number_density.multiply(
            plasma_solver.number_density.index.values,
            axis=0,
        )
        .sum()
        .values
    )
    electron_fraction = (
        plasma_solver.electron_densities / max_electron_number_density
    ).values

    guess = np.zeros(2 * len(plasma_solver.link_t_rad_t_electron))
    guess[::2] = electron_fraction
    guess[1::2] = plasma_solver.link_t_rad_t_electron

    return guess, max_electron_number_density


def test_thermal_balance_solver(
    iip_regression_path,
    type_iip_workflow,
    iip_plasma_after_mc,
    regression_data,
):

    type_iip_workflow.plasma_solver = deepcopy(iip_plasma_after_mc)
    initial_guess, max_electron_number_density = thermal_balance_guess(
        type_iip_workflow.plasma_solver
    )
    initial_residual = type_iip_workflow.thermal_balance_iteration(
        initial_guess,
        max_electron_number_density,
    )
    assert_regression_dataframe(
        regression_data,
        "thermal_balance_iteration_initial_residual",
        initial_residual,
        atol=3e-14,  # values near zero
    )

    type_iip_workflow.plasma_solver = iip_plasma_after_mc
    type_iip_workflow.solve_thermal_balance()

    tau_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_tau_sobolevs_after_tb.h5",
        key="data",
    )

    for attr in PLASMA_SOLVER_REGRESSION_OUTPUTS:
        _assert_regression_dataframe(
            regression_data,
            f"after_mc_{attr}",
            getattr(iip_plasma_after_mc, attr),
            rtol=4e-8,
        )

    for attr, rtol in POST_MONTE_CARLO_CTARDIS_PORT_RTOLS.items():
        _assert_values_allclose(
            getattr(equilibrium_iip_plasma_after_mc, attr),
            getattr(iip_plasma_after_mc, attr),
            rtol=rtol,
        )
    _assert_values_allclose(
        equilibrium_iip_plasma_after_mc.t_rad,
        iip_plasma_after_mc.t_rad,
        rtol=0.0,
    )
    _assert_values_allclose(
        equilibrium_iip_plasma_after_mc.w,
        iip_plasma_after_mc.w,
        rtol=0.0,
    )

    level_totals = equilibrium_iip_plasma_after_mc.level_number_density.groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    _assert_values_allclose(
        level_totals,
        equilibrium_iip_plasma_after_mc.ion_number_density,
        rtol=1e-12,
        atol=0.0,
    )

    assert_regression_dataframe(
        regression_data,
        "after_thermal_balance_fractional_heating",
        iip_plasma_after_mc.fractional_heating,
        atol=2e-13,  # values near zero
    )

    for attr in PLASMA_SOLVER_REGRESSION_OUTPUTS:
        assert_regression_dataframe(
            regression_data,
            f"after_thermal_balance_{attr}",
            getattr(type_iip_workflow.plasma_solver, attr),
            rtol=3e-11,
        )

    final_guess, max_electron_number_density = thermal_balance_guess(
        type_iip_workflow.plasma_solver
    )
    residual = type_iip_workflow.thermal_balance_iteration(
        final_guess,
        max_electron_number_density,
    )
    assert_regression_dataframe(
        regression_data,
        "thermal_balance_iteration_residual",
        residual,
        atol=1e-13,  # values near zero
    )


def test_thermal_balance_solver(type_iip_workflow, thermal_balance_guess):
    plasma = type_iip_workflow.plasma_solver
    estimator_index = plasma.photo_ion_index
    estimator_columns = plasma.electron_densities.index
    zero_photo_estimator = pd.DataFrame(
        0.0, index=estimator_index, columns=estimator_columns
    )
    zero_line_estimator = pd.DataFrame(
        0.0, index=plasma.lines.index, columns=estimator_columns
    )
    plasma.update(
        photoionization_rate_estimator=zero_photo_estimator,
        stimulated_recombination_rate_estimator=zero_photo_estimator,
        collisional_ionization_rate_coefficient=zero_photo_estimator,
        collisional_excitation_rate_coefficient=zero_line_estimator,
        collisional_deexcitation_rate_coefficient=zero_line_estimator,
        free_free_heating_estimator=pd.Series(0.0, index=estimator_columns),
        bound_free_heating_estimator=zero_photo_estimator,
        stimulated_recombination_cooling_estimator=zero_photo_estimator,
    )
    initial_guess, max_electron_number_density = thermal_balance_guess(plasma)
    initial_residual = type_iip_workflow.thermal_balance_iteration(
        initial_guess,
        max_electron_number_density,
    )
    assert np.all(np.isfinite(initial_residual))

    for attr in (
        "electron_densities",
        "t_electrons",
        "ion_number_density",
        "level_number_density",
        "fractional_heating",
    ):
        value = getattr(plasma, attr)
        values = getattr(value, "values", value)
        assert np.all(np.isfinite(values))


def test_solve_montecarlo(type_iip_workflow, regression_data):
    opacity_states = type_iip_workflow.solve_opacity()
    type_iip_workflow.solve_montecarlo(opacity_states, 1000)

    type_iip_workflow.initialize_spectrum_solver()
    real_packets = type_iip_workflow.spectrum_solver.spectrum_real_packets
    expected_lum_dens = regression_data.sync_ndarray(
        real_packets.luminosity_density_lambda.value
    )

    np.testing.assert_allclose(
        real_packets.luminosity_density_lambda.value,
        expected_lum_dens,
        atol=0,
        rtol=4e-6,
    )
