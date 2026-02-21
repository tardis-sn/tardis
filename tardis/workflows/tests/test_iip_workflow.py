import numpy as np
import pandas as pd
import pytest

from tardis.iip_plasma.continuum.base_continuum_data import ContinuumData
from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow


def _max_rel_diff(actual, expected):
    actual_vals = actual.values
    expected_vals = expected.values

    relative_difference = np.abs(actual_vals - expected_vals) / np.abs(
        expected_vals
    )

    return float(np.nanmax(relative_difference))


@pytest.fixture
def iip_regression_path(tardis_regression_path):
    return tardis_regression_path / "tardis" / "workflows" / "tests"


@pytest.fixture
def iip_atom_data(tardis_regression_path):
    # identical atomic data to that used by C Vogl
    atom_data = AtomData.from_hdf(
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )

    # need to set up macroatom and continuum data for ctardis plasma
    atom_data.prepare_atom_data([1], "macroatom", [(1, 0)], [(1, 0)])

    atom_data.continuum_data = ContinuumData(
        atom_data, selected_continuum_species=[(1, 0)]
    )

    # matching prep in workflow
    atom_data.continuum_data.photoionization_data.loc[(1, 0, 0), "x_sect"] *= (
        0.0
    )

    atom_data.yg_data.columns = list(atom_data.collision_data_temperatures)

    atom_data.nlte_data._init_indices()

    atom_data.has_collision_data = False

    return atom_data


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
def iip_plasma(iip_atom_data, elemental_number_density, ctardis_compare_config):
    plasma = LegacyPlasmaArray(
        elemental_number_density,
        iip_atom_data,
        ctardis_compare_config.supernova.time_explosion.to("s").value,
        nlte_config=ctardis_compare_config.plasma.nlte,
        delta_treatment=None,
        ionization_mode="nlte",
        excitation_mode="dilute-lte",
        line_interaction_type=ctardis_compare_config.plasma.line_interaction_type,
        link_t_rad_t_electron=1.0 * np.ones(24),
        # link_t_rad_t_electron=self.ws**0.25,
        helium_treatment="none",
        heating_rate_data_file=None,
        v_inner=None,
        v_outer=None,
        continuum_treatment=True,
    )

    return plasma


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
        **{},
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
        ctardis_compare_config.plasma.nlte,
        initialize_nlte=False,
        n_e_convergence_threshold=0.05,
        **continuum_estimators,
    )

    return iip_plasma_nlte_init


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
    transition_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_transition_probabilities_init_nlte.h5",
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

    print(
        "init transition_probabilities max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.transition_probabilities,
                transition_probabilities_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_nlte_init.transition_probabilities,
        transition_probabilities_ctardis,
        rtol=3e-8,
        atol=0,
        check_dtype=False,
    )

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


# comparison of plasma after the Monte Carlo calculations have been performed
def test_iip_plasma_after_mc(
    iip_regression_path,
    iip_plasma_after_mc,
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
    transition_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_transition_probabilities_after_mc.h5",
        key="data",
    )

    # tolerances are much worse than after init

    print(
        "after MC transition_probabilities max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.transition_probabilities,
                transition_probabilities_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_after_mc.transition_probabilities,
        transition_probabilities_ctardis,
        rtol=7e-6,
        atol=0,
        check_dtype=False,
    )

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
    )


def test_thermal_balance_solver(
    iip_regression_path, type_iip_workflow, iip_plasma_after_mc
):
    type_iip_workflow.plasma_solver = iip_plasma_after_mc
    type_iip_workflow.solve_thermal_balance()

    tau_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_tau_sobolevs_after_tb.h5",
        key="data",
    )

    beta_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_beta_sobolevs_after_tb.h5",
        key="data",
    )
    ion_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_ion_density_after_tb.h5",
        key="data",
    )
    level_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_level_number_density_after_tb.h5",
        key="data",
    )
    transition_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_transition_probabilities_after_tb.h5",
        key="data",
    )

    print(
        "after thermal balance transition_probabilities max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.transition_probabilities,
                transition_probabilities_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.transition_probabilities,
        transition_probabilities_ctardis,
        rtol=7e-7,
        atol=0,
        check_dtype=False,
    )

    print(
        "after thermal balance ion_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.ion_number_density,
                ion_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.ion_number_density,
        ion_number_density_ctardis,
        rtol=6e-7,
        atol=0,
        check_dtype=False,
    )

    # Sobolev values are stored differently between codes, so comparing raw data instead
    print(
        "after thermal balance tau_sobolevs max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.tau_sobolevs,
                tau_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        type_iip_workflow.plasma_solver.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=7e-7,
        atol=0,
    )

    print(
        "after thermal balance beta_sobolev max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.beta_sobolev,
                beta_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        type_iip_workflow.plasma_solver.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=7e-7,
        atol=0,
    )

    print(
        "after thermal balance level_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.level_number_density,
                level_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.level_number_density,
        level_number_density_ctardis,
        rtol=6e-7,
        atol=0,
        check_dtype=False,
    )


def test_solve_continuum_state_after_nlte_init(
    iip_regression_path, type_iip_workflow, iip_plasma_nlte_init
):
    type_iip_workflow.plasma_solver = iip_plasma_nlte_init
    type_iip_workflow.solve_continuum_state(None)

    recombination_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "continuum_recomb_transition_prob_init_nlte.h5",
        key="data",
    )

    pd.testing.assert_frame_equal(
        type_iip_workflow.base_continuum.recombination_transition_probabilities.dataframe,
        recombination_probabilities_ctardis,
        rtol=4e-3,
    )


def test_solve_montecarlo(type_iip_workflow):
    opacity_states = type_iip_workflow.solve_opacity()
    type_iip_workflow.solve_montecarlo(opacity_states, 1000)
