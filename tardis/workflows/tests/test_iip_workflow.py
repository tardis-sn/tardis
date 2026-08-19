from copy import copy, deepcopy
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from astropy import units as u
from scipy.optimize import root

from tardis import constants as const
from tardis.conftest import assert_regression_dataframe
from tardis.iip_plasma.properties.ion_population import NLTEIonNumberDensity
from tardis.iip_plasma.properties.partition_function import (
    PartitionFunction as IIPPartitionFunction,
)
from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.configuration.config_reader import Configuration
from tardis.opacities.continuum.macro_atom_state import (
    ContinuumMacroAtomState,
)
from tardis.opacities.tau_sobolev import SOBOLEV_COEFFICIENT
from tardis.plasma import BasePlasma
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluator,
    calculate_lte_populations,
)
from tardis.plasma.equilibrium.inputs import (
    NumberDensityPerShell,
    SobolevInputs,
)
from tardis.plasma.equilibrium.ion_populations import (
    IonPopulationSolver,
)
from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import (
    EstimatedIonRateMatrix,
    RateMatrix,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
    CollisionalIonizationRateSolver,
    CollisionalIonizationSeaton,
    EstimatedPhotoionizationRateSolver,
    SpontaneousRecombinationCoeffSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.equilibrium.thermal_balance import ThermalBalanceSolver
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.transport.montecarlo.estimators import init_estimators_continuum
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow


@dataclass(frozen=True)
class ContinuumComparisonState:
    plasma: LegacyPlasmaArray
    continuum: ContinuumMacroAtomState
    photoionization_data: pd.DataFrame
    photoionization_index: pd.MultiIndex
    upper_ion_index: pd.MultiIndex
    radiation_field: DilutePlanckianRadiationField
    electron_temperature: u.Quantity
    electron_distribution: ThermalElectronEnergyDistribution
    level_to_ion_population_factor: pd.DataFrame
    lte_ion_population: pd.DataFrame
    lte_level_population: pd.DataFrame


@dataclass(frozen=True)
class CollisionalBoundRates:
    excitation: pd.DataFrame
    deexcitation: pd.DataFrame
    excitation_index: pd.MultiIndex
    deexcitation_index: pd.MultiIndex


def _max_rel_diff(actual, expected):
    """Helper function to print relative diffs for checking broken tests"""
    actual_vals = actual.values
    expected_vals = expected.values

    relative_difference = np.abs(actual_vals - expected_vals) / np.abs(
        expected_vals
    )

    return float(np.nanmax(relative_difference))


PLASMA_SOLVER_REGRESSION_OUTPUTS = (
    "electron_densities",
    "t_electrons",
    "link_t_rad_t_electron",
    "p_fb_deactivation",
    "chi_bf",
    "stimulated_emission_factor",
    "b",
    "j_blues",
)

STANDARD_PLASMA_SOLVER_REGRESSION_OUTPUTS = tuple(
    output for output in PLASMA_SOLVER_REGRESSION_OUTPUTS if output != "b"
)


INITIAL_PLASMA_SOLVER_REGRESSION_OUTPUTS = (
    "ion_number_density",
    "tau_sobolevs",
    "beta_sobolev",
    "level_number_density",
    *STANDARD_PLASMA_SOLVER_REGRESSION_OUTPUTS,
)


@pytest.fixture
def iip_regression_path(tardis_regression_path):
    return tardis_regression_path / "tardis" / "workflows" / "tests"


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

    config.plasma.nlte.species = [
        (1, 0)
    ]  # Configure the hydrogen NLTE state used by the ctardis comparison.
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


@pytest.fixture
def iip_plasma_after_thermal_balance(
    iip_regression_path: Path,
    iip_plasma_after_mc: LegacyPlasmaArray,
) -> LegacyPlasmaArray:
    """Rebuild the stored accepted legacy state for evaluator comparisons."""
    regression_file = (
        iip_regression_path
        / "test_iip_workflow"
        / "test_thermal_balance_solver.h5"
    )
    electron_densities = pd.read_hdf(
        regression_file, key="after_thermal_balance_electron_densities"
    )["value"].to_numpy()
    link_t_rad_t_electron = pd.read_hdf(
        regression_file, key="after_thermal_balance_link_t_rad_t_electron"
    )["value"].to_numpy()
    electron_temperatures = pd.read_hdf(
        regression_file, key="after_thermal_balance_t_electrons"
    )["value"].to_numpy()
    plasma = deepcopy(iip_plasma_after_mc)
    plasma.update(
        previous_ion_number_density=plasma.ion_number_density.copy(),
        previous_electron_densities=electron_densities,
        previous_beta_sobolev=plasma.beta_sobolev.copy(),
        link_t_rad_t_electron=link_t_rad_t_electron,
        previous_b=plasma.b,
        previous_t_electrons=electron_temperatures,
    )
    return plasma


@pytest.fixture
def iip_charge_conserving_rate_matrix(
    iip_plasma_after_mc: LegacyPlasmaArray,
) -> SimpleNamespace:
    """Adapt full-dataset IIP rate tables to the new solver interface."""
    plasma = iip_plasma_after_mc
    legacy_solver = NLTEIonNumberDensity(plasma)
    partition_function = IIPPartitionFunction.calculate(
        plasma.level_boltzmann_factor
    )
    level_population_fraction = legacy_solver._calculate_level_pop_fractions(
        plasma.level_boltzmann_factor,
        partition_function,
    )
    ion_index = partition_function.index
    radiative_recombination = legacy_solver._calculate_alpha_tot(
        plasma.alpha_stim,
        plasma.alpha_sp,
        ion_index,
    ).loc[(1, 0)]
    radiative_ionization = legacy_solver._calculate_tot_ion_rate(
        plasma.gamma,
        level_population_fraction,
        ion_index,
    ).loc[(1, 0)]
    collisional_ionization = legacy_solver._calculate_tot_ion_rate(
        plasma.coll_ion_coeff,
        level_population_fraction,
        ion_index,
    ).loc[(1, 0)]
    collisional_recombination = legacy_solver._calculate_coll_recomb_tot(
        plasma.coll_recomb_coeff,
        ion_index,
    ).loc[(1, 0)]
    radiative_ionization_values = radiative_ionization.to_numpy()
    radiative_recombination_values = radiative_recombination.to_numpy()
    collisional_ionization_values = collisional_ionization.to_numpy()
    collisional_recombination_values = collisional_recombination.to_numpy()

    def solve(
        radiation_field: object,
        electron_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        level_to_continuum_saha_factor: pd.DataFrame,
    ) -> pd.DataFrame:
        """Build the IIP two-stage hydrogen matrices at trial densities."""
        electron_number_density = electron_distribution.number_density.value
        ionization_rate = (
            radiative_ionization_values
            + collisional_ionization_values * electron_number_density
        )
        recombination_rate = (
            radiative_recombination_values * electron_number_density
            + collisional_recombination_values * electron_number_density**2
        )
        rate_matrices = np.empty((len(electron_number_density), 2, 2))
        rate_matrices[:, 0, 0] = -ionization_rate
        rate_matrices[:, 0, 1] = recombination_rate
        rate_matrices[:, 1, :] = 1.0
        rate_matrix_array = np.empty((1, len(rate_matrices)), dtype=object)
        rate_matrix_array[0] = list(rate_matrices)
        return pd.DataFrame(
            rate_matrix_array,
            index=pd.Index([1], name="atomic_number"),
            columns=radiative_recombination.index,
        )

    return SimpleNamespace(
        ion_population_index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
        solve=solve,
    )


def test_charge_conserving_solver_matches_iip_with_full_atomic_data(
    iip_plasma_after_mc: LegacyPlasmaArray,
    iip_charge_conserving_rate_matrix: SimpleNamespace,
) -> None:
    """Match IIP NLTE populations with identical post-Monte-Carlo rates."""
    plasma = iip_plasma_after_mc
    legacy_parent = SimpleNamespace(
        nlte_species=plasma.nlte_species,
        previous_ion_number_density=None,
        previous_electron_densities=None,
    )
    expected_ion_population, expected_electron_density = NLTEIonNumberDensity(
        legacy_parent
    ).calculate(
        plasma.phi,
        plasma.alpha_stim,
        plasma.alpha_sp,
        plasma.gamma,
        plasma.coll_ion_coeff,
        plasma.coll_recomb_coeff,
        plasma.number_density,
        plasma.level_boltzmann_factor,
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        plasma.t_electrons * u.K,
        plasma.electron_densities.to_numpy() / u.cm**3,
    )
    actual_ion_population, actual_electron_density = IonPopulationSolver(
        iip_charge_conserving_rate_matrix
    ).solve(
        None,
        electron_distribution,
        plasma.number_density,
        plasma.lte_level_number_density,
        plasma.level_number_density,
        plasma.lte_ion_number_density,
        plasma.ion_number_density,
        plasma.partition_function,
        plasma.level_boltzmann_factor,
        level_to_continuum_saha_factor=plasma.phi_lucy,
    )

    pd.testing.assert_frame_equal(
        actual_ion_population,
        expected_ion_population,
        check_dtype=False,
        check_names=False,
        rtol=3e-10,
        atol=0.0,
    )
    pd.testing.assert_series_equal(
        actual_electron_density,
        expected_electron_density,
        check_dtype=False,
        check_names=False,
        rtol=3e-10,
        atol=0.0,
    )


def test_charge_conserving_solver_only_resolves_unconverged_shells(
    iip_plasma_after_mc: LegacyPlasmaArray,
    iip_charge_conserving_rate_matrix: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep converged shell charge roots fixed during lagged iterations."""
    plasma = iip_plasma_after_mc
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        plasma.t_electrons * u.K,
        plasma.electron_densities.to_numpy() / u.cm**3,
    )
    estimated_ion_population = plasma.ion_number_density.copy()
    estimated_ion_population.iloc[:, 0] = 0.0
    solver = IonPopulationSolver(iip_charge_conserving_rate_matrix)
    solve_shell_charge = solver.solve_shell_charge
    solved_shell_indices = []

    def record_solve_shell_charge(
        shell_idx: int, *args: object, **kwargs: object
    ) -> float:
        solved_shell_indices.append(shell_idx)
        return solve_shell_charge(shell_idx, *args, **kwargs)

    monkeypatch.setattr(solver, "solve_shell_charge", record_solve_shell_charge)
    solver.solve(
        None,
        electron_distribution,
        plasma.number_density,
        plasma.lte_level_number_density,
        plasma.level_number_density,
        plasma.lte_ion_number_density,
        estimated_ion_population,
        plasma.partition_function,
        plasma.level_boltzmann_factor,
        # This scheduling test starts from the legacy fixed point. The adjacent
        # parity test bounds its difference from the standard owner at 3e-10,
        # so 1e-9 treats those unperturbed shells as converged.
        tolerance=1e-9,
        level_to_continuum_saha_factor=plasma.phi_lucy,
    )

    assert [
        solved_shell_indices.count(shell_idx)
        for shell_idx in range(len(plasma.number_density.columns))
    ] == [2] + [1] * (len(plasma.number_density.columns) - 1)


def test_type_iip_workflow_initial_plasma_regression(
    type_iip_workflow,
    regression_data,
):
    """Compare initial IIP plasma outputs with regression references."""
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
    for attr in INITIAL_PLASMA_SOLVER_REGRESSION_OUTPUTS:
        assert_regression_dataframe(
            regression_data,
            f"workflow_init_{attr}",
            outputs[attr],
            # The standard dilute-LTE bootstrap and legacy IIP initialization
            # use different nonlinear owners before the first MC estimator
            # snapshot. Preserve observable parity without requiring identical
            # intermediate iterates.
            rtol=1e-4,
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
        check_names=False,
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


@pytest.fixture(scope="module")
def continuum_comparison_state(
    tardis_regression_path: Path,
) -> ContinuumComparisonState:
    comparison_config = Configuration.from_yaml(
        Path(__file__).parent / "data" / "ctardis_compare.yml"
    )
    comparison_config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    comparison_config.plasma.nlte.species = [(1, 0)]
    workflow = TypeIIPWorkflow(comparison_config)
    plasma = workflow.plasma_solver
    continuum = workflow.continuum_macro_atom_state
    photoionization_data = plasma.photo_ion_cross_sections
    photoionization_index = photoionization_data.index.unique()
    upper_ion_index = pd.MultiIndex.from_arrays(
        [
            photoionization_index.get_level_values("atomic_number"),
            photoionization_index.get_level_values("ion_number") + 1,
        ],
        names=["atomic_number", "ion_number"],
    ).unique()
    radiation_field = DilutePlanckianRadiationField(
        np.asarray(plasma.t_rad) * u.K,
        np.asarray(plasma.w),
    )
    electron_temperature = np.asarray(plasma.t_electrons) * u.K
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        electron_temperature,
        plasma.electron_densities.to_numpy() / u.cm**3,
    )
    maximum_electron_density = (
        plasma.number_density.multiply(
            plasma.number_density.index.to_numpy(), axis=0
        )
        .sum()
        .to_numpy()
    )
    evaluator = workflow._build_thermal_balance_evaluator(
        maximum_electron_density, analytic=True
    )
    _, level_to_ion_population_factor, partition_function, boltzmann_factor = (
        evaluator.calculate_continuum_coefficients(plasma.t_electrons)
    )
    level_to_ion_population_factor = level_to_ion_population_factor.loc[
        photoionization_index
    ]
    lte_ion_population, lte_level_population = calculate_lte_populations(
        plasma.thermal_phi_lte,
        partition_function,
        plasma.number_density,
        plasma.electron_densities,
        boltzmann_factor,
        plasma.atomic_data.levels.loc[plasma.level_number_density.index],
    )

    return ContinuumComparisonState(
        plasma=plasma,
        continuum=continuum,
        photoionization_data=photoionization_data,
        photoionization_index=photoionization_index,
        upper_ion_index=upper_ion_index,
        radiation_field=radiation_field,
        electron_temperature=electron_temperature,
        electron_distribution=electron_distribution,
        level_to_ion_population_factor=level_to_ion_population_factor,
        lte_ion_population=lte_ion_population,
        lte_level_population=lte_level_population,
    )


@pytest.fixture(scope="module")
def collisional_bound_rates(
    continuum_comparison_state: ContinuumComparisonState,
) -> CollisionalBoundRates:
    collisional_rates = ThermalCollisionalRateSolver(
        continuum_comparison_state.plasma.atomic_data.levels,
        continuum_comparison_state.plasma.atomic_data.lines,
        continuum_comparison_state.plasma.atomic_data.collision_data_temperatures,
        continuum_comparison_state.plasma.atomic_data.yg_data,
        collision_strengths_type="cmfgen",
    ).solve(continuum_comparison_state.electron_temperature)
    source_levels = collisional_rates.index.get_level_values(
        "level_number_source"
    )
    destination_levels = collisional_rates.index.get_level_values(
        "level_number_destination"
    )
    excitation = collisional_rates[source_levels < destination_levels]
    deexcitation = collisional_rates[source_levels > destination_levels]
    excitation_index = excitation.index.droplevel(
        ["ion_number_source", "ion_number_destination"]
    ).rename(
        {
            "level_number_source": "level_number_lower",
            "level_number_destination": "level_number_upper",
        }
    )
    deexcitation_index = (
        deexcitation.index.swaplevel(
            "level_number_source", "level_number_destination"
        )
        .droplevel(["ion_number_source", "ion_number_destination"])
        .rename(
            {
                "level_number_destination": "level_number_lower",
                "level_number_source": "level_number_upper",
            }
        )
    )
    return CollisionalBoundRates(
        excitation,
        deexcitation,
        excitation_index,
        deexcitation_index,
    )


@pytest.fixture(scope="module")
def collisional_ionization_rate(
    continuum_comparison_state: ContinuumComparisonState,
) -> pd.DataFrame:
    return CollisionalIonizationSeaton(
        continuum_comparison_state.photoionization_data
    ).solve(continuum_comparison_state.electron_temperature)


@pytest.fixture(scope="module")
def equilibrium_cooling_channels(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
    collisional_ionization_rate: pd.DataFrame,
) -> npt.NDArray[np.float64]:
    bound_free_cooling = BoundFreeThermalRates(
        continuum_comparison_state.photoionization_data
    ).solve(
        continuum_comparison_state.plasma.level_number_density,
        continuum_comparison_state.plasma.ion_number_density,
        continuum_comparison_state.electron_distribution,
        continuum_comparison_state.level_to_ion_population_factor,
        continuum_comparison_state.radiation_field,
    )[1]
    free_free_cooling = FreeFreeThermalRates().solve(
        pd.Series(
            0.0,
            index=continuum_comparison_state.plasma.electron_densities.index,
        ),
        continuum_comparison_state.electron_distribution,
        continuum_comparison_state.plasma.ion_number_density,
    )[1]
    collisional_ionization_cooling = CollisionalIonizationThermalRates(
        continuum_comparison_state.photoionization_data
    ).solve(
        continuum_comparison_state.electron_distribution.number_density,
        continuum_comparison_state.plasma.ion_number_density,
        continuum_comparison_state.plasma.level_number_density,
        collisional_ionization_rate,
        continuum_comparison_state.level_to_ion_population_factor,
    )[1]
    collisional_bound_cooling = CollisionalBoundThermalRates(
        continuum_comparison_state.plasma.atomic_data.lines.loc[
            collisional_bound_rates.excitation_index
        ]
    ).solve(
        continuum_comparison_state.electron_distribution.number_density,
        collisional_bound_rates.deexcitation.set_axis(
            collisional_bound_rates.deexcitation_index
        ),
        collisional_bound_rates.excitation.set_axis(
            collisional_bound_rates.excitation_index
        ),
        continuum_comparison_state.plasma.level_number_density,
    )[1]

    return np.vstack(
        [
            collisional_bound_cooling.to_numpy(),
            collisional_ionization_cooling.to_numpy(),
            bound_free_cooling.to_numpy(),
            free_free_cooling.to_numpy(),
        ]
    )


def test_radiative_ionization_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
) -> None:

    radiative_ionization_rate = AnalyticCorrectedPhotoionizationCoeffSolver(
        continuum_comparison_state.photoionization_data
    ).solve(
        continuum_comparison_state.radiation_field,
        continuum_comparison_state.electron_temperature,
        continuum_comparison_state.lte_level_population.loc[
            continuum_comparison_state.photoionization_index
        ],
        continuum_comparison_state.plasma.level_number_density.loc[
            continuum_comparison_state.photoionization_index
        ],
        continuum_comparison_state.lte_ion_population.loc[
            continuum_comparison_state.upper_ion_index
        ],
        continuum_comparison_state.plasma.ion_number_density.loc[
            continuum_comparison_state.upper_ion_index
        ],
    )
    pd.testing.assert_index_equal(
        radiative_ionization_rate.index,
        continuum_comparison_state.continuum.radiative_ionization_rate.index,
    )
    np.testing.assert_allclose(
        radiative_ionization_rate.to_numpy(),
        continuum_comparison_state.continuum.radiative_ionization_rate.to_numpy(),
        rtol=2e-6,
        atol=0.0,
    )


def test_radiative_recombination_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
) -> None:
    radiative_recombination_rate = (
        SpontaneousRecombinationCoeffSolver(
            continuum_comparison_state.photoionization_data
        ).solve(continuum_comparison_state.electron_temperature)
        * continuum_comparison_state.level_to_ion_population_factor
    )
    pd.testing.assert_index_equal(
        radiative_recombination_rate.index,
        continuum_comparison_state.continuum.radiative_recombination_rate.index,
    )
    np.testing.assert_allclose(
        radiative_recombination_rate.to_numpy(),
        continuum_comparison_state.continuum.radiative_recombination_rate.to_numpy(),
        rtol=2e-4,
        atol=0.0,
    )


def test_collisional_excitation_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
) -> None:
    np.testing.assert_allclose(
        collisional_bound_rates.excitation.to_numpy(),
        continuum_comparison_state.continuum.collisional_excitation_rate.loc[
            collisional_bound_rates.excitation_index
        ].to_numpy(),
        # The standard and IIP implementations use different cgs constant
        # sources for the complete transition table.
        rtol=2e-5,
        atol=0.0,
    )


def test_collisional_deexcitation_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
) -> None:
    np.testing.assert_allclose(
        collisional_bound_rates.deexcitation.to_numpy(),
        continuum_comparison_state.continuum.collisional_deexcitation_rate.loc[
            collisional_bound_rates.deexcitation_index
        ].to_numpy(),
        rtol=2e-5,
        atol=0.0,
    )


def test_collisional_ionization_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_ionization_rate: pd.DataFrame,
) -> None:
    pd.testing.assert_frame_equal(
        collisional_ionization_rate,
        continuum_comparison_state.continuum.collisional_ionization_rate.loc[
            collisional_ionization_rate.index
        ],
        check_names=False,
        check_column_type=False,
    )


def test_collisional_recombination_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_ionization_rate: pd.DataFrame,
) -> None:
    collisional_recombination_rate = (
        collisional_ionization_rate
        * continuum_comparison_state.level_to_ion_population_factor
    )
    np.testing.assert_allclose(
        collisional_recombination_rate.to_numpy(),
        continuum_comparison_state.continuum.collisional_recombination_rate.loc[
            collisional_recombination_rate.index
        ].to_numpy(),
        rtol=2e-5,
        atol=0.0,
    )


def test_cooling_channel_probabilities_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    equilibrium_cooling_channels: npt.NDArray[np.float64],
) -> None:
    actual = equilibrium_cooling_channels / equilibrium_cooling_channels.sum(
        axis=0
    )
    expected = np.vstack(
        [
            continuum_comparison_state.continuum.collisional_excitation_cooling_probability,
            continuum_comparison_state.continuum.collisional_ionization_cooling_probability,
            continuum_comparison_state.continuum.radiative_recombination_cooling_probability,
            continuum_comparison_state.continuum.free_free_cooling_probability,
        ]
    )
    np.testing.assert_allclose(
        actual,
        expected,
        rtol=3e-4,
        atol=0.0,
    )


@pytest.mark.parametrize(
    "process_name",
    [
        "collisional_excitation",
        "collisional_ionization",
        "radiative_recombination",
    ],
)
def test_iip_process_probabilities_normalize_per_shell(
    continuum_comparison_state: ContinuumComparisonState,
    process_name: str,
) -> None:
    probabilities = getattr(
        continuum_comparison_state.continuum,
        f"{process_name}_cooling_array",
    )
    assert probabilities.shape[0] == len(
        continuum_comparison_state.plasma.t_electrons
    )
    np.testing.assert_allclose(
        probabilities.sum(axis=1),
        1.0,
        rtol=1e-12,
        atol=0.0,
    )


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "The standalone bound-bound solver does not include IIP's coupled "
        "continuum rates, ion-stage ratios, or population-dependent Sobolev "
        "escape probabilities."
    ),
)
def test_standalone_level_populations_equality_with_iip_plasma(
    iip_plasma_nlte_init: LegacyPlasmaArray,
    iip_atom_data: object,
) -> None:
    """Characterize the known standalone-equilibrium/IIP population gap."""
    atom_data = iip_atom_data
    line_index = atom_data.lines.index
    hydrogen_lines = atom_data.lines.loc[
        (line_index.get_level_values("atomic_number") == 1)
        & (line_index.get_level_values("ion_number") == 0)
    ]
    standard_rate_matrix = RateMatrix(
        RadiativeRatesSolver(hydrogen_lines),
        ThermalCollisionalRateSolver(
            atom_data.levels,
            hydrogen_lines,
            atom_data.collision_data_temperatures,
            atom_data.yg_data.loc[(1, 0, slice(None), slice(None)), :],
            "cmfgen",
        ),
        atom_data.levels,
    )
    radiation_field = DilutePlanckianRadiationField(
        np.asarray(iip_plasma_nlte_init.t_rad) * u.K,
        np.asarray(iip_plasma_nlte_init.w),
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        np.asarray(iip_plasma_nlte_init.t_electrons) * u.K,
        np.asarray(iip_plasma_nlte_init.electron_densities) / u.cm**3,
    )
    standard_matrices = standard_rate_matrix.solve(
        radiation_field,
        electron_distribution,
    )
    standard_populations = (
        LevelPopulationSolver(standard_matrices, atom_data.levels)
        .solve()
        .loc[(1, 0)]
    )
    iip_populations = iip_plasma_nlte_init.level_number_density.loc[(1, 0)]
    iip_ion_population = iip_plasma_nlte_init.ion_number_density.loc[(1, 0)]
    iip_populations = iip_populations.divide(iip_ion_population, axis=1)

    common_levels = standard_populations.index.intersection(
        iip_populations.index
    )
    np.testing.assert_allclose(
        standard_populations.loc[common_levels].to_numpy(),
        iip_populations.loc[common_levels].to_numpy(),
        rtol=1e-6,
        atol=0,
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
        check_names=False,
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
    plasma_solver: BasePlasma,
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


def test_thermal_balance_iteration_delegates_to_evaluator() -> None:
    """Map the outer candidate to one pure evaluator residual."""
    workflow = TypeIIPWorkflow.__new__(TypeIIPWorkflow)
    workflow.plasma_solver = SimpleNamespace(t_rad=np.array([1.0e4, 2.0e4]))
    level_seed = pd.DataFrame([[0.6, 0.7], [0.4, 0.3]])
    workflow._thermal_balance_level_seed = level_seed

    class RecordingEvaluator:
        def __init__(self) -> None:
            self.call_count = 0
            self.arguments: tuple[object, ...] | None = None
            self.result = SimpleNamespace(
                electron_residual=pd.Series([0.1, -0.2]),
                fractional_heating=pd.Series([0.3, -0.4]),
            )

        def evaluate(
            self,
            trial_electron_density: npt.ArrayLike,
            electron_temperature: npt.ArrayLike,
            candidate_level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            self.call_count += 1
            self.arguments = (
                np.asarray(trial_electron_density),
                np.asarray(electron_temperature),
                candidate_level_seed,
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
    trial_density, temperature, actual_seed = evaluator.arguments
    np.testing.assert_allclose(trial_density, [1.0e9, 3.0e9])
    np.testing.assert_allclose(temperature, [8.0e3, 2.2e4])
    pd.testing.assert_frame_equal(actual_seed, level_seed)
    np.testing.assert_allclose(residual, [0.1, 0.3, -0.2, -0.4])
    assert not hasattr(workflow, "_thermal_balance_evaluation")


def test_nlte_beta_sobolev_calculation_matches_plasma_property(
    iip_plasma_after_mc,
):
    """Compare optimized NLTE beta Sobolev values with the plasma property."""
    nlte_property = iip_plasma_after_mc.plasma_properties_dict[
        "LevelBoltzmannFactorNLTE"
    ]
    beta_sobolev = nlte_property._calculate_beta_sobolevs(
        iip_plasma_after_mc.level_number_density[0].to_numpy()
    )

    np.testing.assert_allclose(
        beta_sobolev,
        iip_plasma_after_mc.beta_sobolev.values[:, [0]],
        rtol=5e-13,  # AVX-512 tolerance
        atol=0.0,
    )


def test_iip_augmented_and_reduced_level_equations_match(
    iip_plasma_after_mc: LegacyPlasmaArray,
    iip_plasma_after_thermal_balance: LegacyPlasmaArray,
) -> None:
    """Freeze the fixed-density augmented/reduced H I equation map."""
    post_mc_candidate, maximum_electron_density = thermal_balance_guess(
        iip_plasma_after_mc
    )
    accepted_candidate, _ = thermal_balance_guess(
        iip_plasma_after_thermal_balance
    )

    shell_indices = np.array([0, 2, 3, 8, 23])
    number_of_shells = len(iip_plasma_after_mc.t_rad)
    minimum_link = 1500.0 / np.min(iip_plasma_after_mc.t_rad)
    lower_bound = np.tile([0.0, minimum_link], number_of_shells)
    upper_bound = np.tile([1.0, 1.5], number_of_shells)
    midpoint_candidate = (lower_bound + upper_bound) / 2.0
    candidates = {
        "accepted": accepted_candidate,
        "post_mc": post_mc_candidate,
        "midpoint": midpoint_candidate,
        "near_lower": lower_bound + 0.01 * (upper_bound - lower_bound),
        "near_upper": upper_bound - 0.01 * (upper_bound - lower_bound),
    }

    plasma = iip_plasma_after_mc
    species = (1, 0)
    nlte_excitation = plasma.plasma_properties_dict["LevelBoltzmannFactorNLTE"]
    original_calculate_beta = nlte_excitation._calculate_beta_sobolevs
    j_blues = pd.DataFrame(
        plasma.j_blues,
        index=plasma.lines.index,
        columns=range(number_of_shells),
    )

    for candidate in candidates.values():
        candidate_electron_density = candidate[::2] * maximum_electron_density
        candidate_electron_temperature = (
            np.asarray(plasma.t_rad) * candidate[1::2]
        )
        for shell_idx in shell_indices:
            electron_density = np.array([candidate_electron_density[shell_idx]])
            electron_temperature = np.array(
                [candidate_electron_temperature[shell_idx]]
            )
            assert electron_density[0] > 0.0

            (
                _rates_matrix,
                _normalization_vector,
                lines_idx,
                radiative_upward_rates,
                radiative_downward_rates,
                radiative_upward_index,
                radiative_downward_index,
            ) = nlte_excitation._setup_bb_rates(
                species,
                plasma.atomic_data,
                plasma.nlte_data,
                electron_temperature,
                j_blues.iloc[:, [shell_idx]],
                np.ones((len(plasma.lines), 1)),
                pd.Series(electron_density),
            )

            gamma = plasma.gamma.loc[species].iloc[:, shell_idx].to_numpy()
            alpha = (
                plasma.alpha_stim.loc[species].iloc[:, shell_idx].to_numpy()
                + plasma.alpha_sp.loc[species].iloc[:, shell_idx].to_numpy()
            )
            collisional_ionization = (
                plasma.coll_ion_coeff.loc[species].iloc[:, shell_idx].to_numpy()
            )
            collisional_recombination = (
                plasma.coll_recomb_coeff.loc[species]
                .iloc[:, shell_idx]
                .to_numpy()
            )
            collision_rates = nlte_excitation._setup_collision_rate_matrices(
                plasma.coll_exc_coeff.loc[species].iloc[:, [shell_idx]],
                plasma.coll_deexc_coeff.loc[species].iloc[:, [shell_idx]],
                plasma.atomic_data.levels.energy.loc[species].count(),
                electron_density,
            )[:, :, 0]
            gamma_vector = gamma + collisional_ionization * electron_density[0]
            recombination_vector = (
                alpha + collisional_recombination * electron_density[0]
            ) * electron_density[0]
            recombination_vector[0] = 0.0
            remaining_rates = collision_rates - np.diag(gamma_vector)
            phis = (
                plasma.phi.loc[(species[0],), shell_idx].to_numpy()
                / electron_density[0]
            )
            base_level_density = plasma.level_number_density.iloc[
                :, shell_idx
            ].to_numpy(copy=True)
            species_level_positions = (
                nlte_excitation._get_species_level_positions(species)
            )
            number_density = plasma.number_density.loc[species[0], shell_idx]
            arguments = (
                nlte_excitation,
                recombination_vector,
                remaining_rates,
                lines_idx,
                radiative_upward_rates[:, :, 0],
                radiative_downward_rates[:, :, 0],
                radiative_upward_index,
                radiative_downward_index,
                number_density,
                gamma_vector,
                species,
                phis,
                base_level_density,
                species_level_positions,
            )

            def calculate_q(level_fractions: npt.NDArray[np.float64]) -> float:
                return float(
                    np.dot(gamma_vector, level_fractions)
                    / recombination_vector.sum()
                )

            def calculate_augmented_residual(
                values: npt.NDArray[np.float64],
            ) -> npt.NDArray[np.float64]:
                return nlte_excitation._rate_equations(
                    values.copy(),
                    arguments[0],
                    arguments[1].copy(),
                    arguments[2],
                    arguments[3],
                    arguments[4].copy(),
                    arguments[5].copy(),
                    arguments[6],
                    arguments[7],
                    arguments[8],
                    arguments[9],
                    arguments[10],
                    arguments[11],
                    arguments[12],
                    arguments[13],
                )

            def calculate_reduced_residual(
                level_fractions: npt.NDArray[np.float64],
            ) -> npt.NDArray[np.float64]:
                return calculate_augmented_residual(
                    np.hstack([level_fractions, calculate_q(level_fractions)])
                )[:-1]

            initial_fractions = (
                plasma.level_number_density.iloc[:, shell_idx]
                .loc[species]
                .to_numpy(copy=True)
            )
            initial_fractions /= initial_fractions.sum()
            trial_fractions = np.full_like(
                initial_fractions, 1.0 / len(initial_fractions)
            )

            beta_inputs = []
            calculate_beta = original_calculate_beta

            def record_beta(
                level_density_values: npt.NDArray[np.float64],
                line_indices: npt.NDArray[np.int64],
            ) -> npt.NDArray[np.float64]:
                beta_inputs.append(level_density_values.copy())
                return calculate_beta(level_density_values, line_indices)

            nlte_excitation._calculate_beta_sobolevs = record_beta
            augmented_residual = calculate_augmented_residual(
                np.hstack([trial_fractions, calculate_q(trial_fractions)])
            )
            perturbed_fractions = trial_fractions.copy()
            perturbed_fractions[1] *= 1.1
            perturbed_fractions /= perturbed_fractions.sum()
            calculate_augmented_residual(
                np.hstack(
                    [perturbed_fractions, calculate_q(perturbed_fractions)]
                )
            )
            assert not np.array_equal(beta_inputs[-1], beta_inputs[-2])
            reduced_residual = calculate_reduced_residual(trial_fractions)
            np.testing.assert_allclose(
                augmented_residual[:-1],
                reduced_residual,
                rtol=1e-12,
                atol=1e-12,
            )
            assert (
                abs(augmented_residual[-1])
                / (
                    abs(np.dot(gamma_vector, trial_fractions))
                    + abs(
                        recombination_vector.sum()
                        * calculate_q(trial_fractions)
                    )
                )
                < 1e-12
            )
            augmented_solution = root(
                calculate_augmented_residual,
                np.hstack([initial_fractions, calculate_q(initial_fractions)]),
                options={"xtol": 1e-12},
            )
            reduced_solution = root(
                calculate_reduced_residual,
                initial_fractions,
                options={"xtol": 1e-12},
            )
            residual_scale = max(
                1.0,
                np.max(np.abs(gamma_vector)),
                np.max(np.abs(recombination_vector)),
            )
            assert (
                np.max(np.abs(augmented_solution.fun)) / residual_scale < 1e-10
            )
            assert np.max(np.abs(reduced_solution.fun)) / residual_scale < 1e-10
            solved_fractions = augmented_solution.x[:-1]
            reduced_fractions = reduced_solution.x
            assert np.isfinite(solved_fractions).all()
            assert np.isfinite(reduced_fractions).all()
            assert np.all(solved_fractions >= 0.0)
            assert np.all(reduced_fractions >= 0.0)
            np.testing.assert_allclose(solved_fractions.sum(), 1.0)
            np.testing.assert_allclose(reduced_fractions.sum(), 1.0)
            np.testing.assert_allclose(
                reduced_fractions,
                solved_fractions,
                rtol=3e-10,
                atol=0.0,
            )
            np.testing.assert_allclose(
                calculate_q(solved_fractions),
                augmented_solution.x[-1],
                rtol=1e-10,
                atol=0.0,
            )
            nlte_excitation._calculate_beta_sobolevs = original_calculate_beta

            assert np.isfinite(initial_fractions).all()
            assert np.all(initial_fractions >= 0.0)
            np.testing.assert_allclose(initial_fractions.sum(), 1.0)


@pytest.fixture
def iip_equilibrium_evaluator(
    iip_plasma_after_thermal_balance: LegacyPlasmaArray,
) -> PlasmaEquilibriumEvaluator:
    """Build the evaluator used by IIP workflow parity tests."""
    plasma = iip_plasma_after_thermal_balance
    _, maximum_electron_density = thermal_balance_guess(plasma)

    time_simulation = 2.0e5 * u.s
    volume = 3.0e30 * u.cm**3
    estimator_scale = (
        time_simulation.to_value(u.s)
        * volume.to_value(u.cm**3)
        * const.h.cgs.value
    )
    estimators = init_estimators_continuum(
        plasma.photo_ion_estimator.shape, len(plasma.number_density.columns)
    )
    estimators.photo_ion_estimator[:] = (
        np.asarray(plasma.photo_ion_estimator) * estimator_scale
    )
    estimators.stim_recomb_estimator[:] = (
        np.asarray(plasma.stim_recomb_estimator) * estimator_scale
    )
    estimators.bf_heating_estimator[:] = np.asarray(plasma.bf_heating_coeff)
    estimators.stim_recomb_cooling_estimator[:] = np.asarray(
        plasma.stim_recomb_cooling_coeff
    )
    estimators.ff_heating_estimator[:] = np.asarray(plasma.ff_heating_estimator)

    continuum_index = plasma.atomic_data.continuum_data.multi_index_nu_sorted
    equilibrium_levels = plasma.atomic_data.levels.loc[
        plasma.level_number_density.index
    ]
    level2continuum_edge_idx = pd.Series(
        np.arange(len(continuum_index), dtype=np.int64),
        index=continuum_index,
        name="continuum_idx",
    )
    photoionization_data = (
        plasma.atomic_data.continuum_data.photoionization_data
    )
    level_index = plasma.level_number_density.index
    hydrogen_level_positions = np.flatnonzero(
        (
            level_index.get_level_values("atomic_number")
            == plasma.nlte_species[0][0]
        )
        & (
            level_index.get_level_values("ion_number")
            == plasma.nlte_species[0][1]
        )
    )
    population_geometries = tuple(
        NumberDensityPerShell(
            plasma.number_density.loc[1, shell],
            plasma.level_number_density[shell].to_numpy(dtype=np.float64),
            hydrogen_level_positions,
        )
        for shell in plasma.number_density.columns
    )

    line_index = plasma.lines.index
    line_species_index = line_index.droplevel(
        ["level_number_lower", "level_number_upper"]
    )
    nlte_lines_mask = np.asarray(
        line_species_index.isin(plasma.nlte_species), dtype=bool
    )
    time_explosion_seconds = plasma.time_explosion
    if isinstance(time_explosion_seconds, u.Quantity):
        time_explosion_seconds = time_explosion_seconds.to_value("s")
    tau_coefficient = (
        plasma.lines.wavelength_cm.to_numpy()
        * plasma.lines.f_lu.to_numpy()
        * SOBOLEV_COEFFICIENT
        * time_explosion_seconds
    )
    sobolev_input = SobolevInputs(
        plasma.lines_lower_level_index,
        plasma.lines_upper_level_index,
        plasma.g.iloc[plasma.lines_lower_level_index].to_numpy(),
        plasma.g.iloc[plasma.lines_upper_level_index].to_numpy(),
        plasma.metastability.iloc[plasma.lines_upper_level_index].to_numpy(),
        nlte_lines_mask,
        tau_coefficient,
        np.arange(len(line_index), dtype=np.int64),
        line_index,
    )

    return PlasmaEquilibriumEvaluator(
        photoionization_data,
        level2continuum_edge_idx,
        estimators,
        time_simulation,
        volume,
        equilibrium_levels,
        plasma.ionization_data,
        RateMatrix(
            RadiativeRatesSolver(plasma.lines),
            ThermalCollisionalRateSolver(
                equilibrium_levels,
                plasma.lines,
                plasma.atomic_data.collision_data_temperatures,
                plasma.atomic_data.yg_data,
                collision_strengths_type="cmfgen",
            ),
            equilibrium_levels,
        ),
        pd.DataFrame(
            plasma.j_blues,
            index=plasma.lines.index,
            columns=plasma.number_density.columns,
        ),
        population_geometries,
        tuple(sobolev_input for _ in plasma.number_density.columns),
        plasma.level_number_density.index,
        plasma.nlte_species[0],
        plasma.number_density,
        maximum_electron_density,
        ion_population_solver=IonPopulationSolver(
            EstimatedIonRateMatrix(
                EstimatedPhotoionizationRateSolver(
                    photoionization_data,
                    level2continuum_edge_idx,
                    estimators,
                    time_simulation,
                    volume,
                ),
                CollisionalIonizationRateSolver(photoionization_data),
                plasma.phi,
            )
        ),
        ion_population_arguments={
            "radiation_field": None,
            "elemental_number_density": plasma.number_density,
            "lte_level_population": plasma.lte_level_number_density,
            "lte_ion_population": plasma.lte_ion_number_density,
            "estimated_ion_population": plasma.ion_number_density,
            "partition_function": plasma.partition_function,
            "boltzmann_factor": plasma.level_boltzmann_factor,
            "level_to_continuum_saha_factor": plasma.phi_lucy,
        },
        thermal_balance_solver=ThermalBalanceSolver(
            BoundFreeThermalRates(photoionization_data),
            FreeFreeThermalRates(),
            CollisionalIonizationThermalRates(photoionization_data),
            CollisionalBoundThermalRates(
                pd.DataFrame({"nu": np.asarray(plasma.nu_lines_coll)})
            ),
        ),
        thermal_balance_arguments={
            "collisional_ionization_rate_coefficient": plasma.coll_ion_coeff,
            "collisional_deexcitation_rate_coefficient": plasma.coll_deexc_coeff,
            "collisional_excitation_rate_coefficient": plasma.coll_exc_coeff,
            "free_free_heating_estimator": plasma.ff_heating_estimator,
            "level_population_ratio": plasma.phi_lucy,
            "bound_free_heating_estimator": plasma.bf_heating_coeff,
            "stimulated_recombination_estimator": plasma.stim_recomb_cooling_coeff,
        },
        reference_electron_temperature=plasma.t_electrons * u.K,
    )


def test_evaluator_matches_iip_five_shell_path(
    iip_plasma_after_mc: LegacyPlasmaArray,
    iip_plasma_after_thermal_balance: LegacyPlasmaArray,
    type_iip_workflow: TypeIIPWorkflow,
    iip_equilibrium_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Compare the real evaluator composition with accepted IIP shells."""
    plasma = iip_plasma_after_thermal_balance
    shell_indices = pd.Index([0, 2, 3, 8, 23])
    evaluator = iip_equilibrium_evaluator
    expected_normalized_levels = plasma.level_number_density.loc[
        plasma.nlte_species[0]
    ].divide(plasma.ion_number_density.loc[plasma.nlte_species[0]], axis=1)

    def calculate_iip_total_heating(
        state: LegacyPlasmaArray,
    ) -> npt.NDArray[np.float64]:
        thermal_balance = state.outputs_dict["fractional_heating"]
        photoionization_data = (
            state.atomic_data.continuum_data.photoionization_data
        )
        return np.array(
            [
                thermal_balance.heating_function(
                    state.t_electrons[shell],
                    state.ff_heating_estimator,
                    state.bf_heating_coeff,
                    state.stim_recomb_cooling_coeff,
                    None,
                    state.electron_densities,
                    state.ion_number_density,
                    state.level_number_density,
                    state.excitation_energy,
                    state.g,
                    state.levels,
                    state.ionization_data,
                    photoionization_data,
                    state.lines,
                    shell,
                    state.t_rad,
                    state.w,
                    state.time_explosion,
                    state.b,
                    state.previous_t_electrons,
                    state.coll_exc_cooling,
                    state.coll_deexc_heating,
                    phi_lucy=state.phi_lucy[shell],
                )[0]
                for shell in shell_indices
            ]
        )

    result = evaluator.evaluate(
        plasma.electron_densities.to_numpy(),
        plasma.t_electrons,
        expected_normalized_levels,
    )

    np.testing.assert_allclose(
        result.normalized_population.loc[:, shell_indices].to_numpy(),
        expected_normalized_levels.loc[:, shell_indices].to_numpy(),
        rtol=3e-7,
        atol=0.0,
    )
    np.testing.assert_allclose(
        result.diagnostic_ion_ratio.loc[shell_indices].to_numpy(),
        plasma.ion_ratio[shell_indices],
        rtol=3e-7,
        atol=0.0,
    )
    pd.testing.assert_frame_equal(
        result.ion_population.loc[:, shell_indices],
        plasma.ion_number_density.loc[:, shell_indices],
        rtol=1e-5,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_series_equal(
        result.charge_solved_electron_density.loc[shell_indices],
        plasma.electron_densities.loc[shell_indices],
        rtol=2e-8,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        result.absolute_level_population.loc[:, shell_indices],
        plasma.level_number_density.loc[:, shell_indices],
        rtol=1e-5,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    np.testing.assert_allclose(
        result.tau_sobolev.loc[:, shell_indices].to_numpy(),
        plasma.tau_sobolevs.loc[:, shell_indices].to_numpy(),
        rtol=1e-5,
        atol=0.0,
    )
    np.testing.assert_allclose(
        result.beta_sobolev.loc[:, shell_indices].to_numpy(),
        plasma.beta_sobolev.loc[:, shell_indices].to_numpy(),
        rtol=1e-5,
        atol=0.0,
    )
    np.testing.assert_allclose(
        result.fractional_heating.loc[shell_indices].to_numpy(),
        plasma.fractional_heating[shell_indices],
        rtol=0.0,
        atol=2e-7,
    )
    np.testing.assert_allclose(
        result.total_heating.loc[shell_indices].to_numpy(),
        calculate_iip_total_heating(plasma),
        rtol=0.0,
        atol=5e-13,
    )
    np.testing.assert_allclose(
        result.charge_residual.loc[shell_indices].to_numpy(),
        np.zeros(len(shell_indices)),
        atol=1e-10,
    )
    np.testing.assert_allclose(
        result.electron_residual.loc[shell_indices].to_numpy(),
        np.zeros(len(shell_indices)),
        atol=2e-8,
    )
    np.testing.assert_allclose(
        result.fractional_heating.loc[shell_indices].to_numpy(),
        np.zeros(len(shell_indices)),
        atol=2e-7,
    )

    off_root_candidate, off_root_maximum_density = thermal_balance_guess(
        iip_plasma_after_mc
    )
    type_iip_workflow._thermal_balance_evaluator = evaluator
    type_iip_workflow._thermal_balance_radiation_temperature = np.asarray(
        iip_plasma_after_mc.t_rad
    ).copy()
    type_iip_workflow._thermal_balance_level_seed = expected_normalized_levels
    actual_outer_residual = type_iip_workflow.thermal_balance_iteration(
        off_root_candidate, off_root_maximum_density
    )
    off_root_result = evaluator.evaluate(
        off_root_candidate[::2] * off_root_maximum_density,
        np.asarray(iip_plasma_after_mc.t_rad) * off_root_candidate[1::2],
        expected_normalized_levels,
    )
    type_iip_workflow._thermal_balance_evaluation = off_root_result
    type_iip_workflow._publish_thermal_balance_state(off_root_candidate)
    off_root_plasma = type_iip_workflow.plasma_solver
    pd.testing.assert_frame_equal(
        off_root_result.ion_population.loc[:, shell_indices],
        off_root_plasma.ion_number_density.loc[:, shell_indices],
        rtol=1e-5,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_series_equal(
        off_root_result.charge_solved_electron_density.loc[shell_indices],
        off_root_plasma.electron_densities.loc[shell_indices],
        rtol=1e-5,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    np.testing.assert_allclose(
        off_root_result.charge_residual.loc[shell_indices].to_numpy(),
        np.zeros(len(shell_indices)),
        atol=1e-10,
    )
    np.testing.assert_allclose(
        off_root_result.electron_residual.loc[shell_indices].to_numpy(),
        actual_outer_residual[::2][shell_indices],
        rtol=1e-5,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        off_root_result.fractional_heating.loc[shell_indices].to_numpy(),
        actual_outer_residual[1::2][shell_indices],
        rtol=1e-5,
        atol=0.0,
    )
    np.testing.assert_allclose(
        off_root_result.total_heating.loc[shell_indices].to_numpy(),
        calculate_iip_total_heating(iip_plasma_after_mc),
        rtol=1e-5,
        atol=0.0,
    )


def test_evaluator_is_seed_independent_for_fixed_candidates(
    iip_plasma_after_mc: LegacyPlasmaArray,
    iip_plasma_after_thermal_balance: LegacyPlasmaArray,
    iip_equilibrium_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Verify fixed-candidate results are stable for compatible level seeds.

    The accepted and post-Monte-Carlo seeds are compared away from the lower
    corner. At the lower candidate, the accepted-state seed and physical-root
    acceptance policy intentionally match legacy iip_plasma behavior.
    """
    plasma = iip_plasma_after_thermal_balance
    _, maximum_electron_density = thermal_balance_guess(plasma)

    shell_indices = pd.Index([0, 2, 3, 8, 23])
    evaluator = iip_equilibrium_evaluator

    post_mc_candidate, _ = thermal_balance_guess(iip_plasma_after_mc)
    accepted_candidate, _ = thermal_balance_guess(plasma)
    number_of_shells = len(plasma.t_rad)
    minimum_link = 1500.0 / np.min(plasma.t_rad)
    lower_bound = np.tile([0.0, minimum_link], number_of_shells)
    upper_bound = np.tile([1.0, 1.5], number_of_shells)
    candidates = {
        "near_lower": lower_bound + 0.01 * (upper_bound - lower_bound),
        "accepted": accepted_candidate,
        "post_mc": post_mc_candidate,
        "midpoint": (lower_bound + upper_bound) / 2.0,
        "near_upper": upper_bound - 0.01 * (upper_bound - lower_bound),
    }
    accepted_seed = plasma.level_number_density.loc[
        plasma.nlte_species[0]
    ].divide(plasma.ion_number_density.loc[plasma.nlte_species[0]], axis=1)
    post_mc_seed = iip_plasma_after_mc.level_number_density.loc[
        iip_plasma_after_mc.nlte_species[0]
    ].divide(
        iip_plasma_after_mc.ion_number_density.loc[
            iip_plasma_after_mc.nlte_species[0]
        ],
        axis=1,
    )

    result_fields = {
        "normalized_population",
        "diagnostic_ion_ratio",
        "charge_solved_electron_density",
        "ion_population",
        "charge_residual",
        "electron_residual",
        "total_heating",
        "fractional_heating",
    }
    closure_tolerances = {
        "charge_residual": 1e-10,
        "electron_residual": 2e-8,
        "total_heating": 5e-13,
        "fractional_heating": 2e-7,
    }
    for candidate_name, candidate in candidates.items():
        trial_density = candidate[::2] * maximum_electron_density
        temperature = np.asarray(plasma.t_rad) * candidate[1::2]
        seed_results = []
        # Legacy iip_plasma starts this lower-corner solve from the accepted
        # population and accepts a finite, nonnegative root iterate even when
        # SciPy reports failure. Phase 3 freezes that compatibility behavior
        # instead of requiring strict closure or seed independence there.
        seeds = (
            (("accepted", accepted_seed),)
            if candidate_name == "near_lower"
            else (
                ("accepted", accepted_seed),
                ("post_mc", post_mc_seed),
            )
        )
        for seed_name, level_seed in seeds:
            try:
                first_result = evaluator.evaluate(
                    trial_density, temperature, level_seed
                )
                second_result = (
                    copy(evaluator).evaluate(
                        trial_density, temperature, level_seed
                    )
                    if seed_name == "accepted"
                    else None
                )
            except ValueError as error:
                pytest.fail(
                    f"{candidate_name} candidate failed for {seed_name} seed: "
                    f"{error}"
                )
            for field in result_fields:
                first_value = getattr(first_result, field)
                if second_result is None:
                    continue
                second_value = getattr(second_result, field)
                if first_value is None:
                    assert second_value is None
                    continue
                np.testing.assert_allclose(
                    first_value.loc[:, shell_indices].to_numpy()
                    if isinstance(first_value, pd.DataFrame)
                    else first_value.loc[shell_indices].to_numpy(),
                    second_value.loc[:, shell_indices].to_numpy()
                    if isinstance(second_value, pd.DataFrame)
                    else second_value.loc[shell_indices].to_numpy(),
                    rtol=1e-5,
                    atol=closure_tolerances.get(field, 0.0),
                    err_msg=field,
                )
            seed_results.append(first_result)
        if len(seed_results) == 1:
            continue
        for field in result_fields:
            first_value = getattr(seed_results[0], field)
            second_value = getattr(seed_results[1], field)
            if first_value is None:
                assert second_value is None
                continue
            np.testing.assert_allclose(
                first_value.loc[:, shell_indices].to_numpy()
                if isinstance(first_value, pd.DataFrame)
                else first_value.loc[shell_indices].to_numpy(),
                second_value.loc[:, shell_indices].to_numpy()
                if isinstance(second_value, pd.DataFrame)
                else second_value.loc[shell_indices].to_numpy(),
                rtol=1e-5,
                atol=closure_tolerances.get(field, 0.0),
                err_msg=f"{candidate_name}: {field}",
            )


def test_standard_thermal_rates_match_iip_plasma_after_mc(
    iip_plasma_after_mc: LegacyPlasmaArray,
) -> None:
    """Compare each standard thermal rate with the IIP plasma value."""
    plasma = iip_plasma_after_mc
    photoionization_data = (
        plasma.atomic_data.continuum_data.photoionization_data
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        plasma.t_electrons * u.K,
        plasma.electron_densities.to_numpy() * u.cm**-3,
    )
    bound_free_solver = BoundFreeThermalRates(photoionization_data)
    free_free_solver = FreeFreeThermalRates()
    collisional_ionization_solver = CollisionalIonizationThermalRates(
        photoionization_data
    )
    collisional_bound_solver = CollisionalBoundThermalRates(
        pd.DataFrame({"nu": np.asarray(plasma.nu_lines_coll)})
    )

    standard_rates = {
        "bound_free": bound_free_solver.solve(
            plasma.level_number_density,
            plasma.ion_number_density,
            electron_distribution,
            plasma.phi_lucy,
            bound_free_heating_estimator=plasma.bf_heating_coeff,
            stimulated_recombination_estimator=(
                plasma.stim_recomb_cooling_coeff
            ),
        ),
        "free_free": free_free_solver.solve(
            plasma.ff_heating_estimator,
            electron_distribution,
            plasma.ion_number_density,
        ),
        "collisional_ionization": collisional_ionization_solver.solve(
            electron_distribution.number_density,
            plasma.ion_number_density,
            plasma.level_number_density,
            plasma.coll_ion_coeff,
            plasma.phi_lucy,
        ),
        "collisional_bound": collisional_bound_solver.solve(
            electron_distribution.number_density,
            plasma.coll_deexc_coeff,
            plasma.coll_exc_coeff,
            plasma.level_number_density,
        ),
    }

    thermal_balance = plasma.outputs_dict["fractional_heating"]
    legacy_rates = {
        name: np.empty((2, len(plasma.t_electrons))) for name in standard_rates
    }
    for shell, temperature in enumerate(plasma.t_electrons):
        legacy_rates["bound_free"][:, shell] = (
            thermal_balance._calculate_bf_heating_rate(
                plasma.bf_heating_coeff,
                plasma.level_number_density,
                shell,
                temperature,
                photoionization_data,
                plasma.b,
                plasma.previous_t_electrons,
            ),
            thermal_balance._calculate_fb_cooling_rate(
                temperature,
                plasma.stim_recomb_cooling_coeff,
                plasma.phi_lucy[shell],
                plasma.electron_densities,
                plasma.ion_number_density,
                photoionization_data,
                shell,
            )[0],
        )
        legacy_rates["free_free"][:, shell] = (
            thermal_balance._calculate_ff_heating_balance(
                temperature,
                plasma.ff_heating_estimator,
                plasma.electron_densities,
                plasma.ion_number_density,
                shell,
            )
        )
        legacy_rates["collisional_ionization"][:, shell] = (
            thermal_balance._calculate_coll_ion_heating_balance(
                temperature,
                photoionization_data,
                plasma.phi_lucy[shell],
                plasma.electron_densities,
                plasma.level_number_density,
                plasma.ion_number_density,
                shell,
            )
        )
        legacy_rates["collisional_bound"][:, shell] = (
            plasma.coll_deexc_heating[shell],
            plasma.coll_exc_cooling[shell],
        )

    for process, (heating, cooling) in standard_rates.items():
        relative_tolerance = {
            "bound_free": 2e-6,
            "collisional_ionization": 1e-7,
            "collisional_bound": 1e-7,
        }.get(process, 1e-12)
        np.testing.assert_allclose(
            np.vstack([heating, cooling]),
            legacy_rates[process],
            rtol=relative_tolerance,
            atol=0.0,
            err_msg=f"{process} thermal rates differ",
        )

    standard_total = ThermalBalanceSolver(
        bound_free_solver,
        free_free_solver,
        collisional_ionization_solver,
        collisional_bound_solver,
    ).solve(
        electron_distribution,
        plasma.level_number_density,
        plasma.ion_number_density,
        plasma.coll_ion_coeff,
        plasma.coll_deexc_coeff,
        plasma.coll_exc_coeff,
        plasma.ff_heating_estimator,
        plasma.phi_lucy,
        bound_free_heating_estimator=plasma.bf_heating_coeff,
        stimulated_recombination_estimator=(plasma.stim_recomb_cooling_coeff),
    )
    legacy_heating = sum(rates[0] for rates in legacy_rates.values())
    legacy_cooling = sum(rates[1] for rates in legacy_rates.values())
    np.testing.assert_allclose(
        np.vstack(standard_total),
        np.vstack(
            [
                legacy_heating - legacy_cooling,
                plasma.fractional_heating,
            ]
        ),
        rtol=3e-6,
        atol=0.0,
    )


def test_thermal_balance_solver(
    iip_regression_path,
    type_iip_workflow,
    iip_plasma_after_mc,
    regression_data,
):
    plasma = type_iip_workflow.plasma_solver
    continuum_estimators = {
        "photo_ion_estimator": iip_plasma_after_mc.photo_ion_estimator,
        "stim_recomb_estimator": iip_plasma_after_mc.stim_recomb_estimator,
        "bf_heating_estimator": iip_plasma_after_mc.bf_heating_coeff,
        "stim_recomb_cooling_estimator": (
            iip_plasma_after_mc.stim_recomb_cooling_coeff
        ),
        "ff_heating_estimator": iip_plasma_after_mc.ff_heating_estimator,
    }
    type_iip_workflow.simulation_state.t_radiative = (
        np.asarray(iip_plasma_after_mc.t_rad) * u.K
    )
    type_iip_workflow.simulation_state.dilution_factor = np.asarray(
        iip_plasma_after_mc.w
    )
    plasma.update(
        electron_densities=iip_plasma_after_mc.electron_densities,
        ion_number_density=iip_plasma_after_mc.ion_number_density,
        level_number_density=iip_plasma_after_mc.level_number_density,
        link_t_rad_t_electron=iip_plasma_after_mc.link_t_rad_t_electron,
    )
    j_blues = pd.DataFrame(
        np.asarray(iip_plasma_after_mc.j_blues),
        index=plasma.lines.index,
        columns=plasma.number_density.columns,
    )
    type_iip_workflow.solve_plasma(continuum_estimators, j_blues)
    initial_guess, max_electron_number_density = thermal_balance_guess(plasma)
    type_iip_workflow._initialize_thermal_balance_evaluator(
        max_electron_number_density
    )
    initial_residual = type_iip_workflow.thermal_balance_iteration(
        initial_guess,
        max_electron_number_density,
    )
    assert_regression_dataframe(
        regression_data,
        "thermal_balance_iteration_initial_residual",
        pd.DataFrame({"value": initial_residual}),
        rtol=1e-5,
        atol=2e-7,
    )

    type_iip_workflow.solve_thermal_balance()
    final_evaluation = type_iip_workflow._thermal_balance_evaluation
    np.testing.assert_allclose(
        final_evaluation.normalized_population.sum(axis=0),
        1.0,
        rtol=1e-12,
    )
    assert (final_evaluation.normalized_population >= 0.0).all().all()
    np.testing.assert_allclose(
        final_evaluation.trial_level_residual, 0.0, atol=1e-10
    )
    np.testing.assert_allclose(final_evaluation.level_residual, 0.0, atol=1e-10)
    np.testing.assert_allclose(
        final_evaluation.electron_residual, 0.0, atol=2e-8
    )
    np.testing.assert_allclose(
        final_evaluation.fractional_heating, 0.0, atol=2e-7
    )
    np.testing.assert_allclose(final_evaluation.total_heating, 0.0, atol=5e-13)
    np.testing.assert_allclose(
        final_evaluation.charge_residual, 0.0, atol=1e-10
    )
    type_iip_workflow.solve_continuum_state(continuum_estimators)
    type_iip_workflow.solve_opacity()

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
    print(
        "after thermal balance ion_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.ion_number_density,
                ion_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        plasma.ion_number_density,
        ion_number_density_ctardis,
        rtol=6e-7,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    # Sobolev values are stored differently between codes, so comparing raw data instead
    print(
        "after thermal balance tau_sobolevs max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow._tau_sobolev,
                tau_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        type_iip_workflow._tau_sobolev.values,
        tau_sobolevs_ctardis.values,
        # The independently closed standard level root feeds Sobolev opacity;
        # compare it with the plan-wide legacy-parity contract.
        rtol=1e-5,
        atol=0,
    )

    print(
        "after thermal balance beta_sobolev max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow._beta_sobolev,
                beta_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        type_iip_workflow._beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=1e-5,
        atol=0,
    )

    print(
        "after thermal balance level_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                plasma.level_number_density,
                level_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        plasma.level_number_density,
        level_number_density_ctardis,
        rtol=1e-5,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    assert_regression_dataframe(
        regression_data,
        "after_thermal_balance_fractional_heating",
        final_evaluation.fractional_heating,
        atol=2e-7,  # Legacy publication at the standard accepted root.
    )

    # These outputs now come from an independently closed standard state.
    # Apply one plan-wide legacy-parity contract rather than thresholds tuned
    # to individual arrays. J-blues remains an unchanged explicit input.
    regression_tolerances = {"j_blues": 1e-12}
    outputs = {
        "electron_densities": plasma.electron_densities,
        "t_electrons": plasma.t_electrons,
        "link_t_rad_t_electron": plasma.link_t_rad_t_electron,
        "p_fb_deactivation": (
            type_iip_workflow.continuum_opacity_state.p_fb_deactivation
        ),
        "chi_bf": type_iip_workflow.continuum_opacity_state.chi_bf,
        "stimulated_emission_factor": plasma.stimulated_emission_factor,
        "j_blues": plasma.j_blues.reset_index(drop=True),
    }
    for attr in STANDARD_PLASMA_SOLVER_REGRESSION_OUTPUTS:
        assert_regression_dataframe(
            regression_data,
            f"after_thermal_balance_{attr}",
            outputs[attr],
            rtol=regression_tolerances.get(attr, 1e-5),
        )

    final_guess, max_electron_number_density = thermal_balance_guess(plasma)
    residual = type_iip_workflow.thermal_balance_iteration(
        final_guess,
        max_electron_number_density,
    )
    assert_regression_dataframe(
        regression_data,
        "thermal_balance_iteration_residual",
        pd.DataFrame({"value": residual}),
        rtol=1e-5,
        atol=2e-7,  # Legacy-published and standard roots differ slightly.
    )


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
        # The standard-plasma bootstrap is not bitwise identical to the
        # removed legacy IIP owner. Use the plan's general legacy-parity
        # contract rather than tuning this threshold to one packet sample.
        rtol=1e-5,
    )


def test_iip_outer_shell_population_cutoff_second_iteration_opacity(
    tardis_regression_path: Path,
) -> None:
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
        continuum_coefficients,
        level_to_continuum_saha_factor,
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
        continuum_coefficients,
        level_to_continuum_saha_factor,
    )
    assert lte_ion_population.loc[(1, 1)].iloc[-1] == 0.0

    workflow.completed_iterations = 1
    second_iteration_opacity_states = workflow.solve_opacity()

    continuum_state = second_iteration_opacity_states[
        "opacity_state"
    ].continuum_state
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
    workflow.solve_montecarlo(second_iteration_opacity_states, 10)
    assert len(workflow.transport_state.packet_collection.output_energies) == 10
