from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

from tardis.conftest import assert_regression_dataframe
from tardis.iip_plasma.continuum.base_continuum import BaseContinuum
from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.equilibrium.continuum_state import (
    EquilibriumContinuumState,
)
from tardis.plasma.properties.ion_population import IonNumberDensity
from tardis.plasma.properties.level_population import LevelNumberDensity
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow


def _max_rel_diff(actual, expected):
    """Helper function to print relative diffs for checking broken tests"""
    actual_vals = actual.values
    expected_vals = expected.values

    relative_difference = np.abs(actual_vals - expected_vals) / np.abs(
        expected_vals
    )

    return float(np.nanmax(relative_difference))


def _as_regression_dataframe(value: object) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, pd.Series):
        return value.to_frame("value")

    values = np.asarray(value)
    if values.ndim == 1:
        return pd.DataFrame({"value": values})
    return pd.DataFrame(values)


def _assert_regression_dataframe(
    regression_data: object,
    key: str,
    actual: object,
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


@pytest.fixture
def ctardis_after_mc_continuum_estimators(
    iip_regression_path: Path,
) -> dict[str, object]:
    # Replay the stored post-Monte-Carlo inputs so the parity diagnostic
    # isolates the plasma fixed point from transport differences.
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

    return {
        "photoionization_rate_estimator": photo_ion_estimator,
        "stimulated_recombination_rate_estimator": stim_recomb_estimator,
        "bound_free_heating_estimator": bf_heating_estimator,
        "stimulated_recombination_cooling_estimator": (
            stim_recomb_cooling_estimator
        ),
        # C-TARDIS stores these 24 per-shell free-free heating values only in
        # the legacy fixture, so replay them verbatim for the parity input.
        "free_free_heating_estimator": [
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
        ],
    }


LEGACY_PLASMA_REGRESSION_OUTPUTS = (
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


INITIAL_PARITY_XFAIL = pytest.mark.xfail(
    strict=True,
    reason="Step 2 population/Sobolev parity gate is not open yet.",
)

WORKFLOW_INITIAL_REGRESSION_CASES = (
    pytest.param("transition_probabilities", marks=INITIAL_PARITY_XFAIL),
    pytest.param("ion_number_density", marks=INITIAL_PARITY_XFAIL),
    pytest.param("tau_sobolevs", marks=INITIAL_PARITY_XFAIL),
    pytest.param("beta_sobolev", marks=INITIAL_PARITY_XFAIL),
    pytest.param("level_number_density", marks=INITIAL_PARITY_XFAIL),
    pytest.param("electron_densities", marks=INITIAL_PARITY_XFAIL),
    pytest.param("p_fb_deactivation", marks=INITIAL_PARITY_XFAIL),
    pytest.param("chi_bf", marks=INITIAL_PARITY_XFAIL),
    pytest.param("stimulated_emission_factor", marks=INITIAL_PARITY_XFAIL),
    "ion_ratio",
    "t_electrons",
    "link_t_rad_t_electron",
    "j_blues",
)


# The IIP oracle uses the rounded ``8.629e-6`` collisional-rate prefactor,
# whereas the standard solver derives it from ``tardis.constants``.  This
# produces a maximum relative difference of 1.85e-5 for this fixed state.
LEGACY_COLLISIONAL_RATE_RTOL = 2e-5
# Normalizing cooling channels propagates that rate-coefficient difference; the
# maximum relative difference of the fixed legacy state is 2.13e-4.
LEGACY_COOLING_PROBABILITY_RTOL = 2.2e-4


@pytest.fixture
def iip_regression_path(tardis_regression_path):
    return tardis_regression_path / "tardis" / "workflows" / "tests"


@pytest.fixture
def thermal_balance_guess() -> Callable[
    [object], tuple[np.ndarray, np.ndarray]
]:
    def build_guess(plasma_solver: object) -> tuple[np.ndarray, np.ndarray]:
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


@pytest.fixture(scope="module")
def ctardis_compare_config(tardis_regression_path: object) -> Configuration:
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )

    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )

    # Force the NLTE configuration required by the ctardis comparison.
    config.plasma.nlte.species = [(1, 0)]
    return config


@pytest.fixture
def type_iip_workflow(
    ctardis_compare_config: Configuration,
) -> TypeIIPWorkflow:
    return TypeIIPWorkflow(ctardis_compare_config)


@pytest.fixture(scope="module")
def initial_type_iip_workflow(
    ctardis_compare_config: Configuration,
) -> TypeIIPWorkflow:
    """Return one read-only workflow for granular initialization checks."""
    return TypeIIPWorkflow(ctardis_compare_config)


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
def iip_plasma(
    iip_atom_data: object,
    elemental_number_density: pd.DataFrame,
    ctardis_compare_config: Configuration,
) -> LegacyPlasmaArray:
    """Return the restored coupled IIP plasma used as the C-TARDIS oracle."""
    return LegacyPlasmaArray(
        elemental_number_density,
        iip_atom_data,
        ctardis_compare_config.supernova.time_explosion.to("s").value,
        nlte_config=ctardis_compare_config.plasma.nlte,
        delta_treatment=None,
        ionization_mode="nlte",
        excitation_mode="dilute-lte",
        line_interaction_type=(
            ctardis_compare_config.plasma.line_interaction_type
        ),
        link_t_rad_t_electron=np.ones(24),
        helium_treatment="none",
        heating_rate_data_file=None,
        v_inner=None,
        v_outer=None,
        continuum_treatment=True,
    )


# "NLTE init" is the first call to update_radiationfield to set up the plasma
@pytest.fixture
def iip_plasma_nlte_init(
    iip_regression_path: object,
    iip_plasma: LegacyPlasmaArray,
    ctardis_compare_config: Configuration,
) -> LegacyPlasmaArray:
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
def initial_continuum_numerical_contract(
    iip_plasma_nlte_init: LegacyPlasmaArray,
) -> tuple[EquilibriumContinuumState, BaseContinuum, LegacyPlasmaArray]:
    """Build standard and legacy continuum states from identical populations."""
    plasma = iip_plasma_nlte_init
    legacy_atomic_data = plasma.atomic_data
    photoionization_data = legacy_atomic_data.continuum_data.photoionization_data
    atomic_data = SimpleNamespace(
        photoionization_data=photoionization_data,
        levels=legacy_atomic_data.levels,
        lines=legacy_atomic_data.lines,
        collision_data_temperatures=(
            legacy_atomic_data.collision_data_temperatures
        ),
        yg_data=legacy_atomic_data.yg_data,
        ionization_data=legacy_atomic_data.ionization_data,
    )
    standard_inputs = SimpleNamespace(
        atomic_data=atomic_data,
        t_electrons=plasma.t_electrons,
        electron_densities=plasma.electron_densities,
        dilute_planckian_radiation_field=DilutePlanckianRadiationField(
            plasma.t_rad * u.K,
            plasma.w,
        ),
        lte_level_number_density=plasma.lte_level_number_density,
        level_number_density=plasma.level_number_density,
        lte_ion_number_density=plasma.lte_ion_number_density,
        ion_number_density=plasma.ion_number_density,
    )
    standard_state = EquilibriumContinuumState.from_plasma(standard_inputs)
    legacy_continuum = BaseContinuum(
        atom_data=legacy_atomic_data,
        plasma_array=plasma,
        ws=plasma.w,
        radiative_transition_probabilities=plasma.transition_probabilities,
        estimators=None,
    )
    return standard_state, legacy_continuum, plasma


@pytest.mark.parametrize(
    ("standard_name", "legacy_name", "relative_tolerance"),
    [
        ("radiative_ionization_rate", "radiative_ionization", 2e-6),
        ("radiative_recombination_rate", "radiative_recombination", 1e-5),
    ],
)
def test_standard_radiative_continuum_rates_match_numerical_contract(
    initial_continuum_numerical_contract: tuple[
        EquilibriumContinuumState,
        BaseContinuum,
        LegacyPlasmaArray,
    ],
    standard_name: str,
    legacy_name: str,
    relative_tolerance: float,
) -> None:
    """Protect radiative continuum coefficients with an independent solver."""
    standard_state, legacy_continuum, _ = initial_continuum_numerical_contract
    actual = getattr(standard_state, standard_name)
    expected = getattr(legacy_continuum, legacy_name).rate_coefficient
    pd.testing.assert_index_equal(actual.index, expected.index)
    np.testing.assert_allclose(
        actual.to_numpy(),
        expected.to_numpy(),
        rtol=relative_tolerance,
        atol=0.0,
    )


def test_standard_collisional_continuum_rates_match_numerical_contract(
    initial_continuum_numerical_contract: tuple[
        EquilibriumContinuumState,
        BaseContinuum,
        LegacyPlasmaArray,
    ],
) -> None:
    """Protect all collisional continuum coefficients numerically."""
    standard_state, _, plasma = initial_continuum_numerical_contract
    np.testing.assert_allclose(
        standard_state.collisional_excitation_rate,
        plasma.coll_exc_coeff.loc[
            standard_state.collisional_excitation_rate.index
        ],
        rtol=LEGACY_COLLISIONAL_RATE_RTOL,
        atol=0.0,
    )
    np.testing.assert_allclose(
        standard_state.collisional_deexcitation_rate,
        plasma.coll_deexc_coeff.loc[
            standard_state.collisional_deexcitation_rate.index
        ],
        rtol=LEGACY_COLLISIONAL_RATE_RTOL,
        atol=0.0,
    )
    pd.testing.assert_frame_equal(
        standard_state.collisional_ionization_rate,
        plasma.coll_ion_coeff.loc[
            standard_state.collisional_ionization_rate.index
        ],
        rtol=LEGACY_COLLISIONAL_RATE_RTOL,
        check_names=False,
        check_column_type=False,
    )
    np.testing.assert_allclose(
        standard_state.collisional_recombination_rate,
        plasma.coll_recomb_coeff.loc[
            standard_state.collisional_recombination_rate.index
        ],
        rtol=LEGACY_COLLISIONAL_RATE_RTOL,
        atol=0.0,
    )


def test_structured_cooling_numerical_contract_replaces_legacy_output(
    initial_continuum_numerical_contract: tuple[
        EquilibriumContinuumState,
        BaseContinuum,
        LegacyPlasmaArray,
    ],
) -> None:
    """Protect cooling physics without exposing ``sp_fb_cooling_rates``."""
    standard_state, legacy_continuum, _ = initial_continuum_numerical_contract
    expected = legacy_continuum.cooling_rates
    np.testing.assert_allclose(
        np.vstack(
            [
                standard_state.collisional_excitation_cooling_probability,
                standard_state.collisional_ionization_cooling_probability,
                standard_state.radiative_recombination_cooling_probability,
                standard_state.free_free_cooling_probability,
            ]
        ),
        np.vstack(
            [
                expected.collisional_excitation_probability,
                expected.collisional_ionization_probability,
                expected.radiative_recombination_probability,
                expected.free_free_probability,
            ]
        ),
        rtol=LEGACY_COOLING_PROBABILITY_RTOL,
        atol=0.0,
    )
    np.testing.assert_allclose(
        standard_state.radiative_recombination_cooling_array,
        expected.radiative_recombination.probabilities_array,
        rtol=LEGACY_COOLING_PROBABILITY_RTOL,
        atol=0.0,
    )


@pytest.fixture
def iip_plasma_after_mc(
    iip_regression_path: object,
    iip_plasma_nlte_init: LegacyPlasmaArray,
    ctardis_compare_config: Configuration,
) -> LegacyPlasmaArray:
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


@pytest.mark.parametrize("attr", WORKFLOW_INITIAL_REGRESSION_CASES)
def test_type_iip_workflow_initial_plasma_regression(
    initial_type_iip_workflow: TypeIIPWorkflow,
    tardis_regression_path: Path,
    attr: str,
) -> None:
    expected = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "test_iip_workflow"
        / "test_type_iip_workflow_initial_plasma_regression.h5",
        key=f"workflow_init_{attr}",
    )
    if attr == "transition_probabilities":
        macro_atom_state = initial_type_iip_workflow.solve_macro_atom_state()
        bound_bound = macro_atom_state.transition_metadata.transition_type.isin(
            [-1, 0, 1]
        )
        # The historical regression key contains only the three bound-bound
        # channels.  The structured state intentionally retains the
        # continuum and collisional channels required by transport.
        raw_actual = macro_atom_state.transition_probabilities.loc[
            bound_bound
        ].reset_index(drop=True)
    elif attr == "p_fb_deactivation":
        raw_actual = (
            initial_type_iip_workflow.continuum_opacity_state.p_fb_deactivation
        )
    elif attr == "chi_bf":
        raw_actual = initial_type_iip_workflow.continuum_opacity_state.chi_bf
    else:
        raw_actual = getattr(initial_type_iip_workflow.plasma_solver, attr)
    if attr == "link_t_rad_t_electron":
        actual = pd.DataFrame(
            np.full(expected.shape, raw_actual),
            index=expected.index,
            columns=expected.columns,
        )
    else:
        actual = _as_regression_dataframe(raw_actual)
        if attr == "j_blues":
            actual = actual.reset_index(drop=True)
        elif attr == "transition_probabilities":
            actual.index = expected.index
            actual.columns = expected.columns
    pd.testing.assert_frame_equal(
        actual,
        expected,
        rtol=1e-10,  # Mac ARM64 tolerance
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )


def test_standard_initial_continuum_factor_matches_legacy_fixture(
    initial_type_iip_workflow: TypeIIPWorkflow,
    iip_plasma_nlte_init: LegacyPlasmaArray,
    iip_atom_data: object,
) -> None:
    """Compare the standard continuum factor with the legacy oracle."""
    assert iip_atom_data.has_collision_data is False
    expected = iip_plasma_nlte_init.get_value(
        "general_level_boltzmann_factor"
    ).loc[(1, 0)]
    actual = initial_type_iip_workflow.plasma_solver.level_boltzmann_factor.loc[
        (1, 0)
    ]
    np.testing.assert_allclose(expected.loc[0], 2.0)
    np.testing.assert_allclose(actual.loc[0], expected.loc[0])
    pd.testing.assert_frame_equal(
        actual,
        expected,
        check_names=False,
        check_dtype=False,
        rtol=1e-5,
        atol=0.0,
    )


def test_initial_continuum_level_factor_uses_dilute_lte_without_estimators(
    initial_type_iip_workflow: TypeIIPWorkflow,
) -> None:
    """Verify the initial no-estimator continuum plasma contract."""
    standard_plasma = initial_type_iip_workflow.plasma_solver
    # With no estimators, continuum excitation is configured dilute LTE rather
    # than the bound-free statistical-equilibrium update.
    # These 1e-12 comparisons are algebraic property-graph identities, not
    # integrated legacy-parity tolerances.
    pd.testing.assert_frame_equal(
        standard_plasma.level_boltzmann_factor.loc[(1, 0)],
        standard_plasma.general_level_boltzmann_factor.loc[(1, 0)],
        check_dtype=False,
        check_names=False,
        rtol=1e-12,
        atol=0.0,
    )
    expected_ion_number_density, _ = IonNumberDensity(None).calculate(
        standard_plasma.phi,
        standard_plasma.hydrogen_continuum_partition_function,
        standard_plasma.number_density,
    )
    charges = expected_ion_number_density.index.get_level_values("ion_number")
    # Charge neutrality is reconstructed from all ion stages, not just H II.
    expected_electron_densities = expected_ion_number_density.multiply(
        charges, axis=0
    ).sum(axis=0)
    pd.testing.assert_frame_equal(
        standard_plasma.ion_number_density,
        expected_ion_number_density,
        check_dtype=False,
        check_names=False,
        rtol=1e-12,
        atol=0.0,
    )
    pd.testing.assert_series_equal(
        standard_plasma.electron_densities,
        expected_electron_densities,
        check_dtype=False,
        check_names=False,
        rtol=1e-12,
        atol=0.0,
    )
    expected_level_number_density = LevelNumberDensity(None).calculate(
        standard_plasma.hydrogen_continuum_level_boltzmann_factor,
        expected_ion_number_density,
        standard_plasma.levels,
        standard_plasma.hydrogen_continuum_partition_function,
    )
    pd.testing.assert_frame_equal(
        standard_plasma.level_number_density,
        expected_level_number_density,
        check_dtype=False,
        check_names=False,
        rtol=1e-12,
        atol=0.0,
    )


def test_iip_plasma_initialization(
    iip_plasma_nlte_init, iip_regression_path
):
    """Retain the legacy initial state as numerical parity evidence."""
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


# comparison of plasma after the Monte Carlo calculations have been performed
def test_iip_plasma_after_mc(
    iip_regression_path,
    iip_plasma_after_mc,
    regression_data,
):
    """Retain the legacy post-transport state as numerical parity evidence."""
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

    for attr in LEGACY_PLASMA_REGRESSION_OUTPUTS:
        assert_regression_dataframe(
            regression_data,
            f"after_mc_{attr}",
            getattr(iip_plasma_after_mc, attr),
            rtol=4e-8,
        )


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Step 2 needs the independent linear/Picard/scalar-charge solver "
        "before the fixed post-MC parity target can be opened."
    ),
)
def test_standard_post_mc_fixed_point_parity_diagnostic(
    type_iip_workflow: TypeIIPWorkflow,
    ctardis_after_mc_continuum_estimators: dict[str, object],
    iip_plasma_after_mc: LegacyPlasmaArray,
) -> None:
    """Compare standard and legacy states from identical post-MC estimators."""
    workflow = type_iip_workflow
    workflow.simulation_state.t_radiative = iip_plasma_after_mc.t_rad * u.K
    workflow.simulation_state.dilution_factor = iip_plasma_after_mc.w

    workflow.solve_plasma(
        ctardis_after_mc_continuum_estimators,
        iip_plasma_after_mc.j_blues,
    )
    plasma = workflow.plasma_solver

    pd.testing.assert_frame_equal(
        plasma.ion_number_density,
        iip_plasma_after_mc.ion_number_density,
        rtol=1e-5,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        plasma.level_number_density,
        iip_plasma_after_mc.level_number_density,
        rtol=1e-5,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_series_equal(
        plasma.electron_densities,
        iip_plasma_after_mc.electron_densities,
        rtol=1e-5,
        atol=0.0,
        check_dtype=False,
        check_names=False,
    )
    np.testing.assert_allclose(
        plasma.tau_sobolevs.to_numpy(),
        iip_plasma_after_mc.tau_sobolevs.to_numpy(),
        rtol=1e-5,
        atol=0.0,
    )
    np.testing.assert_allclose(
        plasma.beta_sobolev.to_numpy(),
        iip_plasma_after_mc.beta_sobolev.to_numpy(),
        rtol=1e-5,
        atol=0.0,
    )


def test_standard_beta_sobolev_satisfies_escape_probability_identity(
    initial_type_iip_workflow: TypeIIPWorkflow,
) -> None:
    """Check the standard plasma's Sobolev escape probabilities analytically."""
    plasma = initial_type_iip_workflow.plasma_solver
    tau = plasma.tau_sobolevs.to_numpy()
    expected = np.divide(
        -np.expm1(-tau),
        tau,
        out=np.ones_like(tau),
        where=tau != 0.0,
    )
    np.testing.assert_allclose(
        plasma.beta_sobolev.to_numpy(),
        expected,
        rtol=2e-9,
        atol=0.0,
    )


def test_thermal_balance_solver(
    iip_regression_path: object,
    type_iip_workflow: TypeIIPWorkflow,
    thermal_balance_guess: Callable[[object], tuple[np.ndarray, np.ndarray]],
) -> None:
    """Compare the standard post-thermal-balance state to fixed references."""
    opacity_states = type_iip_workflow.solve_opacity()
    type_iip_workflow.solve_montecarlo(opacity_states, 1000)
    continuum_estimators, j_blues = type_iip_workflow.update_estimators()
    type_iip_workflow.solve_plasma(continuum_estimators, j_blues)
    type_iip_workflow.solve_continuum_state(continuum_estimators)

    initial_guess, max_electron_number_density = thermal_balance_guess(
        type_iip_workflow.plasma_solver
    )
    initial_residual = type_iip_workflow.thermal_balance_iteration(
        initial_guess,
        max_electron_number_density,
    )
    assert np.all(np.isfinite(initial_residual))

    type_iip_workflow.solve_thermal_balance()
    type_iip_workflow.solve_continuum_state(continuum_estimators)

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

    macro_atom_state = type_iip_workflow.solve_macro_atom_state()
    bound_bound = macro_atom_state.transition_metadata.transition_type.isin(
        [-1, 0, 1]
    )
    transition_probabilities = macro_atom_state.transition_probabilities.loc[
        bound_bound
    ].reset_index(drop=True)
    transition_probabilities.index = transition_probabilities_ctardis.index
    transition_probabilities.columns = transition_probabilities_ctardis.columns
    pd.testing.assert_frame_equal(
        transition_probabilities,
        transition_probabilities_ctardis,
        rtol=7e-7,
        atol=0,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.ion_number_density,
        ion_number_density_ctardis,
        rtol=6e-7,
        atol=0,
        check_dtype=False,
        check_names=False,
    )
    np.testing.assert_allclose(
        type_iip_workflow.plasma_solver.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=7e-7,
        atol=0,
    )
    np.testing.assert_allclose(
        type_iip_workflow.plasma_solver.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=7e-7,
        atol=0,
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.level_number_density,
        level_number_density_ctardis,
        rtol=6e-7,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    level_totals = type_iip_workflow.plasma_solver.level_number_density.groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    np.testing.assert_allclose(
        level_totals,
        type_iip_workflow.plasma_solver.ion_number_density,
        rtol=1e-12,
        atol=0.0,
    )

    final_guess, max_electron_number_density = thermal_balance_guess(
        type_iip_workflow.plasma_solver
    )
    residual = type_iip_workflow.thermal_balance_iteration(
        final_guess,
        max_electron_number_density,
    )
    assert np.all(np.isfinite(residual))


def test_thermal_balance_iteration_stays_finite(
    type_iip_workflow: object,
    thermal_balance_guess: Callable[[object], tuple[np.ndarray, np.ndarray]],
) -> None:
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


def test_solve_opacity_exposes_structured_macro_atom_state(
    type_iip_workflow: TypeIIPWorkflow,
) -> None:
    opacity_states = type_iip_workflow.solve_opacity()
    continuum_macro_atom_state = opacity_states["macro_atom_state"]
    assert continuum_macro_atom_state.transition_probabilities.shape[0] > 0
    assert continuum_macro_atom_state.absorbing_probability_matrix is not None


def test_standard_plasma_states_conserve_populations(
    initial_type_iip_workflow: TypeIIPWorkflow,
) -> None:
    """Check population identities on the standard structured plasma state."""
    plasma = initial_type_iip_workflow.plasma_solver
    level_totals = plasma.level_number_density.groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    np.testing.assert_allclose(
        level_totals,
        plasma.ion_number_density,
        rtol=1e-12,
        atol=0.0,
    )

    charges = plasma.ion_number_density.index.get_level_values(
        "ion_number"
    ).to_numpy()
    expected_electron_densities = plasma.ion_number_density.multiply(
        charges, axis=0
    ).sum(axis=0)
    np.testing.assert_allclose(
        plasma.electron_densities.to_numpy(),
        expected_electron_densities.to_numpy(),
        rtol=1e-12,
        atol=0.0,
    )


def test_standard_workflow_completes_two_owned_iterations(
    type_iip_workflow: TypeIIPWorkflow,
) -> None:
    """Exercise two complete standard workflow iterations with real states."""
    workflow = type_iip_workflow

    for _iteration in range(2):
        opacity_states = workflow.solve_opacity()
        assert opacity_states["macro_atom_state"].absorbing_probability_matrix is not None

        workflow.solve_montecarlo(opacity_states, 256)
        assert np.all(
            np.isfinite(
                workflow.transport_state.packet_collection.output_energies
            )
        )

        estimated_values, _ = workflow.get_convergence_estimates()
        workflow.solve_simulation_state(estimated_values)
        continuum_estimators, j_blues = workflow.update_estimators()
        workflow.solve_plasma(continuum_estimators, j_blues)
        workflow.solve_thermal_balance()
        workflow.solve_continuum_state(continuum_estimators)

        plasma = workflow.plasma_solver
        assert np.all(np.isfinite(plasma.electron_densities.to_numpy()))
        assert np.all(np.isfinite(plasma.level_number_density.to_numpy()))
        assert np.all(
            np.isfinite(
                workflow.continuum_state.radiative_recombination_rate.to_numpy()
            )
        )
        assert np.all(
            np.isfinite(workflow.continuum_opacity_state.chi_bf.to_numpy())
        )

        workflow.completed_iterations += 1

    assert workflow.completed_iterations == 2


def test_solve_montecarlo_completes(
    type_iip_workflow: TypeIIPWorkflow,
) -> None:
    opacity_states = type_iip_workflow.solve_opacity()
    type_iip_workflow.solve_montecarlo(opacity_states, 1000)

    output_energies = (
        type_iip_workflow.transport_state.packet_collection.output_energies
    )
    assert len(output_energies) == 1000
    assert np.all(np.isfinite(output_energies))


def test_solve_montecarlo(
    type_iip_workflow: TypeIIPWorkflow,
    regression_data: object,
) -> None:
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
