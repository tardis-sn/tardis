import copy

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis import constants as const
from tardis.io.atom_data import AtomData
from tardis.plasma.equilibrium.rates import ThermalCollisionalRateSolver
from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)
from tardis.util.base import intensity_black_body

REFERENCE_RADIATION_TEMPERATURES = np.array([10000.0, 12000.0]) * u.K
REFERENCE_RADIATION_DILUTION_FACTORS = np.array([0.3, 0.7])
REFERENCE_SINGLE_RADIATION_TEMPERATURE = np.array([10000.0]) * u.K
REFERENCE_LOW_DILUTION_FACTOR = np.array([0.25])
REFERENCE_HIGH_DILUTION_FACTOR = np.array([0.5])
REFERENCE_COLLISION_TEMPERATURES = np.array([12500.0, 17500.0]) * u.K


@pytest.fixture
def transition_index() -> pd.MultiIndex:
    return pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )


@pytest.fixture
def real_einstein_data(
    new_chianti_atomic_dataset: AtomData,
) -> pd.DataFrame:
    einstein = (
        new_chianti_atomic_dataset.lines.loc[(1, 0, slice(None), slice(None))]
        .iloc[:1][["A_ul", "B_ul", "B_lu", "nu"]]
        .copy()
    )
    einstein.index = pd.MultiIndex.from_tuples(
        [(1, 0, *index) for index in einstein.index],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )
    return einstein


def test_radiative_rates_match_einstein_relations(
    real_einstein_data: pd.DataFrame,
) -> None:
    radiation_field = DilutePlanckianRadiationField(
        REFERENCE_RADIATION_TEMPERATURES, REFERENCE_RADIATION_DILUTION_FACTORS
    )
    j_nu = radiation_field.calculate_mean_intensity(real_einstein_data.nu.values)
    rates = RadiativeRatesSolver(real_einstein_data).solve(
        pd.DataFrame(j_nu, index=real_einstein_data.index)
    )
    upward = rates.loc[(1, 0, 0, 0, 0, 1)].to_numpy()
    downward = rates.loc[(1, 0, 0, 0, 1, 0)].to_numpy()
    # This is the direct Einstein relation, so only machine round-off is
    # expected.
    npt.assert_allclose(
        upward, real_einstein_data.B_lu.iloc[0] * j_nu[0], rtol=1e-12
    )
    npt.assert_allclose(
        downward,
        real_einstein_data.A_ul.iloc[0]
        + real_einstein_data.B_ul.iloc[0] * j_nu[0],
        rtol=1e-12,
    )
    # The independent Planck evaluation confirms that the rate comparison is
    # using the same cgs radiation intensity as the IIP line-rate path.
    npt.assert_allclose(
        j_nu[0],
        REFERENCE_RADIATION_DILUTION_FACTORS
        * intensity_black_body(
            real_einstein_data.nu.values * u.Hz, REFERENCE_RADIATION_TEMPERATURES
        ),
    )


def test_radiative_rates_use_fixed_j_blues_and_candidate_beta(
    transition_index: pd.MultiIndex,
) -> None:
    """Apply candidate escape probabilities to both Einstein directions."""
    einstein_data = pd.DataFrame(
        {"A_ul": [7.0], "B_ul": [3.0], "B_lu": [5.0], "nu": [1.0e15]},
        index=transition_index,
    )
    shell_index = pd.Index([3, 7], name="shell")
    j_blues = pd.DataFrame(
        [[2.0e8, 5.0e8]], index=transition_index, columns=shell_index
    )
    beta_sobolev = pd.DataFrame(
        [[0.25, 0.6]], index=transition_index, columns=shell_index
    )
    rates = RadiativeRatesSolver(einstein_data).solve(j_blues, beta_sobolev)
    npt.assert_allclose(
        rates.loc[(1, 0, 0, 0, 0, 1)].to_numpy(),
        einstein_data.B_lu.iloc[0]
        * j_blues.iloc[0].to_numpy()
        * beta_sobolev.iloc[0].to_numpy(),
    )
    npt.assert_allclose(
        rates.loc[(1, 0, 0, 0, 1, 0)].to_numpy(),
        (
            einstein_data.A_ul.iloc[0]
            + einstein_data.B_ul.iloc[0] * j_blues.iloc[0].to_numpy()
        )
        * beta_sobolev.iloc[0].to_numpy(),
    )


def test_radiative_rates_scale_linearly_with_dilution(
    real_einstein_data: pd.DataFrame,
) -> None:
    einstein = real_einstein_data.copy()
    einstein["A_ul"] = 0.0
    rates_a = RadiativeRatesSolver(einstein).solve(
        pd.DataFrame(
            DilutePlanckianRadiationField(
                REFERENCE_SINGLE_RADIATION_TEMPERATURE,
                REFERENCE_LOW_DILUTION_FACTOR,
            ).calculate_mean_intensity(einstein.nu.values),
            index=einstein.index,
        )
    )
    rates_b = RadiativeRatesSolver(einstein).solve(
        pd.DataFrame(
            DilutePlanckianRadiationField(
                REFERENCE_SINGLE_RADIATION_TEMPERATURE,
                REFERENCE_HIGH_DILUTION_FACTOR,
            ).calculate_mean_intensity(einstein.nu.values),
            index=einstein.index,
        )
    )
    # With A_ul=0, doubling dilution doubles both stimulated rates exactly.
    npt.assert_allclose(rates_b.to_numpy(), 2.0 * rates_a.to_numpy(), rtol=1e-12)


def test_tabulated_collision_strength_interpolation_matches_iip(
    nlte_atomic_dataset: AtomData,
    transition_index: pd.MultiIndex,
    regression_data: RegressionData,
) -> None:
    atom_data = copy.deepcopy(nlte_atomic_dataset)
    lines = atom_data.lines.loc[transition_index].copy()
    temperatures = atom_data.collision_data_temperatures
    supplied_strengths = pd.DataFrame(
        atom_data.yg_data.loc[transition_index].to_numpy(),
        index=transition_index,
        columns=temperatures,
    )
    actual = ThermalCollisionalRateSolver(
        atom_data.levels,
        lines,
        temperatures,
        supplied_strengths,
        collision_strengths_type="cmfgen",
    ).calculate_collision_strengths(REFERENCE_COLLISION_TEMPERATURES)
    actual = actual.loc[transition_index]
    expected = regression_data.sync_dataframe(actual, key="frame_0")
    # The reference is the IIP interpolation at these tabulated strengths.
    pdt.assert_frame_equal(
        actual,
        expected,
        check_names=False,
        check_column_type=False,
        rtol=1e-12,
        atol=0.0,
    )


def test_collisional_coefficients_satisfy_detailed_balance_and_temperature_scaling(
    nlte_atomic_dataset: AtomData,
    transition_index: pd.MultiIndex,
    regression_data: RegressionData,
) -> None:
    atom_data = copy.deepcopy(nlte_atomic_dataset)
    lines = atom_data.lines.loc[transition_index].copy()
    temperatures = atom_data.collision_data_temperatures[[2, 3]]
    supplied_strengths = pd.DataFrame(
        atom_data.yg_data.loc[transition_index, [2, 3]].to_numpy(),
        index=transition_index,
        columns=temperatures,
    )
    rates = ThermalCollisionalRateSolver(
        atom_data.levels,
        lines,
        temperatures,
        supplied_strengths,
        collision_strengths_type="cmfgen",
    ).solve(temperatures * u.K)
    downward = rates.loc[(1, 0, 0, 0, 1, 0)].to_numpy()
    expected = regression_data.sync_dataframe(
        pd.DataFrame({"value": downward}), key="allclose_0"
    ).to_numpy().ravel()
    npt.assert_allclose(downward, expected, rtol=1e-12)
    upward = rates.loc[(1, 0, 0, 0, 0, 1)].to_numpy()
    delta_energy = (
        atom_data.levels.loc[(1, 0, 1), "energy"]
        - atom_data.levels.loc[(1, 0, 0), "energy"]
    )
    expected_ratio = (
        atom_data.levels.loc[(1, 0, 1), "g"]
        / atom_data.levels.loc[(1, 0, 0), "g"]
    ) * np.exp(-delta_energy / (const.k_B.cgs.value * temperatures))
    # Detailed balance combines statistical weights and the Boltzmann factor.
    npt.assert_allclose(upward / downward, expected_ratio, rtol=1e-10)
