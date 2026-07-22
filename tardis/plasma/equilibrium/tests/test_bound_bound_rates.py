import copy

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.iip_plasma.properties.continuum import (
    CollDeexcRateCoeff as IIPCollDeexcRateCoeff,
)
from tardis.iip_plasma.properties.continuum import (
    Yg as IIPYg,
)
from tardis.iip_plasma.properties.continuum import (
    YgInterpolator as IIPYgInterpolator,
)
from tardis.io.atom_data import AtomData
from tardis.plasma.equilibrium.rates import (
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)
from tardis.util.base import intensity_black_body


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
    einstein = real_einstein_data
    radiation_field = DilutePlanckianRadiationField(
        np.array([10000.0, 12000.0]) * u.K, np.array([0.3, 0.7])
    )
    rates = RadiativeRatesSolver(einstein).solve(radiation_field)
    j_nu = radiation_field.calculate_mean_intensity(einstein.nu.values)

    upward = rates.loc[(1, 0, 0, 0, 0, 1)].to_numpy()
    downward = rates.loc[(1, 0, 0, 0, 1, 0)].to_numpy()
    # This is the direct Einstein relation, so only machine round-off is
    # expected.
    npt.assert_allclose(upward, einstein.B_lu.iloc[0] * j_nu[0], rtol=1e-12)
    npt.assert_allclose(
        downward,
        einstein.A_ul.iloc[0] + einstein.B_ul.iloc[0] * j_nu[0],
        rtol=1e-12,
    )

    # The independent Planck evaluation confirms that the rate comparison is
    # using the same cgs radiation intensity as the IIP line-rate path.
    npt.assert_allclose(
        j_nu[0],
        np.array([0.3, 0.7])
        * intensity_black_body(
            einstein.nu.values * u.Hz, np.array([10000.0, 12000.0]) * u.K
        ),
    )


def test_radiative_rates_scale_linearly_with_dilution(
    real_einstein_data: pd.DataFrame,
) -> None:
    einstein = real_einstein_data
    einstein["A_ul"] = 0.0
    field_a = DilutePlanckianRadiationField(
        np.array([10000.0]) * u.K, np.array([0.25])
    )
    field_b = DilutePlanckianRadiationField(
        np.array([10000.0]) * u.K, np.array([0.5])
    )
    rates_a = RadiativeRatesSolver(einstein).solve(field_a)
    rates_b = RadiativeRatesSolver(einstein).solve(field_b)
    # With A_ul=0, doubling dilution doubles both stimulated rates exactly.
    npt.assert_allclose(
        rates_b.to_numpy(), 2.0 * rates_a.to_numpy(), rtol=1e-12
    )


def test_tabulated_collision_strength_interpolation_matches_iip(
    nlte_atomic_dataset: AtomData,
    transition_index: pd.MultiIndex,
) -> None:
    standard_atom_data = copy.deepcopy(nlte_atomic_dataset)
    iip_atom_data = copy.deepcopy(nlte_atomic_dataset)
    lines = standard_atom_data.lines.loc[transition_index].copy()
    iip_lines = iip_atom_data.lines.loc[transition_index].copy()
    temperatures = standard_atom_data.collision_data_temperatures
    supplied_strengths = pd.DataFrame(
        standard_atom_data.yg_data.loc[transition_index].to_numpy(),
        index=transition_index,
        columns=temperatures,
    )
    test_temperatures = np.array([12500.0, 17500.0]) * u.K

    standard_solver = ThermalCollisionalRateSolver(
        standard_atom_data.levels,
        lines,
        temperatures,
        supplied_strengths,
        collision_strengths_type="cmfgen",
    )
    standard_strengths = standard_solver.calculate_collision_strengths(
        test_temperatures
    )

    iip_interpolator, allowed_index, forbidden_index = IIPYgInterpolator(
        None
    ).calculate(supplied_strengths, temperatures, iip_lines.index)
    iip_strengths = IIPYg(None).calculate(
        iip_interpolator, supplied_strengths.index, test_temperatures.value
    )

    pd.testing.assert_frame_equal(
        standard_strengths.loc[transition_index],
        iip_strengths.loc[transition_index],
        check_names=False,
        check_column_type=False,
        # Both interpolators evaluate the same tabulated values at the same
        # temperatures; this is a numerical interpolation parity check.
        rtol=1e-12,
        atol=0,
    )
    assert allowed_index.equals(transition_index)
    assert forbidden_index.empty


def test_collisional_coefficients_satisfy_detailed_balance_and_temperature_scaling(
    nlte_atomic_dataset: AtomData,
    transition_index: pd.MultiIndex,
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

    upward = rates.loc[(1, 0, 0, 0, 0, 1)].to_numpy()
    downward = rates.loc[(1, 0, 0, 0, 1, 0)].to_numpy()
    delta_energy = (
        atom_data.levels.loc[(1, 0, 1), "energy"]
        - atom_data.levels.loc[(1, 0, 0), "energy"]
    )
    g_lower = atom_data.levels.loc[(1, 0, 0), "g"]
    g_upper = atom_data.levels.loc[(1, 0, 1), "g"]
    expected_ratio = (g_upper / g_lower) * np.exp(
        -delta_energy / (const.k_B.cgs.value * temperatures)
    )
    # Use the IIP de-excitation property with the same thermal Boltzmann ratio
    # to validate the downward coefficient convention.
    lte_bf = pd.DataFrame(
        [
            [
                g_lower,
                g_upper * np.exp(-delta_energy / (const.k_B.cgs.value * t)),
            ]
            for t in temperatures
        ],
        index=temperatures,
    ).T
    lte_bf.index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    upward_frame = pd.DataFrame(
        [upward], index=transition_index, columns=temperatures
    )
    iip_downward = IIPCollDeexcRateCoeff(None).calculate(lte_bf, upward_frame)
    # The IIP de-excitation property is evaluated from the same LTE ratio;
    # differences should be limited to floating-point round-off.
    npt.assert_allclose(downward, iip_downward.to_numpy().ravel(), rtol=1e-12)
    # Detailed balance combines the statistical-weight and Boltzmann factors;
    # the tolerance allows arithmetic round-off in the exponential.
    npt.assert_allclose(upward / downward, expected_ratio, rtol=1e-10)
    assert np.all(upward > 0)
    assert np.all(downward > 0)
    strength_ratio = (
        supplied_strengths.iloc[0, 1] / supplied_strengths.iloc[0, 0]
    )
    expected_temperature_scaling = (
        strength_ratio
        * np.sqrt(temperatures[0] / temperatures[1])
        * np.exp(
            delta_energy
            / const.k_B.cgs.value
            * (1 / temperatures[0] - 1 / temperatures[1])
        )
    )
    # This is the analytic sqrt(T)^-1 and threshold-exponential scaling.
    npt.assert_allclose(
        upward[1] / upward[0], expected_temperature_scaling, rtol=1e-12
    )
