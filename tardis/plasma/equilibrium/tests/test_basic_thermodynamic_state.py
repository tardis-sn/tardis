from typing import Any

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis import constants as const
from tardis.plasma.properties.atomic import IonizationData, Levels
from tardis.plasma.properties.general import (
    BetaElectron,
    BetaRadiation,
    ElectronTemperature,
    GElectron,
)
from tardis.util.base import intensity_black_body


@pytest.fixture
def atomic_property_values(
    basic_thermodynamic_state: dict[str, Any],
) -> dict[str, Any]:
    state = basic_thermodynamic_state
    atom_data = state["atomic_data"]
    selected_atoms = state["selected_atoms"]
    return {
        "levels": Levels(None).calculate(atom_data, selected_atoms),
        "ionization": IonizationData(None).calculate(atom_data, selected_atoms),
        "mass": atom_data.atom_data.loc[selected_atoms, "mass"],
    }


def test_atomic_levels_match_iip(
    atomic_property_values: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    levels = atomic_property_values["levels"]
    actual_index = levels[0].to_frame(index=False)
    expected_index = regression_data.sync_dataframe(actual_index, key="index_0")
    pdt.assert_frame_equal(actual_index, expected_index, check_names=False)

    for assertion_idx, actual in enumerate(levels[1:]):
        actual_frame = actual.to_frame("value")
        expected = regression_data.sync_dataframe(
            actual_frame, key=f"series_{assertion_idx}"
        )
        pdt.assert_frame_equal(actual_frame, expected, check_names=False)


def test_ionization_data_matches_iip(
    atomic_property_values: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    actual = atomic_property_values["ionization"].to_frame("value")
    expected = regression_data.sync_dataframe(actual, key="series_0")
    pdt.assert_frame_equal(actual, expected, check_names=False)


def test_atomic_mass_matches_iip(
    atomic_property_values: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    actual = atomic_property_values["mass"].to_frame("value")
    expected = regression_data.sync_dataframe(actual, key="series_0")
    pdt.assert_frame_equal(actual, expected, check_names=False)


def test_number_density_and_mass_reconstruct_density(
    basic_thermodynamic_state: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    state = basic_thermodynamic_state
    masses = state["atomic_data"].atom_data.loc[state["selected_atoms"], "mass"]
    number_density = (
        state["abundance"].mul(state["density"], axis=1).div(masses, axis=0)
    )
    expected = regression_data.sync_dataframe(
        number_density, key="number_density"
    )
    pdt.assert_frame_equal(
        number_density,
        expected,
        check_names=False,
    )
    reconstructed_density = number_density.mul(masses, axis=0).sum(axis=0)
    npt.assert_allclose(
        reconstructed_density.to_numpy(), state["density"].to_numpy()
    )
    assert (number_density >= 0).all().all()


@pytest.fixture
def thermodynamic_property_values(
    basic_thermodynamic_state: dict[str, Any],
) -> dict[str, Any]:
    state = basic_thermodynamic_state
    t_rad = state["t_rad"].to_numpy()
    link = state["link_t_rad_t_electron"]
    t_electrons = ElectronTemperature(None).calculate(t_rad, link)
    beta_rad = BetaRadiation(None).calculate(t_rad)
    return {
        "t_rad": t_rad,
        "link": link,
        "t_electrons": t_electrons,
        "beta_rad": beta_rad,
        "beta_electron": BetaElectron(None).calculate(t_electrons),
        "g_electron": GElectron(None).calculate(beta_rad),
    }


def expected_array(
    regression_data: RegressionData,
    actual: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    expected = regression_data.sync_dataframe(
        pd.DataFrame({"value": actual}), key="allclose_0"
    )
    return expected.to_numpy().ravel()


def test_electron_temperature_matches_iip(
    thermodynamic_property_values: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    values = thermodynamic_property_values
    actual = values["t_electrons"]
    npt.assert_allclose(actual, expected_array(regression_data, actual))
    npt.assert_allclose(actual, values["link"] * values["t_rad"])


def test_radiation_beta_matches_iip(
    thermodynamic_property_values: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    actual = thermodynamic_property_values["beta_rad"]
    npt.assert_allclose(
        actual, expected_array(regression_data, actual), rtol=3e-7
    )


def test_electron_beta_matches_iip(
    thermodynamic_property_values: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    actual = thermodynamic_property_values["beta_electron"]
    npt.assert_allclose(
        actual, expected_array(regression_data, actual), rtol=3e-7
    )


def test_electron_statistical_factor_matches_iip(
    thermodynamic_property_values: dict[str, Any],
    regression_data: RegressionData,
) -> None:
    actual = thermodynamic_property_values["g_electron"]
    #  iip_plasma uses raw astropy constants not tardis.constants
    npt.assert_allclose(
        actual, expected_array(regression_data, actual), rtol=5e-7
    )
    expected_g = (
        2
        * np.pi
        * const.m_e.cgs.value
        / thermodynamic_property_values["beta_rad"]
        / const.h.cgs.value**2
    ) ** 1.5
    npt.assert_allclose(actual, expected_g)


def test_dilute_planckian_mean_intensity_matches_analytic_planck_function(
    basic_thermodynamic_state: dict[str, Any],
) -> None:
    state = basic_thermodynamic_state
    frequencies = (
        state["atomic_data"]
        .lines.loc[(1, 0, slice(None), slice(None)), "nu"]
        .iloc[:3]
        .to_numpy()
        * u.Hz
    )
    actual = state["radiation_field"].calculate_mean_intensity(frequencies)
    expected = state["dilution_factor"].to_numpy() * intensity_black_body(
        frequencies[np.newaxis].T, state["t_rad"].to_numpy() * u.K
    )
    # ``intensity_black_body`` and the radiation-field API return cgs values
    # without an Astropy unit wrapper.
    npt.assert_allclose(actual, expected)
