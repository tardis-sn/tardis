from pathlib import Path

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest

from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)

REGRESSION_DATA_FPATH = (
    Path(__file__).parent
    / "data"
    / "regression"
    / "dilute_planckian_radiation_field.h5"
)


@pytest.fixture
def radiation_field_regression_fpath():
    assert REGRESSION_DATA_FPATH.exists(), (
        f"Missing regression file: {REGRESSION_DATA_FPATH}"
    )
    return REGRESSION_DATA_FPATH


@pytest.fixture
def temperature(radiation_field_regression_fpath):
    return np.asarray(
        pd.read_hdf(radiation_field_regression_fpath, key="inputs/temperature")
    ).squeeze()


@pytest.fixture
def dilution_factor(radiation_field_regression_fpath):
    return np.asarray(
        pd.read_hdf(
            radiation_field_regression_fpath,
            key="inputs/dilution_factor",
        )
    ).squeeze()


@pytest.fixture
def nu(radiation_field_regression_fpath):
    return np.asarray(
        pd.read_hdf(radiation_field_regression_fpath, key="inputs/nu")
    ).squeeze()


@pytest.fixture
def expected_mean_intensity(radiation_field_regression_fpath):
    return np.asarray(
        pd.read_hdf(
            radiation_field_regression_fpath,
            key="outputs/mean_intensity",
        )
    )


@pytest.fixture
def radiation_field(temperature, dilution_factor):
    return DilutePlanckianRadiationField(
        temperature=temperature * u.K,
        dilution_factor=dilution_factor,
    )


def test_dilute_planckian_calculate_mean_intensity_regression(
    radiation_field,
    nu,
    expected_mean_intensity,
):
    actual = radiation_field.calculate_mean_intensity(nu)

    npt.assert_allclose(actual, expected_mean_intensity, rtol=1e-15, atol=0)


def test_dilute_planckian_temperature_kelvin_regression(
    radiation_field,
    temperature,
):
    actual = radiation_field.temperature_kelvin

    npt.assert_allclose(actual, temperature, rtol=1e-15, atol=0)
