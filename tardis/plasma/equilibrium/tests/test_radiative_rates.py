from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)

REGRESSION_DATA_FPATH = (
    Path(__file__).parent / "data" / "regression" / "radiative_rates_solver.h5"
)


@pytest.fixture
def radiative_rates_regression_fpath():
    assert REGRESSION_DATA_FPATH.exists(), (
        f"Missing regression file: {REGRESSION_DATA_FPATH}"
    )
    return REGRESSION_DATA_FPATH


@pytest.fixture
def einstein_coefficients(radiative_rates_regression_fpath):
    return pd.read_hdf(
        radiative_rates_regression_fpath,
        key="inputs/einstein_coefficients",
    )


@pytest.fixture
def radiation_field(radiative_rates_regression_fpath):
    temperature = np.asarray(
        pd.read_hdf(radiative_rates_regression_fpath, key="inputs/temperature")
    ).squeeze()
    dilution_factor = np.asarray(
        pd.read_hdf(
            radiative_rates_regression_fpath,
            key="inputs/dilution_factor",
        )
    ).squeeze()

    return DilutePlanckianRadiationField(
        temperature=temperature * u.K,
        dilution_factor=dilution_factor,
    )


@pytest.fixture
def expected_rates(radiative_rates_regression_fpath):
    return pd.read_hdf(radiative_rates_regression_fpath, key="outputs/rates")


def test_radiative_rates_solver_solve_regression(
    einstein_coefficients,
    radiation_field,
    expected_rates,
):
    solver = RadiativeRatesSolver(einstein_coefficients)

    actual = solver.solve(radiation_field)

    pdt.assert_frame_equal(actual, expected_rates, atol=0, rtol=1e-15)


def test_radiative_rates_solver_missing_column_raises(einstein_coefficients):
    invalid_einstein_coefficients = einstein_coefficients.drop(columns=["B_lu"])

    with pytest.raises(AssertionError):
        RadiativeRatesSolver(invalid_einstein_coefficients)
