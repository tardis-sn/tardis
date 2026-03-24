import numpy as np
import pytest
import astropy.units as u

from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)


def test_length_mismatch_raises():
    """Temperature and dilution factor lengths must match."""
    temperature = np.array([5000, 6000]) * u.K
    dilution_factor = np.array([0.5])

    with pytest.raises(AssertionError):
        DilutePlanckianRadiationField(temperature, dilution_factor)


def test_negative_temperature_raises():
    """Temperature must be positive."""
    temperature = np.array([-5000, 6000]) * u.K
    dilution_factor = np.array([0.5, 0.6])

    with pytest.raises(AssertionError):
        DilutePlanckianRadiationField(temperature, dilution_factor)


def test_temperature_kelvin_property():
    """temperature_kelvin should return values in Kelvin."""
    temperature = np.array([5000, 6000]) * u.K
    dilution_factor = np.array([0.5, 0.6])

    radiation_field = DilutePlanckianRadiationField(temperature, dilution_factor)

    result = radiation_field.temperature_kelvin

    assert isinstance(result, np.ndarray)
    assert np.allclose(result, np.array([5000, 6000]))


def test_calculate_mean_intensity_regression():
    """
    Regression test for calculate_mean_intensity.

    Uses controlled input and compares against reference values
    generated from a trusted run of the solver.
    """
    temperature = np.array([5000, 6000]) * u.K
    dilution_factor = np.array([0.5, 0.6])

    radiation_field = DilutePlanckianRadiationField(
        temperature, dilution_factor
    )

    nu = np.array([1e14, 2e14])

    intensity = radiation_field.calculate_mean_intensity(nu)

    # getting reference data
    print(intensity)

    # --- Custom reference data (regression values) ---
    #replaced these values with actual output from a trusted run (intensity values)
    expected = np.array([
                [4.57549236e-06, 7.22050336e-06],
                [1.01359416e-05, 1.79098807e-05],
                ])

    # --- Assertions ---
    assert intensity.shape == expected.shape
    assert np.all(np.isfinite(intensity))
    assert np.all(intensity >= 0)

    # Regression check
    assert np.allclose(intensity, expected, rtol=1e-5)