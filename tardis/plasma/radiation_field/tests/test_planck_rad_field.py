import numpy as np
import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
    PlanckianRadiationField,
)


class TestDilutePlanckianRadiationField:
    """Direct tests for DilutePlanckianRadiationField.

    This class is tested indirectly via simulation-level regression
    tests elsewhere, but had no focused unit tests before this PR.
    """

    @pytest.fixture(scope="class")
    def dilute_rad_field(self):
        return DilutePlanckianRadiationField(
            temperature=np.array([5000.0, 9000.0]) * u.K,
            dilution_factor=np.array([0.3, 0.7]),
        )

    def test_temperature_kelvin(self, dilute_rad_field):
        npt.assert_allclose(
            dilute_rad_field.temperature_kelvin,
            np.array([5000.0, 9000.0]),
        )

    def test_mean_intensity_shape(self, dilute_rad_field):
        nu = np.array([1e14, 2e14, 5e14]) * u.Hz
        result = dilute_rad_field.calculate_mean_intensity(nu)
        if hasattr(result, "value"):
            result = result.value
        assert result.shape == (3, 2)

    def test_calculate_mean_intensity(self, regression_data, dilute_rad_field):
        nu = np.array([1e14, 2e14, 5e14]) * u.Hz
        actual = dilute_rad_field.calculate_mean_intensity(nu)
        if hasattr(actual, "value"):
            actual = actual.value
        expected = regression_data.sync_ndarray(actual)
        npt.assert_allclose(actual, expected)

    def test_dilution_factor_scaling(self, dilute_rad_field):
        """Verify that mean intensity scales linearly with W."""
        nu = np.array([3e14]) * u.Hz
        intensity = dilute_rad_field.calculate_mean_intensity(nu)
        if hasattr(intensity, "value"):
            intensity = intensity.value

        planckian = PlanckianRadiationField(
            temperature=dilute_rad_field.temperature,
        )
        undiluted = planckian.calculate_mean_intensity(nu)
        if hasattr(undiluted, "value"):
            undiluted = undiluted.value

        npt.assert_allclose(
            intensity,
            dilute_rad_field.dilution_factor * undiluted,
        )

    def test_to_planckian_radiation_field(self, dilute_rad_field):
        planckian = dilute_rad_field.to_planckian_radiation_field()
        assert isinstance(planckian, PlanckianRadiationField)
        npt.assert_allclose(
            planckian.temperature.value,
            dilute_rad_field.temperature.value,
        )

    def test_negative_temperature_raises(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(
                temperature=np.array([-5000.0]) * u.K,
                dilution_factor=np.array([0.3]),
            )

    def test_mismatched_lengths_raises(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(
                temperature=np.array([5000.0, 9000.0]) * u.K,
                dilution_factor=np.array([0.3]),
            )


class TestPlanckianRadiationField:
    """Direct tests for PlanckianRadiationField."""

    @pytest.fixture(scope="class")
    def planck_rad_field(self):
        return PlanckianRadiationField(
            temperature=np.array([5000.0, 9000.0]) * u.K,
        )

    def test_temperature_kelvin(self, planck_rad_field):
        npt.assert_allclose(
            planck_rad_field.temperature_kelvin,
            np.array([5000.0, 9000.0]),
        )

    def test_mean_intensity_shape(self, planck_rad_field):
        nu = np.array([1e14, 2e14, 5e14]) * u.Hz
        result = planck_rad_field.calculate_mean_intensity(nu)
        if hasattr(result, "value"):
            result = result.value
        assert result.shape == (3, 2)

    def test_calculate_mean_intensity(self, regression_data, planck_rad_field):
        nu = np.array([1e14, 2e14, 5e14]) * u.Hz
        actual = planck_rad_field.calculate_mean_intensity(nu)
        if hasattr(actual, "value"):
            actual = actual.value
        expected = regression_data.sync_ndarray(actual)
        npt.assert_allclose(actual, expected)

    def test_negative_temperature_raises(self):
        with pytest.raises(AssertionError):
            PlanckianRadiationField(
                temperature=np.array([-5000.0]) * u.K,
            )
