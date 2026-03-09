import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
    PlanckianRadiationField,
)
from tardis.util.base import intensity_black_body


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def dilute_rad_field():
    """A simple DilutePlanckianRadiationField with 3 shells."""
    temperature = [10000, 12000, 15000] * u.K
    dilution_factor = np.array([0.5, 0.3, 0.1])
    return DilutePlanckianRadiationField(temperature, dilution_factor)


@pytest.fixture
def sample_nu():
    """A small array of frequencies covering optical/UV range."""
    return np.array([5e14, 1e15, 2e15]) * u.Hz


# ---------------------------------------------------------------------------
# Construction / Validation Tests
# ---------------------------------------------------------------------------


class TestDilutePlanckianRadiationFieldInit:
    """Tests for __init__ validation logic."""

    def test_init_valid(self, dilute_rad_field):
        """Object constructs successfully with valid inputs."""
        assert dilute_rad_field is not None
        assert len(dilute_rad_field.temperature) == 3
        assert len(dilute_rad_field.dilution_factor) == 3

    def test_init_stores_temperature(self, dilute_rad_field):
        """Temperature is stored on the object."""
        npt.assert_array_equal(
            dilute_rad_field.temperature.to(u.K).value,
            [10000, 12000, 15000],
        )

    def test_init_stores_dilution_factor(self, dilute_rad_field):
        """Dilution factor array is stored on the object."""
        npt.assert_array_equal(
            dilute_rad_field.dilution_factor, [0.5, 0.3, 0.1]
        )

    def test_init_mismatched_lengths_raises(self):
        """AssertionError when len(temperature) != len(dilution_factor)."""
        temperature = [10000, 12000] * u.K
        dilution_factor = np.array([0.5, 0.3, 0.1])  # length 3, not 2
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(temperature, dilution_factor)

    def test_init_zero_temperature_raises(self):
        """AssertionError when any temperature is 0 K."""
        temperature = [0, 12000, 15000] * u.K
        dilution_factor = np.array([0.5, 0.3, 0.1])
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(temperature, dilution_factor)

    def test_init_negative_temperature_raises(self):
        """AssertionError when any temperature is negative."""
        temperature = [-5000, 12000, 15000] * u.K
        dilution_factor = np.array([0.5, 0.3, 0.1])
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(temperature, dilution_factor)

    def test_init_negative_dilution_raises(self):
        """AssertionError when any dilution factor is negative."""
        temperature = [10000, 12000, 15000] * u.K
        dilution_factor = np.array([0.5, -0.1, 0.1])
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(temperature, dilution_factor)

    def test_init_zero_dilution_allowed(self):
        """Dilution factor of exactly 0 is permitted (radiation-free shell)."""
        temperature = [10000, 12000, 15000] * u.K
        dilution_factor = np.array([0.0, 0.3, 0.1])
        rad_field = DilutePlanckianRadiationField(temperature, dilution_factor)
        assert rad_field.dilution_factor[0] == 0.0


# ---------------------------------------------------------------------------
# Property Tests
# ---------------------------------------------------------------------------


class TestDilutePlanckianRadiationFieldProperties:
    """Tests for class properties."""

    def test_temperature_kelvin_returns_floats(self, dilute_rad_field):
        """temperature_kelvin returns a plain numpy array of floats, not a Quantity."""
        T_K = dilute_rad_field.temperature_kelvin
        assert isinstance(T_K, np.ndarray)
        npt.assert_array_equal(T_K, [10000.0, 12000.0, 15000.0])

    def test_temperature_kelvin_unit_conversion(self):
        """temperature_kelvin correctly converts from non-Kelvin input."""
        temperature = [0, 100, 200] * u.Celsius + 273.15 * u.K
        # Construct a simple valid field (all-positive after conversion)
        temperature_K = np.array([273.15, 373.15, 473.15]) * u.K
        dilution_factor = np.array([0.1, 0.2, 0.3])
        rad_field = DilutePlanckianRadiationField(temperature_K, dilution_factor)
        npt.assert_allclose(rad_field.temperature_kelvin, [273.15, 373.15, 473.15])


# ---------------------------------------------------------------------------
# calculate_mean_intensity Tests
# ---------------------------------------------------------------------------


class TestDilutePlanckianRadiationFieldMeanIntensity:
    """Tests for calculate_mean_intensity method."""

    def test_mean_intensity_shape(self, dilute_rad_field, sample_nu):
        """Output shape is (n_frequencies, n_shells)."""
        result = dilute_rad_field.calculate_mean_intensity(sample_nu)
        # 3 frequencies × 3 shells
        assert result.shape == (3, 3)

    def test_mean_intensity_equals_dilution_times_blackbody(
        self, dilute_rad_field, sample_nu
    ):
        """
        J_nu = W * B_nu(T).

        The diluted Planckian intensity must equal the dilution factor
        multiplied by the undiluted Planck function, element-wise.
        """
        result = dilute_rad_field.calculate_mean_intensity(sample_nu)

        nu_cgs = sample_nu.to(u.Hz).value
        T_K = dilute_rad_field.temperature_kelvin
        W = dilute_rad_field.dilution_factor

        expected = W * intensity_black_body(nu_cgs[:, np.newaxis], T_K)

        npt.assert_allclose(result, expected, rtol=1e-10)

    def test_zero_dilution_gives_zero_intensity(self, sample_nu):
        """
        Shells with W=0 must have zero mean intensity.

        This represents a radiation-free shell — the radiation field is
        completely suppressed by geometry.
        """
        temperature = [10000, 12000, 15000] * u.K
        dilution_factor = np.array([0.0, 0.0, 0.0])
        rad_field = DilutePlanckianRadiationField(temperature, dilution_factor)

        result = rad_field.calculate_mean_intensity(sample_nu)
        npt.assert_array_equal(result, 0.0)

    def test_unit_dilution_equals_pure_blackbody(self, sample_nu):
        """
        With W=1, DilutePlanckian must equal pure PlanckianRadiationField.

        This is the physical limit: the dilution factor of 1 means the
        observer is inside the photosphere, seeing the full blackbody.
        """
        temperature = [10000, 12000, 15000] * u.K
        dilution_factor = np.ones(3)

        dilute = DilutePlanckianRadiationField(temperature, dilution_factor)
        planckian = PlanckianRadiationField(temperature)

        result_dilute = dilute.calculate_mean_intensity(sample_nu)
        result_planck = planckian.calculate_mean_intensity(sample_nu)

        npt.assert_allclose(result_dilute, result_planck, rtol=1e-10)

    def test_mean_intensity_positive(self, dilute_rad_field, sample_nu):
        """All intensities must be positive for positive temperatures and W>0."""
        result = dilute_rad_field.calculate_mean_intensity(sample_nu)
        assert np.all(result > 0)

    def test_mean_intensity_scales_linearly_with_dilution(self, sample_nu):
        """
        Doubling the dilution factor must double the mean intensity.

        This directly verifies the linear scaling J = W * B_nu(T).
        """
        temperature = [10000] * u.K
        W1 = np.array([0.2])
        W2 = np.array([0.4])  # double

        rad_field_1 = DilutePlanckianRadiationField(temperature, W1)
        rad_field_2 = DilutePlanckianRadiationField(temperature, W2)

        result_1 = rad_field_1.calculate_mean_intensity(sample_nu)
        result_2 = rad_field_2.calculate_mean_intensity(sample_nu)

        npt.assert_allclose(result_2, 2 * result_1, rtol=1e-10)

    def test_mean_intensity_regression(self, dilute_rad_field, sample_nu, regression_data):
        """
        Pin calculate_mean_intensity output against regression data.

        This test ensures that any change to the radiation field calculation
        (Planck formula, dilution scaling, unit handling) is immediately caught.
        """
        result = dilute_rad_field.calculate_mean_intensity(sample_nu)
        result_df = pd.DataFrame(
            result,
            index=pd.Index(sample_nu.value, name="nu"),
        )
        expected_df = regression_data.sync_dataframe(
            result_df, key="dilute_planckian_mean_intensity"
        )
        npt.assert_allclose(
            result_df.values, expected_df.values, rtol=1e-10
        )


# ---------------------------------------------------------------------------
# Conversion Test
# ---------------------------------------------------------------------------


class TestDilutePlanckianRadiationFieldConversion:
    """Tests for to_planckian_radiation_field conversion."""

    def test_to_planckian_returns_correct_type(self, dilute_rad_field):
        """Converted object is a PlanckianRadiationField."""
        planckian = dilute_rad_field.to_planckian_radiation_field()
        assert isinstance(planckian, PlanckianRadiationField)

    def test_to_planckian_preserves_temperature(self, dilute_rad_field):
        """Converted object has the same temperature as the original."""
        planckian = dilute_rad_field.to_planckian_radiation_field()
        npt.assert_array_equal(
            planckian.temperature.to(u.K).value,
            dilute_rad_field.temperature.to(u.K).value,
        )

    def test_to_planckian_has_no_dilution(self, dilute_rad_field, sample_nu):
        """
        PlanckianRadiationField produced by conversion gives W=1 intensity.

        This is the mathematical equivalent of removing the dilution factor —
        the converted field represents the full blackbody.
        """
        planckian = dilute_rad_field.to_planckian_radiation_field()
        full_intensity = planckian.calculate_mean_intensity(sample_nu)

        # Compare to manually constructed W=1 dilute field
        unit_dilute = DilutePlanckianRadiationField(
            dilute_rad_field.temperature, np.ones(3)
        )
        unit_intensity = unit_dilute.calculate_mean_intensity(sample_nu)

        npt.assert_allclose(full_intensity, unit_intensity, rtol=1e-10)
