from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_array_almost_equal

from tardis.model import SimulationState
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
    PlanckianRadiationField,
)


def test_dilution_factor_setter_accepts_read_only_initial_array() -> None:
    geometry = SimpleNamespace(
        v_inner_boundary_idx=1,
        v_outer_boundary_idx=3,
        no_of_shells_active=2,
    )
    dilution_factor = np.array([0.1, 0.2, 0.3, 0.4])
    dilution_factor.flags.writeable = False
    radiation_field_state = DilutePlanckianRadiationField(
        np.ones(4) * 10000 * u.K,
        dilution_factor,
        geometry,
    )
    simulation_state = SimulationState(
        geometry=geometry,
        composition=None,
        radiation_field_state=radiation_field_state,
        time_explosion=None,
        packet_source=None,
    )

    simulation_state.dilution_factor = np.array([0.5, 0.6])

    assert_array_almost_equal(
        simulation_state.radiation_field_state.dilution_factor,
        np.array([0.1, 0.5, 0.6, 0.4]),
    )

def test_planckian_field_temperature_kelvin() -> None:
    temperature = np.array([8000, 6000, 4000]) * u.K
    radiation_field = PlanckianRadiationField(temperature)

    assert_array_almost_equal(
        radiation_field.temperature_kelvin,
        np.array([8000, 6000, 4000]),
    )


def test_planckian_field_calculate_mean_intensity(
    atomic_dataset, regression_data
):
    temperature = np.array([8000, 6000, 4000]) * u.K
    radiation_field = PlanckianRadiationField(temperature)

    nu = atomic_dataset.lines["nu"].values
    actual_intensities = np.asarray(
        radiation_field.calculate_mean_intensity(nu)
    )

    expected_intensities = regression_data.sync_ndarray(actual_intensities)

    np.testing.assert_allclose(
        actual_intensities,
        expected_intensities,
        atol=0,
        rtol=2e-14,
    )


def test_planckian_field_rejects_negative_temperature() -> None:
    with pytest.raises(AssertionError):
        PlanckianRadiationField(np.array([-1000, 6000, 4000]) * u.K)


def test_planckian_field_rejects_zero_temperature() -> None:
    with pytest.raises(AssertionError):
        PlanckianRadiationField(np.array([0, 6000, 4000]) * u.K)


def test_planckian_field_rejects_temperature_without_units() -> None:
    with pytest.raises(u.UnitConversionError):
        PlanckianRadiationField(np.array([8000, 6000, 4000]))


def test_dilute_field_temperature_and_dilution_length_mismatch() -> None:
    with pytest.raises(AssertionError):
        DilutePlanckianRadiationField(
            np.array([8000, 6000]) * u.K,
            np.array([1.0, 1.0, 1.0]),
        )


def test_dilute_field_rejects_negative_dilution() -> None:
    with pytest.raises(AssertionError):
        DilutePlanckianRadiationField(
            np.array([8000, 6000, 4000]) * u.K,
            np.array([1.0, -0.1, 1.0]),
        )


def test_dilute_field_zero_dilution_gives_zero_intensity() -> None:
    radiation_field = DilutePlanckianRadiationField(
        np.array([8000, 6000, 4000]) * u.K,
        np.zeros(3),
    )
    frequencies = np.array([1e14, 2e14, 3e14])

    intensity = radiation_field.calculate_mean_intensity(frequencies)

    assert_array_almost_equal(intensity, np.zeros((3, 3)))


def test_dilute_field_converts_to_planckian_field() -> None:
    temperature = np.array([8000, 6000, 4000]) * u.K
    radiation_field = DilutePlanckianRadiationField(
        temperature,
        np.array([1.0, 0.5, 0.2]),
    )

    planckian_field = radiation_field.to_planckian_radiation_field()

    assert isinstance(planckian_field, PlanckianRadiationField)
    assert_array_almost_equal(
        planckian_field.temperature.to(u.K).value,
        temperature.to(u.K).value,
    )
