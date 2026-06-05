from pathlib import Path

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.model.base import SimulationState
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


UNIFORM_10KK = [10000, 10000, 10000] * u.K
UNIFORM_ONE_DILUTION_FACTOR = np.array([1, 1, 1])
N_SHELLS = 3

negative_temps = [-10, -10, -10] * u.K
zero_temps = [0, 0, 0] * u.K
no_unit_temps = np.array([10000, 10000, 10000])

CONFIG_PATHS = [
    Path("tardis") / "plasma" / "tests" / "data" / "config_init_trad.yml"
]


@pytest.fixture(scope="class", params=CONFIG_PATHS, ids=["trad_init"])
def config(request):
    return Configuration.from_yaml(request.param)


@pytest.fixture(scope="class", params=CONFIG_PATHS, ids=["trad_init"])
def simulation_state(request, atomic_dataset):
    config = Configuration.from_yaml(request.param)
    return SimulationState.from_config(config, atom_data=atomic_dataset)


@pytest.fixture(
    scope="class",
    params=[
        pytest.param(
            (UNIFORM_10KK, UNIFORM_ONE_DILUTION_FACTOR), id="T_10kK_dilution_1"
        ),
        pytest.param(
            (UNIFORM_10KK, np.array([0.8, 0.6, 0.4])),
            id="T_10kK_dilution_0.8_0.6_0.4",
        ),
        pytest.param((UNIFORM_10KK, np.array([0, 0, 0])), id="T_10kK_dilution_0"),
        pytest.param(
            ([8000, 6000, 4000] * u.K, UNIFORM_ONE_DILUTION_FACTOR),
            id="T_8k_6k_4kK",
        ),
    ],
)
def valid_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature, dilution)


@pytest.fixture(scope="class")
def simulation_state_rad_field(simulation_state):
    return DilutePlanckianRadiationField(
        simulation_state.t_radiative,
        dilution_factor=np.zeros_like(simulation_state.t_radiative),
    )


class TestValidFields:
    def test_temperature_and_dilution_lengths(self, valid_rad_field):
        assert len(valid_rad_field.temperature) == N_SHELLS
        assert len(valid_rad_field.dilution_factor) == N_SHELLS
        assert len(valid_rad_field.dilution_factor) == len(
            valid_rad_field.temperature
        )

    def test_calculate_mean_intensity(
        self, valid_rad_field, atomic_dataset, regression_data
    ):
        nu = atomic_dataset.lines["nu"].values
        actual_intensities = valid_rad_field.calculate_mean_intensity(nu)
        expected_intensities = regression_data.sync_ndarray(actual_intensities)
        npt.assert_array_equal(actual_intensities, expected_intensities)

    def test_temperature_and_dilution_lengths_from_simulation_state(
        self, simulation_state_rad_field, simulation_state
    ):
        n_shells = len(simulation_state.t_radiative)
        assert len(simulation_state_rad_field.temperature) == n_shells
        assert len(simulation_state_rad_field.dilution_factor) == n_shells
        assert len(simulation_state_rad_field.dilution_factor) == len(
            simulation_state_rad_field.temperature
        )

    def test_calculate_mean_intensity_from_simulation_state(
        self, simulation_state_rad_field, atomic_dataset, regression_data
    ):
        nu = atomic_dataset.lines["nu"].values
        actual_intensities = np.array(
            simulation_state_rad_field.calculate_mean_intensity(nu)
        )
        expected_intensities = regression_data.sync_ndarray(actual_intensities)
        npt.assert_array_equal(actual_intensities, expected_intensities)




class TestInvalidFields:
    def test_negative_temperature(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(negative_temps, UNIFORM_ONE_DILUTION_FACTOR)

    def test_dilution_factors_negative(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(
                UNIFORM_10KK, np.array([-1, -1, -1])
            )

    def test_zero_temperature(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(zero_temps, UNIFORM_ONE_DILUTION_FACTOR)

    def test_no_units(self):
        with pytest.raises(u.UnitConversionError):
            DilutePlanckianRadiationField(no_unit_temps, UNIFORM_ONE_DILUTION_FACTOR)

    def test_dilution_factors_no_numpy(self):
        with pytest.raises(TypeError):
            DilutePlanckianRadiationField(UNIFORM_10KK, [1, 1, 1])
