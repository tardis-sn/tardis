import itertools
from pathlib import Path

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pytest

from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.io.model import parse_geometry_configuration
from tardis.model.base import SimulationState
from tardis.plasma.radiation_field import DilutePlanckianRadiationField

DILUTION_FACTORS_CONSTANT = [np.array([1, 1, 1])]
TEMPERATURE_CONSTANT = [[10000, 10000, 10000] * u.K]

valid_dilution_factors = [
    np.array([1, 1, 1]),
    np.array([0.8, 0.6, 0.4]),
    np.array([0, 0, 0]),
]
valid_temps = [[10000, 10000, 10000] * u.K, [8000, 6000, 4000] * u.K]

negative_temps = [-10, -10, -10] * u.K
zero_temps = [0, 0, 0] * u.K
no_unit_temps = np.array([10000, 10000, 10000])

CONFIG_PATHS = [
    Path("tardis") / "plasma" / "tests" / "data" / "config_init_trad.yml"
]


@pytest.fixture(scope="class")
def atom_dataset(tardis_regression_path):
    atomic_data_fname = (
        tardis_regression_path
        / "atom_data"
        / "kurucz_cd23_chianti_H_He_latest.h5"
    )
    return AtomData.from_hdf(atomic_data_fname)


@pytest.fixture(scope="class", params=CONFIG_PATHS, ids=["trad_init"])
def config(request):
    return Configuration.from_yaml(request.param)


@pytest.fixture(scope="class", params=CONFIG_PATHS, ids=["trad_init"])
def geometry(config):
    time_explosion = config.supernova.time_explosion.cgs
    return parse_geometry_configuration.parse_geometry_from_config(
        config, time_explosion
    )


@pytest.fixture(scope="class", params=CONFIG_PATHS, ids=["trad_init"])
def simulation_state(request, atom_dataset):
    config = Configuration.from_yaml(request.param)
    return SimulationState.from_config(config, atom_data=atom_dataset)


@pytest.fixture(
    scope="class",
    params=list(
        itertools.product(TEMPERATURE_CONSTANT, valid_dilution_factors)
    ),
)
def valid_dilute_factors_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature, dilution)


@pytest.fixture(
    scope="class",
    params=list(itertools.product(valid_temps, DILUTION_FACTORS_CONSTANT)),
)
def valid_temperature_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature, dilution)


@pytest.fixture(scope="class")
def simulation_state_rad_field(simulation_state):
    return DilutePlanckianRadiationField(
        simulation_state.t_radiative,
        dilution_factor=np.zeros_like(simulation_state.t_radiative),
    )


@pytest.fixture(scope="class", params=CONFIG_PATHS, ids=["trad_init"])
def geometry_rad_field(config, geometry):
    temperature = np.ones(geometry.no_of_shells) * config.plasma.initial_t_rad
    dilution = np.ones(geometry.no_of_shells)
    return DilutePlanckianRadiationField(temperature, dilution, geometry)


class TestValidFields:
    def test_temp_len_dilute_focus(self, valid_dilute_factors_rad_field):
        assert len(valid_dilute_factors_rad_field.temperature) == 3

    def test_dilute_factors_len_dilute_focus(
        self, valid_dilute_factors_rad_field
    ):
        assert len(valid_dilute_factors_rad_field.dilution_factor) == 3

    def test_dilute_factors_len_equals_temp_len_dilute_focus(
        self, valid_dilute_factors_rad_field
    ):
        assert len(valid_dilute_factors_rad_field.dilution_factor) == len(
            valid_dilute_factors_rad_field.temperature
        )

    def test_calculate_mean_intensity_dilute_focus(
        self, valid_dilute_factors_rad_field, atom_dataset, regression_data
    ):
        nu = atom_dataset.lines["nu"].values
        actual_intensities = (
            valid_dilute_factors_rad_field.calculate_mean_intensity(nu)
        )
        expected_intensities = regression_data.sync_ndarray(actual_intensities)
        npt.assert_array_equal(actual_intensities, expected_intensities)

    def test_temp_len_temp_focus(self, valid_temperature_rad_field):
        assert len(valid_temperature_rad_field.temperature) == 3

    def test_dilute_factors_len_temp_focus(self, valid_temperature_rad_field):
        assert len(valid_temperature_rad_field.dilution_factor) == 3

    def test_dilute_factors_len_equals_temp_len_temp_focus(
        self, valid_temperature_rad_field
    ):
        assert len(valid_temperature_rad_field.dilution_factor) == len(
            valid_temperature_rad_field.temperature
        )

    def test_calculate_mean_intensity_temp_focus(
        self, valid_temperature_rad_field, atom_dataset, regression_data
    ):
        nu = atom_dataset.lines["nu"].values
        actual_intensities = (
            valid_temperature_rad_field.calculate_mean_intensity(nu)
        )
        expected_intensities = regression_data.sync_ndarray(actual_intensities)
        npt.assert_array_equal(actual_intensities, expected_intensities)

    def test_temp_len_sim_state(
        self, simulation_state_rad_field, regression_data
    ):
        actual_len = np.array([len(simulation_state_rad_field.temperature)])
        expected_len = regression_data.sync_ndarray(actual_len)
        npt.assert_array_equal(actual_len, expected_len)

    def test_dilute_factors_len_sim_state(
        self, simulation_state_rad_field, regression_data
    ):
        actual_len = np.array([len(simulation_state_rad_field.dilution_factor)])
        expected_len = regression_data.sync_ndarray(actual_len)
        npt.assert_array_equal(actual_len, expected_len)

    def test_dilute_factors_len_equals_temp_len_sim_state(
        self, simulation_state_rad_field
    ):
        assert len(simulation_state_rad_field.dilution_factor) == len(
            simulation_state_rad_field.temperature
        )

    def test_calculate_mean_intensity_sim_state(
        self, simulation_state_rad_field, atom_dataset, regression_data
    ):
        nu = atom_dataset.lines["nu"].values
        actual_intensities = np.array(
            simulation_state_rad_field.calculate_mean_intensity(nu)
        )
        expected_intensities = regression_data.sync_ndarray(actual_intensities)
        npt.assert_array_equal(actual_intensities, expected_intensities)

    def test_temp_len_geometry(self, geometry_rad_field, geometry):
        assert len(geometry_rad_field.temperature) == geometry.no_of_shells

    def test_dilute_factors_len_geometry(self, geometry_rad_field, geometry):
        assert len(geometry_rad_field.dilution_factor) == geometry.no_of_shells

    def test_dilute_factors_len_equals_temp_len_geometry(
        self, geometry_rad_field
    ):
        assert len(geometry_rad_field.dilution_factor) == len(
            geometry_rad_field.temperature
        )

    def test_calculate_mean_intensity_geometry(
        self, geometry_rad_field, atom_dataset, regression_data
    ):
        nu = atom_dataset.lines["nu"].values
        actual_intensities = geometry_rad_field.calculate_mean_intensity(nu)
        expected_intensities = regression_data.sync_ndarray(actual_intensities)
        npt.assert_array_equal(actual_intensities, expected_intensities)


class TestInvalidFields:
    def test_negative_temperature(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(
                negative_temps, DILUTION_FACTORS_CONSTANT[0]
            )

    def test_dilution_factors_negative(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(
                TEMPERATURE_CONSTANT[0], np.array([-1, -1, -1])
            )

    def test_zero_temperature(self):
        with pytest.raises(AssertionError):
            DilutePlanckianRadiationField(
                zero_temps, DILUTION_FACTORS_CONSTANT[0]
            )

    def test_no_units(self):
        with pytest.raises(u.UnitConversionError):
            DilutePlanckianRadiationField(
                no_unit_temps, DILUTION_FACTORS_CONSTANT[0]
            )

    def test_dilution_factors_no_numpy(self):
        with pytest.raises(TypeError):
            DilutePlanckianRadiationField(TEMPERATURE_CONSTANT[0], [1, 1, 1])
