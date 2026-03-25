from pathlib import Path
import itertools
import pytest

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.io.configuration.config_reader import Configuration
from tardis.model.base import SimulationState
from tardis.io.atom_data import AtomData

# non np array input


DILUTION_FACTORS_CONSTANT = [[1,1,1]]
TEMPERATURE_CONSTANT = [[10000,10000,10000]] * u.K 


# create radiation fields with 3 shells of valid dilution parameters and constant temperature of 100K: uniform density, varied density, uniform 0 density, 
valid_dilution_factors = [[1,1,1],[0.8,0.6,0.4],[0,0,0]] 

# create radiation fields with 3 shells of valid temperature parameters and constant dilution factor of 1: uniform temperature, varied temperature
valid_temps = [[10000,10000,10000] * u.K, [8000,6000,4000] * u.K]

# create radiation field with 3 shells of negative temperatures: 
negative_temps = [[-10000,-10000,-10000]] * u.K 
# create radiation field with 3 shells of zero temperatures: 
zero_temps = [[0,0,0] * u.K]
# create radiation field with 3 shells of no unit specified temperatures: 
no_unit_temps = [[10000,10000,10000]]

# load atomic data todo: change tardis_regression_path
@pytest.fixture
def kurucz_atom_pure_simple_dataset(tardis_regression_path):
    atomic_data_fname = (
        tardis_regression_path / "atom_data" / "kurucz_atom_pure_simple.h5"
    )
    return AtomData.from_hdf(atomic_data_fname)

# define a SimulationState configured by plasma_base_test_config.yml
@pytest.fixture(scope="class")
def simulation_state():
    config = Configuration.from_yaml(
        Path("tardis")
        / "plasma"
        / "tests"
        / "data"
        / "plasma_base_test_config.yml"
    )
    return SimulationState.from_config(config, atom_data=kurucz_atom_pure_simple_dataset)

# create radiation field from very simple valid geometry 

# create radiation field from invalid geometry 

# create radiation field from csvy model - temp defined but no dilution factor

# create radition field from csvy model -  no temp dilution factor defined



@pytest.fixture(scope="class",params=list(itertools.product(TEMPERATURE_CONSTANT,valid_dilution_factors)))
def valid_dilute_factors_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(valid_temps,DILUTION_FACTORS_CONSTANT)))
def valid_temperature_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(negative_temps,DILUTION_FACTORS_CONSTANT)))
def negative_temperature_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(zero_temps,DILUTION_FACTORS_CONSTANT)))
def zero_temperature_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(no_unit_temps,DILUTION_FACTORS_CONSTANT)))
def no_units_temperature_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class")
def simulation_state_rad_field():
    return DilutePlanckianRadiationField(simulation_state.t_radiative,dilution_factor=np.zeros_like(simulation_state.t_radiative))

class TestValidDilutionFactors:
    def test_temp_len(valid_dilute_factors_rad_field):
        assert len(valid_dilute_factors_rad_field.temperature) == 3

    def test_dilute_factors_len(valid_dilute_factors_rad_field):
        assert len(valid_dilute_factors_rad_field.dilution_factor) == 3

    def test_dilute_factors_len_equals_temp_len(valid_dilute_factors_rad_field):
        assert len(valid_dilute_factors_rad_field.dilution_factor) == len(valid_dilute_factors_rad_field.temperature)
    
    @pytest.mark.parametrize("nu",kurucz_atom_pure_simple_dataset.lines["nu"].values)
    def test_calculate_mean_intensity_single_nu(nu):
        pass 
    
    def test_calculate_mean_intensity_nu_array(kurucz_atom_pure_simple_dataset):
        nu = kurucz_atom_pure_simple_dataset.lines["nu"].values
        pass 

class TestValidTemperatures:
    def test_temp_len(valid_temperature_rad_field):
        assert len(valid_temperature_rad_field.temperature) == 3

    def test_dilute_factors_len(valid_temperature_rad_field):
        assert len(valid_temperature_rad_field.dilution_factor) == 3

    def test_dilute_factors_len_equals_temp_len(valid_temperature_rad_field):
        assert len(valid_temperature_rad_field.dilution_factor) == len(valid_temperature_rad_field.temperature)
    
    @pytest.mark.parametrize("nu",kurucz_atom_pure_simple_dataset.lines["nu"].values)
    def test_calculate_mean_intensity_single_nu(nu):
        pass 
    
    def test_calculate_mean_intensity_nu_array(kurucz_atom_pure_simple_dataset):
        nu = kurucz_atom_pure_simple_dataset.lines["nu"].values
        actual_intensity = valid_temperature_rad_field.calculate_mean_intensity(nu)
        #expected_intensity_frame = 
        pass 
