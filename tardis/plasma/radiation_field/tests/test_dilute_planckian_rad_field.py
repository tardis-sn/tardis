from pathlib import Path
import itertools
import pytest

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.io.configuration.config_reader import Configuration
from tardis.model.base import SimulationState


# intensity tests
# test mean_intensity for all the fields

# test the conversion for all the fields

DILUTION_FACTORS_CONSTANT = [1,1,1]
TEMPERATURE_CONSTANT = [10000,10000,10000] * u.k 


# create radiation fields with 3 shells of valid dilution parameters and constant temperature of 100K: uniform density, varied density, uniform 0 density, 
valid_dilution_factors = [[1,1,1],[0.8,0.6,0.4],[0,0,0]] 

# create radiation fields with 3 shells of valid temperature parameters and constant dilution factor of 1: uniform temperature, varied temperature
valid_temps = [[10000,10000,10000] * u.k, [8000,6000,4000] * u.k]

# create radiation field with 3 shells of negative temperatures: 
negative_temps = [-10000,-10000,-10000] * u.k 
# create radiation field with 3 shells of zero temperatures: 
zero_temps = [0,0,0] * u.k 
# create radiation field with 3 shells of no unit specified temperatures: 
no_unit_temps = [10000,10000,10000]

# create radiation field from SimulationState configured from plasma_base_test_config.yml
# this simulation state is used in test_level_populations so it seems like a good idea to test it directly
def simulation_state(new_chianti_atomic_dataset_si):
    config = Configuration.from_yaml(
        Path("tardis")
        / "plasma"
        / "tests"
        / "data"
        / "plasma_base_test_config.yml"
    )
    return SimulationState.from_config(
        config, atom_data=new_chianti_atomic_dataset_si
    )

# create radiation field from very simple valid geometry 


# create radiation field from invalid geometry 

# create radiation field from csvy model - temp defined but no dilution factor

# create radition field from csvy model -  no temp dilution factor defined


@pytest.fixture(scope="class",params=list(itertools.product(TEMPERATURE_CONSTANT,valid_dilution_factors)))
def valid_dilute_factors_test_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(DILUTION_FACTORS_CONSTANT,valid_temps)))
def valid_temperature_test_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(DILUTION_FACTORS_CONSTANT,negative_temps)))
def negative_temperature_test_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(DILUTION_FACTORS_CONSTANT,zero_temps)))
def zero_temperature_test_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(DILUTION_FACTORS_CONSTANT,no_unit_temps)))
def no_units_temperature_test_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class",params=list(itertools.product(DILUTION_FACTORS_CONSTANT,no_unit_temps)))
def no_units_temperature_test_rad_field(request):
    temperature, dilution = request.param
    return DilutePlanckianRadiationField(temperature,dilution)

@pytest.fixture(scope="class")
def simulation_state_rad_field():
    return DilutePlanckianRadiationField(simulation_state.t_radiative,dilution_factor=np.zeros_like(simulation_state.t_radiative))