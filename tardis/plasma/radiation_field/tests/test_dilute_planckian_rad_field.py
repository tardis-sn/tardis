from numpy import testing as npt
from pandas import testing as pdt

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)
def test_
def test_dilute_planckian_rad_field_collisional_sim_state(self,collisional_simulation_state): 
    expected_rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(
            collisional_simulation_state.t_radiative
        ),
    )