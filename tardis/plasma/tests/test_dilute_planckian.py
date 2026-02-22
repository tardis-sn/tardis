import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy import units as u
from tardis.plasma.radiation_field.planck_rad_field import DilutePlanckianRadiationField

def test_dilute_planckian_radiation_field(regression_data):
    # 1. Inputs with Astropy units
    t_rad = np.array([10000.0, 12000.0]) * u.K
    w = np.array([0.5, 0.4])
    nu = np.array([1e14, 2e14]) * u.Hz

    # 2. Class initialize
    radiation_field_prop = DilutePlanckianRadiationField(temperature=t_rad, dilution_factor=w)

    # 3. Calculate mean intensity
    actual_output = radiation_field_prop.calculate_mean_intensity(nu=nu)

    # 4. Save and compare regression data (FIXED: removed .value)
    expected_output = regression_data.sync_ndarray(actual_output)
    assert_allclose(actual_output, expected_output)
