import pytest
import numpy as np
from pathlib import Path
from astropy import units as u

from tardis.io.model.readers.snec.input_profiles import read_snec_input_profile



def test_read_snec_input_profile(minimal_profile_file):
    profile = read_snec_input_profile(minimal_profile_file)
    
    assert len(profile.enclosed_mass) == 3
    assert len(profile.radius) == 3
    assert profile.isotope_mass_fraction.shape == (3, 4)
    
    assert profile.enclosed_mass.unit == u.g
    assert profile.radius.unit == u.cm
    
    assert profile.enclosed_mass[0].value == 1.0e10
    assert profile.radius[1].value == 1.5e8
    
    assert list(profile.isotope_mass_fraction.columns.get_level_values('element_number')) == [6, 8, 14, 28]
    assert list(profile.isotope_mass_fraction.columns.get_level_values('mass_number')) == [12, 16, 28, 56]
    
    assert np.isclose(profile.isotope_mass_fraction.iloc[2, 3], 0.4)
