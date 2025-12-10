import astropy
import pytest
from astropy import units as u

from tardis.io.model.readers.snec.snec_input import (
    SNECIsotopeProfile,
    read_snec_isotope_profile,
)


@pytest.fixture
def regression_input_profile_file(regression_data):
    return regression_data.regression_data_path / "testdata" / "MESA_STIR_MESA_SNEC" / "input" / "profile8.data.iso.dat"


def test_read_snec_input_profile(regression_input_profile_file):
    profile = read_snec_isotope_profile(regression_input_profile_file)
    assert isinstance(profile, SNECIsotopeProfile)
    assert all(isinstance(p, astropy.units.Quantity) for p in [profile.radius, profile.enclosed_mass])
    assert profile.radius.unit == u.cm
    assert profile.enclosed_mass.unit == u.g
