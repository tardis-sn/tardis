from pathlib import Path

import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.io.model.csvy.data import CSVYData
from tardis.io.model.csvy.readers import load_csvy


@pytest.fixture
def csvy_data() -> CSVYData:
    return load_csvy(Path(__file__).parent / "data" / "csvy_geometry.csvy")


def test_to_geometry_splits_velocity_boundaries_into_shells(csvy_data):
    geometry = csvy_data.to_geometry(time_explosion=13 * u.day)

    npt.assert_allclose(geometry.v_inner.to("km/s").value, [10000, 11000])
    npt.assert_allclose(geometry.v_outer.to("km/s").value, [11000, 12000])


def test_to_nonhomologous_geometry_radii_follow_homologous_expansion(
    csvy_data,
):
    time_explosion = 13 * u.day
    geometry = csvy_data.to_nonhomologous_geometry(
        time_explosion=time_explosion
    )

    expected_r_inner = ([10000, 11000] * u.km / u.s * time_explosion).cgs
    npt.assert_allclose(geometry.r_inner.cgs.value, expected_r_inner.value)
