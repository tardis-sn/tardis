import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.io.model.csvy.data import CSVYData
from tardis.io.model.csvy.readers import load_csvy

CSVY_TEXT = """---
name: csvy_geometry
model_density_time_0: 1 day
model_isotope_time_0: 0 day
description: Three velocity boundaries, two shells.
tardis_model_config_version: v1.0
datatype:
  fields:
    - name: velocity
      unit: km/s
      desc: velocities of shell boundaries.
    - name: H
      desc: fractional H abundance

density:
  type: uniform
  value: 1e-10 g/cm^3
---
velocity,H
10000,0.0
11000,1.0
12000,1.0
"""


@pytest.fixture
def csvy_data(tmp_path) -> CSVYData:
    fname = tmp_path / "geometry.csvy"
    fname.write_text(CSVY_TEXT)
    return load_csvy(fname)


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
