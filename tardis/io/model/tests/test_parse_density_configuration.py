import numpy as np
import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.io.configuration.config_reader import ConfigurationNameSpace
from tardis.io.model.csvy.readers import load_csvy
from tardis.io.model.parse_density_configuration import (
    calculate_density_after_time,
    parse_density_from_csvy,
    parse_density_section_config,
)
from tardis.io.model.parse_geometry_configuration import (
    parse_velocity_from_csvy,
)

CSVY_WITH_DENSITY_SECTION_AND_CSV_VELOCITY = """---
name: csvy_density_section_csv_velocity
model_density_time_0: 1 day
model_isotope_time_0: 0 day
description: Density section with velocity given in the CSV section.
tardis_model_config_version: v1.0
datatype:
  fields:
    - name: velocity
      unit: km/s
      desc: velocities of shell outer boundaries.
    - name: H
      desc: fractional H abundance
    - name: He
      desc: fractional He abundance

density:
  type: branch85_w7
---
velocity,H,He
10000,0.3,0.7
11000,0.6,0.4
18000,0.6,0.4
19000,0.6,0.4
"""

CSVY_WITH_DENSITY_SECTION_AND_YAML_VELOCITY = """---
name: csvy_density_section_yaml_velocity
model_density_time_0: 1 day
model_isotope_time_0: 0 day
description: Density section with velocity given in the YAML section.
tardis_model_config_version: v1.0

velocity:
  start: 9000 km/s
  stop: 12000 km/s
  num: 5
density:
  type: branch85_w7
abundance:
  type: uniform
  H: 1.0
---
"""


def test_parse_density_from_csvy_density_section_with_csv_velocity(tmp_path):
    csvy_path = tmp_path / "density_section.csvy"
    csvy_path.write_text(CSVY_WITH_DENSITY_SECTION_AND_CSV_VELOCITY)
    csvy_data = load_csvy(csvy_path)
    time_explosion = 13 * u.day

    density = parse_density_from_csvy(
        csvy_data.model_config, csvy_data.raw_csv_data, time_explosion
    )

    velocity = parse_velocity_from_csvy(
        csvy_data.model_config, csvy_data.raw_csv_data
    )
    v_middle = 0.5 * (velocity[1:] + velocity[:-1])
    density_0, time_0 = parse_density_section_config(
        csvy_data.model_config.density, v_middle, time_explosion
    )
    expected = calculate_density_after_time(density_0, time_0, time_explosion)
    npt.assert_allclose(density.to("g/cm^3").value, expected.to("g/cm^3").value)
    assert len(density) == len(csvy_data.raw_csv_data) - 1


def test_parse_density_from_csvy_density_section_with_yaml_velocity(
    tmp_path,
):
    csvy_path = tmp_path / "density_section_yaml_velocity.csvy"
    csvy_path.write_text(CSVY_WITH_DENSITY_SECTION_AND_YAML_VELOCITY)
    csvy_data = load_csvy(csvy_path)

    density = parse_density_from_csvy(
        csvy_data.model_config, csvy_data.raw_csv_data, 13 * u.day
    )

    assert len(density) == csvy_data.model_config.velocity.num


V_MIDDLE = [1e9, 2e9] * u.cm / u.s


@pytest.mark.parametrize(
    ("density_configuration", "expected_density_0"),
    [
        (
            {"type": "uniform", "value": 2e-10 * u.g / u.cm**3},
            [2e-10, 2e-10],
        ),
        (
            {
                "type": "power_law",
                "v_0": 1e9 * u.cm / u.s,
                "rho_0": 1e-10 * u.g / u.cm**3,
                "exponent": -2,
            },
            [1e-10, 2.5e-11],
        ),
        (
            {
                "type": "exponential",
                "v_0": 1e9 * u.cm / u.s,
                "rho_0": 1e-10 * u.g / u.cm**3,
            },
            [1e-10 * np.exp(-1), 1e-10 * np.exp(-2)],
        ),
    ],
)
def test_parse_density_section_config_analytic_profiles(
    density_configuration, expected_density_0
):
    time_explosion = 13 * u.day

    density_0, time_0 = parse_density_section_config(
        ConfigurationNameSpace(density_configuration), V_MIDDLE, time_explosion
    )

    npt.assert_allclose(density_0.to("g/cm^3").value, expected_density_0)
    assert time_0 == time_explosion


def test_parse_density_section_config_rejects_unknown_type():
    with pytest.raises(ValueError, match="Unrecognized density type 'foo'"):
        parse_density_section_config(
            ConfigurationNameSpace({"type": "foo"}), V_MIDDLE, 13 * u.day
        )


def test_calculate_density_after_time_scales_with_inverse_cube_of_time():
    density = calculate_density_after_time(
        [8.0, 16.0] * u.g / u.cm**3, 1 * u.day, 2 * u.day
    )

    npt.assert_allclose(density.to("g/cm^3").value, [1.0, 2.0])
