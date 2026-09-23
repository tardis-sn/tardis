import numpy as np
import numpy.testing as npt
import pandas as pd
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.parse_geometry_configuration import (
    parse_homologous_geometry_from_config,
    parse_nonhomologous_geometry_from_config,
    parse_nonhomologous_geometry_from_csvy,
)
from tardis.model.geometry.radial1d import Radial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    HomologousRadial1DGeometry,
)


def make_specific_structure_config(
    radius: dict[str, u.Quantity] | None = None,
) -> Configuration:
    """Construct a minimal specific-structure configuration."""
    structure = {
        "type": "specific",
        "velocity": {
            "start": 1_000 * u.km / u.s,
            "stop": 4_000 * u.km / u.s,
            "num": 3,
        },
    }
    if radius is not None:
        structure["radius"] = radius
    return Configuration({"model": {"structure": structure}})


def test_homologous_config_derives_radius_from_velocity() -> None:
    """A structure without radius boundaries uses homologous expansion."""
    config = make_specific_structure_config()
    time_explosion = 10 * u.day

    geometry = parse_homologous_geometry_from_config(config, time_explosion)

    assert isinstance(geometry, HomologousRadial1DGeometry)
    npt.assert_allclose(
        geometry.r_inner.to_value(u.cm),
        (geometry.v_inner * time_explosion).to_value(u.cm),
    )
    npt.assert_allclose(
        geometry.r_outer.to_value(u.cm),
        (geometry.v_outer * time_explosion).to_value(u.cm),
    )


def test_nonhomologous_config_preserves_radius_and_velocity() -> None:
    """Explicit radius and velocity boundaries remain independent."""
    config = make_specific_structure_config(
        radius={
            "start": 1.0e14 * u.cm,
            "stop": 7.0e14 * u.cm,
        }
    )

    geometry = parse_nonhomologous_geometry_from_config(config)

    assert type(geometry) is Radial1DGeometry
    npt.assert_allclose(
        np.concatenate((geometry.r_inner[:1], geometry.r_outer)).to_value(u.cm),
        np.linspace(1.0e14, 7.0e14, 4),
    )
    npt.assert_allclose(
        np.concatenate((geometry.v_inner[:1], geometry.v_outer)).to_value(
            u.km / u.s
        ),
        np.linspace(1_000, 4_000, 4),
    )


def test_nonhomologous_csvy_accepts_yaml_velocity() -> None:
    """A CSVY radius column can accompany YAML velocity boundaries."""
    config = Configuration({"model": {}})
    csvy_model_config = Configuration(
        {
            "velocity": {
                "start": 1_000 * u.km / u.s,
                "stop": 3_000 * u.km / u.s,
                "num": 2,
            },
            "datatype": {"fields": [{"name": "radius", "unit": "cm"}]},
        }
    )
    csvy_model_data = pd.DataFrame({"radius": [1.0e14, 4.0e14, 9.0e14]})

    geometry = parse_nonhomologous_geometry_from_csvy(
        config, csvy_model_config, csvy_model_data
    )

    assert type(geometry) is Radial1DGeometry
    npt.assert_allclose(
        np.concatenate((geometry.r_inner[:1], geometry.r_outer)).to_value(u.cm),
        [1.0e14, 4.0e14, 9.0e14],
    )
    npt.assert_allclose(
        np.concatenate((geometry.v_inner[:1], geometry.v_outer)).to_value(
            u.km / u.s
        ),
        [1_000, 2_000, 3_000],
    )


def test_nonhomologous_csvy_preserves_boundary_columns() -> None:
    """CSVY radius and velocity columns remain independent."""
    config = Configuration({"model": {}})
    csvy_model_config = Configuration(
        {
            "datatype": {
                "fields": [
                    {"name": "radius", "unit": "cm"},
                    {"name": "velocity", "unit": "km/s"},
                ]
            }
        }
    )
    csvy_model_data = pd.DataFrame(
        {
            "radius": [1.0e14, 4.0e14, 9.0e14],
            "velocity": [1_000, 2_000, 3_000],
        }
    )

    geometry = parse_nonhomologous_geometry_from_csvy(
        config, csvy_model_config, csvy_model_data
    )

    assert type(geometry) is Radial1DGeometry
    npt.assert_allclose(
        np.concatenate((geometry.r_inner[:1], geometry.r_outer)).to_value(u.cm),
        csvy_model_data["radius"],
    )
    npt.assert_allclose(
        np.concatenate((geometry.v_inner[:1], geometry.v_outer)).to_value(
            u.km / u.s
        ),
        csvy_model_data["velocity"],
    )
