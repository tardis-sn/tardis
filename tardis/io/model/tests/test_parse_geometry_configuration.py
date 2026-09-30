from collections.abc import Callable

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
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


@pytest.fixture
def make_specific_structure_config() -> Callable[..., Configuration]:
    """Return a builder for minimal specific-structure configurations."""

    def make_config(
        radius: dict[str, u.Quantity] | None = None,
    ) -> Configuration:
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

    return make_config


def test_homologous_config_derives_radius_from_velocity(
    make_specific_structure_config: Callable[..., Configuration],
) -> None:
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


def test_nonhomologous_config_preserves_radius_and_velocity(
    make_specific_structure_config: Callable[..., Configuration],
) -> None:
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


def test_nonhomologous_config_requires_radius(
    make_specific_structure_config: Callable[..., Configuration],
) -> None:
    """Nonhomologous configuration requires explicit radius boundaries."""
    config = make_specific_structure_config()

    with pytest.raises(
        ValueError,
        match=r"^Nonhomologous geometry requires explicit radius boundaries\.$",
    ):
        parse_nonhomologous_geometry_from_config(config)


@pytest.mark.parametrize(
    (
        "velocity_start",
        "velocity_stop",
        "radius_unit",
        "radius_boundaries",
        "expected_radius_cm",
        "expected_velocity_km_per_s",
    ),
    [
        pytest.param(
            500 * u.km / u.s,
            1_500 * u.km / u.s,
            "cm",
            [2.0e13, 8.0e13],
            [2.0e13, 8.0e13],
            [500, 1_500],
            id="one-shell",
        ),
        pytest.param(
            1_000 * u.km / u.s,
            3_000 * u.km / u.s,
            "cm",
            [1.0e14, 4.0e14, 9.0e14],
            [1.0e14, 4.0e14, 9.0e14],
            [1_000, 2_000, 3_000],
            id="two-shell-nonuniform-radius",
        ),
        pytest.param(
            4.0e7 * u.cm / u.s,
            1.6e8 * u.cm / u.s,
            "km",
            [2.0e9, 3.0e9, 5.0e9, 8.0e9],
            [2.0e14, 3.0e14, 5.0e14, 8.0e14],
            [400, 800, 1_200, 1_600],
            id="three-shell-unit-conversion",
        ),
    ],
)
def test_nonhomologous_csvy_accepts_yaml_velocity(
    velocity_start: u.Quantity,
    velocity_stop: u.Quantity,
    radius_unit: str,
    radius_boundaries: list[float],
    expected_radius_cm: list[float],
    expected_velocity_km_per_s: list[float],
) -> None:
    """A CSVY radius column can accompany YAML velocity boundaries."""
    config = Configuration({"model": {}})
    csvy_model_config = Configuration(
        {
            "velocity": {
                "start": velocity_start,
                "stop": velocity_stop,
                "num": len(radius_boundaries) - 1,
            },
            "datatype": {"fields": [{"name": "radius", "unit": radius_unit}]},
        }
    )
    csvy_model_data = pd.DataFrame({"radius": radius_boundaries})

    geometry = parse_nonhomologous_geometry_from_csvy(
        config, csvy_model_config, csvy_model_data
    )

    assert type(geometry) is Radial1DGeometry
    npt.assert_allclose(
        np.concatenate((geometry.r_inner[:1], geometry.r_outer)).to_value(u.cm),
        expected_radius_cm,
        rtol=1e-12,
        atol=0,
    )
    npt.assert_allclose(
        np.concatenate((geometry.v_inner[:1], geometry.v_outer)).to_value(
            u.km / u.s
        ),
        expected_velocity_km_per_s,
        rtol=1e-12,
        atol=0,
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
