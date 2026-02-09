"""Shared fixtures for density tests."""

import numpy as np
import pytest
from astropy import units as u

from tardis.model.mesh import HomologousRadial1DMesh


@pytest.fixture
def simple_mesh():
    """Create a simple test mesh."""
    velocity_interfaces = np.array([1000, 2000, 3000, 4000, 5000]) * u.km / u.s
    time_explosion = 1.0 * u.day
    return HomologousRadial1DMesh.from_velocity_interfaces(
        velocity=velocity_interfaces, time_explosion=time_explosion
    )
