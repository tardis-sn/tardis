"""Integration tests for density solvers with realistic meshes."""

import numpy as np
from astropy import units as u
from numpy.testing import assert_allclose

from tardis.io.model.classic.parse_density import (
    ExponentialDensitySolver,
    PowerLawDensitySolver,
    UniformDensitySolver,
)
from tardis.model.mesh import HomologousRadial1DMesh


class TestSolverIntegration:
    """Integration tests for density solvers with realistic meshes."""

    def test_uniform_density_integration(self):
        """Test uniform density solver with realistic mesh."""
        velocity = np.linspace(10000, 20000, 21) * u.km / u.s
        time_explosion = 10.0 * u.day
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        density_value = 5e-14 * u.g / u.cm**3
        solver = UniformDensitySolver(density_value=density_value)
        result = solver.solve(mesh)

        assert len(result.data) == 20
        assert_allclose(
            result.data.value, np.ones(20) * density_value.value, rtol=1e-10
        )

    def test_power_law_density_integration(self):
        """Test power law density solver with realistic mesh."""
        velocity = np.linspace(10000, 20000, 21) * u.km / u.s
        time_explosion = 10.0 * u.day
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 10000 * u.km / u.s
        exponent = 3.0
        solver = PowerLawDensitySolver(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )
        result = solver.solve(mesh)

        assert len(result.data) == 20
        # Density should decrease monotonically
        assert np.all(np.diff(result.data.value) < 0)

    def test_exponential_density_integration(self):
        """Test exponential density solver with realistic mesh."""
        velocity = np.linspace(10000, 20000, 21) * u.km / u.s
        time_explosion = 10.0 * u.day
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 5000 * u.km / u.s
        solver = ExponentialDensitySolver(
            density_0=density_0, velocity_0=velocity_0
        )
        result = solver.solve(mesh)

        assert len(result.data) == 20
        # Density should decrease monotonically
        assert np.all(np.diff(result.data.value) < 0)
