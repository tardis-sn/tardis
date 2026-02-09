"""Unit tests for individual density solver classes."""

import numpy as np
from astropy import units as u
from numpy.testing import assert_allclose

from tardis.io.model.classic.parse_density import (
    ExponentialDensitySolver,
    PowerLawDensitySolver,
    UniformDensitySolver,
    W7DensitySolver,
)


class TestUniformDensitySolver:
    """Test UniformDensitySolver class."""

    def test_init(self):
        """Test initialization of UniformDensitySolver."""
        density_value = 1e-14 * u.g / u.cm**3
        solver = UniformDensitySolver(density_value=density_value)
        assert solver.density_value == density_value

    def test_solve(self, simple_mesh):
        """Test solve method returns correct uniform density."""
        density_value = 1e-14 * u.g / u.cm**3
        solver = UniformDensitySolver(density_value=density_value)

        result = solver.solve(simple_mesh)

        # Check shape
        n_cells = len(simple_mesh.velocity.data) - 1
        assert len(result.data) == n_cells

        # Check all values are uniform
        expected = np.ones(n_cells) * density_value
        assert_allclose(result.data.value, expected.value)
        assert result.data.unit == expected.unit

    def test_solve_preserves_units(self, simple_mesh):
        """Test that solve preserves astropy units."""
        density_value = 5e-15 * u.g / u.cm**3
        solver = UniformDensitySolver(density_value=density_value)

        result = solver.solve(simple_mesh)

        assert result.data.unit == u.g / u.cm**3


class TestPowerLawDensitySolver:
    """Test PowerLawDensitySolver class."""

    def test_init(self):
        """Test initialization of PowerLawDensitySolver."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        exponent = 2.0
        solver = PowerLawDensitySolver(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )
        assert solver.density_0 == density_0
        assert solver.velocity_0 == velocity_0
        assert solver.exponent == exponent

    def test_solve(self, simple_mesh):
        """Test solve method returns correct power law density."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        exponent = 2.0
        solver = PowerLawDensitySolver(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )

        result = solver.solve(simple_mesh)

        # Check shape
        n_cells = len(simple_mesh.velocity.data) - 1
        assert len(result.data) == n_cells

        # Check power law behavior: rho(v) = rho_0 * (v / v_0)^(-exponent)
        velocity_middle = 0.5 * (
            simple_mesh.velocity.data[:-1] + simple_mesh.velocity.data[1:]
        )
        expected = density_0 * (velocity_middle / velocity_0) ** (-exponent)
        assert_allclose(result.data.value, expected.value, rtol=1e-10)

    def test_solve_decreasing_density(self, simple_mesh):
        """Test that power law produces decreasing density with positive exponent."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        exponent = 2.0
        solver = PowerLawDensitySolver(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )

        result = solver.solve(simple_mesh)

        # Density should decrease with increasing velocity
        assert np.all(np.diff(result.data.value) < 0)

    def test_solve_preserves_units(self, simple_mesh):
        """Test that solve preserves astropy units."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        exponent = 2.0
        solver = PowerLawDensitySolver(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )

        result = solver.solve(simple_mesh)

        assert result.data.unit == u.g / u.cm**3


class TestExponentialDensitySolver:
    """Test ExponentialDensitySolver class."""

    def test_init(self):
        """Test initialization of ExponentialDensitySolver."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        solver = ExponentialDensitySolver(
            density_0=density_0, velocity_0=velocity_0
        )
        assert solver.density_0 == density_0
        assert solver.velocity_0 == velocity_0

    def test_solve(self, simple_mesh):
        """Test solve method returns correct exponential density."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        solver = ExponentialDensitySolver(
            density_0=density_0, velocity_0=velocity_0
        )

        result = solver.solve(simple_mesh)

        # Check shape
        n_cells = len(simple_mesh.velocity.data) - 1
        assert len(result.data) == n_cells

        # Check exponential behavior: rho(v) = rho_0 * exp(-v / v_0)
        velocity_middle = 0.5 * (
            simple_mesh.velocity.data[:-1] + simple_mesh.velocity.data[1:]
        )
        expected = density_0 * np.exp(-velocity_middle / velocity_0)
        assert_allclose(result.data.value, expected.value, rtol=1e-10)

    def test_solve_decreasing_density(self, simple_mesh):
        """Test that exponential produces decreasing density."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        solver = ExponentialDensitySolver(
            density_0=density_0, velocity_0=velocity_0
        )

        result = solver.solve(simple_mesh)

        # Density should decrease with increasing velocity
        assert np.all(np.diff(result.data.value) < 0)

    def test_solve_preserves_units(self, simple_mesh):
        """Test that solve preserves astropy units."""
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 1000 * u.km / u.s
        solver = ExponentialDensitySolver(
            density_0=density_0, velocity_0=velocity_0
        )

        result = solver.solve(simple_mesh)

        assert result.data.unit == u.g / u.cm**3


class TestW7DensitySolver:
    """Test W7DensitySolver class."""

    def test_init(self):
        """Test initialization of W7DensitySolver."""
        time_0 = 0.01 * u.day
        density_0 = 1e-13 * u.g / u.cm**3
        solver = W7DensitySolver(time_0=time_0, density_0=density_0)
        assert solver.time_0 == time_0
        assert solver.density_0 == density_0

    def test_solve(self, simple_mesh):
        """Test solve method returns W7 density profile."""
        time_0 = 0.01 * u.day
        density_0 = 1e-13 * u.g / u.cm**3
        solver = W7DensitySolver(time_0=time_0, density_0=density_0)

        result = solver.solve(simple_mesh)

        # Check shape
        n_cells = len(simple_mesh.velocity.data) - 1
        assert len(result.data) == n_cells

        # Check that densities are positive
        assert np.all(result.data.value > 0)

        # Check units
        assert result.data.unit == u.g / u.cm**3

    def test_solve_reasonable_values(self, simple_mesh):
        """Test that W7 density values are in reasonable range."""
        time_0 = 0.01 * u.day
        density_0 = 1e-13 * u.g / u.cm**3
        solver = W7DensitySolver(time_0=time_0, density_0=density_0)

        result = solver.solve(simple_mesh)

        # W7 densities should be in typical supernova range
        assert np.all(result.data.value < 1e-10)
        assert np.all(result.data.value > 1e-20)
