"""Tests for parse_density_solver_from_density_config function."""

import pytest
from astropy import units as u

from tardis.io.model.classic.parse_density import (
    ExponentialDensitySolver,
    PowerLawDensitySolver,
    UniformDensitySolver,
    parse_density_solver_from_density_config,
)


class TestParseDensitySolverFromDensityConfig:
    """Test parse_density_solver_from_density_config function."""

    def test_uniform_density(self):
        """Test parsing uniform density configuration."""

        class DensityConfig:
            type = "uniform"
            value = 1e-14 * u.g / u.cm**3

        solver = parse_density_solver_from_density_config(DensityConfig())

        assert isinstance(solver, UniformDensitySolver)
        assert solver.density_value.value == 1e-14
        assert solver.density_value.unit == u.g / u.cm**3

    def test_power_law_density(self):
        """Test parsing power law density configuration."""

        class DensityConfig:
            type = "power_law"
            rho_0 = 1e-13 * u.g / u.cm**3
            v_0 = 1000 * u.km / u.s
            exponent = 2.0

        solver = parse_density_solver_from_density_config(DensityConfig())

        assert isinstance(solver, PowerLawDensitySolver)
        assert solver.density_0.value == 1e-13
        assert solver.velocity_0.value == 1000
        assert solver.exponent == 2.0

    def test_exponential_density(self):
        """Test parsing exponential density configuration."""

        class DensityConfig:
            type = "exponential"
            rho_0 = 1e-13 * u.g / u.cm**3
            v_0 = 1000 * u.km / u.s

        solver = parse_density_solver_from_density_config(DensityConfig())

        assert isinstance(solver, ExponentialDensitySolver)
        assert solver.density_0.value == 1e-13
        assert solver.velocity_0.value == 1000

    def test_w7_density(self):
        """Test parsing W7 density configuration."""

        class DensityConfig:
            type = "branch85_w7"
            w7_time_0 = 0.01 * u.day
            w7_v_0 = 1000.0 * u.km / u.s
            w7_rho_0 = 1e-13 * u.g / u.cm**3

        solver = parse_density_solver_from_density_config(DensityConfig())

        # W7 is implemented as a PowerLawDensitySolver with exponent -7
        assert isinstance(solver, PowerLawDensitySolver)
        assert solver.density_0.value == 1e-13
        assert solver.velocity_0.value == 1000.0
        assert solver.exponent == -7

    def test_invalid_type_raises_error(self):
        """Test that invalid density type raises ValueError."""

        class DensityConfig:
            type = "invalid_type"

        with pytest.raises(ValueError, match="Unrecognized density type"):
            parse_density_solver_from_density_config(DensityConfig())
