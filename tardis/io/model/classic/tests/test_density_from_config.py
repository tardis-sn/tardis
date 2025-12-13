"""Tests for density solvers using actual YAML configuration files."""

import numpy as np
from astropy import units as u
from numpy.testing import assert_allclose

from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.classic.parse_density import (
    parse_density_solver_from_density_config,
)
from tardis.model.mesh import HomologousRadial1DMesh


class TestDensityFromConfig:
    """Test density solvers using actual YAML configuration files."""

    def test_w7_density_from_config(self, example_configuration_dir):
        """Test W7 (branch85_w7) density profile from configuration file."""
        config = Configuration.from_yaml(
            example_configuration_dir / "paper1_tardis_configv1.yml"
        )

        # Extract density config and create mesh
        density_config = config.model.structure.density
        velocity_config = config.model.structure.velocity
        time_explosion = config.supernova.time_explosion

        # Create solver from config
        solver = parse_density_solver_from_density_config(density_config)

        # Create mesh matching config
        velocity = np.linspace(
            velocity_config.start.to(u.km / u.s).value,
            velocity_config.stop.to(u.km / u.s).value,
            velocity_config.num + 1
        ) * u.km / u.s
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        # Solve for density
        result = solver.solve(mesh)

        # Verify density was calculated
        assert len(result.data) == velocity_config.num
        assert result.data.unit == u.g / u.cm**3

        # Check against known values from tardis/model/tests/test_base.py
        # TestModelFromPaper1Config.test_densities() (lines 45-52)
        # First cell density
        assert_allclose(
            result.data[0].cgs.value,
            (7.542803599143591e-14 * u.Unit("g/cm^3")).value,
        )
        # Last cell density
        assert_allclose(
            result.data[-1].cgs.value,
            (1.432259798833509e-15 * u.Unit("g/cm^3")).value,
        )

    def test_uniform_density_from_config(self, example_configuration_dir):
        """Test uniform density profile from configuration file."""
        config = Configuration.from_yaml(
            example_configuration_dir / "tardis_configv1_uniform_density.yml"
        )

        # Extract density config and create mesh
        density_config = config.model.structure.density
        velocity_config = config.model.structure.velocity
        time_explosion = config.supernova.time_explosion

        # Create solver from config
        solver = parse_density_solver_from_density_config(density_config)

        # Create mesh matching config
        velocity = np.linspace(
            velocity_config.start.to(u.km / u.s).value,
            velocity_config.stop.to(u.km / u.s).value,
            velocity_config.num + 1
        ) * u.km / u.s
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        # Solve for density
        result = solver.solve(mesh)

        # Verify uniform density
        assert len(result.data) == velocity_config.num
        expected_density = 1e-14 * u.g / u.cm**3
        assert_allclose(
            result.data.value,
            np.ones(velocity_config.num) * expected_density.value,
            rtol=1e-10,
        )

    def test_power_law_density_from_config(self, example_configuration_dir):
        """Test W7 density profile from verysimple YAML configuration."""
        # Use verysimple config which has branch85_w7
        config = Configuration.from_yaml(
            example_configuration_dir / "tardis_configv1_verysimple.yml"
        )

        # Extract density config and create mesh
        density_config = config.model.structure.density
        velocity_config = config.model.structure.velocity
        time_explosion = config.supernova.time_explosion

        # Create solver from config
        solver = parse_density_solver_from_density_config(density_config)

        # Create mesh matching config
        velocity = np.linspace(
            velocity_config.start.to(u.km / u.s).value,
            velocity_config.stop.to(u.km / u.s).value,
            velocity_config.num + 1
        ) * u.km / u.s
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        # Solve for density
        result = solver.solve(mesh)

        # Verify it returns W7DensitySolver and produces valid results
        from tardis.io.model.classic.parse_density import W7DensitySolver
        
        assert isinstance(solver, W7DensitySolver)
        assert len(result.data) == velocity_config.num
        assert result.data.unit == u.g / u.cm**3
        # Check that density values are positive and in reasonable range
        assert np.all(result.data.value > 0)
        assert np.all(result.data.value < 1e-10)
        assert np.all(result.data.value > 1e-16)
