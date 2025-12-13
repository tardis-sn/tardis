"""Tests for density solver classes."""

from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.classic.parse_density import (
    ExponentialDensitySolver,
    PowerLawDensitySolver,
    UniformDensitySolver,
    W7DensitySolver,
    parse_density_solver_from_density_config,
)
from tardis.model import SimulationState
from tardis.model.mesh import HomologousRadial1DMesh


@pytest.fixture(scope="module")
def example_configuration_dir():
    """Path to example configuration files."""
    return Path("tardis/io/configuration/tests/data")


@pytest.fixture
def simple_mesh():
    """Create a simple test mesh."""
    velocity_interfaces = np.array([1000, 2000, 3000, 4000, 5000]) * u.km / u.s
    time_explosion = 1.0 * u.day
    return HomologousRadial1DMesh.from_velocity_interfaces(
        velocity=velocity_interfaces, time_explosion=time_explosion
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


class TestDensityFromConfig:
    """Test density solvers using actual YAML configuration files."""

    @pytest.fixture
    def example_configuration_dir(self):
        """Path to example configuration directory."""
        return Path(__file__).parent.parent.parent.parent / "io" / "configuration" / "tests" / "data"

    @pytest.fixture
    def atomic_dataset(self):
        """Load atomic dataset for tests."""
        atom_data_path = pytest.config.getoption("--atomic-dataset")
        if atom_data_path is not None:
            return AtomData.from_hdf(atom_data_path)
        pytest.skip("--atomic-dataset not specified")

    def test_w7_density_from_config(self, example_configuration_dir, atomic_dataset):
        """Test W7 (branch85_w7) density profile from configuration file."""
        config = Configuration.from_yaml(
            example_configuration_dir / "paper1_tardis_configv1.yml"
        )

        simulation_state = SimulationState.from_config(
            config, atom_data=atomic_dataset
        )

        # Verify density was calculated
        assert len(simulation_state.density) == 20
        assert simulation_state.density.unit == u.g / u.cm**3

        # Check expected density values for W7 profile
        # First cell should have higher density than last
        assert simulation_state.density[0] > simulation_state.density[-1]

        # Check against known values from paper1 config
        assert_allclose(
            simulation_state.density[0].cgs.value,
            7.542803599143591e-14,
            rtol=1e-5,
        )

    def test_uniform_density_from_config(self, example_configuration_dir, atomic_dataset):
        """Test uniform density profile from configuration file."""
        config = Configuration.from_yaml(
            example_configuration_dir / "tardis_configv1_uniform_density.yml"
        )

        simulation_state = SimulationState.from_config(
            config, atom_data=atomic_dataset
        )

        # Verify uniform density
        assert len(simulation_state.density) == 20
        expected_density = 1e-14 * u.g / u.cm**3
        assert_allclose(
            simulation_state.density.value,
            np.ones(20) * expected_density.value,
            rtol=1e-10,
        )

    def test_power_law_density_from_config(self, example_configuration_dir):
        """Test power law density profile using solver directly with config-like object."""
        # Simulate config values
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 10000 * u.km / u.s
        exponent = 3.0

        solver = PowerLawDensitySolver(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )

        # Create test mesh
        velocity = np.linspace(10000, 20000, 21) * u.km / u.s
        time_explosion = 13.0 * u.day
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        result = solver.solve(mesh)

        # Verify power law behavior
        assert len(result.data) == 20
        # Density should decrease monotonically for positive exponent
        assert np.all(np.diff(result.data.value) < 0)

        # Check first value matches expected power law
        v_mid_0 = 0.5 * (velocity[0] + velocity[1])
        expected_0 = density_0 * (v_mid_0 / velocity_0) ** (-exponent)
        assert_allclose(result.data[0].value, expected_0.value, rtol=1e-10)

    def test_exponential_density_from_config(self, example_configuration_dir):
        """Test exponential density profile using solver directly with config-like object."""
        # Simulate config values
        density_0 = 1e-13 * u.g / u.cm**3
        velocity_0 = 5000 * u.km / u.s

        solver = ExponentialDensitySolver(
            density_0=density_0, velocity_0=velocity_0
        )

        # Create test mesh
        velocity = np.linspace(10000, 20000, 21) * u.km / u.s
        time_explosion = 13.0 * u.day
        mesh = HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=velocity, time_explosion=time_explosion
        )

        result = solver.solve(mesh)

        # Verify exponential behavior
        assert len(result.data) == 20
        # Density should decrease monotonically
        assert np.all(np.diff(result.data.value) < 0)

        # Check first value matches expected exponential
        v_mid_0 = 0.5 * (velocity[0] + velocity[1])
        expected_0 = density_0 * np.exp(-v_mid_0 / velocity_0)
        assert_allclose(result.data[0].value, expected_0.value, rtol=1e-10)
