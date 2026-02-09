"""Density profile solver classes for TARDIS classic models."""

from dataclasses import dataclass

import numpy as np
from astropy import units as u

from tardis.io.model.parse_density_configuration import (
    calculate_density_after_time,
)
from tardis.model.mesh import HomologousRadial1DMesh, MeshLocation, MeshQuantity


@dataclass
class UniformDensitySolver:
    """Uniform density profile solver.

    Parameters
    ----------
    density_value : u.Quantity
        The uniform density value.
    """

    density_value: u.Quantity

    @classmethod
    def from_config(cls, config):
        """Create UniformDensitySolver from a TARDIS config.

        Parameters
        ----------
        config : object
            Configuration object with model.structure.density section.

        Returns
        -------
        UniformDensitySolver
            The uniform density profile solver.
        """
        density_config = config.model.structure.density
        density_value = density_config.value * u.Unit(density_config.unit)
        return cls(density_value=density_value)

    def solve(self, mesh: HomologousRadial1DMesh) -> MeshQuantity:
        """Calculate density for the given mesh.

        Parameters
        ----------
        mesh : HomologousRadial1DMesh
            The mesh to calculate density for.

        Returns
        -------
        MeshQuantity
            Density at cell volumes.
        """
        n_cells = len(mesh.velocity.data) - 1
        density_data = np.ones(n_cells) * self.density_value

        return MeshQuantity(
            data=density_data, defined_at=MeshLocation.VOLUME, name="density"
        )


@dataclass
class PowerLawDensitySolver:
    """Power law density profile solver: rho(v) = rho_0 * (v / velocity_0)^(-n).

    Parameters
    ----------
    density_0 : u.Quantity
        Reference density at reference velocity.
    velocity_0 : u.Quantity
        Reference velocity.
    exponent : float
        Power law exponent (positive).
    """

    density_0: u.Quantity
    velocity_0: u.Quantity
    exponent: float

    @classmethod
    def from_config(cls, config):
        """Create PowerLawDensitySolver from a TARDIS config.

        Parameters
        ----------
        config : object
            Configuration object with model.structure.density section.

        Returns
        -------
        PowerLawDensitySolver
            The power law density profile solver.
        """
        density_config = config.model.structure.density
        density_0 = density_config.rho_0 * u.Unit(density_config.unit)
        velocity_0 = density_config.v_0 * u.Unit(density_config.v_unit)
        exponent = density_config.exponent
        return cls(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )

    def solve(self, mesh: HomologousRadial1DMesh) -> MeshQuantity:
        """Calculate density for the given mesh.

        Parameters
        ----------
        mesh : HomologousRadial1DMesh
            The mesh to calculate density for.

        Returns
        -------
        MeshQuantity
            Density at cell volumes.
        """
        velocity_middle = 0.5 * (mesh.velocity.data[:-1] + mesh.velocity.data[1:])
        density_data = self.density_0 * (velocity_middle / self.velocity_0) ** (
            -self.exponent
        )

        return MeshQuantity(
            data=density_data, defined_at=MeshLocation.VOLUME, name="density"
        )


@dataclass
class ExponentialDensitySolver:
    """Exponential density profile solver: rho(v) = density_0 * exp(-v / velocity_0).

    Parameters
    ----------
    density_0 : u.Quantity
        Reference density.
    velocity_0 : u.Quantity
        Characteristic velocity scale.
    """

    density_0: u.Quantity
    velocity_0: u.Quantity

    @classmethod
    def from_config(cls, config):
        """Create ExponentialDensitySolver from a TARDIS config.

        Parameters
        ----------
        config : object
            Configuration object with model.structure.density section.

        Returns
        -------
        ExponentialDensitySolver
            The exponential density profile solver.
        """
        density_config = config.model.structure.density
        density_0 = density_config.rho_0 * u.Unit(density_config.unit)
        velocity_0 = density_config.v_0 * u.Unit(density_config.v_unit)
        return cls(density_0=density_0, velocity_0=velocity_0)

    def solve(self, mesh: HomologousRadial1DMesh) -> MeshQuantity:
        """Calculate density for the given mesh.

        Parameters
        ----------
        mesh : HomologousRadial1DMesh
            The mesh to calculate density for.

        Returns
        -------
        MeshQuantity
            Density at cell volumes.
        """
        v_mid = 0.5 * (mesh.velocity.data[:-1] + mesh.velocity.data[1:])
        density_data = self.density_0 * np.exp(
            -(v_mid / self.velocity_0).decompose()
        )

        return MeshQuantity(
            data=density_data, defined_at=MeshLocation.VOLUME, name="density"
        )


@dataclass
class W7DensitySolver:
    """W7 Type Ia supernova density profile solver (Branch85 W7).

    This implements the Branch85 W7 model using a power law density profile
    with exponent -7, scaled by the time of explosion.

    Parameters
    ----------
    time_0 : u.Quantity
        Reference time for the W7 model (default: 0.000231481 day).
    density_0 : u.Quantity
        Reference density normalization (default: 3e29 g/cm^3).
    velocity_0 : u.Quantity
        Reference velocity (default: 1 km/s).
    """

    time_0: u.Quantity
    density_0: u.Quantity
    velocity_0: u.Quantity

    @classmethod
    def from_config(cls, config):
        """Create W7DensitySolver from a TARDIS config.

        Parameters
        ----------
        config : object
            Configuration object with model.structure.density section.

        Returns
        -------
        W7DensitySolver
            The W7 density profile solver.
        """
        density_config = config.model.structure.density
        time_0 = density_config.time_0 * u.Unit(density_config.time_unit)
        density_0 = density_config.rho_0 * u.Unit(density_config.unit)
        velocity_0 = density_config.v_0 * u.Unit(density_config.v_unit)
        return cls(time_0=time_0, density_0=density_0, velocity_0=velocity_0)

    def solve(self, mesh: HomologousRadial1DMesh) -> MeshQuantity:
        """Calculate density for the given mesh.

        Parameters
        ----------
        mesh : HomologousRadial1DMesh
            The mesh to calculate density for.

        Returns
        -------
        MeshQuantity
            Density at cell volumes.
        """
        velocity_middle = 0.5 * (mesh.velocity.data[:-1] + mesh.velocity.data[1:])
        
        # Calculate power law density: rho = density_0 * (v / v_0)^(-7)
        density_at_time_0 = self.density_0 * (velocity_middle / self.velocity_0) ** (-7)
        
        # Scale density from time_0 to time_explosion using calculate_density_after_time
        density_data = calculate_density_after_time(
            density_at_time_0, self.time_0, mesh.time_explosion
        )

        return MeshQuantity(
            data=density_data, defined_at=MeshLocation.VOLUME, name="density"
        )


def parse_density_solver_from_density_config(density_config):
    """Create a density solver from density configuration section.

    This is a transition layer that works with the density configuration object
    directly (as used in parse_density_configuration.py).

    Parameters
    ----------
    density_config : ConfigurationNameSpace
        Density configuration section (e.g., config.model.structure.density).

    Returns
    -------
    solver : UniformDensitySolver, PowerLawDensitySolver, ExponentialDensitySolver, or W7DensitySolver
        The appropriate density solver based on configuration type.

    Raises
    ------
    ValueError
        If density type is not recognized.
    """
    density_type = density_config.type

    if density_type == "uniform":
        density_value = density_config.value
        return UniformDensitySolver(density_value=density_value)
    if density_type == "power_law":
        density_0 = density_config.rho_0
        velocity_0 = density_config.v_0
        exponent = density_config.exponent
        return PowerLawDensitySolver(
            density_0=density_0, velocity_0=velocity_0, exponent=exponent
        )
    if density_type == "exponential":
        density_0 = density_config.rho_0
        velocity_0 = density_config.v_0
        return ExponentialDensitySolver(density_0=density_0, velocity_0=velocity_0)
    if density_type == "branch85_w7":
        time_0 = density_config.w7_time_0
        velocity_0 = density_config.w7_v_0
        rho_0 = density_config.w7_rho_0
        return W7DensitySolver(
            time_0=time_0, density_0=rho_0, velocity_0=velocity_0
        )
    raise ValueError(f"Unrecognized density type '{density_type}'")
