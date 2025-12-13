.. _build_model:

**************************************
Building TARDIS Models from Components
**************************************

This guide shows you how to build TARDIS supernova models using the modular architecture. You'll learn to create meshes, calculate physical quantities, and assemble complete models.

.. contents:: Contents
   :local:
   :depth: 2

Prerequisites
=============

You should understand:

- Basic TARDIS concepts (supernova ejecta, 1D models)
- Python programming basics
- The :doc:`/concepts/modular_architecture` (recommended)

Quick Start
===========

Three-Step Model Building
--------------------------

Building a TARDIS model follows three steps:

.. code-block:: python

    from tardis.model.mesh import HomologousRadial1DMesh
    from tardis.io.model.classic.parse_density import UniformDensitySolver
    from astropy import units as u
    import numpy as np

    # Step 1: Create mesh (define spatial structure)
    velocities = np.linspace(10000, 20000, 21) * u.km/u.s
    mesh = HomologousRadial1DMesh.from_velocity_interfaces(
        velocity=velocities,
        time_explosion=10 * u.day
    )

    # Step 2: Calculate density (add physics)
    solver = UniformDensitySolver(density_value=1e-14 * u.g/u.cm**3)
    density = solver.solve(mesh)

    # Step 3: Add more physics and create simulation
    # (abundances, temperature, etc.)

Step 1: Creating Meshes
========================

From Velocity Structure (Homologous Expansion)
-----------------------------------------------

Most supernova models use homologous expansion where :math:`r = v \times t`:

.. code-block:: python

    from tardis.model.mesh import HomologousRadial1DMesh
    from astropy import units as u
    import numpy as np

    # Define velocity boundaries (N+1 interfaces for N cells)
    N_SHELLS = 20
    V_INNER = 8000  # km/s
    V_OUTER = 20000  # km/s
    
    velocity_interfaces = np.linspace(
        V_INNER, V_OUTER, N_SHELLS + 1
    ) * u.km/u.s
    
    # Create mesh
    mesh = HomologousRadial1DMesh.from_velocity_interfaces(
        velocity=velocity_interfaces,
        time_explosion=10 * u.day
    )
    
    print(f"Created mesh with {mesh.n_cells} cells")
    print(f"Velocity range: {mesh.velocity.data[0]} to {mesh.velocity.data[-1]}")

**Key Points**:
   - N+1 velocity interfaces define N cells
   - Velocities are at INTERFACE locations
   - Time explosion is required for :math:`r = v \times t`

From Spatial Structure (Radial Mesh)
-------------------------------------

If you have radii directly (e.g., from hydro simulations):

.. code-block:: python

    from tardis.model.mesh import Radial1DMesh
    from astropy import units as u

    # Define radii at cell interfaces
    radii = [1e8, 2e8, 3e8, 4e8, 5e8] * u.cm
    
    mesh = Radial1DMesh.from_cell_interfaces(radii)
    
    print(f"Cells: {mesh.n_cells}")
    print(f"Radii: {mesh.radius.data}")

Step 2: Calculating Density
============================

Using Density Solvers
---------------------

TARDIS provides several density solvers:

**Uniform Density**

.. code-block:: python

    from tardis.io.model.classic.parse_density import UniformDensitySolver
    
    solver = UniformDensitySolver(density_value=1e-14 * u.g/u.cm**3)
    density = solver.solve(mesh)
    
    # density is a MeshQuantity with N values (one per cell)
    print(f"Density shape: {density.data.shape}")
    print(f"Defined at: {density.defined_at}")  # VOLUME

**Power Law Density**

.. code-block:: python

    from tardis.io.model.classic.parse_density import PowerLawDensitySolver
    
    solver = PowerLawDensitySolver(
        density_0=3e-13 * u.g/u.cm**3,
        velocity_0=10000 * u.km/u.s,
        exponent=2  # rho ~ v^(-2)
    )
    density = solver.solve(mesh)

**W7 Type Ia Density (Branch85)**

.. code-block:: python

    from tardis.io.model.classic.parse_density import W7DensitySolver
    
    solver = W7DensitySolver(
        time_0=0.000231481 * u.day,
        density_0=3e29 * u.g/u.cm**3,
        velocity_0=1 * u.km/u.s
    )
    density = solver.solve(mesh)

Direct Assignment
-----------------

For custom density profiles:

.. code-block:: python

    from tardis.model.mesh import MeshQuantity, MeshLocation
    
    # Create custom density profile (N values for N cells)
    density_data = np.array([1e-13, 8e-14, 5e-14, 3e-14]) * u.g/u.cm**3
    
    density = MeshQuantity(
        data=density_data,
        defined_at=MeshLocation.VOLUME,
        name="density"
    )
    
    # Verify it matches mesh
    assert len(density.data) == mesh.n_cells

Step 3: Working with Configurations
====================================

From YAML Config Files
----------------------

TARDIS configs work with the modular architecture:

.. code-block:: python

    from tardis.io.configuration.config_reader import Configuration
    from tardis.io.model.classic.parse_density import (
        parse_density_solver_from_density_config
    )
    
    # Load configuration
    config = Configuration.from_yaml("my_config.yml")
    
    # Extract mesh parameters
    v_config = config.model.structure.velocity
    velocities = np.linspace(
        v_config.start.to(u.km/u.s).value,
        v_config.stop.to(u.km/u.s).value,
        v_config.num + 1
    ) * u.km/u.s
    
    # Create mesh
    mesh = HomologousRadial1DMesh.from_velocity_interfaces(
        velocity=velocities,
        time_explosion=config.supernova.time_explosion
    )
    
    # Parse density configuration and solve
    density_config = config.model.structure.density
    solver = parse_density_solver_from_density_config(density_config)
    density = solver.solve(mesh)

Common Patterns
===============

Pattern 1: Creating Test Meshes
--------------------------------

For unit tests:

.. code-block:: python

    import pytest
    
    @pytest.fixture
    def simple_mesh():
        """Simple test mesh with 5 cells."""
        v = np.array([1000, 2000, 3000, 4000, 5000, 6000]) * u.km/u.s
        return HomologousRadial1DMesh.from_velocity_interfaces(
            velocity=v,
            time_explosion=1 * u.day
        )
    
    def test_density_solver(simple_mesh):
        solver = UniformDensitySolver(density_value=1e-14 * u.g/u.cm**3)
        result = solver.solve(simple_mesh)
        
        assert len(result.data) == simple_mesh.n_cells
        assert result.defined_at == MeshLocation.VOLUME

Pattern 2: Verifying Consistency
---------------------------------

Always check that quantities match mesh structure:

.. code-block:: python

    def verify_mesh_quantity(quantity, mesh, expected_location):
        """Verify MeshQuantity is consistent with mesh."""
        if expected_location == MeshLocation.VOLUME:
            assert len(quantity.data) == mesh.n_cells
        elif expected_location == MeshLocation.INTERFACE:
            assert len(quantity.data) == mesh.n_cells + 1
        
        assert quantity.defined_at == expected_location
    
    # Usage
    verify_mesh_quantity(density, mesh, MeshLocation.VOLUME)

Pattern 3: Converting Between Mesh Types
-----------------------------------------

Convert homologous to spatial mesh:

.. code-block:: python

    # Create homologous mesh
    hom_mesh = HomologousRadial1DMesh.from_velocity_interfaces(
        velocity=velocities,
        time_explosion=time_exp
    )
    
    # Convert to spatial mesh
    spatial_mesh, velocity_field = hom_mesh.to_spatial_mesh()
    
    # spatial_mesh is a Radial1DMesh
    # velocity_field preserves the velocity data
    print(f"Radii: {spatial_mesh.radius.data}")
    print(f"Velocities preserved: {velocity_field.data}")

Complete Example: Type Ia Model
================================

Here's a complete example building a W7-like Type Ia model:

.. code-block:: python

    from tardis.model.mesh import HomologousRadial1DMesh, MeshLocation
    from tardis.io.model.classic.parse_density import W7DensitySolver
    from astropy import units as u
    import numpy as np

    # 1. Define model parameters
    N_SHELLS = 20
    V_INNER = 11000  # km/s
    V_OUTER = 20000  # km/s
    T_EXPLOSION = 13  # days

    # 2. Create velocity mesh
    velocity_interfaces = np.linspace(
        V_INNER, V_OUTER, N_SHELLS + 1
    ) * u.km/u.s
    
    mesh = HomologousRadial1DMesh.from_velocity_interfaces(
        velocity=velocity_interfaces,
        time_explosion=T_EXPLOSION * u.day
    )

    # 3. Calculate W7 density profile
    w7_solver = W7DensitySolver(
        time_0=0.000231481 * u.day,
        density_0=3e29 * u.g/u.cm**3,
        velocity_0=1 * u.km/u.s
    )
    
    density = w7_solver.solve(mesh)

    # 4. Verify results
    print(f"Model created with {mesh.n_cells} cells")
    print(f"Density range: {density.data.min()} to {density.data.max()}")
    print(f"First cell density: {density.data[0]}")
    print(f"Last cell density: {density.data[-1]}")
    
    # Density should decrease with radius for W7
    assert density.data[0] < density.data[-1]

Troubleshooting
===============

Common Errors
-------------

**"MeshQuantity has wrong length"**

.. code-block:: python

    # Wrong: N+1 values for volume quantity
    density_data = np.ones(mesh.n_cells + 1) * u.g/u.cm**3  # ❌
    
    # Correct: N values for volume quantity
    density_data = np.ones(mesh.n_cells) * u.g/u.cm**3  # ✓

**"Units not compatible"**

.. code-block:: python

    # Always include units
    velocities = np.array([1000, 2000, 3000])  # ❌ No units
    velocities = np.array([1000, 2000, 3000]) * u.km/u.s  # ✓

**"MeshLocation mismatch"**

.. code-block:: python

    # Check where your quantity should live
    # Density is always at VOLUME, not INTERFACE
    density = MeshQuantity(
        data=rho_data,
        defined_at=MeshLocation.INTERFACE,  # ❌
        name="density"
    )

Next Steps
==========

**Learn More**:
   - :doc:`/io/model/geometry_mesh_classes` - Interactive notebook examples
   - :doc:`/concepts/modular_architecture` - Architectural concepts
   - :doc:`/io/configuration/index` - Configuration file format

**API Documentation**:
   - :class:`~tardis.model.mesh.HomologousRadial1DMesh`
   - :class:`~tardis.model.mesh.Radial1DMesh`
   - :class:`~tardis.model.mesh.MeshQuantity`
   - :mod:`tardis.io.model.classic.parse_density` - Density solvers
