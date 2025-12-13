.. _modular_architecture:

*****************************
TARDIS Modular Architecture
*****************************

TARDIS uses a modular architecture that separates topology, physics, and simulation into distinct, composable layers. This guide explains the conceptual foundations and design philosophy.

.. contents:: Contents
   :local:
   :depth: 2

Core Architectural Principles
==============================

Separation of Concerns
----------------------

TARDIS separates three distinct aspects of supernova modeling:

**Topology (Spatial Structure)**
   The mesh defines the geometric structure without any physics attached.
   
   - Cell boundaries (interfaces)
   - Cell centers (volumes)
   - Coordinate systems (radial, velocity-based)

**Physics (Physical Quantities)**
   Physical data bound to specific mesh locations via ``MeshQuantity``.
   
   - Density, temperature, abundances
   - Each quantity knows where it lives (interface or volume)
   - Type-safe binding prevents location errors

**Simulation (High-Level Components)**
   Simulation components use mesh + physics together.
   
   - Monte Carlo transport
   - Plasma calculations
   - Convergence strategies

The Three-Layer Architecture
-----------------------------

.. code-block:: none

    ┌─────────────────────────────────────────────────────┐
    │                 SIMULATION LAYER                    │
    │  (SimulationState, Transport, Plasma)              │
    └─────────────────────────────────────────────────────┘
                            ▲
                            │ uses
                            │
    ┌─────────────────────────────────────────────────────┐
    │                  PHYSICS LAYER                      │
    │  (MeshQuantity: density, temperature, abundances)   │
    └─────────────────────────────────────────────────────┘
                            ▲
                            │ bound to
                            │
    ┌─────────────────────────────────────────────────────┐
    │                 TOPOLOGY LAYER                      │
    │  (Radial1DMesh, HomologousRadial1DMesh)            │
    └─────────────────────────────────────────────────────┘

Each layer can be developed, tested, and modified independently.

Why This Architecture?
=======================

Benefits
--------

**Modularity**
   Components have single, well-defined responsibilities. Change one layer without affecting others.

**Testability**
   Test topology, physics, and simulation independently. No need for full simulation setup to test mesh creation or density solvers.

**Physical Clarity**
   Matches how physicists think: define spatial structure, specify conditions, solve equations.

**Flexibility**
   Easy to extend with new mesh types, physical models, or simulation strategies.

**Type Safety**
   ``MeshLocation`` enum prevents accidentally using interface values where volume values are expected (and vice versa).

Comparison with Legacy Approach
--------------------------------

**Legacy (Configuration-Driven)**:
   Everything mixed in YAML - spatial structure, physics parameters, solver configs all together. Implicit dependencies, hard to test components in isolation.

**Modular Architecture**:
   Explicit separation. Create mesh → calculate physics → build simulation. Each step is clear, testable, and composable.

Topology vs. Physics
=====================

Understanding MeshLocation
--------------------------

Physical quantities exist at specific locations on the mesh:

**INTERFACE**
   Values at cell boundaries (N+1 values for N cells)
   
   .. code-block:: none

       Interface:  |-------|-------|-------|
                   0       1       2       3
       
       4 interfaces define 3 cells
   
   *Examples*: radii, interface velocities

**VOLUME**
   Values at cell centers (N values for N cells)
   
   .. code-block:: none

       Volume:     |   0   |   1   |   2   |
       Interface:  0       1       2       3
       
       3 volumes between 4 interfaces
   
   *Examples*: density, temperature, abundances

Why Location Matters
--------------------

**Physical Correctness**
   Different quantities naturally live at different locations:
   
   - Boundary conditions (velocity) → interfaces
   - Bulk properties (density) → volumes
   - Fluxes → interfaces

**Numerical Stability**
   Proper location prevents interpolation errors and maintains conservation laws.

**Explicit Semantics**
   Code clearly states where data lives:

   .. code-block:: python

       density = MeshQuantity(
           data=rho_data,
           defined_at=MeshLocation.VOLUME,
           name="density"
       )

Mesh Types
==========

TARDIS provides mesh types for different use cases:

Radial1DMesh
------------

**Purpose**: Direct spatial mesh with physical radii

**Use Case**: When spatial structure comes from external sources (e.g., hydrodynamic simulations)

**Features**:
   - N+1 interface radii define N cells
   - No time dependence
   - Pure spatial structure

HomologousRadial1DMesh
----------------------

**Purpose**: Velocity-based mesh with homologous expansion (:math:`r = v \times t`)

**Use Case**: Standard supernova models

**Features**:
   - Define using velocity interfaces
   - Requires explosion time
   - Converts to spatial mesh when needed
   - Preserves velocity field

Design Decisions
================

Why Explicit MeshLocation?
---------------------------

**Alternative**: Infer location from array shape

**Why Explicit**:
   - No ambiguity
   - Self-documenting code
   - Catches errors at creation time
   - Easy to extend (EDGE, CORNER, etc.)

Why Immutable MeshQuantity?
----------------------------

**Design**: ``MeshQuantity`` is a frozen dataclass

**Rationale**:
   - Prevents accidental modification
   - Thread-safe
   - Clear separation: "calculate" not "modify"

**Pattern**: Create new objects rather than modifying existing ones

Why Not Unstructured Meshes?
-----------------------------

**Current**: Radial 1D with regular structure

**Rationale**:
   - Matches TARDIS's 1D spherical symmetry
   - Simpler and more efficient
   - Can extend to 2D/3D later without breaking API

Extensibility
=============

Future Directions
-----------------

The architecture supports future extensions:

**Multi-Dimensional Meshes**
   Architecture allows 2D/3D meshes without breaking existing code.

**Adaptive Meshes**
   Can add adaptive refinement while keeping same ``MeshQuantity`` interface.

**New Physical Models**
   New solvers plug into existing mesh framework.

**Alternative Geometries**
   Cartesian, cylindrical, etc. can coexist with radial meshes.

Migration Path
--------------

Existing YAML configs still work - they're parsed to create mesh and physics objects explicitly.

New code can use modular approach directly:

.. code-block:: python

    # Create mesh
    mesh = HomologousRadial1DMesh.from_velocity_interfaces(...)
    
    # Calculate physics
    solver = UniformDensitySolver(density_value=...)
    density = solver.solve(mesh)
    
    # Build simulation
    sim = Simulation.from_mesh_and_quantities(mesh, density, ...)

See Also
========

**How-To Guides**:
   - :doc:`/how_to/build_model` - Practical model building
   - :doc:`/io/model/geometry_mesh_classes` - Mesh usage examples

**API Reference**:
   - :mod:`tardis.model.mesh` - Mesh classes
   - :class:`~tardis.model.mesh.MeshQuantity`
   - :class:`~tardis.model.mesh.MeshLocation`
