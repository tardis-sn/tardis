.. _model_io:

*******************
TARDIS Model Input
*******************

This section covers how TARDIS constructs and manages supernova models.

.. contents:: Contents
   :local:
   :depth: 2

Overview
========

TARDIS uses a **mesh-first architecture** where the simulation space is defined first,
and physical quantities are bound to specific locations on that mesh. This provides
a clean separation between topology (mesh structure) and physics (density, temperature, etc.).

Architecture Guide
==================

See :doc:`/concepts/modular_architecture` for the complete architectural guide explaining
the philosophical foundations and design principles of TARDIS's modular architecture.
**Start there** to understand the conceptual framework before diving into practical usage.

Practical Guide
===============

See :doc:`/how_to/build_model` for a hands-on guide to building TARDIS models from components,
with code examples and common patterns.

Practical Guides
================

.. toctree::
   :maxdepth: 1

   geometry_mesh_classes

The :doc:`geometry_mesh_classes` notebook provides hands-on examples of creating
and manipulating meshes in TARDIS.

Additional Examples
===================

.. toctree::
   :maxdepth: 1

   parse_snec_to_tardis

The :doc:`parse_snec_to_tardis` notebook demonstrates how to import data from
hydrodynamic simulations into TARDIS's mesh framework.

Key Concepts
============

Mesh Types
----------

**Radial1DMesh**
   Spatial mesh with radii in physical units (e.g., cm)

**HomologousRadial1DMesh**
   Velocity-based mesh with homologous expansion (:math:`r = v \times t`)

Mesh Locations
--------------

**INTERFACE**
   Values at cell boundaries (N+1 values for N cells)
   
   Example: radii, interface velocities

**VOLUME**
   Values at cell centers (N values for N cells)
   
   Example: density, temperature, abundances

See Also
========

- :doc:`/io/configuration/index` - Configuration file format
- :doc:`/physics_walkthrough/setup/index` - Physics setup walkthrough
- API: :mod:`tardis.model.mesh`
