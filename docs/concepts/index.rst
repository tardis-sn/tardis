.. _concepts:

******************
TARDIS Concepts
******************

This section explains the fundamental concepts and architectural principles that underlie TARDIS.

.. toctree::
   :maxdepth: 2

   modular_architecture

Core Architectural Concepts
============================

**Modular Architecture**
   TARDIS separates topology, physics, and simulation into distinct, composable layers.
   This enables modularity, testability, and clear separation of concerns.
   See :doc:`modular_architecture` for the complete guide.

**Separation of Concerns**
   TARDIS separates topology (mesh structure) from physics (density, temperature)
   from simulation (Monte Carlo transport, plasma calculations).

**Type Safety**
   Using ``MeshQuantity`` and ``MeshLocation`` enums to ensure physical quantities
   are bound to the correct mesh locations.
