.. _plasma_module:

*********************
TARDIS Plasma Module
*********************

Plasma Structure
================

NetworkX - Plasma Graphs
-------------------------

The plasma is structured as a network of calculations. Each quantity
(temperature, ionization, optical depth, etc.) is a "property." Properties
depend on one another, and TARDIS connects them into a graph so that when
one value changes, only the downstream values get recalculated in the
correct order.

The plasma lives in the class ``BasePlasma``, defined in ``plasma/base.py``.

There are two types of properties:

* **Input properties** — values fed in from the model or config. They have
  outputs but no inputs, so they are the starting points of the graph.

* **Processing properties** — values that get calculated. TARDIS reads the
  argument names of their ``calculate`` function to learn what inputs
  they need.

For example, ``PhiSahaLTE`` (defined in
``plasma/properties/ion_population.py``) has the function::

   calculate(g_electron, beta_rad, partition_function, ionization_data)

Its inputs are those four names and its output is ``phi``. Connections are
made by matching names — an output called ``phi`` automatically links to
any calculation that has an argument called ``phi``.

Building the Graph
~~~~~~~~~~~~~~~~~~

TARDIS builds the graph inside ``_build_graph`` in three steps:

1. Add one node for every property.
2. Mark the input properties as starting points (they have no inputs).
3. For each calculated property, find what it needs, locate the property
   that produces each name, and draw an arrow from producer to consumer.
   If a property needs something that nothing produces, TARDIS raises an
   error immediately, catching broken or incomplete setups.

When an input changes (for example, the temperature), TARDIS finds
everything downstream of that change and recalculates in topological
order, so each property is updated only after the values it depends
on are fresh.

For a guide on how to display a plasma graph, see
:doc:`../../../how_to/output/how_to_plasma_graph`.

Plasma Solver Factory
---------------------

The ``PlasmaSolverFactory`` class, defined in ``plasma/assembly/base.py``,
initializes, configures, and assembles the plasma used in TARDIS runs. It
has two key methods:

* ``prepare_factory`` — accepts specific property collections, which are
  collections of classes that define plasma properties.
* ``assemble`` — returns an instance of the ``BasePlasma`` class.

There is also an ``assemble_plasma`` function defined in two locations:

* ``plasma/assembly/legacy_assembly.py``
* ``plasma/standard_plasmas.py``

The active ``assemble_plasma`` function used by the TARDIS simulation
class is the one defined in ``plasma/assembly/legacy_assembly.py``. It
returns an assembled plasma by:

1. Defining a ``PlasmaSolverFactory`` instance.
2. Calling the ``prepare_factory`` method.
3. Calling the ``assemble`` method.

Plasma Properties
-----------------

The plasma properties used throughout the plasma module are defined in
``legacy_property_collections.py``. A separate ``property_collections.py``
file exists but is not imported by the main plasma solver factory and is
therefore not active during standard TARDIS runs. The standard plasma solver
imports its property collections from ``legacy_property_collections.py``. 

Property Categories
~~~~~~~~~~~~~~~~~~~

The following describes the property categories defined in
``legacy_property_collections.py`` and the integration test coverage
status of each class:

.. note::

   * **Green** — Covered in integration tests, with unit tests written
     where deemed useful.
   * **Yellow** — Referenced in the TARDIS codebase but not exercised
     by integration tests. These classes are not built during standard
     TARDIS runs.
   * **Red** — Not referenced outside their class definition.
   * **Unhighlighted** — See inline notes for explanation.

Adiabatic Cooling Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``AdiabaticCoolingRate``

Basic Inputs
^^^^^^^^^^^^

.. note::

   For Basic Inputs, green indicates the input is used by plasma
   configurations that are active in standard TARDIS runs. Red indicates
   the input is not activated by any current simulation configuration.
   Input classes do not require unit tests.

* ``DilutePlanckianRadField``
* ``NumberDensity``
* ``TimeExplosion``
* ``AtomicData``
* ``JBlues``
* ``LinkTRadTElectron``
* ``HeliumTreatment``
* ``ContinuumInteractionSpecies`` — Not exercised by any current
  integration test. The associated continuum interaction properties
  are scoped to Type II supernova plasmas and are not activated in
  standard Type Ia simulation configs.
* ``NLTEIonizationSpecies``
* ``NLTEExcitationSpecies``

Basic Properties
^^^^^^^^^^^^^^^^

* ``BetaRadiation``
* ``DilutionFactor``
* ``ElectronTemperature``
* ``GElectron``
* ``IonizationData``
* ``Levels``
* ``Lines``
* ``LinesLowerLevelIndex``
* ``LinesUpperLevelIndex``
* ``PartitionFunction``
* ``SelectedAtoms``
* ``StimulatedEmission``
* ``TRadiative``
* ``TauSobolev``

Continuum Interaction Inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``PhotoIonRateCoeff``
* ``StimRecombRateFactor``
* ``BfHeatingRateCoeffEstimator``
* ``StimRecombCoolingRateCoeffEstimator``
* ``YgData`` — Tested directly in
  ``plasma/equilibrium/tests/test_collisional_transitions.py`` and
  ``plasma/tests/test_plasma_continuum.py``, but not wired into the
  main plasma solver factory and therefore not built during standard
  TARDIS runs.

Continuum Interaction Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``BoundFreeOpacity``
* ``CollDeexcRateCoeff`` — Defined in two files:
  ``plasma/equilibrium/rates/collision_strengths.py`` and
  ``plasma/properties/continuum_processes/rates.py``. The active class
  is the one imported into the property collection from ``rates.py``.
  The class in ``collision_strengths.py`` is not covered by the test
  suite.
* ``CollExcRateCoeff`` — Same as above.
* ``CollIonRateCoeffSeaton``
* ``CollRecombRateCoeff``
* ``ContinuumInteractionHandler``
* ``CorrPhotoIonRateCoeff``
* ``FreeBoundCoolingRate``
* ``FreeBoundEmissionCDF``
* ``FreeFreeCoolingRate``
* ``LevelIdxs2LineIdx``
* ``LevelIdxs2TransitionIdx``
* ``LevelNumberDensityLTE``
* ``MarkovChainIndex``
* ``MarkovChainTransProbs``
* ``MarkovChainTransProbsCollector``
* ``MonteCarloTransProbs``
* ``NonContinuumTransProbsMask``
* ``NonMarkovChainTransitionProbabilities``
* ``PhotoIonBoltzmannFactor``
* ``PhotoIonizationData`` — Registered in the plasma solver factory
  but not exercised by any current integration test configuration.
* ``RawCollIonTransProbs``
* ``RawCollisionTransProbs``
* ``RawPhotoIonTransProbs``
* ``RawRadBoundBoundTransProbs``
* ``RawRecombTransProbs``
* ``SahaFactor``
* ``SpontRecombCoolingRateCoeff``
* ``SpontRecombRateCoeff`` — Defined in the property collection but
  not referenced by the main plasma solver factory. Developers
  extending continuum interaction support should note this class
  requires explicit wiring to become active.
* ``StimRecombRateCoeff``
* ``ThermalGElectron``
* ``ThermalLTEPartitionFunction``
* ``ThermalLevelBoltzmannFactorLTE``
* ``ThermalPhiSahaLTE``
* ``YgInterpolator``

Dilute LTE Excitation Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``LevelBoltzmannFactorDiluteLTE``

Helium LTE Properties
^^^^^^^^^^^^^^^^^^^^^^

* ``LevelNumberDensity``
* ``IonNumberDensity``

Helium NLTE Properties
^^^^^^^^^^^^^^^^^^^^^^^

* ``HeliumNLTE``
* ``RadiationFieldCorrection``
* ``ZetaData``
* ``BetaElectron``
* ``LevelNumberDensityHeNLTE``
* ``IonNumberDensityHeNLTE``

Helium Numerical NLTE Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``HeliumNumericalNLTE`` — Requires a specific numerical NLTE solver
  and a specific atomic dataset. Neither are distributed with TARDIS.
  This class is not activated in standard simulation runs.

LTE Excitation Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``LevelBoltzmannFactorLTE``

LTE Ionization Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``PhiSahaLTE``

Nebular Ionization Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``PhiSahaNebular``
* ``ZetaData``
* ``BetaElectron``
* ``RadiationFieldCorrection``

NLTE LU Solver Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``NLTEIndexHelper``
* ``NLTEPopulationSolverLU`` — Tested in ``test_nlte_solver`` but not
  exercised in full integration tests. This solver is not activated
  in standard simulation configurations.

NLTE Properties
^^^^^^^^^^^^^^^^

* ``LevelBoltzmannFactorNLTE``
* ``NLTEData``
* ``PreviousElectronDensities``
* ``PreviousBetaSobolev``
* ``BetaSobolev``

NLTE Root Solver Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``NLTEIndexHelper``
* ``NLTEPopulationSolverRoot`` — Tested in ``test_nlte_solver`` but
  not exercised in full integration tests. This solver is not activated
  in standard simulation configurations.

Non NLTE Properties
^^^^^^^^^^^^^^^^^^^^

* ``LevelBoltzmannFactorNoNLTE``

Two Photon Properties
^^^^^^^^^^^^^^^^^^^^^^

* ``RawTwoPhotonTransProbs``
* ``TwoPhotonData``
* ``TwoPhotonEmissionCDF``
* ``TwoPhotonFrequencySampler``

Equilibrium
-----------

The ``tardis/plasma/equilibrium/rates/`` folder contains an alternative
plasma implementation based on ``RateSolver`` classes. This implementation
is independent of the NetworkX-based plasma module used in standard TARDIS
simulations and is not invoked during normal runs. The rate solver classes
have dedicated unit tests and are retained for future development.

Developer Notes: Inactive and Partially Covered Modules
========================================================

The following plasma module files are not invoked during standard TARDIS
simulation runs. Developers modifying the plasma module should be aware
of their status.

``tardis/plasma/standard_plasmas.py``
   Not covered by the plasma test suite and not called during standard
   TARDIS runs. This module is used exclusively by
   ``grotrian_mockup.ipynb``. The active plasma assembly entry point
   is ``assemble_plasma`` in ``plasma/assembly/legacy_assembly.py``.

``tardis/plasma/properties/property_collections.py``
   An alternative property collection used by ``standard_plasmas.py``.
   Not imported by the main plasma solver factory and therefore inactive
   during normal simulation runs. Developers should use
   ``legacy_property_collections.py`` as the authoritative reference
   for active plasma properties.

``tardis/plasma/properties/nlte_excitation_data.py``
   Defines ``NLTEExcitationData``, which is not referenced elsewhere in 
   the active codebase and is the existing implementation relevant to 
   future NLTE excitation work.

``tardis/plasma/properties/hydrogen_continuum.py``
   Not activated by any current simulation configuration and not
   covered by the plasma test suite.

``tardis/plasma/properties/continuum_processes/fast_array_util.py``
   Provides two Numba-accelerated integration utilities:

   * ``numba_cumulative_trapezoid`` — Used by ``TwoPhotonEmissionCDF``
     in ``tardis/plasma/properties/continuum_processes/rates.py``.
   * ``cumulative_integrate_array_by_blocks`` — Used by the same class.

   Both functions are not exercised by the plasma test suite, as
   ``TwoPhotonEmissionCDF`` is not activated in standard simulation
   configurations.

Developer Reference: Test Coverage and Plasma Graphs
=====================================================

Test Coverage
-------------

The plasma module test suite provides a practical way to identify which
properties are exercised by standard plasma configurations.
Running the plasma module tests with coverage enabled allows developers
to verify which property ``calculate`` methods are exercised by the
current set of simulation configurations.

Key architectural notes for developers:

* ``plasma/standard_plasmas.py`` is not covered by the plasma test
  suite. The active assembly entry point is
  ``plasma/assembly/legacy_assembly.py``.
* Continuum interaction properties in ``PlasmaSolverFactory``
  (``plasma/assembly/base.py``) are not activated by any current
  integration test configuration. Continuum interactions are not
  used in standard Type Ia supernova simulations.

.. TODO: Add screenshot of PlasmaSolverFactory continuum interactions code here

* ``hydrogen_continuum.py`` is not activated by any current simulation
  configuration and is not covered by the plasma test suite.

To determine whether a specific property class is active, inspect its
``calculate`` method under coverage. For example:

* ``LinesLowerLevelIndex`` has full coverage, confirming it is
  calculated in every standard plasma configuration built by
  ``test_complete_plasmas``.

.. TODO: Add screenshot of LinesLowerLevelIndex coverage here

* ``SpontRecombRateCoeff``, defined in
  ``continuum_interactions_properties``, has no coverage in
  ``test_complete_plasmas``, confirming it is not built by any
  current standard simulation configuration.

.. TODO: Add screenshot of SpontRecombRateCoeff coverage here

Coverage analysis of ``test_complete_plasmas`` can be used to:

* Identify which property classes require additional unit tests.
* Determine which features are scoped to Type II plasma and should
  not be assumed available in Type Ia runs.
* Identify property classes that are candidates for removal or
  consolidation.

Plasma Graph Generation
-----------------------

TARDIS can generate diagrams of the plasma property graph, showing
which properties depend on one another. Generating graphs across
multiple simulation configurations is a practical method for
identifying which plasma properties are relevant to standard
simulation runs.

Useful resources:

* TARDIS configs repository:
  `tardis-sn/tardis-configs <https://github.com/tardis-sn/tardis-configs>`_
* How to draw plasma graphs:
  :doc:`../../../how_to/output/how_to_plasma_graph`