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

The ``assemble_plasma`` function used by the simulation class is the one
defined in ``plasma/assembly``. It returns an assembled plasma by:

1. Defining a ``PlasmaSolverFactory`` instance.
2. Calling the ``prepare_factory`` method.
3. Calling the ``assemble`` method.

Plasma Properties
-----------------

The plasma properties used throughout the plasma module are defined in
``legacy_property_collections.py``. There is also a
``property_collections.py`` file which contains fewer classes. The plasma
built throughout the TARDIS codebase imports from
``legacy_property_collections.py``. The ``property_collections.py`` file
is not referenced elsewhere in the TARDIS codebase for plasma that gets
built during actual runs.

Property Categories
~~~~~~~~~~~~~~~~~~~

The following describes the property categories defined in
``legacy_property_collections.py`` and whether each class is built or
referenced outside its definition:

.. note::

   * **Green** — Covered in integration tests, with unit tests written
     where deemed useful.
   * **Yellow** — Referenced in the TARDIS codebase but never run in
     integration tests. These classes are not usually built during
     real TARDIS runs.
   * **Red** — Never referenced outside their class definition.
   * **Unhighlighted** — See inline notes for explanation.

Adiabatic Cooling Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``AdiabaticCoolingRate``

Basic Inputs
^^^^^^^^^^^^

.. note::

   For Basic Inputs, green indicates the input is used for classes built
   by plasma that TARDIS users actually use. Red indicates the input is
   usually never used. No unit tests are needed for Inputs.

* ``DilutePlanckianRadField``
* ``NumberDensity``
* ``TimeExplosion``
* ``AtomicData``
* ``JBlues``
* ``LinkTRadTElectron``
* ``HeliumTreatment``
* ``ContinuumInteractionSpecies`` — Never used as an input in any
  integration test. Many associated properties seem to be part of an
  unfinished refactor attempt.
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
* ``YgData`` — Referenced in
  ``plasma/equilibrium/tests/test_collisional_transitions.py`` and
  ``plasma/tests/test_plasma_continuum.py`` but not referenced in the
  plasma solver factory.

Continuum Interaction Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``BoundFreeOpacity``
* ``CollDeexcRateCoeff`` — Imported into ``collision_strengths.py``.
  Defined in two files: ``plasma/equilibrium/rates/collision_strengths.py``
  and ``plasma/properties/continuum_processes/rates.py``. The class
  defined in ``collision_strengths.py`` is fully uncovered.
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
* ``PhotoIonizationData`` — Referenced in the main plasma solver factory
  but never run.
* ``RawCollIonTransProbs``
* ``RawCollisionTransProbs``
* ``RawPhotoIonTransProbs``
* ``RawRadBoundBoundTransProbs``
* ``RawRecombTransProbs``
* ``SahaFactor``
* ``SpontRecombCoolingRateCoeff``
* ``SpontRecombRateCoeff`` — Only defined in the class and put into
  property collections but never referenced in the main plasma solver
  factory.
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

* ``HeliumNumericalNLTE`` — Requires a specific numerical NLTE solver and
  a specific atomic dataset, neither of which are distributed with TARDIS.

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
* ``NLTEPopulationSolverLU`` — Referenced in ``test_nlte_solver`` but
  never run in a full integration test.

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
* ``NLTEPopulationSolverRoot`` — Referenced in ``test_nlte_solver`` but
  never run in a full integration test.

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

In ``tardis/plasma/equilibrium`` there are files under a ``rates/`` folder
that represent another attempt at a plasma module refactor using
``RateSolver`` classes. These are not used in the actual implemented plasma
module but do have direct tests.

Uncovered Plasma Module Files
==============================

The following files in the plasma module are either fully uncovered or have
very little coverage, beyond what was already noted in the property
sections above.

``tardis/plasma/standard_plasmas.py``
   Zero percent coverage. Appears to be an attempt at refactoring how
   ``assemble_plasma`` is used by ``PlasmaSolverFactory``. Referenced
   only in ``grotarian_mockup.ipynb``.

``tardis/plasma/properties/property_collections.py``
   Property collection for a plasma module refactor that was never
   completed. Appears to be imported into
   ``tardis/plasma/standard_plasmas.py``.

``tardis/plasma/properties/nlte_excitation_data.py``
   Defines ``NLTEExcitationData``, which is not used anywhere else in
   the codebase.

``tardis/plasma/properties/hydrogen_continuum.py``
   Zero percent coverage by plasma tests.

``tardis/plasma/properties/continuum_processes/fast_array_util.py``
   Defines two numba functions:

   * ``numba_cumulative_trapezoid`` — Used in the uncovered class
     ``TwoPhotonEmissionCDF`` in
     ``tardis/plasma/properties/continuum_processes/rates.py``.
   * ``cumulative_integrate_array_by_blocks`` — Same as above.

Tools to Understand the Plasma Module
======================================

Test Coverage Reports
---------------------

Generating coverage reports of the plasma module test cases is a useful
way to understand what parts of the plasma module are actually used by
researchers running TARDIS. Running exclusively the tests in the plasma
module and examining coverage results can reveal, for example:

* ``plasma/standard_plasmas.py`` is fully uncovered. The
  ``assemble_plasma`` function defined there is never used in the actual
  TARDIS simulation. The ``assemble_plasma`` that TARDIS actually uses is
  in ``plasma/assembly/legacy_assembly.py``.
* The continuum interactions setup in ``PlasmaSolverFactory``
  (``plasma/assembly/base.py``) is never triggered because none of the
  tested plasma configs set up continuum interactions.

.. TODO: Add screenshot of PlasmaSolverFactory continuum interactions code here

* ``hydrogen_continuum.py`` is not used anywhere and has zero percent
  coverage.

The best way to find which classes are actually being called is to go
through ``legacy_property_collections.py``, click on each class, and
check whether the code inside the ``calculate`` function is covered.

For example, ``LinesLowerLevelIndex`` is in legacy properties and has
full coverage, meaning it is actually being calculated when
``test_full_plasmas`` builds its plasma configurations.

.. TODO: Add screenshot of LinesLowerLevelIndex coverage here

On the other hand, ``SpontRecombRateCoeff`` is in
``continuum_interactions_properties`` and its ``calculate`` method is
never run by any config in ``test_complete_plasmas``.

.. TODO: Add screenshot of SpontRecombRateCoeff coverage here

Going through coverage results for ``test_complete_plasma`` is a good
way to:

* Guide what is worth writing unit tests for.
* Identify features unused by the plasma module that should be migrated
  to Type II plasma exclusively.
* Find features that should be removed completely.
* Find features that need to be fleshed out further.

Plasma Graph Generation
-----------------------

As described in the :ref:`NetworkX - Plasma Graphs <plasma_module>`
section above, TARDIS can generate diagrams showing how plasma graph
properties connect to one another. Generating diagrams across multiple
configs is a good way to understand which plasma properties are relevant
to general TARDIS users.

Useful resources:

* TARDIS configs repository:
  `tardis-sn/tardis-configs <https://github.com/tardis-sn/tardis-configs>`_
* How to draw plasma graphs:
  :doc:`../../../how_to/output/how_to_plasma_graph`