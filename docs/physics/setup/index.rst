.. _setup:

*************************
Setting Up the Simulation
*************************

The first step executed when TARDIS runs is to call an instance of the ``Simulation`` class. This sets up a lot of
things that TARDIS will need during its run. The main things that are set up are the supernova model (a ``Radial1DModel`` object), the
initial plasma state (a ``BasePlasma`` object, which may be updated throughout the simulation, see :doc:`../update_and_conv/update_and_conv`),
and a ``MonteCarloRunner`` object. The pages linked below explain how the former two are calculated (the latter is used mainly in the
:doc:`next step of the calculation <../montecarlo/index>`), as well as showing this in action by calling an instance of the ``Simulation``
class.

.. toctree::
    setup_example
    model
    plasma/index