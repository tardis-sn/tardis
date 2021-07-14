.. _setup:

*************************
Setting Up the Simulation
*************************

The first step executed when TARDIS runs is to call an instance of the ``Simulation`` class. This sets up a lot of
things that TARDIS will need during its run. The two main things that are set up are the supernova model and the
initial plasma state (which may be updated throughout the simulation, see :ref:`est_and_conv`). The pages linked
below explain how these are calculated, as well as showing this in action by calling an instance of the
``Simulation`` class.

.. toctree::
    model
    plasma/index