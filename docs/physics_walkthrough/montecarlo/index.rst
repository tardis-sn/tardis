.. _montecarlo:

*********************
Monte Carlo Iteration
*********************

After setting up the simulation, TARDIS runs the simulation using the ``.run_convergence()`` method. This runs several Monte
Carlo iterations (which will be described in the links below), corresponding to one less than the number of iterations specified
in the :ref:`Monte Carlo Configuration <montecarlo-config>`. As will be described in :doc:`../update_and_conv/update_and_conv` and
:ref:`spectrum`, these iterations will eventually be used to calculate the steady-state plasma properties.
TARDIS runs the last iteration of the simulation using the ``.run_final()`` method to determine the spectrum.

The following pages provide a very basic introduction to Monte Carlo radiative
transfer techniques as they are used in TARDIS. All the information listed here
can also be found in various papers by L. Lucy and in the main TARDIS publication
(cf. :cite:`Abbott1985`, :cite:`Mazzali1993`, :cite:`Lucy1999`,
:cite:`Long2002`, :cite:`Lucy2002`, :cite:`Lucy2003`, :cite:`Lucy2005`,
:cite:`Kerzendorf2014`).

.. toctree::
    :maxdepth: 2

    basicprinciples
    initialization
    propagation
    lineinteraction
    estimators