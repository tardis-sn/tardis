.. _est_and_conv:

**************************
Estimators and Convergence
**************************

As light travels through a real plasma, it has effects on the properties of the plasma due to light-matter
interactions as well as the presence of extra energy from the light. Additionally, as :ref:`previously discussed
<propagation>`, properties of the plasma affect how light travels through it. Things like this (where two things
both affect each other on different ways) frequently occur in physics. The solution is finding a steady-state for
the plasma properties; that is, the actual plasma will be in a state such that the plasma state will not change as
light propagates through it, because the effects of the light on the plasma and the effects of the plasma on the
light have been "balenced out."

One of TARDIS's main goals is to determine this plasma state (as we need the actual plasma properties in order to
get an accurate spectrum). This is done in an iterative process. After each :ref:`Monte Carlo iteration
<montecarlo>` (which sends light through the supernova ejecta), TARDIS uses objects called estimators to determine
how the propagating light affects the plasma state, after which the plasma state is updated (as will be demonstrated
in the eatimators page linked below). We do this many times, and attempt to have the plasma state converge
to the steady-state we are looking for. In fact, all but the last Monte Carlo iteration is used for this purpose
(after which TARDIS will have the needed plasma state for its last iteration which calculates the spectrum).

.. note::
    For all but the last iteration, TARDIS uses the number of Monte Carlo packets specified in the
    :ref:`configuration file <montecarlo-config>` under the ``no_of_packets`` argument. This is because
    a different number of packets may be necessary to calculate the spectrum as opposed to calculate the
    plasma state.

.. toctree::

    estimators
    convergence