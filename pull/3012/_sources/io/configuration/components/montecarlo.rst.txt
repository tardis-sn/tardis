.. _montecarlo-config:

*************************
Monte Carlo Configuration
*************************

The ``montecarlo`` section describes the parameters for the Monte Carlo radiation transport and convergence criteria:

.. jsonschema:: schemas/montecarlo.yml

The ``seed`` parameter seeds the random number generator first for the creation of the packets
(:math:`\nu` and :math:`\mu`) and then the interactions in the actual Monte Carlo process.
The ``no_of_packets`` parameter can take a float number for input convenience and gives the number of packets normally
used in each Monte Carlo loop. The parameters ``last_no_of_packets`` and ``no_of_virtual_packets`` influence the last run
of the Monte Carlo loop (which calculates the final spectrum!) when the radiation field should have converged. ``last_no_of_packets`` is normally higher than
``no_of_packets`` to create a less noisy output spectrum. ``no_of_virtual_packets`` can also be set to greater than 0 (a useful number is 3) to
use the Virtual Packet formalism. Increasing this number drastically increases computational costs (and memory requirements if they are logged).
The ``iterations`` parameter describes the maximum number of Monte Carlo loops executed in a simulation before it ends. Convergence criteria can be used to make the simulation stop
sooner when the convergence threshold has been reached (see :doc:`../../../physics/update_and_conv/update_and_conv`).
 
.. _conv-config:

Convergence Strategy
--------------------

The ``convergence_criteria`` section has a ``type`` keyword. Currently, one type is allowed: ``damped``.
All convergence criteria can be specified separately for the three variables for which convergence can be checked
(``t_inner``, ``t_rad``, ``ws``) by specifying subsections in the ``convergence_criteria`` of the same name. These then
override the defaults. Two more schemas are presented that further explain the ``damped`` and
``custom`` convergence strategies:

.. _damped-config:

.. jsonschema:: schemas/montecarlo_definitions.yml#/definitions/convergence_strategy/damped

``damped`` only has one parameter ``damping-constant`` and does not check for convergence. This can be used to fix the
temperature of the inner boundary.

.. jsonschema:: schemas/montecarlo_definitions.yml#/definitions/convergence_strategy/custom


.. _parallelization:

Parallel Execution with Numba
-----------------------------
The ``montecarlo`` section of the Configuration file accepts the parameter ``nthreads`` which sets the number of
threads to be used for parallelisation. Setting the value of the parameter between 1 and the environment variable
``NUMBA_NUM_THREADS`` (which is, by default, the number of CPU cores on your system) will automatically invoke Numba
to parallelise the code. (See :ref:`config-file` section).