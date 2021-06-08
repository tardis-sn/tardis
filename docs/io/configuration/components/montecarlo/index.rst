*************************
Monte Carlo Configuration
*************************

.. jsonschema:: ../schemas/montecarlo.yml

.. jsonschema:: ../schemas/montecarlo_definitions.yml#/definitions/convergence_strategy/damped

.. jsonschema:: ../schemas/montecarlo_definitions.yml#/definitions/convergence_strategy/custom

Parallel Execution with Numba
-----------------------------
The ``montecarlo`` section of the Configuration file accepts the parameter ``nthreads`` which sets the number of
threads to be used for parallelisation. Setting the value of the parameter between 1 and the environment variable
``NUMBA_NUM_THREADS`` (which is, by default, the number of CPU cores on your system) will automatically invoke Numba
to parallelise the code. (See :ref:`config-file` section).