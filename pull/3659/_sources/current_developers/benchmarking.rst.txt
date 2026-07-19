Benchmarks
----------

.. _explanation-benchmarking:

Explanation: Benchmarking
~~~~~~~~~~~~~~~~~~~~~~~~~

The benchmarking system detects performance regressions in TARDIS. It lets
developers visually check whether performance has improved or worsened.

TARDIS uses AirSpeed Velocity, or ASV. ASV is designed to run benchmarks on
random servers, such as GitHub-hosted runners, and reduce noise caused by
technical differences between servers. ASV produces graphs that indicate
whether a regression occurred and can identify commits that affected performance
in specific functions.

Benchmark files live in the ``benchmarks/`` directory. Results are stored under
``.asv/``.

.. _reference-benchmark-command-reference:

Reference: Benchmark Command Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Common ASV commands:

.. code-block:: shell

   asv setup
   asv machine --yes
   asv run
   asv publish
   asv preview


Benchmark classes live under ``benchmarks/``; for example,
``benchmarks/spectrum_formal_integral.py`` defines
``BenchmarkTransportMontecarloFormalIntegral`` with ASV methods such as
``time_intensity_black_body``.

.. _how-to-guide-run-benchmarks:

How-To Guide: Run Benchmarks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS uses AirSpeed Velocity, or ASV, for benchmarks.

Install ASV and its required environment tooling. ASV needs Conda or Miniconda
and Mamba. Mini-forge includes these installers and can simplify configuration.

Create the benchmark environment:

.. code-block:: shell

   export MAMBA_ENV_NAME="benchmark"
   mamba create --yes --name "${MAMBA_ENV_NAME}" python asv mamba
   mamba init


Set up ASV for TARDIS:

.. code-block:: shell

   cd tardis
   export MAMBA_ENV_NAME="benchmark"
   mamba activate "${MAMBA_ENV_NAME}"
   asv setup
   asv machine --yes


Run and publish benchmarks:

.. code-block:: shell

   cd tardis
   export MAMBA_ENV_NAME="benchmark"
   mamba activate "${MAMBA_ENV_NAME}"
   asv run
   asv publish


Preview benchmark output:

.. code-block:: shell

   asv preview


You can also view the generated data with a local web server of your choice.

When iterating on one benchmark file, inspect the benchmark name in
``benchmarks/spectrum_formal_integral.py``, then run a matching ASV benchmark
locally before publishing results:

.. code-block:: shell

   asv run -b time_intensity_black_body
   asv preview
