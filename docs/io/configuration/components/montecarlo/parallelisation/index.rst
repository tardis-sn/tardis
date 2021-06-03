.. _parallelisation:

*****************************************
Parallel execution with Numba and OpenMP
*****************************************

Enabling parallel execution with OpenMP
========================================

Manually, cloning the repository enables other options such as running the code in parallel (enabling OpenMP).
In general, we encourage downloading the compilers from `conda` as we then can ensure that they work with TARDIS.
Within the TARDIS conda environment do::

    conda install -c conda-forge compilers

For macOS::

    conda install -c conda-forge llvm-openmp

For Linux::

    conda install -c conda-forge openmp

To compile TARDIS for parallel execution::

    python setup.py install --with-openmp



Numba Usage Guide
=================

The ``montecarlo`` section of the Configuration file accepts the parameter ``nthreads`` which sets the number of
threads to be used for parallelisation. Setting the value of the parameter between 1 and the environment variable
``NUMBA_NUM_THREADS`` (which is, by default, the number of CPU cores on your system) will automatically invoke Numba
to parallelise the code. (See :ref:`config-file` section).