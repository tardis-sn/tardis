******************************
Parallel execution with OpenMP
******************************

.. note:: TARDIS no longer uses C, and thus OpenMP is not relevant to TARDIS.

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