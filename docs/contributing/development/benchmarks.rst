.. _benchmarks:

**********
Benchmarks
**********

The objective of the benchmarking system is to detect regressions
that affect the performance of the TARDIS. This means we can visually
check whether there is a positive or negative spike in the TARDIS' performance.


AirSpeed Velocity (``ASV``)
===========================

TARDIS bases its benchmarking system on
`AirSpeed Velocity <https://asv.readthedocs.io/en/latest/index.html/>`_ (``ASV``).
Since it has a great advantage, which is that it is designed to run benchmarks on
`random servers <https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners/about-github-hosted-runners#standard-github-hosted-runners-for-public-repositories>`_
such as those provided by
`GitHub Actions <https://docs.github.com/en/actions>`_.
ASV eliminates the noise due to the technical differences between the servers
and produces graphs that indicate whether there is a regression or not.
It indicates the commit for the changes added or removed that affected performance
in some functions.


Installation
============

The complete guide is on the
`official ASV site <https://asv.readthedocs.io/en/latest/installing.html#installing-airspeed-velocity>`_,
however, here is detailed and summarized information to configure TARDIS with ASV.

ASV needs `Conda <https://docs.conda.io/en/latest/>`_
(or `Miniconda <https://docs.anaconda.com/free/miniconda/>`_)
and `Mamba <https://github.com/mamba-org/mamba>`_.
To make configuration easier, you can use
`Mini-forge <https://github.com/conda-forge/miniforge>`_,
which includes the installers mentioned above.
In this step, Mamba installs Python, ASV, and Mamba;
however, this step does not configure Mamba with the TARDIS.

.. code-block:: shell

    > export MAMBA_ENV_NAME="benchmark"
    > mamba create --yes --name "${MAMBA_ENV_NAME}" python asv mamba
    > mamba init


Set up
======

In this step, ASV configures TARDIS through Mamba.
Packages that use TARIDS are downloaded here.
These packages are mainly found in this ``tardis_env3.yml`` file.
The environment is also configured for ASV to execute benchmarks
and store the results through the ``asv.conf.json`` file.

.. code-block:: shell

    > cd tardis
    > export MAMBA_ENV_NAME="benchmark"
    > mamba activate "${MAMBA_ENV_NAME}"
    > asv setup
    > asv machine --yes


Execution
=========

ASV commands are used for execution. The first ``run`` command execute
the benchmarks found in the Python files that are in the ``benchmarks/``
folder. Subsequently, the data and information are stored in the ``.asv/`` folder.

.. code-block:: shell

    > cd tardis
    > export MAMBA_ENV_NAME="benchmark"
    > mamba activate "${MAMBA_ENV_NAME}"
    > asv run
    > asv publish


Visualization
=============

There are two ways to view the data. The simplest thing is
to execute the ``asv preview`` command, creating a local web server.
The second is to run a local web server of your choice.
