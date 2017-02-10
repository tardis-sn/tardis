.. _installation:

************
Installation
************

Before installing TARDIS, please check its :ref:`requirements
<requirements_label>`. We provide instructions for installing TARDIS using
    Anaconda :ref:`guide
<anaconda_inst_label>`. If you encounter problems, consult the
:ref:`troubleshooting <troubleshooting_inst_label>` section. Once you have
installed TARDIS, check out :doc:`running` for instructions of how to perform
simple TARDIS calculations.

.. _requirements_label:

Requirements
============

.. warning::
    TARDIS only works on 64-bit python installations.

You can see a list of all the requirements of TARDIS in the `environment definition file <https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env27.yml>`_.

TARDIS is using astropy's excellent installation helpers and thus uses similar
instructions to astropy.

.. _anaconda_inst_label:

Installing TARDIS with Anaconda
===============================

We highly recommend using the Anaconda python environment to install TARDIS (or
any other scientific packages for that matter). Anaconda has the advantage of
being an isolated environment that can be set to be the default one, but by no
means will mess with your other environments. It will also work on computers
where ``root``-rights are not available. Use these `instructions
<http://docs.continuum.io/anaconda/install.html>`_ to install Anaconda on your
machine. The next step is to create an environment for TARDIS that contains all
of the necessary packages (this ensures that TARDIS requirements won't clash
with any other python installs on disc::

First, download the `environment definition file <https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env27.yml>`_ from::

    https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env27.yml

To create the environment, change to the directory that you downloaded the environment definition file and run::

    conda-env create -f tardis_env27.yml

Then to activate this environment simply do::

    source activate tardis

and after you are done with TARDIS you can deactivate::

    source deactivate

One does not need to recreate the environment, but simply activate it every time
TARDIS is used.

Since TARDIS has reached a mature state, we recommend always installing the latest development version::

    pip install git+https://github.com/tardis-sn/tardis


To install the latest stable version of TARDIS simply do (usually outdated)::

    pip install tardis-sn




.. _troubleshooting_inst_label:

Installation Troubles (FAQ)
===========================

We highly encourage with any installation problems to try the recommended install
method because this often fix problems. Here are some common problems when
installing and their fixes:

**Problem:** While building tardis via ``python setup.py`` build you
may encounter the following error::

    error: tardis/montecarlo/montecarlo.c: Could not find C file tardis/montecarlo/montecarlo.c for Cython file tardis/montecarlo/montecarlo.pyx when building extension tardis.montecarlo.montecarlo. Cython must be installed to build from a git checkout.


**Solution:** There are several solutions to this problem. A clean checkout will
help. To clean up your repository please try ``python setup.py clean`` and
then ``git clean -dfx`` (**WARNING** will delete any non tardis file in that directory)
This will often clean this problem. If it still persists:

Go into the tardis/montecarlo directory and build montecarlo.c by hand::

    cython montecarlo.pyx

Then, ``python setup.py build`` should run without problems.
