.. _installation:

************
Installation
************


Requirements
============

.. warning::
    Currently TARDIS only works on 64-bit python installations. We're working on making it work on 32-bit python
    distributions.


TARDIS has the following requirements:

- `Python <http://www.python.org/>`_ 2.7

- `Numpy <http://www.numpy.org/>`_ 1.7 or later

- `Scipy <http://www.scipy.org/>`_ 0.10 or later

- `Astropy <http://www.astropy.org/>`_ 0.4 or later

- `h5py <http://www.h5py.org/>`_ 2.0.0 or later

- `pandas <http://pandas.pydata.org/>`_ 0.12.0 or later

- `pyyaml <http://pyyaml.org/>`_ 3.0 or later

Most of these requirements are easy to install using package managers like
OS X's macports or linux package managers or using the Anaconda python
distribution.

TARDIS is using astropy's excellent installation helpers and thus uses similar
instructions to astropy.

If you are interested in doing some development for TARDIS please read
:doc:`workflow/development_workflow`.


Once you are done you can run the simple example listed in :doc:`running`.




Installing TARDIS with Anaconda (recommended)
=============================================


We highly recommend using the Anaconda python environment to install TARDIS (or
any other scientific packages for that matter). Anaconda has the advantage of
being an isolated environment that can be set to be the default one, but
by no means will mess with your other environments. It will also work on
computers where ``root``-rights are not available. Use these
`instructions <http://docs.continuum.io/anaconda/install.html>`_ to install
Anaconda on your machine. The next step is to create an environment for tardis
that contains all of the necessary packages (this ensures that TARDIS
requirements won't clash with any other python installs on disc::

    conda create -n tardis --file https://raw.githubusercontent.com/tardis-sn/tardis/master/pip-requirements python=2 pip

Then to activate this environment simply do::

    source activate tardis

and after you are done with TARDIS you can deactivate::

    source deactivate

One does not need to recreate the environment, but simply activate it everytime
TARDIS is used.

To install the latest stable version of TARDIS simply do::

    pip install tardis-sn

or to use the development version::

    pip install git+https://github.com/tardis-sn/tardis

Installing TARDIS with virtualenvs
==================================


A virtual environment is python's way to ensure that the versions of third-party libraries
that TARDIS requires do not interfere with the system-wide installation. This
is also the way that the majority of core developers for TARDIS operate.

It is nevertheless recommended to install a number of python packages using the
system packagemanager. This ensures that third-party non-python libraries like
`libhdf5`, `lapack`, etc. are installed.

For OS X we recommend the `macports <http://www.macports.org/install.php>`_ package
manager::

    sudo port install python27 py27-scipy py27-numpy py27-virtualenv py27-astropy py27-h5py py27-yaml py27-pandas py27-pip py27-tables

In Ubuntu 14.04 the pre-requesite install works with this::

    sudo apt-get install python-virtualenv python-numpy python-pandas python-scipy python-h5py python-yaml ipython python-matplotlib cython git

First use the following guide to install virtualenv and the TARDIS requirements
:doc:`workflow/python_environment`.

After the virtualenv and the requirements are installed, there are a few options
of how to proceed.

If one just wants to use TARDIS you can install the latest stable version::

    pip install tardis-sn

the latest development version can be installed using::

    pip install git+https://github.com/tardis-sn/tardis

If you are interested in doing some development for TARDIS please read
:doc:`workflow/development_workflow`.


Once you are done you can run the simple example listed in :doc:`running`.



Installing TARDIS system-wide (not recommended)
===============================================

On Ubuntu (14.04)
-----------------

We use a clean install of Ubuntu 14.04 as one of our testing grounds. Here's how we get TARDIS to run::

    sudo apt-get install python-dev python-pip python-numpy python-scipy python-h5py python-pandas python-yaml

We now need to install the newest astropy and we will install everything into our users directory::

    pip install astropy --user
    
Once astropy is installed, install TARDIS::

    pip install tardis-sn

.. note::
    pip often tries to take care of many of the dependencies, this might be annoying as they already exist.
     Adding `--no-deps` will help with this problem.


On MAC OS X (10.8.5)
--------------------

On a clean install of Mountain Lion, here's how we get TARDIS running::

First install `macports <http://www.macports.org/install.php>`_

Use macports to install::

    sudo port install python27 py27-astropy py27-h5py py27-yaml py27-pandas py27-pip

Then install TARDIS::

    pip-2.7 install tardis-sn --user --pre

Before running, ensure that the directory ~/Library/Python/2.7/bin is in the appropriate path.

.. note::
    This has also been successfully tested on a clean MAC OS 10.9.1 (Mavericks) install.


