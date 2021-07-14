.. _installation:

************
Installation
************

Before installing TARDIS, please check its :ref:`requirements
<requirements_label>`. We provide :ref:`instructions <anaconda_inst_label>` for installing TARDIS using 
Anaconda. If you encounter problems, consult the
:ref:`troubleshooting <troubleshooting_inst_label>` section. Once you have
installed TARDIS, check out the "Using TARDIS" section (in sidebar) for instructions regarding how to perform simple TARDIS calculations.

.. _requirements_label:


.. note::

    - TARDIS is only compatible with Python >=3.6
    - TARDIS only supports MacOS and Linux.
    - We strongly recommend installing TARDIS within an Anaconda environment and to always use the latest GitHub development version.


Requirements
============

You can see a list of all the requirements of TARDIS in the `environment definition file <https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env3.yml>`_.

TARDIS is using Astropy's excellent installation helpers and thus uses similar
instructions to Astropy.

.. _anaconda_inst_label:

Installing TARDIS with Anaconda
===============================

We highly recommend using the `Anaconda <https://www.anaconda.com/>`_ Python environment to install TARDIS (or
any other scientific packages for that matter). Anaconda has the advantage of
being an isolated environment that can be set to be the default one, but by no
means will mess with your other environments. It will also work on computers
where ``root``-rights are not available. Use these `instructions
<http://docs.continuum.io/anaconda/install.html>`_ to install Anaconda on your
machine. The next step is to create an environment for TARDIS that contains all
of the necessary packages (this ensures that TARDIS requirements won't clash
with any other Python installs on disc):

First, download the `environment definition file <https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env3.yml>`_ from:
::

    https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env3.yml

To create the environment, change to the directory that you downloaded the environment definition file and run:
::

    conda env create -f tardis_env3.yml

Then to activate this environment simply do:
::

    source activate tardis

or the new method:
::

    conda activate tardis

and after you are done with TARDIS you can deactivate:
::

    conda deactivate

One does not need to recreate the environment, but simply activate it every time
TARDIS is used.

For TARDIS development purposes please follow the steps :ref:`here <forking>`
until the step to install TARDIS in the development mode
``python setup.py develop``. Development guidelines for
TARDIS can be found `here <https://tardis-sn.github.io/tardis/development/index.html>`_.

To install TARDIS, it is recommended to first clone our repository and
then install TARDIS, as follows:
::

    git clone https://github.com/tardis-sn/tardis.git
    cd tardis
    python setup.py install



