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
    - TARDIS only supports Linux and MacOS.
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
with any other Python installs on disc)::

First, download the `environment definition file <https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env3.yml>`_ from::

    https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env3.yml

To create the environment, change to the directory that you downloaded the environment definition file and run::

    conda env create -f tardis_env3.yml

Then to activate this environment simply do::

    source activate tardis

or the new method::

    conda activate tardis

and after you are done with TARDIS you can deactivate::

    source deactivate

One does not need to recreate the environment, but simply activate it every time
TARDIS is used.

For TARDIS development purposes please follow the steps :ref:`here <forking>`
until the step to install TARDIS in the development mode
``python setup.py develop``. Development guidelines for
TARDIS can be found `here <https://tardis-sn.github.io/tardis/development/index.html>`_.

To install TARDIS, it is recommended to first clone our repository and
then install TARDIS, as follows::

    git clone https://github.com/tardis-sn/tardis.git
    cd tardis
    python setup.py install


.. _install_openmp:

Enabling parallel execution with OpenMP
---------------------------------------


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


.. _troubleshooting_inst_label:

Installation Troubles (FAQ)
===========================

We highly encourage with any installation problems to try the recommended install
method because this often fixes problems. Here are some common problems when
installing and their fixes:

**Problem:** While building TARDIS via ``python setup.py`` build you
may encounter the following error::

    error: tardis/montecarlo/montecarlo.c: Could not find C file tardis/montecarlo/montecarlo.c for Cython file tardis/montecarlo/montecarlo.pyx when building extension tardis.montecarlo.montecarlo. Cython must be installed to build from a git checkout.


**Solution:** There are several solutions to this problem. A clean checkout will
help. To clean up your repository please try ``python setup.py clean`` and
then ``git clean -dfx`` (**WARNING** will delete any non-TARDIS file in that directory)
This will often clean this problem. If it still persists:

Go into the tardis/montecarlo directory and build montecarlo.c by hand::

    cython montecarlo.pyx

Then, ``python setup.py build`` should run without problems.


**Problem:** when trying to set up CC=gcc python setup.py develop --with-openmp the following error popped up: 
from tardis/_compiler.c:1: /Users/yssavo/miniconda2/envs/tardis-show2/lib/gcc/x86_64-apple-darwin13.4.0/5.2.0/include-fixed/limits.h:168:61: fatal error: limits.h: No such file or directory 
        
**Solution:** Run on terminal: 

    open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

**Problem:** Symbol not found: _GOMP_parallel when compiling with `--with-openmp`

**Solution:** Install gcc8 from macports and then install with these flags: `link_args = ['-fopenmp','-Wl,-rpath,/opt/local/lib/gcc8/']`

**Problem:** While building TARDIS (via python 2.7) via ``python setup.py`` build you
may encounter the following error::

     TypeError: super() argument 1 must be type, not None
    
    ----------------------------------------
    Command "python setup.py egg_info" failed with error code 1 in /tmp/pip-req-build-wPB39p/


**Solution:** The cause for this problem is Sphinx or Sphinx version. It can be easily solved by installing Sphinx 1.5.6.
              The command for the same is :

    pip install sphinx==1.5.6
    
    or
    
    conda install sphinx==1.5.6

Then, ``python setup.py build install`` should run without problems.
