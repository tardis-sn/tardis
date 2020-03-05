Installing TARDIS
=================

Welcome to the Installation Guide for **TARDIS**. Please check the
`requirements <file:///home/harpreet/Tardis_sn/tardis/docs/_build/html/installation.html#requirements>`__
before the installation. You can read out the
`Troubleshooting <file:///home/harpreet/Tardis_sn/tardis/docs/_build/html/installation.html#troubleshooting-installation-faq>`__
section the **FAQ** as well for some common issues you might
encounter during installation.

.. note:: Consider these points before proceeding furthur:

1. TARDIS is only compatbile with Python Version above 3.6
2. Anaconda environment is strongly recommended to install TARDIS.
3. Try to use the latest github development version of Anaconda.
4. Linux is strongly recommended as TARDIS is not supported on Windows.

Requirements
------------

**Anaconda** is preffered for installing Tardis, You can follow these
`instructions <https://docs.continuum.io/anaconda/install/>`__ to
install anaconda.

The **environment definition file** contains the list of all the
requirements necessary for TARDIS. You can find the file
`here <https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env3.yml>`__.

TARDIS is using Astropy’s excellent installation helpers and thus uses
similar instructions to **Astropy**.

Installing TARDIS with Anaconda
-------------------------------

Anaconda python environment is strongly recommended for the installation
(or any other scientific packages for that matter). Anaconda benifits as
an isolated environment that can be set to be the default one, but by no
means will mess with your other environments. It will also work on
computers where ``root-``\ rights are not available. You can use these
instructions in order to install Anaconda on your local machine.

Step.1: Download the environment defination file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can get the environment defination file from
`here <https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis_env3.yml>`__.

Step.2: Create the environment:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Change to the directory that you downloaded the environment definition
file and run:

::

   conda env create -f tardis_env3.yml

Sit back and wait till the system prompts you.

Step.3: Activating the environment:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run either of the following commands to activate the environment

::

   source activate tardis

or use the new command

::

   conda activate tardis

Step.4: Deactivating the environment:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

when you feel that you are done with TARDIS the deactivation can be done
easily by running the following command

::

   source deactivate

You don’t need to recreate the environment next time you want to use
Tardis. Just activate it again using the same command:

::

   conda activate tardis    

That’s all with the installation of TARDIS. Click Next for furthur
instruction to **Develop TARDIS**

Troubleshooting Installation(FAQ)
---------------------------------

This section includes some of the common problems most prominent to be
encountered during installation and their solutions. It’s adviced to use
the recommended install method for any installation problem,

Problem: Could not find C file tardis/montecarlo/montecarlo.c
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While building tardis via ``python setup.py`` build you may encounter
the following error:

::

   error: tardis/montecarlo/montecarlo.c: Could not find C file tardis/montecarlo/montecarlo.c for Cython file tardis/montecarlo/montecarlo.pyx when building extension tardis.montecarlo.montecarlo. Cython must be installed to build from a git checkout.

Solution:
~~~~~~~~~

The most obvious solution to this is a clean checkout. Try running the
following commands:

::

   python setup.py clean

::

   git clean -dfx

**WARNING** This will delete any non tardis file in that directory and
will often clean this problem. If it still persists ``cd`` to
``tardis/montecarlo`` and buils montecarlo.c manually:

::

   cython montecarlo.pyx

Now ``python setup.py build`` will run perfectly.

Problem: limits.h: No such file or directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

when trying to set up CC=gcc python setup.py develop –with-openmp the
following error popped up: ``from tardis/_compiler.c:1``:

::

   /Users/user-name/miniconda2/envs/tardis-show2/lib/gcc/x86_64-apple-darwin13.4.0/5.2.0/include-fixed/limits.h:168:61: fatal error: limits.h: No such file or directory

Solution:
~~~~~~~~~

Try this Command:

::

   open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

Problem: Symbol not found: \_GOMP_parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Symbol not found: ``_GOMP_parallel`` when compiling with
–with-``openmp``

Solution:
~~~~~~~~~

Install ``gcc8`` from ``macports`` and then install with these flags:
``link_args = [‘-fopenmp’,’-Wl,-rpath,/opt/local/lib/gcc8/’]``

Problem: “python setup.py egg_info” failed with error code 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While building tardis(via python 2.7) via ``python setup.py build`` you
may encounter this error:

::

   TypeError: super() argument 1 must be type, not None

   ----------------------------------------
   Command "python setup.py egg_info" failed with error code 1 in /tmp/pip-req-build-wPB39p/

Solution:
~~~~~~~~~

The cause for this problem is sphinx , or sphinx version . It can be
easily solved by installing sphinx 1.5.6. try these commands:

::

   pip install sphinx==1.5.6

or

::

   conda install sphinx==1.5.6

You won’t face any errors in ``python setup.py build`` install now.
