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

- `Numpy <http://www.numpy.org/>`_ 1.5 or later

- `Scipy <http://www.scipy.org/>`_ 0.10 or later

- `Astropy <http://www.astropy.org/>`_ 0.2.4 or later

- `h5py <http://www.h5py.org/>`_ 2.0.0 or later

- `pandas <http://pandas.pydata.org/>`_ 0.12.0 or later

- `pyyaml <http://pyyaml.org/>`_ 3.0 or later

Most of these requirements are easy to install using package managers like OS X's macports or normal linux package managers.

TARDIS is using astropy's excellent installation helpers and thus uses similar instructions to astropy.


Installing TARDIS
=================

On Ubuntu (13.10)
-----------------

We use a clean install of Ubuntu 13.10 as one of our testing grounds. Here's how we get TARDIS to run::

    sudo apt-get install python-dev python-pip python-numpy python-scipy python-h5py python-pandas python-yaml

We now need to install the newest astropy and we will install everything into our users directory::

    pip install astropy --user
    
Once astropy is installed, install TARDIS::

    pip install tardis-sn --user --pre

Add a `--pre` to install the latest development version (currently no stable version is available).


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


Building from source
====================

Prerequisites
-------------

You will need a compiler suite and the development headers for Python and
Numpy in order to build TARDIS. On Linux, using the package manager for your
distribution will usually be the easiest route, while on MacOS X you will
need the XCode command line tools.

The `instructions for building Numpy from source
<http://docs.scipy.org/doc/numpy/user/install.html>`_ are also a good
resource for setting up your environment to build Python packages.

You will also need `Cython <http://cython.org/>`_ installed to build
from source.

.. note:: If you are using MacOS X, you will need to the XCode command line
          tools.  One way to get them is to install `XCode
          <https://developer.apple.com/xcode/>`_. If you are using OS X 10.7
          (Lion) or later, you must also explicitly install the command line
          tools. You can do this by opening the XCode application, going to
          **Preferences**, then **Downloads**, and then under **Components**,
          click on the Install button to the right of **Command Line Tools**.
          Alternatively, on 10.7 (Lion) or later, you do not need to install
          XCode, you can download just the command line tools from
          https://developer.apple.com/downloads/index.action (requires an Apple
          developer account).

