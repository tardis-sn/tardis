************
Installation
************

Requirements
============

TARDIS has the following strict requirements:

- `Python <http://www.python.org/>`_ 2.6, 2.7, 3.1 or 3.2

- `Numpy <http://www.numpy.org/>`_ 1.4 or later

- `Astropy <http://www.astropy.org/>`_ 0.2.4 or later



Installing TARDIS
=================

Using pip
---------

To install TARDIS with `pip`, first install astropy::

    pip install astropy

Once astropy is installed, install TARDIS::

    pip install tardis-sn

Add a `--pre` to install the latest development version.


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
from source, unless you are installing a numbered release. (The releases
packages have the necessary C files packaged with them, and hence do not
require Cython.)

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

