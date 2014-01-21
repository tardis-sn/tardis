TARDIS
******

TARDIS is a Monte Carlo radiative-transfer spectral synthesis code for 1D models of supernova ejecta. It is designed for
rapid spectral modelling of supernovae. The code is described in Kerzendorf & Sim 2014.

The code is built on a few principles:

 * **open** - the code is fully open source and we invite usage and contributions from the community
 * **modular** - the code has different microphysics modules and can be easily extended
 * **fast** - the code is geared towards rapid spectral synthesis to fit supernovae and other transients
 * **easy** - the code is designed to be easily installed and run as well as a detailed documentation


We encourage you to subscribe to `tardis-sn-user <http://groups.google.com/forum/#!forum/tardis-sn-users>`_, which also contains the developers.

If you use this code for any publications or presentations please acknowledge it accordingly. For this first version
please mention the website and cite Kerzendorf & Sim 2014.

The current stable version of TARDIS is 0.9 and can be downloaded `here <https://pypi.python.org/pypi/tardis-sn>`_, further installation instructions are
available here :ref:`installation`.

A package containing an example file and a atomic database can be found in the section :ref:`running`

.. warning::
    Currently TARDIS only works on 64-bit python installations. We're working on making it work on 32-bit python
    distributions.

.. toctree::
    :maxdepth: 1

    installation
    running
    gui
    atomic
    plasma
    montecarlo
    glossary
    zreferences
    credits

..    configuration







