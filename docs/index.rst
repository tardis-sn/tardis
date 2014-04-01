TARDIS
******

TARDIS is a Monte Carlo radiative-transfer spectral synthesis code for 1D models of supernova ejecta. It is designed for
rapid spectral modelling of supernovae. The code is described in this documentation and Kerzendorf & Sim 2014.

.. note::
    This documentation is currently under construction and does not describe all of the modes of operations available for TARDIS.


.. toctree::
    :maxdepth: 1

    installation
    running
    uses
    workflow/development_workflow
    examples/examples
    testing
    atomic/atomic_data
    plasma
    montecarlo
    changelog
    glossary
    zreferences
    credits

The code is built on a few principles:

 * **open** - the code is fully open source and we invite usage and contributions from the community
 * **modular** - the code has different microphysics modules and can be easily extended
 * **fast** - the code is geared towards rapid spectral synthesis to fit supernovae and other transients
 * **easy** - the code is designed to be easily installed and run as well as a detailed documentation


We encourage you to subscribe to `tardis-sn-user <http://groups.google.com/forum/#!forum/tardis-sn-users>`_ to ask questions about TARDIS.

If you use this code for any publications or presentations please acknowledge it accordingly. For this first version
please mention the website and cite Kerzendorf & Sim 2014.

User modifications and additions that lead to publications need to be handed back to the community by incorporating them
into this publicly available version of TARDIS.

The current stable version of TARDIS is 0.9 and can be downloaded `here <https://pypi.python.org/pypi/tardis-sn>`_, further installation instructions are
available here :ref:`installation`.

A file containing an example configuration file and an atomic database can be found in the section :ref:`running`

If you're interested in contributing to the code, either contact us or you can contribute directly via github.
We are using Astropy's excellent workflow - more details can be found at `<http://astropy.readthedocs.org/en/latest/development/workflow/index.html>`_.




.. warning::
    Currently TARDIS only works on 64-bit python installations. We're working on making it work on 32-bit python
    distributions.



..    configuration
..    gui







