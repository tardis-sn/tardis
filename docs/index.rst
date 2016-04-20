..
  .. image:: graphics/tardis_logo.jpg

.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

##################################
Tardis Core Package Documentation
##################################

.. |logo_svg| image:: graphics/tardis_banner.svg

.. raw:: html

   <img src="_images/tardis_banner.svg" onerror="this.src='_images/tardis_banner.svg'; this.onerror=null;" width="458"/>

.. _tardis:

TARDIS is an open-source Monte Carlo radiative-transfer spectral synthesis code
for 1D models of supernova ejecta. It is designed for rapid spectral modelling
of supernovae. If you use this code for any publications or presentations please
acknowledge it by citing Kerzendorf & Sim 2014.

============
Using Tardis
============

* :doc:`Installation <installation>`
* :doc:`Quickstart guide <running>`
* Examples
* Credit & Publication Policies

======================
Looking under the hood
======================

* Physics behind Tardis
* Monte Carlo primer

=================
Developing Tardis
=================

* Reporting issues
* Contributing


..
  .. note::
      This documentation is currently under construction and does not describe all of the modes of operations available for TARDIS.
  
  
  .. toctree::
      :maxdepth: 1
  
      installation
      running
      uses
      bugs
      configuration/index
      examples/examples
      testing
      atomic/atomic_data
      workflow/development_workflow
      physics/index
      changelog
      glossary
      zreferences
      credits
      license
  
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
  
  The current stable version of TARDIS is 0.9.2 and can be downloaded `here <https://pypi.python.org/pypi/tardis-sn>`_, further installation instructions are
  available here :ref:`installation`.
  
  A file containing an example configuration file and an atomic database can be found in the section :ref:`running`
  
  If you're interested in contributing to the code, either contact us or you can contribute directly via github.
  We are using Astropy's excellent workflow - more details can be found at `<http://astropy.readthedocs.org/en/latest/development/workflow/maintainer_workflow.html>`_.
  
  We encourage you to subscribe to `tardis-sn-user <http://groups.google.com/forum/#!forum/tardis-sn-users>`_ to ask questions about TARDIS.
  
  .. warning::
      Currently TARDIS only works on 64-bit python installations. We're working on making it work on 32-bit python
      distributions.
  
  
  
  ..    configuration
  ..    gui
