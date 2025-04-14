.. _in-progress:
  
********************
Features In-Progress
********************

This is a current status summary of ongoing development projects in the TARDIS collaboration. Brief summaries of the projects are listed below. Other projects and their descriptions can be found in the `GitHub Projects <https://github.com/tardis-sn/tardis/projects>`_ or the `TARDIS Enhancement Proposals (TEPS) repository <https://github.com/tardis-sn/tep>`_.

Additionally, documentation for features that are not yet implemented can be found here.


Google Summer of Code (GSOC)
============================

Many new TARDIS developments are made through GSOC students. To learn more about potential GSOC projects, see the `GSOC ideas page <https://tardis-sn.github.io/gsoc/tardis-gsoc-ideas.html>`_.



Restructure
===========

Restructuring TARDIS to be more modularized in order to facilitate the addition of new physics modules and improve maintainability.


New Physics
===========

TARDIS is currently capable of simulating spherically symmetric type Ia and stripped envelope supernovae at a single point in time in the earlier stages of the explosion. To expand to more types of supernovae, such as type Ib, Ic, and IIP, and to model a larger variety of times after the explosion, new physics capabilities are necessary.


Simulate Supernovae in the Nebular Phase
----------------------------------------

* Currently, TARDIS assumes an inner boundary approximation to the photosphere. New physics is required to handle the breakdown of the photosphere at late times. This involves a detailed treatment of the decay products of the radioactive isotopes in SN ejecta, as well as their deposition in various processes such as Compton scattering, photoabsorption and pair production.
* Deposited energy needs to be thermalized by solving the Spencer-Fano equations, resulting in fractions of the energy going into heating, non-thermal excitation and non-thermal ionization.
* At late times when the densities are low collisions become too infrequent to quickly de-excite metastable energy levels. Forbidden lines arising from these levels need to be included. Transitions between levels have to include ion-electron collisions.

.. toctree::
    nebular_phase/gammaray_deposition
    nebular_phase/positronium


Implement More Continuum Interactions
-------------------------------------

* Currently, TARDIS considers Thomson electron scattering and line interactions when performing radiative transfer. In many cases, other continuum interactions have nontrivial effects of the propagation of light.
* Treatment of bound-free and free-free processes have been developed and are being integrated into TARDIS.


Simulate Non-Homologous Expansion
---------------------------------

* Currently, TARDIS models supernovae as expanding homologously. This is not always the case in many types of supernovae, and much of the simulation needs to be modified in order to drop this assumption.


Simulate Asymmetric Supernovae
------------------------------

* Currently, TARDIS assumes spherical symmetry in its calculations. A framework is being developed to allow for dropping this assumption and modeling 2D and 3D radiative transfer for spherically asymmetric supernovae.


Adding Time Dependence
----------------------

* Currently, TARDIS performs time-independent radiative transfer. Groundwork is being laid for the eventual goal of having TARDIS perform *time-dependent* radiative transfer.


Interoperability
================

Read Different Model Data
-------------------------

Adding parsers to generate TARDIS models produced by other codes.


STARDIS
=======

Developing a radiative transfer code using the TARDIS infrastructure that produces synthetic stellar spectra. See https://github.com/tardis-sn/stardis.
