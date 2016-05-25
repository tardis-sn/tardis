**************************
Ejecta Model Configuration
**************************

Overview
========

TARDIS requires information about the ejecta structure. In particular, the
ejecta density and its composition have to be specified on a velocity grid.

.. note:: 

   Since homology is implicitly assumed, the velocity grid immediately translates into a radial mesh

The density and composition may be specified independently in the following ways

Specifying the density
======================

* :doc:`using a power-law density profile <densitypow>`
* :doc:`using an exponential density profile <densityexp>`
* :doc:`using a custom density profile <densitycust>`

Specifying the composition
==========================

* :doc:`using a uniform composition <abundanceuni>`
* :doc:`using a custom stratified composition <abundancecust>`
