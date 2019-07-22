.. _modelconfig:

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

.. toctree::
    densitypow
    densityexp
    densitycust


Specifying the composition
==========================

.. toctree::
    abundanceuni
    abundancecust
