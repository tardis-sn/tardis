**************
Example Models
**************

Basic Model specifications
==========================

TARDIS requires information about the ejecta structure. In particular, the
ejecta density and its composition have to be specified on a velocity grid.

.. note:: 

   Since homology is implicitly assumed, the velocity grid immediately translates into a radial mesh

The density and composition may be specified independently in the following ways

Specifying the density
----------------------

* :doc:`using a power-law density profile <densitypow>`
* :doc:`using an exponential density profile <densityexp>`
* :doc:`using a custom density profile <densitycust>`

Specifying the composition
--------------------------

* using a uniform composition
* using a custom stratified composition

Examples
========

Ejecta with Exponential density profile and uniform composition - tardis_example
--------------------------------------------------------------------------------


Detailed explosion model - mapping a delayed-detonation into TARDIS
-------------------------------------------------------------------


..
  Old Version
  ===========
  
  Here's an overview of some of the different modes of operation and models for TARDIS.
  
  
  .. toctree::
      :maxdepth: 1
  
      profileexp
  
      profilepl
  
      profileuniform
  
      profilemodel
  
  
  .. note:: 
  
     Increasing the number of virtual packets will improve the
     signal-to-noise of TARDIS spectra. You may wish to consider using a
     filter (e.g. Savitzky–Golay) to suppress the Monte Carlo noise for
     some applications.
  
  .. warning::
  
     The usefulness of any TARDIS calculations depends on the
     quality of the input atomic data. For further information on the atomic data -
     including details of how to develop your own dataset to suit your
     needs - please contact us.
