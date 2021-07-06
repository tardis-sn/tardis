.. _basic-spectrum:

*************************
Basic Spectrum Generation
*************************

The first and most simple way in which TARDIS calculates spectra calculates it directly from the Monte Carlo
packets after the final :ref:`Monte Carlo iteration <montecarlo>`. This simply requires knowledge of the each 
packet's energy and frequency in the lab frame (see :ref:`referenceframes`) at the end of the iteration. The only
other quantity needed is the time duration of the simulation :math:`\Delta t`, which isc alculated based off of the 
luminosity of the supernova's photosphere (see :ref:`initialization`).

.. note:: 
    
    The only packets which are used for this calculation are the packets which escape the outer boundary of the
    computational domain -- those reabsorbed into the photosphere are not included (see :ref:`propagation`).

The spectrum calculation is very straightforward. A packet of energy :math:`E_\mathrm{packet}` contributes a
luminosity

.. math:: L_\mathrm{packet} = \frac{E_\mathrm{packet}}{\Delta t}

to the spectrum at its frequency.