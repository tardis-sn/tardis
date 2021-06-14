.. _spectrum-config:

**********************
Spectrum Configuration
**********************

Finally, the spectrum configuration tells TARDIS information needed for spectrum generation (see :ref:`spectrum`):

.. jsonschema:: schemas/spectrum.yml

``Start`` and ``end`` are given as Quantities with units. If they are given in
frequency space they are switched around if necessary. The number of bins is
just an integer. Finally, the method option selects the final spectral synthesis mode. Currently, there are three options:
 
* real: construct spectrum from the real packet population alone
* virtual: use the :ref:`virtual packet scheme <virtual_packets>` for spectral synthesis
* integrated: use the :ref:`formal integral method <formal_integral>` of Lucy 1999
 
.. warning::
    Currently, the "integrated" mode only works with the downbranching line
    interaction mode. Note also the limitations listed at the bottom of the
    dedicated page.
