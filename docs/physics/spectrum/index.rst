.. _spectrum:

*******************
Spectrum Generation
*******************

During the final :ref:`Monte Carlo iteration <montecarlo>`, TARDIS calculates the emitted spectrum. We currently 
employ three diiferent methods for doing this: a basic spectrum generation directly from the Monte Carlo packets,
a method using so-called "virtual packets," and by the formal intrgral method. These are all detailed in the links
below.

.. note::
    For this last iteration, TARDIS uses the number of Monte Carlo packets specified in the :ref:`configuration file
    <montecarlo-config>` under the ``last_no_of_packets`` argument. Because this last iteration is used for
    calculating the actual spectrum, users may want to employ more packets than are used in the previous iterations.

.. toctree::
    basic
    virtualpackets
    sourceintegration