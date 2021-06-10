**************************************************
Spectral Synthesis with the Formal Integral Method
**************************************************

In addition to generating the final Monte Carlo spectrum from the population of
Monte Carlo packets and the implemented variant of the "peeling-off" technique
(see :ref:`virtual_packets`), TARDIS supports spectral synthesis with
so-called formal integral method by :cite:`Lucy1999a` (see a detailed
description of the method at :ref:`formal_integral`)

Using the Formal Integral Method
================================

The formal integral spectral synthesis mode is activated by setting the

.. code-block:: none

    method: integrated 

configuration option in the ``spectrum`` block of the yml file.

Note that the integrated spectrum (``simulation.runner.spectrum_integrated``)
is a lazy property so it will be only generated (and then cached) once it is
accessed. The spectrum integration routine is Open-MP parallelized, so this
process may be significantly sped-up if TARDIS is built with the
``--with-openmp`` option and more then one thread is used. More instructions on
how to enable parallelisation of the code is given in the :ref:`parallelization`
section.
