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

The three methods can be specified when running the simulation (see the :doc:`quickstart guide <../../../quickstart/quickstart>` for an example of this).

The following example shows how to edit variables for the integrated method. 

.. code-block:: yaml

        spectrum:
                start: 500 angstrom
                stop: 20000 angstrom
                num: 1000
                method: integrated
                integrated:
                        points: 2000
                        interpolate_shells: 100

Since the default method is set to virtual, the following shows a different example for editing the virtual method parameters.

.. code-block:: yaml

        spectrum:
                start: 5e2 angstrom
                stop: 2e4 angstrom
                num: 1000
                virtual:
                        tau_russian: 15
                        survival_probability: 0.1
                        enable_biasing: True
                        virtual_packet_logging: True

One can also change these parameters as they wish by reading in the configuration file and editing them before running the simulation (see :doc:`Reading a Configuration <../read_configuration>`).

.. warning::
    Currently, the "integrated" mode only works with the downbranching line
    interaction mode. Note also the limitations listed at the bottom of the
    dedicated page.
