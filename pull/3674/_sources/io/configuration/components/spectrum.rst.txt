.. _spectrum-config:

**********************
Spectrum Configuration
**********************

Finally, the spectrum configuration tells TARDIS information needed for spectrum generation (see :ref:`spectrum`):

.. jsonschema:: schemas/spectrum.yml

``start`` and ``stop`` are given as wavelength values with units. They define
the boundaries of the frequency grid used for the real- and virtual-packet histograms and for
the generated spectrum. The histogram range is independent of the range that
controls which R-packets may spawn virtual packets.
``num`` specifies the number of bins used to build the spectrum and must be given as an integer. 
TARDIS produces the spectrum via three different methods. For more information on these methods, visit the
pages below:
 
* real: :doc:`Basic Spectrum Generation <../../../physics_walkthrough/spectrum/basic>`
* virtual: :ref:`Virtual Packet Scheme <virtual_packets>`
* integrated: :ref:`Formal Integral Method <formal_integral>`

The three methods can be specified when plotting the spectrum (see the :doc:`quickstart guide <../../../quickstart>` for an example of this).

The following example shows how to edit variables for the different methods. 

.. code-block:: yaml

        spectrum:
                start: 500 angstrom
                stop: 20000 angstrom
                num: 1000
                method: integrated
                integrated:
                        points: 2000
                        interpolate_shells: 100
                virtual:
                        tau_russian: 15
                        survival_probability: 0.1
                        enable_biasing: True
                        virtual_packet_logging: True

The virtual-packet spawn range is configured separately in the ``montecarlo``
section. It limits the wavelength of an R-packet when a virtual-packet volley
is considered; it does not change the histogram bin edges. The limits are
converted to frequency internally, so the lower wavelength limit corresponds
to the upper frequency limit::

        montecarlo:
                virtual_spectrum_spawn_range:
                        start: 1000 angstrom
                        end: 10000 angstrom
 

One can also change these parameters as they wish by reading in the configuration file and editing them before running the simulation (see :doc:`Reading a Configuration <../tutorial_read_configuration>`).


.. warning::
    As of now, the `method` argument serves no purpose other than adding 
    the integrated spectrum to the HDF output when "integrated" is used as the method
    (see :doc:`How to Store Simulations to HDF <../../../how_to/output/how_to_to_hdf>`). 
