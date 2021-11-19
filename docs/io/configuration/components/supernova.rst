.. _supernova-config:

***********************
Supernova Configuration
***********************

The supernova component of the configuration file contains some key information about the supernova being modeled, namelt the time since the supernova and the luminosity the user wishes TARDIS should output:

.. jsonschema:: schemas/supernova.yml

As luminosity (in units of energy/s) is computed by integrating over the spectral luminosity (in units of energy/s/wavelength), TARDIS sums over all discrete energy packets to compute luminosity when ran, hence why we set `luminosity_wavelength_start` to 0 and `luminosity_wavelength_end` to infinity.

During a run of TARDIS, we attempt to converge the output spectrum to match the requested luminosity
(see :ref:`est_and_conv`).
