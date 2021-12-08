.. _supernova-config:

***********************
Supernova Configuration
***********************

The supernova component of the configuration file contains some key information about the supernova being modeled, namely the time since the supernova and the luminosity the user wishes TARDIS should output:

.. jsonschema:: schemas/supernova.yml

As luminosity (in units of energy/s) is computed by integrating over the spectral luminosity (in units of energy/s/wavelength), TARDIS sums over all discrete energy packets to compute luminosity when ran, attempting to converge the output spectrum to match the luminosity requested (see :ref:`est_and_conv`). The range over which TARDIS sums these energy packets is set from 0 to infinity via the `luminosity_wavelength_start` and `luminosity_wavelength_end`, respectively, to generate a spectrum whose luminosity is known across the entire spectrum by default. However, if in the event only the luminosity within a certain range of wavelengths is known, then `luminosity_wavelength_start` and `luminosity_wavelength_end` can be changed as necessary to reflect this, allowing TARDIS to attempt to create a spectrum whose luminosity within the set range will converge to the value defined in `luminosity_requested`.

As an example, here is a sample code which will generate a specturm where only the luminosity of the visible light portion of the spectrum is given.

.. code-block:: yaml
        
    supernova:
      luminosity_requested: 9.44 log_lsun
      time_explosion: 13 day
      luminosity_wavelength_start: 400 nm
      luminosity_wavelength_end: 700 nm

.. warning::
    If `luminosity_wavelength_start` and `luminosity_wavelength_end` are given in terms of frequency, 
    the larger frequency should be set as the start and the lower frequency should be set as the end
    since TARDIS will convert these into wavelengths.

