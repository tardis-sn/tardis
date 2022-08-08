.. _supernova-config:

***********************
Supernova Configuration
***********************

The supernova component of the configuration file contains some key information about the supernova being modeled, namely the time since the supernova (which is used throughout TARDIS calculations) and the luminosity the user wishes TARDIS should output:

.. jsonschema:: schemas/supernova.yml

As luminosity (in units of energy/s) is computed by integrating over the spectral luminosity (in units of energy/s/wavelength), TARDIS sums over all discrete energy packets to compute luminosity when ran, attempting to converge the output spectrum to match `luminosity_requested` (see :doc:`../../../physics/update_and_conv/update_and_conv`). `luminosity_requested` can be given in standard units, such as erg/s or J/s, or in logarithmic units such as log_lsun. The range over which TARDIS sums these energy packets is set by default from 0 to infinity via `luminosity_wavelength_start` and `luminosity_wavelength_end`, respectively, so as to generate a spectrum whose total luminosity across the entire spectrum is `luminosity_requested`. However, if in the event only the luminosity within a certain range of wavelengths is known, then `luminosity_wavelength_start` and `luminosity_wavelength_end` can be changed as necessary to reflect this, allowing TARDIS to attempt to create a spectrum whose luminosity within the set range will converge to the value defined in `luminosity_requested`.

As an example, here is a sample code which will generate a specturm where only the luminosity of the visible light portion of the spectrum is given. Here, the output spectrum will have a luminosity of approximately :math:`10^{9.44}L_{sun}` within the visible range.  

.. code-block:: yaml
        
    supernova:
      luminosity_requested: 9.44 log_lsun
      time_explosion: 13 day
      luminosity_wavelength_start: 400 nm
      luminosity_wavelength_end: 700 nm


