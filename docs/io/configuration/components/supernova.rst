.. _supernova-config:

***********************
Supernova Configuration
***********************

The supernova component of the configuration file contains some key information about the supernova being modeled, namely the time since the supernova and the luminosity the user wishes TARDIS should output:

.. jsonschema:: schemas/supernova.yml

As luminosity (in units of energy/s) is computed by integrating over the spectral luminosity (in units of energy/s/wavelength), TARDIS sums over all discrete energy packets to compute luminosity when ran, attempting to converge the output spectrum to match the luminosity requested (see :ref:`est_and_conv`), hence why `luminosity_wavelength_start` is set to 0 and `luminosity_wavelength_end` to infinity as a defualt. These values, however, can be changed depending on which wavelengths one has a luminosity. For example, if the user only wishes to see luminosities within the visible light range, then one can change `luminosity_wavelength_start` to 400 nm and `luminosity_wavelength_end` to 700 nm to obtain a spectrum within the visible light range. If the user wishes to see the entire spectrum, however, it is recommended to leave these parameters at their default value.   

