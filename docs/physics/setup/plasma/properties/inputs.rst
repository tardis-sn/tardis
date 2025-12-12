******
Inputs
******

Before beginning plasma calculations, TARDIS will require a set of initial inputs as defined within the configuration file. Some inputs can be automatically calculated or are set to a default unchangable value by TARDIS, while others need an input specified by the user.

Required User Inputs
====================

Abundance
------------------
Denoted as: abundance
TARDIS will require the Abundance of elements and isotopes within the supernova. This is defined within the abundance section of the model configuration (see Model Configuration).


Atomic Data
-----------------
Denoted as: atomic-data
In addition to the mass abundances of elements and isotopes, TARDIS requires the Atomic Data to be able to define the properties of each element. This can be done defining the file containing the data within the configuration file via the :code:`atom_data` property and having the file within your working directory.


Density
-----------------
Denoted as: :math:`\rho`
The Density of each shell is necessary to perform calculations as well. This is defined within the density section of the Model Configuration (see the Model page for more information on density).


Time Since Explosion
----------------
Denoted as: :math:`t_{exp}`
TARDIS will need the time since the explosion in units of time primarily to calculate :math:`\tau_{\textrm{sobolev}}`. This time value can be defined within the supernova section of the model configuration (see Supernova Configuration).


Optional User Inputs
==================

Radiative Temperature
-----------------
Denoted as: :math:`T_{rad}`
An Initial Radiative Temperature (in Kelvin) is required by TARDIS, which can be defined within the plasma configuration (see Plasma Configuration). This value will be updated throughout the simulation (see the Convergence Section). If no initial radiative temperature is defined, TARDIS will assign a default value of -1 to automatically calculate the initial temperatures. For more information on the role of radiative temperature, see the Model page.

Dilution Factor
-----------------
Denoted as: :math:`W`
TARDIS requires, as well, an initial fractional energy density at the observer compared to the blackbody known as the Dilution Factor, which can either be defined per shell within a CSVY model configuration (see CSVY Model) or can be left to be automatically calculated by TARDIS (see Model for more details). The dilution factor will be updated throughout the simulation (see the Convergence Section). 


LinkTRadTElectron
------------------
Denoted as: :math:`T_{electron}/T_{rad}`
Within its calculations, TARDIS will compute the Electron Temperature by using the Radiative Temperature and the ratio between the radiative temperature (also known as the photospheric effective temperature) and the electron temperature, denoted as :code:`LinkTRadTElectron`. TARDIS sets this ratio to a default value of 0.9 when assembling the plasma states to assume an isothermal flow of stellar winds as per :cite:`MazzaliLucy93`. 

.. note::
    This parameter cannot be edited or changed by the user.

    
Helium Treatment
----------------
Denoted as: HeliumTreatment
The method TARDIS uses to handle how helium is treated within the simulation depends on the specified treatment type. As a default, TARDIS will treat helium as it does with all other elements, however within the plasma configuration, one can specify TARDIS use NLTE approximations when handling helium (see Plasma Configuration).

.. note::
    Specifying any type of approximation will require a 
    heating rate/light curve data file within the working directory
    

Continuum Interaction Species
-----------------
Denoted as: ContinuumInteractionSpecies
Certain isotopes can be requested to be treated with ionization and recombination if a list of the specified isotopes are passed through the species subsection under the continuum_interaction section of the plasma configuration (see Plasma Configuration). 

.. warning::
    This feature as of now does not run correctly under inputted species and
    should be left in its default state. 


.. note::
    Both :code:`HeliumTreatment` and :code:`Continuum Interaction Species` are used 
    for treatments of specific elements within the span of the simulation 
    and thus do not contribute to any one property directly.