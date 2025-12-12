************************
Basic General Properties
************************

After reading in the inputs, TARDIS begins its calculations of the plasma states by first calculating some basic properties that sets a foundation for more advanced calculations later. These properties generally involve properties of atoms and elements. 

Selected Atoms
==============
Dentoted as: selected-atoms
Calculated from: abundance

In order to obtain any calculations regarding an atom, its atomic number must be known, that is, the number of protons within the nucleus of the atom. When the atomic data file is read in via the :code:`atom-data` parameter in the configuration, TARDIS is made aware of the atomic number of every element within the file, so long as the file is formatted correctly, and must chose which atomic numbers correspond to the elements and isotopes stated within the :code:`abundace` field of the configuration. To do this, TARDIS simply in the index of the specified abundance list, giving it the atomic numbers of each element and isotope present within the simulation, while saving these numbers as a `Pandas <https://pandas.pydata.org/docs/>`_ index named :code:`selected-atoms`.

Atomic Mass
===========
Denoted as: atomic-mass
Calculated from: atomic-data, selected-atoms

After reading in the atomic number of each atom, TARDIS will parse through the atomic data to find the atomic mass (the amount of matter within the nucleus) in grams given its corresponding atomic number, saving these masses into the :code:`atomic-mass` Pandas series.


Ionization Data
===============
Denoted as: ionization-data
Calculated from: atomic-data, selected-atoms

The ionization energy of an atom is the energy necessary to remove a single electron from the outermost orbital of an atom. This removal will cause the atom to become a cation, and since the atom is now positively charged, then the next subsequent ionization energy to remove another electron will increase due to Coulombic forces. These energies are refered to by ordinal numbers, such as the First Ionization Energy, Second Ionization Energy, and so on, while the ordinal numbers themselves can be refered to as the Ion Number of the atom. Within a supernova, the thermal energy emitted will be enough to cause certain atoms to ionize due to the energy exceeding the atom's ionization energy, and it is vital for TARDIS to keep track of these energies to factor in how and when the atoms ionize during a simulation. Using the :code:`selected-atoms` list, TARDIS will parse through the atom data file which contains the ionization energies in ergs for each ion number of an atom and will select the energies corresponding to the atomic numbers it has saved, saving these energies and ion numbers for each selected atom within the :code:`ionization-data` Pandas series. 

Levels
======
Denoted as: levels
Calculated from: atomic-data, selected-atoms

In order to keep the atomic masses and ionization data for each atom within the simulation, TARDIS will place them into a single Pandas series to keep track of the number of energy levels and ion numbers per atomic number. In future calculations, TARDIS will parse through this series to acquire any necessary parameters.


Number Density
==============
Denoted as: :math:`N_i`
Calculated from: atomic-mass, :math:`\rho`

The number density, in the case of TARDIS, refers to the number of given atoms per unit volume. Within the configuration file, the density of the atoms is is make known to TARDIS via the specified density profile (see the Model Configuration).

Lines
=====
Denoted as: lines
Calculated from: atomic-data, selected-atoms

Within an atom lie many different energy levels which electrons can transition into, either when an atom absorbs a photon which results in a low-to-high transition or when an atom emits a photon resulting in a high-to-low transition. Due to quantum mechanical behaviors, however, these transitions do not always occur, so these transitions are measured via their probability of occuring. TARDIS, in creating and updating the plasma state, must take these transitions into account to create a spectrum with absoprtion and emission lines, so TARDIS will parse through the atomic data and the atoms within the abdunances section and find the following atomic number and ion numbers of each atom. TARDIS will then go through every possible energy level within each atom and record the following:

* A Numeric Identification for the energy level
* The wavelength of light associated with that energy level in angstrom and in centimeters 
* The probability for an electron in that level to transition to a higher (:code:`f_lu`) or lower (:code:`f_ul`) state
* The frequency of each level (denoted by :math:`\nu`) calculated from :math:`\dfrac{\lambda}{c}`
* The probabilities for a stimulated absorption and emission from atom-photon interactions (the B Einstein Coefficients)
* The probability for spontaneous emission (the A Einstein Coefficient)

TARDIS will then save all this information into a single Pandas dataframe entitled :code:`lines`, with :code:`nu`, :code:`f_lu`, :code:`wavelength_cm` as attributes to access them easily. 