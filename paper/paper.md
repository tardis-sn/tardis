---
title: 'TARDIS: a framework for transient radiative transfer'

tags:
- Python
- radiative transfer
- astrophysics

authors:
- name: Wolfgang E. Kerzendorf
  orcid: 0000-0002-0479-7235
  affiliation: 1

authors:
- name: and others
  orcid: XXXX
  affiliation: 1

affiliations:
- name: Department of Physics and Astronomy
  index: 1

bibliography: paper.bib
---

# Summary

TARDIS is a radiative-transfer framework using Monte-Carlo techniques to calculate synthetic spectra transients (with a specific focus). The code calculates a synthetic spectrum from a 1D model input assuming homologous expansion. The input model is structured into different shells with the required parameters of density and isotopic/elemental mass fractions. The current version of TARDIS then calculates the spectrum assuming a steady-state solution at a given time after explosion. The code scales the model (density and radii) and injects energy approximated by an inner boundary black-body of a given temperature.

# Detailed physics description

TARDIS uses a 1D input model (spherically symmetric shells; called cells in our documentation) that requires at the very least a density an elemental/isotopic mass fraction. The model can include more information such as radiative temperature, electron gas temperature, and ??? that will be taken into account at time of computation. Following the parsing of the input model, the code scales the densities and radii of the model to the appropriate time after explosion. We use the PYNE toolkit to calculate the decay chains for any isotopes that are potentially given and then generate a dataframe with elemental abundances for the given time after explosion. The team is currently working on an extension to include energy from nuclear decay into account (requiring gamma-ray deposition). 

TARDIS can calculate the plasma state using a variety of different approximation ranging from simple pure thermal approximation via semi-analytical to full Non-LTE calculations (currently only for excitation). The final step in the plasma calculation phase is to calculate optical depths (assuming a Sobolev approximation) for the radiative transfer step.

The Monte Carlo radiative transfer part uses an indivisible packet scheme in which each packet is described by a frequency, an energy (can be thought of as a collection of photons), a direction, and position in the ejecta. The packet can interact through several physical mechanisms (line interaction, thomson scattering) and can change direction on the path throught the ejecta. The final outcome is that the packet escapes the outer boundary is counted towards the spectrum or moves to the inner boundary and is reabsorbed and lost. The code tracks several packet statistics that are used to calculate the state of the radiation field (e.g. temperature).  

The updated radiation field parameters are then used in the plasma part to calculate a new plasma state. The MonteCarlo calculation and the plasma calculation are then performed iteratively until the plasma parameters converge. The code returns a converged model and several spectral outputs (specifically a formal integral solution).



# Acknowledgements

XXXX

# References