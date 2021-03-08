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

affiliations:
- name: Department of Physics and Astronomy
  index: 1
bibliography: tardis.bib
---

# Summary

TARDIS is a radiative-transfer framework using Monte-Carlo techniques to calculate synthetic spectra for transients. The code is an open-science implementation of the indivisible energy-packet Monte Carlo formalism [@Abbott1985].  The code calculates a synthetic spectrum from a 1D model input assuming homologous expansion. The input model is structured into different shells with the required parameters of density and isotopic/elemental mass fractions. The current version of TARDIS then calculates the spectrum assuming a steady-state solution at a given time after explosion. The code scales the model (density and radii) and injects energy approximated by an inner boundary black-body of a given temperature.

# Statement of Need

Supernovae are critical components of the Universe and impact large fields of astronomy. The community has invested considerable resources to acquire high-quality data of many of these transient phenomena. However, the analysis of these data (especially by comparison to highly complex simulations) has been lacking as only very few scientists can run these complex codes. 

Answering the question of the connection between progenitor stars, explosion mechanisms, and the wide range of observed SNe requires models that allow us to investigate the physical processes that are involved. Stellar evolution codes simulate the process across long time-scales from stellar birth to the point of explosions. Spectral synthesis is a critical component of the aforementioned link between SNe and the processes that produce them. Specifically, there are large collections of supernova spectra available from dedicated repositories (cite wise; cite sne.space).

TARDIS aims to provide a community-developed open-science based code designed to synthesize spectra for transients from parameterized abundance and density models. The code is structured in a modular way that allows to use different microphysics approximations for different applications (e.g. an analytical NLTE solution for computational expedience, a numerical NLTE solution for accuracy). The TARDIS framework also provides an open reference implementation of the widely used indivisible energy packet Monte Carlo method which is used in several other radiative transfer codes. 

# Acknowledgements

The TARDIS community would like to acknowledge the continued support from the Google Summer of Code initiative and the Python Software Foundation. We have received support from NumFOCUS. The PI was supported by an ESO Fellowship for part of the project.



# References