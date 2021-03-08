.. _roadmap:
  
*******************
Development Roadmap
*******************

This is a current status summary of ongoing development projects
in the TARDIS collaboration. Brief summaries of the projects as
well as timelines are listed below. More detailed descriptions
of the projects can be found in the Github Projects or the TARDIS
Enhancement Proposals (TEPS) repository under the TARDIS-SN
collaboration github page.

Physics
#######

* Simulate Type IIP Supernovae  (Fall 2020)

 - TARDIS is currently capable of simulating type Ia and stripped
   envelope supernovae, but Type IIP require new physics.
 - NLTE ionization and excitation treatment of Hydrogen, relativistic
   effects, free-free and bound-free processes have all been implemented
   and are being integrated into the main TARDIS branch.
   
* Simulate Supernovae in the Nebular Phase (Summer 2021)

 - Currently TARDIS assumes an inner boundary approximation to the
   photosphere.  New physics is required to handle the breakdown of the
   photosphere at late times. This involves a detailed treatment of the
   decay products of the radioactive isotopes in SN ejecta, as well as
   their deposition in various processes such as Compton scattering,
   photoabsorption and pair production.
 - Deposited energy needs to be thermalized by solving the Spencer-Fano
   equations, resulting in fractions of the energy going into heating,
   non-thermal excitation and non-thermal ionization.
 - At late times when the densities are low collisions become too
   infrequent to quickly de-excite metastable energy levels. Forbidden
   lines arising from these levels need to be included. Transitions between
   levels have to include ion-electron collisions.

Website and Documentation
#########################

* Notebookify (Fall 2021)

 - TARDIS currently has detailed code and physics documentation
   (https://tardis-sn.github.io/tardis/)
 - In order to help new members join the TARDIS development team, we want
   to create interactive jupyter notebooks that are incorporated into the
   documentation to better illustrate how the code is used and show how
   the physics of TARDIS works.
 - Apply to Google Summer of Docs to get a dedicated, professional
   documentation writer.

Analysis
########

* Visualization (Fall 2020)

 - TARDIS is used to explore potential models for observed supernovae,
   and it is often difficult to find optimal parameter values.
 - We are working on an interactive Graphical User Interface for
   TARDIS that will allow users to visually inspect their models and
   obtain physical information from the TARDIS simulation easily.
 - Most of the analysis will be used via Jupyter Notebooks.
 - This is a Google Summer of Code project
   (see https://github.com/tardis-sn/tardis/projects/5 ) with time to
   completion aimed Fall 2020.

Maintenance
###########

* TARDIS Code Restructure (Fall 2020)

 - TARDIS needs to be more modularized to facilitate the addition of
   new physics modules (see https://github.com/tardis-sn/tardis/projects/11 )
   and improve maintainability.
 - Google Summer of Code project is currently underway to restructure the
   model classes that TARDIS uses.
 - Cython and C code for monte carlo simulation is currently being replaced
   with Numba to facilitate modularity and debugging.

Performance
###########

* Numba (Spring 2021)

 - In the past, TARDIS has used Cython and C for performance-critical parts
   of the code. However, this requires significant code overhead and is difficult
   to maintain given that many of our developers are junior researchers who may
   be less familiar with C code.
 - As one of TARDIS’ key focuses is its use as an educational tool, we would like
   to minimize the complexity of the code so that it can be understood by a wider
   audience.
 - We are currently in the process of converting this Cython and C code to Numba,
   which we hope will allow us to maintain the performance of Cython but with code
   that is now exclusively written in (numba-fied) Python.
   
* GPU (Spring 20201)

 - As part of our work to optimize TARDIS’ performance, we would like to explore
   GPUs and see if there are parts of the code which would benefit from GPU
   acceleration.
 - We are looking into using CUDA-accelerated Numba for this (does TARDIS use
   matrices a lot? CuPy might also be useful here).

Interoperability
################

* Carsus (Fall 2020)

 - There is evidence that the choice of atomic data has a significant impact on
   model results. As part of Google Summer of Code a version control system for
   atomic data is being implemented, so that previously published results become
   reproducible.
 - Carsus will include readers for the biggest atomic data compilations such as
   CMFGEN, Chianti and Kurucz, and output the data in the Tardis format.
 - Updated properties of the atomic data will be automatically collected and
   ingested by Carsus. In addition, a comparison spectrum will be computed,
   highlighting the effects of the new atomic data.
   
* Read different model data (Fall 2020)

 - Add readers for more types of initial model data.