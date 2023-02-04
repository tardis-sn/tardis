"""
Faciliating the MonteCarlo iterations. 

During a simulation run, a number of MonteCarlo iterations specified
in the configuration are run using the numba compiler.
Most of the iterations are used to calculate the steady-state plasma 
properties and with the last iteration, the spectrum is determined.
"""

from tardis.montecarlo.base import *
