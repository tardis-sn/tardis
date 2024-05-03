"""
Faciliating the MonteCarlo iterations. 

During a simulation run, a number of MonteCarlo iterations specified
in the configuration are run using the numba compiler.
Most of the iterations are used to calculate the steady-state plasma 
properties and with the last iteration, the spectrum is determined.
"""

"""
Implements the main loop of the MonteCarlo routine.
"""

from llvmlite import binding

binding.set_option("tmp", "-non-global-value-max-name-size=2048")
njit_dict = {"fastmath": True, "error_model": "numpy", "parallel": True}
njit_dict_no_parallel = {
    "fastmath": True,
    "error_model": "numpy",
    "parallel": False,
}

from tardis.transport.montecarlo.r_packet import RPacket
from tardis.transport.montecarlo.montecarlo_main_loop import (
    montecarlo_main_loop,
)
from tardis.transport.montecarlo.packet_collections import (
    PacketCollection,
)
