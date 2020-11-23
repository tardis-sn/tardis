from llvmlite import binding

binding.set_option("tmp", "-non-global-value-max-name-size=2048")
njit_dict = {"fastmath": True, "error_model": "numpy"}

from tardis.montecarlo.montecarlo_numba.r_packet import RPacket
from tardis.montecarlo.montecarlo_numba.base import montecarlo_radial1d
from tardis.montecarlo.montecarlo_numba.numba_interface import PacketCollection
