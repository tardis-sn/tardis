from llvmlite import binding

binding.set_option("tmp", "-non-global-value-max-name-size=2048")
njit_dict = {"fastmath": True, "error_model": "numpy", "parallel": True}
njit_dict_no_parallel = {
    "fastmath": True,
    "error_model": "numpy",
    "parallel": False,
}

from tardis.montecarlo.transport.r_packet import RPacket
from tardis.montecarlo.transport.base import montecarlo_radial1d
from tardis.montecarlo.transport.numba_interface import PacketCollection
