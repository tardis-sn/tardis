from llvmlite import binding

binding.set_option("tmp", "-non-global-value-max-name-size=2048")
<<<<<<< HEAD
njit_dict = {"fastmath": True, "error_model": "numpy", "parallel": True, 'nogil':True}
=======
njit_dict = {"fastmath": True, "error_model": "numpy", "parallel": False, 'nogil':True}
>>>>>>> 68942e53f1421b2d9edb3dfabd485f7515b96e65
njit_dict_no_parallel = {"fastmath": True, "error_model": "numpy", "parallel": False, 'nogil':True}

from tardis.montecarlo.montecarlo_numba.r_packet import RPacket
from tardis.montecarlo.montecarlo_numba.base import montecarlo_radial1d
from tardis.montecarlo.montecarlo_numba.numba_interface import PacketCollection
