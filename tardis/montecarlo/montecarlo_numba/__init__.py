from llvmlite import binding
binding.set_option("tmp", "-non-global-value-max-name-size=2048")
njit_dict = {'fastmath': False}

from tardis.montecarlo.montecarlo_numba.rpacket import RPacket
from tardis.montecarlo.montecarlo_numba.base import montecarlo_radial1d