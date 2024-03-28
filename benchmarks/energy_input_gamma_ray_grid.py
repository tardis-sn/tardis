"""
Basic TARDIS Benchmark.
"""
import numpy as np
from asv_runner.benchmarks.mark import skip_benchmark

from tardis.energy_input import H_CGS_KEV
from tardis.energy_input.GXPacket import GXPacketStatus, GXPacket
from tardis.energy_input.gamma_ray_grid import move_packet


# @skip_benchmark
class BenchmarkEnergyInputGammaRayGrid:
    """
    Class to benchmark the gamma ray grid function.
    """

    def time_move_packet(self):
        packet = GXPacket(
            location=np.array([1.36375693e13, 4.10589818e14, 9.11718168e14]),
            direction=np.array([-0.97113853, 0.23134328, -0.05805379]),
            energy_rf=1e52,
            energy_cmf=1e52,
            nu_rf=1000.0e3 / H_CGS_KEV,
            nu_cmf=1000.0e3 / H_CGS_KEV,
            status=GXPacketStatus.IN_PROCESS,
            shell=1,
            time_current=1000,
        )
        distance = 1.0e15

        move_packet(packet, distance)
