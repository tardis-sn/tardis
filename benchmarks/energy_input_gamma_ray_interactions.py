"""
Basic TARDIS Benchmark.
"""
import numpy as np
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.energy_input import H_CGS_KEV, ELECTRON_MASS_ENERGY_KEV
from tardis.energy_input.GXPacket import GXPacketStatus, GXPacket
from tardis.energy_input.gamma_ray_interactions import pair_creation_packet, scatter_type


# @skip_benchmark
class BenchmarkEnergyInputGammaRayInteractions(BenchmarkBase):
    """
    Class to benchmark the gamma ray interactions function.
    """

    def __init__(self):
        self.basic_gamma_ray = GXPacket(
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

    def time_pair_creation(self):
        np.random.seed(2)
        self.basic_gamma_ray.nu_cmf = 2 * ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV
        pair_creation_packet(self.basic_gamma_ray)

    @parameterize({"Scatter type": [
        {
            "compton_opacity": 1,
            "photoabsorption_opacity": 0,
            "total_opacity": 1,
        },
        {
            "compton_opacity": 0,
            "photoabsorption_opacity": 1,
            "total_opacity": 1,
        },
        {
            "compton_opacity": 0,
            "photoabsorption_opacity": 0,
            "total_opacity": 1,
        },
    ]})
    def time_scatter_type(self, values: dict):
        scatter_type(
            values.get('compton_opacity'), values.get('photoabsorption_opacity'), values.get('total_opacity')
        )
