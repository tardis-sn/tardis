"""
Basic TARDIS Benchmark.
"""

import functools

from numba import config


import tardis.spectrum.formal_integral as formal_integral
from benchmarks.benchmark_base import BenchmarkBase

config.THREADING_LAYER = "workqueue"


class BenchmarkTransportMontecarloFormalIntegral(BenchmarkBase):
    """
    Class to benchmark the numba formal integral function.
    """

    @functools.cache
    def setup(self):
        self.sim = self.simulation_verysimple
        self.FormalIntegrator = formal_integral.FormalIntegrator(
            self.sim.simulation_state, self.sim.plasma, self.sim.transport
        )

    # Bencmark for intensity black body function
    def time_intensity_black_body(self):
        nu = 1e14
        temperature = 1e4
        formal_integral.intensity_black_body(nu, temperature)

    # Benchmark for functions in FormalIntegrator class
    def time_FormalIntegrator_functions(self):
        self.FormalIntegrator.calculate_spectrum(
            self.sim.spectrum_solver.spectrum_real_packets.frequency
        )
        self.FormalIntegrator.make_source_function()
        self.FormalIntegrator.generate_numba_objects()
        self.FormalIntegrator.formal_integral(
            self.sim.spectrum_solver.spectrum_real_packets.frequency, 1000
        )
