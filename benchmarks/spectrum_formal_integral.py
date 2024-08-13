"""
Basic TARDIS Benchmark.
"""

import functools

from numba import config

from benchmarks.benchmark_base import BenchmarkBase
from tardis.spectrum import formal_integral

config.THREADING_LAYER = "workqueue"


class BenchmarkTransportMontecarloFormalIntegral(BenchmarkBase):
    """
    Class to benchmark the numba formal integral function.
    """

    repeat = 2

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
            self.sim.spectrum_solver.spectrum_frequency_grid
        )
        self.FormalIntegrator.make_source_function()
        self.FormalIntegrator.generate_numba_objects()
        self.FormalIntegrator.formal_integral(
            self.sim.spectrum_solver.spectrum_frequency_grid, 1000
        )
