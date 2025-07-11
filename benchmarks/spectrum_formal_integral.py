"""
Basic TARDIS Benchmark.
"""

import functools

from numba import config

from benchmarks.benchmark_base import BenchmarkBase
from tardis.spectrum.formal_integral.base import make_source_function
from tardis.spectrum.formal_integral.formal_integral import FormalIntegrator
from tardis.spectrum.formal_integral.formal_integral_numba import intensity_black_body

config.THREADING_LAYER = "workqueue"


class BenchmarkTransportMontecarloFormalIntegral(BenchmarkBase):
    """
    Class to benchmark the numba formal integral function.
    """

    repeat = 2

    @functools.cache
    def setup(self):
        self.sim = self.simulation_verysimple
        self.FormalIntegrator = FormalIntegrator(
            self.sim.simulation_state, self.sim.plasma, self.sim.transport
        )

    # Benchmark for intensity black body function
    def time_intensity_black_body(self):
        nu = 1e14
        temperature = 1e4
        intensity_black_body(nu, temperature)

    # Benchmark for functions in FormalIntegrator class
    def time_FormalIntegrator_functions(self):
        self.FormalIntegrator.calculate_spectrum(
            self.sim.spectrum_solver.spectrum_frequency_grid
        )
        make_source_function(
            self.FormalIntegrator.simulation_state, 
            self.FormalIntegrator.opacity_state,
            self.FormalIntegrator.transport,  
            self.FormalIntegrator.plasma, 
            self.FormalIntegrator.interpolate_shells
        )
        self.FormalIntegrator.generate_numba_objects()
        self.FormalIntegrator.formal_integral(
            self.sim.spectrum_solver.spectrum_frequency_grid, 1000
        )
