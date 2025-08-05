"""
Basic TARDIS Benchmark.
"""

import functools

from numba import config

from benchmarks.benchmark_base import BenchmarkBase
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
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
        self.formal_integrator = FormalIntegrator(
            self.sim.simulation_state, self.sim.plasma, self.sim.transport
        )

    # Benchmark for intensity black body function
    def time_intensity_black_body(self):
        nu = 1e14
        temperature = 1e4
        intensity_black_body(nu, temperature)

    # Benchmark for functions in FormalIntegrator class
    def time_FormalIntegrator_functions(self):
        self.formal_integrator.calculate_spectrum(
            self.sim.spectrum_solver.spectrum_frequency_grid
        )

        source_function_solver = SourceFunctionSolver(self.formal_integrator.transport.line_interaction_type)
        source_function_solver.solve(
            self.formal_integrator.simulation_state, 
            self.formal_integrator.opacity_state,  
            self.formal_integrator.transport.transport_state, 
            self.formal_integrator.plasma,
            self.formal_integrator.atomic_data)

        self.formal_integrator.generate_numba_objects()
        self.formal_integrator.formal_integral(
            self.sim.spectrum_solver.spectrum_frequency_grid, 1000
        )
