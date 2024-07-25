"""
Basic TARDIS Benchmark.
"""

from asv_runner.benchmarks.mark import parameterize
from numba import config

import tardis.spectrum.formal_integral as formal_integral
from benchmarks.benchmark_base import BenchmarkBase

config.THREADING_LAYER='workqueue'

class BenchmarkTransportMontecarloFormalIntegral(BenchmarkBase):
    """
    Class to benchmark the numba formal integral function.
    """

    @parameterize(
        {
            "Parameters": [
                {
                    "nu": 1e14,
                    "temperature": 1e4,
                }
            ]
        }
    )
    def time_intensity_black_body(self, parameters):
        nu = parameters["nu"]
        temperature = parameters["temperature"]
        formal_integral.intensity_black_body(nu, temperature)

    # Benchmark for functions in FormalIntegrator class
    def time_FormalIntegrator_functions(self):
        FormalIntegrator = formal_integral.FormalIntegrator(
            self.simulation_verysimple.simulation_state, self.simulation_verysimple.plasma, self.simulation_verysimple.transport
        )
        FormalIntegrator.calculate_spectrum(self.simulation_verysimple.spectrum_solver.spectrum_real_packets.frequency)
        FormalIntegrator.make_source_function()
        FormalIntegrator.generate_numba_objects()
        FormalIntegrator.formal_integral(
            self.simulation_verysimple.spectrum_solver.spectrum_real_packets.frequency,
            1000
        )
