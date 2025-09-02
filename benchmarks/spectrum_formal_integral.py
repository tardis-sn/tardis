"""
Basic TARDIS Benchmark.
"""

from numba import config

from benchmarks.benchmark_base import BenchmarkBase
from tardis.spectrum.formal_integral.base import intensity_black_body
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)

config.THREADING_LAYER = "workqueue"


class BenchmarkTransportMontecarloFormalIntegral(BenchmarkBase):
    """
    Class to benchmark the numba formal integral function.
    """

    repeat = 2

    def setup(self):
        self.sim = self.simulation_verysimple
        integrator_settings = self.sim.spectrum_solver.integrator_settings
        self.formal_integral_solver = FormalIntegralSolver(
            integrator_settings.points,
            integrator_settings.interpolate_shells,
            getattr(integrator_settings, "method", None),
        )

    # Benchmark for intensity black body function
    def time_intensity_black_body(self):
        nu = 1e14
        temperature = 1e4
        intensity_black_body(nu, temperature)

    # Benchmark for functions in FormalIntegrator class
    def time_FormalIntegrator_functions(self):
        sim_state = self.sim.simulation_state
        transport_solver = self.sim.transport
        plasma = self.sim.plasma
        nu = self.sim.spectrum_solver.spectrum_frequency_grid[:-1]

        # Solve the formal integral - setup is called internally
        self.formal_integral_solver.solve(
            nu,
            self.sim.simulation_state,
            self.sim.transport,
            self.sim.opacity_state,
            self.sim.plasma.atomic_data,
            self.sim.plasma.electron_densities,
            self.sim.macro_atom_state
        )
