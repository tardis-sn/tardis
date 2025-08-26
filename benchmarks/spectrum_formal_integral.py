"""
Basic TARDIS Benchmark.
"""

import functools

from numba import config

from benchmarks.benchmark_base import BenchmarkBase
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)
from tardis.spectrum.formal_integral.formal_integral_numba import (
    intensity_black_body,
)

config.THREADING_LAYER = "workqueue"


class BenchmarkTransportMontecarloFormalIntegral(BenchmarkBase):
    """
    Class to benchmark the numba formal integral function.
    """

    repeat = 2

    @functools.cache
    def setup(self):
        self.sim = self.simulation_verysimple
        self.formal_integral_solver = FormalIntegralSolver(
            self.sim.spectrum_solver.integrator_settings
        )

    # Benchmark for intensity black body function
    def time_intensity_black_body(self):
        nu = 1e14
        temperature = 1e4
        intensity_black_body(nu, temperature)

    # Benchmark for functions in FormalIntegrator class
    def time_FormalIntegrator_functions(self):
        sim_state = self.sim.simulation_state
        transport = self.sim.transport
        plasma = self.sim.plasma
        nu = self.sim.spectrum_solver.spectrum_frequency_grid[:-1]

        self.formal_integral_solver.solve(  # does the work of calculate spectrum and formal_integral
            nu, sim_state, transport, plasma
        )

        atomic_data, opacity_state = self.formal_integral_solver.setup(
            transport, plasma
        )
        source_function_solver = SourceFunctionSolver(
            transport.line_interaction_type
        )
        source_function_state = source_function_solver.solve(
            sim_state, opacity_state, transport.transport_state, atomic_data
        )

        interpolate_shells = (
            self.formal_integral_solver.integrator_settings.interpolate_shells
        )
        (
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            _,
            r_inner_itp,
            r_outer_itp,
            tau_sobolevs_integ,
            electron_densities_integ,
        ) = self.formal_integral_solver.get_interpolated_quantities(
            source_function_state,
            interpolate_shells,
            sim_state,
            transport,
            opacity_state,
            plasma,
        )

        att_S_ul = att_S_ul.flatten(order="F")
        Jred_lu = Jred_lu.flatten(order="F")
        Jblue_lu = Jblue_lu.flatten(order="F")

        self.formal_integral_solver.setup_integrator(
            opacity_state, sim_state.time_explosion, r_inner_itp, r_outer_itp
        )

        self.formal_integral_solver.integrator.formal_integral(
            sim_state.t_inner,
            nu,
            nu.shape[0],
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolevs_integ,
            electron_densities_integ,
            self.formal_integral_solver.integrator_settings.points,
        )
