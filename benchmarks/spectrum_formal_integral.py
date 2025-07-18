"""
Basic TARDIS Benchmark.
"""

import functools

from numba import config

from benchmarks.benchmark_base import BenchmarkBase
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
from tardis.spectrum.formal_integral.formal_integral import FormalIntegrator
from tardis.spectrum.formal_integral.formal_integral_solver import FormalIntegralSolver
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
        self.FormalIntegrator = FormalIntegralSolver(
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

        self.FormalIntegrator.solve( # does the work of calculate spectrum and formal_integral
            nu,
            sim_state,
            transport,
            plasma
        )

        atomic_data, levels, opacity_state = self.FormalIntegrator.setup(transport, plasma)
        sourceFunction = SourceFunctionSolver(transport.line_interaction_type, atomic_data)
        res = sourceFunction.solve(
            sim_state, 
            opacity_state,  
            transport.transport_state, 
            levels)
        
        interpolate_shells = self.integrator_settings.interpolate_shells
        att_S_ul, Jred_lu, Jblue_lu, _, r_inner_i, r_outer_i, tau_sobolevs_integ, electron_densities_integ = self.get_interpolated_quantities(
            res, interpolate_shells, sim_state, transport, opacity_state, plasma
        ) 

        # self.FormalIntegrator.generate_numba_objects()
        # self.FormalIntegrator.formal_integral(
        #     self.sim.spectrum_solver.spectrum_frequency_grid, 1000
        # )
        self.FormalIntegrator.setup_integrator(opacity_state, sim_state.time_explosion, r_inner_i, r_outer_i)
        self.FormalIntegrator.integrator.formal_integral(
            sim_state.t_inner,
            nu,
            nu.shape[0],
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolevs_integ,
            electron_densities_integ,
            self.FormalIntegrator.integrator_settings.points
        )

