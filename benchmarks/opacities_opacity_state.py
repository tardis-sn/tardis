"""
Basic TARDIS Benchmark.
"""

from asv_runner.benchmarks.mark import parameterize

from tardis.opacities.opacity_state import opacity_state_initialize
from benchmarks.benchmark_base import BenchmarkBase
import functools

@parameterize({"Input params": ["scatter", "macroatom"]})
class BenchmarkOpacitiesOpacityState(BenchmarkBase):
    """
    Class to benchmark the numba interface function.
    """

    @functools.cache
    def setup(self, input_params):
        self.sim = self.nb_simulation_verysimple

    def time_opacity_state_initialize(self, input_params):
        line_interaction_type = input_params
        plasma = self.sim.plasma
        opacity_state_initialize(
            plasma,
            line_interaction_type,
            True
        )
