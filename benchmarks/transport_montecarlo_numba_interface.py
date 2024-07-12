"""
Basic TARDIS Benchmark.
"""

import numpy as np
from asv_runner.benchmarks.mark import parameterize

from tardis.opacities.opacity_state import opacity_state_initialize
from benchmarks.benchmark_base import BenchmarkBase


class BenchmarkMontecarloMontecarloNumbaNumbaInterface(BenchmarkBase):
    """
    Class to benchmark the numba interface function.
    """

    @parameterize({"Input params": ["scatter", "macroatom", "downbranch"]})
    def time_opacity_state_initialize(self, input_params):
        line_interaction_type = input_params
        plasma = self.nb_simulation_verysimple.plasma
        opacity_state_initialize(
            plasma,
            line_interaction_type,
            self.verysimple_disable_line_scattering,
            self.verysimple_continuum_processes_enabled,
        )