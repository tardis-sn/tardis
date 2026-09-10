"""
Basic TARDIS Benchmark.
"""

import functools

from asv_runner.benchmarks.mark import parameterize

from benchmarks.benchmark_base import BenchmarkBase


@parameterize({"Input params": ["scatter", "macroatom"]})
class BenchmarkOpacitiesOpacityState(BenchmarkBase):
    """
    Class to benchmark the opacity-state conversion.
    """

    repeat = 2

    @functools.cache
    def setup(self, input_params):
        self.sim = self.nb_simulation_verysimple

    def time_opacity_state_to_numba(self, input_params: str) -> None:
        """Time the opacity state to_numba method"""
        macro_atom_state = (
            None if input_params == "scatter" else self.sim.macro_atom_state
        )
        self.sim.opacity_state.to_numba(macro_atom_state, input_params)
