"""
Basic TARDIS Benchmark.
"""

from benchmarks.benchmark_base import BenchmarkBase


class BenchmarkMacroAtomSolver(BenchmarkBase):
    """
    Class to benchmark the macro_atom_solver.
    """

    repeat = 2

    def setup(self):
        self.sim = self.nb_simulation_verysimple
        self.macro_atom = self.sim.macro_atom
        self.plasma = self.sim.plasma

    def time_macro_atom_solve(self):
        """Time the macro atom transition probability computation."""
        self.macro_atom.solve(
            self.plasma.j_blues,
            self.plasma.beta_sobolev,
            self.plasma.stimulated_emission_factor,
        )
