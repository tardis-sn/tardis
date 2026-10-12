"""
Benchmarks for plasma and macro atom helper kernels.
"""

import numpy as np
import pandas as pd

from benchmarks.benchmark_base import BenchmarkBase
from tardis.opacities.macro_atom.macroatom_solver import (
    BoundBoundMacroAtomSolver,
)
from tardis.plasma.array_util import (
    cumulative_integrate_array_by_blocks,
)


class BenchmarkPlasma(BenchmarkBase):
    """
    Benchmarks for helper functions used by the plasma workflow.
    """

    repeat = 3

    def setup(self):
        n_blocks = 720
        block_size = 120
        n_shells = 24
        block_offsets = np.repeat(np.arange(n_blocks) * 2.0, block_size)
        block_grid = np.tile(np.linspace(1.0, 2.0, block_size), n_blocks)
        self.frequency = block_grid + block_offsets
        self.block_references = np.arange(n_blocks + 1) * block_size
        shell_scale = 1000.0 + np.arange(n_shells)
        self.alpha_sp_energy = np.exp(
            -self.frequency[:, np.newaxis] / shell_scale
        )
        transition_blocks = 12000
        transitions_per_block = 4
        probabilities = np.random.default_rng(1963).random(
            (transition_blocks * transitions_per_block, n_shells)
        )
        self.transition_probabilities = pd.DataFrame(probabilities)
        self.transition_metadata = pd.DataFrame(
            {
                "source": [
                    (1, 0, source_level)
                    for source_level in np.repeat(
                        np.arange(transition_blocks), transitions_per_block
                    )
                ]
            }
        )
        cumulative_integrate_array_by_blocks(
            self.alpha_sp_energy,
            self.frequency,
            self.block_references,
        )
        self.normalize_transition_probabilities()

    def time_cumulative_integrate_array_by_blocks(self):
        cumulative_integrate_array_by_blocks(
            self.alpha_sp_energy,
            self.frequency,
            self.block_references,
        )

    def time_normalize_transition_probabilities(self):
        self.normalize_transition_probabilities()

    def normalize_transition_probabilities(self):
        BoundBoundMacroAtomSolver.normalize_transition_probabilities(
            self.transition_probabilities.copy(),
            self.transition_metadata,
        )
