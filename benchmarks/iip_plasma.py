"""
Benchmarks for IIP plasma helper kernels.
"""

import numpy as np
import pandas as pd

from benchmarks.benchmark_base import BenchmarkBase
from tardis.iip_plasma.continuum.base import ContinuumProcess
from tardis.iip_plasma.properties.continuum import (
    integrate_array_by_level_groups,
)
from tardis.plasma.array_util import (
    cumulative_integrate_array_by_blocks,
)


class BenchmarkIIPPlasma(BenchmarkBase):
    """
    Benchmarks for helper functions used by the IIP plasma workflow.
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
        level_index = pd.MultiIndex.from_arrays(
            [
                np.repeat(np.arange(n_blocks), block_size),
                np.zeros(n_blocks * block_size, dtype=int),
                np.zeros(n_blocks * block_size, dtype=int),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        self.frequency_series = pd.Series(self.frequency, index=level_index)
        transition_blocks = 12000
        transitions_per_block = 4
        source_level_idx = np.repeat(
            np.arange(transition_blocks), transitions_per_block
        )
        destination_level_idx = np.tile(
            np.arange(transitions_per_block), transition_blocks
        )
        index = pd.MultiIndex.from_arrays(
            [source_level_idx, destination_level_idx]
        )
        probabilities = np.random.default_rng(1963).random(
            (transition_blocks * transitions_per_block, n_shells)
        )
        self.transition_probabilities = pd.DataFrame(probabilities, index=index)
        self.transition_probabilities.insert(0, "lines_idx", -1)
        self.transition_probabilities.insert(0, "transition_type", 1)
        cumulative_integrate_array_by_blocks(
            self.alpha_sp_energy,
            self.frequency,
            self.block_references,
        )
        ContinuumProcess._normalize_transition_probabilities(
            self.transition_probabilities,
            no_ref_columns=2,
        )
        integrate_array_by_level_groups(
            self.alpha_sp_energy,
            self.frequency_series,
        )

    def time_cumulative_integrate_array_by_blocks(self):
        cumulative_integrate_array_by_blocks(
            self.alpha_sp_energy,
            self.frequency,
            self.block_references,
        )

    def time_normalize_transition_probabilities(self):
        ContinuumProcess._normalize_transition_probabilities(
            self.transition_probabilities,
            no_ref_columns=2,
        )

    def time_integrate_array_by_level_groups(self):
        integrate_array_by_level_groups(
            self.alpha_sp_energy,
            self.frequency_series,
        )
