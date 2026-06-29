"""
Benchmarks for IIP plasma helper kernels.
"""

import numpy as np
import pandas as pd

from benchmarks.benchmark_base import BenchmarkBase
from tardis.iip_plasma.continuum.base import ContinuumProcess
from tardis.iip_plasma.properties.continuum import (
    integrate_array_by_prepared_level_groups,
    prepare_level_group_integration,
)
from tardis.plasma.properties.continuum_processes.fast_array_util import (
    cumulative_integrate_array_by_blocks,
)

N_BLOCKS = 720
BLOCK_SIZE = 120
N_SHELLS = 24
N_TRANSITION_BLOCKS = 12000
N_TRANSITIONS_PER_BLOCK = 4
RNG_SEED = 1963


class BenchmarkIIPPlasma(BenchmarkBase):
    """
    Benchmarks for helper functions used by the IIP plasma workflow.
    """

    repeat = 3

    def setup(self):
        self._setup_frequency_grid()
        self._setup_transition_probabilities()

        cumulative_integrate_array_by_blocks(
            self.alpha_sp_energy,
            self.frequency,
            self.block_references,
        )
        ContinuumProcess._normalize_transition_probabilities(
            self.transition_probabilities,
            no_ref_columns=2,
        )
        integrate_array_by_prepared_level_groups(
            self.alpha_sp_energy[:, 0],
            *self.integration_data,
        )

    def _setup_frequency_grid(self):
        block_offsets = np.repeat(
            np.arange(N_BLOCKS) * 2.0,
            BLOCK_SIZE,
        )
        block_grid = np.tile(np.linspace(1.0, 2.0, BLOCK_SIZE), N_BLOCKS)
        self.frequency = block_grid + block_offsets
        self.block_references = np.arange(N_BLOCKS + 1) * BLOCK_SIZE

        shell_scale = 1000.0 + np.arange(N_SHELLS)
        self.alpha_sp_energy = np.exp(
            -self.frequency[:, np.newaxis] / shell_scale
        )

        level_index = pd.MultiIndex.from_arrays(
            [
                np.repeat(np.arange(N_BLOCKS), BLOCK_SIZE),
                np.zeros(N_BLOCKS * BLOCK_SIZE, dtype=int),
                np.zeros(N_BLOCKS * BLOCK_SIZE, dtype=int),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        frequency_series = pd.Series(self.frequency, index=level_index)
        self.integration_data = prepare_level_group_integration(
            frequency_series
        )

    def _setup_transition_probabilities(self):
        source_level_idx = np.repeat(
            np.arange(N_TRANSITION_BLOCKS),
            N_TRANSITIONS_PER_BLOCK,
        )
        destination_level_idx = np.tile(
            np.arange(N_TRANSITIONS_PER_BLOCK),
            N_TRANSITION_BLOCKS,
        )
        index = pd.MultiIndex.from_arrays(
            [source_level_idx, destination_level_idx]
        )

        probabilities = np.random.default_rng(RNG_SEED).random(
            (N_TRANSITION_BLOCKS * N_TRANSITIONS_PER_BLOCK, N_SHELLS)
        )
        self.transition_probabilities = pd.DataFrame(
            probabilities,
            index=index,
        )
        self.transition_probabilities.insert(0, "lines_idx", -1)
        self.transition_probabilities.insert(0, "transition_type", 1)

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

    def time_integrate_array_by_prepared_level_groups(self):
        integrate_array_by_prepared_level_groups(
            self.alpha_sp_energy[:, 0],
            *self.integration_data,
        )
