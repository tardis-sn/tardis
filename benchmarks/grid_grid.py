"""
Basic TARDIS Benchmark.
"""
import numpy as np
import pandas as pd
from asv_runner.benchmarks.mark import skip_benchmark

import tardis.grid as grid
from benchmarks.benchmark_base import BenchmarkBase


# @skip_benchmark
class BenchmarkGridGrid(BenchmarkBase):
    """
    Class to benchmark the grid function.
    """

    def __init__(self):
        self.DATA_PATH = self.get_absolute_path('tardis/grid/tests/data')

    def time_grid(self):
        dfpath = f"{self.DATA_PATH}/example_grid.txt"
        ymlpath = f"{self.DATA_PATH}/example.yml"
        axesdict = {
            "model.structure.velocity.start": np.arange(10000, 15000, 1000),
            "model.abundances.He": np.arange(0, 1, 0.1),
            "model.abundances.H": np.arange(0, 1, 0.25),
        }

        df = pd.read_csv(dfpath)
        g = grid.tardisGrid(configFile=ymlpath, gridFrame=df)
        grid.tardisGrid.from_axes(configFile=ymlpath, axesdict=axesdict)

        # Check that grid attribute has the right shape
        ax_len = 1
        for key in axesdict:
            ax_len *= len(axesdict[key])

        g.grid_row_to_config(row_index=0)

        # Verify that a model can be returned.
        g.grid_row_to_simulation_state(
            row_index=0, atomic_data=self.atomic_dataset
        )
