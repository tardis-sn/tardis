"""
Basic TARDIS Benchmark.
"""
import pandas as pd
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.model.matter.decay import IsotopicMassFraction


# @skip_benchmark
class BenchmarkIoDecay(BenchmarkBase):
    """
    Class to benchmark the decay function.
    """

    def __init__(self):
        pass

    @property
    def simple_abundance_model(self):
        index = pd.MultiIndex.from_tuples(
            [(28, 56)], names=["atomic_number", "mass_number"]
        )
        return IsotopicMassFraction([[1.0, 1.0]], index=index)

    @property
    def raw_abundance_simple(self):
        abundances = pd.DataFrame([[0.2, 0.2], [0.1, 0.1]], index=[28, 30])
        abundances.index.rename("atomic_number", inplace=True)
        return abundances

    def time_simple_decay(self):
        self.simple_abundance_model.decay(100)

    def time_abundance_merge(self):
        decayed_df = self.simple_abundance_model.decay(100)
        decayed_df.as_atoms()
        decayed_df.merge(self.raw_abundance_simple, normalize=False)
