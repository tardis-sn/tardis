"""
Basic TARDIS Benchmark.
"""
from asv_runner.benchmarks.mark import parameterize, skip_benchmark
from numpy.testing import assert_almost_equal

import tardis.montecarlo.montecarlo_numba.nonhomologous_grid as nonhomologous_grid
from benchmarks.benchmark_base import BenchmarkBase


# @skip_benchmark
class BenchmarkMontecarlo_nonhomologous(BenchmarkBase):
    """
    Class to benchmark the nonhomologous function.
    """

    def __init__(self):
        pass

    @parameterize({"Parameters": [
        {
            "a": 0.0,
            "b": 0.0,
            "c": 0.0,
            "d": 2.0,
            "e": -1.0,
            "threshold": 0.0,
            "expected_roots": [0.5],
        },
        {
            "a": 1.0,
            "b": 2.0,
            "c": 0.0,
            "d": 2.0,
            "e": 0.0,
            "threshold": 0.0,
            "expected_roots": [],
        },
        {
            "a": 1.0,
            "b": -14.0,
            "c": 71.0,
            "d": -154.0,
            "e": 120.0,
            "threshold": 2.5,
            "expected_roots": [3, 4, 5],
        },
    ]})
    def time_quartic_roots(self, parameters):
        a = parameters['a']
        b = parameters['b']
        c = parameters['c']
        d = parameters['d']
        e = parameters['e']
        threshold = parameters['threshold']
        expected_roots = parameters['expected_roots']
        obtained_roots = nonhomologous_grid.quartic_roots(a, b, c, d, e, threshold)
        obtained_roots.sort()

        assert_almost_equal(obtained_roots, expected_roots)
