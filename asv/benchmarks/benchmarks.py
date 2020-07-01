# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.

import numpy as np
from tardis.tests import montecarlo_test_wrappers as montecarlo

LINE_SIZE = 10000000


class TimeSuite:
    """
    An example benchmark that times the performance of various kinds
    of iterating over dictionaries in Python.
    """

    def setup(self):
        self.line = np.arange(LINE_SIZE, 1, -1).astype(np.float64)

    def time_binarysearch(self):
        for _ in range(LINE_SIZE):
            montecarlo.binary_search_wrapper(
                self.line, np.random.random() * LINE_SIZE, 0, LINE_SIZE - 1
            )

    def time_compute_distance2outer(self):
        for _ in range(1000000):
            montecarlo.compute_distance2outer_wrapper(0.0, 0.5, 1.0)
            montecarlo.compute_distance2outer_wrapper(1.0, 0.5, 1.0)
            montecarlo.compute_distance2outer_wrapper(0.3, 1.0, 1.0)
            montecarlo.compute_distance2outer_wrapper(0.3, -1.0, 1.0)
            montecarlo.compute_distance2outer_wrapper(0.5, 0.0, 1.0)

    def time_compute_distance2inner(self):
        for _ in range(1000000):
            montecarlo.compute_distance2inner_wrapper(1.5, -1.0, 1.0)
            montecarlo.compute_distance2inner_wrapper(0.0, 0.0, 0.0)
            montecarlo.compute_distance2inner_wrapper(1.2, -0.7, 1.0)

    def time_compute_distance2line(self):
        for _ in range(1000000):
            montecarlo.compute_distance2line_wrapper(
                2.20866912e15,
                -0.251699059004,
                1.05581082105e15,
                1.06020910733e15,
                1693440.0,
                5.90513983371e-07,
                1.0602263591e15,
                1.06011723237e15,
                2,
            )
            montecarlo.compute_distance2line_wrapper(
                2.23434667994e15,
                -0.291130548401,
                1.05581082105e15,
                1.06733618121e15,
                1693440.0,
                5.90513983371e-07,
                1.06738407486e15,
                1.06732933961e15,
                3,
            )

    def time_compute_distance2electron(self):
        for _ in range(1000000):
            montecarlo.compute_distance2electron_wrapper(0.0, 0.0, 2.0, 2.0)
