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
            montecarlo.binary_search_wrapper(self.line, np.random.random() * LINE_SIZE, 0, LINE_SIZE - 1)
