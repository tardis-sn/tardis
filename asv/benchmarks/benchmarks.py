# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
from tardis.util import base
import numpy as np

class Suite:
    def time_number2element(self):
        for i in np.arange(1,100):
            base.atomic_number2element_symbol(i)




