from tardis.montecarlo import montecarlo
import os
from ctypes import CDLL, byref, c_int64, c_double, c_uint
clib = CDLL(os.path.join(montecarlo.__file__))


class montecarlo_test_wrappers:

        @staticmethod
        def binary_search_wrapper(x, x_insert, imin, imax):

                x = (c_double * (imax - imin + 1))(*x)
                x_insert = c_double(x_insert)
                imin = c_int64(imin)
                imax = c_int64(imax)
                obtained_result = c_int64(0)
                clib.binary_search.restype = c_uint
                clib.binary_search(byref(x), x_insert, imin, imax, byref(obtained_result))
