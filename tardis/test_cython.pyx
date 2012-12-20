# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: checknone=False


import numpy as np
cimport numpy as np

#cimport packet_source


ctypedef np.float64_t float_type_t
ctypedef np.int64_t int_type_t


def test_list(np.ndarray[float_type_t, ndim=1] packet_input, np.ndarray[float_type_t, ndim=1] output):
    cdef int i
    cdef np.ndarray[float_type_t, ndim=1] input_single
    input_single = packet_input
    for j in range(10000):


        for i in range(len(input_single)):
            output[i] = input_single[i]**2