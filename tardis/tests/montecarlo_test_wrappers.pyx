include "../montecarlo.pyx"

def binary_search_wrapper(np.ndarray x, float_type_t x_insert, int_type_t imin, int_type_t imax):
    cdef float_type_t* x_pointer
    x_pointer = <float_type_t*> x.data
    return binary_search(x_pointer, x_insert, imin, imax)

def line_search_wrapper(np.ndarray nu, float_type_t nu_insert,
                        int_type_t number_of_lines):
    cdef float_type_t* nu_pointer
    nu_pointer = <float_type_t*> nu.data
    return line_search(nu_pointer, nu_insert, number_of_lines)