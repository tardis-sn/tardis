include "../montecarlo.pyx"

def get_miss_distance():
    return miss_distance

def binary_search_wrapper(np.ndarray x, float_type_t x_insert, int_type_t imin, int_type_t imax):
    cdef float_type_t* x_pointer
    x_pointer = <float_type_t*> x.data
    return binary_search(x_pointer, x_insert, imin, imax)

def line_search_wrapper(np.ndarray nu, float_type_t nu_insert,
                        int_type_t number_of_lines):
    cdef float_type_t* nu_pointer
    nu_pointer = <float_type_t*> nu.data
    return line_search(nu_pointer, nu_insert, number_of_lines)

def compute_distance2outer_wrapper(float_type_t r, float_type_t mu, float_type_t r_outer):
    return compute_distance2outer(r, mu, r_outer)

def compute_distance2inner_wrapper(float_type_t r, float_type_t mu, float_type_t r_inner):
    return compute_distance2inner(r, mu, r_inner)

def compute_distance2line_wrapper(float_type_t r, float_type_t mu, float_type_t nu, float_type_t nu_line, float_type_t t_exp, float_type_t inverse_t_exp, float_type_t last_line, float_type_t next_line, int_type_t cur_zone_id):
    return compute_distance2line(r, mu, nu, nu_line, t_exp, inverse_t_exp, last_line, next_line, cur_zone_id)

def compute_distance2electron_wrapper(float_type_t r, float_type_t mu, float_type_t tau_event, float_type_t inverse_ne):
    return compute_distance2electron(r, mu, tau_event, inverse_ne)