# cython: profile=False
# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
cimport numpy as np


cdef extern from "math.h":
    double log(double)
    double sqrt(double)
    double abs(double)
    double exp(double)
    
    
ctypedef np.float64_t float_type_t
ctypedef np.int64_t int_type_t


def calculate_beta_sobolev(np.ndarray[float_type_t, ndim=1] tau_sobolevs, np.ndarray[float_type_t, ndim=1] beta_sobolevs):
    cdef float_type_t beta_sobolev
    cdef float_type_t tau_sobolev
    cdef int i
    
    for i in range(len(tau_sobolevs)):
        tau_sobolev = tau_sobolevs[i]
        
        if tau_sobolev > 1e3:
            beta_sobolev = 1 / tau_sobolev
        elif tau_sobolev < 1e-4:
            beta_sobolev = 1 - 0.5 * tau_sobolev
        else:
            beta_sobolev = (1 - exp(-tau_sobolev)) / tau_sobolev
        beta_sobolevs[i] = beta_sobolev
        
        #beta_sobolev = vec_safe_beta_sobolev(tau_sobolev[self.target_line_total])
        
        
def normalize_transition_probabilities(np.ndarray[float_type_t, ndim=1] p_transition,
                                        np.ndarray[int_type_t, ndim=1] reference_levels,
                                        np.ndarray[int_type_t, ndim=1] count_total):
        cdef int_type_t i,j = 0
        cdef float_type_t norm_factor = 0.0
        cdef int_type_t start_id = 0
        cdef int_type_t end_id = 0
        for i in range(len(reference_levels)):
            start_id = reference_levels[i]
            end_id = start_id + count_total[i]
            norm_factor = 0.0
            for j in range(start_id, end_id):
                norm_factor += p_transition[j]
            if norm_factor == 0.0:
                continue
            for j in range(start_id, end_id):
                p_transition[j] = p_transition[j] / norm_factor
            
            
            