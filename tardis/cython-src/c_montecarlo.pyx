# cython: profile=False
# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
from cython.parallel import prange
cimport numpy as np

cdef extern from "math.h":
    double log(double)
    double sqrt(double)
    double abs(double)

cdef extern from "../randomkit/randomkit.h":
    ctypedef struct rk_state:
        unsigned long key[624]
        int pos
        int has_gauss
        double gauss

    ctypedef enum rk_error:
        RK_NOERR = 0
        RK_ENODEV = 1
        RK_ERR_MAX = 2
    
    void rk_seed(unsigned long seed, rk_state *state)
    double rk_double(rk_state *state)

    

ctypedef np.float64_t float_type_t

cdef rk_state mt_state
rk_seed(250819801106, &mt_state)


#np.random.seed(250819801106)

#constants
cdef float_type_t miss_distance = 1e99
cdef float_type_t c = 2.99792458e10 # cm/s
cdef float_type_t inverse_c = 1 / c
cdef float_type_t sigma_thomson = 6.652486e-25 #cm^(-2)

cdef float_type_t inverse_sigma_thomson = 1 / sigma_thomson

#variables are restframe if not specified by prefix comov_


cdef float_type_t move_packet(float_type_t* r, float_type_t* mu, float_type_t nu, float_type_t energy, float_type_t distance, float_type_t* js, float_type_t* nubars, float_type_t inverse_t_exp):
    
    cdef float_type_t new_r, doppler_factor, comov_energy, comov_nu
    
    doppler_factor = (1 - (mu[0] * r[0] * inverse_t_exp * inverse_c))
    comov_energy = energy * doppler_factor
    comov_nu = nu * doppler_factor
    js[0] += comov_energy * distance
    
    nubars[0] += comov_energy * distance * comov_nu 
    
    new_r = sqrt(r[0]**2 + distance**2 + 2 * r[0] * distance * mu[0])
    mu[0] = (distance**2 + new_r**2 - r[0]**2) / (2*distance*new_r)
    r[0] = new_r
    return doppler_factor
    

cdef float_type_t compute_distance2outer(float_type_t r, float_type_t  mu, float_type_t r_outer):
    return sqrt(r_outer**2 + ((mu**2 - 1.) * r**2)) - (r * mu)


cdef float_type_t compute_distance2inner(float_type_t r, float_type_t mu, float_type_t r_inner):
    #compute distance to the inner layer
    #check if intersection is possible?
    cdef float_type_t check
    
    check = r_inner**2 + (r**2 * (mu**2 - 1.))
    if check < 0:
        return miss_distance
    else:
        if mu < 0:
           return -r * mu - sqrt(check)
        else:
            return miss_distance

cdef float_type_t compute_distance2line(float_type_t r, float_type_t mu,
                                    float_type_t nu, float_type_t nu_line,
                                    float_type_t t_exp, float_type_t inverse_t_exp):
        #computing distance to line
        cdef float_type_t comov_nu, doppler_factor
        doppler_factor = (1. - (mu * r * inverse_t_exp * inverse_c))
        comov_nu = nu * doppler_factor

        if comov_nu < nu_line:
            #TODO raise exception
            print "WARNING comoving nu less than nu_line shouldn't happen", comov_nu, nu_line, (comov_nu - nu_line)/nu_line
        else:
            return ((comov_nu - nu_line) / nu) * c * t_exp


cdef float_type_t compute_distance2electron(float_type_t r, float_type_t mu, float_type_t tau_event, float_type_t inverse_ne):
    return tau_event * inverse_ne * inverse_sigma_thomson



cdef float_type_t get_r_sobolev(float_type_t r, float_type_t mu, float_type_t d_line):
    return sqrt(r**2 + d_line**2 + 2 * r * d_line * mu)


def run_simple_oned(np.ndarray[float_type_t, ndim=1] packets,
                np.ndarray[float_type_t, ndim=1] mus,
                np.ndarray[float_type_t, ndim=1] line_list_nu,
                np.ndarray[float_type_t, ndim=1] tau_lines,
                float_type_t r_inner,
                float_type_t r_outer,
                float_type_t v_inner,
                float_type_t ne,
                float_type_t packet_energy):
    
    cdef float_type_t test_random=0.
    
    cdef float_type_t t_exp = r_inner / v_inner
    cdef float_type_t inverse_t_exp = 1 / t_exp
    cdef float_type_t inverse_ne = 1 / ne
    
    cdef int no_of_packets = packets.shape[0]
    cdef int no_of_lines = line_list_nu.shape[0]

    #outputs
    cdef np.ndarray[float_type_t, ndim=1] nus = np.zeros(no_of_packets, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] energies = np.zeros(no_of_packets, dtype=np.float64)
    
    cdef np.ndarray[float_type_t, ndim=1] nus_reabsorbed = np.zeros(no_of_packets, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] energies_reabsorbed = np.zeros(no_of_packets, dtype=np.float64)
    
    
    cdef float_type_t js = 0
    cdef float_type_t nubars = 0
    #cdef np.ndarray[float_type_t, ndim=1] js = np.zeros(1, dtype=np.float64)
    #cdef np.ndarray[float_type_t, ndim=1] nubars = np.zeros(1, dtype=np.float64)
    
    cdef float_type_t nu_line = 0.0
    cdef float_type_t nu_electron = 0.0
    cdef float_type_t current_r = 0.0
    cdef float_type_t prev_r = 0.0
    cdef float_type_t current_mu = 0.0
    cdef float_type_t current_nu = 0.0
    cdef float_type_t comov_current_nu = 0.0
    cdef float_type_t comov_nu = 0.0
    cdef float_type_t comov_energy = 0.0
    cdef float_type_t comov_current_energy = 0.0
    cdef float_type_t current_energy = 0.0
    cdef float_type_t energy_electron = 0.0
    
    #doppler factor definition
    cdef float_type_t doppler_factor = 0.0
    cdef float_type_t old_doppler_factor = 0.0
    cdef float_type_t inverse_doppler_factor = 0.0
    
    cdef float_type_t tau_line = 0.0
    cdef float_type_t tau_electron = 0.0
    cdef float_type_t tau_combined = 0.0
    cdef float_type_t tau_event = 0.0
    #indices
    cdef int cur_line_id
    
    #defining distances
    cdef float_type_t d_inner = 0.0
    cdef float_type_t d_outer = 0.0
    cdef float_type_t d_line = 0.0
    cdef float_type_t d_electron = 0.0
    
    #Flags for close lines and last line, etc
    cdef int last_line = 0
    cdef int close_line = 0
    cdef int reabsorbed = 0
    
    cdef int i=0
    
    for i in range(no_of_packets):
        
        if i % 1000 == 0: print "@packet %d" % i
        
        current_nu = packets[i]
        current_energy = packet_energy
        
        current_mu = mus[i]
        current_r = r_inner
        
        tau_event = -log(rk_double(&mt_state))
        
        comov_current_nu = current_nu * (1. - (current_mu * v_inner * inverse_c))
        
        cur_line_id = line_list_nu.size - line_list_nu[::-1].searchsorted(comov_current_nu)
        
        if cur_line_id == line_list_nu.size: last_line=1
        else: last_line = 0
        
        #numerical problem: when the particles are sitting on the r_inner and one calculatese d_inner numerical instabilities can make it negative (depending on mu)
        #Solution we give the packets a little nudge (choosing r_inner*1e-8 as moving distance)
        
        move_packet(&current_r, &current_mu, current_nu, current_energy, r_inner*1e-8, &js, &nubars, inverse_t_exp)
        
        while True:
            #check if we are at the end of linelist
            if last_line == 0:
                nu_line  = line_list_nu[cur_line_id]
            
            
            if close_line == 1:
                d_line = 0.0
                close_line = 0
                
                #CHECK if 3 lines in a row work
                
            else:
                d_inner = compute_distance2inner(current_r, current_mu, r_inner)
                d_outer = compute_distance2outer(current_r, current_mu, r_outer)
                if last_line == 1:
                    d_line = 1e99
                else:
                    d_line = compute_distance2line(current_r, current_mu, current_nu, nu_line, t_exp, inverse_t_exp)
                d_electron = compute_distance2electron(current_r, current_mu, tau_event, inverse_ne)
                
                if current_r <0:
                    print " aha. r=",current_r, 
                if d_line < 0:
                    print "d-line negative"
                if d_electron < 0:
                    print "d-electron negative"
                if d_inner < 0:
                    print "d-inner negative", current_r, r_inner, d_inner, d_inner/current_r
                
                if d_outer < 0:
                    print "d-outer negative"
            if (d_outer < d_inner) and (d_outer < d_electron) and (d_outer < d_line):
                #escaped
                reabsorbed = 0
                #print "That one got away"
                move_packet(&current_r, &current_mu, current_nu, current_energy, d_outer, &js, &nubars, inverse_t_exp)
                break
            
            #packet reabsorbing into core
            
            elif (d_inner < d_outer) and (d_inner < d_electron) and (d_inner < d_line):
                #reabsorbed
                reabsorbed = 1
                #print "another one bites the dust"
                move_packet(&current_r, &current_mu, current_nu, current_energy, d_inner, &js, &nubars, inverse_t_exp)
                break
            
            elif (d_electron < d_outer) and (d_electron < d_inner) and (d_electron < d_line):
            #electron scattering
                doppler_factor = move_packet(&current_r, &current_mu, current_nu, current_energy, d_electron, &js, &nubars, inverse_t_exp)
                
                
                comov_nu = current_nu * doppler_factor
                comov_energy = current_energy * doppler_factor
                
                #new mu chosen
                current_mu = 2*rk_double(&mt_state) - 1
                inverse_doppler_factor = 1/(1 - (current_mu * current_r * inverse_t_exp * inverse_c))
                current_nu = comov_nu * inverse_doppler_factor
                current_energy = comov_energy * inverse_doppler_factor
                tau_event = -log(rk_double(&mt_state))
            
            elif (d_line < d_outer) and (d_line < d_inner) and (d_line < d_electron):
            #Line scattering
                #It has a chance to hit the line
                tau_line = tau_lines[cur_line_id]
                tau_electron = sigma_thomson * ne * d_line
                tau_combined = tau_line + tau_electron
                prev_r = current_r
                

                
                
                cur_line_id += 1
                
                #check for last line
                if cur_line_id >= no_of_lines:
                    cur_line_id = no_of_lines
                    last_line = 1
                
                #check for same line        
                if last_line == 0:
                    if abs(line_list_nu[cur_line_id] - nu_line)/nu_line < 1e-7:
                        close_line = 1
                    
                #Check for line interaction
                if tau_event < tau_combined:
                    #line event happens - move and scatter packet
                    #choose new mu
                    old_doppler_factor = move_packet(&current_r, &current_mu, current_nu, current_energy, d_line, &js, &nubars, inverse_t_exp)
                    comov_current_energy = current_energy * old_doppler_factor
                    
                    current_mu = 2*rk_double(&mt_state) - 1
                    inverse_doppler_factor = 1 / (1 - (current_mu * current_r * inverse_t_exp * inverse_c))
                    current_nu = nu_line * inverse_doppler_factor
                    current_energy = comov_current_energy * inverse_doppler_factor
                    #current_energy = 0.0
                    tau_event = -log(rk_double(&mt_state))
                else:
                    tau_event -= tau_line
                if tau_event < 0:
                    print 'ola, what happened here'
            
        if current_energy < 0:
            1/0
        if reabsorbed == 1:
        #TODO bin them right away
            nus_reabsorbed[i] = current_nu
            energies_reabsorbed[i] = current_energy
            
        elif reabsorbed == 0:
            nus[i] = current_nu
            energies[i] = current_energy
            
        else:
            1/0

    return nus, energies, nus_reabsorbed, energies_reabsorbed, js, nubars