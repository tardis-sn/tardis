
# cython: profile=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True


import numpy as np
cimport numpy as np
np.import_array()
from cython.parallel import prange

cdef extern from "math.h":
    double log(double)
    double sqrt(double)
    double abs(double)
    bint isnan(double x)

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
ctypedef np.int64_t int_type_t

cdef rk_state mt_state
rk_seed(250819801106, &mt_state)




#constants
cdef float_type_t miss_distance = 1e99
cdef float_type_t c = 2.99792458e10 # cm/s
cdef float_type_t inverse_c = 1 / c
cdef float_type_t sigma_thomson = 6.652486e-25 #cm^(-2)

cdef float_type_t inverse_sigma_thomson = 1 / sigma_thomson



cdef int_type_t binary_search(np.ndarray[float_type_t, ndim=1] nu, float_type_t nu_insert, int_type_t imin, int_type_t imax):
    #continually narrow search until just one element remains
    cdef int_type_t imid
    while imax - imin > 2:
        imid = (imin + imax) / 2
        
        #code must guarantee the interval is reduced at each iteration
        #assert(imid < imax);
        # note: 0 <= imin < imax implies imid will always be less than imax
   
        # reduce the search
        if (nu[imid] < nu_insert):
            imax = imid + 1
        else:
            imin = imid
        #print imin, imax, imid, imax - imin
    return imin+1
    
#variables are restframe if not specified by prefix comov_
cdef int_type_t macro_atom(int_type_t activate_level,
                            np.ndarray[float_type_t, ndim=2] p_transition,
                            np.ndarray[int_type_t, ndim=1] type_transition,
                            np.ndarray[int_type_t, ndim=1] target_level_id,
                            np.ndarray[int_type_t, ndim=1] target_line_id,
                            np.ndarray[int_type_t, ndim=1] unroll_reference,
			    int_type_t cur_zone_id):
    cdef int_type_t emit, i = 0
    cdef float_type_t p = 0.0
    
    while True:
        event_random = rk_double(&mt_state)
        i = unroll_reference[activate_level]
        p = 0.0
        while True:
            p = p + p_transition[cur_zone_id, i]
            if p > event_random:
                emit = type_transition[i]
                activate_level = target_level_id[i]
                break
            i += 1
        if emit == 1:
            #print "emitted with target_line_id %d" % target_line_id[i]
            return target_line_id[i]
        
        

cdef float_type_t move_packet(float_type_t* r,
			float_type_t* mu,
			float_type_t nu,
			float_type_t energy,
			float_type_t distance,
			np.ndarray[float_type_t, ndim=1] js,
			np.ndarray[float_type_t, ndim=1] nubars,
			float_type_t inverse_t_exp,
			int_type_t cur_zone_id):
    
    cdef float_type_t new_r, doppler_factor, comov_energy, comov_nu
    doppler_factor = (1 - (mu[0] * r[0] * inverse_t_exp * inverse_c))
    
    if distance <= 0:
        return doppler_factor
    
    comov_energy = energy * doppler_factor
    comov_nu = nu * doppler_factor
    js[cur_zone_id] += comov_energy * distance
    
    nubars[cur_zone_id] += comov_energy * distance * comov_nu 
    
    new_r = sqrt(r[0]**2 + distance**2 + 2 * r[0] * distance * mu[0])
    #print "move packet before mu", mu[0], distance, new_r, r[0]
#    if distance/new_r > 1e-6:
#        mu[0] = (distance**2 + new_r**2 - r[0]**2) / (2*distance*new_r)
    mu[0] = (mu[0] * r[0] + distance) / new_r
        
    if mu[0] == 0.0:
        print "-------- move packet: mu turned 0.0"
        print distance, new_r, r[0], new_r, cur_zone_id
        print distance/new_r
    #print "move packet after mu", mu[0]
    r[0] = new_r
    return doppler_factor
    

cdef float_type_t compute_distance2outer(float_type_t r, float_type_t  mu, float_type_t r_outer):
       d_outer = sqrt(r_outer**2 + ((mu**2 - 1.) * r**2)) - (r * mu)
       return d_outer

cdef float_type_t compute_distance2inner(float_type_t r, float_type_t mu, float_type_t r_inner):
    #compute distance to the inner layer
    #check if intersection is possible?
    cdef float_type_t check, d_inner
    check = r_inner**2 + (r**2 * (mu**2 - 1.))
    if check < 0:
        return miss_distance
    else:
        if mu < 0:
           d_inner = -r * mu - sqrt(check)
           return d_inner
        else:
            return miss_distance

cdef float_type_t compute_distance2line(float_type_t r, float_type_t mu,
                                    float_type_t nu, float_type_t nu_line,
                                    float_type_t t_exp, float_type_t inverse_t_exp,
                                    float_type_t last_line, float_type_t next_line, int_type_t cur_zone_id):
        #computing distance to line
        cdef float_type_t comov_nu, doppler_factor
        doppler_factor = (1. - (mu * r * inverse_t_exp * inverse_c))
        comov_nu = nu * doppler_factor

        if comov_nu < nu_line:
            #TODO raise exception
            print "WARNING comoving nu less than nu_line shouldn't happen:"
            print "comov_nu = ", comov_nu
            print "nu_line", nu_line
            print "(comov_nu - nu_line) nu_lines", (comov_nu - nu_line)/nu_line
            print "last_line", last_line
            print "next_line", next_line
            print "r", r
            print "mu", mu
            print "nu", nu
            print "doppler_factor", doppler_factor
            print "cur_zone_id", cur_zone_id
            #raise Exception('wrong')

        return ((comov_nu - nu_line) / nu) * c * t_exp


cdef float_type_t compute_distance2electron(float_type_t r, float_type_t mu, float_type_t tau_event, float_type_t inverse_ne):
    return tau_event * inverse_ne * inverse_sigma_thomson



cdef float_type_t get_r_sobolev(float_type_t r, float_type_t mu, float_type_t d_line):
    return sqrt(r**2 + d_line**2 + 2 * r * d_line * mu)


def run_simple_oned(np.ndarray[float_type_t, ndim=1] packets,
                np.ndarray[float_type_t, ndim=1] mus,
                np.ndarray[float_type_t, ndim=1] line_list_nu,
                np.ndarray[float_type_t, ndim=2] tau_lines,
                np.ndarray[float_type_t, ndim=1] r_inner,
                np.ndarray[float_type_t, ndim=1] r_outer,
                np.ndarray[float_type_t, ndim=1] v_inner,
                np.ndarray[float_type_t, ndim=1] ne,
                float_type_t packet_energy,
                np.ndarray[float_type_t, ndim=2] p_transition,
                np.ndarray[int_type_t, ndim=1] type_transition,
                np.ndarray[int_type_t, ndim=1] target_level_id,
                np.ndarray[int_type_t, ndim=1] target_line_id,
                np.ndarray[int_type_t, ndim=1] unroll_reference,
                np.ndarray[int_type_t, ndim=1] line2level
                ):
    
    cdef int_type_t no_of_zones = len(r_inner)
    cdef float_type_t t_exp = r_inner[0] / v_inner[0]
    cdef float_type_t inverse_t_exp = 1 / t_exp
    cdef np.ndarray[float_type_t, ndim=1] inverse_ne = 1 / ne
    
    cdef int no_of_packets = packets.size
    cdef int no_of_lines = line_list_nu.size

    #outputs
    cdef np.ndarray[float_type_t, ndim=1] nus = np.zeros(no_of_packets, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] energies = np.zeros(no_of_packets, dtype=np.float64)
    
    cdef np.ndarray[float_type_t, ndim=1] nus_reabsorbed = np.zeros(no_of_packets, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] energies_reabsorbed = np.zeros(no_of_packets, dtype=np.float64)
    
    
    cdef np.ndarray[float_type_t, ndim=1] js = np.zeros(no_of_zones, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] nubars = np.zeros(no_of_zones, dtype=np.float64)
    #cdef np.ndarray[float_type_t, ndim=1] js = np.zeros(1, dtype=np.float64)
    #cdef np.ndarray[float_type_t, ndim=1] nubars = np.zeros(1, dtype=np.float64)
    
    cdef float_type_t nu_line = 0.0
    cdef float_type_t nu_electron = 0.0
    cdef float_type_t current_r = 0.0
    cdef float_type_t current_mu = 0.0
    cdef float_type_t current_nu = 0.0
    cdef float_type_t comov_current_nu = 0.0
    cdef float_type_t comov_nu = 0.0
    cdef float_type_t comov_energy = 0.0
    cdef float_type_t comov_current_energy = 0.0
    cdef float_type_t current_energy = 0.0
    cdef float_type_t energy_electron = 0.0
    cdef int_type_t emission_line_id = 0

    #doppler factor definition
    cdef float_type_t doppler_factor = 0.0
    cdef float_type_t old_doppler_factor = 0.0
    cdef float_type_t inverse_doppler_factor = 0.0
    
    cdef float_type_t tau_line = 0.0
    cdef float_type_t tau_electron = 0.0
    cdef float_type_t tau_combined = 0.0
    cdef float_type_t tau_event = 0.0
    #indices
    cdef int_type_t current_line_id = 0
    cdef int_type_t current_zone_id = 0
    cdef int_type_t current_line_list_id = 0
    #defining distances
    cdef float_type_t d_inner = 0.0
    cdef float_type_t d_outer = 0.0
    cdef float_type_t d_line = 0.0
    cdef float_type_t d_electron = 0.0
    
    #Flags for close lines and last line, etc
    cdef int_type_t last_line = 0
    cdef int_type_t close_line = 0
    cdef int_type_t reabsorbed = 0
    cdef int_type_t recently_crossed_boundary = 0    

    cdef int i=0
    
    for i in range(no_of_packets):
        
        if i % 1000 == 0: print "@packet %d" % i
        
        current_nu = packets[i]
        current_energy = packet_energy
        
        current_mu = mus[i]
        #print "initial_mu" , current_mu
        current_r = r_inner[0]
        current_zone_id = 0

        recently_crossed_boundary = 1
 
        tau_event = -log(rk_double(&mt_state))
        

        
        comov_current_nu = current_nu * (1 - (current_mu * current_r * inverse_t_exp * inverse_c))
        
        #cur_line_id = line_list_nu.size - line_list_nu[::-1].searchsorted(comov_current_nu)

        current_line_id = binary_search(line_list_nu, comov_current_nu, 0, no_of_lines)
        current_line_list_id = current_line_id
        #print "cur_line_id binary_search %d" % test_line_id
        if current_line_id == line_list_nu.size: last_line=1
        else: last_line = 0

        # ---- main loop stars
        while True:
                

            #check if we are at the end of linelist
            if last_line == 0:
                nu_line = line_list_nu[current_line_id]
            
            
            if close_line == 1:
                d_line = 0.0
                close_line = 0
                
                #CHECK if 3 lines in a row work
    
            else:
                if recently_crossed_boundary == 1:
		    #if the packet just crossed the inner boundary it will not intersect again unless it interacts. So skip
		    #calculation of d_inner
                    d_inner = miss_distance
                else:
                    d_inner = compute_distance2inner(current_r, current_mu, r_inner[current_zone_id])
                    
                d_outer = compute_distance2outer(current_r, current_mu, r_outer[current_zone_id])
                if last_line == 1:
                    d_line = miss_distance
                else:
                    d_line = compute_distance2line(current_r, current_mu, current_nu, nu_line, t_exp, inverse_t_exp, line_list_nu[current_line_id-1], line_list_nu[current_line_id+1], current_zone_id)
                d_electron = compute_distance2electron(current_r, current_mu, tau_event, inverse_ne[current_zone_id])
                
# ------------------------------------------------------------------------------------------            
            if isnan(d_outer):
                print '-------- d_outer nan'
                print current_mu, d_inner, d_outer, current_zone_id
                print current_r, current_mu, r_inner[current_zone_id], r_outer[current_zone_id]
                print recently_crossed_boundary
                return 0
            if current_mu == 0.0:
                print '-------- mu==0.0'
                print current_mu, d_inner, d_outer, current_zone_id
                print current_r, current_mu, r_inner[current_zone_id], r_outer[current_zone_id]
                print recently_crossed_boundary
                return 0
#            print current_mu, d_inner, d_outer, cur_zone_id
            # ---- calculating distances #    
# ------------------------------------------------------------------------------------------            
            if (d_outer < d_inner) and (d_outer < d_electron) and (d_outer < d_line):
	        #moving one zone outwards. If it's already in the outermost one this is escaped. Otherwise just move, change the zone index
                #and flag as an outwards propagating packet
                move_packet(&current_r, &current_mu, current_nu, current_energy, d_outer, js, nubars, inverse_t_exp, current_zone_id)
                if (current_zone_id < no_of_zones - 1):
                    current_zone_id += 1
                    recently_crossed_boundary = 1
#                    print "----zone change outwards"
#                    print cur_zone_id, current_r, r_inner[cur_zone_id], r_outer[cur_zone_id], current_mu
                else:
                    #escaped
                    reabsorbed = 0
                    #print "That one got away"
                    break
                    #      
            elif (d_inner < d_outer) and (d_inner < d_electron) and (d_inner < d_line):
    	        #moving one zone inwards. If it's already in the innermost zone this is a reabsorption
                move_packet(&current_r, &current_mu, current_nu, current_energy, d_inner, js, nubars, inverse_t_exp, current_zone_id)
                if current_zone_id > 0:
                    current_zone_id -= 1
#                    print "----zone change inwards"
#                    print cur_zone_id, current_r, r_inner[cur_zone_id], r_outer[cur_zone_id], current_mu
                    recently_crossed_boundary = -1
                else:
                    #reabsorbed
                    reabsorbed = 1
                    #print "another one bites the dust"
                    break
            
            elif (d_electron < d_outer) and (d_electron < d_inner) and (d_electron < d_line):
            #electron scattering
#                print "electron scattering happened"
                doppler_factor = move_packet(&current_r, &current_mu, current_nu, current_energy, d_electron, js, nubars, inverse_t_exp, current_zone_id)
                
                
                comov_nu = current_nu * doppler_factor
                comov_energy = current_energy * doppler_factor
                
                #new mu chosen
                current_mu = 2*rk_double(&mt_state) - 1
                inverse_doppler_factor = 1/(1 - (current_mu * current_r * inverse_t_exp * inverse_c))
                current_nu = comov_nu * inverse_doppler_factor
                current_energy = comov_energy * inverse_doppler_factor
                #currrent_energy = 0.0
                tau_event = -log(rk_double(&mt_state))
                #scattered so can re-cross a boundary now
                recently_crossed_boundary = 0
            
            elif (d_line < d_outer) and (d_line < d_inner) and (d_line < d_electron):
            #Line scattering
                #It has a chance to hit the line
                tau_line = tau_lines[current_zone_id, current_line_id]
                tau_electron = sigma_thomson * ne[current_zone_id] * d_line
                tau_combined = tau_line + tau_electron
                prev_r = current_r
                

                
                
                current_line_id += 1
                
                #check for last line
                if current_line_id >= no_of_lines:
                    current_line_id = no_of_lines
                    last_line = 1
                
                #check for same line        
                    
                #Check for line interaction
                if tau_event < tau_combined:
#                    print "line event happened"
                    #line event happens - move and scatter packet
                    #choose new mu
                    old_doppler_factor = move_packet(&current_r, &current_mu, current_nu, current_energy, d_line, js, nubars, inverse_t_exp, current_zone_id)
                    comov_current_energy = current_energy * old_doppler_factor
                    
                    current_mu = 2*rk_double(&mt_state) - 1
                    
                    inverse_doppler_factor = 1 / (1 - (current_mu * current_r * inverse_t_exp * inverse_c))
                    
                    #here comes the macro atom
                    activate_level_id = line2level[current_line_id]
                    emission_line_id = macro_atom(activate_level_id,
                                                p_transition,
                                                type_transition,
                                                target_level_id,
                                                target_line_id,
                                                unroll_reference,
		    				 current_zone_id)
                    #emission_line_id = cur_line_id - 1
                    current_nu = line_list_nu[emission_line_id] * inverse_doppler_factor
                    nu_line = line_list_nu[emission_line_id]
                    current_line_id = emission_line_id + 1
                    ### end of macro_atom
                    
                    current_energy = comov_current_energy * inverse_doppler_factor
                    
                    #print "current_energy after scattering (energy, nu, nu_line, mu)"
                    #print current_energy, current_nu, nu_line, current_mu
                    #print '-----'
                    tau_event = -log(rk_double(&mt_state))
        		    #it has scattered and so can now recross a boundary
                    recently_crossed_boundary = 0
                else:
                    tau_event -= tau_line
                if tau_event < 0:
                    print 'ola, what happened here'
                
                if last_line == 0:
                    if abs(line_list_nu[current_line_id] - nu_line)/nu_line < 1e-7:
                        close_line = 1

        if current_energy < 0:
            #TODO logging warning
            print "current_energy less than 0"
        if reabsorbed == 1:
        #TODO bin them right away
            nus_reabsorbed[i] = current_nu
            energies_reabsorbed[i] = current_energy
            
        elif reabsorbed == 0:
            nus[i] = current_nu
            #print "emitted packet energy", current_energy
            energies[i] = current_energy

    return nus, energies, nus_reabsorbed, energies_reabsorbed, js, nubars