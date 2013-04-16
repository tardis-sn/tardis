# cython: profile=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True


import numpy as np
import logging
import time

from astropy import constants

cimport numpy as np
np.import_array()

from cython.parallel import *


ctypedef np.float64_t float_type_t
ctypedef np.int64_t int_type_t


cdef extern from "math.h":
    float_type_t log(float_type_t)
    float_type_t sqrt(float_type_t)
    float_type_t exp(float_type_t)
    int_type_t floor(float_type_t)
    bint isnan(double x)



cdef extern from "randomkit/randomkit.h":
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
    float_type_t rk_double(rk_state *state)
cdef rk_state mt_state

cdef np.ndarray x
cdef class StorageModel:
    """
    Class for storing the arrays in a cythonized way (as pointers). This ensures fast access during the calculations.
    """

    cdef np.ndarray packet_nus_a
    cdef float_type_t*packet_nus
    cdef np.ndarray packet_mus_a
    cdef float_type_t*packet_mus
    cdef np.ndarray packet_energies_a
    cdef float_type_t*packet_energies
    cdef int_type_t no_of_packets
    cdef int_type_t no_of_shells
    cdef np.ndarray r_inner_a
    cdef float_type_t*r_inner
    cdef np.ndarray r_outer_a
    cdef float_type_t*r_outer
    cdef np.ndarray v_inner_a
    cdef float_type_t*v_inner
    cdef float_type_t time_explosion
    cdef float_type_t inverse_time_explosion
    cdef np.ndarray electron_density_a
    cdef float_type_t*electron_density
    cdef np.ndarray inverse_electron_density_a
    cdef float_type_t*inverse_electron_density
    cdef np.ndarray line_list_nu_a
    cdef float_type_t*line_list_nu
    cdef np.ndarray line_lists_tau_sobolevs_a
    cdef float_type_t*line_lists_tau_sobolevs
    cdef int_type_t line_lists_tau_sobolevs_nd

    #J_BLUES initialize
    cdef np.ndarray line_lists_j_blues_a
    cdef float_type_t*line_lists_j_blues
    cdef int_type_t line_lists_j_blues_nd

    cdef int_type_t no_of_lines
    cdef int_type_t line_interaction_id
    cdef np.ndarray transition_probabilities_a
    cdef float_type_t*transition_probabilities
    cdef int_type_t transition_probabilities_nd
    cdef np.ndarray line2macro_level_upper_a
    cdef int_type_t*line2macro_level_upper
    cdef np.ndarray macro_block_references_a
    cdef int_type_t*macro_block_references
    cdef np.ndarray transition_type_a
    cdef int_type_t*transition_type
    cdef np.ndarray destination_level_id_a
    cdef int_type_t*destination_level_id
    cdef np.ndarray transition_line_id_a
    cdef int_type_t*transition_line_id
    cdef np.ndarray js_a
    cdef float_type_t*js
    cdef np.ndarray nubars_a
    cdef float_type_t*nubars
    cdef float_type_t spectrum_start_nu
    cdef float_type_t spectrum_delta_nu
    cdef float_type_t spectrum_end_nu
    cdef float_type_t*spectrum_virt_nu

    cdef float_type_t sigma_thomson
    cdef float_type_t inverse_sigma_thomson

    def __init__(self, model):

        rk_seed(250819801106, &mt_state)

        cdef np.ndarray[float_type_t, ndim=1] packet_nus = model.packet_src.packet_nus
        self.packet_nus_a = packet_nus
        self.packet_nus = <float_type_t*> self.packet_nus_a.data
        #
        cdef np.ndarray[float_type_t, ndim=1] packet_mus = model.packet_src.packet_mus
        self.packet_mus_a = packet_mus
        self.packet_mus = <float_type_t*> self.packet_mus_a.data
        #
        cdef np.ndarray[float_type_t, ndim=1] packet_energies = model.packet_src.packet_energies
        self.packet_energies_a = packet_energies
        self.packet_energies = <float_type_t*> self.packet_energies_a.data
        #
        self.no_of_packets = packet_nus.size
        #@@@ Setup of Geometry @@@
        self.no_of_shells = model.no_of_shells
        cdef np.ndarray[float_type_t, ndim=1] r_inner = model.r_inner
        self.r_inner_a = r_inner
        self.r_inner = <float_type_t*> self.r_inner_a.data
        #
        cdef np.ndarray[float_type_t, ndim=1] r_outer = model.r_outer
        self.r_outer_a = r_outer
        self.r_outer = <float_type_t*> self.r_outer_a.data
        #
        cdef np.ndarray[float_type_t, ndim=1] v_inner = model.v_inner
        self.v_inner_a = v_inner
        self.v_inner = <float_type_t*> self.v_inner_a.data
        #@@@ Setup the rest @@@
        #times
        self.time_explosion = model.time_explosion
        self.inverse_time_explosion = 1 / model.time_explosion
        #electron density
        cdef np.ndarray[float_type_t, ndim=1] electron_density = model.electron_density
        self.electron_density_a = electron_density
        self.electron_density = <float_type_t*> self.electron_density_a.data
        #
        cdef np.ndarray[float_type_t, ndim=1] inverse_electron_density = 1 / electron_density
        self.inverse_electron_density_a = inverse_electron_density
        self.inverse_electron_density = <float_type_t*> self.inverse_electron_density_a.data
        #Line lists
        cdef np.ndarray[float_type_t, ndim=1] line_list_nu = model.line_list_nu.values
        self.line_list_nu_a = line_list_nu
        self.line_list_nu = <float_type_t*> self.line_list_nu_a.data
        #
        self.no_of_lines = line_list_nu.size

        cdef np.ndarray[float_type_t, ndim=2] line_lists_tau_sobolevs = model.tau_sobolevs
        self.line_lists_tau_sobolevs_a = line_lists_tau_sobolevs
        self.line_lists_tau_sobolevs = <float_type_t*> self.line_lists_tau_sobolevs_a.data
        self.line_lists_tau_sobolevs_nd = self.line_lists_tau_sobolevs_a.shape[1]

        cdef np.ndarray[float_type_t, ndim=2] line_lists_j_blues = model.j_blues
        model.j_blues[:] = 0.0
        self.line_lists_j_blues_a = line_lists_j_blues
        self.line_lists_j_blues = <float_type_t*> self.line_lists_j_blues_a.data
        self.line_lists_j_blues_nd = self.line_lists_j_blues_a.shape[1]

        #
        self.line_interaction_id = model.line_interaction_id
        #macro atom & downbranch
        cdef np.ndarray[float_type_t, ndim=2] transition_probabilities
        cdef np.ndarray[int_type_t, ndim=1] line2macro_level_upper
        cdef np.ndarray[int_type_t, ndim=1] macro_block_references
        cdef np.ndarray[int_type_t, ndim=1] transition_type
        cdef np.ndarray[int_type_t, ndim=1] destination_level_id
        cdef np.ndarray[int_type_t, ndim=1] transition_line_id
        if model.line_interaction_id >= 1:
            transition_probabilities = model.transition_probabilities
            self.transition_probabilities_a = transition_probabilities
            self.transition_probabilities = <float_type_t*> self.transition_probabilities_a.data
            self.transition_probabilities_nd = self.transition_probabilities_a.shape[1]
            #
            line2macro_level_upper = model.atom_data.lines_upper2macro_reference_idx
            self.line2macro_level_upper_a = line2macro_level_upper
            self.line2macro_level_upper = <int_type_t*> self.line2macro_level_upper_a.data
            macro_block_references = model.atom_data.macro_atom_references['block_references'].values
            self.macro_block_references_a = macro_block_references
            self.macro_block_references = <int_type_t*> self.macro_block_references_a.data
            transition_type = model.atom_data.macro_atom_data['transition_type'].values
            self.transition_type_a = transition_type
            self.transition_type = <int_type_t*> self.transition_type_a.data

            #Destination level is not needed and/or generated for downbranch
            destination_level_id = model.atom_data.macro_atom_data['destination_level_idx'].values
            #print "preparing stuff " , destination_level_id[0], destination_level_id[5], destination_level_id[10]
            self.destination_level_id_a = destination_level_id
            #print "preparing stuff(2) " , self.destination_level_id_a[0], self.destination_level_id_a[5], self.destination_level_id_a[10]
            self.destination_level_id = <int_type_t*> self.destination_level_id_a.data

            #print "preparing stuff(3) " , self.destination_level_id[0], self.destination_level_id[5], self.destination_level_id[10]
            transition_line_id = model.atom_data.macro_atom_data['lines_idx'].values
            self.transition_line_id_a = transition_line_id
            self.transition_line_id = <int_type_t*> self.transition_line_id_a.data

        cdef np.ndarray[float_type_t, ndim=1] js = np.zeros(model.no_of_shells, dtype=np.float64)
        cdef np.ndarray[float_type_t, ndim=1] nubars = np.zeros(model.no_of_shells, dtype=np.float64)
        self.js_a = js
        self.js = <float_type_t*> self.js_a.data
        self.nubars_a = nubars
        self.nubars = <float_type_t*> self.nubars_a.data
        self.spectrum_start_nu = model.spec_nu_bins.min()
        self.spectrum_end_nu = model.spec_nu_bins.max()
        self.spectrum_delta_nu = model.spec_nu_bins[1] - model.spec_nu_bins[0]

        cdef np.ndarray[float_type_t, ndim=1] spectrum_virt_nu = model.spec_virtual_flux_nu
        self.spectrum_virt_nu = <float_type_t*> spectrum_virt_nu.data

        if model.sigma_thomson is None:
            self.sigma_thomson = 6.652486e-25 #cm^(-2)
        else:
            self.sigma_thomson = model.sigma_thomson

        self.inverse_sigma_thomson = 1 / self.sigma_thomson
        #
        #
        #

DEF packet_logging = False
IF packet_logging == True:
    packet_logger = logging.getLogger('tardis_packet_logger')
#DEF packet_logging = False

logger = logging.getLogger(__name__)
#Log Level
#cdef int_type_t loglevel = logging.DEBUG








#constants
cdef float_type_t miss_distance = 1e99
cdef float_type_t c = constants.c.cgs.value # cm/s
cdef float_type_t inverse_c = 1 / c
#DEBUG STATEMENT TAKE OUT




cdef int_type_t binary_search(float_type_t*nu, float_type_t nu_insert, int_type_t imin,
                              int_type_t imax):
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
    return imin + 1

#variables are restframe if not specified by prefix comov_
cdef inline int_type_t macro_atom(int_type_t activate_level,
                                  float_type_t*p_transition,
                                  int_type_t p_transition_nd,
                                  int_type_t*type_transition,
                                  int_type_t*target_level_id,
                                  int_type_t*target_line_id,
                                  int_type_t*unroll_reference,
                                  int_type_t cur_zone_id):
    cdef int_type_t emit, i = 0
    cdef float_type_t p, event_random = 0.0
    #print "Activating Level %d" % activate_level
    while True:
        event_random = rk_double(&mt_state)
        #activate_level = 7
        i = unroll_reference[activate_level]
        #i_end = unroll_reference[activate_level + 1]
        #print "checking tprops %.5g" % (p_transition[cur_zone_id * p_transition_nd +i: cur_zone_id * p_transition_nd +i_end)
        #        print "jumping to block_references %d" % i
        p = 0.0
        #cur_zone_id=10

        #print "bunch of numbers %g %g %g %g %g %g %g %g %g" % ( p_transition[i], p_transition[i+1], p_transition[i+2], p_transition[i+3], p_transition[i+4], p_transition[i+5], p_transition[i+6], p_transition[i+7], p_transition[i+8])
        #print "new bunch of numbers %g %g %g %g %g %g %g %g %g" % ( p_transition[cur_zone_id * p_transition_nd + i], p_transition[cur_zone_id * p_transition_nd + i+1], p_transition[cur_zone_id * p_transition_nd + i+2], p_transition[cur_zone_id * p_transition_nd + i+3], p_transition[cur_zone_id * p_transition_nd + i+4], p_transition[cur_zone_id * p_transition_nd + i+5], p_transition[cur_zone_id * p_transition_nd + i+6], p_transition[cur_zone_id * p_transition_nd + i+7], p_transition[cur_zone_id * p_transition_nd + i+8])
        #print "p_transition_nd %g" % p_transition_nd


        while True:
            p = p + p_transition[cur_zone_id * p_transition_nd + i]
            #print " p %g %g %g" % (p, cur_zone_id, i)
            if p > event_random:
                emit = type_transition[i]
                #                print "assuming transition_id %d" % emit
                activate_level = target_level_id[i]
                #print "reActivating Level %g %g" % (activate_level, i)
                break
            i += 1
            if i == unroll_reference[activate_level + 1]:
                logger.warn("overrolling!")
                logger.warn(
                    "activate level %g unroll_reference[activate_level] %g unroll_reference[activate_level+1] %g i=%g event_random=%g current_p=%g" % (
                        activate_level, unroll_reference[activate_level], unroll_reference[activate_level + 1], i,
                        event_random,
                        p))
                time.sleep(10000000)

                #print "I just broke " , emit
        if emit == -1:
            IF packet_logging == True:
                packet_logger.debug('Emitting in level %d', activate_level + 1)

            return target_line_id[i]

cdef float_type_t move_packet(float_type_t*r,
                              float_type_t*mu,
                              float_type_t nu,
                              float_type_t energy,
                              float_type_t distance,
                              float_type_t*js,
                              float_type_t*nubars,
                              float_type_t inverse_t_exp,
                              int_type_t cur_zone_id,
                              int_type_t virtual_packet):
    cdef float_type_t new_r, doppler_factor, comov_energy, comov_nu
    doppler_factor = (1 - (mu[0] * r[0] * inverse_t_exp * inverse_c))
    IF packet_logging == True:
        if distance < 0:
            packet_logger.warn('Distance in move packets less than 0.')

    if distance <= 0:
        return doppler_factor

    new_r = sqrt(r[0] ** 2 + distance ** 2 + 2 * r[0] * distance * mu[0])
    mu[0] = (mu[0] * r[0] + distance) / new_r

    if mu[0] == 0.0:
        print "-------- move packet: mu turned 0.0"
        print distance, new_r, r[0], new_r, cur_zone_id
        print distance / new_r
        #print "move packet after mu", mu[0]
    r[0] = new_r

    if (virtual_packet > 0):
        return doppler_factor

    comov_energy = energy * doppler_factor
    comov_nu = nu * doppler_factor
    js[cur_zone_id] += comov_energy * distance

    nubars[cur_zone_id] += comov_energy * distance * comov_nu

    #print "move packet before mu", mu[0], distance, new_r, r[0]
    #    if distance/new_r > 1e-6:
    #        mu[0] = (distance**2 + new_r**2 - r[0]**2) / (2*distance*new_r)
    #if ((((mu[0] * r[0] + distance) / new_r) < 0.0) and (mu[0] > 0.0)):

    return doppler_factor

cdef void increment_j_blue_estimator(int_type_t*current_line_id, float_type_t*current_nu, float_type_t*current_energy,
                                     float_type_t*mu, float_type_t*r, float_type_t d_line, int_type_t j_blue_idx,
                                     StorageModel storage):
    cdef float_type_t comov_energy, comov_nu, r_interaction, mu_interaction, distance, doppler_factor

    distance = d_line

    r_interaction = sqrt(r[0] ** 2 + distance ** 2 + 2 * r[0] * distance * mu[0])
    mu_interaction = (mu[0] * r[0] + distance) / r_interaction

    doppler_factor = (1 - (mu_interaction * r_interaction * storage.inverse_time_explosion * inverse_c))

    comov_energy = current_energy[0] * doppler_factor
    comov_nu = current_nu[0] * doppler_factor

    storage.line_lists_j_blues[j_blue_idx] += (comov_energy / current_nu[0])
    #print "incrementing j_blues = %g" % storage.line_lists_j_blues[j_blue_idx]

cdef float_type_t compute_distance2outer(float_type_t r, float_type_t  mu, float_type_t r_outer):
    cdef float_type_t d_outer
    d_outer = sqrt(r_outer ** 2 + ((mu ** 2 - 1.) * r ** 2)) - (r * mu)
    return d_outer

cdef float_type_t compute_distance2inner(float_type_t r, float_type_t mu, float_type_t r_inner):
    #compute distance to the inner layer
    #check if intersection is possible?
    cdef float_type_t check, d_inner
    check = r_inner ** 2 + (r ** 2 * (mu ** 2 - 1.))
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
        print "(comov_nu - nu_line) nu_lines", (comov_nu - nu_line) / nu_line
        print "last_line", last_line
        print "next_line", next_line
        print "r", r
        print "mu", mu
        print "nu", nu
        print "doppler_factor", doppler_factor
        print "cur_zone_id", cur_zone_id
        #raise Exception('wrong')

    return ((comov_nu - nu_line) / nu) * c * t_exp

cdef float_type_t compute_distance2electron(float_type_t r, float_type_t mu, float_type_t tau_event,
                                            float_type_t inverse_ne):
    return tau_event * inverse_ne # * inverse_sigma_thomson folded into inverse_ne

cdef inline float_type_t get_r_sobolev(float_type_t r, float_type_t mu, float_type_t d_line):
    return sqrt(r ** 2 + d_line ** 2 + 2 * r * d_line * mu)

def montecarlo_radial1d(model, int_type_t virtual_packet_flag=0):
    """
    Parameters
    ---------


    model : `tardis.model_radial_oned.ModelRadial1D`
        complete model

    param photon_packets : PacketSource object
        photon packets

    Returns
    -------

    output_nus : `numpy.ndarray`

    output_energies : `numpy.ndarray`



    TODO
                    np.ndarray[float_type_t, ndim=1] line_list_nu,
                    np.ndarray[float_type_t, ndim=2] tau_lines,
                    np.ndarray[float_type_t, ndim=1] ne,
                    float_type_t packet_energy,
                    np.ndarray[float_type_t, ndim=2] p_transition,
                    np.ndarray[int_type_t, ndim=1] type_transition,
                    np.ndarray[int_type_t, ndim=1] target_level_id,
                    np.ndarray[int_type_t, ndim=1] target_line_id,
                    np.ndarray[int_type_t, ndim=1] unroll_reference,
                    np.ndarray[int_type_t, ndim=1] line2level,
                    int_type_t log_packets,
                    int_type_t do_scatter

    """

    storage = StorageModel(model)

    ######## Setting up the output ########
    cdef np.ndarray[float_type_t, ndim=1] output_nus = np.zeros(storage.no_of_packets, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] output_energies = np.zeros(storage.no_of_packets, dtype=np.float64)

    ######## Setting up the running variable ########
    cdef float_type_t nu_line = 0.0
    cdef float_type_t current_r = 0.0
    cdef float_type_t current_mu = 0.0
    cdef float_type_t current_nu = 0.0
    cdef float_type_t comov_current_nu = 0.0
    cdef float_type_t current_energy = 0.0

    #indices
    cdef int_type_t current_line_id = 0
    cdef int_type_t current_shell_id = 0

    #Flags for close lines and last line, etc
    cdef int_type_t last_line = 0
    cdef int_type_t close_line = 0
    cdef int_type_t reabsorbed = 0
    cdef int_type_t recently_crossed_boundary = 0
    cdef int i = 0

    for i in range(storage.no_of_packets):
        if i % (storage.no_of_packets / 5) == 0:
            logger.info("At packet %d of %d", i, storage.no_of_packets)


        #setting up the properties of the packet
        current_nu = storage.packet_nus[i]
        current_energy = storage.packet_energies[i]
        current_mu = storage.packet_mus[i]

        #these have been drawn for the comoving frame so we want to convert them        
        comov_current_nu = current_nu

        #Location of the packet
        current_shell_id = 0
        current_r = storage.r_inner[0]
        current_nu = current_nu / (1 - (current_mu * current_r * storage.inverse_time_explosion * inverse_c))
        current_energy = current_energy / (1 - (current_mu * current_r * storage.inverse_time_explosion * inverse_c))

        #linelists
        current_line_id = binary_search(storage.line_list_nu, comov_current_nu, 0, storage.no_of_lines)

        if current_line_id == storage.no_of_lines:
            #setting flag that the packet is off the red end of the line list
            last_line = 1
        else:
            last_line = 0
            #print "storage.no_of_lines %d" % (storage.no_of_lines)
            #print "current_line_id %d" % (current_line_id)
            #print "storage.line_list_nu[current_line_id] %g" % (storage.line_list_nu[current_line_id])
            #print "current_nu %g comov_current_nu %g" % (current_nu, comov_current_nu)
            #return

        #### FLAGS ####
        #Packet recently crossed the inner boundary
        recently_crossed_boundary = 1

        if (virtual_packet_flag > 0):
        #this is a run for which we want the virtual packet spectrum. So first thing we need to do is spawn virtual packets to track the input packet
        #print "I'M STARTING A VIRTUAL FOR A NEW PACKET"
            reabsorbed = montecarlo_one_packet(storage, &current_nu, &current_energy, &current_mu, &current_shell_id,
                                               &current_r, &current_line_id, &last_line, &close_line,
                                               &recently_crossed_boundary, virtual_packet_flag,
                                               -1)

            #print "I'M ABOUT TO START A REAL PACKET"
        #Now can do the propagation of the real packet
        reabsorbed = montecarlo_one_packet(storage, &current_nu, &current_energy, &current_mu, &current_shell_id,
                                           &current_r, &current_line_id, &last_line, &close_line,
                                           &recently_crossed_boundary, virtual_packet_flag, 0)

        if reabsorbed == 1: #reabsorbed
            output_nus[i] = -current_nu
            output_energies[i] = -current_energy

        elif reabsorbed == 0: #emitted
            output_nus[i] = current_nu
            output_energies[i] = current_energy

            #^^^^^^^^^^^^^^^^^^^^^^^^ RESTART MAINLOOP ^^^^^^^^^^^^^^^^^^^^^^^^^

    return output_nus, output_energies, storage.js_a, storage.nubars_a

#
#
#

#Virtual_mode controls whether the packet being propagated is a real packet or a virtual one, to be used to construct a synthetic spectrum. Setting
#this to zero makes it real. Setting it to != 0 (integer) specifies that we want to do virtual particles.
#
#When this routine is called, it is always sent properties of a REAL packet. The issue is whether we are extracting that packet or really propagating it.

cdef int_type_t montecarlo_one_packet(StorageModel storage, float_type_t*current_nu, float_type_t*current_energy,
                                      float_type_t*current_mu, int_type_t*current_shell_id, float_type_t*current_r,
                                      int_type_t*current_line_id, int_type_t*last_line, int_type_t*close_line,
                                      int_type_t*recently_crossed_boundary, int_type_t virtual_packet_flag,
                                      int_type_t virtual_mode):
    cdef int_type_t i
    cdef float_type_t current_nu_virt
    cdef float_type_t current_energy_virt
    cdef float_type_t current_mu_virt
    cdef int_type_t current_shell_id_virt
    cdef float_type_t current_r_virt
    cdef int_type_t current_line_id_virt
    cdef int_type_t last_line_virt
    cdef int_type_t close_line_virt
    cdef int_type_t recently_crossed_boundary_virt
    cdef float_type_t mu_bin
    cdef float_type_t mu_min
    cdef float_type_t doppler_factor_ratio
    cdef float_type_t weight

    if (virtual_mode == 0):
        #print "I'M DOING A REAL PACKET"
        reabsorbed = montecarlo_one_packet_loop(storage, current_nu, current_energy, current_mu, current_shell_id,
                                                current_r, current_line_id, last_line, close_line,
                                                recently_crossed_boundary, virtual_packet_flag, 0)
    else:
        #print "I'M DOING AN EXTRACT"
        for i in range(virtual_packet_flag):
            recently_crossed_boundary_virt = recently_crossed_boundary[0]
            close_line_virt = close_line[0]
            last_line_virt = last_line[0]
            current_line_id_virt = current_line_id[0]
            #comov_current_nu_virt = comov_current_nu[0]
            current_r_virt = current_r[0]
            current_shell_id_virt = current_shell_id[0]

            #choose a direction for the extract packet. We don't want any directions that will hit the inner boundary. So this sets a minimum value for the packet mu
            mu_min = -1. * sqrt(1.0 - ( storage.r_inner[0] / current_r_virt) ** 2)
            mu_bin = (1 - mu_min) / virtual_packet_flag
            current_mu_virt = mu_min + ((i + rk_double(&mt_state)) * mu_bin)

            if (virtual_mode < 0):
                #this is a virtual packet calculation based on a newly born packet - so the weights are more subtle than for a isotropic emission process
                weight = 2. * current_mu_virt / virtual_packet_flag
            else:
                #isotropic emission case ("normal case") for a source in the ejecta
                weight = (1 - mu_min) / 2. / virtual_packet_flag

            #the virtual packets are spawned with known comoving frame energy and frequency

            doppler_factor_ratio = (1 - (current_mu[0] * current_r[0] * storage.inverse_time_explosion * inverse_c)) / (
                1 - (current_mu_virt * current_r_virt * storage.inverse_time_explosion * inverse_c))

            current_energy_virt = current_energy[0] * doppler_factor_ratio
            current_nu_virt = current_nu[0] * doppler_factor_ratio


            #print "DOING i=%g" % i 
            #print "current_nu_virt %g current_energy_virt %g current_mu_virt %g current_shell_id_virt %g current_r_virt %g comov_current_nu_virt %g current_line_id_virt %g last_line_virt %g close_line_virt %g recently_crossed_boundary_virt %g virtual_packet_flag %g " % (current_nu_virt , current_energy_virt , current_mu_virt , current_shell_id_virt , current_r_virt, comov_current_nu_virt, current_line_id_virt, last_line_virt, close_line_virt, recently_crossed_boundary_virt, virtual_packet_flag)
            #print "r_inner %g r_outer %g" % (storage.r_inner[current_shell_id[0]], storage.r_outer[current_shell_id[0]])
            reabsorbed = montecarlo_one_packet_loop(storage, &current_nu_virt, &current_energy_virt, &current_mu_virt,
                                                    &current_shell_id_virt, &current_r_virt, &current_line_id_virt,
                                                    &last_line_virt, &close_line_virt,
                                                    &recently_crossed_boundary_virt, virtual_packet_flag, 1)


            #Putting the virtual nu into the output spectrum
            if (current_nu_virt < storage.spectrum_end_nu) and (current_nu_virt > storage.spectrum_start_nu):
                virt_id_nu = floor(( current_nu_virt - storage.spectrum_start_nu) / storage.spectrum_delta_nu)
                storage.spectrum_virt_nu[virt_id_nu] += current_energy_virt * weight
                #

                #print "I FINISHED DOING A VIRTUAL PACKET"

    return reabsorbed


#
#
cdef int_type_t montecarlo_one_packet_loop(StorageModel storage, float_type_t*current_nu, float_type_t*current_energy,
                                           float_type_t*current_mu, int_type_t*current_shell_id, float_type_t*current_r,
                                           int_type_t*current_line_id, int_type_t*last_line, int_type_t*close_line,
                                           int_type_t*recently_crossed_boundary, int_type_t virtual_packet_flag,
                                           int_type_t virtual_packet):
    cdef float_type_t nu_electron = 0.0
    cdef float_type_t comov_nu = 0.0
    cdef float_type_t comov_energy = 0.0
    cdef float_type_t energy_electron = 0.0
    cdef int_type_t emission_line_id = 0
    cdef int_type_t activate_level_id = 0

    #doppler factor definition
    cdef float_type_t doppler_factor = 0.0
    cdef float_type_t old_doppler_factor = 0.0
    cdef float_type_t inverse_doppler_factor = 0.0

    cdef float_type_t tau_line = 0.0
    cdef float_type_t tau_electron = 0.0
    cdef float_type_t tau_combined = 0.0
    cdef float_type_t tau_event = 0.0


    #defining distances
    cdef float_type_t d_inner = 0.0
    cdef float_type_t d_outer = 0.0
    cdef float_type_t d_line = 0.0
    cdef float_type_t d_electron = 0.0

    cdef int_type_t reabsorbed = 0
    cdef float_type_t nu_line = 0.0

    cdef int_type_t virtual_close_line = 0
    cdef int_type_t j_blue_idx = -1

    #Initializing tau_event if it's a real packet
    if (virtual_packet == 0):
        tau_event = -log(rk_double(&mt_state))

    #For a virtual packet tau_event is the sum of all the tau's that the packet passes.


    #@@@ MAIN LOOP START @@@
    #-----------------------

    while True:
        #check if we are at the end of linelist
        if last_line[0] == 0:
            nu_line = storage.line_list_nu[current_line_id[0]]

        #check if the last line was the same nu as the current line
        if close_line[0] == 1:
            #if yes set the distance to the line to 0.0
            d_line = 0.0
            #reset close_line
            close_line[0] = 0

            #CHECK if 3 lines in a row work
            #print "CLOSE LINE WAS 1"
        else:# -- if close line didn't happen start calculating the the distances
        #print "CLOSE LINE WASN'T 1"
            # ------------------ INNER DISTANCE CALCULATION ---------------------
            if recently_crossed_boundary[0] == 1:
                #if the packet just crossed the inner boundary it will not intersect again unless it interacts. So skip
                #calculation of d_inner
                d_inner = miss_distance
            else:
                #compute distance to the inner shell
                d_inner = compute_distance2inner(current_r[0], current_mu[0], storage.r_inner[current_shell_id[0]])
                # ^^^^^^^^^^^^^^^^^^ INNER DISTANCE CALCULATION ^^^^^^^^^^^^^^^^^^^^^

            # ------------------ OUTER DISTANCE CALCULATION ---------------------
            #computer distance to the outer shell basically always possible
            d_outer = compute_distance2outer(current_r[0], current_mu[0], storage.r_outer[current_shell_id[0]])
            # ^^^^^^^^^^^^^^^^^^ OUTER DISTANCE CALCULATION ^^^^^^^^^^^^^^^^^^^^^

            # ------------------ LINE DISTANCE CALCULATION ---------------------
            if last_line[0] == 1:
                d_line = miss_distance
            else:
                d_line = compute_distance2line(current_r[0], current_mu[0], current_nu[0], nu_line,
                                               storage.time_explosion,
                                               storage.inverse_time_explosion,
                                               storage.line_list_nu[current_line_id[0] - 1],
                                               storage.line_list_nu[current_line_id[0] + 1],
                                               current_shell_id[0])
                # ^^^^^^^^^^^^^^^^^^ LINE DISTANCE CALCULATION ^^^^^^^^^^^^^^^^^^^^^

            # ------------------ ELECTRON DISTANCE CALCULATION ---------------------
            # a virtual packet should never be stopped by continuum processes
            if (virtual_packet > 0):
                d_electron = miss_distance
            else:
                d_electron = compute_distance2electron(current_r[0], current_mu[0], tau_event,
                                                       storage.inverse_electron_density[current_shell_id[0]] * \
                                                       storage.inverse_sigma_thomson)
                # ^^^^^^^^^^^^^^^^^^ ELECTRON DISTANCE CALCULATION ^^^^^^^^^^^^^^^^^^^^^


        # ------------------------------ LOGGING ---------------------- (with precompiler IF)
        IF packet_logging == True:
            packet_logger.debug('%s\nCurrent packet state:\n'
                                'current_mu=%s\n'
                                'current_nu=%s\n'
                                'current_energy=%s\n'
                                'd_inner=%s\n'
                                'd_outer=%s\n'
                                'd_electron=%s\n'
                                'd_line=%s\n%s',
                                '-' * 80,
                                current_mu[0],
                                current_nu[0],
                                current_energy[0],
                                d_inner,
                                d_outer,
                                d_electron,
                                d_line,
                                '-' * 80)

            if isnan(d_inner) or d_inner < 0:
                packet_logger.warning('d_inner is nan or less than 0')
            if isnan(d_outer) or d_outer < 0:
                packet_logger.warning('d_outer is nan or less than 0')
            if isnan(d_electron) or d_electron < 0:
                packet_logger.warning('d_electron is nan or less than 0')
            if isnan(d_line) or d_line < 0:
                packet_logger.warning('d_line is nan or less than 0')

                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^


                #if ((abs(storage.r_inner[current_shell_id[0]] - current_r[0]) < 10.0) and (current_mu[0] < 0.0)):
                #print "d_outer %g, d_inner %g, d_line %g, d_electron %g" % (d_outer, d_inner, d_line, d_electron)
                #print "nu_line %g current_nu[0] %g current_line_id[0] %g" % (nu_line, current_nu[0], current_line_id[0])
                #print "%g " % current_shell_id[0]

        # ------------------------ PROPAGATING OUTWARDS ---------------------------
        if (d_outer <= d_inner) and (d_outer <= d_electron) and (d_outer < d_line):
            #moving one zone outwards. If it's already in the outermost one this is escaped. Otherwise just move, change the zone index
            #and flag as an outwards propagating packet
            move_packet(current_r, current_mu, current_nu[0], current_energy[0], d_outer, storage.js, storage.nubars,
                        storage.inverse_time_explosion,
                        current_shell_id[0], virtual_packet)
            #for a virtual packet, add on the opacity contribution from the continuum
            if (virtual_packet > 0):
                tau_event += (d_outer * storage.electron_density[current_shell_id[0]] * storage.sigma_thomson)
            else:
                tau_event = -log(rk_double(&mt_state))

            if (current_shell_id[0] < storage.no_of_shells - 1): # jump to next shell
                current_shell_id[0] += 1
                recently_crossed_boundary[0] = 1



            else:
                # ------------------------------ LOGGING ---------------------- (with precompiler IF)
                IF packet_logging == True:
                    packet_logger.debug(
                        'Packet has left the simulation through the outer boundary nu=%s mu=%s energy=%s',
                        current_nu[0], current_mu[0], current_energy[0])

                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^

                reabsorbed = 0
                #print "That one got away"
                break
                #
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^ PROPAGATING OUTWARDS ^^^^^^^^^^^^^^^^^^^^^^^^^^


        # ------------------------ PROPAGATING inwards ---------------------------
        elif (d_inner <= d_outer) and (d_inner <= d_electron) and (d_inner < d_line):
            #moving one zone inwards. If it's already in the innermost zone this is a reabsorption
            move_packet(current_r, current_mu, current_nu[0], current_energy[0], d_inner, storage.js, storage.nubars,
                        storage.inverse_time_explosion,
                        current_shell_id[0], virtual_packet)

            #for a virtual packet, add on the opacity contribution from the continuum
            if (virtual_packet > 0):
                tau_event += (d_inner * storage.electron_density[current_shell_id[0]] * storage.sigma_thomson)
            else:
                tau_event = -log(rk_double(&mt_state))

            if current_shell_id[0] > 0:
                current_shell_id[0] -= 1
                recently_crossed_boundary[0] = -1



            else:
                # ------------------------------ LOGGING ---------------------- (with precompiler IF)
                IF packet_logging == True:
                    packet_logger.debug(
                        'Packet has left the simulation through the inner boundary nu=%s mu=%s energy=%s',
                        current_nu[0], current_mu[0], current_energy[0])
                reabsorbed = 1
                break
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^

        # ^^^^^^^^^^^^^^^^^^^^^^^^^^ PROPAGATING INWARDS ^^^^^^^^^^^^^^^^^^^^^^^^^^


        # ------------------------ ELECTRON SCATTER EVENT ELECTRON ---------------------------
        elif (d_electron <= d_outer) and (d_electron <= d_inner) and (d_electron < d_line):
            # we should never enter this branch for a virtual packet
            # ------------------------------ LOGGING ----------------------
            IF packet_logging == True:
                packet_logger.debug('%s\nElectron scattering occuring\n'
                                    'current_nu=%s\n'
                                    'current_mu=%s\n'
                                    'current_energy=%s\n',
                                    '-' * 80,
                                    current_nu[0],
                                    current_mu[0],
                                    current_energy[0])

            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^

            doppler_factor = move_packet(current_r, current_mu, current_nu[0], current_energy[0], d_electron,
                                         storage.js, storage.nubars
                , storage.inverse_time_explosion, current_shell_id[0], virtual_packet)

            comov_nu = current_nu[0] * doppler_factor
            comov_energy = current_energy[0] * doppler_factor

            #new mu chosen
            current_mu[0] = 2 * rk_double(&mt_state) - 1
            inverse_doppler_factor = 1 / (
                1 - (current_mu[0] * current_r[0] * storage.inverse_time_explosion * inverse_c))
            current_nu[0] = comov_nu * inverse_doppler_factor
            current_energy[0] = comov_energy * inverse_doppler_factor
            # ------------------------------ LOGGING ----------------------
            IF packet_logging == True:
                packet_logger.debug('Electron scattering occured\n'
                                    'current_nu=%s\n'
                                    'current_mu=%s\n'
                                    'current_energy=%s\n%s',
                                    current_nu[0],
                                    current_mu[0],
                                    current_energy[0],
                                    '-' * 80)
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^

            tau_event = -log(rk_double(&mt_state))

            #scattered so can re-cross a boundary now
            recently_crossed_boundary[0] = 0
            #We've had an electron scattering event in the SN. This corresponds to a source term - we need to spawn virtual packets now
            if (virtual_packet_flag > 0):
                #print "AN ELECTRON SCATTERING HAPPENED: CALLING VIRTUAL PARTICLES!!!!!!"
                montecarlo_one_packet(storage, current_nu, current_energy, current_mu, current_shell_id, current_r,
                                      current_line_id, last_line, close_line, recently_crossed_boundary,
                                      virtual_packet_flag, 1)

        # ^^^^^^^^^^^^^^^^^^^^^^^^^ SCATTER EVENT LINE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        # ------------------------ LINE SCATTER EVENT  ---------------------------
        elif (d_line <= d_outer) and (d_line <= d_inner) and (d_line <= d_electron):
        #Line scattering
            #It has a chance to hit the line
            if virtual_packet == 0:
                j_blue_idx = current_shell_id[0] * storage.line_lists_j_blues_nd + current_line_id[0]
                increment_j_blue_estimator(current_line_id, current_nu, current_energy, current_mu, current_r, d_line,
                                           j_blue_idx, storage)

            tau_line = storage.line_lists_tau_sobolevs[
                current_shell_id[0] * storage.line_lists_tau_sobolevs_nd + current_line_id[0]]

            tau_electron = storage.sigma_thomson * storage.electron_density[current_shell_id[0]] * d_line
            tau_combined = tau_line + tau_electron

            # ------------------------------ LOGGING ----------------------
            IF packet_logging == True:
                packet_logger.debug('%s\nEntering line scattering routine\n'
                                    'Scattering at line %d (nu=%s)\n'
                                    'tau_line=%s\n'
                                    'tau_electron=%s\n'
                                    'tau_combined=%s\n'
                                    'tau_event=%s\n',
                                    '-' * 80,
                                    current_line_id[0] + 1,
                                    storage.line_list_nu[current_line_id[0]],
                                    tau_line,
                                    tau_electron,
                                    tau_combined,
                                    tau_event)
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^

            #Advancing to next line
            current_line_id[0] += 1

            #check for last line
            if current_line_id[0] >= storage.no_of_lines:
                current_line_id[0] = storage.no_of_lines
                last_line[0] = 1

            if (virtual_packet > 0):
                #here we should never stop the packet - only account for its tau increment from the line
                tau_event += tau_line

            else:
                #Check for line interaction
                if tau_event < tau_combined:
                #line event happens - move and scatter packet
                #choose new mu
                # ------------------------------ LOGGING ----------------------
                #IF packet_logging == True:
                #      packet_logger.debug('Line interaction happening. Activating macro atom at level %d',
                #          line2macro_level_upper[current_line_id] + 1)
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^



                    old_doppler_factor = move_packet(current_r, current_mu, current_nu[0], current_energy[0], d_line,
                                                     storage.js,
                                                     storage.nubars, storage.inverse_time_explosion,
                                                     current_shell_id[0], virtual_packet)

                    current_mu[0] = 2 * rk_double(&mt_state) - 1

                    inverse_doppler_factor = 1 / (
                        1 - (current_mu[0] * current_r[0] * storage.inverse_time_explosion * inverse_c))

                    comov_nu = current_nu[0] * old_doppler_factor
                    comov_energy = current_energy[0] * old_doppler_factor

                    #new mu chosen
                    current_energy[0] = comov_energy * inverse_doppler_factor



                    #here comes the macro atom

                    #print "A LINE EVENT HAPPENED AT R = %g in shell %g which has range %g to %g" % (current_r[0], current_shell_id[0], storage.r_inner[current_shell_id[0]], storage.r_outer[current_shell_id[0]])
                    #print "d_outer %g, d_inner %g, d_line %g, d_electron %g" % (d_outer, d_inner, d_line, d_electron)
                    #print "nu_line_1 %g nu_line_2 %g" % (storage.line_list_nu[current_line_id[0]-1], storage.line_list_nu[current_line_id[0]])
                    #print "nu_line_1 %g nu_line_2 %g" % (current_line_id[0]-1, current_line_id[0])
                    #print "last_line %g" % (last_line[0])




                    if storage.line_interaction_id == 0: #scatter
                        emission_line_id = current_line_id[0] - 1
                    elif storage.line_interaction_id >= 1:# downbranch & macro
                        activate_level_id = storage.line2macro_level_upper[current_line_id[0] - 1]
                        #print "HERE %g %g" %(current_line_id[0], activate_level_id)
                        #print "DEST " , (storage.destination_level_id[0], storage.destination_level_id[5], storage.destination_level_id[10])
                        #print "DEST " , (storage.macro_block_references[0], storage.macro_block_references[5], storage.macro_block_references[10])
                        emission_line_id = macro_atom(activate_level_id,
                                                      storage.transition_probabilities,
                                                      storage.transition_probabilities_nd,
                                                      storage.transition_type,
                                                      storage.destination_level_id,
                                                      storage.transition_line_id,
                                                      storage.macro_block_references,
                                                      current_shell_id[0])

                    current_nu[0] = storage.line_list_nu[emission_line_id] * inverse_doppler_factor
                    nu_line = storage.line_list_nu[emission_line_id]
                    current_line_id[0] = emission_line_id + 1

                    IF packet_logging == True:
                        packet_logger.debug('Line interaction over. New Line %d (nu=%s; rest)', emission_line_id + 1,
                                            storage.line_list_nu[emission_line_id])




                    # getting new tau_event
                    tau_event = -log(rk_double(&mt_state))

                    # reseting recently crossed boundary - can intersect with inner boundary again
                    recently_crossed_boundary[0] = 0



                    #We've had a line event in the SN. This corresponds to a source term - we need to spawn virtual packets now
                    if (virtual_packet_flag > 0):
                        #print "LINE EVENT HAPPENED: TRIGGERING A VIRTUAL PACKET CALCULATION %g!!!!!!" % current_line_id[0]
                        #print "nu_line %g current_nu[0] %g storage.line_list_nu[emission_line_id] %g storage.line_list_nu[emission_line_id-1] %g emission_line_id %g" % (nu_line, current_nu[0],  storage.line_list_nu[emission_line_id] , storage.line_list_nu[emission_line_id-1] , emission_line_id)

                        virtual_close_line = 0
                        if last_line[
                            0] == 0: #Next line is basically the same just making sure we take this into account
                            if abs(storage.line_list_nu[current_line_id[0]] - nu_line) / nu_line < 1e-7:
                                virtual_close_line = 1

                        montecarlo_one_packet(storage, current_nu, current_energy, current_mu, current_shell_id,
                                              current_r, current_line_id, last_line, &virtual_close_line,
                                              recently_crossed_boundary,
                                              virtual_packet_flag, 1)
                        virtual_close_line = 0

                else: #tau_event > tau_line no interaction so far
                    #reducing event tau to take this probability into account
                    tau_event -= tau_line

                    # ------------------------------ LOGGING ----------------------
                    IF packet_logging == True:
                        packet_logger.debug('No line interaction happened. Tau_event decreasing %s\n%s', tau_event,
                                            '-' * 80)
                        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^

            # ------------------------------ LOGGING ----------------------
            IF packet_logging == True: #THIS SHOULD NEVER HAPPEN
                if tau_event < 0:
                    logging.warn('tau_event less than 0: %s', tau_event)
                    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^ LOGGING # ^^^^^^^^^^^^^^^^^^^^^^^^^^^

            if last_line[0] == 0: #Next line is basically the same just making sure we take this into account
                if abs(storage.line_list_nu[current_line_id[0]] - nu_line) / nu_line < 1e-7:
                    close_line[0] = 1
                    #print "CLOSE LINE: line_list_nu[current_line_id[0]] %g nu_line %g abs(storage.line_list_nu[current_line_id[0]] - nu_line) %g abs(storage.line_list_nu[current_line_id[0]] - nu_line) / nu_line %g" % (storage.line_list_nu[current_line_id[0]], nu_line, abs(storage.line_list_nu[current_line_id[0]] - nu_line), abs(storage.line_list_nu[current_line_id[0]] - nu_line) / nu_line)
                    #print "here %g %g %g" % (sqrt(2.), abs(-2.), log(2.7))
                    # ^^^^^^^^^^^^^^^^^^^^^^^^^ SCATTER EVENT LINE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        if (virtual_packet > 0):
            if (tau_event > 10.0):
                tau_event = 100.0
                reabsorbed = 0
                break


    # ------------------------------ LOGGING ----------------------
    IF packet_logging == True: #SHOULD NEVER HAPPEN
        if current_energy[0] < 0:
            logging.warn('Outcoming energy less than 0: %s', current_energy)
            # ^^^^^^^^^^^^^^^^^^^^^^^^^ SCATTER EVENT LINE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    if (virtual_packet > 0):
        current_energy[0] = current_energy[0] * exp(-1. * tau_event)


        #if (reabsorbed == 1):
        #print "reabsorbed mu = %g" % (current_mu[0])

    if (current_energy[0] < 0):
        print "negative energy"

    return reabsorbed
