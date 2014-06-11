# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True


import logging
import time

import numpy as np
cimport numpy as np
from astropy import constants

np.import_array()

ctypedef np.float64_t float_type_t
ctypedef np.int64_t int_type_t

ctypedef struct rpacket_t:
    float_type_t nu
    float_type_t mu
    float_type_t energy
    float_type_t r
    int_type_t current_shell_id
    int_type_t next_line_id
    int_type_t last_line
    int_type_t close_line
    int_type_t recently_crossed_boundary
    int_type_t virtual_packet_flag

ctypedef struct storage_model_t:
    float_type_t *packet_nus
    float_type_t *packet_mus
    float_type_t *packet_energies
    float_type_t *output_nus
    float_type_t *output_energies
    int_type_t *last_line_interaction_in_id
    int_type_t *last_line_interaction_out_id
    int_type_t *last_line_interaction_shell_id
    int_type_t *last_interaction_type
    int_type_t no_of_packets
    int_type_t no_of_shells
    float_type_t *r_inner
    float_type_t *r_outer
    float_type_t *v_inner
    float_type_t time_explosion
    float_type_t inverse_time_explosion
    float_type_t *electron_densities
    float_type_t *inverse_electron_densities
    float_type_t *line_list_nu
    float_type_t *line_lists_tau_sobolevs
    int_type_t line_lists_tau_sobolevs_nd
    float_type_t *line_lists_j_blues
    int_type_t line_lists_j_blues_nd
    int_type_t no_of_lines
    int_type_t line_interaction_id
    float_type_t *transition_probabilities
    int_type_t transition_probabilities_nd
    int_type_t *line2macro_level_upper
    int_type_t *macro_block_references
    int_type_t *transition_type
    int_type_t *destination_level_id
    int_type_t *transition_line_id
    float_type_t *js
    float_type_t *nubars
    float_type_t spectrum_start_nu
    float_type_t spectrum_delta_nu
    float_type_t spectrum_end_nu
    float_type_t *spectrum_virt_nu
    float_type_t sigma_thomson
    float_type_t inverse_sigma_thomson
    float_type_t inner_boundary_albedo
    int_type_t reflective_inner_boundary
    int_type_t current_packet_id


cdef extern int_type_t line_search(float_type_t *nu, float_type_t nu_insert, int_type_t number_of_lines) except -1
cdef extern int_type_t binary_search(float_type_t *x, float_type_t x_insert, int_type_t imin, int_type_t imax) except -1
cdef extern float_type_t compute_distance2outer(float_type_t r, float_type_t mu, float_type_t r_outer)
cdef extern float_type_t compute_distance2inner(float_type_t r, float_type_t mu, float_type_t r_inner)
cdef extern float_type_t compute_distance2line(float_type_t r, float_type_t mu, float_type_t nu, float_type_t nu_line, float_type_t t_exp, float_type_t inverse_t_exp, float_type_t last_line, float_type_t next_line, int_type_t cur_zone_id) except? 0
cdef extern float_type_t compute_distance2electron(float_type_t r, float_type_t mu, float_type_t tau_event, float_type_t inverse_ne)
cdef extern int_type_t macro_atom(int_type_t activate_level, float_type_t *p_transition, int_type_t p_transition_nd, int_type_t *type_transition, int_type_t *target_level_id, int_type_t *target_line_id, int_type_t *unroll_reference, int_type_t cur_zone_id)
cdef extern float_type_t move_packet(float_type_t *r, float_type_t *mu, float_type_t nu, float_type_t energy, float_type_t distance, float_type_t *js, float_type_t *nubars, float_type_t inverse_t_exp, int_type_t cur_zone_id, int_type_t virtual_packet)
cdef extern void increment_j_blue_estimator(int_type_t *current_line_id, float_type_t *current_nu, float_type_t *current_energy, float_type_t *mu, float_type_t *r, float_type_t d_line, int_type_t j_blue_idx, float_type_t inverse_time_explosion, float_type_t *line_lists_j_blues)
cdef extern int_type_t montecarlo_one_packet(storage_model_t *storage, float_type_t *current_nu, float_type_t *current_energy, float_type_t *current_mu, int_type_t *current_shell_id, float_type_t *current_r, int_type_t *current_line_id, int_type_t *last_line, int_type_t *close_line, int_type_t *recently_crossed_boundary, int_type_t virtual_packet_flag, int_type_t virtual_mode)
cdef extern int_type_t montecarlo_one_packet(storage_model_t *storage, float_type_t *current_nu, float_type_t *current_energy, float_type_t *current_mu, int_type_t *current_shell_id, float_type_t *current_r, int_type_t *current_line_id, int_type_t *last_line, int_type_t *close_line, int_type_t *recently_crossed_boundary, int_type_t virtual_packet_flag, int_type_t virtual_mode)
cdef extern void rpacket_init(rpacket_t *packet, storage_model_t *storage, float_type_t nu, float_type_t mu, float_type_t energy, int_type_t virtual_packet)

cdef extern from "math.h":
    float_type_t log(float_type_t)
    float_type_t sqrt(float_type_t)
    float_type_t exp(float_type_t)
    int_type_t floor(float_type_t)
    bint isnan(double x)

cdef extern from "randomkit.h":
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

cdef extern rk_state mt_state

#constants
cdef float_type_t miss_distance = 1e99
cdef float_type_t c = constants.c.cgs.value # cm/s
cdef float_type_t inverse_c = 1 / c
#DEBUG STATEMENT TAKE OUT

def montecarlo_radial1d(model, int_type_t virtual_packet_flag=0):
    """
    Parameters
    ----------
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
    cdef storage_model_t storage
    cdef rpacket_t packet
    rk_seed(model.tardis_config.montecarlo.seed, &mt_state)
    cdef np.ndarray[float_type_t, ndim=1] packet_nus = model.packet_src.packet_nus
    storage.packet_nus = <float_type_t*> packet_nus.data
    cdef np.ndarray[float_type_t, ndim=1] packet_mus = model.packet_src.packet_mus
    storage.packet_mus = <float_type_t*> packet_mus.data
    cdef np.ndarray[float_type_t, ndim=1] packet_energies = model.packet_src.packet_energies
    storage.packet_energies = <float_type_t*> packet_energies.data
    storage.no_of_packets = packet_nus.size
    # Setup of structure
    structure = model.tardis_config.structure
    storage.no_of_shells = structure.no_of_shells
    cdef np.ndarray[float_type_t, ndim=1] r_inner = structure.r_inner.to('cm').value
    storage.r_inner = <float_type_t*> r_inner.data
    cdef np.ndarray[float_type_t, ndim=1] r_outer = structure.r_outer.to('cm').value
    storage.r_outer = <float_type_t*> r_outer.data
    cdef np.ndarray[float_type_t, ndim=1] v_inner = structure.v_inner.to('cm/s').value
    storage.v_inner = <float_type_t*> v_inner.data
    # Setup the rest
    # times
    storage.time_explosion = model.tardis_config.supernova.time_explosion.to('s').value
    storage.inverse_time_explosion = 1.0 / storage.time_explosion
    #electron density
    cdef np.ndarray[float_type_t, ndim=1] electron_densities = model.plasma_array.electron_densities.values
    storage.electron_densities = <float_type_t*> electron_densities.data
    cdef np.ndarray[float_type_t, ndim=1] inverse_electron_densities = 1.0 / electron_densities
    storage.inverse_electron_densities = <float_type_t*> inverse_electron_densities.data
    # Line lists
    cdef np.ndarray[float_type_t, ndim=1] line_list_nu = model.atom_data.lines.nu.values
    storage.line_list_nu = <float_type_t*> line_list_nu.data
    storage.no_of_lines = line_list_nu.size
    cdef np.ndarray[float_type_t, ndim=2] line_lists_tau_sobolevs = model.plasma_array.tau_sobolevs.values.transpose()
    storage.line_lists_tau_sobolevs = <float_type_t*> line_lists_tau_sobolevs.data
    storage.line_lists_tau_sobolevs_nd = line_lists_tau_sobolevs.shape[1]
    cdef np.ndarray[float_type_t, ndim=2] line_lists_j_blues = model.j_blue_estimators
    storage.line_lists_j_blues = <float_type_t*> line_lists_j_blues.data
    storage.line_lists_j_blues_nd = line_lists_j_blues.shape[1]
    line_interaction_type = model.tardis_config.plasma.line_interaction_type
    if line_interaction_type == 'scatter':
        storage.line_interaction_id = 0
    elif line_interaction_type == 'downbranch':
        storage.line_interaction_id = 1
    elif line_interaction_type == 'macroatom':
        storage.line_interaction_id = 2
    else:
        storage.line_interaction_id = -99
    # macro atom & downbranch
    cdef np.ndarray[float_type_t, ndim=2] transition_probabilities
    cdef np.ndarray[int_type_t, ndim=1] line2macro_level_upper
    cdef np.ndarray[int_type_t, ndim=1] macro_block_references
    cdef np.ndarray[int_type_t, ndim=1] transition_type
    cdef np.ndarray[int_type_t, ndim=1] destination_level_id
    cdef np.ndarray[int_type_t, ndim=1] transition_line_id
    if storage.line_interaction_id >= 1:
        transition_probabilities = model.transition_probabilities.values.transpose()
        storage.transition_probabilities = <float_type_t*> transition_probabilities.data
        storage.transition_probabilities_nd = transition_probabilities.shape[1]
        line2macro_level_upper = model.atom_data.lines_upper2macro_reference_idx
        storage.line2macro_level_upper = <int_type_t*> line2macro_level_upper.data
        macro_block_references = model.atom_data.macro_atom_references['block_references'].values
        storage.macro_block_references = <int_type_t*> macro_block_references.data
        transition_type = model.atom_data.macro_atom_data['transition_type'].values
        storage.transition_type = <int_type_t*> transition_type.data
        # Destination level is not needed and/or generated for downbranch
        destination_level_id = model.atom_data.macro_atom_data['destination_level_idx'].values
        storage.destination_level_id = <int_type_t*> destination_level_id.data
        transition_line_id = model.atom_data.macro_atom_data['lines_idx'].values
        storage.transition_line_id = <int_type_t*> transition_line_id.data
    cdef np.ndarray[float_type_t, ndim=1] output_nus = np.zeros(storage.no_of_packets, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] output_energies = np.zeros(storage.no_of_packets, dtype=np.float64)
    storage.output_nus = <float_type_t*> output_nus.data
    storage.output_energies = <float_type_t*> output_energies.data
    cdef np.ndarray[int_type_t, ndim=1] last_line_interaction_in_id = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] last_line_interaction_out_id = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] last_line_interaction_shell_id = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] last_interaction_type = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    storage.last_line_interaction_in_id = <int_type_t*> last_line_interaction_in_id.data
    storage.last_line_interaction_out_id = <int_type_t*> last_line_interaction_out_id.data
    storage.last_line_interaction_shell_id = <int_type_t*> last_line_interaction_shell_id.data
    storage.last_interaction_type = <int_type_t*> last_interaction_type.data
    cdef np.ndarray[float_type_t, ndim=1] js = np.zeros(storage.no_of_shells, dtype=np.float64)
    cdef np.ndarray[float_type_t, ndim=1] nubars = np.zeros(storage.no_of_shells, dtype=np.float64)
    storage.js = <float_type_t*> js.data
    storage.nubars = <float_type_t*> nubars.data
    storage.spectrum_start_nu = model.tardis_config.spectrum.frequency.value.min()
    storage.spectrum_end_nu = model.tardis_config.spectrum.frequency.value.max()
    storage.spectrum_delta_nu = model.tardis_config.spectrum.frequency.value[1] - model.tardis_config.spectrum.frequency.value[0]
    cdef np.ndarray[float_type_t, ndim=1] spectrum_virt_nu = model.montecarlo_virtual_luminosity
    storage.spectrum_virt_nu = <float_type_t*> spectrum_virt_nu.data
    storage.sigma_thomson = model.tardis_config.montecarlo.sigma_thomson.to('1/cm^2').value
    storage.inverse_sigma_thomson = 1.0 / storage.sigma_thomson
    storage.reflective_inner_boundary = model.tardis_config.montecarlo.enable_reflective_inner_boundary
    storage.inner_boundary_albedo = model.tardis_config.montecarlo.inner_boundary_albedo
    storage.current_packet_id = -1
    ######## Setting up the output ########
    #cdef np.ndarray[float_type_t, ndim=1] output_nus = np.zeros(storage.no_of_packets, dtype=np.float64)
    #cdef np.ndarray[float_type_t, ndim=1] output_energies = np.zeros(storage.no_of_packets, dtype=np.float64)
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
        storage.current_packet_id = i
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
        current_line_id = line_search(storage.line_list_nu, comov_current_nu, storage.no_of_lines)
        if current_line_id == storage.no_of_lines:
            #setting flag that the packet is off the red end of the line list
            last_line = 1
        else:
            last_line = 0
        #### FLAGS ####
        #Packet recently crossed the inner boundary
        recently_crossed_boundary = 1
        rpacket_init(&packet, &storage, current_nu, current_mu, current_energy, virtual_packet_flag)
        if (virtual_packet_flag > 0):
            #this is a run for which we want the virtual packet spectrum. So first thing we need to do is spawn virtual packets to track the input packet
            reabsorbed = montecarlo_one_packet(&storage, &current_nu, &current_energy, &current_mu, &current_shell_id,
                                               &current_r, &current_line_id, &last_line, &close_line,
                                               &recently_crossed_boundary, virtual_packet_flag,
                                               -1)
        #Now can do the propagation of the real packet
        reabsorbed = montecarlo_one_packet(&storage, &current_nu, &current_energy, &current_mu, &current_shell_id,
                                           &current_r, &current_line_id, &last_line, &close_line,
                                           &recently_crossed_boundary, virtual_packet_flag, 0)
        storage.output_nus[i] = current_nu
        storage.output_energies[i] = -current_energy if reabsorbed == 1 else current_energy
    return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id
