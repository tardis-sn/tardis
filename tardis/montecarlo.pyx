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

ctypedef np.int64_t int_type_t

cdef extern from "cmontecarlo.h":
    ctypedef enum rpacket_status_t:
        TARDIS_PACKET_STATUS_IN_PROCESS = 1 << 3
        TARDIS_PACKET_STATUS_EMITTED = 1 << 2
        TARDIS_PACKET_STATUS_REABSORBED = 1 << 1
        TARDIS_PACKET_STATUS_DISABLED = 1 << 0

    ctypedef enum packet_status_t:
        TARDIS_R_PACKET_IN_PROCESS = 1 << 6 | 1 << 3
        TARDIS_R_PACKET_STATUS_EMITTED = 1 << 6 | 1 << 2
        TARDIS_R_PACKET_STATUS_REABSORBED = 1 << 6 | 1 << 1
        TARDIS_R_PACKET_STATUS_DISABLED = 1 << 6 | 1 << 0
        TARDIS_K_PACKET_IN_PROCESS = 1 << 5 | 1 << 3
        TARDIS_K_PACKET_STATUS_DISABLED = 1 << 5 | 1 << 0
        TARDIS_I_PACKET_IN_PROCESS = 1 << 4 | 1 << 3
        TARDIS_I_PACKET_STATUS_DISABLED = 1 << 4 | 1 << 0

    ctypedef struct rpacket_t:
        double nu
        double mu
        double energy
        double comov_nu
        double comov_energy
        double r
        double tau_event
        double nu_line
        int_type_t current_shell_id
        int_type_t next_line_id
        int_type_t last_line
        int_type_t close_line
        int_type_t recently_crossed_boundary
        int_type_t virtual_packet_flag
        int_type_t virtual_packet
        double d_line
        double d_electron
        double d_boundary
        rpacket_status_t next_shell_id
        packet_status_t status
        double chi_bf
        double chi_th
        double chi_ff
        double chi_cont
        double d_bf
        double d_th
        double d_ff
        double d_cont
        double last_bf_edge
        double *chi_bf_tmp_partial
        int_type_t chi_bf_tmp_partial_last_shell_id
        double chi_bf_tmp_partial_last_nu
        double Cr_fb_max
        double Cr_ff_max
        double Cr_bb_max
        double Cr_ion_max


    ctypedef struct storage_model_t:
        double *packet_nus
        double *packet_mus
        double *packet_energies
        double *output_nus
        double *output_energies
        int_type_t *last_line_interaction_in_id
        int_type_t *last_line_interaction_out_id
        int_type_t *last_line_interaction_shell_id
        int_type_t *last_interaction_type
        int_type_t no_of_packets
        int_type_t no_of_shells
        double *r_inner
        double *r_outer
        double *v_inner
        double time_explosion
        double inverse_time_explosion
        double *electron_densities
        double *inverse_electron_densities
        double *line_list_nu
        double *line_lists_tau_sobolevs
        int_type_t line_lists_tau_sobolevs_nd
        double *line_lists_j_blues
        int_type_t line_lists_j_blues_nd
        int_type_t no_of_lines
        int_type_t line_interaction_id
        double *transition_probabilities
        int_type_t transition_probabilities_nd
        int_type_t *line2macro_level_upper
        int_type_t *macro_block_references
        int_type_t *transition_type
        int_type_t *destination_level_id
        int_type_t *transition_line_id
        double *js
        double *nubars
        double spectrum_start_nu
        double spectrum_delta_nu
        double spectrum_end_nu
        double *spectrum_virt_nu
        double sigma_thomson
        double inverse_sigma_thomson
        double inner_boundary_albedo
        int_type_t reflective_inner_boundary
        int_type_t current_packet_id

        int_type_t *chi_bf_index_to_level
        int_type_t chi_bf_index_to_level_nrow
        int_type_t chi_bf_index_to_level_ncolum

        double *bf_level_population
        int_type_t bf_level_population_nrow
        int_type_t bf_level_population_ncolum

        double *bf_lpopulation_ratio
        int_type_t bf_lpopulation_ratio_nrow
        int_type_t bf_lpopulation_ratio_ncolum

        double *bf_lpopulation_ratio_nlte_lte
        int_type_t bf_lpopulation_ratio_nlte_lte_nrow
        int_type_t bf_lpopulation_ratio_nlte_lte_ncolum

        double *bf_cross_sections

        double *bound_free_th_frequency

        double *t_electrons
        double *chi_bf_tmp_partial

        double *Cr_fb_ijk_all
        int_type_t Cr_fb_ijk_all_nrow
        int_type_t Cr_fb_ijk_all_ncolum

        double *Cr_fb_ijk_cumsum_all
        int_type_t Cr_fb_ijk_cumsum_all_nrow
        int_type_t Cr_fb_ijk_cumsum_all_ncolum
        int_type_t *Cr_fb_ijk_index

        double *Cr_fb_ijk_th_frequency
        int_type_t Cr_fb_ijk_th_frequency_nrow
        int_type_t Cr_fb_ijk_th_frequency_ncolum

        double *Cr_ff_jk_all
        int_type_t Cr_ff_jk_all_nrow
        int_type_t Cr_ff_jk_all_ncolum

        double *Cr_ff_jk_cumsum_all
        int_type_t Cr_ff_jk_cumsum_all_nrow
        int_type_t Cr_ff_jk_cumsum_all_ncolum
        int_type_t *Cr_ff_jk_index

        double *Cr_bb_ijk_all
        int_type_t Cr_bb_ijk_all_nrow
        int_type_t Cr_bb_ijk_all_ncolum

        double *Cr_bb_ijk_cumsum_all
        int_type_t Cr_bb_ijk_cumsum_all_nrow
        int_type_t Cr_bb_ijk_cumsum_all_ncolum

        int_type_t *Cr_bb_ijk_index

        double *Cr_ion_ijk_all
        int_type_t Cr_ion_ijk_all_nrow
        int_type_t Cr_ion_ijk_all_ncolum

        double *Cr_ion_ijk_cumsum_all
        int_type_t Cr_ion_ijk_cumsum_all_nrow
        int_type_t Cr_ion_ijk_cumsum_all_ncolum

        int_type_t *Cr_ion_ijk_index

    #double kB

    int_type_t montecarlo_one_packet(storage_model_t *storage, rpacket_t *packet, int_type_t virtual_mode)
    int rpacket_init(rpacket_t *packet, storage_model_t *storage, int packet_index, int virtual_packet_flag)
    double rpacket_get_nu(rpacket_t *packet)
    double rpacket_get_energy(rpacket_t *packet)
    void initialize_random_kit(unsigned long seed)




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
                    np.ndarray[double, ndim=1] line_list_nu,
                    np.ndarray[double, ndim=2] tau_lines,
                    np.ndarray[double, ndim=1] ne,
                    double packet_energy,
                    np.ndarray[double, ndim=2] p_transition,
                    np.ndarray[int_type_t, ndim=1] type_transition,
                    np.ndarray[int_type_t, ndim=1] target_level_id,
                    np.ndarray[int_type_t, ndim=1] target_line_id,
                    np.ndarray[int_type_t, ndim=1] unroll_reference,
                    np.ndarray[int_type_t, ndim=1] line2level,
                    int_type_t log_packets,
                    int_type_t do_scatter
    """
    print("Start montecarlo_radial1d")
    cdef storage_model_t storage
    cdef rpacket_t packet
    initialize_random_kit(model.tardis_config.montecarlo.seed)
    cdef np.ndarray[double, ndim=1] packet_nus = model.packet_src.packet_nus
    storage.packet_nus = <double*> packet_nus.data
    cdef np.ndarray[double, ndim=1] packet_mus = model.packet_src.packet_mus
    storage.packet_mus = <double*> packet_mus.data
    cdef np.ndarray[double, ndim=1] packet_energies = model.packet_src.packet_energies
    storage.packet_energies = <double*> packet_energies.data
    storage.no_of_packets = packet_nus.size
    # Setup of structure
    structure = model.tardis_config.structure
    storage.no_of_shells = structure.no_of_shells
    cdef np.ndarray[double, ndim=1] r_inner = structure.r_inner.to('cm').value
    storage.r_inner = <double*> r_inner.data
    cdef np.ndarray[double, ndim=1] r_outer = structure.r_outer.to('cm').value
    storage.r_outer = <double*> r_outer.data
    cdef np.ndarray[double, ndim=1] v_inner = structure.v_inner.to('cm/s').value
    storage.v_inner = <double*> v_inner.data
    # Setup the rest
    # times
    storage.time_explosion = model.tardis_config.supernova.time_explosion.to('s').value
    storage.inverse_time_explosion = 1.0 / storage.time_explosion
    #electron density
    cdef np.ndarray[double, ndim=1] electron_densities = model.plasma_array.electron_densities.values
    storage.electron_densities = <double*> electron_densities.data
    cdef np.ndarray[double, ndim=1] inverse_electron_densities = 1.0 / electron_densities
    storage.inverse_electron_densities = <double*> inverse_electron_densities.data
    # Line lists
    cdef np.ndarray[double, ndim=1] line_list_nu = model.atom_data.lines.nu.values
    storage.line_list_nu = <double*> line_list_nu.data
    storage.no_of_lines = line_list_nu.size
    cdef np.ndarray[double, ndim=2] line_lists_tau_sobolevs = model.plasma_array.tau_sobolevs.values.transpose()
    storage.line_lists_tau_sobolevs = <double*> line_lists_tau_sobolevs.data
    storage.line_lists_tau_sobolevs_nd = line_lists_tau_sobolevs.shape[1]
    cdef np.ndarray[double, ndim=2] line_lists_j_blues = model.j_blue_estimators
    storage.line_lists_j_blues = <double*> line_lists_j_blues.data
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
    cdef np.ndarray[double, ndim=2] transition_probabilities
    cdef np.ndarray[int_type_t, ndim=1] line2macro_level_upper
    cdef np.ndarray[int_type_t, ndim=1] macro_block_references
    cdef np.ndarray[int_type_t, ndim=1] transition_type
    cdef np.ndarray[int_type_t, ndim=1] destination_level_id
    cdef np.ndarray[int_type_t, ndim=1] transition_line_id
    if storage.line_interaction_id >= 1:
        transition_probabilities = model.transition_probabilities.values.transpose()
        storage.transition_probabilities = <double*> transition_probabilities.data
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
    cdef np.ndarray[double, ndim=1] output_nus = np.zeros(storage.no_of_packets, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] output_energies = np.zeros(storage.no_of_packets, dtype=np.float64)
    storage.output_nus = <double*> output_nus.data
    storage.output_energies = <double*> output_energies.data
    cdef np.ndarray[int_type_t, ndim=1] last_line_interaction_in_id = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] last_line_interaction_out_id = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] last_line_interaction_shell_id = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] last_interaction_type = -1 * np.ones(storage.no_of_packets, dtype=np.int64)
    storage.last_line_interaction_in_id = <int_type_t*> last_line_interaction_in_id.data
    storage.last_line_interaction_out_id = <int_type_t*> last_line_interaction_out_id.data
    storage.last_line_interaction_shell_id = <int_type_t*> last_line_interaction_shell_id.data
    storage.last_interaction_type = <int_type_t*> last_interaction_type.data
    cdef np.ndarray[double, ndim=1] js = np.zeros(storage.no_of_shells, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] nubars = np.zeros(storage.no_of_shells, dtype=np.float64)
    storage.js = <double*> js.data
    storage.nubars = <double*> nubars.data
    storage.spectrum_start_nu = model.tardis_config.spectrum.frequency.value.min()
    storage.spectrum_end_nu = model.tardis_config.spectrum.frequency.value.max()
    storage.spectrum_delta_nu = model.tardis_config.spectrum.frequency.value[1] - model.tardis_config.spectrum.frequency.value[0]
    cdef np.ndarray[double, ndim=1] spectrum_virt_nu = model.montecarlo_virtual_luminosity
    storage.spectrum_virt_nu = <double*> spectrum_virt_nu.data
    storage.sigma_thomson = model.tardis_config.montecarlo.sigma_thomson.to('1/cm^2').value
    storage.inverse_sigma_thomson = 1.0 / storage.sigma_thomson
    storage.reflective_inner_boundary = model.tardis_config.montecarlo.enable_reflective_inner_boundary
    storage.inner_boundary_albedo = model.tardis_config.montecarlo.inner_boundary_albedo
    storage.current_packet_id = -1


    cdef np.ndarray[double, ndim=1, mode='c'] chi_bf_tmp_partial = np.zeros(2 * storage.bf_level_population_nrow, dtype=np.double)
    storage.chi_bf_tmp_partial = <double *> chi_bf_tmp_partial.data

    ######## Setting up the output ########
    #cdef np.ndarray[double, ndim=1] output_nus = np.zeros(storage.no_of_packets, dtype=np.float64)
    #cdef np.ndarray[double, ndim=1] output_energies = np.zeros(storage.no_of_packets, dtype=np.float64)
    cdef int_type_t reabsorbed = 0

    for packet_index in range(storage.no_of_packets):
        if not packet_index % (storage.no_of_packets/20):
            print(packet_index)
        storage.current_packet_id = packet_index
        rpacket_init(&packet, &storage, packet_index, virtual_packet_flag)
        if (virtual_packet_flag > 0):
            #this is a run for which we want the virtual packet spectrum. So first thing we need to do is spawn virtual packets to track the input packet
            reabsorbed = montecarlo_one_packet(&storage, &packet, -1)
        #Now can do the propagation of the real packet
        reabsorbed = montecarlo_one_packet(&storage, &packet, 0)
        storage.output_nus[packet_index] = rpacket_get_nu(&packet)
        storage.output_energies[packet_index] = -rpacket_get_energy(&packet) if reabsorbed == 1 else rpacket_get_energy(
            &packet)
    return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id

