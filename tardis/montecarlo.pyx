# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

IF OPENMP:
    from cython.parallel cimport parallel, prange, threadid
    cimport openmp
    cdef int num_threads
    openmp.omp_set_dynamic(1)


np.import_array()

ctypedef np.int64_t int_type_t

cdef extern from "cmontecarlo.h":
    ctypedef enum rpacket_status_t:
        TARDIS_PACKET_STATUS_IN_PROCESS = 0
        TARDIS_PACKET_STATUS_EMITTED = 1
        TARDIS_PACKET_STATUS_REABSORBED = 2

    ctypedef struct rpacket_t:
        double nu
        double mu
        double energy
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

    int_type_t montecarlo_one_packet(storage_model_t *storage, rpacket_t *packet, int_type_t virtual_mode) nogil
    int rpacket_init(rpacket_t *packet, storage_model_t *storage, int packet_index, int virtual_packet_flag) nogil
    double rpacket_get_nu(rpacket_t *packet) nogil
    double rpacket_get_energy(rpacket_t *packet)  nogil
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
    cdef rpacket_t* packet
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
    ######## Setting up the output ########
    #cdef np.ndarray[double, ndim=1] output_nus = np.zeros(storage.no_of_packets, dtype=np.float64)
    #cdef np.ndarray[double, ndim=1] output_energies = np.zeros(storage.no_of_packets, dtype=np.float64)
    cdef int_type_t reabsorbed = 0
    cdef int_type_t no_of_packets = storage.no_of_packets
    cdef int_type_t packet_index = 0
    cdef int_type_t num_threads


    IF OPENMP:
        with nogil, parallel():
            num_threads = openmp.omp_get_num_threads()
            for packet_index in prange(no_of_packets):
                packet = <rpacket_t *> malloc(sizeof(rpacket_t))
                if not packet_index % (no_of_packets/20):
                    printf("%d\n",packet_index)
                storage.current_packet_id = packet_index
                rpacket_init(packet, &storage, packet_index, virtual_packet_flag)
                if (virtual_packet_flag > 0):
                #this is a run for which we want the virtual packet spectrum. So first thing we need to do is spawn virtual packets to track the input packet
                    reabsorbed = montecarlo_one_packet(&storage, packet, -1)
            #Now can do the propagation of the real packet
                reabsorbed = montecarlo_one_packet(&storage, packet, 0)
                storage.output_nus[packet_index] = rpacket_get_nu(packet)
                if reabsorbed ==1 :
                    storage.output_energies[packet_index] = -rpacket_get_energy(packet)
                else:
                    storage.output_energies[packet_index] = rpacket_get_energy(packet)
                free(packet)
    ELSE:
        for packet_index in range(no_of_packets):
            packet = <rpacket_t *> malloc(sizeof(rpacket_t))
            if not packet_index % (no_of_packets/20):
                printf("%d\n",packet_index)
            storage.current_packet_id = packet_index
            rpacket_init(packet, &storage, packet_index, virtual_packet_flag)
            if (virtual_packet_flag > 0):
                #this is a run for which we want the virtual packet spectrum. So first thing we need to do is spawn virtual packets to track the input packet
                reabsorbed = montecarlo_one_packet(&storage, packet, -1)
                #Now can do the propagation of the real packet
            reabsorbed = montecarlo_one_packet(&storage, packet, 0)
            storage.output_nus[packet_index] = rpacket_get_nu(packet)
            if reabsorbed ==1 :
                storage.output_energies[packet_index] = -rpacket_get_energy(packet)
            else:
                storage.output_energies[packet_index] = rpacket_get_energy(packet)
            free(packet)


    return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id

