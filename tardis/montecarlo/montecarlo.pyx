# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True


import logging
import time

import numpy as np
cimport numpy as np
from astropy import constants
from astropy import units
from libc.stdlib cimport free

np.import_array()

ctypedef np.int64_t int_type_t

cdef extern from "src/cmontecarlo.h":
    ctypedef enum ContinuumProcessesStatus:
        CONTINUUM_OFF = 0
        CONTINUUM_ON = 1

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
        double *continuum_list_nu
        int_type_t line_lists_tau_sobolevs_nd
        double *line_lists_j_blues
        int_type_t line_lists_j_blues_nd
        int_type_t no_of_lines
        int_type_t no_of_edges
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
        double spectrum_virt_start_nu
        double spectrum_virt_end_nu
        double spectrum_start_nu
        double spectrum_delta_nu
        double spectrum_end_nu
        double *spectrum_virt_nu
        double sigma_thomson
        double inverse_sigma_thomson
        double inner_boundary_albedo
        int_type_t reflective_inner_boundary
        double *chi_bf_tmp_partial
        double *t_electrons
        double *l_pop
        double *l_pop_r
        ContinuumProcessesStatus cont_status
        double *virt_packet_nus
        double *virt_packet_energies
        int_type_t virt_packet_count
        int_type_t virt_array_size

    void montecarlo_main_loop(storage_model_t * storage, int_type_t virtual_packet_flag, int nthreads, unsigned long seed)

def montecarlo_radial1d(model, runner, int_type_t virtual_packet_flag=0,
                        int nthreads=4):
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
    cdef storage_model_t storage
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
    # Switch for continuum processes
    storage.cont_status = CONTINUUM_OFF
    # Continuum data
    cdef np.ndarray[double, ndim=1] continuum_list_nu
    cdef np.ndarray[double, ndim =1] chi_bf_tmp_partial
    cdef np.ndarray[double, ndim=1] l_pop
    cdef np.ndarray[double, ndim=1] l_pop_r
    if storage.cont_status == CONTINUUM_ON:
        continuum_list_nu = np.array([9.0e14, 8.223e14, 6.0e14, 3.5e14, 3.0e14])  # sorted list of threshold frequencies
        storage.continuum_list_nu = <double*> continuum_list_nu.data
        storage.no_of_edges = continuum_list_nu.size
        chi_bf_tmp_partial = np.zeros(continuum_list_nu.size)
        storage.chi_bf_tmp_partial = <double*> chi_bf_tmp_partial.data
        l_pop = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)
        storage.l_pop = <double*> l_pop.data
        l_pop_r = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)
        storage.l_pop_r = <double*> l_pop_r.data
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
    cdef double [:, :] transition_probabilities
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
    storage.spectrum_virt_start_nu = model.tardis_config.montecarlo.virtual_spectrum_range.end.to('Hz', units.spectral()).value
    storage.spectrum_virt_end_nu = model.tardis_config.montecarlo.virtual_spectrum_range.start.to('Hz', units.spectral()).value
    storage.spectrum_delta_nu = model.tardis_config.spectrum.frequency.value[1] - model.tardis_config.spectrum.frequency.value[0]
    cdef np.ndarray[double, ndim=1] spectrum_virt_nu = model.montecarlo_virtual_luminosity
    storage.spectrum_virt_nu = <double*> spectrum_virt_nu.data
    storage.sigma_thomson = model.tardis_config.montecarlo.sigma_thomson.to('1/cm^2').value
    storage.inverse_sigma_thomson = 1.0 / storage.sigma_thomson
    storage.reflective_inner_boundary = model.tardis_config.montecarlo.enable_reflective_inner_boundary
    storage.inner_boundary_albedo = model.tardis_config.montecarlo.inner_boundary_albedo
    # Data for continuum implementation
    cdef np.ndarray[double, ndim=1] t_electrons = model.plasma_array.t_electrons
    storage.t_electrons = <double*> t_electrons.data
    ######## Setting up the output ########
    #cdef np.ndarray[double, ndim=1] output_nus = np.zeros(storage.no_of_packets, dtype=np.float64)
    #cdef np.ndarray[double, ndim=1] output_energies = np.zeros(storage.no_of_packets, dtype=np.float64)
    montecarlo_main_loop(&storage, virtual_packet_flag, nthreads, model.tardis_config.montecarlo.seed)

    cdef np.ndarray[double, ndim=1] virt_packet_nus = np.zeros(storage.virt_packet_count, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] virt_packet_energies = np.zeros(storage.virt_packet_count, dtype=np.float64)
    for i in range(storage.virt_packet_count):
        virt_packet_nus[i] = storage.virt_packet_nus[i]
        virt_packet_energies[i] = storage.virt_packet_energies[i]
    free(<void *>storage.virt_packet_nus)
    free(<void *>storage.virt_packet_energies)
    runner._packet_nu = output_nus
    runner._packet_energy = output_energies
    runner.j_estimator = js
    runner.nu_bar_estimator = nubars
    runner.last_line_interaction_in_id = last_line_interaction_in_id
    runner.last_line_interaction_out_id = last_line_interaction_out_id
    runner.last_interaction_type = last_interaction_type
    runner.last_line_interaction_shell_id = last_line_interaction_shell_id
    runner.virt_packet_nus = virt_packet_nus
    runner.virt_packet_energies = virt_packet_energies
    
    #return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id, virt_packet_nus, virt_packet_energies


