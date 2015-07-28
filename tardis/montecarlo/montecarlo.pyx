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
        double *transition_probabilities_continuum
        int_type_t transition_probabilities_nd_continuum
        int_type_t *cont_edge2macro_continuum  # continuum equivalent to line2macro_level_upper
        int_type_t *macro_block_references_continuum
        int_type_t *transition_type_continuum
        int_type_t *destination_level_id_continuum
        int_type_t *transition_continuum_id  # connects index i in transition_probabilities_nd_continuum to continuum_id
        # of emission
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

    void montecarlo_main_loop(storage_model_t * storage, int_type_t virtual_packet_flag, int nthreads, unsigned long seed)

def montecarlo_radial1d(model, int_type_t virtual_packet_flag=0, int nthreads=4):
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
    cdef np.ndarray[double, ndim=1] transition_probabilities_continuum
    cdef np.ndarray[int_type_t, ndim=1] cont_edge2macro_continuum
    cdef np.ndarray[int_type_t, ndim=1] macro_block_references_continuum
    cdef np.ndarray[int_type_t, ndim=1] transition_type_continuum
    cdef np.ndarray[int_type_t, ndim=1] destination_level_id_continuum
    cdef np.ndarray[int_type_t, ndim=1] transition_continuum_id
    if storage.cont_status == CONTINUUM_ON:
        I_H = 13.5984 * units.eV.to(units.erg)
        #I_HeI = 24.5874 * units.eV.to(units.erg)
        #I_HeII = 54.417760 * units.eV.to(units.erg)
        energy_list_nu_H = I_H - model.atom_data.levels.ix[(1, 0)].energy.values
        cont_list_nu_H = (units.erg).to(units.Hz, equivalencies=units.spectral()) * energy_list_nu_H
        continuum_list_nu = cont_list_nu_H
        #print continuum_list_nu , type(continuum_list_nu)
        #time.sleep(20)
        #cont_list_nu_HeI = I_HeI - model.atom_data.levels.ix[(2, 0)].energy.values
        #cont_list_nu_HeII = I_HeII - model.atom_data.levels.ix[(2, 1)].energy.values
        #continuum_list_nu = np.concatenate((cont_list_nu_H, cont_list_nu_HeI, cont_list_nu_HeII))
        #sorting_order = np.argsort(continuum_list_nu)[::-1]
        #continuum_list_nu = continuum_list_nu[sorting_order]
        #l_pop = np.zeros(0)
        #for i in range(0, storage.no_of_shells):
        #    l_pop_H = model.plasma_array.level_populations[0].ix[(1, 0)].values
        #    l_pop_HeI = model.plasma_array.level_populations[0].ix[(2, 0)].values
        #    l_pop_HeII = model.plasma_array.level_populations[0].ix[(2, 1)].values
        #    l_pop_shell = np.concatenate((l_pop_H, l_pop_HeI, l_pop_HeII))
        #    l_pop_shell = l_pop_shell[sorting_order]
        #    l_pop = np.append(l_pop, l_pop_shell)
        # Prepare continuum data for macro_atom_new
        edge2cont_HI = np.zeros(len(cont_list_nu_H), dtype=np.int64)
        cont_edge2macro_continuum = edge2cont_HI
        macro_block_references_continuum = np.zeros(1, dtype=np.int64)
        transition_probabilities_HI_cont = np.ones(2 * len(cont_list_nu_H) * storage.no_of_shells, dtype=np.float64)
        transition_probabilities_continuum = transition_probabilities_HI_cont
        #
        transition_type_continuum = np.ones(2 * len(cont_list_nu_H), dtype=np.int64) * (-3)
        destination_level_id_continuum = np.ones(2 * len(cont_list_nu_H), dtype=np.int64) * 15
        transition_continuum_id = np.ones(2 * len(cont_list_nu_H), dtype=np.int64) * 10  # not needed atm

        #edge2cont_HeI = np.ones(len(cont_list_nu_HeI), dtype = np.int64)
        #edge2cont_HeII = np.ones(len(cont_list_nu_HeII), dtype = np.int64) * 2
        #cont_edge2macro_continuum = np.concatenate((edge2cont_HI, edge2cont_HeI, edge2cont_HeII))
        #cont_edge2macro_continuum = cont_edge2macro_continuum[sorting_order]

        storage.transition_probabilities_continuum = <double*> transition_probabilities_continuum.data
        storage.transition_probabilities_nd_continuum = 2 * len(cont_list_nu_H)
        storage.cont_edge2macro_continuum = <int_type_t*> cont_edge2macro_continuum.data
        storage.macro_block_references_continuum = <int_type_t*> macro_block_references_continuum.data
        storage.transition_type_continuum = <int_type_t*> transition_type_continuum.data
        storage.destination_level_id_continuum = <int_type_t*> destination_level_id_continuum.data
        storage.transition_continuum_id = <int_type_t*> transition_continuum_id.data
        #continuum_list_nu = np.array(
        #    [9.0e14, 8.223e14, 6.0e14, 3.5e14, 3.0e14])  # sorted list of threshold frequencies
        #l_pop = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)

        l_pop = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)

        storage.continuum_list_nu = <double*> continuum_list_nu.data
        storage.no_of_edges = continuum_list_nu.size
        chi_bf_tmp_partial = np.zeros(continuum_list_nu.size)
        storage.chi_bf_tmp_partial = <double*> chi_bf_tmp_partial.data
        storage.l_pop = <double*> l_pop.data
        l_pop_r = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)
        storage.l_pop_r = <double*> l_pop_r.data
        #storage.cont_edge2macro_continuum = <int_type_t*> cont_edge2macro_continuum.data
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
    return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id

