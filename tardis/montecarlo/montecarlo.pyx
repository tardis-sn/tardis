# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True



import numpy as np
cimport numpy as np
from numpy cimport PyArray_DATA
from astropy import constants
from astropy import units

from libc.stdlib cimport malloc, free

np.import_array()



ctypedef np.int64_t int_type_t

cdef extern from "src/cmontecarlo.h":
    ctypedef enum ContinuumProcessesStatus:
        CONTINUUM_OFF = 0
        CONTINUUM_ON = 1

    ctypedef enum FreeFreeStatus:
        FREE_FREE_OFF = 0
        FREE_FREE_ON = 1

    ctypedef struct photo_xsect_1level:
        double *nu
        double *x_sect
        int_type_t no_of_points

    ctypedef struct storage_model_t:
        double *packet_nus
        double *packet_mus
        double *packet_energies
        double *output_nus
        double *output_energies
        double *last_interaction_in_nu
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
        int_type_t *ion_charge
        double *ion_population
        int_type_t no_of_ions
        ContinuumProcessesStatus cont_status
        double *virt_packet_nus
        double *virt_packet_energies
        double *virt_packet_last_interaction_in_nu
        int_type_t *virt_packet_last_interaction_type
        int_type_t *virt_packet_last_line_interaction_in_id
        int_type_t *virt_packet_last_line_interaction_out_id
        int_type_t virt_packet_count
        int_type_t virt_array_size
        FreeFreeStatus ff_status
        photo_xsect_1level ** photo_xsect

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

    # Setting up the storage struct
    cdef storage_model_t storage

    storage.no_of_packets = model.packet_src.packet_nus.size
    storage.packet_nus = <double*> PyArray_DATA(model.packet_src.packet_nus)
    storage.packet_mus = <double*> PyArray_DATA(model.packet_src.packet_mus)
    storage.packet_energies = <double*> PyArray_DATA(
        model.packet_src.packet_energies)

    # Setup of structure
    structure = model.tardis_config.structure
    storage.no_of_shells = structure.no_of_shells


    storage.r_inner = <double*> PyArray_DATA(runner.r_inner_cgs)
    storage.r_outer = <double*> PyArray_DATA(runner.r_outer_cgs)
    storage.v_inner = <double*> PyArray_DATA(runner.v_inner_cgs)

    # Setup the rest
    # times
    storage.time_explosion = model.tardis_config.supernova.time_explosion.to(
        's').value
    storage.inverse_time_explosion = 1.0 / storage.time_explosion
    #electron density
    storage.electron_densities = <double*> PyArray_DATA(
        model.plasma_array.electron_densities.values)
    storage.inverse_electron_densities = <double*> PyArray_DATA(
        1.0 / model.plasma_array.electron_densities.values)

    # Switch for continuum processes
    if model.tardis_config.plasma['continuum_treatment'] == True:
        storage.cont_status = CONTINUUM_ON
    else:
        storage.cont_status = CONTINUUM_OFF
    # Continuum data

    # Switch for ff processes
    storage.ff_status = FREE_FREE_OFF
    # ff-data
    cdef np.ndarray[int_type_t, ndim=1] ion_charge
    cdef np.ndarray[double, ndim=2] ion_population
    if storage.ff_status == FREE_FREE_ON:
        ion_charge = model.plasma_array.ion_populations[0].index.get_level_values(1).values
        storage.ion_charge = <int_type_t*> ion_charge.data
        ion_population = model.plasma_array.ion_populations.values.transpose()
        storage.ion_population = <double *> ion_population.data
        storage.no_of_ions = ion_population.shape[1]

    # Line lists
    storage.no_of_lines = model.atom_data.lines.nu.values.size
    storage.line_list_nu = <double*> PyArray_DATA(model.atom_data.lines.nu.values)
    storage.line_lists_tau_sobolevs = <double*> PyArray_DATA(
        model.plasma_array.tau_sobolevs.values)
    storage.line_lists_j_blues = <double*> PyArray_DATA(runner.j_blue_estimator)

    storage.line_interaction_id = runner.get_line_interaction_id(
        model.tardis_config.plasma.line_interaction_type)

    # macro atom & downbranch
    if storage.line_interaction_id >= 1:
        storage.transition_probabilities = <double*> PyArray_DATA(model.transition_probabilities.values)
        storage.line2macro_level_upper = <int_type_t*> PyArray_DATA(
            model.atom_data.lines_upper2macro_reference_idx)
        storage.macro_block_references = <int_type_t*> PyArray_DATA(
            model.atom_data.macro_atom_references['block_references'].values)
        storage.transition_type = <int_type_t*> PyArray_DATA(
            model.atom_data.macro_atom_data['transition_type'].values)

        # Destination level is not needed and/or generated for downbranch
        storage.destination_level_id = <int_type_t*> PyArray_DATA(
            model.atom_data.macro_atom_data['destination_level_idx'].values)
        storage.transition_line_id = <int_type_t*> PyArray_DATA(
            model.atom_data.macro_atom_data['lines_idx'].values)

    storage.output_nus = <double*> PyArray_DATA(runner._packet_nu)
    storage.output_energies = <double*> PyArray_DATA(runner._packet_energy)

    storage.last_line_interaction_in_id = <int_type_t*> PyArray_DATA(
        runner.last_line_interaction_in_id)
    storage.last_line_interaction_out_id = <int_type_t*> PyArray_DATA(
        runner.last_line_interaction_out_id)
    storage.last_line_interaction_shell_id = <int_type_t*> PyArray_DATA(
        runner.last_line_interaction_shell_id)
    storage.last_interaction_type = <int_type_t*> PyArray_DATA(
        runner.last_interaction_type)
    storage.last_interaction_in_nu = <double*> PyArray_DATA(
        runner.last_interaction_in_nu)

    storage.js = <double*> PyArray_DATA(runner.j_estimator)
    storage.nubars = <double*> PyArray_DATA(runner.nu_bar_estimator)

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

    cdef np.ndarray[double, ndim=1] continuum_list_nu
    cdef np.ndarray[double, ndim =1] chi_bf_tmp_partial
    cdef np.ndarray[double, ndim=2] l_pop
    cdef np.ndarray[double, ndim=1] l_pop_r
    cdef np.ndarray[int_type_t, ndim=1] cont_edge2macro_continuum
    cdef np.ndarray[int_type_t, ndim=1] macro_block_references_continuum
    cdef np.ndarray[int_type_t, ndim=1] transition_type_continuum
    cdef np.ndarray[int_type_t, ndim=1] destination_level_id_continuum
    cdef np.ndarray[int_type_t, ndim=1] transition_continuum_id
    cdef np.ndarray[double, ndim=2] transition_probabilities_continuum

    # kind of redundant atm
    cdef int no_levels_with_photdata = model.atom_data.continuum_data.no_levels_with_photdata
    # Temporary
    cdef photo_xsect_1level ** photo_xsect = <photo_xsect_1level **> malloc(
        no_levels_with_photdata * sizeof(photo_xsect_1level *))

    if storage.cont_status == CONTINUUM_ON:
        for i in range(no_levels_with_photdata):
            photo_xsect[i] = <photo_xsect_1level *> malloc(sizeof(photo_xsect_1level))
            phot_table_xsect = model.atom_data.continuum_data.get_phot_table_xsect(i)
            phot_table_nu = model.atom_data.continuum_data.get_phot_table_nu(i)
            photo_xsect[i].no_of_points = len(phot_table_xsect)
            photo_xsect[i].nu = <double *> (<np.ndarray[double, ndim =1]> (phot_table_nu)).data
            photo_xsect[i].x_sect = <double *> (<np.ndarray[double, ndim =1]> (phot_table_xsect)).data

    if storage.cont_status == CONTINUUM_ON:
        continuum_list_nu = model.atom_data.continuum_data.continuum_edges_list
        transition_probabilities_continuum = model.transition_probabilities_continuum.data.ix[:, 3:].values.transpose()
        transition_probabilities_nd_continuum = len(model.transition_probabilities_continuum.data)
        macro_block_references_continuum = model.atom_data.continuum_data.continuum_references[
            'block_references'].values
        cont_edge2macro_continuum = model.atom_data.continuum_data.cont_edge2macro_continuum
        destination_level_id_continuum = model.transition_probabilities_continuum.data['destination_level_idx'].values
        transition_type_continuum = model.transition_probabilities_continuum.data['transition_type'].values
        transition_continuum_id = model.transition_probabilities_continuum.data['continuum_edge_idx'].values
        l_pop = model.atom_data.continuum_data.level_number_density
        chi_bf_tmp_partial = np.zeros(continuum_list_nu.size)
        l_pop_r = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)

        storage.l_pop = <double*> l_pop.data
        storage.transition_probabilities_continuum = <double*> transition_probabilities_continuum.data
        storage.continuum_list_nu = <double*> continuum_list_nu.data
        storage.destination_level_id_continuum = <int_type_t*> destination_level_id_continuum.data
        storage.transition_type_continuum = <int_type_t*> transition_type_continuum.data
        storage.transition_continuum_id = <int_type_t*> transition_continuum_id.data
        storage.macro_block_references_continuum = <int_type_t*> macro_block_references_continuum.data
        storage.chi_bf_tmp_partial = <double*> chi_bf_tmp_partial.data
        storage.l_pop_r = <double*> l_pop_r.data
        storage.cont_edge2macro_continuum = <int_type_t*> cont_edge2macro_continuum.data
        storage.transition_probabilities_nd_continuum = transition_probabilities_nd_continuum
        storage.no_of_edges = continuum_list_nu.size
        storage.photo_xsect = photo_xsect

    montecarlo_main_loop(&storage, virtual_packet_flag, nthreads,
                         model.tardis_config.montecarlo.seed)

    cdef np.ndarray[double, ndim=1] virt_packet_nus = np.zeros(storage.virt_packet_count, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] virt_packet_energies = np.zeros(storage.virt_packet_count, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] virt_packet_last_interaction_in_nu = np.zeros(storage.virt_packet_count, dtype=np.float64)
    cdef np.ndarray[int_type_t, ndim=1] virt_packet_last_interaction_type = np.zeros(storage.virt_packet_count, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] virt_packet_last_line_interaction_in_id = np.zeros(storage.virt_packet_count, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] virt_packet_last_line_interaction_out_id = np.zeros(storage.virt_packet_count, dtype=np.int64)

    for i in range(storage.virt_packet_count):
        virt_packet_nus[i] = storage.virt_packet_nus[i]
        virt_packet_energies[i] = storage.virt_packet_energies[i]
        virt_packet_last_interaction_in_nu[i] = storage.virt_packet_last_interaction_in_nu[i]
        virt_packet_last_interaction_type[i] = storage.virt_packet_last_interaction_type[i]
        virt_packet_last_line_interaction_in_id[i] = storage.virt_packet_last_line_interaction_in_id[i]
        virt_packet_last_line_interaction_out_id[i] = storage.virt_packet_last_line_interaction_out_id[i]
    free(<void *>storage.virt_packet_nus)
    free(<void *>storage.virt_packet_energies)
    free(<void *>storage.virt_packet_last_interaction_in_nu)
    free(<void *>storage.virt_packet_last_interaction_type)
    free(<void *>storage.virt_packet_last_line_interaction_in_id)
    free(<void *>storage.virt_packet_last_line_interaction_out_id)

    # Necessary?
    if storage.cont_status == CONTINUUM_ON:
        for i in range(storage.no_of_edges):
            free(<photo_xsect_1level *> storage.photo_xsect[i])

    runner.virt_packet_nus = virt_packet_nus
    runner.virt_packet_energies = virt_packet_energies
    runner.virt_packet_last_interaction_in_nu = virt_packet_last_interaction_in_nu
    runner.virt_packet_last_interaction_type = virt_packet_last_interaction_type
    runner.virt_packet_last_line_interaction_in_id = virt_packet_last_line_interaction_in_id
    runner.virt_packet_last_line_interaction_out_id = virt_packet_last_line_interaction_out_id
    
    #return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id, virt_packet_nus, virt_packet_energies


