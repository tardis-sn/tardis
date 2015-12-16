# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True



import numpy as np
cimport numpy as np
from numpy cimport PyArray_DATA
from astropy import constants
from astropy import units
from libc.stdlib cimport free

from tardis.continuum.exceptions import ContinuumBuildError

from libc.stdlib cimport malloc, free

np.import_array()



ctypedef np.int64_t int_type_t

cdef extern from "src/cmontecarlo.h":
    ctypedef enum ContinuumProcessesStatus:
        CONTINUUM_OFF = 0
        CONTINUUM_ON = 1

    cdef int LOG_VPACKETS

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
        double *chi_ff_factor
        double *fb_cooling_prob
        double *ff_cooling_prob
        double *coll_exc_cooling_prob
        double *coll_ion_cooling_prob
        double *fb_cooling_prob_individual
        double *coll_exc_cooling_prob_individual
        double *coll_ion_cooling_prob_individual
        int_type_t *fb_cooling_references
        int_type_t *coll_ion_cooling_references
        int_type_t *coll_exc_cooling_references
        int_type_t fb_cooling_prob_nd
        int_type_t coll_ion_cooling_prob_nd
        int_type_t coll_exc_cooling_prob_nd

    void montecarlo_main_loop(storage_model_t * storage, int_type_t virtual_packet_flag, int nthreads, unsigned long seed)




cdef initialize_storage_model(model, runner, storage_model_t *storage):
    """
    Initializing the storage struct.

    """

    storage.no_of_packets = runner.input_nu.size
    storage.packet_nus = <double*> PyArray_DATA(runner.input_nu)
    storage.packet_mus = <double*> PyArray_DATA(runner.input_mu)
    storage.packet_energies = <double*> PyArray_DATA(runner.input_energy)

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

    runner.inverse_electron_densities = (
        1.0 / model.plasma_array.electron_densities.values)
    storage.inverse_electron_densities = <double*> PyArray_DATA(
        runner.inverse_electron_densities)

    # Switch for continuum processes
    if model.tardis_config.plasma['continuum_treatment'] == True:
        if WITH_CONTINUUM:
            storage.cont_status = CONTINUUM_ON
            storage.ff_status = FREE_FREE_ON
        else:
            raise ContinuumBuildError
    else:
        storage.cont_status = CONTINUUM_OFF
        storage.ff_status = FREE_FREE_OFF

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
        storage.line2macro_level_upper = <int_type_t*> PyArray_DATA(
            model.atom_data.lines_upper2macro_reference_idx)

        if model.tardis_config.plasma['continuum_treatment'] == False:
            storage.transition_probabilities = <double*> PyArray_DATA(
                model.plasma_array.transition_probabilities.values)
            storage.macro_block_references = <int_type_t*> PyArray_DATA(
                model.atom_data.macro_atom_references['block_references'].values)
            storage.transition_type = <int_type_t*> PyArray_DATA(
                model.atom_data.macro_atom_data['transition_type'].values)

            # Destination level is not needed and/or generated for downbranch
            storage.destination_level_id = <int_type_t*> PyArray_DATA(
                model.atom_data.macro_atom_data['destination_level_idx'].values)
            storage.transition_line_id = <int_type_t*> PyArray_DATA(
                model.atom_data.macro_atom_data['lines_idx'].values)
        else:
            storage.transition_probabilities = <double*> PyArray_DATA(
                model.base_continuum.transition_probabilities.data_array)
            storage.macro_block_references = <int_type_t*> PyArray_DATA(
                model.base_continuum.transition_probabilities.block_references)
            storage.destination_level_id = <int_type_t*> PyArray_DATA(
                model.base_continuum.transition_probabilities.destination_level_id)
            storage.transition_type = <int_type_t*> PyArray_DATA(
                model.base_continuum.transition_probabilities.transition_type)
            storage.transition_line_id = <int_type_t*> PyArray_DATA(
                model.base_continuum.transition_probabilities.transition_line_id)

    storage.output_nus = <double*> PyArray_DATA(runner._output_nu)
    storage.output_energies = <double*> PyArray_DATA(runner._output_energy)

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

    storage.spectrum_virt_nu = <double*> PyArray_DATA(
        runner.legacy_montecarlo_virtual_luminosity)

    storage.sigma_thomson = (
        model.tardis_config.montecarlo.sigma_thomson.cgs.value)
    storage.inverse_sigma_thomson = 1.0 / storage.sigma_thomson
    storage.reflective_inner_boundary = model.tardis_config.montecarlo.enable_reflective_inner_boundary
    storage.inner_boundary_albedo = model.tardis_config.montecarlo.inner_boundary_albedo

    # Data for continuum implementation
    storage.t_electrons = <double*> PyArray_DATA(model.plasma_array.t_electrons)

    #cdef int no_levels_with_photdata
    cdef photo_xsect_1level ** photo_xsect

    if storage.cont_status == CONTINUUM_ON:
        # Photoionization data
        no_levels_with_photdata = model.atom_data.continuum_data.no_levels_with_photdata
        photo_xsect = <photo_xsect_1level **> malloc(no_levels_with_photdata * sizeof(photo_xsect_1level *))
        for i in range(no_levels_with_photdata):
            photo_xsect[i] = <photo_xsect_1level *> malloc(sizeof(photo_xsect_1level))
            phot_table_xsect = model.atom_data.continuum_data.get_phot_table_xsect(i)
            photo_xsect[i].no_of_points = len(phot_table_xsect)
            photo_xsect[i].x_sect = <double*> PyArray_DATA(phot_table_xsect)
            photo_xsect[i].nu = <double*> PyArray_DATA(model.atom_data.continuum_data.get_phot_table_nu(i))

        storage.photo_xsect = photo_xsect

        storage.chi_ff_factor = <double *> PyArray_DATA(model.base_continuum.free_free.chi_ff_factor)

        # Bound-free data
        storage.l_pop = <double*> PyArray_DATA(model.atom_data.continuum_data.level_number_density)
        storage.l_pop_r = <double*> PyArray_DATA(model.atom_data.continuum_data.level_number_density_ratio)
        storage.continuum_list_nu = <double*> PyArray_DATA(model.atom_data.continuum_data.continuum_edges_list)
        storage.no_of_edges = model.atom_data.continuum_data.continuum_edges_list.size
        storage.cont_edge2macro_continuum = <int_type_t*> PyArray_DATA(
            model.atom_data.continuum_data.cont_edge2macro_continuum)

        transition_probabilities_continuum = model.base_continuum.recombination_transition_probabilities
        storage.transition_probabilities_nd_continuum = transition_probabilities_continuum.data_array_nd
        storage.transition_probabilities_continuum = <double *> PyArray_DATA(
            transition_probabilities_continuum.data_array)
        storage.macro_block_references_continuum = <int_type_t*> PyArray_DATA(
            transition_probabilities_continuum.block_references)
        storage.transition_type_continuum = <int_type_t*> PyArray_DATA(
            transition_probabilities_continuum.dataframe['transition_type'].values)
        storage.transition_continuum_id = <int_type_t*> PyArray_DATA(
            transition_probabilities_continuum.dataframe['continuum_edge_idx'].values)
        storage.destination_level_id_continuum = <int_type_t*> PyArray_DATA(
            transition_probabilities_continuum.dataframe['destination_level_idx'].values)

        cooling_rates = model.base_continuum.cooling_rates

        for process_name in cooling_rates.cooling_processes:
            if process_name == 'free_free':
                storage.ff_cooling_prob = <double*> PyArray_DATA(cooling_rates.free_free_probability)
            else:
                cooling_data = getattr(cooling_rates, process_name)
                cooling_prob = <double*> PyArray_DATA(cooling_data.cooling_probability)
                cooling_prob_individual = <double *> PyArray_DATA(cooling_data.probabilities_array)
                references = <int_type_t*> PyArray_DATA(cooling_data.references)
                prob_array_nd = cooling_data.prob_array_nd
                if process_name == 'collisional_excitation':
                    storage.coll_exc_cooling_prob = cooling_prob
                    storage.coll_exc_cooling_prob_individual = cooling_prob_individual
                    storage.coll_exc_cooling_references = references
                    storage.coll_exc_cooling_prob_nd = prob_array_nd
                elif process_name == 'radiative_recombination':
                    storage.fb_cooling_prob = cooling_prob
                    storage.fb_cooling_prob_individual = cooling_prob_individual
                    storage.fb_cooling_references = references
                    storage.fb_cooling_prob_nd = prob_array_nd
                elif process_name == 'collisional_ionization':
                    storage.coll_ion_cooling_prob = cooling_prob
                    storage.coll_ion_cooling_prob_individual = cooling_prob_individual
                    storage.coll_ion_cooling_references = references
                    storage.coll_ion_cooling_prob_nd = prob_array_nd

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

    initialize_storage_model(model, runner, &storage)

    montecarlo_main_loop(&storage, virtual_packet_flag, nthreads,
                         model.tardis_config.montecarlo.seed)
    cdef np.ndarray[double, ndim=1] virt_packet_nus = np.zeros(storage.virt_packet_count, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] virt_packet_energies = np.zeros(storage.virt_packet_count, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] virt_packet_last_interaction_in_nu = np.zeros(storage.virt_packet_count, dtype=np.float64)
    cdef np.ndarray[int_type_t, ndim=1] virt_packet_last_interaction_type = np.zeros(storage.virt_packet_count, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] virt_packet_last_line_interaction_in_id = np.zeros(storage.virt_packet_count, dtype=np.int64)
    cdef np.ndarray[int_type_t, ndim=1] virt_packet_last_line_interaction_out_id = np.zeros(storage.virt_packet_count, dtype=np.int64)
    runner.virt_logging = LOG_VPACKETS
    if LOG_VPACKETS != 0:
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
        runner.virt_packet_nus = virt_packet_nus
        runner.virt_packet_energies = virt_packet_energies
        runner.virt_packet_last_interaction_in_nu = virt_packet_last_interaction_in_nu
        runner.virt_packet_last_interaction_type = virt_packet_last_interaction_type
        runner.virt_packet_last_line_interaction_in_id = virt_packet_last_line_interaction_in_id
        runner.virt_packet_last_line_interaction_out_id = virt_packet_last_line_interaction_out_id
    else:
        runner.virt_packet_nus = None
        runner.virt_packet_energies = None
        runner.virt_packet_last_interaction_in_nu = None
        runner.virt_packet_last_interaction_type = None
        runner.virt_packet_last_line_interaction_in_id = None
        runner.virt_packet_last_line_interaction_out_id = None
    #return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id, virt_packet_nus, virt_packet_energies

    # Necessary?
    if storage.cont_status == CONTINUUM_ON:
        for i in range(storage.no_of_edges):
            free(<void *> storage.photo_xsect[i])
