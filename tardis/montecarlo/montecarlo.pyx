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

np.import_array()



ctypedef np.int64_t int_type_t

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef c_array_to_numpy(void *ptr, int dtype, np.npy_intp N):
    cdef np.ndarray arr = np.PyArray_SimpleNewFromData(1, &N, dtype, ptr)
    PyArray_ENABLEFLAGS(arr, np.NPY_OWNDATA)
    return arr

cdef extern from "src/cmontecarlo.h":
    ctypedef enum ContinuumProcessesStatus:
        CONTINUUM_OFF = 0
        CONTINUUM_ON = 1

    cdef int LOG_VPACKETS

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
        int_type_t *last_interaction_out_type
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
        double *line_lists_Edotlu
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
        photo_xsect_1level ** photo_xsect
        double *chi_ff_factor
        double *t_electrons
        double *l_pop
        double *l_pop_r
        ContinuumProcessesStatus cont_status
        double *virt_packet_nus
        double *virt_packet_energies
        double *virt_packet_last_interaction_in_nu
        int_type_t *virt_packet_last_interaction_type
        int_type_t *virt_packet_last_line_interaction_in_id
        int_type_t *virt_packet_last_line_interaction_out_id
        int_type_t virt_packet_count
        int_type_t virt_array_size
        int_type_t kpacket2macro_level
        int_type_t *cont_edge2macro_level
        double *photo_ion_estimator
        double *stim_recomb_estimator
        int_type_t *photo_ion_estimator_statistics
        double *bf_heating_estimator
        double *ff_heating_estimator
        double *stim_recomb_cooling_estimator

    void montecarlo_main_loop(storage_model_t * storage, int_type_t virtual_packet_flag, int nthreads, unsigned long seed)




cdef initialize_storage_model(model, plasma, runner, storage_model_t *storage):
    """
    Initializing the storage struct.

    """

    storage.no_of_packets = runner.input_nu.size
    storage.packet_nus = <double*> PyArray_DATA(runner.input_nu)
    storage.packet_mus = <double*> PyArray_DATA(runner.input_mu)
    storage.packet_energies = <double*> PyArray_DATA(runner.input_energy)

    # Setup of structure
    storage.no_of_shells = model.no_of_shells


    storage.r_inner = <double*> PyArray_DATA(runner.r_inner_cgs)
    storage.r_outer = <double*> PyArray_DATA(runner.r_outer_cgs)
    storage.v_inner = <double*> PyArray_DATA(runner.v_inner_cgs)

    # Setup the rest
    # times
    storage.time_explosion = model.time_explosion.to('s').value
    storage.inverse_time_explosion = 1.0 / storage.time_explosion
    #electron density
    storage.electron_densities = <double*> PyArray_DATA(
        plasma.electron_densities.values)

    runner.inverse_electron_densities = (
        1.0 / plasma.electron_densities.values)
    storage.inverse_electron_densities = <double*> PyArray_DATA(
        runner.inverse_electron_densities)
    # Switch for continuum processes
    storage.cont_status = CONTINUUM_OFF
    # Continuum data
    cdef np.ndarray[double, ndim=1] continuum_list_nu
    cdef np.ndarray[double, ndim=1] l_pop
    cdef np.ndarray[double, ndim=1] l_pop_r

    if storage.cont_status == CONTINUUM_ON:
        continuum_list_nu = np.array([9.0e14, 8.223e14, 6.0e14, 3.5e14, 3.0e14])  # sorted list of threshold frequencies
        storage.continuum_list_nu = <double*> continuum_list_nu.data
        storage.no_of_edges = continuum_list_nu.size
        l_pop = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)
        storage.l_pop = <double*> l_pop.data
        l_pop_r = np.ones(storage.no_of_shells * continuum_list_nu.size, dtype=np.float64)
        storage.l_pop_r = <double*> l_pop_r.data

    # Line lists
    storage.no_of_lines = plasma.atomic_data.lines.nu.values.size
    storage.line_list_nu = <double*> PyArray_DATA(plasma.atomic_data.lines.nu.values)
    runner.line_lists_tau_sobolevs = (
            plasma.tau_sobolevs.values.flatten(order='F')
            )
    storage.line_lists_tau_sobolevs = <double*> PyArray_DATA(
            runner.line_lists_tau_sobolevs
            )
    storage.line_lists_j_blues = <double*> PyArray_DATA(
            runner.j_blue_estimator)

    storage.line_lists_Edotlu = <double*> PyArray_DATA(
            runner.Edotlu_estimator)

    storage.line_interaction_id = runner.get_line_interaction_id(
        runner.line_interaction_type)

    # macro atom & downbranch
    if storage.line_interaction_id >= 1:
        runner.transition_probabilities = (
                plasma.transition_probabilities.values.flatten(order='F')
        )
        storage.transition_probabilities = <double*> PyArray_DATA(
                runner.transition_probabilities
                )
        storage.transition_probabilities_nd = (
        plasma.transition_probabilities.values.shape[0])
        storage.line2macro_level_upper = <int_type_t*> PyArray_DATA(
            plasma.atomic_data.lines_upper2macro_reference_idx)
        storage.macro_block_references = <int_type_t*> PyArray_DATA(
            plasma.atomic_data.macro_atom_references['block_references'].values)
        storage.transition_type = <int_type_t*> PyArray_DATA(
            plasma.atomic_data.macro_atom_data['transition_type'].values)

        # Destination level is not needed and/or generated for downbranch
        storage.destination_level_id = <int_type_t*> PyArray_DATA(
            plasma.atomic_data.macro_atom_data['destination_level_idx'].values)
        storage.transition_line_id = <int_type_t*> PyArray_DATA(
            plasma.atomic_data.macro_atom_data['lines_idx'].values)

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

    storage.spectrum_start_nu = runner.spectrum_frequency.to('Hz').value.min()
    storage.spectrum_end_nu = runner.spectrum_frequency.to('Hz').value.max()
    # TODO: Linspace handling for virtual_spectrum_range
    storage.spectrum_virt_start_nu = runner.virtual_spectrum_range.stop.to('Hz', units.spectral()).value
    storage.spectrum_virt_end_nu = runner.virtual_spectrum_range.start.to('Hz', units.spectral()).value
    storage.spectrum_delta_nu = runner.spectrum_frequency.to('Hz').value[1] - runner.spectrum_frequency.to('Hz').value[0]

    storage.spectrum_virt_nu = <double*> PyArray_DATA(
        runner.legacy_montecarlo_virtual_luminosity)

    storage.sigma_thomson = runner.sigma_thomson.cgs.value
    storage.inverse_sigma_thomson = 1.0 / storage.sigma_thomson
    storage.reflective_inner_boundary = runner.enable_reflective_inner_boundary
    storage.inner_boundary_albedo = runner.inner_boundary_albedo
    # Data for continuum implementation
    cdef np.ndarray[double, ndim=1] t_electrons = plasma.t_electrons
    storage.t_electrons = <double*> t_electrons.data

def montecarlo_radial1d(model, plasma, runner, int_type_t virtual_packet_flag=0,
                        int nthreads=4,last_run=False):
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

    initialize_storage_model(model, plasma, runner, &storage)

    montecarlo_main_loop(&storage, virtual_packet_flag, nthreads, runner.seed)
    runner.virt_logging = LOG_VPACKETS
    if LOG_VPACKETS != 0:
        runner.virt_packet_nus = c_array_to_numpy(storage.virt_packet_nus, np.NPY_DOUBLE, storage.virt_packet_count)
        runner.virt_packet_energies = c_array_to_numpy(storage.virt_packet_energies, np.NPY_DOUBLE, storage.virt_packet_count)
        runner.virt_packet_last_interaction_in_nu = c_array_to_numpy(storage.virt_packet_last_interaction_in_nu, np.NPY_DOUBLE, storage.virt_packet_count)
        runner.virt_packet_last_interaction_type = c_array_to_numpy(storage.virt_packet_last_interaction_type, np.NPY_INT64, storage.virt_packet_count)
        runner.virt_packet_last_line_interaction_in_id = c_array_to_numpy(storage.virt_packet_last_line_interaction_in_id, np.NPY_INT64,
                                                                          storage.virt_packet_count)
        runner.virt_packet_last_line_interaction_out_id = c_array_to_numpy(storage.virt_packet_last_line_interaction_out_id, np.NPY_INT64,
                                                                           storage.virt_packet_count)
    else:
        runner.virt_packet_nus = np.zeros(0)
        runner.virt_packet_energies = np.zeros(0)
        runner.virt_packet_last_interaction_in_nu = np.zeros(0)
        runner.virt_packet_last_interaction_type = np.zeros(0)
        runner.virt_packet_last_line_interaction_in_id = np.zeros(0)
        runner.virt_packet_last_line_interaction_out_id = np.zeros(0)
