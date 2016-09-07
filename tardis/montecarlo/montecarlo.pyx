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
from tardis.util import intensity_black_body

np.import_array()



ctypedef np.int64_t int_type_t
ctypedef np.float64_t float_t

cdef extern from "src/cmontecarlo.h":
    ctypedef enum ContinuumProcessesStatus:
        CONTINUUM_OFF = 0
        CONTINUUM_ON = 1

    cdef int LOG_VPACKETS

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
        double *chi_bf_tmp_partial
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


cdef extern from "src/integrator.h":
    void montecarlo_main_loop(storage_model_t * storage, int_type_t virtual_packet_flag, int nthreads, unsigned long seed)
    void debug_print_arg(double* arg, int len)
    void debug_print_2d_arg(double* arg, int len1, int len2)
    void integrate_source_functions(double* L_nu, double* line_nu, double* taus, double* att_S_ul, double* I_BB, double* nus, 
              double* ps, double* Rs, double R_ph, double inv_ct, int_type_t* lens)

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
        model.plasma.electron_densities.values)

    runner.inverse_electron_densities = (
        1.0 / model.plasma.electron_densities.values)
    storage.inverse_electron_densities = <double*> PyArray_DATA(
        runner.inverse_electron_densities)
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
    storage.no_of_lines = model.atom_data.lines.nu.values.size
    storage.line_list_nu = <double*> PyArray_DATA(model.atom_data.lines.nu.values)
    runner.line_lists_tau_sobolevs = (
            model.plasma.tau_sobolevs.values.flatten(order='F')
            )
    storage.line_lists_tau_sobolevs = <double*> PyArray_DATA(
            runner.line_lists_tau_sobolevs
            )
    storage.line_lists_j_blues = <double*> PyArray_DATA(
            runner.j_blue_estimator)

    storage.line_lists_Edotlu = <double*> PyArray_DATA(
            runner.Edotlu_estimator)

    storage.line_interaction_id = runner.get_line_interaction_id(
        model.tardis_config.plasma.line_interaction_type)

    # macro atom & downbranch
    if storage.line_interaction_id >= 1:
        runner.transition_probabilities = (
                model.plasma.transition_probabilities.values.flatten(order='F')
        )
        storage.transition_probabilities = <double*> PyArray_DATA(
                runner.transition_probabilities
                )
        storage.transition_probabilities_nd = (
        model.plasma.transition_probabilities.values.shape[0])
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
    cdef np.ndarray[double, ndim=1] t_electrons = model.plasma.t_electrons
    storage.t_electrons = <double*> t_electrons.data

def montecarlo_radial1d(model, runner, int_type_t virtual_packet_flag=0,
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
    #return output_nus, output_energies, js, nubars, last_line_interaction_in_id, last_line_interaction_out_id, last_interaction_type, last_line_interaction_shell_id, virt_packet_nus, virt_packet_energies
    else:
        runner.virt_packet_nus = np.zeros(0)
        runner.virt_packet_energies = np.zeros(0)
        runner.virt_packet_last_interaction_in_nu = np.zeros(0)
        runner.virt_packet_last_interaction_type = np.zeros(0)
        runner.virt_packet_last_line_interaction_in_id = np.zeros(0)
        runner.virt_packet_last_line_interaction_out_id = np.zeros(0)

def integrate(np.ndarray[float_t] L_nu,np.ndarray[float_t] line_nu, np.ndarray[float_t, ndim=2] taus, 
              np.ndarray[float_t, ndim=2] att_S_ul, float R_max, float T, np.ndarray[float_t] nus, 
              np.ndarray[float_t] ps_outer, np.ndarray[float_t] ps_inner,
              np.ndarray[float_t, ndim=2] z_ct_outer, np.ndarray[float_t, ndim=2] z_ct_inner,
              int num_shell, np.ndarray[int_type_t] n_shell_p_outer):   

    cdef int k, nu_idx, p_idx,  z_len
    cdef float nu_start, nu_end, p,
    cdef np.ndarray I_outer, I_inner, z_cross_p, shell_idx, ks, ps, I_p

    for nu_idx,nu in enumerate(nus.value):
        I_outer = np.zeros(len(ps_outer))
        for p_idx,p in enumerate(ps_outer):
            z_cross_p = z_ct_outer[z_ct_outer[:,p_idx] > 0,p_idx]
            if len(z_cross_p) == 1:
                z_cross_p = np.hstack((-z_cross_p,z_cross_p))
            else:
                z_cross_p = np.hstack((-z_cross_p,z_cross_p[::-1][1:],0)) # Zero ensures empty ks in last step below
                                                                          # 1: avoids double counting center
            shell_idx = (num_shell-1) - np.arange(n_shell_p_outer[p_idx]) # -1 for 0-based indexing
            shell_idx = np.hstack((shell_idx,shell_idx[::-1][1:]))
            
            z_len = len(z_cross_p)
            for idx,z_cross in enumerate(z_cross_p[0:z_len-1]):
                if z_cross_p[idx+1] == 0:
                    continue
                nu_start = nu * (1 - z_cross) 
                nu_end   = nu * (1 - z_cross_p[idx+1])
                shell = shell_idx[idx]
                # Note the direction of the comparisons
                ks, = np.where( (line_nu < nu_start) & (line_nu >= nu_end) )

                if len(ks) < 2:
                    ks = np.array(ks)

                for k in ks:
                    I_outer[p_idx] = I_outer[p_idx] * np.exp(-taus[k,shell]) + att_S_ul[k,shell]

        I_inner = np.zeros(len(ps_inner))
        for p_idx,p in enumerate(ps_inner):
            z_cross_p = z_ct_inner[z_ct_inner[:,p_idx] > 0,p_idx]
            z_cross_p = np.hstack((z_cross_p[::-1],0)) # Zero ensures empty ks in last step below

            shell_idx = np.hstack(( np.arange(num_shell), 0 ))
            I_inner[p_idx] = intensity_black_body(nu,T)
            z_len = len(z_cross_p)

            for idx,z_cross in enumerate(z_cross_p[0:z_len-1]):
                if z_cross_p[idx+1] == 0:
                    continue
                nu_start = nu * (1 - z_cross) 
                nu_end   = nu * (1 - z_cross_p[idx+1])
                shell = shell_idx[idx]
                # Note the direction of the comparisons
                ks, = np.where( (line_nu < nu_start) & (line_nu >= nu_end) )

                if len(ks) < 2:
                    ks = np.array(ks)

                for k in ks:
                    I_inner[p_idx] = I_inner[p_idx] * np.exp(-taus[k,shell]) + att_S_ul[k,shell]


        if ( nu_idx % 200 ) == 0:
            print "{:3.0f} %".format( 100*float(nu_idx)/len(nus))
            print I_outer, I_inner
        ps = np.hstack((ps_outer,ps_inner))*R_max
        I_p = np.hstack((I_outer,I_inner))*ps
        L_nu[nu_idx] = 8 * np.pi**2 *  np.trapz(y = I_p[::-1],x = ps[::-1])
    return  L_nu


def print_c_version(arg,twod):
    if twod is not None:
        debug_print_2d_arg(<double*> PyArray_DATA(arg),  arg.shape[0], arg.shape[1])
    else:
        debug_print_arg(<double*> PyArray_DATA(arg), len(arg))

def c_source_integrate(np.ndarray[float_t] L_nu, np.ndarray[float_t] line_nu, np.ndarray[float_t, ndim=2] taus, 
              np.ndarray[float_t, ndim=2] att_S_ul, np.ndarray[float_t] I_BB, np.ndarray[float_t] nus,
              np.ndarray[float_t] ps, np.ndarray[float_t] Rs, int num_shell, double R_ph, double inv_ct):
              
    lens = np.array([0,0,0,0])
    lens[0] = len(L_nu)
    lens[1] = len(line_nu)
    lens[2] = len(ps)
    lens[3] = num_shell

    integrate_source_functions(<double*> PyArray_DATA(L_nu), <double*> PyArray_DATA(line_nu), <double*> PyArray_DATA(taus),
            <double*> PyArray_DATA(att_S_ul), <double*> PyArray_DATA(I_BB), <double*> PyArray_DATA(nus), <double*> PyArray_DATA(ps), 
            <double*> PyArray_DATA(Rs), <double>R_ph, <double>inv_ct, <int_type_t*> PyArray_DATA(lens))
    return L_nu

