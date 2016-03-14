#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cmontecarlo.h"


rpacket_t * rp;
storage_model_t * sm;

double TIME_EXPLOSION =  5.2e7; /* 10 days(in seconds)   ~      51840000.0 */
double R_INNER_VALUE =  6.2e11; /* 12,000xTIME_EXPLOSION ~  622080000000.0 */
double R_OUTER_VALUE =  7.8e12; /* 15,000xTIME_EXPLOSION ~ 7776600000000.0 */
double SIGMA_THOMSON = 6.652486e-25;

double test_compute_distance2boundary(void);
double test_compute_distance2line(void);
double test_compute_distance2continuum(void);
double test_rpacket_doppler_factor(void);
//double test_move_packet(void);
bool test_move_packet_across_shell_boundary(void);
int64_t test_montecarlo_one_packet(void);
int64_t test_montecarlo_one_packet_loop(void);
bool test_montecarlo_line_scatter(void);
double test_increment_j_blue_estimator(void);
bool test_montecarlo_thomson_scatter(void);
//bool test_macro_atom(void);
double test_calculate_chi_bf(void);
bool test_montecarlo_bound_free_scatter(void);
double test_bf_cross_section(void);
int64_t test_montecarlo_free_free_scatter(void);

/* initialise RPacket */
static void init_rpacket(rpacket_t *rp){
	double MU = 0.3;
	double R = 7.5e14;
	double ENERGY = 0.9;
	int NEXT_LINE_ID = 1;
	double NU = 0.4;
	double NU_LINE = 0.2;
	int CURRENT_SHELL_ID = 0;

	double TAU_EVENT = 2.9e13;

	rpacket_set_current_shell_id(rp, CURRENT_SHELL_ID);
	rpacket_set_next_shell_id(rp, CURRENT_SHELL_ID+1);
	rpacket_set_mu(rp, MU);
	rpacket_set_nu(rp, NU);
	rpacket_set_r(rp, R);
	rpacket_set_last_line(rp, false);
	rpacket_set_close_line(rp, false);
	rpacket_set_nu_line(rp, NU_LINE);

	rpacket_set_next_line_id(rp, NEXT_LINE_ID);

	rpacket_set_tau_event(rp, TAU_EVENT);
	rpacket_set_virtual_packet(rp, 0);
	rpacket_set_energy(rp, ENERGY);
	rpacket_set_virtual_packet_flag(rp, true);
	rpacket_set_status(rp, TARDIS_PACKET_STATUS_IN_PROCESS);
	rpacket_set_id(rp, 0);

	rpacket_set_current_continuum_id(rp, 1);
}

/* initialise storage model */
static void init_storage_model(storage_model_t *sm){
	int NUMBER_OF_SHELLS = 2;

	sm->time_explosion = TIME_EXPLOSION;
	sm->inverse_time_explosion = 1.0/TIME_EXPLOSION;
	sm->inverse_sigma_thomson = 1.0/SIGMA_THOMSON;

	/* R_OUTER = {8.64e14, 1.0368e15} */
	sm->r_outer = (double *) malloc(sizeof(double)*NUMBER_OF_SHELLS);
	sm->r_outer[0] = 8.64e14;
	sm->r_outer[1] = 1.0368e15;

	/* R_INNER = {6.912e14, 8.64e14} */
	sm->r_inner = (double *) malloc(sizeof(double)*NUMBER_OF_SHELLS);
	sm->r_inner[0] = 6.912e14;
	sm->r_inner[1] = 8.64e14;
	sm->no_of_lines = 2;

	/* LINE_LIST_NU = { 1.26318289e+16,  1.26318289e+16,
						1.23357675e+16,  1.23357675e+16,
						1.16961598e+16 };
	*/
	sm->line_list_nu = (double *) malloc(sizeof(double)*5);
	sm->line_list_nu[0] = 1.26318289e+16;
	sm->line_list_nu[1] = 1.26318289e+16;
	sm->line_list_nu[2] = 1.23357675e+16;
	sm->line_list_nu[3] = 1.23357675e+16;
	sm->line_list_nu[4] = 1.16961598e+16;

	/* INVERSE_ELECTRON_DENSITIES = {} */
	sm->inverse_electron_densities = (double *) malloc(sizeof(double)*NUMBER_OF_SHELLS);
	sm->inverse_electron_densities[0] = 1e-9;
	sm->inverse_electron_densities[1] = 1e-9;
	sm->electron_densities = (double *) malloc(sizeof(double )*2);
	sm->electron_densities[0] = 1e9;
	sm->electron_densities[1] = 1e9;


	sm->js = (double *) malloc(sizeof(double)*NUMBER_OF_SHELLS);
	sm->js[0] = 0;
	sm->js[1] = 0;

	sm->nubars = (double *) malloc(sizeof(double)*NUMBER_OF_SHELLS);
	sm->nubars[0] = 0;
	sm->nubars[1] = 0;

	sm->last_line_interaction_in_id = (int64_t *) malloc(sizeof(int64_t)*NUMBER_OF_SHELLS);
	sm->last_line_interaction_in_id[0] = 0;
	sm->last_line_interaction_in_id[1] = 0;

	sm->last_line_interaction_shell_id = (int64_t *) malloc(sizeof(int64_t)*NUMBER_OF_SHELLS);
	sm->last_line_interaction_shell_id[0] = 0;
	sm->last_line_interaction_shell_id[1] = 0;

	sm->last_interaction_type = (int64_t *) malloc(sizeof(int64_t)*1);
	sm->last_interaction_type[0] = 2;

	sm->line_interaction_id = 0;

	sm->line_lists_j_blues_nd = 0;

	sm->line_lists_j_blues = (double *) malloc(sizeof(double )*2);
	sm->line_lists_j_blues[0] = 1e-10;
	sm->line_lists_j_blues[1] = 1e-10;

	sm->line_lists_tau_sobolevs = (double *) malloc(sizeof(double )*2);
	sm->line_lists_tau_sobolevs[0] = 1e-5;
	sm->line_lists_tau_sobolevs[1] = 1000;

	sm->line2macro_level_upper = (int64_t *) malloc(sizeof(int64_t)*2);
	sm->line2macro_level_upper[0] = 0;
	sm->line2macro_level_upper[1] = 0;

	sm->reflective_inner_boundary = false;
	sm->inner_boundary_albedo = 0.0;
	sm->no_of_shells = NUMBER_OF_SHELLS;

	sm->spectrum_start_nu = 1.e14;
	sm->spectrum_delta_nu = 293796608840.0;
	sm->spectrum_end_nu = 6.e15;
/* FIXME: copy paste error below?
	sm->spectrum_virt_start_nu = 1.e14;
	sm->spectrum_virt_end_nu = 293796608840.0;
	sm->spectrum_delta_nu = 6.e15;
    */

	sm->spectrum_virt_start_nu = 1.e14;
	sm->spectrum_virt_end_nu = 6.e15;

	sm->spectrum_virt_nu = (double *) calloc(20000,sizeof(double ));

	/*
	*  Initialising the below values to 0 untill
	*  I get the real data !
	*/

	sm->t_electrons = (double *) malloc(sizeof(double )*2);
	sm->t_electrons[0] = 2;
	sm->t_electrons[1] = 2;

	sm->l_pop = (double *) malloc(sizeof(double )*20000);
        for (size_t i=0; i<20000; ++i)
          sm->l_pop[i]=2;

	sm->l_pop_r = (double *) malloc(sizeof(double )* 20000);
        for (size_t i=0; i<20000; ++i)
          sm->l_pop_r[i]=3;

	sm->continuum_list_nu = (double *) malloc(sizeof(double )* 20000);
        for (size_t i=0; i<20000; ++i)
          sm->continuum_list_nu[i]=1e13;

	sm->no_of_edges = 100;

	sm->chi_bf_tmp_partial = (double *) malloc(sizeof(double )* 20000);
        for (size_t i=0; i<20000; ++i)
          sm->chi_bf_tmp_partial[i]=160;

}
static void dealloc_storage_model(storage_model_t *sm){
        free(sm->r_outer);
        free(sm->r_inner);
        free(sm->line_list_nu);
        free(sm->inverse_electron_densities);
        free(sm->electron_densities);
        free(sm->js);
        free(sm->nubars);
        free(sm->last_line_interaction_in_id);
        free(sm->last_line_interaction_shell_id);
        free(sm->last_interaction_type);
        free(sm->line_lists_j_blues);
        free(sm->line_lists_tau_sobolevs);
        free(sm->line2macro_level_upper);
        free(sm->spectrum_virt_nu);
        free(sm->t_electrons);
        free(sm->l_pop);
        free(sm->l_pop_r);
        free(sm->continuum_list_nu);
        free(sm->chi_bf_tmp_partial);
}

static void irandom(rk_state *mt_state)
{
memset(mt_state,0,sizeof(rk_state));
}

double
test_compute_distance2boundary(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        compute_distance2boundary(&rp, &sm);
        double D_BOUNDARY = rpacket_get_d_boundary(&rp);
        dealloc_storage_model(&sm);
        return D_BOUNDARY;
}
double
test_compute_distance2line(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        // FIXME MR: return status of compute_distance2line() is ignored
        compute_distance2line(&rp, &sm);
        double D_LINE = rpacket_get_d_line(&rp);
        dealloc_storage_model(&sm);
        return D_LINE;
}

double
test_compute_distance2continuum(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
	compute_distance2continuum(&rp, &sm);
        dealloc_storage_model(&sm);
	return rp.d_cont;
}

double
test_rpacket_doppler_factor(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
	double res= rpacket_doppler_factor(&rp, &sm);
        dealloc_storage_model(&sm);
        return res;
}

//double
//test_move_packet(void){
//	double DISTANCE = 1e13;
//        rpacket_t rp;
//        storage_model_t sm;
//        init_rpacket(&rp);
//        init_storage_model(&sm);
//	double res = move_packet(&rp, &sm, DISTANCE);
//        dealloc_storage_model(&sm);
//        return res;
//}

double
test_increment_j_blue_estimator(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        int64_t j_blue_idx = 0;
        compute_distance2line(&rp, &sm);
        move_packet(&rp, &sm, 1e13);
        double d_line = rpacket_get_d_line(&rp);
        increment_j_blue_estimator(&rp, &sm, d_line, j_blue_idx);
        double res = sm.line_lists_j_blues[j_blue_idx];
        dealloc_storage_model(&sm);
        return res;
}

bool
test_montecarlo_line_scatter(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
	double DISTANCE = 1e13;
	montecarlo_line_scatter(&rp, &sm, DISTANCE, &mt_state);
        dealloc_storage_model(&sm);
	return true;
}

bool
test_montecarlo_thomson_scatter(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
	double DISTANCE = 1e13;
	montecarlo_thomson_scatter(&rp, &sm, DISTANCE, &mt_state);
        dealloc_storage_model(&sm);
	return true;
}

bool
test_move_packet_across_shell_boundary(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
	double DISTANCE = 0.95e13;
// MR: wrong: move_packet_across_shell_boundary() returns void
        move_packet_across_shell_boundary(&rp, &sm, DISTANCE, &mt_state);
        dealloc_storage_model(&sm);
        return true;
}


int64_t
test_montecarlo_one_packet(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
        int64_t res = montecarlo_one_packet(&sm, &rp, 1, &mt_state);
        dealloc_storage_model(&sm);
	return res;
}

int64_t
test_montecarlo_one_packet_loop(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
	int64_t res= montecarlo_one_packet_loop(&sm, &rp, 1, &mt_state);
        dealloc_storage_model(&sm);
        return res;
}

//bool
//test_macro_atom(void){
//	macro_atom(rp, sm);
//	return true;
//}

double
test_calculate_chi_bf(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
        int64_t j_blue_idx = 0;
        compute_distance2line(&rp, &sm);
        move_packet(&rp, &sm, 1e13);
        double d_line = rpacket_get_d_line(&rp);
        increment_j_blue_estimator(&rp, &sm, d_line, j_blue_idx);
        double DISTANCE = 1e13;
        montecarlo_line_scatter(&rp, &sm, DISTANCE, &mt_state);
        DISTANCE = 0.95e13;
// MR: wrong: move_packet_across_shell_boundary() returns void
        move_packet_across_shell_boundary(&rp, &sm, DISTANCE, &mt_state);
        montecarlo_one_packet(&sm, &rp, 1, &mt_state);
        montecarlo_one_packet_loop(&sm, &rp, 1, &mt_state);
        DISTANCE = 1e13;
        montecarlo_thomson_scatter(&rp, &sm, DISTANCE, &mt_state);
        calculate_chi_bf(&rp, &sm);
	double res = rpacket_doppler_factor (&rp, &sm);
        dealloc_storage_model(&sm);
        return res;
}

bool
test_montecarlo_bound_free_scatter(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
	double DISTANCE = 1e13;
	montecarlo_bound_free_scatter(&rp, &sm, DISTANCE, &mt_state);
        dealloc_storage_model(&sm);
	return (rpacket_get_status(&rp)!=0);
}

double
test_bf_cross_section(void){
        storage_model_t sm;
        init_storage_model(&sm);
	double CONV_MU = 0.4;
	double res = bf_cross_section(&sm, 1, CONV_MU);
        dealloc_storage_model(&sm);
        return res;
}

int64_t
test_montecarlo_free_free_scatter(void){
        rpacket_t rp;
        storage_model_t sm;
        init_rpacket(&rp);
        init_storage_model(&sm);
        rk_state mt_state;
        irandom(&mt_state);
	double DISTANCE = 1e13;
	montecarlo_free_free_scatter(&rp, &sm, DISTANCE,&mt_state);
        dealloc_storage_model(&sm);
	return rpacket_get_status(&rp);
}
