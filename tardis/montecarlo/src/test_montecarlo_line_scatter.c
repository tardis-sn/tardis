#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#include "cmontecarlo.h"


rpacket_t * rp;
storage_model_t * sm;

double TIME_EXPLOSION =     51840000.0;
double R_INNER_VALUE =  622080000000.0; /* 12,000xTIME_EXPLOSION */
double R_OUTER_VALUE = 7776600000000.0; /* 15,000xTIME_EXPLOSION */

void init_rpacket(void);
void init_storage_model(void);
double test_compute_distance2boundary(void);

/* initialise RPacket */
void
init_rpacket(void){
	rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	
	double MU = 0.3;
	double R = 673920000000.0; /* 13,000xTIME_EXPLOSION */
	double ENERGY = 0.9;

	rpacket_set_current_shell_id(rp, 0);
	rpacket_set_mu(rp, MU);
	rpacket_set_r(rp, R);
	rpacket_set_recently_crossed_boundary(rp, 1);

	/*
	rpacket_set_virtual_packet(rp, 1);
	rpacket_set_current_shell_id(rp, 1);
	rpacket_set_tau_event(rp, TAU_EVENT);
	rpacket_set_nu_line(rp, NU_LINE);
	rpacket_set_last_line(rp, true);
	rpacket_set_energy(rp, ENERGY);
	rpacket_set_status(rp, TARDIS_PACKET_STATUS_IN_PROCESS);
	rpacket_set_d_electron(rp, D_ELECTRON);
	rpacket_set_d_line(rp, D_LINE);
	rpacket_set_next_line_id(rp, NEXT_LINE_ID);
	*/

}

/* initialise storage model */
void
init_storage_model(void){
	sm = (storage_model_t *) malloc(sizeof(storage_model_t));
	double R_OUTER[2] = {R_OUTER_VALUE, R_OUTER_VALUE};
	double R_INNER[2] = {R_INNER_VALUE, R_INNER_VALUE};
	sm->r_outer = R_OUTER;
	sm->r_inner = R_INNER;
	
	/*
	double LINE_LIST_NU[2] = {1.5, 3.2};
	double LINE_LISTS_J_BLUES[2] = {5.6, 4.4};

	sm->no_of_packets = 2;
	sm->time_explosion = TIME_EXPLOSION;
	sm->inverse_time_explosion = 1.0/TIME_EXPLOSION;
	sm->r_outer = R_OUTER;
	sm->r_inner = R_INNER;
	sm->line_list_nu = LINE_LIST_NU;
	sm->line_lists_j_blues = LINE_LISTS_J_BLUES;
	*/
}

double
test_compute_distance2boundary(){
	/*
	* this returns nan
	*/
	double D_BOUNDARY = compute_distance2boundary(rp, sm);
	rpacket_set_d_boundary(rp, D_BOUNDARY);
	return D_BOUNDARY;
}