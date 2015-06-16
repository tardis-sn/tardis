#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#include "cmontecarlo.h"


rpacket_t * rp;
storage_model_t * sm;

bool test_montecarlo_line_scatter(void);

/* initialise RPacket */
void
init_rpacket(void){
	rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	
	double MU = 1.5*pow(10, 5);
	double R = 3.3*pow(10, 3);
	double TAU_EVENT = 3.2*pow(10, 4);
	double NU_LINE = 1.5;
	double ENERGY = 10.0;
	double D_BOUNDARY = 5.8;
	double D_ELECTRON = 4.3;
	double D_LINE = 3.9;
	unsigned int NEXT_LINE_ID = 1;

	rpacket_set_mu(rp, MU);
	rpacket_set_r(rp, R);
	rpacket_set_virtual_packet(rp, 1);
	rpacket_set_current_shell_id(rp, 1);
	rpacket_set_tau_event(rp, TAU_EVENT);
	rpacket_set_recently_crossed_boundary(rp, 1);
	rpacket_set_nu_line(rp, NU_LINE);
	rpacket_set_last_line(rp, true);
	rpacket_set_energy(rp, ENERGY);
	rpacket_set_status(rp, TARDIS_PACKET_STATUS_IN_PROCESS);
	rpacket_set_d_boundary(rp, D_BOUNDARY);
	rpacket_set_d_electron(rp, D_ELECTRON);
	rpacket_set_d_line(rp, D_LINE);
	rpacket_set_next_line_id(rp, NEXT_LINE_ID);
}

/* initialise storage model */
void
init_storage_model(void){
	sm = (storage_model_t *) malloc(sizeof(storage_model_t));
	double R_OUTER[2] = {1.2, 2.6};
	double R_INNER[2] = {2.3, 3.4};
	double TIME_EXPLOSION = 36.0;
	double LINE_LIST_NU[2] = {1.5, 3.2};
	double LINE_LISTS_J_BLUES[2] = {5.6, 4.4};

	sm->no_of_packets = 2;
	sm->time_explosion = TIME_EXPLOSION;
	sm->inverse_time_explosion = 1.0/TIME_EXPLOSION;
	sm->r_outer = R_OUTER;
	sm->r_inner = R_INNER;
	sm->line_list_nu = LINE_LIST_NU;
	sm->line_lists_j_blues = LINE_LISTS_J_BLUES;
}

bool
test_montecarlo_line_scatter(){
	double DISTANCE = 10.5;
	montecarlo_line_scatter(rp, sm, DISTANCE);

	/* Return true if the test passes */
	return true;
}