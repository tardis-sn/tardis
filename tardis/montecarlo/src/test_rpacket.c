#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "rpacket.h"
#include "status.h"

bool test_rpacket_get_nu(double);
bool test_rpacket_get_mu(double);
bool test_rpacket_get_energy(double);
bool test_rpacket_get_r(double);
bool test_rpacket_get_tau_event(double);
bool test_rpacket_get_nu_line(double);
bool test_rpacket_get_d_boundary(double);
bool test_rpacket_get_d_electron(double);
bool test_rpacket_get_d_line(double);


bool test_rpacket_get_current_shell_id(unsigned int);
bool test_rpacket_get_next_line_id(unsigned int);

bool test_rpacket_get_recently_crossed_boundary(int);
bool test_rpacket_get_virtual_packet_flag(int);
bool test_rpacket_get_virtual_packet(int);
bool test_rpacket_get_next_shell_id(int);

bool test_rpacket_get_last_line(void);
bool test_rpacket_get_close_line(void);
bool test_rpacket_get_status(void);

bool
test_rpacket_get_nu(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_nu(rp, value);
	if( value != rpacket_get_nu(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_mu(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_mu(rp, value);
	if( value != rpacket_get_mu(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_energy(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_energy(rp, value);
	if( value != rpacket_get_energy(rp) ){
		return false;
	}
	return true;	
}

bool
test_rpacket_get_r(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_r(rp, value);
	if( value != rpacket_get_r(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_tau_event(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_tau_event(rp, value);
	if( value != rpacket_get_tau_event(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_nu_line(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_nu_line(rp, value);
	if( value != rpacket_get_nu_line(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_current_shell_id(unsigned int value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_nu_line(rp, value);
	if( value != rpacket_get_nu_line(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_next_line_id(unsigned int value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_next_line_id(rp, value);
	if( value != rpacket_get_next_line_id(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_last_line(void){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_last_line(rp, true);
	return rpacket_get_last_line(rp);	
}

bool
test_rpacket_get_close_line(void){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_last_line(rp, true);
	return rpacket_get_last_line(rp);
}

bool
test_rpacket_get_recently_crossed_boundary(int value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_recently_crossed_boundary(rp, value);
	if( value != rpacket_get_recently_crossed_boundary(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_virtual_packet_flag(int value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_virtual_packet_flag(rp, value);
	if( value != rpacket_get_virtual_packet_flag(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_virtual_packet(int value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_virtual_packet(rp, value);
	if( value != rpacket_get_virtual_packet(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_d_boundary(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_d_boundary(rp, value);
	if( value != rpacket_get_d_boundary(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_d_electron(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_d_electron(rp, value);
	if( value != rpacket_get_d_electron(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_d_line(double value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_d_line(rp, value);
	if( value != rpacket_get_d_line(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_next_shell_id(int value){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_next_shell_id(rp, value);
	if( value != rpacket_get_next_shell_id(rp) ){
		return false;
	}
	return true;
}

bool
test_rpacket_get_status(void){
	rpacket_status_t inProcess = TARDIS_PACKET_STATUS_IN_PROCESS;
	rpacket_status_t emitted = TARDIS_PACKET_STATUS_EMITTED;
	rpacket_status_t reabsorbed = TARDIS_PACKET_STATUS_REABSORBED;

	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	rpacket_set_status(rp, inProcess);
	if( inProcess != rpacket_get_status(rp) ){
		return false;
	}
	rpacket_set_status(rp, emitted);
	if( emitted != rpacket_get_status(rp) ){
		return false;
	}
	rpacket_set_status(rp, reabsorbed);
	if( reabsorbed != rpacket_get_status(rp) ){
		return false;
	}
	return true;
}