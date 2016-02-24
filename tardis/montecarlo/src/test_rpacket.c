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
bool test_rpacket_get_id(void);

bool
test_rpacket_get_nu(double value){
	rpacket_t rp;
	rpacket_set_nu(&rp, value);
        return value==rpacket_get_nu(&rp);
}

bool
test_rpacket_get_mu(double value){
	rpacket_t rp;
	rpacket_set_mu(&rp, value);
        return value==rpacket_get_mu(&rp);
}

bool
test_rpacket_get_energy(double value){
	rpacket_t rp;
	rpacket_set_energy(&rp, value);
        return value==rpacket_get_energy(&rp);
}

bool
test_rpacket_get_r(double value){
	rpacket_t rp;
	rpacket_set_r(&rp, value);
        return value==rpacket_get_r(&rp);
}

bool
test_rpacket_get_tau_event(double value){
	rpacket_t rp;
	rpacket_set_tau_event(&rp, value);
        return value==rpacket_get_tau_event(&rp); 
}

bool
test_rpacket_get_nu_line(double value){
	rpacket_t rp;
	rpacket_set_nu_line(&rp, value);
        return value==rpacket_get_nu_line(&rp); 
}

bool
test_rpacket_get_current_shell_id(unsigned int value){
	rpacket_t rp;
	rpacket_set_current_shell_id(&rp, value);
	return value==rpacket_get_current_shell_id(&rp);
}

bool
test_rpacket_get_next_line_id(unsigned int value){
	rpacket_t rp;
	rpacket_set_next_line_id(&rp, value);
	return value==rpacket_get_next_line_id(&rp);
}

bool
test_rpacket_get_last_line(void){
	rpacket_t rp;
	rpacket_set_last_line(&rp, true);
	return rpacket_get_last_line(&rp);	
}

bool
test_rpacket_get_close_line(void){
	rpacket_t rp;
	rpacket_set_last_line(&rp, true);
	return rpacket_get_last_line(&rp);
}

bool
test_rpacket_get_recently_crossed_boundary(int value){
	rpacket_t rp;
	rpacket_set_recently_crossed_boundary(&rp, value);
	return value==rpacket_get_recently_crossed_boundary(&rp);
}

bool
test_rpacket_get_virtual_packet_flag(int value){
	rpacket_t rp;
	rpacket_set_virtual_packet_flag(&rp, value);
	return value==rpacket_get_virtual_packet_flag(&rp);
}

bool
test_rpacket_get_virtual_packet(int value){
	rpacket_t rp;
	rpacket_set_virtual_packet(&rp, value);
	return value==rpacket_get_virtual_packet(&rp);
}

bool
test_rpacket_get_d_boundary(double value){
	rpacket_t rp;
	rpacket_set_d_boundary(&rp, value);
	return value==rpacket_get_d_boundary(&rp);
}

bool
test_rpacket_get_d_electron(double value){
	rpacket_t rp;
	rpacket_set_d_electron(&rp, value);
	return value==rpacket_get_d_electron(&rp);
}

bool
test_rpacket_get_d_line(double value){
	rpacket_t rp;
	rpacket_set_d_line(&rp, value);
	return value==rpacket_get_d_line(&rp);
}

bool
test_rpacket_get_next_shell_id(int value){
	rpacket_t rp;
	rpacket_set_next_shell_id(&rp, value);
	return value==rpacket_get_next_shell_id(&rp);
}

bool
test_rpacket_get_status(void){
	rpacket_status_t inProcess = TARDIS_PACKET_STATUS_IN_PROCESS;
	rpacket_status_t emitted = TARDIS_PACKET_STATUS_EMITTED;
	rpacket_status_t reabsorbed = TARDIS_PACKET_STATUS_REABSORBED;

	rpacket_t rp;
	rpacket_set_status(&rp, inProcess);
	bool res=  inProcess==rpacket_get_status(&rp);
        rpacket_set_status(&rp, emitted);
	res &= emitted==rpacket_get_status(&rp);
	rpacket_set_status(&rp, reabsorbed);
        res &= reabsorbed==rpacket_get_status(&rp);
	return res;
}

bool
test_rpacket_get_id(void){
	rpacket_t rp;
	rpacket_set_id(&rp, 2);
	return rpacket_get_id(&rp) == 2;
}
