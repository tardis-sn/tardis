#include "cmontecarlo.h"
#include "stdlib.h"

int rpacket_get_nu_wrapper(int, int);
int rpacket_set_nu_wrapper(int);
int rpacket_get_mu_wrapper(int, int);
int rpacket_set_mu_wrapper(int);
int rpacket_get_energy_wrapper(int, int);
int rpacket_set_energy_wrapper(int);
int rpacket_get_nu_line_wrapper(int, int);
int rpacket_set_nu_line_wrapper(int);

int rpacket_get_nu_wrapper(int nu, int expected) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	packet->nu = nu;
	double actual_nu;	
 	actual_nu = rpacket_get_nu(packet);
	free(packet);
	return (int)actual_nu == expected;
}

int rpacket_set_nu_wrapper(int nu) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	rpacket_set_nu(packet, (double)nu);
	double actual_nu;
	actual_nu = rpacket_get_nu(packet);
	free(packet);
	return (int)actual_nu == nu;
}

int rpacket_get_mu_wrapper(int mu, int expected) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	packet->mu = mu;
	double actual_mu;	
 	actual_mu = rpacket_get_mu(packet);
	free(packet);
	return (int)actual_mu == expected;
}

int rpacket_set_mu_wrapper(int mu) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	rpacket_set_mu(packet, (double)mu);
	double actual_mu;
	actual_mu = rpacket_get_mu(packet);
	free(packet);
	return (int)actual_mu == mu;
}

int rpacket_get_energy_wrapper(int energy, int expected) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	packet->energy = energy;
	double actual_energy;	
 	actual_energy = rpacket_get_energy(packet);
	free(packet);
	return (int)actual_energy == expected;
}

int rpacket_set_energy_wrapper(int energy) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	rpacket_set_energy(packet, (double)energy);
	double actual_energy;
	actual_energy = rpacket_get_energy(packet);
	free(packet);
	return (int)actual_energy == energy;
}

int rpacket_get_nu_line_wrapper(int nu_line, int expected) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	packet->nu_line = nu_line;
	double actual_nu_line;	
 	actual_nu_line = rpacket_get_nu_line(packet);
	free(packet);
	return (int)actual_nu_line == expected;
}

int rpacket_set_nu_line_wrapper(int nu_line) {
	rpacket_t *packet;
	packet = (rpacket_t *) malloc(sizeof(packet));
	rpacket_set_nu_line(packet, (double)nu_line);
	double actual_nu_line;
	actual_nu_line = rpacket_get_nu_line(packet);
	return (int)actual_nu_line == nu_line;
}
