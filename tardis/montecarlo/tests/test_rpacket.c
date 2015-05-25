#include "../src/rpacket.h"
#include <assert.h>
#include <stdlib.h>

void test_rpacket_get_nu(){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	double value = 10.2;
	rpacket_set_nu(rp, value);
	assert(value == rpacket_get_nu(rp));
}

int main(){

	test_rpacket_get_nu();
	return 0;
}