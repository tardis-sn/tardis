#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "../src/rpacket.h"

int testing_rpacket_get_nu();

int
testing_rpacket_get_nu(){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	double value = 10.2;
	/*rpacket_set_nu(rp, value);
	if( value != rpacket_get_nu(rp) ){
		return 0;
	}*/
	return true;
}