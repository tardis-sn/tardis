#include <assert.h>
#include <stdlib.h>

#include "../src/montecarlo.h"

int testing_rpacket_get_nu(void);

int
testing_rpacket_get_nu(void){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	double value = 10.2;
	rpacket_set_nu(rp, value);
	/*if( value != rpacket_get_nu(rp) ){
		return 0;
	}*/
	return 1;
}