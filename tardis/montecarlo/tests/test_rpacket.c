#include "../src/rpacket.h"
#include <assert.h>
#include <stdlib.h>
#include <Python.h>

void test_rpacket_get_nu(){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	double value = 10.2;
	rpacket_set_nu(rp, value);
	if( value != rpacket_get_nu(rp) ){
		PyErr_Format(PyExc_AssertionError, "value %g not %g", value, rpacket_get_nu(rp));
	}
}

int main(){

	test_rpacket_get_nu();
	return 0;
}