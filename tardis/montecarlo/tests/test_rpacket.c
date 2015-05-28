#include <assert.h>
#include <stdlib.h>
#include <Python.h>

#include "../src/rpacket.h"


static PyObject *
test_rpacket_get_nu(){
	rpacket_t * rp = (rpacket_t *) malloc(sizeof(rpacket_t));
	double value = 10.2;
	rpacket_set_nu(rp, value);
	if( value != rpacket_get_nu(rp) ){
		PyErr_Format(PyExc_AssertionError, "value %g not %g", value, rpacket_get_nu(rp));
	}
}

static PyMethodDef TestMethods[] = {
	{ "test_rpacket_get_nu", test_rpacket_get_nu, METH_VARARGS },
	{ NULL, NULL, 0, NULL }
};


/*
 * Python calls this to let us initialize our module
 */
void inittest_cmontecarlo()
{
  (void) Py_InitModule("test_cmontecarlo", TestMethods);
}
