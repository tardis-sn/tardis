#include "cmontecarlo.h"
#include "test_macros.h"
#include "stdlib.h"


void test_rpacket_get_set()
{
	double nu = 10.3;
	rpacket_t * packet = (rpacket_t * ) malloc(sizeof(rpacket_t));
	rpacket_set_nu(packet, nu);
	double ret = rpacket_get_nu( packet );
	ASSERT_EQUAL_FLOAT(ret, nu);
}

int main()
{

	test_rpacket_get_set();
	return 0;
}
