#include "test_macros.h"
#include "../montecarlo/src/cmontecarlo.h"
#include <stdio.h>

void test_c_value()
{
    ASSERT_EQUAL_FLOAT(29979245800.0, C);
}
int main()
{
    test_c_value();
    return 0;
}
