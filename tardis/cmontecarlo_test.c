#include "cmontecarlo.h"
#include "CuTest.h"

/*
  Running tests by hand on Linux:
  make
  ./cmontecarlo_test
*/

void example_test_success(CuTest *tc)
{
  CuAssertTrue(tc, 0 == 0);
}

void example_test_failure(CuTest *tc)
{
  CuAssertTrue(tc, 0 == 1);
}

int main(void)
{
  CuString *output = CuStringNew();
  CuSuite *suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, example_test_success);
  SUITE_ADD_TEST(suite, example_test_failure);
  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("%s\n", output->buffer);
}
