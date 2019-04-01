/** @example ErrorHandlingSimulation.cpp
 * Description: This example simulates MOAB's enhanced error handling in parallel. \n
 * All of the errors are contrived, used for simulation purpose only. \n
 *
 * Note: We do not need a moab instance for this example
 *
 * <b>To run</b>: mpiexec -np 4 ./ErrorHandlingSimulation <test_case_num(1 to 4)> \n
 */

#include "moab/MOABConfig.h"
#include "moab/ErrorHandler.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <iostream>
#include <stdlib.h>

using namespace moab;
using namespace std;

// Functions that create and handle contrived errors
// Call hierarchy: A calls B, and B calls C
ErrorCode FunctionC(int test_case_num, int rank)
{
  switch (test_case_num) {
    case 1:
      // Simulate a global fatal error MB_NOT_IMPLEMENTED on all processors
      // Note, it is printed by root processor 0 only
      MB_SET_GLB_ERR(MB_NOT_IMPLEMENTED, "A contrived global error MB_NOT_IMPLEMENTED");
      break;
    case 2:
      // Simulate a per-processor relevant error MB_INDEX_OUT_OF_RANGE on all processors
      // Note, it is printed by all processors
      MB_SET_ERR(MB_INDEX_OUT_OF_RANGE, "A contrived error MB_INDEX_OUT_OF_RANGE on processor " << rank);
      break;
    case 3:
      // Simulate a per-processor relevant error MB_TYPE_OUT_OF_RANGE on all processors except root
      // Note, it is printed by all non-root processors
      if (0 != rank)
        MB_SET_ERR(MB_TYPE_OUT_OF_RANGE, "A contrived error MB_TYPE_OUT_OF_RANGE on processor " << rank);
      break;
    case 4:
      // Simulate a per-processor relevant error MB_INDEX_OUT_OF_RANGE on processor 1
      // Note, it is printed by processor 1 only
      if (1 == rank)
        MB_SET_ERR(MB_INDEX_OUT_OF_RANGE, "A contrived error MB_INDEX_OUT_OF_RANGE on processor 1");

      // Simulate a per-processor relevant error MB_TYPE_OUT_OF_RANGE on processor 3
      // Note, it is printed by processor 3 only
      if (3 == rank)
        MB_SET_ERR(MB_TYPE_OUT_OF_RANGE, "A contrived error MB_TYPE_OUT_OF_RANGE on processor 3");
      break;
    default:
      break;
  }

  return MB_SUCCESS;
}

ErrorCode FunctionB(int test_case_num, int rank)
{
  ErrorCode err_code = FunctionC(test_case_num, rank);MB_CHK_ERR(err_code);

  return MB_SUCCESS;
}

ErrorCode FunctionA(int test_case_num, int rank)
{
  ErrorCode err_code = FunctionB(test_case_num, rank);MB_CHK_ERR(err_code);

  return MB_SUCCESS;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <test_case_num(1 to 4)>" << endl;
    return 0;
  }

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Initialize error handler, required for this example (not using a moab instance)
  MBErrorHandler_Init();

  int test_case_num = atoi(argv[1]);
  int rank = 0;
#ifdef MOAB_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  ErrorCode rval = FunctionA(test_case_num, rank);MB_CHK_ERR(rval);

  // Finalize error handler, required for this example (not using a moab instance)
  MBErrorHandler_Finalize();

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
