/** @example TestErrorHandlingPar.cpp \n
 * Description: This example tests MOAB's trace back error handler in parallel.\n
 *
 * <b>To run</b>: mpiexec -np <n> ./TestErrorHandlingPar <test_case_num(1 to 2)> \n
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <iostream>

using namespace moab;
using namespace std;

// In this test case, a global fatal error MB_NOT_IMPLEMENTED is returned by MOAB
// Note, it is printed by root processor 0 only
ErrorCode TestErrorHandlingPar_1()
{
  Core moab;
  Interface& mb = moab;

  string opts = ";;";
#ifdef MOAB_HAVE_MPI
  // Use parallel options
  opts += "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ";
#endif

  // Load a CAM-FV file and read a variable on edges (not supported yet)
  string test_file = string(MESH_DIR) + string("/io/fv3x46x72.t.3.nc");
  opts += ";VARIABLE=US";
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, opts.c_str());MB_CHK_ERR(rval);

  return MB_SUCCESS;
}

// In this test case, a per-processor relevant error MB_FAILURE is returned by MOAB
// Note, it is printed by all processors
ErrorCode TestErrorHandlingPar_2()
{
  Core moab;
  Interface& mb = moab;

  string opts = ";;";
#ifdef MOAB_HAVE_MPI
  // Use parallel options
  opts += "PARALLEL=READ_PART;PARTITION_METHOD=UNKNOWN";
#endif

  // Load a CAM-FV file with an unknown partition method specified
  string test_file = string(MESH_DIR) + string("/io/fv3x46x72.t.3.nc");
  opts += ";VARIABLE=T";
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, opts.c_str());MB_CHK_ERR(rval);

  return MB_SUCCESS;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <test_case_num(1 to 2>" << endl;
    return 0;
  }

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Initialize error handler, optional for this example (using moab instances)
  MBErrorHandler_Init();

  ErrorCode rval = MB_SUCCESS;

  int test_case_num = atoi(argv[1]);
  switch (test_case_num) {
    case 1:
      rval = TestErrorHandlingPar_1();MB_CHK_ERR(rval);
      break;
    case 2:
      rval = TestErrorHandlingPar_2();MB_CHK_ERR(rval);
      break;
    default:
      break;
  }

  // Finalize error handler, optional for this example (using moab instances)
  MBErrorHandler_Finalize();

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
