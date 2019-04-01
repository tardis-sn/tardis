/** @example TestErrorHandling.cpp \n
 * Description: This example tests MOAB's trace back error handler in serial. \n
 *
 * <b>To run</b>: ./TestErrorHandling <test_case_num(1 to 4)> \n
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <iostream>

using namespace moab;
using namespace std;

// In this test case, an error MB_NOT_IMPLEMENTED is returned by MOAB
ErrorCode TestErrorHandling_1()
{
  Core moab;
  Interface& mb = moab;

  // Load a CAM-FV file and read a variable on edges (not supported yet)
  string test_file = string(MESH_DIR) + string("/io/fv3x46x72.t.3.nc");
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, "VARIABLE=US");MB_CHK_ERR(rval);

  return MB_SUCCESS;
}

// In this test case, an error MB_TYPE_OUT_OF_RANGE is returned by MOAB
ErrorCode TestErrorHandling_2()
{
  Core moab;
  Interface& mb = moab;

  // Load a HOMME file with an invalid GATHER_SET option
  string test_file = string(MESH_DIR) + string("/io/homme3x3458.t.3.nc");
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, "VARIABLE=T;GATHER_SET=0.1");MB_CHK_ERR(rval);

  return MB_SUCCESS;
}

// In this test case, an error MB_FAILURE is returned by MOAB
ErrorCode TestErrorHandling_3()
{
  Core moab;
  Interface& mb = moab;

  // Load a CAM-FV file with NOMESH option and a NULL file set
  string test_file = string(MESH_DIR) + string("/io/fv3x46x72.t.3.nc");
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, "NOMESH;VARIABLE=");MB_CHK_ERR(rval);

  return MB_SUCCESS;
}

// In this test case, an error MB_VARIABLE_DATA_LENGTH is returned by MOAB
ErrorCode TestErrorHandling_4()
{
  Core moab;
  Interface& mb = moab;

  // Create 100 vertices
  const int NUM_VTX = 100;
  vector<double> coords(3 * NUM_VTX);
  Range verts;
  ErrorCode rval = mb.create_vertices(&coords[0], NUM_VTX, verts);MB_CHK_SET_ERR(rval, "Failed to create vertices");

  // Create a variable-length dense tag
  Tag tag;
  rval = mb.tag_get_handle("var_len_den", 1, MB_TYPE_INTEGER, tag,
                          MB_TAG_VARLEN | MB_TAG_DENSE | MB_TAG_CREAT);MB_CHK_SET_ERR(rval, "Failed to create a tag");

  // Attempt to iterate over a variable-length tag, which will never be possible
  void* ptr = NULL;
  int count = 0;
  rval = mb.tag_iterate(tag, verts.begin(), verts.end(),
                        count, ptr);MB_CHK_SET_ERR(rval, "Failed to iterate over tag on " << NUM_VTX << " vertices");

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

  // Initialize error handler, optional for this example (using moab instances)
  MBErrorHandler_Init();

  ErrorCode rval = MB_SUCCESS;

  int test_case_num = atoi(argv[1]);
  switch (test_case_num) {
    case 1:
      rval = TestErrorHandling_1();MB_CHK_ERR(rval);
      break;
    case 2:
      rval = TestErrorHandling_2();MB_CHK_ERR(rval);
      break;
    case 3:
      rval = TestErrorHandling_3();MB_CHK_ERR(rval);
      break;
    case 4:
      rval = TestErrorHandling_4();MB_CHK_ERR(rval);
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
