/** @example HelloMOAB.cpp
 * Description: read a mesh, get the entities.\n
 * HelloMOAB is a simple test file which is used to read meshes from VTK file and test how many entities there are.\n
 *
 * To run: ./HelloMOAB [meshfile]\n
 * (default values can run if users don't specify a mesh file)
 */


#include "moab/Core.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string(MESH_DIR) + string("/3k-tri-sphere.vtk");

int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }

  // Load the mesh from vtk file
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);

  // Get verts entities, by type
  Range verts;
  rval = mb->get_entities_by_type(0, MBVERTEX, verts);MB_CHK_ERR(rval);

  // Get edge entities, by type
  Range edges;
  rval = mb->get_entities_by_type(0, MBEDGE, edges);MB_CHK_ERR(rval);

  // Get faces, by dimension, so we stay generic to entity type
  Range faces;
  rval = mb->get_entities_by_dimension(0, 2, faces);MB_CHK_ERR(rval);

  // Get regions, by dimension, so we stay generic to entity type
  Range elems;
  rval = mb->get_entities_by_dimension(0, 3, elems);MB_CHK_ERR(rval);

  // Output the number of entities
  cout << "Number of vertices is " << verts.size() << endl;
  cout << "Number of edges is " << edges.size() << endl;
  cout << "Number of faces is " << faces.size() << endl;
  cout << "Number of elements is " << elems.size() << endl;

  delete mb;

  return 0;
}
