/** @example GetEntities.cpp
 * Description: Get entities and report non-vertex entity connectivity and vertex adjacencies.\n
 * This example shows how to get connectivity and adjacencies.\n
 *
 * To run: ./GetEntities [meshfile]\n
 * (default values can run if users don't specify a mesh file)
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string(MESH_DIR) + string("/hex01.vtk");

int main(int argc, char **argv)
{
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }

  // Instantiate & load a mesh from a file
  Core* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);

  Range ents;

  // Get all entities in the database
  rval = mb->get_entities_by_handle(0, ents);MB_CHK_ERR(rval);

  for (Range::iterator it = ents.begin(); it != ents.end(); ++it) {
    if (MBVERTEX == mb->type_from_handle(*it)) {
      Range adjs;
      rval = mb->get_adjacencies(&(*it), 1, 3, false, adjs);MB_CHK_ERR(rval);
      cout << "Vertex " << mb->id_from_handle(*it) << " adjacencies:" << endl;
      adjs.print();
    }
    else if (mb->type_from_handle(*it) < MBENTITYSET) {
      const EntityHandle *connect;
      int num_connect;
      rval = mb->get_connectivity(*it, connect, num_connect);MB_CHK_ERR(rval);
      cout << CN::EntityTypeName(mb->type_from_handle(*it)) << " " << mb->id_from_handle(*it)
           << " vertex connectivity is: ";
      for (int i = 0; i < num_connect; i++)
        cout << mb->id_from_handle(connect[i]) << " ";
      cout << endl;
    }
  }

  delete mb;

  return 0;
}
