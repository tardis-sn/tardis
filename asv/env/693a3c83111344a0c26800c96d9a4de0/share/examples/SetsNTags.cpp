/** @example SetsNTags.cpp
 * Description: Get the sets representing materials and Dirichlet/Neumann boundary conditions and list their contents.\n
 * This example shows how to get entity sets, and tags on those sets.
 *
 * To run: ./SetsNTags [meshfile]\n
 * (default values can run if users don't specify a mesh file)
 */

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"

#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string(MESH_DIR) + string("/hex01.vtk");

// Tag names for these conventional tags come from MBTagConventions.hpp
const char *tag_nms[] = {MATERIAL_SET_TAG_NAME, DIRICHLET_SET_TAG_NAME, NEUMANN_SET_TAG_NAME};

int main(int argc, char **argv)
{
  // Get the material set tag handle
  Tag mtag;
  ErrorCode rval;
  Range sets, set_ents;

  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }

  // Load a file
  rval = mb->load_file(test_file_name.c_str());MB_CHK_ERR(rval);

  // Loop over set types
  for (int i = 0; i < 3; i++) {
    // Get the tag handle for this tag name; tag should already exist (it was created during file read)
    rval = mb->tag_get_handle(tag_nms[i], 1, MB_TYPE_INTEGER, mtag);MB_CHK_ERR(rval);

    // Get all the sets having that tag (with any value for that tag)
    sets.clear();
    rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, NULL, 1, sets);MB_CHK_ERR(rval);

    // Iterate over each set, getting the entities and printing them
    Range::iterator set_it;
    for (set_it = sets.begin(); set_it != sets.end(); ++set_it) {
      // Get the id for this set
      int set_id;
      rval = mb->tag_get_data(mtag, &(*set_it), 1, &set_id);MB_CHK_ERR(rval);

      // Get the entities in the set, recursively
      rval = mb->get_entities_by_handle(*set_it, set_ents, true);MB_CHK_ERR(rval);

      cout << tag_nms[i] << " " << set_id << " has " << set_ents.size() << " entities:" << endl;
      set_ents.print("   ");
      set_ents.clear();
    }
  }

  delete mb;

  return 0;
}
