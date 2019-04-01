/** @example TestExodusII.cpp
 * This example demonstrates how to retrieve material, dirichlet and neumann sets
 * from an ExodusII file. \n
 * Sets in MOAB contain entities and have a tag on them to give meaning to the entities.
 * Tag names: MATERIAL_SET", "DIRICHLET_SET" and "NEUMANN_SET" are reserved and
 * are associated with their corresponding entity sets.
 * Sets are traversed to find out the type number of entities contained in each set. \n
 *
 * <b>Steps in this example </b>:
 *    -# Instantiate MOAB
 *    -# Get input mesh file name and load it.
 *    -# Loop over the three sets: material, dirichlet and neumann
 *      -# Get TagHandle and EntitySet(corresponding to the TagHandle)
 *      -# Loop thru all the EntitySet's
 *        -# Get the set id and entities in this set
 *    -# Destroy the MOAB instance
 *
 *
 * <b> To compile: </b>
 *    make TestExodusII MOAB_DIR=<installdir> \n
 *
 * <b> To run: </b>
 *    -# TestExodusII <mesh-file> \n
 *    -# TestExodusII (This uses the default <mesh-file>: <MOAB_SRC_DIR>/MeshFiles/unittest/mbtest2.g)
 */
#include <iostream>

// Include header for MOAB instance and range
#include "moab/Core.hpp"

using namespace moab;
using namespace std;

string test_file_name = string(MESH_DIR) + string("/mbtest2.g");
int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_NETCDF
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Get the material set tag handle
  Tag mtag;
  ErrorCode rval;
  const char *tag_nms[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  Range sets, set_ents;

  // Load a file
  if (argc == 1) {
    cout << "Running default case, loading " << test_file_name << endl;
    cout << "Usage: " << argv[0] << " <filename>\n" << endl;
    rval = mb->load_file(test_file_name.c_str());MB_CHK_ERR(rval);
  }
  else {
    rval = mb->load_file(argv[argc - 1]);MB_CHK_ERR(rval);
    cout << "Loaded mesh file: " << argv[argc - 1] << endl;
  }

  // Loop over set types
  for (int i = 0; i < 3; i++) {
    rval = mb->tag_get_handle(tag_nms[i], 1, MB_TYPE_INTEGER, mtag);MB_CHK_ERR(rval);

    // Get all the sets of that type in the mesh
    sets.clear();
    rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, NULL, 1, sets);MB_CHK_ERR(rval);

    // Iterate over each set, getting entities
    Range::iterator set_it;
    for (set_it = sets.begin(); set_it != sets.end(); ++set_it) {
      EntityHandle this_set = *set_it;

      // Get the id for this set
      int set_id;
      rval = mb->tag_get_data(mtag, &this_set, 1, &set_id);MB_CHK_ERR(rval);

      // Get the entities in the set, recursively
      rval = mb->get_entities_by_handle(this_set, set_ents, true);MB_CHK_ERR(rval);

      cout << tag_nms[i] << " " << set_id << " has " << set_ents.size() << " entities:" << endl;

      // Print the entities contained in this set
      set_ents.print("   ");
      set_ents.clear();
    }
  }

  delete mb;
#else
  cout<< " This test needs moab configured with netcdf \n"; 
#endif

  return 0;
}
