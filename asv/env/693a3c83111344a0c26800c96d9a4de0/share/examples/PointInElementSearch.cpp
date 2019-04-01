/** \brief This test shows how to perform local point-in-element searches with MOAB's new tree searching functionality.  
 *
 * MOAB's SpatialLocator functionality performs point-in-element searches over a local or parallel mesh.
 * SpatialLocator is flexible as to what kind of tree is used and what kind of element basis functions are 
 * used to localize elements and interpolate local fields.
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/LinearHex.hpp"
#include "moab/CN.hpp"
#include "moab/SpatialLocator.hpp"

using namespace moab;
using namespace std;

#ifdef MOAB_HAVE_HDF5
string test_file_name = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");
#else
string test_file_name = string(MESH_DIR) + string("/mbtest1.vtk");
#endif
int main(int argc, char **argv)
{
  int num_queries = 1000000;

  if (argc > 3) {
    cout << "Usage: " << argv[0] << " <filename> [num_queries]" << endl;
    return 0;
  }
  else if (argc == 3) {
    test_file_name = argv[1];
    num_queries = atoi(argv[2]);
  }
  else {
    num_queries = 100;
  }

  // Instantiate
  Core mb;

  // Load the file
  ErrorCode rval = mb.load_file(test_file_name.c_str());MB_CHK_SET_ERR(rval, "Error loading file");

  // Get all 3d elements in the file
  Range elems;
  rval = mb.get_entities_by_dimension(0, 3, elems);MB_CHK_SET_ERR(rval, "Error getting 3d elements");

  // Create a tree to use for the location service
  AdaptiveKDTree tree(&mb);

  // Specify an evaluator based on linear hexes
  ElemEvaluator el_eval(&mb);

  // Build the SpatialLocator
  SpatialLocator sl(&mb, elems, &tree);
  
  // Get the box extents
  CartVect box_extents, pos;
  BoundBox box = sl.local_box();
  box_extents = box.bMax - box.bMin;

  // Query at random places in the tree
  CartVect params;
  int is_inside = 0;
  int num_inside = 0;
  EntityHandle elem;
  for (int i = 0; i < num_queries; i++) {
    pos = box.bMin + CartVect(box_extents[0] * .01 * (rand() % 100), box_extents[1] * .01 * (rand() % 100),
        box_extents[2] * .01 * (rand() % 100));
    rval = sl.locate_point(pos.array(), elem, params.array(), &is_inside, 0.0, 0.0);MB_CHK_ERR(rval);
    if (is_inside) num_inside++;
  }

  cout << "Mesh contains " << elems.size() << " elements of type "
            << CN::EntityTypeName(mb.type_from_handle(*elems.begin())) << endl;
  cout << "Bounding box min-max = (" << box.bMin[0] << "," << box.bMin[1] << "," << box.bMin[2] << ")-("
            << box.bMax[0] << "," << box.bMax[1] << "," << box.bMax[2] << ")" << endl;
  cout << "Queries inside box = " << num_inside << "/" << num_queries << " = "
            << 100.0*((double)num_inside) / num_queries << "%" << endl;

  return 0;
}
