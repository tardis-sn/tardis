/** @example LoadPartial.cpp \n
 * \brief Load a part of a file  \n
 * <b>To run</b>: LoadPartial <file> <tag_name> <val1> <val2> ...\n
 *
 * In this example, it is shown how to load only a part of one file; the file must be organized in sets.
 * (cherry-picking only the sets we want)
 * The sets to load are identified by a tag name and the tag values for the sets of interest.
 * This procedure is used  when reading in parallel, as each processor will load only
 * its part of the file, identified either by partition or by material/block sets
 *  by default, this example will load parallel partition sets
 *  with values 1, 2, and 5 from ../MeshFiles/unittest/64bricks_1khex.h5m
 *  The example will always write the output to a file name part.h5m
 */

#include <iostream>
#include <vector>

// Include header for MOAB instance and tag conventions for
#include "moab/Core.hpp" 
#include "MBTagConventions.hpp"

using namespace moab;
using namespace std;

string test_file_name = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");
int main(int argc, char **argv)
{
  // Get MOAB instance
#ifdef MOAB_HAVE_HDF5
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  ErrorCode rval;
  if (argc <= 1) {
    // The default file to load
    int set_tag_values[] = {1, 2, 5};
    int num_set_tag_values = 3;
    // This file is in the mesh files directory
    rval = mb->load_file(test_file_name.c_str(),
            0, 0, PARALLEL_PARTITION_TAG_NAME, set_tag_values, num_set_tag_values);MB_CHK_SET_ERR(rval, "Failed to read");
  }
  else {
    // First arg is input file, second is tag name, then are the tag values
    if (argc < 4) {
      cout << " usage is " << argv[0] << " <file> <tag_name> <value1> <value2> .. \n";
      return 0;
    }
    else {
      vector<int> vals(argc - 3); // The first 3 args are exe, file, tag name; the rest are values
      for (int i = 3; i < argc; i++)
        vals[i - 3] = atoi(argv[i]);
      rval = mb->load_file(argv[1], 0, 0, argv[2], &vals[0], (int)vals.size());MB_CHK_SET_ERR(rval, "Failed to read");
    }
  }

  // If HANDLEID tag present, convert to long, and see what we read from file
  Tag handleid_tag;
  rval = mb->tag_get_handle("HANDLEID", handleid_tag);
  if (MB_SUCCESS == rval) {
    // Convert a few values for a few vertices
    Range verts;
    rval = mb->get_entities_by_type(0, MBVERTEX, verts);MB_CHK_SET_ERR(rval, "Failed to get vertices");
    vector<long> valsTag(verts.size());
    rval = mb->tag_get_data(handleid_tag, verts, &valsTag[0]);
    if (MB_SUCCESS == rval)
      cout << "First 2 long values recovered: " << valsTag[0] << " " << valsTag[1] << "\n";
  }

  rval = mb->write_file("part.h5m");MB_CHK_SET_ERR(rval, "Failed to write partial file");
  cout << "Wrote successfully part.h5m.\n";

  delete mb;
#else
  std::cout <<"configure MOAB with hdf5 for this to work\n";
#endif

  return 0;
}
