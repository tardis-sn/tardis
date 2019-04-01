/** @example VisTags.cpp \n
 * \brief tool for visualizing multi level tags  \n
 * <b>To run</b>: VisTags  <inp_file>  <outfile> -O <read_opts> -t <tags> -l <levels>  -d <dim> \n
 *
 * In this example, it is shown how to create some simple tags for those tags that come from 
 *  climate data, multiple levels.
 *  you can read directly nc data, or *.h5m file that will have the tag with multi levels
 *   output will be a vtk file with dense tags of form tag_name_<level> 
 * the tag name might contain a time index too, like T0 or U0
 * <tag> is a list of tags, separated by commas, no spaces
 * <levels> is a list of levels, separated by commas, no spaces
 *  dimension of entities with the tags will be specified with -d (default 2)
 *
 * an example of use
 *
 * VisTags gcrm_r3.nc out.vtk -O VARIABLE=u -t u0,u1 -l 0,1,2 -d 2
 * (we knew that it had variable u in the file, that it had 256 levels, that there are 2 time
 *  steps, etc)
 *
 * or
 *  VisTags gcrm_r3.nc out.vtk  -t u0 -l 0,1,2 -d 2
 *  (it will read all variables, but we need to know that u0 will be created as a tag)
 *
 *  the out.vtk file will contain u0_0, u0_1, as simple dense double tags
 */

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

// Include header for MOAB instance and tag conventions
#include "moab/Core.hpp" 
#include "MBTagConventions.hpp"
#include "moab/FileOptions.hpp"

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_NETCDF
  ErrorCode rval;
  string file_input,file_output;
  string read_opts, tags; // Tags to write, separated by commas; it is the name of the tag
  if (argc < 2) {
    file_input = string(MESH_DIR) + string("/io/gcrm_r3.nc");
    file_output = "VisTagsOut.vtk";
  }
  else {
    file_input = argv[1];
    file_output = argv[2];
  }
  read_opts = "";
  tags = "";

  // Instantiate
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  int dimension = 2;
  // In MOAB, it may have index after reading (T0, T1, etc)
  char* levels = NULL; // Levels, separated by commas, no spaces (like 0, 1, 19)
  if (argc > 3) {
    int index = 3;
    while (index < argc) {
      if (!strcmp(argv[index], "-O")) // This is for reading options, optional
        read_opts = argv[++index];
      if (!strcmp(argv[index], "-t"))
        tags = argv[++index];
      if (!strcmp(argv[index], "-l"))
        levels = argv[++index];
      if (!strcmp(argv[index], "-d"))
        dimension = atoi(argv[++index]);
      index++;
    }
  }

  ostringstream opts;
  opts << ";;TAGS=" << tags << ";LEVELS=" << levels << "\0" ;
  FileOptions fo(opts.str().c_str());

  vector<string> tagsNames;
  vector<int> levelsArray;
  fo.get_strs_option("TAGS", tagsNames);
  fo.get_ints_option("LEVELS", levelsArray);

  // Load the input file with the specified options
  rval = mb->load_file(file_input.c_str(), 0, read_opts.c_str());MB_CHK_SET_ERR(rval, "not loading file");

  Range ents;
  rval = mb->get_entities_by_dimension(0, dimension, ents);MB_CHK_SET_ERR(rval, "not getting ents");

  // Now create double tags for entities of dimension
  for (size_t i = 0; i < tagsNames.size(); i++) {
    string tagName = tagsNames[i];
    Tag tagh;
    rval = mb->tag_get_handle(tagName.c_str(), tagh);
    if (MB_SUCCESS != rval) {
      cout << "not getting tag " << tagName.c_str() << "\n";
      continue;
    }

    int len = 0;
    rval = mb->tag_get_length(tagh, len);
    if (MB_SUCCESS != rval) {
      cout << "not getting tag len " << tagName.c_str() << "\n";
      continue;
    }

    DataType type;
    rval = mb->tag_get_data_type(tagh, type) ;
    if (MB_SUCCESS != rval) {
      cout << "not getting tag type " << tagName.c_str() << "\n";
      continue;
    }

    int count;
    void* dataptr; // Assume double tags, for simplicity
    rval = mb->tag_iterate(tagh,
        ents.begin(),
        ents.end(),
        count,
        dataptr);
    if (MB_SUCCESS != rval || count != (int)ents.size()) {
      cout << "not getting tag iterate right " << tagName.c_str() << "\n";
      continue;
    }

    // Now create a new tag, with a new name, concatenated, and copy data there , for each level
    for (size_t j = 0; j < levelsArray.size(); j++) {
      int level = levelsArray[j];
      if (level >= len) {
        cout << "level too big at " << level << "\n";
        continue;
      }

      ostringstream newTagName;
      newTagName << tagName << "_" << level  ;
      Tag newTagh;
      rval = mb->tag_get_handle(newTagName.str().c_str(), 1, type, newTagh,
          MB_TAG_DENSE | MB_TAG_CREAT);
      if (MB_SUCCESS != rval) {
        cout <<"not getting new tag " << newTagName.str() << "\n";
        continue;
      }

      void* newDataPtr;
      rval = mb->tag_iterate(newTagh,
                            ents.begin(),
                            ents.end(),
                            count,
                            newDataPtr);
      if (MB_SUCCESS != rval || count != (int) ents.size()) {
        cout << "not getting new tag iterate " << newTagName.str() << "\n";
        continue;
      }

      if (MB_TYPE_DOUBLE == type) {
        double* ptrD = (double*)newDataPtr;
        double* oldData = (double*)dataptr;
        for (int k = 0; k < count; k++, ptrD++)
          *ptrD = oldData[level + count*k];
      }
    } // for (size_t j = 0; j < levelsArray.size(); j++)

    mb->tag_delete(tagh); // No need for the tag anymore, write it to the new file
  } // for (size_t i = 0; i < tagsNames.size(); i++)

  rval = mb->write_file(file_output.c_str());MB_CHK_SET_ERR(rval, "Can't write file " << file_output);
  cout << "Successfully wrote file " << file_output << "\n";

  delete mb;
#else
  std::cout <<" configure with netcdf for this example to work\n";
#endif
  return 0;
}
