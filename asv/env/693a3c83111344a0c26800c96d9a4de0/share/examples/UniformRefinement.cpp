/*This unit test is for the uniform refinement capability based on AHF datastructures*/
#include <iostream>
#include "moab/Core.hpp"
#include "moab/NestedRefine.hpp"

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

int main(int argc, char *argv[])
{

  Core mb;
  Interface* mbImpl = &mb;
  ErrorCode error;

  if (argc==1)
    {
      std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
      return 1;
    }

  const char *filename = argv[1];
  error = mbImpl->load_file(filename); MB_CHK_ERR(error);

  NestedRefine uref(&mb);

  // Usage: The level degrees array controls number of refinemetns and 
  // the degree of refinement at each level.
  // Example: int level_degrees[4] = {2,3,2,3};
  int level_degrees[2] = {3,2};
  int num_levels = sizeof(level_degrees) / sizeof(int);
  std::vector<EntityHandle> set;

  std::cout<<"Starting hierarchy generation"<<std::endl;

  error = uref.generate_mesh_hierarchy(num_levels, level_degrees, set); MB_CHK_ERR(error);

  std::cout<<"Finished hierarchy generation"<<std::endl;

  std::stringstream file;
  file <<"mesh_hierarchy.h5m";
  error = mbImpl->write_file(file.str().c_str()); MB_CHK_ERR(error);
  return 0;
}

