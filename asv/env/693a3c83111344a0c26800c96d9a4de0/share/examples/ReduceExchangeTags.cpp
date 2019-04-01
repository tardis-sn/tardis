/** @example ReduceExchangeTags.cpp
 * \brief Example program that shows the use case for performing tag data exchange
 * between parallel processors in order to sync data on shared entities. The reduction
 * operation on tag data is also shown where the user can perform any of the actions supported
 * by MPI_Op on data residing on shared entities. \n
 *
 * <b>This example </b>:
 *    -# Initialize MPI and instantiate MOAB
 *    -# Get user options: Input mesh file name, tag name (default: USERTAG), tag value (default: 1.0)
 *    -# Create the root and partition sets
 *    -# Instantiate ParallelComm and read the mesh file in parallel using appropriate options
 *    -# Create two tags: USERTAG_EXC (exchange) and USERTAG_RED (reduction)
 *    -# Set tag data and exchange shared entity information between processors
 *      -# Get entities in all dimensions and set local (current rank, dimension) dependent data for
 *     exchange tag (USERTAG_EXC)
 *      -# Perform exchange of tag data so that data on shared entities are synced via ParallelCommunicator.
 *    -#  Set tag data and reduce shared entity information between processors using MPI_SUM
 *      -#  Get higher dimensional entities in the current partition and set local (current rank)
 *     dependent data for reduce tag (USERTAG_EXC)
 *      -#  Perform the reduction operation (MPI_SUM) on shared entities via ParallelCommunicator.
 *    -#  Destroy the MOAB instance and finalize MPI
 *
 * <b>To run:</b> \n mpiexec -n 2 ./ReduceExchangeTags <mesh_file> <tag_name> <tag_value> \n
 * <b>Example:</b> \n mpiexec -n 2 ./ReduceExchangeTags ../MeshFiles/unittest/64bricks_1khex.h5m USERTAG 100 \n
 *
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif
#include "MBParallelConventions.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace moab;
using namespace std;

// Error routines for use with MPI API
#define MPICHKERR(CODE, MSG) \
  do { \
    if (0 != CODE) { \
      cerr << MSG << endl; \
      MPI_Finalize(); \
    } \
  } while (false)

#define dbgprint(MSG) \
  do { \
      if (!rank) cerr << MSG << endl; \
  } while (false)

#define dbgprintall(MSG) \
  do { \
      cerr << "[" << rank << "]: " << MSG << endl; \
  } while (false)

// Function to parse input parameters
ErrorCode get_file_options(int argc, char **argv,
                           string& filename,
                           string& tagName,
                           double& tagValues)
{
  // Get mesh filename
  if (argc > 1)
    filename = string(argv[1]);
  else
    filename = string(MESH_DIR) + string("/64bricks_1khex.h5m");

  // Get tag selection options
  if (argc > 2)
    tagName = string(argv[2]);
  else
    tagName = "USERTAG";

  if (argc > 3)
    tagValues = atof(argv[3]);
  else
    tagValues = 1.0;

  return MB_SUCCESS;
}

//
// Start of main test program
//
int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_MPI
  ErrorCode err;
  int ierr, rank;
  string filename, tagName;
  double tagValue;
  MPI_Comm comm = MPI_COMM_WORLD;
  /// Parallel Read options:
  ///   PARALLEL = type {READ_PART}
  ///   PARTITION = PARALLEL_PARTITION : Partition as you read
  ///   PARALLEL_RESOLVE_SHARED_ENTS : Communicate to all processors to get the shared adjacencies consistently in parallel
  ///   PARALLEL_GHOSTS : a.b.c
  ///                   : a = 3 - highest dimension of entities
  ///                   : b = 0 -
  ///                   : c = 1 - number of layers
  ///   PARALLEL_COMM = index
  string read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_DISTRIBUTE;PARALLEL_GHOSTS=3.0.1;PARALLEL_COMM=0";

  // Print usage if not enough arguments
  if (argc < 1) {
    cerr << "Usage: ";
    cerr << argv[0] << " <file_name> <tag_name> <tag_value>" << endl;
    cerr << "file_name    : mesh file name" << endl;
    cerr << "tag_name     : name of tag to add to mesh" << endl;
    cerr << "tag_value    : a double valued string to set for highest-dimensional entities in the mesh for the named tag" << endl;

    ierr = MPI_Finalize();
    MPICHKERR(ierr, "MPI_Finalize failed; Aborting");

    return 1;
  }

  // Initialize MPI first
  ierr = MPI_Init(&argc, &argv);
  MPICHKERR(ierr, "MPI_Init failed");

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPICHKERR(ierr, "MPI_Comm_rank failed");

  dbgprint("********** reduce_exchange_tags **********\n");

  // Create the moab instance
  Interface* mbi = new (std::nothrow) Core;
  if (NULL == mbi)
    return 1;

  // Get the input options
  err = get_file_options(argc, argv, filename, tagName, tagValue);MB_CHK_SET_ERR(err, "get_file_options failed");

  // Print out the input parameters
  dbgprint(" Input Parameters - ");
  dbgprint("   Filenames: " << filename);
  dbgprint("   Tag: Name=" << tagName << " Value=" << tagValue << endl);

  // Create root sets for each mesh.  Then pass these
  // to the load_file functions to be populated.
  EntityHandle rootset, partnset;
  err = mbi->create_meshset(MESHSET_SET, rootset);MB_CHK_SET_ERR(err, "Creating root set failed");
  err = mbi->create_meshset(MESHSET_SET, partnset);MB_CHK_SET_ERR(err, "Creating partition set failed");

  // Create the parallel communicator object with the partition handle associated with MOAB
  ParallelComm *parallel_communicator = ParallelComm::get_pcomm( mbi, partnset, &comm );

  // Load the file from disk with given options
  err = mbi->load_file(filename.c_str(), &rootset, read_options.c_str());MB_CHK_SET_ERR(err, "MOAB::load_file failed");

  // Create two tag handles: Exchange and Reduction operations
  dbgprint("-Creating tag handle " << tagName << "...");
  Tag tagReduce, tagExchange;
  {
    stringstream sstr;
    // Create the exchange tag: default name = USERTAG_EXC
    sstr << tagName << "_EXC";
    err = mbi->tag_get_handle(sstr.str().c_str(), 1, MB_TYPE_INTEGER, tagExchange,
        MB_TAG_CREAT|MB_TAG_DENSE, &tagValue);MB_CHK_SET_ERR(err, "Retrieving tag handles failed");

    // Create the exchange tag: default name = USERTAG_RED
    sstr.str(""); sstr << tagName << "_RED";
    err = mbi->tag_get_handle(sstr.str().c_str(), 1, MB_TYPE_DOUBLE, tagReduce,
        MB_TAG_CREAT|MB_TAG_DENSE, &tagValue);MB_CHK_SET_ERR(err, "Retrieving tag handles failed");
  }

  // Perform exchange tag data
  dbgprint("-Exchanging tags between processors ");
  {
    Range partEnts, dimEnts;
    for (int dim = 0; dim <= 3; dim++) {
      // Get all entities of dimension = dim
      err = mbi->get_entities_by_dimension(rootset, dim, dimEnts, false);MB_CHK_ERR(err);

      vector<int> tagValues(dimEnts.size(), static_cast<int>(tagValue)*(rank+1)*(dim+1));
      // Set local tag data for exchange
      err = mbi->tag_set_data(tagExchange, dimEnts,
          &tagValues[0]);MB_CHK_SET_ERR(err, "Setting local tag data failed during exchange phase");
      // Merge entities into parent set
      partEnts.merge(dimEnts);
    }

    // Exchange tags between processors
    err = parallel_communicator->exchange_tags(tagExchange, partEnts);MB_CHK_SET_ERR(err, "Exchanging tags between processors failed");
  }

  // Perform reduction of tag data
  dbgprint("-Reducing tags between processors ");
  {
    Range partEnts;
    // Get all higher dimensional entities belonging to current partition
    err = parallel_communicator->get_part_entities(partEnts);MB_CHK_SET_ERR(err, "ParallelComm::get_part_entities failed");

    // Output what is in current partition sets
    dbgprintall("Number of Partitioned entities: " << partEnts.size());
    MPI_Barrier(comm);

    // Set local tag data for reduction
    vector<double> tagValues(partEnts.size(), tagValue*(rank+1));
    err = mbi->tag_set_data(tagReduce, partEnts,
        &tagValues[0]);MB_CHK_SET_ERR(err, "Setting local tag data failed during reduce phase");

    Range dummy;
    // Reduce tag data using MPI_SUM on the interface between partitions
    err = parallel_communicator->reduce_tags(tagReduce, MPI_SUM,
        dummy/*partEnts*/);MB_CHK_SET_ERR(err, "Reducing tags between processors failed");
  }
  // Write out to output file to visualize reduction/exchange of tag data
  err = mbi->write_file("test.h5m", "H5M", "PARALLEL=WRITE_PART");MB_CHK_ERR(err);

  // Done, cleanup
  delete mbi;

  dbgprint("\n********** reduce_exchange_tags DONE! **********");

  MPI_Finalize();
#else
  std::cout <<" compile with MPI and HDF5 for this example to work \n";
#endif
  return 0;
}
