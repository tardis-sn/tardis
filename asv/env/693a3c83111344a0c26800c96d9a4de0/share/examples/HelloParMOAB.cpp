/** @example HelloParMOAB.cpp \n
 * \brief Read mesh into MOAB and resolve/exchange/report shared and ghosted entities \n
 * <b>To run</b>: mpiexec -np 4 HelloParMOAB [filename]\n
 *
 *  It shows how to load the mesh independently, on multiple
 *  communicators (with second argument, the number of comms)
 *
 *  mpiexec -np 8 HelloParMOAB [filename] [nbComms]
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif
#include "MBParallelConventions.h"
#include <iostream>

using namespace moab;
using namespace std;

string test_file_name = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");

int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);

  string options;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }  

  int nbComms = 1;
  if (argc > 2)
    nbComms = atoi(argv[2]);

  options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";

  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  MPI_Comm comm;
  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_size);

  int color = global_rank % nbComms; // For each angle group a different color
  if (nbComms > 1) {
    // Split the communicator, into ngroups = nbComms
    MPI_Comm_split(MPI_COMM_WORLD, color, global_rank, &comm);
  }
  else
    comm = MPI_COMM_WORLD;

  // Get the ParallelComm instance
  ParallelComm* pcomm = new ParallelComm(mb, comm);
  int nprocs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();
#ifndef NDEBUG
  MPI_Comm rcomm = pcomm->proc_config().proc_comm();
  assert(rcomm == comm);
#endif
  if (0 == global_rank)
    cout << " global rank:" << global_rank << " color:" << color << " rank:" << rank << " of " << nprocs << " processors\n";

  if (1 == global_rank)
    cout << " global rank:" << global_rank << " color:" << color << " rank:" << rank << " of " << nprocs << " processors\n";

  MPI_Barrier(MPI_COMM_WORLD);

  if (0 == global_rank)
    cout << "Reading file " << test_file_name << "\n with options: " << options <<
         "\n on " << nprocs << " processors on " << nbComms << " communicator(s)\n";

  // Read the file with the specified options
  ErrorCode rval = mb->load_file(test_file_name.c_str(), 0, options.c_str());MB_CHK_ERR(rval);

  Range shared_ents;
  // Get entities shared with all other processors
  rval = pcomm->get_shared_entities(-1, shared_ents);MB_CHK_ERR(rval);

  // Filter shared entities with not not_owned, which means owned
  Range owned_entities;
  rval = pcomm->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_entities);MB_CHK_ERR(rval);

  unsigned int nums[4] = {0}; // to store the owned entities per dimension
  for (int i = 0; i < 4; i++)
    nums[i] = (int)owned_entities.num_of_dimension(i);
  vector<int> rbuf(nprocs*4, 0);
  MPI_Gather(nums, 4, MPI_INT, &rbuf[0], 4, MPI_INT, 0, comm);
  // Print the stats gathered:
  if (0 == global_rank) {
    for (int i = 0; i < nprocs; i++)
      cout << " Shared, owned entities on proc " << i << ": " << rbuf[4*i] << " verts, " <<
          rbuf[4*i + 1] << " edges, " << rbuf[4*i + 2] << " faces, " << rbuf[4*i + 3] << " elements" << endl;
  }

  // Now exchange 1 layer of ghost elements, using vertices as bridge
  // (we could have done this as part of reading process, using the PARALLEL_GHOSTS read option)
  rval = pcomm->exchange_ghost_cells(3, // int ghost_dim
                                     0, // int bridge_dim
                                     1, // int num_layers
                                     0, // int addl_ents
                                     true);MB_CHK_ERR(rval); // bool store_remote_handles

  // Repeat the reports, after ghost exchange
  shared_ents.clear();
  owned_entities.clear();
  rval = pcomm->get_shared_entities(-1, shared_ents);MB_CHK_ERR(rval);
  rval = pcomm->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_entities);MB_CHK_ERR(rval);

  // Find out how many shared entities of each dimension are owned on this processor
  for (int i = 0; i < 4; i++)
    nums[i] = (int)owned_entities.num_of_dimension(i);

  // Gather the statistics on processor 0
  MPI_Gather(nums, 4, MPI_INT, &rbuf[0], 4, MPI_INT, 0, comm);
  if (0 == global_rank) {
    cout << " \n\n After exchanging one ghost layer: \n";
    for (int i = 0; i < nprocs; i++) {
      cout << " Shared, owned entities on proc " << i << ": " << rbuf[4*i] << " verts, " <<
          rbuf[4*i + 1] << " edges, " << rbuf[4*i + 2] << " faces, " << rbuf[4*i + 3] << " elements" << endl;
    }
  }

  delete mb;

  MPI_Finalize();
#else
  std::cout<<" compile with MPI and hdf5 for this example to work\n";

#endif
  return 0;
}
