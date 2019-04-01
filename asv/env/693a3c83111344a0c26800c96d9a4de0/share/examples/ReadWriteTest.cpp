/** @example ReadWriteTest.cpp \n
 * \brief Read mesh into MOAB and write some back \n
 *
 * <b>To run</b>: mpiexec -np 4 ReadWriteTest [input] [output] -O <read_opts> -o <write_opts>\n
 *
 * used for stress test of reader/writer
 *  report times to read and write
 *
 *  example ReadWriteTest ../MeshFiles/io/fv3x46x72.t.3.nc out.nc  \
 *  -O PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=T,U;  \
 *  -o PARALLEL=WRITE_PART;VARIABLE=T,U
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif
#include <iostream>
#include <time.h>

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_MPI
  string input_file,output_file,read_opts,write_opts;

  MPI_Init(&argc, &argv);

  // Need option handling here for input filename
  if (argc < 3) {
#ifdef MOAB_HAVE_NETCDF
    input_file = string(MESH_DIR) + string("/io/fv3x46x72.t.3.nc");
    output_file = "ReadWriteTestOut.h5m";
    read_opts = "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=T,U";
    write_opts = "PARALLEL=WRITE_PART";
#else
    cout << "Usage: mpiexec -n $NP ReadWriteTest [input] [output] -O <read_opts> -o <write_opts>\n";
    return 0;
#endif
  }
  else {
    input_file = argv[1];
    output_file = argv[2];
  }

  if (argc > 3) {
    int index = 3;
    while (index < argc) {
      if (!strcmp(argv[index], "-O")) // This is for reading options, optional
        read_opts = argv[++index];
      if (!strcmp(argv[index], "-o"))
        write_opts = argv[++index];
      index++;
    }
  }

  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Get the ParallelComm instance
  ParallelComm* pcomm = new ParallelComm(mb, MPI_COMM_WORLD);
  int nprocs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  EntityHandle set;
  ErrorCode rval = mb->create_meshset(MESHSET_SET, set);MB_CHK_ERR(rval);

  clock_t tt = clock();

  if (0 == rank)
    cout << "Reading file " << input_file << "\n with options: " << read_opts
         << "\n on " << nprocs << " processors\n";

  // Read the file with the specified options
  rval = mb->load_file(input_file.c_str(), &set, read_opts.c_str());MB_CHK_ERR(rval);

  if (0 == rank) {
    cout << "Time: " << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << endl;
    tt = clock();
  }

  // Write back the file with the specified options
  rval = mb->write_file(output_file.c_str(), 0, write_opts.c_str(), &set, 1);MB_CHK_ERR(rval);

  if (0 == rank) {
    cout << "Writing file " << output_file << "\n with options: " << write_opts << endl;
    cout << "Time:  " << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << endl;
    tt = clock();
  }

  delete mb;

  MPI_Finalize();
#else
  std::cout << " compile MOAB with mpi for this example to work\n";
#endif

  return 0;
}
