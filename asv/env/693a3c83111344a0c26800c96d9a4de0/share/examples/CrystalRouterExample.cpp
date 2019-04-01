/*
 * This example will show one of the building blocks of parallel infrastructure in MOAB
 * More exactly, if we have some homogeneous data to communicate from each processor to a list of other
 * processors, how do we do it?
 *
 * introduce the TupleList and crystal router to MOAB users.
 *
 * This technology is used in resolving shared vertices / sets between partitions
 * It is used in the mbcoupler for sending data (target points) to the proper processor, and communicate
 *   back the results.
 * Also, it is used to communicate departure mesh for intersection in parallel
 *
 *  It is a way of doing  MPI_gatheralltoallv(), when the communication matrix is sparse
 *
 *  It is assumed that every proc needs to communicate only with a few of the other processors.
 *  If every processor needs to communicate with all other, then we will have to use paired isend and irecv, the
 *  communication matrix is full
 *
 *  the example needs to be launched in parallel.
 *  Every proc will build a list of tuples, that will be send to a few procs;
 *  In general, we will send to num_comms tasks, and about num_tuples to each task
 *  We vary num_comms and num_tuples for processor
 *
 *  we will send long ints of the form
 *    100000 * send + 1000* rank +j, where j is the index of tuple
 *
 *  after routing, we verify we received
 *    100000 * rank + 1000 * from
 *
 *    For some reportrank we also print the tuples.
 *
 *  after routing, we will see if we received, as expected. Should run on at least 2 processors.
 *
 * Note: We do not need a moab instance for this example
 *
 */

/** @example CrystalRouterExample.cpp \n
 * \brief generalized gather scatter using tuples \n
 * <b>To run</b>: mpiexec -np <n> CrystalRouterExample -r [reportrank] -t [num_tuples] -n [num_comms] \n
 *
 */
//
#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
#include "moab/ProcConfig.hpp"
#endif
#include "moab/TupleList.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ErrorHandler.hpp"
#include <time.h>
#include <iostream>
#include <sstream>

const char BRIEF_DESC[] = "Example of gather scatter with tuple lists \n";
std::ostringstream LONG_DESC;

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);

  // Initialize error handler, required for this example (not using a moab instance)
  MBErrorHandler_Init();

  ProcConfig pc(MPI_COMM_WORLD);
  int size = pc.proc_size();
  int rank = pc.proc_rank();

  // Start copy
  LONG_DESC << "This program does a gather scatter with a list of tuples. \n"
          " It tries to see how much communication costs in terms of time and memory. \n"
          << "It starts with creating a list of tuples to be sent from each processor, \n to a list of other processors.\n" <<
          "The number of tuples and how many tasks to communicate to are controlled by input parameters.\n" <<
          "After communication, we verify locally if we received what we expected. \n";
  ProgOptions opts(LONG_DESC.str(), BRIEF_DESC);

  // How many procs communicate to current proc, on average (we will vary that too)
  int num_comms = 2;
  opts.addOpt<int>("num_comms,n",
       "each task will send to about num_comms other tasks some tuples (default 2)", &num_comms);

  int num_tuples = 4;
  opts.addOpt<int>("num_tuples,t",
        "each task will send to some task about num_tuples tuples (default 4)", &num_tuples);

  int reportrank = size + 1;
  opts.addOpt<int>("reporting_rank,r",
      "this rank will report the tuples sent and the tuples received; it could be higher than num_procs, then no reporting"
      ,&reportrank);

  opts.parseCommandLine(argc, argv);

  if (rank == reportrank || (reportrank >= size && 0 == rank)) {
    cout << " There are " << size << " tasks in example.\n";
    cout << " We will send groups of " << num_tuples << " from each task towards " <<
        num_comms << " other tasks.\n";
  }

  // Send some data from proc i to i + n/2, also to i + n/2 + 1 modulo n, where n is num procs

  gs_data::crystal_data *cd = pc.crystal_router();

  long total_n_tuples = num_comms * num_tuples;

  // Vary the number of tasks to send to, and the number of tuples to send
  if (rank < size/2)
    num_comms--;
  else
    num_comms++;

  if (rank < size/3)
    num_tuples *= 2;
  else if (rank > size - size/3)
    num_tuples /= 2;

  TupleList tl;
  // At most num_tuples* num_comms to send
  // We do a preallocate with this; some tuples on some processors might need more memory, to be able
  // to grow locally; Some tasks might receive more tuples though, and in the process, some might grow more than
  // others. By doing these logP sends/receives, we do not grow local memory too much.
  tl.initialize(1, 1, 0, 1, num_tuples*num_comms);
  tl.enableWriteAccess();
  // Form num_tuples*num_comms tuples, send to various ranks
  unsigned int n = tl.get_n();
  for (int i = 0; i < num_comms; i++) {
    int sendTo = rank + i*size/2 + 1; // Spread out the send to, for a stress-like test
    sendTo = sendTo % size;//
    long intToSend = 1000*rank + 100000*sendTo;
    for (int j = 0; j < num_tuples; j++) {
      n = tl.get_n();
      tl.vi_wr[n] = sendTo;
      tl.vl_wr[n] = intToSend + j;
      tl.vr_wr[n] = 10000.*rank + j;
      tl.inc_n();
    }
  }

  if (rank == reportrank) {
    cout << "rank " << rank << "\n";
    tl.print(" before sending");
  }

  clock_t tt = clock();
  // All communication happens here; no mpi calls for the user
  ErrorCode rval = cd->gs_transfer(1, tl, 0);MB_CHK_SET_ERR(rval, "Error in tuple transfer");

  double secs = 0;
  if (rank == reportrank || (reportrank >= size && 0 == rank)) {
    secs = (clock() - tt) / (double) CLOCKS_PER_SEC;
  }
  if (rank == reportrank) {
    cout << "rank " << rank << "\n";
    tl.print(" after transfer");
  }

  // Check that all tuples received have the form 10000*rank + 100*from
  unsigned int received = tl.get_n();
  for (int i = 0; i < (int)received; i++) {
    int from = tl.vi_rd[i];
    long valrec = tl.vl_rd[i];
    int remainder = valrec - 100000*rank - 1000*from;
    if (remainder < 0 || remainder >= num_tuples*4)
      cout << " error: tuple " << i << " received at proc rank " << rank << " from proc " << from <<
      " has value " << valrec << " remainder " << remainder << "\n";
  }

  if (rank == reportrank || (reportrank >= size && 0 == rank)) {
    cout << "communication of about " << total_n_tuples << " tuples/per proc took " << secs << " seconds" << "\n";
    tt = clock();
  }

  // Finalize error handler, required for this example (not using a moab instance)
  MBErrorHandler_Finalize();

  MPI_Finalize();
#else
  std::cout<<" Build with MPI for this example to work\n";
#endif

  return 0;
}
