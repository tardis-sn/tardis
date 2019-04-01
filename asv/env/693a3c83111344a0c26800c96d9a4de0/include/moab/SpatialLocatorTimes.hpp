/**\file SpatialLocatorTimes.hpp
 * \class moab::SpatialLocatorTimes
 * \brief Statistics for spatial location
 *
 * Class to accumulate statistics on performance of spatial location.  This structure stores
 * only local (single proc) statistics, but provides functions for accumulating max/min/avg
 * time properties for performance reporting.
 *
 * Similar to TreeStats, this class is very lightweight, with most variables publicly accessible.
 * 
 */

#ifndef SPATIALLOCATORTIMES_HPP
#define SPATIALLOCATORTIMES_HPP

#include "moab/Interface.hpp"

#include <iostream>

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

namespace moab 
{
class SpatialLocatorTimes
{
public:
    /* constructor
     */
  SpatialLocatorTimes() {reset();}
  
    /* \brief Reset all stats defined here
     */
  void reset();
  
    /* \brief Output header for stats names
     */
  void output_header(bool print_endl = false) const;

    /* \brief Output stats all on one line
     */
  void output(bool print_head = false, bool print_endl = false) const;

#ifdef MOAB_HAVE_MPI
    /* \brief Accumulate times over all processors into provided array
     * Max and min is accumulated over all processors, onto root, for each stat.  If avg_stats is non-NULL,
     * all stats are gathered to root so average can be taken too.
     *
     * This function must be called collectively on the communicator
     * \param comm MPI communictor over which this accumulation is done
     * \param max_times Array of size NUM_STATS into which max times are inserted
     * \param min_times Array of size NUM_STATS into which min times are inserted
     * \param avg_times Array of size NUM_STATS into which avg times are inserted; if NULL, no avg times are computed.

     */
  ErrorCode accumulate_times(MPI_Comm comm, double *max_times, double *min_times, 
                             double *avg_times = NULL);
#endif

    /* \brief Enumeration to identify fields in performance data
     */
  enum {INTMED_INIT = 0,       // time to compute intermediate partition, incl global bounding box         
        INTMED_SEND,           // time to send search points from target to intermediate parts             
        INTMED_SEARCH,         // time to find candidate src boxes for search points on intermidiate procs 
        SRC_SEND,              // time to send search points to src procs                                  
        SRC_SEARCH,            // time to search local box/elements on src procs                           
        TARG_RETURN,           // time to return point location data to target procs
        TARG_STORE,            // time to store point location into local SpatialLocator object
        NUM_STATS              // number of stats, useful for array sizing and terminating loops over stats
  };
  
  double slTimes[NUM_STATS];
};

inline void SpatialLocatorTimes::reset() 
{
  for (int i = 0; i < NUM_STATS; i++) slTimes[i] = 0.0;
}

#ifdef MOAB_HAVE_MPI
inline ErrorCode SpatialLocatorTimes::accumulate_times(MPI_Comm comm, 
                                                       double *min_times, double *max_times, double *avg_times) 
{
  ErrorCode rval = MB_SUCCESS;
  int success = MPI_Reduce(slTimes, min_times, NUM_STATS, MPI_DOUBLE, MPI_MIN, 0, comm);
  if (!success) rval = MB_FAILURE;
  
  success = MPI_Reduce(slTimes, max_times, NUM_STATS, MPI_DOUBLE, MPI_MAX, 0, comm);
  if (!success) rval = MB_FAILURE;
  
  if (avg_times) {
    int sz, rank;
    MPI_Comm_size(comm, &sz);
    MPI_Comm_rank(comm, &rank);
    std::vector<double> all_times;
    if (!rank) all_times.resize(NUM_STATS*sz+NUM_STATS);
    success = MPI_Gather(slTimes, NUM_STATS, MPI_DOUBLE, (rank ? NULL : &all_times[0]), NUM_STATS, MPI_DOUBLE, 0, comm);
    if (!success) rval = MB_FAILURE;
    if (!rank) {
      std::fill(avg_times, avg_times+NUM_STATS, 0.0);
      for (int p = 0; p < sz; p++) {
        for (int s = 0; s < NUM_STATS; s++)
          avg_times[s] += all_times[p*NUM_STATS+s];
      }
    
      for (int s = 0; s <= NUM_STATS; s++) avg_times[s] /= (double)sz;
    }
  }
  
  return rval;
}
#endif

  /* \brief Output stats all on one line
   */
inline void SpatialLocatorTimes::output_header(bool print_endl) const 
{
  std::cout << "Intmed_init Intmed_send Intmed_search src_send src_search targ_return targ_store";
  if (print_endl) std::cout << std::endl;
}

  /* \brief Output stats all on one line
   */
inline void SpatialLocatorTimes::output(bool print_head, bool print_endl) const 
{
  if (print_head) output_header(true);
  for (int i = 0; i < NUM_STATS; i++) std::cout << slTimes[i] << " ";
  
  if (print_endl) std::cout << std::endl;
}
    

} // namespace moab

#endif
    
  
