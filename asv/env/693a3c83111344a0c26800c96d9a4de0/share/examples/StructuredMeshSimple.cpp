/** @example StructuredMeshSimple.cpp
 * \brief Show creation and query of structured mesh, serial or parallel, through MOAB's structured mesh interface.
 * This is an example showing creation and query of a 3D structured mesh.  In serial, a single N*N*N block of elements
 * is created; in parallel, each proc gets an N*N*N block, with blocks arranged in a 1d column, sharing vertices
 * and faces at their interfaces (proc 0 has no left neighbor and proc P-1 no right neighbor).
 * Each square block of hex elements is then referenced by its ijk parameterization.
 * 1D and 2D examples could be made simply by changing the dimension parameter passed into the MOAB functions. \n
 *
 * <b>This example </b>:
 *    -# Instantiate MOAB and get the structured mesh interface
 *    -# Decide what the local parameters of the mesh will be, based on parallel/serial and rank.
 *    -# Create a N^d structured mesh, which includes (N+1)^d vertices and N^d elements.
 *    -# Get the vertices and elements from moab and check their numbers against (N+1)^d and N^d, resp.
 *    -# Loop over elements in d nested loops over i, j, k; for each (i,j,k):
 *      -# Get the element corresponding to (i,j,k)
 *      -# Get the connectivity of the element
 *      -# Get the coordinates of the vertices comprising that element
 *    -# Release the structured mesh interface and destroy the MOAB instance
 *
 * <b> To run: </b> ./StructuredMeshSimple [d [N] ] \n
 * (default values so can run w/ no user interaction)
 */

#include "moab/Core.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CN.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif
#include <iostream>
#include <vector>

using namespace moab;
using namespace std;

int main(int argc, char **argv) 
{
  int N = 10, dim = 3;

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  ProgOptions opts;
  opts.addOpt<int>(string("dim,d"), string("Dimension of mesh (default=3)"),
                   &dim);
  opts.addOpt<int>(string(",n"), string("Number of elements on a side (default=10)"),
                   &N);
  opts.parseCommandLine(argc, argv);

  // 0. Instantiate MOAB and get the structured mesh interface
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;
  ScdInterface *scdiface;
  ErrorCode rval = mb->query_interface(scdiface);MB_CHK_ERR(rval); // Get a ScdInterface object through moab instance

  // 1. Decide what the local parameters of the mesh will be, based on parallel/serial and rank.
#ifdef MOAB_HAVE_MPI
  int rank = 0, nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int ilow = rank * N, ihigh = ilow + N;
#else
  int rank = 0;
  int ilow = 0, ihigh = N;
#endif  

  // 2. Create a N^d structured mesh, which includes (N+1)^d vertices and N^d elements.
  ScdBox *box;
  rval = scdiface->construct_box(HomCoord(ilow, (dim > 1? 0 : -1), (dim > 2? 0 : -1)), // Use in-line logical tests to handle dimensionality
                                 HomCoord(ihigh, (dim > 1? N : -1), (dim > 2? N : -1)),
                                 NULL, 0, // NULL coords vector and 0 coords (don't specify coords for now)
                                 box);MB_CHK_ERR(rval); // box is the structured box object providing the parametric
                                                     // structured mesh interface for this rectangle of elements

  // 3. Get the vertices and elements from moab and check their numbers against (N+1)^d and N^d, resp.
  Range verts, elems;
  rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_ERR(rval); // First '0' specifies "root set", or entire MOAB instance, second the entity dimension being requested
  rval = mb->get_entities_by_dimension(0, dim, elems);MB_CHK_ERR(rval);

#define MYSTREAM(a) if (!rank) cout << a << endl

  if (pow(N, dim) == (int) elems.size() && pow(N+1, dim) == (int) verts.size()) { // Expected #e and #v are N^d and (N+1)^d, resp.
#ifdef MOAB_HAVE_MPI
    MYSTREAM("Proc 0: ");
#endif
    MYSTREAM("Created " << elems.size() << " " << CN::EntityTypeName(mb->type_from_handle(*elems.begin())) 
             << " elements and " << verts.size() << " vertices." << endl);
  }
  else
    cout << "Created the wrong number of vertices or hexes!" << endl;

  // 4. Loop over elements in 3 nested loops over i, j, k; for each (i,j,k):
  vector<double> coords(3*pow(N+1, dim));
  vector<EntityHandle> connect;
  for (int k = 0; k < (dim > 2? N : 1); k++) {
    for (int j = 0; j < (dim > 1? N : 1); j++) {
      for (int i = 0; i < N - 1; i++) {
        // 4a. Get the element corresponding to (i,j,k)
        EntityHandle ehandle = box->get_element(i, j, k);
        if (0 == ehandle) return MB_FAILURE;
        // 4b. Get the connectivity of the element
        rval = mb->get_connectivity(&ehandle, 1, connect);MB_CHK_ERR(rval); // Get the connectivity, in canonical order
        // 4c. Get the coordinates of the vertices comprising that element
        rval = mb->get_coords(&connect[0], connect.size(), &coords[0]);MB_CHK_ERR(rval); // Get the coordinates of those vertices
      }
    }
  }

  // 5. Release the structured mesh interface and destroy the MOAB instance
  mb->release_interface(scdiface); // Tell MOAB we're done with the ScdInterface

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
