/** @example LloydRelaxation.cpp \n
 * \brief Perform Lloyd relaxation on a mesh and its dual \n
 * <b>To run</b>: mpiexec -np <np> LloydRelaxation [filename]\n
 *
 * Briefly, Lloyd relaxation is a technique to smooth out a mesh.  The centroid of each cell is computed from its
 * vertex positions, then vertices are placed at the average of their connected cells' centroids.
 *
 * In the parallel algorithm, an extra ghost layer of cells is exchanged.  This allows us to compute the centroids
 * for boundary cells on each processor where they appear; this eliminates the need for one round of data exchange
 * (for those centroids) between processors.  New vertex positions must be sent from owning processors to processors
 * sharing those vertices.  Convergence is measured as the maximum distance moved by any vertex.  
 * 
 * In this implementation, a fixed number of iterations is performed.  The final mesh is output to 'lloydfinal.h5m'
 * in the current directory (H5M format must be used since the file is written in parallel).
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#  include "moab/ParallelComm.hpp"
#  include "MBParallelConventions.h"
#endif
#include "moab/Skinner.hpp"
#include "moab/CN.hpp"
#include "moab/CartVect.hpp"
#include <iostream>
#include <sstream>

using namespace moab;
using namespace std;

string test_file_name = string(MESH_DIR) + string("/surfrandomtris-4part.h5m");

ErrorCode perform_lloyd_relaxation(Interface *mb, Range &verts, Range &cells, Tag fixed, 
                                   int num_its, int report_its);

int main(int argc, char **argv)
{
  int num_its = 10;
  int report_its = 1;

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }  

  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

#ifdef MOAB_HAVE_MPI
  // Get the ParallelComm instance
  ParallelComm *pcomm = new ParallelComm(mb, MPI_COMM_WORLD);
  int nprocs = pcomm->size();
#else
  int nprocs = 1;
#endif
  string options;
  if (nprocs > 1) // If reading in parallel, need to tell it how
    options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=2.0.1;DEBUG_IO=0;DEBUG_PIO=0";

  // Load the test file with specified options
  ErrorCode rval = mb->load_file(test_file_name.c_str(), 0, options.c_str());MB_CHK_ERR(rval);

  // Make tag to specify fixed vertices, since it's input to the algorithm; use a default value of non-fixed
  // so we only need to set the fixed tag for skin vertices
  Tag fixed;
  int def_val = 0;
  rval = mb->tag_get_handle("fixed", 1, MB_TYPE_INTEGER, fixed, MB_TAG_CREAT | MB_TAG_DENSE, &def_val);MB_CHK_ERR(rval);

  // Get all vertices and faces
  Range verts, faces, skin_verts;
  rval = mb->get_entities_by_type(0, MBVERTEX, verts);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_dimension(0, 2, faces);MB_CHK_ERR(rval);

  // Get the skin vertices of those faces and mark them as fixed; we don't want to fix the vertices on a
  // part boundary, but since we exchanged a layer of ghost faces, those vertices aren't on the skin locally
  // ok to mark non-owned skin vertices too, I won't move those anyway
  // use MOAB's skinner class to find the skin
  Skinner skinner(mb);
  rval = skinner.find_skin(0, faces, true, skin_verts);MB_CHK_ERR(rval); // 'true' param indicates we want vertices back, not faces

  vector<int> fix_tag(skin_verts.size(), 1); // Initialized to 1 to indicate fixed
  rval = mb->tag_set_data(fixed, skin_verts, &fix_tag[0]);MB_CHK_ERR(rval);

  // Now perform the Lloyd relaxation
  rval = perform_lloyd_relaxation(mb, verts, faces, fixed, num_its, report_its);MB_CHK_ERR(rval);

  // Delete fixed tag, since we created it here
  rval = mb->tag_delete(fixed);MB_CHK_ERR(rval);

  // Output file, using parallel write

#ifdef MOAB_HAVE_MPI
  options = "PARALLEL=WRITE_PART";
#endif

  rval = mb->write_file("lloydfinal.h5m", NULL, options.c_str());MB_CHK_ERR(rval);

  // Delete MOAB instance
  delete mb;

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

ErrorCode perform_lloyd_relaxation(Interface *mb, Range &verts, Range &faces, Tag fixed, 
                                   int num_its, int report_its) 
{
  ErrorCode rval;

#ifdef MOAB_HAVE_MPI
  ParallelComm *pcomm = ParallelComm::get_pcomm(mb, 0);
  int nprocs = pcomm->size();
#else
  int nprocs = 1;
#endif

  // Perform Lloyd relaxation:
  // 1. Setup: set vertex centroids from vertex coords; filter to owned verts; get fixed tags

  // Get all verts coords into tag; don't need to worry about filtering out fixed verts,
  // we'll just be setting to their fixed coords
  vector<double> vcentroids(3*verts.size());
  rval = mb->get_coords(verts, &vcentroids[0]);MB_CHK_ERR(rval);

  Tag centroid;
  rval = mb->tag_get_handle("centroid", 3, MB_TYPE_DOUBLE, centroid, MB_TAG_CREAT | MB_TAG_DENSE);MB_CHK_ERR(rval);
  rval = mb->tag_set_data(centroid, verts, &vcentroids[0]);MB_CHK_ERR(rval);

  // Filter verts down to owned ones and get fixed tag for them
  Range owned_verts, shared_owned_verts;
  if (nprocs > 1) {
#ifdef MOAB_HAVE_MPI
    rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_verts);MB_CHK_ERR(rval);
#endif
  }
  else
    owned_verts = verts;
  vector<int> fix_tag(owned_verts.size());
  rval = mb->tag_get_data(fixed, owned_verts, &fix_tag[0]);MB_CHK_ERR(rval);

  // Now fill vcentroids array with positions of just owned vertices, since those are the ones
  // we're actually computing
  vcentroids.resize(3*owned_verts.size());
  rval = mb->tag_get_data(centroid, owned_verts, &vcentroids[0]);MB_CHK_ERR(rval);

#ifdef MOAB_HAVE_MPI
  // Get shared owned verts, for exchanging tags
  rval = pcomm->get_shared_entities(-1, shared_owned_verts, 0, false, true);MB_CHK_ERR(rval);
  // Workaround: if no shared owned verts, put a non-shared one in the list, to prevent exchanging tags
  // for all shared entities
  if (shared_owned_verts.empty()) shared_owned_verts.insert(*verts.begin());
#endif

  // Some declarations for later iterations
  vector<double> fcentroids(3*faces.size()); // fcentroids for face centroids
  vector<double> ctag(3*CN::MAX_NODES_PER_ELEMENT); // Temporary coordinate storage for verts bounding a face
  const EntityHandle *conn; // const ptr & size to face connectivity
  int nconn;
  Range::iterator fit, vit; // For iterating over faces, verts
  int f, v; // For indexing into centroid vectors
  vector<EntityHandle> adj_faces; // Used in vertex iteration

  // 2. For num_its iterations:
  for (int nit = 0; nit < num_its; nit++) {
    double mxdelta = 0.0;

    // 2a. For each face: centroid = sum(vertex centroids)/num_verts_in_cell
    for (fit = faces.begin(), f = 0; fit != faces.end(); ++fit, f++) {
      // Get verts for this face
      rval = mb->get_connectivity(*fit, conn, nconn);MB_CHK_ERR(rval);
      // Get centroid tags for those verts
      rval = mb->tag_get_data(centroid, conn, nconn, &ctag[0]);MB_CHK_ERR(rval);
      fcentroids[3*f + 0] = fcentroids[3*f + 1] = fcentroids[3*f + 2] = 0.0;
      for (v = 0; v < nconn; v++) {
        fcentroids[3*f + 0] += ctag[3*v + 0];
        fcentroids[3*f + 1] += ctag[3*v + 1];
        fcentroids[3*f + 2] += ctag[3*v + 2];
      }
      for (v = 0; v < 3; v++) fcentroids[3*f + v] /= nconn;
    }
    rval = mb->tag_set_data(centroid, faces, &fcentroids[0]);MB_CHK_ERR(rval);

    // 2b. For each owned vertex:
    for (vit = owned_verts.begin(), v = 0; vit != owned_verts.end(); ++vit, v++) {
      // if !fixed
      if (fix_tag[v]) continue;
      // vertex centroid = sum(cell centroids)/ncells
      adj_faces.clear();
      rval = mb->get_adjacencies(&(*vit), 1, 2, false, adj_faces);MB_CHK_ERR(rval);
      rval = mb->tag_get_data(centroid, &adj_faces[0], adj_faces.size(), &fcentroids[0]);MB_CHK_ERR(rval);
      double vnew[] = {0.0, 0.0, 0.0};
      for (f = 0; f < (int)adj_faces.size(); f++) {
        vnew[0] += fcentroids[3*f + 0];
        vnew[1] += fcentroids[3*f + 1];
        vnew[2] += fcentroids[3*f + 2];
      }
      for (f = 0; f < 3; f++) vnew[f] /= adj_faces.size();
      double delta = (CartVect(vnew) - CartVect(&vcentroids[3*v])).length();
      mxdelta = std::max(delta, mxdelta);
      for (f = 0; f < 3; f++) vcentroids[3*v + f] = vnew[f];
    }

    // Set the centroid tag; having them only in vcentroids array isn't enough, as vertex centroids are
    // accessed randomly in loop over faces
    rval = mb->tag_set_data(centroid, owned_verts, &vcentroids[0]);MB_CHK_ERR(rval);

    // 2c. Exchange tags on owned verts
    if (nprocs > 1) {
#ifdef MOAB_HAVE_MPI
      rval = pcomm->exchange_tags(centroid, shared_owned_verts);MB_CHK_ERR(rval);
#endif
    }

    if (!(nit%report_its)) {
      // Global reduce for maximum delta, then report it
      double global_max = mxdelta;
      int myrank = 0;
#ifdef MOAB_HAVE_MPI
      if (nprocs > 1)
        MPI_Reduce(&mxdelta, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, pcomm->comm());
      myrank = pcomm->rank();
#endif
      if (1 == nprocs || 0 == myrank)
        cout << "Max delta = " << global_max << endl;
    }
  }

  // Write the tag back onto vertex coordinates
  rval = mb->set_coords(owned_verts, &vcentroids[0]);MB_CHK_ERR(rval);

  // Delete the centroid tag, since we don't need it anymore
  rval = mb->tag_delete(centroid);MB_CHK_ERR(rval);

  return MB_SUCCESS;
}
