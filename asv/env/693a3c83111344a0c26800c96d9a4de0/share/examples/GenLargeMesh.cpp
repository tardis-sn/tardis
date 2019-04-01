/** @example GenLargeMesh.cpp \n
 * \brief Create a large structured mesh, partitioned \n
 *
 *  It shows how to create a mesh on the fly, on multiple processors
 *  Each processor will create its version of a block mesh, partitioned
 *  as AxBxC blocks. Each block will be with blockSize^3 hexahedrons, and
 *  will get a different PARALLEL_PARTITION tag
 *
 *  The number of tasks will be MxNxK, and it must match the mpi size
 *  Each task p will generate its mesh at location (m,n,k), and it is
 *  lexicographic ordering: rank = m + n * M + k * M * N
 *
 *  By default M=1, N=1, K=1, so by default it should be launched on 1 proc
 *  By default, blockSize is 4, and A=2, B=2, C=2, so each task will generate locally
 *  blockSize^3 x A x B x C hexahedrons (value = 64x8 = 512 hexas, in 8 partitions)
 *  (if -t, multiple by 6 for total number of cells/tets)
 *  The total number of partitions will be A*B*C*M*N*K (default 8)
 *
 *  Each part in partition will get a proper tag, and the numbering is in
 *  lexicographic ordering; x direction varies first, then y, then z.
 *  The same principle is used for global id numbering of the nodes and cells.
 *  (x varies first)
 *
 *  The vertices will get a proper global id, which will be used to resolve the
 *  shared entities
 *  The output will be written in parallel, and we will try sizes as big as we can
 *  (up to a billion vertices, because we use int for global ids)
 *
 *  Within each partition, the hexas entity handles will be contiguous, and also the
 *  vertices; The global id will be determined easily, for a vertex, but the entity
 *  handle space will be more interesting to handle, within a partition (we want
 *  contiguous handles within a partition). After merging the vertices, some fragmentation
 *  will occur in the vertices handle space, within each partition.
 *
 *  To run: ./GenLargeMesh
 *
 *  When launched on more procs, you have to make sure
 *  num procs = M*N*K
 *
 *  So you can launch with
 *  mpiexec -np 8 ./GenLargeMesh -M 2 -N 2 -K 2
 *
 *  We also added -q option; it works now only for hexa mesh, it will generate
 *  quadratic hex27 elements
 *
 *  -t option will generate tetrahedrons instead of hexahedra. Each hexahedra is
 *  decomposed into 6 tetrahedrons.
 *
 *  -f option will also generate all edges and faces in the model.
 *  -w will use a newer merging method locally. Merging is necessary to merge
 *  vertices on the local task, and the new method does not use a searching tree,
 *  but rather the global id set on the vertices in a consistent manner
 *
 *  -d and -i options can be used to add some artificial tags on the model;
 *  you can have multiple -d and -i options; -i <tag_name> will set an integer
 *  tag with name tag_name on the vertices; -d < tag_name2> will generate
 *  double tags on cells (3d elements). You can have multiple tags, like
 *  -i tag1 -i tag2 -i tag3 -d tag4
 *
 *  -x, -y, -z options will control the geometric dimensions of the final mesh, in
 *  x, y and z directions.
 *
 *  -o <out_file> controls the name of the output file; it needs to have extension h5m,
 *  because the file is written in parallel.
 *
 *  -k will keep the edges and faces that are generated as part of resolving shared entities
 *  (by default these edges and faces are removed); when -f option is used, the -k option is
 *  enabled too (so no faces and edges are deleted)
 *
 */

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/MergeMesh.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "moab/ParallelMergeMesh.hpp"
#include "MBParallelConventions.h"
#endif
#include "moab/ReadUtilIface.hpp"

#include <time.h>
#include <iostream>
#include <vector>

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
  int A = 2, B = 2, C = 2, M = 1, N = 1, K = 1;
  int blockSize = 4;
  double xsize = 1., ysize = 1., zsize = 1.; // The size of the region
  int GL = 0; // number of ghost layers

  bool newMergeMethod = false;
  bool quadratic = false;
  bool keep_skins = false;
  bool tetra = false;
  bool adjEnts = false;
  bool parmerge = false;
  bool nosave = false;

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  ProgOptions opts;

  opts.addOpt<int>(string("blockSize,b"),
      string("Block size of mesh (default=4)"), &blockSize);
  opts.addOpt<int>(string("xproc,M"),
      string("Number of processors in x dir (default=1)"), &M);
  opts.addOpt<int>(string("yproc,N"),
      string("Number of processors in y dir (default=1)"), &N);
  opts.addOpt<int>(string("zproc,K"),
      string("Number of processors in z dir (default=1)"), &K);

  opts.addOpt<int>(string("xblocks,A"),
      string("Number of blocks on a task in x dir (default=2)"), &A);
  opts.addOpt<int>(string("yblocks,B"),
      string("Number of blocks on a task in y dir (default=2)"), &B);
  opts.addOpt<int>(string("zblocks,C"),
      string("Number of blocks on a task in x dir (default=2)"), &C);

  opts.addOpt<double>(string("xsize,x"),
      string("Total size in x direction (default=1.)"), &xsize);
  opts.addOpt<double>(string("ysize,y"),
      string("Total size in y direction (default=1.)"), &ysize);
  opts.addOpt<double>(string("zsize,z"),
      string("Total size in z direction (default=1.)"), &zsize);

  opts.addOpt<void>("newMerge,w", "use new merging method", &newMergeMethod);

  opts.addOpt<void>("quadratic,q", "use hex 27 elements", &quadratic);

  opts.addOpt<void>("keep_skins,k", "keep skins with shared entities", &keep_skins);

  opts.addOpt<void>("tetrahedrons,t", "generate tetrahedrons", &tetra);

  opts.addOpt<void>("faces_edges,f", "create all faces and edges", &adjEnts);

  opts.addOpt<int>(string("ghost_layers,g"),
  string("Number of ghost layers (default=0)"), &GL);

  vector<string> intTagNames;
  string firstIntTag;
  opts.addOpt<string>("int_tag_vert,i", "add integer tag on vertices", &firstIntTag);

  vector<string> doubleTagNames;
  string firstDoubleTag;
  opts.addOpt<string>("double_tag_cell,d", "add double tag on cells", &firstDoubleTag);

  string outFileName = "GenLargeMesh.h5m";
  opts.addOpt<string>("outFile,o", "Specify the output file name string (default GenLargeMesh.h5m)", &outFileName);

#ifdef MOAB_HAVE_HDF5_PARALLEL
  bool readb = false;
  opts.addOpt<void>("readback,r", "read back the generated mesh", &readb);

  bool readAndGhost = false;
  opts.addOpt<void>("readAndGhost,G", "read back the generated mesh and ghost one layer", &readAndGhost);
#endif

  opts.addOpt<void>("parallel_merge,p", "use parallel mesh merge, not vertex ID based merge", &parmerge);

  opts.addOpt<void>("no_save,n", "do not save the file", &nosave);

  opts.parseCommandLine(argc, argv);

  opts.getOptAllArgs("int_tag_vert,i", intTagNames);
  opts.getOptAllArgs("double_tag_cell,d", doubleTagNames);

  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb) {
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  ReadUtilIface* iface;
  ErrorCode rval = mb->query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

  int rank=0, size=1;

#ifdef MOAB_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if (M*N*K != size) {
    if (0 == rank)
      cout << "M*N*K = " << M*N*K << " != size = " << size << "\n";
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  if (adjEnts)
    keep_skins = true; // Do not delete anything

  // Determine m, n, k for processor rank
  int m, n, k;
  k = rank / (M*N);
  int leftover = rank % (M*N);
  n = leftover / M;
  m = leftover % M;
  if (rank == size - 1)
    cout << "m, n, k for last rank: " << m << " " << n << " " << k << "\n";

  // So there are a total of M * A * blockSize elements in x direction (so M * A * blockSize + 1 verts in x direction)
  // So there are a total of N * B * blockSize elements in y direction (so N * B * blockSize + 1 verts in y direction)
  // So there are a total of K * C * blockSize elements in z direction (so K * C * blockSize + 1 verts in z direction)

  // There are (M * A blockSize)       * (N * B * blockSize)     * (K * C * blockSize)     hexas
  // There are (M * A * blockSize + 1) * (N * B * blockSize + 1) * (K * C * blockSize + 1) vertices
  // x is the first dimension that varies

  clock_t tt = clock();

  // Used for nodes increments
  int q = (quadratic)? 2 : 1;
  // Used for element increments
  int factor = (tetra)? 6 : 1;

  double dx = xsize / (A*M*blockSize*q); // Distance between 2 nodes in x direction
  double dy = ysize / (B*N*blockSize*q); // Distance between 2 nodes in y direction
  double dz = zsize / (C*K*blockSize*q); // Distance between 2 nodes in z direction

  int NX = (q * M * A * blockSize + 1);
  int NY = (q * N * B * blockSize + 1);
  int nex = M * A * blockSize; // Number of elements in x direction, used for global id on element
  int ney = N * B * blockSize; // Number of elements in y direction ...
  // int NZ = (K * C * blockSize + 1); // Not used
  int blockSize1 = q*blockSize + 1; // Used for vertices
  long num_total_verts = (long) NX * NY * (K * C * blockSize + 1);
  if (0 == rank)
  {
    cout << "Generate mesh on " << size << " processors \n";
    std::cout << "Total number of vertices: " << num_total_verts << "\n";
  }
  //int xstride = 1;
  int ystride = blockSize1;

  int zstride = blockSize1 * blockSize1;
  // Generate the block at (a, b, c); it will represent a partition, it will get a partition tag

  Tag global_id_tag;
  rval = mb->tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER,
                             global_id_tag);MB_CHK_SET_ERR(rval, "Can't get global id tag");

  // set global ids
  Tag new_id_tag;
  if (!parmerge)
  {
    rval = mb->tag_get_handle("HANDLEID", sizeof(long), MB_TYPE_OPAQUE,
      new_id_tag, MB_TAG_CREAT|MB_TAG_DENSE);MB_CHK_SET_ERR(rval, "Can't get handle id tag");
  }
  Tag part_tag;
  int dum_id = -1;
  rval = mb->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                             part_tag, MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id);MB_CHK_SET_ERR(rval, "Can't get parallel partition tag");

  // Create tags on vertices and cells, look in the list of options
  vector<Tag> intTags(intTagNames.size());
  vector<Tag> doubleTags(doubleTagNames.size());
  for (size_t i = 0; i < intTagNames.size(); i++) {
    rval = mb->tag_get_handle(intTagNames[i].c_str(), 1, MB_TYPE_INTEGER, intTags[i],
                              MB_TAG_CREAT | MB_TAG_DENSE, &dum_id);MB_CHK_SET_ERR(rval, "Can't create integer tag");
  }

  double defval = 0.;
  for (size_t i = 0; i < doubleTagNames.size(); i++) {
    rval = mb->tag_get_handle(doubleTagNames[i].c_str(), 1, MB_TYPE_DOUBLE, doubleTags[i],
                              MB_TAG_CREAT | MB_TAG_DENSE, &defval);MB_CHK_SET_ERR(rval, "Can't create double tag");
  }
  Range wsets; // write only part sets
  for (int a = 0; a < A; a++) {
    for (int b = 0; b < B; b++) {
      for (int c = 0; c < C; c++) {
        // We will generate (q*block + 1)^3 vertices, and block^3 hexas; q is 1 for linear, 2 for quadratic
        // the global id of the vertices will come from m, n, k, a, b, c
        // x will vary from  m*A*q*block + a*q*block to m*A*q*block + (a+1)*q*block etc;
        int num_nodes = blockSize1 * blockSize1 * blockSize1;

        vector<double*> arrays;
        EntityHandle startv;
        rval = iface->get_node_coords(3, num_nodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

        // Will start with the lower corner:
        int x = m*A*q*blockSize + a*q*blockSize;
        int y = n*B*q*blockSize + b*q*blockSize;
        int z = k*C*q*blockSize + c*q*blockSize;
        int ix = 0;
        vector<int> gids(num_nodes);
        vector<long> lgids(num_nodes);
        Range verts(startv, startv + num_nodes - 1);
        for (int kk = 0; kk < blockSize1; kk++) {
          for (int jj = 0; jj < blockSize1; jj++) {
            for (int ii = 0; ii < blockSize1; ii++) {
              arrays[0][ix] = (x+ii) * dx;
              arrays[1][ix] = (y+jj) * dy;
              arrays[2][ix] = (z+kk) * dz;
              gids[ix] = 1 + (x+ii) + (y+jj) * NX + (z+kk) * (NX*NY);
              if (!parmerge)
                lgids[ix] = 1 + (x+ii) + (y+jj) * NX + (long)(z+kk) * (NX*NY);
              // Set int tags, some nice values?
              EntityHandle v = startv + ix;
              for (size_t i = 0; i < intTags.size(); i++) {
                int valv = gids[ix]/2 + 3 + i*1000;
                rval = mb->tag_set_data(intTags[i], &v, 1, &valv);MB_CHK_SET_ERR(rval, "Can't set integer tag on a vertex");
              }
              ix++;
            }
          }
        }

        rval = mb->tag_set_data(global_id_tag, verts, &gids[0]);MB_CHK_SET_ERR(rval, "Can't set global ids to vertices");
        if (!parmerge)
        {
          rval = mb->tag_set_data(new_id_tag, verts, &lgids[0]);MB_CHK_SET_ERR(rval, "Can't set the new handle id tags");
        }
        int num_hexas = blockSize * blockSize * blockSize;
        int num_el = num_hexas * factor;

        EntityHandle starte; // Connectivity
        EntityHandle* conn;
        int num_v_per_elem = 8;
        if (quadratic) {
          num_v_per_elem = 27;
          rval = iface->get_element_connect(num_el, 27, MBHEX, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
        }
        else if (tetra) {
          num_v_per_elem = 4;
          rval = iface->get_element_connect(num_el, 4, MBTET, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
        }
        else {
          rval = iface->get_element_connect(num_el, 8, MBHEX, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");
        }

        Range cells(starte, starte + num_el - 1); // Should be elements
        // Fill cells
        ix = 0;
        // Identify the elements at the lower corner, for their global ids
        int xe = m*A*blockSize + a*blockSize;
        int ye = n*B*blockSize + b*blockSize;
        int ze = k*C*blockSize + c*blockSize;
        gids.resize(num_el);
        lgids.resize(num_el);
        int ie = 0; // Index now in the elements, for global ids
        for (int kk = 0; kk < blockSize; kk++) {
          for (int jj = 0; jj < blockSize; jj++) {
            for (int ii = 0; ii < blockSize; ii++) {
              EntityHandle corner = startv + q * ii + q * jj * ystride + q * kk * zstride;
              // These could overflow for large numbers
              gids[ie] = 1 + ((xe+ii) + (ye+jj) * nex + (ze+kk) * (nex*ney))*factor ; // 6 more for tetra
              lgids[ie] = 1 + ((xe+ii) + (ye+jj) * nex + (long)(ze+kk) * (nex*ney))*factor ; // 6 more for tetra
              EntityHandle eh = starte + ie;
              for (size_t i = 0; i < doubleTags.size(); i++) {
                double valv = gids[ie]/30. + i*5000.;
                rval = mb->tag_set_data(doubleTags[i], &eh, 1, &valv);MB_CHK_SET_ERR(rval, "Can't set double tag on an element");
              }
              ie++;
              if (quadratic) {
  //                    4   ----- 19   -----  7
  //                .   |                 .   |
  //            16         25         18      |
  //         .          |          .          |
  //      5   ----- 17   -----  6             |
  //      |            12       | 23         15
  //      |                     |             |
  //      |     20      |  26   |     22      |
  //      |                     |             |
  //     13         21  |      14             |
  //      |             0   ----- 11   -----  3
  //      |         .           |         .
  //      |      8         24   |     10
  //      |  .                  |  .
  //      1   -----  9   -----  2
  //
                conn[ix] =    corner;
                conn[ix+1] =  corner + 2;
                conn[ix+2] =  corner + 2 + 2 * ystride;
                conn[ix+3] =  corner +     2 * ystride;
                conn[ix+4] =  corner                   + 2 * zstride;
                conn[ix+5] =  corner + 2               + 2 * zstride;
                conn[ix+6] =  corner + 2 + 2 * ystride + 2 * zstride;
                conn[ix+7] =  corner +     2 * ystride + 2 * zstride;
                conn[ix+8] =  corner + 1;                                           // 0-1
                conn[ix+9] =  corner + 2 +     ystride;                             // 1-2
                conn[ix+10] = corner + 1 + 2 * ystride;                             // 2-3
                conn[ix+11] = corner +         ystride;                             // 3-0
                conn[ix+12] = corner +                       zstride;               // 0-4
                conn[ix+13] = corner + 2 +                   zstride;               // 1-5
                conn[ix+14] = corner + 2 + 2 * ystride +     zstride;               // 2-6
                conn[ix+15] = corner +     2 * ystride +     zstride;               // 3-7
                conn[ix+16] = corner + 1 +               2 * zstride;               // 4-5
                conn[ix+17] = corner + 2 +     ystride + 2 * zstride;               // 5-6
                conn[ix+18] = corner + 1 + 2 * ystride + 2 * zstride;               // 6-7
                conn[ix+19] = corner +         ystride + 2 * zstride;               // 4-7
                conn[ix+20] = corner + 1 +                   zstride;               // 0154
                conn[ix+21] = corner + 2 +     ystride +     zstride;               // 1265
                conn[ix+22] = corner + 1 + 2 * ystride +     zstride;               // 2376
                conn[ix+23] = corner +         ystride +     zstride;               // 0374
                conn[ix+24] = corner + 1 +     ystride;                             // 0123
                conn[ix+25] = corner + 1 +     ystride + 2 * zstride;               // 4567
                conn[ix+26] = corner + 1 +     ystride +     zstride;               // center
                ix += 27;
              }
              else if (tetra) {
                //        E      H
                //     F     G
                //
                //        A     D
                //     B     C
                EntityHandle AA = corner;
                EntityHandle BB = corner + 1;
                EntityHandle CC = corner + 1 + ystride;
                EntityHandle D =  corner +     ystride;
                EntityHandle E =  corner +               zstride;
                EntityHandle F =  corner + 1 +           zstride;
                EntityHandle G =  corner + 1 + ystride + zstride;
                EntityHandle H =  corner +     ystride + zstride;

                // tet EDHG
                conn[ix]    = E;
                conn[ix+1]  = D;
                conn[ix+2]  = H;
                conn[ix+3]  = G;

                // tet ABCF
                conn[ix+4]  = AA;
                conn[ix+5]  = BB;
                conn[ix+6]  = CC;
                conn[ix+7]  = F;

                // tet ADEF
                conn[ix+8]  = AA;
                conn[ix+9]  = D;
                conn[ix+10] = E;
                conn[ix+11] = F;

                // tet CGDF
                conn[ix+12] = CC;
                conn[ix+13] = G;
                conn[ix+14] = D;
                conn[ix+15] = F;

                // tet ACDF
                conn[ix+16] = AA;
                conn[ix+17] = CC;
                conn[ix+18] = D;
                conn[ix+19] = F;

                // tet DGEF
                conn[ix+20] = D;
                conn[ix+21] = G;
                conn[ix+22] = E;
                conn[ix+23] = F;
                ix += 24;
                for (int ff = 0; ff < factor - 1; ff++) {
                  gids[ie] = gids[ie-1] + 1; // 6 more for tetra

                  eh = starte + ie;
                  for (size_t i = 0; i < doubleTags.size(); i++) {
                    double valv = gids[ie]/30. + i*5000.;
                    rval = mb->tag_set_data(doubleTags[i], &eh, 1, &valv);MB_CHK_SET_ERR(rval, "Can't set double tag on an element");
                  }
                  ie++;
                }
              }
              else { // Linear hex
                conn[ix] =   corner;
                conn[ix+1] = corner + 1;
                conn[ix+2] = corner + 1 + ystride;
                conn[ix+3] = corner +     ystride;
                conn[ix+4] = corner +               zstride;
                conn[ix+5] = corner + 1 +           zstride;
                conn[ix+6] = corner + 1 + ystride + zstride;
                conn[ix+7] = corner +     ystride + zstride;
                ix += 8;
              }
            }
          }
        }

        EntityHandle part_set;
        rval = mb->create_meshset(MESHSET_SET, part_set);MB_CHK_SET_ERR(rval, "Can't create mesh set");
        rval = mb->add_entities(part_set, cells);MB_CHK_SET_ERR(rval, "Can't add entities to set");
        // If needed, add all edges and faces
        if (adjEnts) {
          // We need to update adjacencies now, because some elements are new
          rval = iface->update_adjacencies(starte, num_el, num_v_per_elem, conn);MB_CHK_SET_ERR(rval, "Can't update adjacencies");
          // Generate all adj entities dimension 1 and 2 (edges and faces/ tri or qua)
          Range edges, faces;
          rval = mb->get_adjacencies(cells, 1, true, edges,
                                     Interface::UNION);MB_CHK_SET_ERR(rval, "Can't get edges");
          rval = mb->get_adjacencies(cells, 2, true, faces,
                                     Interface::UNION);MB_CHK_SET_ERR(rval, "Can't get faces");
          //rval = mb->add_entities(part_set, edges);MB_CHK_SET_ERR(rval, "Can't add edges to partition set");
          //rval = mb->add_entities(part_set, faces);MB_CHK_SET_ERR(rval, "Can't add faces to partition set");
        }

        rval = mb->tag_set_data(global_id_tag, cells, &gids[0]);MB_CHK_SET_ERR(rval, "Can't set global ids to elements");
        if (!parmerge){
          rval = mb->tag_set_data(new_id_tag, cells, &lgids[0]);MB_CHK_SET_ERR(rval, "Can't set new ids to elements");
        }
        int part_num = a + m*A + (b + n*B)*(M*A) + (c + k*C)*(M*A * N*B);
        rval = mb->tag_set_data(part_tag, &part_set, 1, &part_num);MB_CHK_SET_ERR(rval, "Can't set part tag on set");
        wsets.insert(part_set);
      }
    }
  }


  /*
  // Before merge locally
  rval = mb->write_file("test0.h5m", 0, ";;PARALLEL=WRITE_PART");MB_CHK_SET_ERR(rval, "Can't write in parallel, before merging");
  */
  // After the mesh is generated on each proc, merge the vertices
  MergeMesh mm(mb);
  Range all3dcells;
  rval = mb->get_entities_by_dimension(0, 3, all3dcells);MB_CHK_SET_ERR(rval, "Can't get all 3d cells elements");

  if (0 == rank) {
    cout << "generate local mesh: "
         << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
    tt = clock();
    cout << "number of elements on rank 0: " << all3dcells.size() << endl;
    cout << "Total number of elements " << all3dcells.size()*size << endl;
    cout << "Element type: " << ( tetra ? "MBTET" : "MBHEX") << " order:" << 
          (quadratic? "quadratic" : "linear" ) << endl;
  }
  Range verts;
  rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_SET_ERR(rval, "Can't get all vertices");

#ifdef MOAB_HAVE_MPI
  if (A*B*C != 1) { // Merge needed
    if (newMergeMethod) {
      rval = mm.merge_using_integer_tag(verts, global_id_tag);MB_CHK_SET_ERR(rval, "Can't merge");
    }
    else {
      rval = mm.merge_entities(all3dcells, 0.0001);MB_CHK_SET_ERR(rval, "Can't merge");
    }

    if (0 == rank) {
       cout << "merge locally: "
            << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
       tt = clock();
    }
  }
  // if adjEnts, add now to each set
  if (adjEnts)
  {
    for (Range::iterator wsit =wsets.begin(); wsit!=wsets.end(); ++wsit)
    {
      EntityHandle ws=*wsit;// write set
      Range cells,edges, faces;
      rval = mb->get_entities_by_dimension(ws, 3, cells);MB_CHK_SET_ERR(rval, "Can't get cells");
      rval = mb->get_adjacencies(cells, 1, false, edges,
                               Interface::UNION);MB_CHK_SET_ERR(rval, "Can't get edges");
      rval = mb->get_adjacencies(cells, 2, false, faces,
                               Interface::UNION);MB_CHK_SET_ERR(rval, "Can't get faces");
      rval = mb->add_entities(ws, edges);MB_CHK_SET_ERR(rval, "Can't add edges to partition set");
      rval = mb->add_entities(ws, faces);MB_CHK_SET_ERR(rval, "Can't add faces to partition set");
    }
  }
  if (size > 1) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(mb, 0);
    if (NULL == pcomm)
      pcomm = new ParallelComm(mb, MPI_COMM_WORLD);
    EntityHandle mesh_set;
    rval = mb->create_meshset(MESHSET_SET, mesh_set);MB_CHK_SET_ERR(rval, "Can't create new set");
    mb->add_entities(mesh_set, all3dcells);

    if (parmerge)
    {
      ParallelMergeMesh pm(pcomm, 0.00001);
      rval = pm.merge();MB_CHK_SET_ERR(rval, "Can't resolve shared ents");
      if (0 == rank) {
         cout << "parallel mesh merge: "
              << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
         tt = clock();
      }
    }
    else
    {
      rval = pcomm->resolve_shared_ents(mesh_set, -1, -1, &new_id_tag);MB_CHK_SET_ERR(rval, "Can't resolve shared ents");

      if (0 == rank) {
         cout << "resolve shared entities: "
              << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
         tt = clock();
      }
    }
    if (!keep_skins) { // Default is to delete the 1- and 2-dimensional entities
      // Delete all quads and edges
      Range toDelete;
      rval = mb->get_entities_by_dimension(0, 1, toDelete);MB_CHK_SET_ERR(rval, "Can't get edges");

      rval = mb->get_entities_by_dimension(0, 2, toDelete);MB_CHK_SET_ERR(rval, "Can't get faces");

      rval = pcomm->delete_entities(toDelete);MB_CHK_SET_ERR(rval, "Can't delete entities");

      if (0 == rank) {
        cout << "delete edges and faces, and correct sharedEnts: "
             << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
        tt = clock();
      }
    }
    // do some ghosting if required
    if (GL>0)
    {
      rval = pcomm->exchange_ghost_cells(3, // int ghost_dim
                                         0, // int bridge_dim
                                         GL, // int num_layers
                                         0, // int addl_ents
                                         true);MB_CHK_ERR(rval); // bool store_remote_handles
      if (0 == rank) {
         cout << "exchange  " << GL << " ghost layer(s) :"
              << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
         tt = clock();
      }
    }
  }
#endif

#ifdef MOAB_HAVE_HDF5_PARALLEL
  if (!parmerge)
  {
    rval = mb->tag_delete(new_id_tag); MB_CHK_SET_ERR(rval, "Can't delete new ID tag");
  }
  if (!nosave){
    rval = mb->write_file(outFileName.c_str(), 0, ";;PARALLEL=WRITE_PART;CPUTIME;", wsets);MB_CHK_SET_ERR(rval, "Can't write in parallel");
    if (0 == rank) {
      cout << "write file " << outFileName << " in "
           << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
      tt = clock();
    }
  }
  // delete the mesh that we already have in-memory
  size_t nLocalVerts = verts.size();
  size_t nLocalCells = all3dcells.size();
#else
  rval = mb->write_file("GenLargeMesh.vtk", 0, "", wsets);MB_CHK_SET_ERR(rval, "Can't write in serial");
#endif
  mb->delete_mesh();

#ifdef MOAB_HAVE_HDF5_PARALLEL
  if (!nosave && readb)
  {
    // now recreate a core instance and load the file we just wrote out to verify
    Core mb2;
    std::string read_opts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;CPUTIME;");
    if (readAndGhost)
      read_opts+="PARALLEL_GHOSTS=3.0.1;";
    rval = mb2.load_file(outFileName.c_str(), 0, read_opts.c_str());MB_CHK_SET_ERR(rval, "Can't read in parallel");
    if (0 == rank) {
      cout << "read back file " << outFileName << " with options: \n" << read_opts <<
          " in "  << (clock() - tt) / (double)CLOCKS_PER_SEC << " seconds" << endl;
      tt = clock();
    }
    moab::Range nverts, ncells;
    rval = mb2.get_entities_by_dimension(0, 0, nverts);MB_CHK_SET_ERR(rval, "Can't get all vertices");
    rval = mb2.get_entities_by_dimension(0, 3, ncells);MB_CHK_SET_ERR(rval, "Can't get all 3d cells elements");

    if (readAndGhost && size > 1)
    {
      // filter out the ghost nodes and elements, for comparison with original mesh
      // first get the parallel comm
      ParallelComm* pcomm2 = ParallelComm::get_pcomm(&mb2, 0);
      if (NULL == pcomm2) MB_SET_ERR(MB_FAILURE, "can't get parallel comm.");
      rval = pcomm2->filter_pstatus(nverts, PSTATUS_GHOST, PSTATUS_NOT);MB_CHK_SET_ERR(rval, "Can't filter ghost vertices");
      rval = pcomm2->filter_pstatus(ncells, PSTATUS_GHOST, PSTATUS_NOT);MB_CHK_SET_ERR(rval, "Can't filter ghost cells");
    }
    if (nverts.size() != nLocalVerts && ncells.size() != nLocalCells ) {
      MB_SET_ERR(MB_FAILURE, "Reading back the output file led to inconsistent number of entities.");
    }

    // delete the mesh that we already have in-memory
    mb2.delete_mesh();
  }
#endif 

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
