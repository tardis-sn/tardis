/** @example DirectAccessWithHoles.cpp \n
 * \brief Use direct access to MOAB data to avoid calling through API \n
 *
 * This example creates a 1d row of quad elements, with a user-specified number of "holes" (missing quads) in the row:
 * \verbatim
 *  ----------------------      ----------------------      --------
 *  |      |      |      |      |      |      |      |      |      |       
 *  |      |      |      |(hole)|      |      |      |(hole)|      | ...
 *  |      |      |      |      |      |      |      |      |      |
 *  ----------------------      ----------------------      --------
 * \endverbatim
 * This makes (nholes+1) contiguous runs of quad handles in the handle space
 * This example shows how to use the xxx_iterate functions in MOAB (xxx = coords, connect, tag, adjacencies) to get 
 * direct pointer access to MOAB internal storage, which can be used without calling through the MOAB API.
 *
 *    -#  Initialize MOAB \n
 *    -#  Create a quad mesh with holes, as depicted above
 *    -#  Create 2 dense tags (tag1, tag2) for avg position to assign to quads, and # verts per quad (tag3)
 *    -#  Get connectivity, coordinate, tag1 iterators
 *    -#  Iterate through quads, computing midpoint based on vertex positions, set on quad-based tag1
 *    -#  Set up map from starting quad handle for a chunk to struct of (tag1_ptr, tag2_ptr, tag3_ptr), pointers to
 *        the dense tag storage for those tags for the chunk
 *    -#  Iterate through vertices, summing positions into tag2 on connected quads and incrementing vertex count
 *    -#  Iterate through quads, normalizing tag2 by vertex count and comparing values of tag1 and tag2
 *
 * <b>To compile</b>: \n
 *    make DirectAccessWithHoles MOAB_DIR=<installdir> \n
 * <b>To run</b>: ./DirectAccess [-nquads <# quads>] [-holes <# holes>]\n
 *
 */

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ReadUtilIface.hpp"
#include <map>
#include <iostream>
#include <assert.h>

using namespace moab;
using namespace std;

ErrorCode create_mesh_with_holes(Interface *mbImpl, int nquads, int nholes);

struct tag_struct {double *avg_ptr; int *nv_ptr;};

int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mbImpl = new (std::nothrow) Core;
  if (NULL == mbImpl)
    return 1;

  int nquads = 1000, nholes = 1;

  // Parse options
  ProgOptions opts;
  opts.addOpt<int>(string("nquads,n"), string("Number of quads in the mesh (default = 1000"), &nquads);
  opts.addOpt<int>(string("holes,H"), string("Number of holes in the element handle space (default = 1"), &nholes);
  opts.parseCommandLine(argc, argv);
  if (nholes >= nquads) {
    cerr << "Number of holes needs to be < number of elements." << endl;
    return 1;
  }

  // Create simple structured mesh with hole, but using unstructured representation
  ErrorCode rval = create_mesh_with_holes(mbImpl, nquads, nholes);MB_CHK_SET_ERR(rval, "Trouble creating mesh");

  // Get all vertices and non-vertex entities
  Range verts, elems;
  rval = mbImpl->get_entities_by_handle(0, elems);MB_CHK_SET_ERR(rval, "Trouble getting all entities");
  verts = elems.subset_by_dimension(0);
  elems -= verts;

  // Create tag1 (element-based avg), tag2 (vertex-based avg), tag3 (# connected verts)
  Tag tag1, tag2, tag3;
  rval = mbImpl->tag_get_handle("tag1", 3, MB_TYPE_DOUBLE, tag1, MB_TAG_DENSE | MB_TAG_CREAT);MB_CHK_SET_ERR(rval, "Trouble creating tag1");
  double def_val[3] = {0.0, 0.0, 0.0}; // Need a default value for tag2 because we sum into it
  rval = mbImpl->tag_get_handle("tag2", 3, MB_TYPE_DOUBLE, tag2, MB_TAG_DENSE | MB_TAG_CREAT, def_val);MB_CHK_SET_ERR(rval, "Trouble creating tag2");
  int def_val_int = 0;  // Need a default value for tag3 because we increment it
  rval = mbImpl->tag_get_handle("tag3", 1, MB_TYPE_INTEGER, tag3, MB_TAG_DENSE | MB_TAG_CREAT, &def_val_int);MB_CHK_SET_ERR(rval, "Trouble creating tag3");

  // Get connectivity, coordinate, tag, and adjacency iterators
  EntityHandle *conn_ptr;
  double *x_ptr, *y_ptr, *z_ptr, *tag1_ptr, *tag2_ptr;
  int *tag3_ptr;

  // First vertex is at start of range (ranges are sorted), and is offset for vertex index calculation
  EntityHandle first_vert = *verts.begin();

  // When iterating over elements, each chunk can have a different # vertices; also, count tells you how many
  // elements are in the current chunk
  int vpere, count;

  // Get coordinates iterator, just need this once because we know verts handle space doesn't have holes
  rval = mbImpl->coords_iterate(verts.begin(), verts.end(), x_ptr, y_ptr, z_ptr, count);MB_CHK_SET_ERR(rval, "Error in coords_iterate");
  assert(count == (int) verts.size()); // Should end up with just one contiguous chunk of vertices

  // Iterate through elements, computing midpoint based on vertex positions, set on element-based tag1
  // Control loop by iterator over elem range
  Range::iterator e_it = elems.begin();

  while (e_it != elems.end()) {
    // Get conn_ptr, tag1_ptr for next contiguous chunk of element handles, and coords pointers for all verts
    rval = mbImpl->connect_iterate(e_it, elems.end(), conn_ptr, vpere, count);MB_CHK_SET_ERR(rval, "Error in connect_iterate");
    rval = mbImpl->tag_iterate(tag1, e_it, elems.end(), count, (void*&)tag1_ptr);MB_CHK_SET_ERR(rval, "Error in tag1_iterate");

    // Iterate over elements in this chunk
    for (int i = 0; i < count; i++) {
      tag1_ptr[0] = tag1_ptr[1] = tag1_ptr[2] = 0.0; // Initialize position for this element
      for (int j = 0; j < vpere; j++) { // Loop over vertices in this element
        int v_index = conn_ptr[j] - first_vert; // vert index is just the offset from first vertex
        tag1_ptr[0] += x_ptr[v_index];
        tag1_ptr[1] += y_ptr[v_index]; // Sum vertex positions into tag1...
        tag1_ptr[2] += z_ptr[v_index];
      }
      tag1_ptr[0] /= vpere;
      tag1_ptr[1] /= vpere; // Then normalize
      tag1_ptr[2] /= vpere;

      // Done with this element; advance connect_ptr and tag1_ptr to next element
      conn_ptr += vpere;
      tag1_ptr += 3;
    } // Loop over elements in chunk

    // Done with chunk; advance range iterator by count; will skip over gaps in range
    e_it += count;
  } // While loop over all elements

  // Iterate through vertices, summing positions into tag2 on connected elements and incrementing vertex count
  // Iterate over chunks the same as elements, even though we know we have only one chunk here, just to show
  // how it's done

  // Create a std::map from EntityHandle (first entity handle in chunk) to
  // tag_struct (ptrs to start of avg/#verts tags for that chunk); then for a given entity handle, we can quickly
  // find the chunk it's in using map::lower_bound; could have set up this map in earlier loop over elements, but do
  // it here for clarity

  map<EntityHandle, tag_struct> elem_map;
  e_it = elems.begin();
  while (e_it != elems.end()) {
    tag_struct ts = {NULL, NULL};
    rval = mbImpl->tag_iterate(tag2, e_it, elems.end(), count, (void*&)ts.avg_ptr);MB_CHK_SET_ERR(rval, "Error in tag2_iterate");
    rval = mbImpl->tag_iterate(tag3, e_it, elems.end(), count, (void*&)ts.nv_ptr);MB_CHK_SET_ERR(rval, "Error in tag3_iterate");
    elem_map[*e_it] = ts;
    e_it += count;
  }

  // Call a vertex-quad adjacencies function to generate vertex-element adjacencies in MOAB
  Range::iterator v_it = verts.begin();
  Range dum_range;
  rval = mbImpl->get_adjacencies(&(*v_it), 1, 2, false, dum_range);MB_CHK_SET_ERR(rval, "Error in get_adjacencies");
  const vector<EntityHandle> **adjs_ptr;
  while (v_it != verts.end()) {
    // Get coords ptrs, adjs_ptr; can't set tag2_ptr by direct access, because of hole in element handle space
    rval = mbImpl->coords_iterate(v_it, verts.end(), x_ptr, y_ptr, z_ptr, count);MB_CHK_SET_ERR(rval, "Error in coords_iterate");
    rval = mbImpl->adjacencies_iterate(v_it, verts.end(), adjs_ptr, count);MB_CHK_SET_ERR(rval, "Error in adjacencies_iterate");

    for (int v = 0; v < count; v++) {
      const vector<EntityHandle> *avec = *(adjs_ptr + v);
      for (vector<EntityHandle>::const_iterator ait = avec->begin(); ait != avec->end(); ++ait) {
        // Get chunk that this element resides in; upper_bound points to the first element strictly > key, so get that and decrement
        // (would work the same as lower_bound with an if-test and decrement)
        map<EntityHandle, tag_struct>::iterator mit = elem_map.upper_bound(*ait); --mit;
        // Index of *ait in that chunk
        int a_ind = *ait - (*mit).first;
        tag_struct ts = (*mit).second;
        ts.avg_ptr[3*a_ind + 0] += x_ptr[v]; // Tag on each element is 3 doubles, x/y/z
        ts.avg_ptr[3*a_ind + 1] += y_ptr[v];
        ts.avg_ptr[3*a_ind + 2] += z_ptr[v];
        ts.nv_ptr[a_ind]++; // Increment the vertex count
      }
    }

    v_it += count;
  }

  // Normalize tag2 by vertex count; loop over elements using same approach as before
  // At the same time, compare values of tag1 and tag2
  e_it = elems.begin();
  while (e_it != elems.end()) {
    // Get conn_ptr, tag1_ptr for next contiguous chunk of element handles, and coords pointers for all verts
    rval = mbImpl->tag_iterate(tag1, e_it, elems.end(), count, (void*&)tag1_ptr);MB_CHK_SET_ERR(rval, "Error in tag1_iterate");
    rval = mbImpl->tag_iterate(tag2, e_it, elems.end(), count, (void*&)tag2_ptr);MB_CHK_SET_ERR(rval, "Error in tag2_iterate");
    rval = mbImpl->tag_iterate(tag3, e_it, elems.end(), count, (void*&)tag3_ptr);MB_CHK_SET_ERR(rval, "Error in tag3_iterate");

    // Iterate over elements in this chunk
    for (int i = 0; i < count; i++) {
      for (int j = 0; j < 3; j++)
        tag2_ptr[3*i + j] /= (double)tag3_ptr[i]; // Normalize by # verts
      if (tag1_ptr[3*i] != tag2_ptr[3*i] || tag1_ptr[3*i + 1] != tag2_ptr[3*i + 1] || tag1_ptr[3*i + 2] != tag2_ptr[3*i + 2])
        cout << "Tag1, tag2 disagree for element " << *e_it + i << endl;
    }

    e_it += count;
  }

  // Ok, we're done, shut down MOAB
  delete mbImpl;

  return 0;
}

ErrorCode create_mesh_with_holes(Interface *mbImpl, int nquads, int nholes) 
{
  // First make the mesh, a 1d array of quads with left hand side x = elem_num; vertices are numbered in layers
  ReadUtilIface *read_iface;
  ErrorCode rval = mbImpl->query_interface(read_iface);MB_CHK_SET_ERR(rval, "Error in query_interface");
  vector<double*> coords;
  EntityHandle start_vert, start_elem, *connect;
  // Create verts, num is 4(nquads+1) because they're in a 1d row; will initialize coords in loop over elems later
  rval = read_iface->get_node_coords(3, 2*(nquads + 1), 0, start_vert, coords);MB_CHK_SET_ERR(rval, "Error in get_node_arrays");
  // Create elems
  rval = read_iface->get_element_connect(nquads, 4, MBQUAD, 0, start_elem, connect);MB_CHK_SET_ERR(rval, "Error in get_element_connect");
  for (int i = 0; i < nquads; i++) {
    coords[0][2*i] = coords[0][2*i + 1] = (double) i; // x values are all i
    coords[1][2*i] = 0.0; coords[1][2*i + 1] = 1.0; // y coords
    coords[2][2*i] = coords[2][2*i + 1] = (double) 0.0; // z values, all zero (2d mesh)
    EntityHandle quad_v = start_vert + 2*i;
    for (int j = 0; j < 4; j++) connect[4*i + j] = quad_v + j; // Connectivity of each quad is a sequence starting from quad_v
  }
  // Last two vertices
  // Cppcheck warning (false positive): variable coords is assigned a value that is never used
  coords[0][2*nquads] = coords[0][2*nquads + 1] = (double) nquads;
  coords[1][2*nquads] = 0.0; coords[1][2*nquads + 1] = 1.0; // y coords
  coords[2][2*nquads] = coords[2][2*nquads + 1] = (double) 0.0; // z values, all zero (2d mesh)

  // Now delete nholes elements, spaced approximately equally through mesh, so contiguous size is about nquads/(nholes + 1)
  // reinterpret start_elem as the next element to be deleted
  int de = nquads / (nholes + 1);
  for (int i = 0; i < nholes; i++) {
    start_elem += de;
    rval = mbImpl->delete_entities(&start_elem, 1);MB_CHK_SET_ERR(rval, "Error in delete_entities");
  }

  return MB_SUCCESS;
}
