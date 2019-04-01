/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifndef MOAB_READ_UTIL_IFACE_HPP
#define MOAB_READ_UTIL_IFACE_HPP

#include <vector>
#include <string>
#include "moab/Types.hpp"
#include "moab/Compiler.hpp"

namespace moab {

class Range;

//! Interface implemented in MOAB which provides memory for mesh reading utilities
class ReadUtilIface
{
public:

  //! Constructor
  ReadUtilIface(){}

  //! Destructor
  virtual ~ReadUtilIface(){}

  //! Given a requested number of vertices and number of coordinates, returns
  //! memory space which will be used to store vertex coordinates and information
  //! about what handles those new vertices are assigned; allows direct read of
  //! coordinate data into memory
  //! \param num_arrays Number of node position arrays requested
  //! \param num_nodes Number of nodes
  //! \param preferred_start_id Preferred integer id starting value
  //! \param actual_start_handle Actual starting id value
  //! \param arrays STL vector of double*'s, point to memory storage to be used for
  //!     these vertices
  //! \param sequence_size If specified, allocate this sequence size instead of
  //!      SequenceManager::DEFAULT_VERTEX_SEQUENCE_SIZE
  //! \return status Success/failure of this call
  virtual ErrorCode get_node_coords(const int num_arrays,
                                    const int num_nodes,
                                    const int preferred_start_id,
                                    EntityHandle& actual_start_handle,
                                    std::vector<double*>& arrays,
                                    const int sequence_size = -1) = 0;

  //! Given requested number of elements, element type, and number of
  //! elements, returns pointer to memory space allocated to store connectivity
  //! of those elements; allows direct read of connectivity data into memory
  //! \param num_elements Number of elements being requested
  //! \param verts_per_element Number of vertices per element (incl. higher-order nodes)
  //! \param mdb_type Element type
  //! \param preferred_start_id Preferred integer id for first element
  //! \param actual_start_handle Actual integer id for first element (returned)
  //! \param array Pointer to memory allocated for storing connectivity for these elements
  //! \param sequence_size If specified, allocate this sequence size instead of
  //!      SequenceManager::DEFAULT_VERTEX_SEQUENCE_SIZE
  //! \return status Success/failure of this call
  virtual ErrorCode get_element_connect(const int num_elements,
                                        const int verts_per_element,
                                        const EntityType mdb_type,
                                        const int preferred_start_id,
                                        EntityHandle& actual_start_handle,
                                        EntityHandle*& array,
                                        int sequence_size = -1) = 0;

  /**
   *\brief Gather entities related to those in the partition
   * Gather entities related to those in the input partition. Related
   * means down-adjacent to, contained in, etc.
   * \param partition Entities for which to gather related entities
   * \param related_ents Related entities
   * \param file_set If non-NULL, entity sets contained in this set will be checked;
   *        otherwise, all sets in the instance will be checked
   */
  virtual ErrorCode gather_related_ents(Range &partition,
                                        Range &related_ents,
                                        EntityHandle *file_set = NULL) = 0;

  virtual ErrorCode create_entity_sets(EntityID num_sets,
                                       const unsigned* set_flags,
                                       EntityID preffered_start_id,
                                       EntityHandle& actual_start_handle) = 0;

  //! Update adjacencies
  //! Given information about new elements, adjacency information will be updated
  //! in MOAB. Think of this function as a way of Readers telling MOAB what elements are
  //! new because we aren't using the Interface to create elements.
  //! \param start_handle Handle of first new element
  //! \param number_elements Number of new elements
  //! \param number_vertices_per_element Number of vertices in each new element
  //! \param conn_array Connectivity of new elements
  //! \return status Success/failure of this call
  virtual ErrorCode update_adjacencies(const EntityHandle start_handle,
                                       const int number_elements,
                                       const int number_vertices_per_element,
                                       const EntityHandle* conn_array) = 0;

  /**\brief Re-order incoming element connectivity
   *
   * Permute the connectivity of each element such that the node
   * order is that of MBCN rather than the target file format.
   *\param order The permutation to use.  Must be an array of 'node_per_elem'
   *             integers and be a permutation of the values [0..node_per_elem-1].
   *             Such that for a single element:
   *             mbcn_conn[order[i]] == target_conn[i]
   *\param conn  The connectivity array to re-order
   *\param num_elem  The number of elements in the connectivity array
   *\param node_per_elem The number of nodes in each element's connectivity list.
   */
  static inline void reorder(const int* order, EntityHandle* conn,
                             int num_elem, int node_per_elem);

  //! Given an ordered list of bounding entities and the sense of
  //! those entities, return an ordered list of vertices
  virtual ErrorCode get_ordered_vertices(EntityHandle *bound_ents,
                                         int *sense,
                                         int num_bound,
                                         int dim,
                                         EntityHandle *bound_verts,
                                         EntityType &etype) = 0;

  //! Assign sequential IDS to entities in range and store IDs in tag
  virtual ErrorCode assign_ids(Tag id_tag, const Range& ents,
                               int start = 0 ) = 0;

  //! Assign to each entity in an array the ID that is its position
  //! in the array plus the value of 'start'.  For any non-zero handles
  //! in the array, store the ID value in the passed tag.
  virtual ErrorCode assign_ids(Tag id_tag, const EntityHandle* ents,
                               size_t num_ents, int start = 0) = 0;

  //! Create a new gather set with tag GATHER_SET
  virtual ErrorCode create_gather_set(EntityHandle& gather_set) = 0;

  //! Get entity handle of an existing gather set
  virtual ErrorCode get_gather_set(EntityHandle& gather_set) = 0;
};

inline void ReadUtilIface::reorder(const int* order, EntityHandle* conn,
                                   int num_elem, int node_per_elem)
{
  std::vector<EntityHandle> elem(node_per_elem);
  EntityHandle* const end = conn + num_elem*node_per_elem;
  while (conn != end) {
    std::copy(conn, conn + node_per_elem, elem.begin());
    for (int j = 0; j < node_per_elem; ++j)
      conn[order[j]] = elem[j];
    conn += node_per_elem;
  }
}

} // namespace moab 

#endif
