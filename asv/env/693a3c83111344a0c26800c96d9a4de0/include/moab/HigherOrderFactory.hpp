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


#ifndef MOAB_HIGHER_ORDER_FACTORY_HPP
#define MOAB_HIGHER_ORDER_FACTORY_HPP

#ifndef IS_BUILDING_MB
#error "HigherOrderFactory.hpp isn't supposed to be included into an application"
#endif

#include "moab/Interface.hpp"

namespace moab {

class ElementSequence;
class Core;

/** \class HigherOrderFactory
 * \brief Functions for converting to/from higher-order elements
 *  \authors Clinton Stimpson
 *  \date    11/25/02
 *  \brief   
 *          
 */ 
class HigherOrderFactory
{
public:

  HigherOrderFactory(Core* , Interface::HONodeAddedRemoved* function_object);
  ~HigherOrderFactory();

  ErrorCode convert(const EntityHandle meshset, const bool mid_edge_nodes, 
                       const bool mid_face_nodes, const bool mid_volume_nodes);

  ErrorCode convert( const Range& entities, const bool mid_edge_nodes, 
                       const bool mid_face_nodes, const bool mid_volume_nodes);

  unsigned char mNodeMap[MBMAXTYPE][8][8];

private:

  //static bool mMapInitialized;
  void initialize_map();

  Core* mMB;
  Interface::HONodeAddedRemoved* mHONodeAddedRemoved;

  ErrorCode convert_sequence(ElementSequence* sequence, 
                               EntityHandle sequence_subset_start,
                               EntityHandle sequence_subset_end,
                               bool mid_edge_nodes, 
                               bool mid_face_nodes,
                               bool mid_volume_nodes);
  ErrorCode add_mid_edge_nodes(ElementSequence*);
  ErrorCode add_mid_face_nodes(ElementSequence*);
  ErrorCode add_mid_volume_nodes(ElementSequence*);

  //! returns the handle of the first center node found between the two corner nodes.
  //! returns zero if none found
  //! entities that share those two corner nodes and have space allocated for mid-edge nodes are returned in a vector
  EntityHandle center_node_exist(EntityHandle corner1, EntityHandle corner2,
         std::vector<EntityHandle>& adj_entities);
  
  //! returns the handle of the first center node found between the 3-4 corner nodes.
  //! set the last node to zero if you want only 3 nodes
  //! returns zero if none found
  //! entities that share those corner nodes and have space allocated for mid face nodes are returned in a vector
  EntityHandle center_node_exist(EntityHandle corners[4], std::vector<EntityHandle>& adj_entities);

  //! adds a center node to element between corner nodes, returns success
  bool add_center_node(EntityType type, EntityHandle* element_conn, int conn_size, 
      EntityHandle corner_node1, EntityHandle corner_node2, EntityHandle center_node);


  ErrorCode copy_corner_nodes( ElementSequence* src, ElementSequence* dst );
  ErrorCode copy_mid_edge_nodes( ElementSequence* src, ElementSequence* dst ); 
  ErrorCode copy_mid_face_nodes( ElementSequence* src, ElementSequence* dst ); 
  ErrorCode copy_mid_volume_nodes( ElementSequence* src, ElementSequence* dst ); 
  ErrorCode copy_nodes( ElementSequence* src, 
                          ElementSequence* dst,
                          unsigned nodes_per_elem_to_copy,
                          unsigned src_conn_offset,
                          unsigned dst_conn_offset ); 

  ErrorCode zero_mid_edge_nodes( ElementSequence* dst ); 
  ErrorCode zero_mid_face_nodes(  ElementSequence* dst ); 
  ErrorCode zero_mid_volume_nodes( ElementSequence* dst ); 
  ErrorCode zero_nodes( ElementSequence* dst,
                        unsigned nodes_per_elem_to_zero,
                        unsigned dst_conn_offset ); 

  ErrorCode remove_mid_edge_nodes( ElementSequence* seq, 
                                     EntityHandle start,
                                     EntityHandle stop,
                                     Tag deletable_ndoes );
  ErrorCode remove_mid_face_nodes( ElementSequence* seq, 
                                     EntityHandle start,
                                     EntityHandle stop,
                                     Tag deletable_ndoes );
  ErrorCode remove_mid_volume_nodes( ElementSequence* seq, 
                                       EntityHandle start,
                                       EntityHandle stop,
                                       Tag deletable_ndoes );
  ErrorCode remove_ho_nodes( ElementSequence* sequence,
                               EntityHandle subset_start_handle,
                               EntityHandle subset_end_handle,
                               int nodes_per_elem_to_remove,
                               int elem_conn_offset_to_remove,
                               Tag deletable_nodes );
  bool tag_for_deletion( EntityHandle element_with_node,
                         int node_index_in_elem_connectivity,
                         ElementSequence* sequence );
};

} // namespace moab 

#endif

