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


#ifndef MOAB_MESH_TOPO_UTIL_HPP
#define MOAB_MESH_TOPO_UTIL_HPP

#include "moab/Forward.hpp"

namespace moab {

/*!
 *  \authors Tim Tautges
 *  \date    2/04
 *  \brief   MeshTopoUtil contains general mesh utility functions 
 *          
 */ 
class MeshTopoUtil
{
public:
  MeshTopoUtil(Interface *impl) : mbImpl(impl) {}
  
  ~MeshTopoUtil() {}

    //! generate all the AEntities bounding the vertices
  ErrorCode construct_aentities(const Range &vertices);

    //! given an entity, get its average position (avg vertex locations)
  ErrorCode get_average_position(Range &entities,
                                   double *avg_position);

    //! given an entity, get its average position (avg vertex locations)
  ErrorCode get_average_position(const EntityHandle entity,
                                   double *avg_position);

    //! given a set of entities, get their average position (avg vertex locations)
  ErrorCode get_average_position(const EntityHandle *entities,
                                   const int num_entities,
                                   double *avg_position);

    //! get (target_dim)-dimensional manifold entities connected to star_entity; that is,
    //! the entities with <= 1 connected (target_dim+2)-dimensional adjacent entities;
    //! for target_dim=3, just return all of them
    //! just insert into the list, w/o clearing manifold list first
  ErrorCode get_manifold(const EntityHandle star_entity,
                           const int target_dim,
                           Range &manifold);
  
  //! given an entity, find the entities of next higher dimension around
  //! that entity, ordered by connection through next higher dimension entities; 
  //! if any of the star entities is in only entity of next higher dimension, 
  //! on_boundary is returned true
  ErrorCode star_entities(const EntityHandle star_center,
                            std::vector<EntityHandle> &star_entities,
                            bool &bdy_entity,
                            const EntityHandle starting_star_entity = 0,
                            std::vector<EntityHandle> *star_entities_dp1 = NULL,
                            Range *star_entities_candidates_dp1 = NULL);

    //! Get a series of (d+1)-dimensional stars around a d-dimensional entity, such that
    //! each star is on a (d+2)-manifold containing the d-dimensional entity; each star
    //! is either open or closed, and also defines a (d+2)-star whose entities are bounded by
    //! (d+1)-entities on the star and on the (d+2)-manifold
  ErrorCode star_entities_nonmanifold(const EntityHandle star_entity,
                                        std::vector<std::vector<EntityHandle> > &stars,
                                        std::vector<bool> *bdy_flags = NULL,
                                        std::vector<std::vector<EntityHandle> > *dp2_stars = NULL);

    //! given a star_center, a last_entity (whose dimension should be 1 greater than center)
    //! and last_dp1 (dimension 2 higher than center), returns the next star entity across
    //! last_dp1, and the next dp1 entity sharing next_entity; if star_candidates is non-empty,
    //! star must come from those
  ErrorCode star_next_entity(const EntityHandle star_center,
                               const EntityHandle last_entity,
                               const EntityHandle last_dp1,
                               Range *star_candidates_dp1,
                               EntityHandle &next_entity,
                               EntityHandle &next_dp1);
  
    //! get "bridge" or "2nd order" adjacencies, going through dimension bridge_dim
  ErrorCode get_bridge_adjacencies(Range &from_entities,
                                     int bridge_dim,
                                     int to_dim, 
                                     Range &to_ents,
                                     int num_layers = 1);
  
    //! get "bridge" or "2nd order" adjacencies, going through dimension bridge_dim
  ErrorCode get_bridge_adjacencies(const EntityHandle from_entity,
                                     const int bridge_dim,
                                     const int to_dim,
                                     Range &to_adjs);

    //! return a common entity of the specified dimension, or 0 if there isn't one
  EntityHandle common_entity(const EntityHandle ent1,
                               const EntityHandle ent2,
                               const int dim);
  
  //! return the opposite side entity given a parent and bounding entity.
  //! This function is only defined for certain types of parent/child types;
  //! See MBCN.hpp::OppositeSide for details.
  //!
  //! \param parent The parent element
  //! \param child The child element
  //! \param opposite_element The index of the opposite element
  ErrorCode opposite_entity(const EntityHandle parent,
                              const EntityHandle child,
                              EntityHandle &opposite_element);

    //! split entity which is non-manifold, that is, which has > 2 connected entities
    //! of next higher dimension; assumes that there are >= 2 connected regions of
    //! (d+2)-dimensional entities; a new d-entity is created for each region after the
    //! first, and it's made explicitly-adjacent to the region to which it corresponds
  ErrorCode split_entity_nonmanifold(EntityHandle split_ent,
                                       Range &old_adjs,
                                       Range &new_adjs,
                                       EntityHandle &new_entity);
  
    //! split entities that are manifold (shared by two or less entities of each higher dimension),
    //! optionally creating an entity of next higher dimension to fill the gap
    /**
       \param entities The entities to be split
       \param new_entities New entities, in order of correspondence to that of entities
       \param fill_entities If non-NULL, create an entity of next higher dimension to fill the gap,
                       passing it back in *fill_entities
    */
  ErrorCode split_entities_manifold(Range &entities,
                                      Range &new_entities,
                                      Range *fill_entities);
  
    //! split entities that are manifold (shared by two or less entities of each higher dimension),
    //! optionally creating an entity of next higher dimension to fill the gap
    /**
       \param entities The entities to be split
       \param new_entities New entities, in order of correspondence to that of entities
       \param fill_entities If non-NULL, create an entity of next higher dimension to fill the gap,
                       passing it back in *fill_entities
       \param gowith_ents If non-NULL, each of the new entities will adj to the
                       corresponding gowith entities after the split; this parameter is
                       ignored for boundary split entities; in that case, the split entity
                       remains on the boundary (i.e. not adj to any entity of higher 
                       dimension).  Dimension of gowith_ents must be the same as entities.
    */
  ErrorCode split_entities_manifold(EntityHandle *entities,
                                      const int num_entities,
                                      EntityHandle *new_entities,
                                      Range *fill_entities,
                                      EntityHandle *gowith_ents = NULL);

    //! return whether entity is equivalent to any other of same type and same vertices;
    //! if equivalent entity is found, it's returned in equiv_ents and return value is true,
    //! false otherwise.
  bool equivalent_entities(const EntityHandle entity,
                           Range *equiv_ents = NULL);
  
                                  
private:
  Interface *mbImpl;
  
};

} // namespace moab 

#endif

