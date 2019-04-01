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



#ifndef MOAB_GEOM_TOPO_TOOL_HPP
#define MOAB_GEOM_TOPO_TOOL_HPP

#include "moab/Forward.hpp"
#include "moab/Range.hpp"
#include "moab/OrientedBoxTreeTool.hpp"

#include <map>

namespace moab {

/** \class GeomTopoTool
 * \brief Tool for interpreting geometric topology sets in MOAB database
 * Tool for interpreting geometric topology sets in MOAB database; see MOAB metadata_info
 * document for information on how geometric topology sets are read and represented.
 */
class GeomTopoTool
{
public:
  GeomTopoTool(Interface *impl, bool find_geoments = false, EntityHandle modelRootSet = 0);
  ~GeomTopoTool() {}
  
    //! Restore parent/child links between GEOM_TOPO mesh sets
  ErrorCode restore_topology();

  //! Store sense of entity relative to wrt_entity.
     //!\return MB_MULTIPLE_ENTITIES_FOUND if surface already has a forward volume.
     //!        MB_SUCCESS if successful
     //!        otherwise whatever internal error code occured.
   ErrorCode set_sense( EntityHandle entity,
                        EntityHandle wrt_entity,
                        int sense);

     //! Get the sense of entity with respect to wrt_entity
     //! Returns MB_ENTITY_NOT_FOUND if no relationship found
   ErrorCode get_sense( EntityHandle entity,
                        EntityHandle wrt_entity,
                        int & sense );

  ErrorCode get_senses (EntityHandle entity,
    std::vector<EntityHandle> &wrt_entities,
    std::vector<int> &senses);

  ErrorCode set_senses (EntityHandle entity,
                          std::vector<EntityHandle> &wrt_entities,
                          std::vector<int> &senses);

    /** \brief get the other (d-1)-dimensional entity bounding a set across a (d-2)-dimensional entity
     *
     * Given a d-dimensional entity and one (d-1)-dimensional entity, return the (d-1) dimensional
     * entity across a specified (d-2)-dimensional entity.  For example, given a surface, edge, and vertex,
     * returns the other edge bounding the surface sharing the vertex.  In the case of degenerate results,
     * e.g. two loops bounding a surface and sharing a vertex, tries to step in positively-oriented
     * direction.  This won't always work; in those cases, will return MB_MULTIPLE_ENTITIES_FOUND.
     *
     * In the special case where bounded is a curve, then not_this can be a vertex and across zero.
     * This function returns the other vertex on the curve.
     */
  ErrorCode other_entity(EntityHandle bounded, EntityHandle not_this, EntityHandle across,
                         EntityHandle &other);

    /** \brief return the dimension of the set, or -1 if it's not a geom_dimension set
     */
  int dimension(EntityHandle this_set);
  
  // used mostly for debugging purposes
  int global_id(EntityHandle this_set);
  
  ErrorCode find_geomsets(Range *ranges = NULL);

  ErrorCode construct_obb_trees(bool make_one_vol = false);

  ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle &root);

  EntityHandle get_one_vol_root();

  OrientedBoxTreeTool *obb_tree() {return &obbTree;}

  // this could make the obb tree out of date
  ErrorCode add_geo_set(EntityHandle set, int dimension, int global_id  = 0);

  // will assume no geo sets are defined for this surface
  // will output a mesh_set that contains everything (all sets of interest), for proper output
  ErrorCode geometrize_surface_set(EntityHandle surface, EntityHandle & output);

  // this would be a deep copy, into a new geom topo tool
  // sets will be duplicated, but entities not
  // modelSet will be a new one;
  // will take as input a pointer to a std::vector of gents (surfaces and volumes, usually),
  // which will serve to filter the gents from modelSet (only dependents will be part of the new gtt)
  // if the pointer is null, all gsets in the original modelSet are duplicated

  ErrorCode duplicate_model(GeomTopoTool *& duplicate, std::vector<EntityHandle> * pvGEnts = NULL);

  EntityHandle get_root_model_set() { return modelSet; }

  bool check_model();
  // should be used instead of keeping multiple ranges, for example in FBEngine
  const Range * geoRanges() { return geomRanges ; }
private:
  Interface *mdbImpl;
  Tag sense2Tag;
  Tag senseNEntsTag, senseNSensesTag;
  Tag geomTag;
  Tag gidTag;
  // the model set encompasses a full topological model
  EntityHandle modelSet;
  Range geomRanges[5];// add one more dimension, for set of gentities; by default, they will
                      // have geom_dimension 4
  int maxGlobalId[5]; // one max global id for each dimension
  bool updated;

  OrientedBoxTreeTool obbTree;
  EntityHandle setOffset;
  std::vector<EntityHandle> rootSets;

  bool contiguous;
  std::map<EntityHandle, EntityHandle>  mapRootSets;
  EntityHandle oneVolRootSet;

    //! compute vertices inclusive and put on tag on sets in geom_sets
  ErrorCode construct_vertex_ranges(const Range &geom_sets,
				      const Tag verts_tag);
  
    //! given a range of geom topology sets, separate by dimension
  ErrorCode separate_by_dimension(const Range &geom_sets);

  // verify sense face tag
  ErrorCode check_face_sense_tag(bool create);

  // verify sense edge tags
  ErrorCode check_edge_sense_tags(bool create);

};

// get the root of the obbtree for a given entity
inline ErrorCode GeomTopoTool::get_root(EntityHandle vol_or_surf, EntityHandle &root) 
{
   if(contiguous)
   {
     unsigned int index = vol_or_surf - setOffset;
     root = (index < rootSets.size() ? rootSets[index] : 0);
   }
   else
      root = mapRootSets[vol_or_surf];
   return (root ? MB_SUCCESS : MB_INDEX_OUT_OF_RANGE);
}

inline EntityHandle GeomTopoTool::get_one_vol_root()
{
  return oneVolRootSet;
}

}
 // namespace moab 

#endif

