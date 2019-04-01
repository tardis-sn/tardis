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


#ifndef MOAB_DUAL_TOOL_HPP
#define MOAB_DUAL_TOOL_HPP

#include "moab/Forward.hpp"

namespace moab {

/*!
 *  \authors Tim Tautges
 *  \date    2/04
 *  \brief   Tools for constructing and working with mesh duals (both tet- and hex-based,
 *           though some functions may not make sense for tet duals)
 *          
 */ 
class DualTool
{
public:
    //! tag name for dual surfaces
  static const char *DUAL_SURFACE_TAG_NAME;

    //! tag name for dual curves
  static const char *DUAL_CURVE_TAG_NAME;

    //! tag name for dual cells
  static const char *IS_DUAL_CELL_TAG_NAME;

    //! tag name for dual entitys
  static const char *DUAL_ENTITY_TAG_NAME;

    //! tag name for dual entitys
  static const char *EXTRA_DUAL_ENTITY_TAG_NAME;

    //! tag name for dual entitys
  static const char *DUAL_GRAPHICS_POINT_TAG_NAME;

    //! struct for storing a graphics pt
  class GraphicsPoint 
  {
  public:
    GraphicsPoint() 
      {xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0; id = -1;}

    GraphicsPoint(float xi, float yi, float zi, int idi) 
      {xyz[0] = xi; xyz[1] = yi; xyz[2] = zi; id = idi;}

    GraphicsPoint(float xyzi[3], int idi)
      {xyz[0] = xyzi[0]; xyz[1] = xyzi[1]; xyz[2] = xyzi[2]; id = idi;}

    GraphicsPoint(double xyzi[3], int idi)
      {xyz[0] = xyzi[0]; xyz[1] = xyzi[1]; xyz[2] = xyzi[2]; id = idi;}

    GraphicsPoint(const GraphicsPoint &gp) 
      {xyz[0] = gp.xyz[0]; xyz[1] = gp.xyz[1]; xyz[2] = gp.xyz[2]; id = gp.id;}
    
    float xyz[3];
    int id;
  };
  
  DualTool(Interface *impl);
  
  ~DualTool();

    //! construct the dual entities for the entire mesh
  ErrorCode construct_dual(EntityHandle *entities, 
                             const int num_entities);
  
    //! construct the dual entities for a hex mesh, including dual surfaces & curves
  ErrorCode construct_hex_dual(EntityHandle *entities,
                                 const int num_entities);
  
    //! construct the dual entities for a hex mesh, including dual surfaces & curves
  ErrorCode construct_hex_dual(Range &entities);
  
  //! get the dual entities; if non-null, only dual of entities passed in are returned
  ErrorCode get_dual_entities(const int dim, 
                                EntityHandle *entities, 
                                const int num_entities, 
                                Range &dual_ents);
  
  //! get the dual entities; if non-null, only dual of entities passed in are returned
  ErrorCode get_dual_entities(const int dim, 
                                EntityHandle *entities, 
                                const int num_entities, 
                                std::vector<EntityHandle> &dual_ents);

    //! return the corresponding dual entity
  EntityHandle get_dual_entity(const EntityHandle this_ent) const;
  
    //! return the corresponding extra dual entity
  EntityHandle get_extra_dual_entity(const EntityHandle this_ent);
  
    //! get the d-dimensional hyperplane sets; static 'cuz it's easy to do without an active
    //! dualtool
  static ErrorCode get_dual_hyperplanes(const Interface *impl, const int dim, 
                                          Range &dual_ents);

    //! get the graphics points for single entity (dual_ent CAN'T be a set);
    //! returns multiple facets, each with npts[i] points
  ErrorCode get_graphics_points(EntityHandle dual_ent,
                                  std::vector<int> &npts,
                                  std::vector<GraphicsPoint> &gpoints);
  
    //! get the graphics points for a range of entities or sets (if set, the
    //! entities in those sets); optionally reset ids on points
  ErrorCode get_graphics_points(const Range &in_range,
                                  std::vector<GraphicsPoint> &gpoints,
                                  const bool assign_ids = false,
                                  const int start_id = 0);
  
    //! given a last_v (possibly zero) and this_v, find the next loop vertex on 
    //! this dual surface
  EntityHandle next_loop_vertex(const EntityHandle last_v,
                                  const EntityHandle this_v,
                                  const EntityHandle dual_surf);
  
    //! get/set the tag for dual surfaces
  Tag dualSurface_tag() const;
  ErrorCode dualSurface_tag(const Tag tag);
  
    //! get/set the tag for dual curves
  Tag dualCurve_tag() const;
  ErrorCode dualCurve_tag(const Tag tag);

    //! get/set the tag for dual cells
  Tag isDualCell_tag() const;
  ErrorCode isDualCell_tag(const Tag tag);

    //! get/set the tag for dual entities
  Tag dualEntity_tag() const;
  ErrorCode dualEntity_tag(const Tag tag);

    //! get/set the tag for dual entities
  Tag extraDualEntity_tag() const;
  ErrorCode extraDualEntity_tag(const Tag tag);

    //! get/set the tag for dual entities
  Tag dualGraphicsPoint_tag() const;
  ErrorCode dualGraphicsPoint_tag(const Tag tag);

    //! get/set the global id tag
  Tag globalId_tag() const;
  ErrorCode globalId_tag(const Tag tag);

    //! given an entity, return any dual surface or curve it's in
  EntityHandle get_dual_hyperplane(const EntityHandle ncell);

    //! returns true if first & last vertices are dual to hexes (not faces)
  bool is_blind(const EntityHandle chord);
  
    //! set the dual surface or curve for an entity
  ErrorCode set_dual_surface_or_curve(EntityHandle entity, 
                                        const EntityHandle dual_hyperplane,
                                        const int dimension);
  
    //! effect atomic pillow operation
  ErrorCode atomic_pillow(EntityHandle odedge, EntityHandle &quad1,
                            EntityHandle &quad2);

    //! effect reverse atomic pillow operation
  ErrorCode rev_atomic_pillow(EntityHandle pillow, Range &chords);

    //! effect face shrink operation
  ErrorCode face_shrink(EntityHandle odedge);
  
    //! effect reverse atomic pillow operation
  ErrorCode rev_face_shrink(EntityHandle edge);
  
    //! effect a face open-collapse operation
  ErrorCode face_open_collapse(EntityHandle ocl, EntityHandle ocr);

    //! given the two 1-cells involved in the foc, get entities associated with
    //! the quads being opened/collapsed; see implementation for more details
  ErrorCode foc_get_ents(EntityHandle ocl, 
                           EntityHandle ocr, 
                           EntityHandle *quads, 
                           EntityHandle *split_edges, 
                           EntityHandle *split_nodes, 
                           Range &hexes, 
                           EntityHandle *other_edges, 
                           EntityHandle *other_nodes);
  
    //! given a 1-cell and a chord, return the neighboring vertices on the
    //! chord, in the same order as the 1-cell's vertices
  ErrorCode get_opposite_verts(const EntityHandle middle_edge, 
                                 const EntityHandle chord, 
                                 EntityHandle *verts);

    //! given a dual surface or curve, return the 2-cells, 1-cells, 0-cells, and
    //! loop 0/1-cells, if requested; any of those range pointers can be NULL,
    //! in which case that range isn't returned
  ErrorCode get_dual_entities(const EntityHandle dual_ent,
                                Range *dcells,
                                Range *dedges,
                                Range *dverts,
                                Range *dverts_loop,
                                Range *dedges_loop);
  
  ErrorCode list_entities(const Range &entities) const;
  ErrorCode list_entities(const EntityHandle *entities,
                            const int num_entities) const;
  
    //! delete all the dual data
  ErrorCode delete_whole_dual();
  
    //! check dual-primal adjacencies
  ErrorCode check_dual_adjs();

private:

    //! construct dual vertices for specified regions
  ErrorCode construct_dual_vertices(const Range &all_regions,
                                      Range &new_dual_ents);
  
    //! construct dual edges for specified faces
  ErrorCode construct_dual_edges(const Range &all_faces,
                                      Range &new_dual_ents);
  
    //! construct dual faces for specified edges
  ErrorCode construct_dual_faces(const Range &all_edges,
                                      Range &new_dual_ents);
  
    //! construct dual cells for specified vertices
  ErrorCode construct_dual_cells(const Range &all_verts,
                                   Range &new_dual_ents);
  
    //! traverse dual faces of input dimension, constructing
    //! dual hyperplanes of them in sets as it goes
  ErrorCode construct_dual_hyperplanes(const int dim, 
                                         EntityHandle *entities, 
                                         const int num_entities);

    //! order 1cells on a chord 
  ErrorCode order_chord(EntityHandle chord_set);
  
    //! make a new dual hyperplane with the specified id; if the id specified is -1,
    //! set the new one's id to the max found
  ErrorCode construct_new_hyperplane(const int dim, EntityHandle &new_hyperplane,
                                       int &id);
  
    //! traverse the cells of a dual hyperplane, starting with this_ent (dimension
    //! of this_ent determines hyperplane dimension)
    //! simpler method for traversing hyperplane, using same basic algorithm but
    //! using MeshTopoUtil::get_bridge_adjacencies
  ErrorCode traverse_hyperplane(const Tag hp_tag, 
                                  EntityHandle &this_hp, 
                                  EntityHandle this_ent);
  
    //! connect dual surfaces with dual curves using parent/child connections
  ErrorCode construct_hp_parent_child();
  
  //! given an edge handle, return a list of dual vertices in radial order 
  //! around the edge; also returns whether this edge is on the boundary
  ErrorCode get_radial_dverts(const EntityHandle edge,
                                std::vector<EntityHandle> &rad_verts,
                                bool &bdy_edge);
  
  ErrorCode construct_dual_vertex(EntityHandle entity, 
                                    EntityHandle &dual_ent, 
                                    const bool extra = false,
                                    const bool add_graphics_pt = true);

    //! add a graphics point to an entity (on a tag)
  ErrorCode add_graphics_point(EntityHandle entity,
                                 double *avg_pos = NULL);
  
    //! get points defining facets of a 2cell
  ErrorCode get_cell_points(EntityHandle dual_ent,
                              std::vector<int> &npts,
                              std::vector<GraphicsPoint> &points);

    //! if this_ent is an edge, is a dual entity, and has quads as
    //! its vertices' dual entities, return true, otherwise false
  bool check_1d_loop_edge(EntityHandle this_ent);
  
    //! go through potential dual equivalent edges (edges whose nodes define
    //! multiple edges), and add explicit adjacencies to corrent 2cells
  ErrorCode check_dual_equiv_edges(Range &dual_edges);
  
    //! delete a dual entity; updates primal to no longer point to it
  ErrorCode delete_dual_entities(EntityHandle *entities, const int num_entities);

    //! delete a range of dual entities; updates primal to no longer point to them
  ErrorCode delete_dual_entities(Range &entities);
  
    //! check sense of connect arrays, and reverse/rotate if necessary
  ErrorCode fs_check_quad_sense(EntityHandle hex0,
                                  EntityHandle quad0,
                                  std::vector<EntityHandle> *connects);
  
    //! get the three quads for a face shrink, the two hexes, and the connectivity
    //! of the three quads
  ErrorCode fs_get_quads(EntityHandle odedge, 
                           EntityHandle *quads,
                           EntityHandle *hexes,
                           std::vector<EntityHandle> *connects);
  
    //! get loops of quads around 2 hexes, ordered similarly to vertex loops
  ErrorCode  fs_get_quad_loops(EntityHandle *hexes, 
                                 std::vector<EntityHandle> *connects, 
                                 std::vector<EntityHandle> *side_quads);
  
    //! given connectivity of first 3 quads for reverse face shrink, 
    //! get fourth (outer 4 verts to be shared by two inner hexes) and quads
    //! around the side of the structure
  ErrorCode fsr_get_fourth_quad(std::vector<EntityHandle> *connects,
                                  std::vector<EntityHandle> *side_quads);

/*  
    //! get pairs of entities to be merged as part of foc operation
  ErrorCode foc_get_merge_ents(EntityHandle *quads, EntityHandle *new_quads, 
                                 Range &edge, Range &new_edge,
                                 std::vector<EntityHandle> &merge_ents);
*/

    //! function for deleting dual prior to foc operation; special because in
    //! many cases need to delete a sheet in preparation for merging onto another
  ErrorCode foc_delete_dual(EntityHandle *split_quads,
                              EntityHandle *split_edges,
                              Range &hexes);

    //! split a pair of quads and the edge(s) shared by them
  ErrorCode split_pair_nonmanifold(EntityHandle *split_quads,
                                     EntityHandle *split_edges,
                                     EntityHandle *split_nodes,
                                     std::vector<EntityHandle> *star_dp1,
                                     std::vector<EntityHandle> *star_dp2,
                                     EntityHandle *other_edges,
                                     EntityHandle *other_nodes,
                                     EntityHandle *new_quads,
                                     EntityHandle *new_edges,
                                     EntityHandle *new_nodes);
  
    //! for foc's splitting two shared edges, there might be additional entities
    //! connected to the split node that also have to be updated
  ErrorCode foc_get_addl_ents(std::vector<EntityHandle> *star_dp1, 
                                std::vector<EntityHandle> *star_dp2, 
                                EntityHandle *split_edges,
                                EntityHandle split_node,
                                Range *addl_ents);
  
    //! given the split quads and edges, get the face and hex stars around the
    //! edge(s), separated into halves, each of which goes with the new or old entities
    //! after the split
  ErrorCode foc_get_stars(EntityHandle *split_quads,
                            EntityHandle *split_edges,
                            std::vector<EntityHandle> *star_dp1,
                            std::vector<EntityHandle> *star_dp2);
  
  void print_cell(EntityHandle cell);
  
    //! private copy of interface *
  Interface *mbImpl;

    //! static constant number of points bounding any cell
  enum
  {
    GP_SIZE=20
  };
  
    //! tags used for dual surfaces, curves, cells, entities
  Tag dualCurveTag;
  Tag dualSurfaceTag;
  Tag isDualCellTag;
  Tag dualEntityTag;
  Tag extraDualEntityTag;
  Tag dualGraphicsPointTag;
  Tag categoryTag;
  Tag globalIdTag;

  int maxHexId;
};

inline Tag DualTool::dualSurface_tag() const
{
  return dualSurfaceTag;
}

inline Tag DualTool::dualCurve_tag() const
{
  return dualCurveTag;
}

inline Tag DualTool::isDualCell_tag() const
{
  return isDualCellTag;
}

inline Tag DualTool::dualEntity_tag() const
{
  return dualEntityTag;
}

inline Tag DualTool::extraDualEntity_tag() const
{
  return extraDualEntityTag;
}

inline Tag DualTool::dualGraphicsPoint_tag() const
{
  return dualGraphicsPointTag;
}

inline Tag DualTool::globalId_tag() const
{
  return globalIdTag;
}

  //! get/set the tag for dual surfaces
inline ErrorCode DualTool::dualSurface_tag(const Tag tag) 
{
  ErrorCode result = MB_FAILURE;
  if ((0 == dualSurfaceTag && tag) || dualSurfaceTag != tag) {
    dualSurfaceTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual curves
inline ErrorCode DualTool::dualCurve_tag(const Tag tag)
{
  ErrorCode result = MB_FAILURE;
  if ((0 == dualCurveTag && tag) || dualCurveTag != tag) {
    dualCurveTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual cells
inline ErrorCode DualTool::isDualCell_tag(const Tag tag)
{
  ErrorCode result = MB_FAILURE;
  if ((0 == isDualCellTag && tag) || isDualCellTag != tag) {
    isDualCellTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual entities
inline ErrorCode DualTool::dualEntity_tag(const Tag tag)
{
  ErrorCode result = MB_FAILURE;
  if ((0 == dualEntityTag && tag) || dualEntityTag != tag) {
    dualEntityTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual entities
inline ErrorCode DualTool::extraDualEntity_tag(const Tag tag)
{
  ErrorCode result = MB_FAILURE;
  if ((0 == extraDualEntityTag && tag) || extraDualEntityTag != tag) {
    extraDualEntityTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual entities
inline ErrorCode DualTool::dualGraphicsPoint_tag(const Tag tag)
{
  ErrorCode result = MB_FAILURE;
  if ((0 == dualGraphicsPointTag && tag) || dualGraphicsPointTag != tag) {
    dualGraphicsPointTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual entities
inline ErrorCode DualTool::globalId_tag(const Tag tag)
{
  ErrorCode result = MB_FAILURE;
  if ((0 == globalIdTag && tag) || globalIdTag != tag) {
    globalIdTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}

} // namespace moab 

#endif

