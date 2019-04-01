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



#ifndef MOAB_SKINNER_HPP
#define MOAB_SKINNER_HPP

#include "moab/Forward.hpp"
#include "moab/Range.hpp"
#include <vector>

namespace moab {

class ScdBox;
    
/** \class Skinner
 * \brief Class for constructing and querying skin of a mesh
 * Class for constructing and querying skin of a mesh, defined as the outside lower-dimensional 
 * boundary of a mesh or a given set of entities.  This class provides options for finding the
 * forward- and reverse-oriented members of the skin.  Several methods are available for computing
 * the skin, e.g. using geometric topology sets, vertex-entity adjacencies, or directly from
 * (n-1)-dimensional entities.
 */
class Skinner 
{

  enum direction{FORWARD=1, REVERSE=-1};
protected:
  //! the MB instance that this works with
  Interface* thisMB;

  Tag mDeletableMBTag;
  Tag mAdjTag;
  int mTargetDim;

public:
  //! constructor, takes mdb instance
  Skinner(Interface* mdb) 
    : thisMB(mdb), mDeletableMBTag(0), mAdjTag(0), mTargetDim(0) {}

  //! destructor
  ~Skinner();

  ErrorCode find_geometric_skin(const EntityHandle meshset, Range &forward_target_entities);
  
    /**\brief will accept entities all of one dimension and 
     *        return entities of n-1 dimension; NOTE: get_vertices argument controls whether
     * vertices or entities of n-1 dimension are returned, and only one of these is allowed
     * (i.e. this function returns only vertices or (n-1)-dimensional entities, but not both)
     * \param entities The elements for which to find the skin
     * \param get_vertices If true, vertices on the skin are returned 
     *        in the range, otherwise elements are returned
     * \param output_handles Range holding skin entities returned
     * \param output_reverse_handles Range holding entities on skin which 
     *        are reversed wrt entities
     * \param create_vert_elem_adjs If true, this function will cause 
     *        vertex-element adjacencies to be generated
     * \param create_skin_elements If true, this function will cause creation 
     *        of skin entities, otherwise only skin entities already extant 
     *        will be returned
     */
  ErrorCode find_skin( const EntityHandle meshset,
                       const Range &entities,
                       bool get_vertices,
                       Range &output_handles,
                       Range *output_reverse_handles = 0,
                       bool create_vert_elem_adjs = false,
                       bool create_skin_elements = true,
                       bool look_for_scd = false);

    /**\brief will accept entities all of one dimension and 
     *        return entities of n-1 dimension; NOTE: get_vertices argument controls whether
     * vertices or entities of n-1 dimension are returned, and only one of these is allowed
     * (i.e. this function returns only vertices or (n-1)-dimensional entities, but not both)
     * \param entities Pointer to elements for which to find the skin
     * \param num_entities Number of entities in vector
     * \param get_vertices If true, vertices on the skin are returned 
     *        in the range, otherwise elements are returned
     * \param output_handles Range holding skin entities returned
     * \param output_reverse_handles Range holding entities on skin which 
     *        are reversed wrt entities
     * \param create_vert_elem_adjs If true, this function will cause 
     *        vertex-element adjacencies to be generated
     * \param create_skin_elements If true, this function will cause creation 
     *        of skin entities, otherwise only skin entities already extant 
     *        will be returned
     */
  ErrorCode find_skin( const EntityHandle this_set,
                       const EntityHandle *entities,
                       int num_entities,
                       bool get_vertices,
                       Range &output_handles,
                       Range *output_reverse_handles = 0,
                       bool create_vert_elem_adjs = false,
                       bool create_skin_elements = true,
                       bool look_for_scd = false);

    /**\brief get skin entities of prescribed dimension
     * \param entities The elements for which to find the skin
     * \param dim Dimension of skin entities requested
     * \param skin_entities Range holding skin entities returned
     * \param create_vert_elem_adjs If true, this function will cause 
     *        vertex-element adjacencies to be generated
     */
  ErrorCode find_skin(const EntityHandle this_set,
                      const Range &entities,
                      int dim,
                      Range &skin_entities,
                      bool create_vert_elem_adjs = false,
                      bool create_skin_elements = true);

  ErrorCode classify_2d_boundary( const Range &boundary,
                                     const Range &bar_elements,
                                     EntityHandle boundary_edges,
                                     EntityHandle inferred_edges,
                                     EntityHandle non_manifold_edges,
                                     EntityHandle other_edges,
                                     int &number_boundary_nodes);
  
  //!given a skin of dimension 2, will classify and return edges
  //! as boundary, inferred, and non-manifold, and the rest (other)
  ErrorCode classify_2d_boundary( const Range  &boundary,
                                     const Range  &mesh_1d_elements,
                                     Range  &boundary_edges,
                                     Range  &inferred_edges,
                                     Range  &non_manifold_edges,
                                     Range  &other_edges,
                                     int &number_boundary_nodes);

protected:
  
  ErrorCode initialize();
  
  ErrorCode deinitialize();

  ErrorCode find_skin_noadj( const Range &source_entities,
                               Range &forward_target_entities,
                               Range &reverse_target_entities );

  ErrorCode add_adjacency(EntityHandle entity);
  
  void add_adjacency(EntityHandle entity, const EntityHandle *conn,
                     const int num_nodes);

  ErrorCode remove_adjacency(EntityHandle entity);

  bool entity_deletable(EntityHandle entity);

  void find_match( EntityType type, 
                   const EntityHandle *conn, 
                   const int num_nodes,
                   EntityHandle& match,
                   Skinner::direction &direct);

  bool connectivity_match(const EntityHandle *conn1,
                          const EntityHandle *conn2,
                          const int num_verts,
                          Skinner::direction &direct);

  void find_inferred_edges(Range &skin_boundary,
                           Range &candidate_edges,
                           Range &inferred_edges,
                           double reference_angle_degrees);

  bool has_larger_angle(EntityHandle &entity1,
                       EntityHandle &entity2,
                       double reference_angle_cosine);


    /**\brief Find vertices on the skin of a set of mesh entities.
     *\param entities The elements for which to find the skin.  Range
     *                may NOT contain vertices, polyhedra, or entity sets.
     *                All elements in range must be of the same dimension.
     *\param skin_verts Output: the vertices on the skin.
     *\param skin_elems Optional output: elements representing sides of entities 
     *                    that are on the skin
     *\param create_if_missing If skin_elemts is non-null and this is true, 
     *                    create new elements representing the sides of 
     *                    entities on the skin.  If this is false, skin_elems
     *                    will contain only those skin elements that already
     *                    exist.
     */
  ErrorCode find_skin_vertices( const EntityHandle this_set,
                                const Range& entities,
                                  Range* skin_verts = 0,
                                  Range* skin_elems = 0,
                                  Range* rev_elems = 0,
                                  bool create_if_missing = true,
                                  bool corners_only = false );

  /**\brief Skin edges
   *
   * Return any vertices adjacent to exactly one of the input edges.
   */
  ErrorCode find_skin_vertices_1D( Tag tag,
                                     const Range& edges,
                                     Range& skin_verts );
                                     
  /**\brief Skin faces
   *
   * For the set of face sides (logical edges), return 
   * vertices on such sides and/or edges equivalent to such sides.
   *\param faces  Set of toplogically 2D entities to skin.
   *\param skin_verts If non-NULL, skin vertices will be added to this container.
   *\param skin_edges If non-NULL, skin edges will be added to this container
   *\param reverse_edges If skin_edges is not NULL and this is not NULL, then
   *                  any existing skin edges that are reversed with respect
   *                  to the skin side will be placed in this range instead of
   *                  skin_edges.  Note: this argument is ignored if skin_edges
   *                  is NULL.
   *\param create_edges If true, edges equivalent to face sides on the skin
   *                  that don't already exist will be created.  Note: this
   *                  parameter is honored regardless of whether or not skin
   *                  edges or vertices are returned.
   *\param corners_only If true, only skin vertices that correspond to the
   *                  corners of sides will be returned (i.e. no higher-order
   *                  nodes.)  This argument is ignored if skin_verts is NULL.
   */
  ErrorCode find_skin_vertices_2D(const EntityHandle this_set, Tag tag,
                                     const Range& faces,
                                     Range* skin_verts = 0,
                                     Range* skin_edges = 0,
                                     Range* reverse_edges = 0,
                                     bool create_edges = false,
                                     bool corners_only = false );
                                     
  /**\brief Skin volume mesh
   *
   * For the set of element sides (logical faces), return 
   * vertices on such sides and/or faces equivalent to such sides.
   *\param entities  Set of toplogically 3D entities to skin.
   *\param skin_verts If non-NULL, skin vertices will be added to this container.
   *\param skin_faces If non-NULL, skin faces will be added to this container
   *\param reverse_faces If skin_faces is not NULL and this is not NULL, then
   *                  any existing skin faces that are reversed with respect
   *                  to the skin side will be placed in this range instead of
   *                  skin_faces.  Note: this argument is ignored if skin_faces
   *                  is NULL.
   *\param create_faces If true, face equivalent to sides on the skin
   *                  that don't already exist will be created.  Note: this
   *                  parameter is honored regardless of whether or not skin
   *                  faces or vertices are returned.
   *\param corners_only If true, only skin vertices that correspond to the
   *                  corners of sides will be returned (i.e. no higher-order
   *                  nodes.)  This argument is ignored if skin_verts is NULL.
   */
  ErrorCode find_skin_vertices_3D(const EntityHandle this_set, Tag tag,
                                     const Range& entities,
                                     Range* skin_verts = 0,
                                     Range* skin_faces = 0,
                                     Range* reverse_faces = 0,
                                     bool create_faces = false,
                                     bool corners_only = false );

  ErrorCode create_side(const EntityHandle this_set, EntityHandle element,
                           EntityType side_type,
                           const EntityHandle* side_corners,
                           EntityHandle& side_elem_handle_out );
                           
  bool edge_reversed( EntityHandle face, const EntityHandle edge_ends[2] );
  bool face_reversed( EntityHandle region, const EntityHandle* face_conn, 
                      EntityType face_type );

    //! look for structured box comprising source_entities, and if one is found use
    //! structured information to find the skin
  ErrorCode find_skin_scd(const Range& source_entities,
                          bool get_vertices,
                          Range& output_handles,
                          bool create_skin_elements);
  
    //! skin a structured box, taking advantage of structured information
  ErrorCode skin_box(ScdBox *box, bool get_vertices, Range &output_handles, 
                     bool create_skin_elements);
};

inline ErrorCode Skinner::find_skin( const EntityHandle this_set,
                                     const EntityHandle *entities,
                                     int num_entities,
                                     bool get_vertices,
                                     Range &output_handles,
                                     Range *output_reverse_handles,
                                     bool create_vert_elem_adjs,
                                     bool create_skin_elements,
                                     bool look_for_scd) 
{
  Range ents;
  std::copy(entities, entities+num_entities, range_inserter(ents));
  return find_skin(this_set, ents, get_vertices, output_handles, output_reverse_handles,
                   create_vert_elem_adjs, create_skin_elements, look_for_scd);
}
      
} // namespace moab 

#endif

