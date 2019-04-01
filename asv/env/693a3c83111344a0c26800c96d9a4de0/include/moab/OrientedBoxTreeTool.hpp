/*
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

/**\file OrientedBoxTreeTool.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 */

#ifndef MOAB_ORIENTED_BOX_TREE_TOOL_HPP
#define MOAB_ORIENTED_BOX_TREE_TOOL_HPP

#include "moab/Forward.hpp"

#include <iosfwd>
#include <list>
#include <vector>

namespace moab {

class Range;
class OrientedBox;
class StatData;
class CartVect;

/** \class OrientedBoxTreeTool
 * \brief Class for constructing and querying Hierarchical Oriented Bounding Box trees
 */
class OrientedBoxTreeTool
{
  public:
  
    /**\brief Misc. knobs controlling tree subdivision
     *
     * Available settings for controlling when and how nodes in the tree
     * are split.  The constructor will initialize to the default
     * settings.  All settings except best_split_ratio control when
     * a node is subdivided.  best_split_ratio influences the choice
     * of how the node is subdivided.
     *
     * A calculated ratio is used in the determination of when and how
     * to split a node.  The ratio is calculated as:
     * - \f$max(\frac{|n_L - n_R|}{n_L+n_R}, f*\frac{n_I}{n_L+n_R})\f$
     * - \f$n_L\f$ : num entities to be placed in left child
     * - \f$n_R\f$ : num entities to be placed in right child
     * - \f$f\f$ : Settings::intersect_ratio_factor
     * - \f$n_I\f$: num entities intersecting split plane
     *
     * ALL of the following conditions must be met for a node to be further
     * subdivied:
     *  - Depth must be less than max_depth
     *  - Node must contain more than max_leaf_entities entities.
     *  - The 'ratio' must be less than worst_split_ratio
     *
     * The node will be subdivided using a plane normal to one of the
     * box axis and containing the box center.  The planes are tested 
     * beginning with the one orthogonal to the longest box axis and
     * finishing with the one orthogonal to the shortest box axis.  The
     * search will stop at the first plane for which the 'ratio' is
     * at least Settings::best_split_ratio .  Giving Settings::best_split_ratio
     * a non-zero value gives preference to a split orthogonal to larger
     * box dimensions.
     */
    struct Settings {
      public:
        Settings();              //!< set defaults
        int max_leaf_entities;   //!< Average number of entities per leaf
        int max_depth;           //!< Maximum tree depth - 0->no limit
        //! Must be in [best_split_ratio,1.0]
        //! A tree node will not be split if the ratio of children
        //! in the child nodes is greater than this value.
        double worst_split_ratio;
        //! Must be in [0.0,worst_split_ratio]
        //! The search for an optimal split plane for splitting a node
        //! will stop if at least this ratio is achieved for the number of
        //! entities on each side of the split plane.
        double best_split_ratio;
        //! Flags used to create entity sets representing tree nodes
        unsigned int set_options;
        //! Check if settings are valid.
        bool valid() const;
    };
  
    OrientedBoxTreeTool( Interface* i, 
                           const char* tag_name = 0,
                           bool destroy_created_trees = false ) ;
  
    ~OrientedBoxTreeTool();
  
    /**\brief Build oriented bounding box tree
     *
     * Build an oriented bounding box tree.  
     *\param entities A list of either vertices or 2-D elements (not both)
     *                for which to build a tree.
     *\param set_handle_out A handle for the entity set representing the
     *                root of the tree.
     */
    ErrorCode build( const Range& entities, 
                       EntityHandle& set_handle_out,
                       const Settings* settings = 0 );
     
    /**\brief Build a tree of sets, where each set contains triangles.
     *
     * Build a tree of sets.  Each set must contain at least one triangle
     * to define its geometry.  Each passed set will become a leaf of
     * the OBB tree.  Settings controlling tree depth are ignored by
     * this method.  The tree will be as deep as it needs to be for each
     * input set to be a leaf.
     *
     * To build a tree representing the surfaces of a geometric volume,
     * 1) Build an OBB tree for each surface using the 'build' method
     * 2) Add each surface to the contents of the resulting OBB tree root set
     * 3) Build a tree from all the surface OBB tree root sets using this
     *    method to get a combined tree for the volume.
     */
    ErrorCode join_trees( const Range& tree_roots,
                            EntityHandle& root_set_out,
                            const Settings* settings = 0 );


    /**\brief Traversal statistics structure
     *
     * Structure to accumulate statistics on traversal performance. Passed optionally
     * to query functions, this structure contains the count of nodes visited at each
     * level in a tree, and the count of traversals that ended at each level.
     * One TrvStats structure can be used with multiple OBB trees or multiple queries,
     * or used on only a single tree or a single query.
     *
     * Note that these traversal statistics are not related to the stats() query below,
     * which calculates static information about a tree.  These statistics relate
     * to a tree's dynamic behavior on particular operations.
     */
    class TrvStats{ 
      public:

        //! return counts of nodes visited, indexed by tree depth.  
        //! the counts include both leaves and interior nodes
        const std::vector< unsigned >& nodes_visited() const
        { return nodes_visited_count; }
        //! return counts of tree leaves visited, indexed by tree depth
        const std::vector< unsigned >& leaves_visited() const
        { return leaves_visited_count; }
        //! return counts of traversals ended, indexed by tree depth
        const std::vector< unsigned >& traversals_ended() const 
        { return traversals_ended_count; }        
        //! return total number of ray-triangle intersection tests performed
        //! in calls made with this TrvStats
        unsigned int ray_tri_tests() const
        { return ray_tri_tests_count; }
        //! reset all counters on this structure
        void reset();
        //! print the contents of this structure to given stream
        void print( std::ostream& str ) const ;

        TrvStats() : ray_tri_tests_count(0) {}

      private: 
      
        std::vector< unsigned > nodes_visited_count;
        std::vector< unsigned > leaves_visited_count;
        std::vector< unsigned > traversals_ended_count;
        unsigned int ray_tri_tests_count;

        void increment( unsigned depth );
        void increment_leaf( unsigned depth );
        void end_traversal( unsigned depth );

      friend class OrientedBoxTreeTool;

    };


    /**\brief Intersect a ray with the triangles contained within the tree
     *
     * Intersect a ray with the triangles contained in the tree and return
     * the distance at which the intersection occured.
     *\param distances_out The output list of intersection points on the ray.
     *\param facets_out    Handles of intersected triangles corresponding to distances_out 
     *\param root_set      The MBENTITYSET representing the root of the tree.
     *\param tolerance     The tolerance to use in intersection checks.
     *\param ray_point     The base point of the ray.
     *\param unit_ray_dir  The ray direction vector (must be unit length)
     *\param ray_length    Optional ray length (intersect segment instead of ray.)
     */
    ErrorCode ray_intersect_triangles( std::vector<double>& distances_out,
                                         std::vector<EntityHandle>& facets_out,
                                         EntityHandle root_set,
                                         double tolerance,
                                         const double ray_point[3],
                                         const double unit_ray_dir[3],
                                         const double* ray_length = 0,
                                         TrvStats* accum = 0 );
    
    /**\brief Intersect ray with tree
     *
     * Return the tree nodes (as MBENTITYSET handles) for the leaf boxes
     * of the tree intersected by a ray.
     *\param boxes_out    The boxes intersected by the ray.
     *\param tolerance     The tolerance to use in intersection checks.
     *\param ray_point     The base point of the ray.
     *\param unit_ray_dir  The ray direction vector (must be unit length)
     *\param ray_length    Optional ray length (intersect segment instead of ray.)
     */
    ErrorCode ray_intersect_boxes( Range& boxes_out,
                                     EntityHandle root_set,
                                     double tolerance,
                                     const double ray_point[3],
                                     const double unit_ray_dir[3],
                                     const double* ray_length = 0,
                                     TrvStats* accum = 0 );

    /**\brief Intersect ray with triangles contained in passed MBENTITYSETs 
     * 
     * \param raytri_test_count    If non-NULL, count of ray-triangle intersect tests
     *                             will be added to the value at which this points.
     */
    ErrorCode ray_intersect_triangles( 
                          std::vector<double>& intersection_distances_out,
                          std::vector<EntityHandle>& intersection_facets_out,
                          const Range& leaf_boxes_containing_tris,
                          double tolerance,
                          const double ray_point[3],
                          const double unit_ray_dir[3],
                          const double* ray_length = 0, 
                          unsigned int* raytri_test_count = 0);
                          

    /**\brief Intersect a ray with the triangles contained within the tree
     *
     * Intersect a ray with the triangles contained in the tree and return
     * the distance at which the intersection occured.
     *\param distances_out The output list of intersection points on the ray.
     *\param sets_out      The contained set encountered during the tree traversal
     *                     (see 'set_build').  For the most common use, this is the
     *                     set corresponding to the geometric surface containing the
     *                     intersected triangle.
     *\param facets_out    Handles of intersected triangles corresponding to distances_out 
     *\param root_set      The MBENTITYSET representing the root of the tree.
     *\param min_tolerance_intersections This method returns all intersections
     *                     within 'tolerance' of the start of the ray and if 
     *                     the number of intersections within the 'tolerance' of the
     *                     ray start point is less than this number, the next closest
     *                     intersection.  If the desired result is only the closest
     *                     intersection, pass zero for this argument.
     *                     This function will return all intersections, regardless
     *                     of distance from the start of the ray, if this value
     *                     is negative.
     *\param tolerance     The tolerance to use in intersection checks.
     *\param ray_point     The base point of the ray.
     *\param unit_ray_dir  The ray direction vector (must be unit length)
     *\param nonneg_ray_len Optional ray length ahead of the ray_point (intersect 
     *                     segment instead of ray.)
     *\param accum         Optional class for tree traversal statistics.
     *\param neg_ray_len   Optional ray length behind the ray_point to search for 
     *                     intersections.
     *\param geom_vol      Optional handle of the geometry set being searched. When
     *                     used, glancing intersections are rejected. Must be used
     *                     used with sense_tag.
     *\param sense_tag     Must be used if geom_vol is used. Saves >4% of execution
     *                     time by avoiding tag_get_handle call.
     *\param desired_orient Optional ptr used to screen intersections by orientation.
     *                     Pass 1 to keep intersections with surface normals in the
     *                     same direction as the ray. Pass -1 for opposite orientation.
     *                     Requires use of geom_vol.
     *\param prev_facets   Optional vector of triangles that cannot be returned
     *                     as intersections.
     */
    ErrorCode ray_intersect_sets( std::vector<double>&       distances_out,
                                  std::vector<EntityHandle>& sets_out,
                                  std::vector<EntityHandle>& facets_out,
                                  EntityHandle               root_set,
                                  double                     tolerance,
                                  int                        min_tolerace_intersections,
                                  const double               ray_point[3],
                                  const double               unit_ray_dir[3],
                                  const double*              nonneg_ray_len = 0,
                                  TrvStats*                  accum          = 0,
                                  const double*              neg_ray_len    = 0,
                                  const EntityHandle*        geom_vol       = 0,
                                  const Tag*                 sense_tag      = 0,
                                  const int*                 desired_orient = 0,
                                  const std::vector<EntityHandle>* prev_facets = 0 );
    
    /**\brief Find closest surface, facet in surface, and location on facet
     *
     * Find the closest location in the tree to the specified location.
     *\param point Location to search from
     *\param point_out Closest location on closest facet
     *\param facet_out Closest 2D element to input position
     *\param set_out Set containing closest facet.  0 if tree was not 
     *               constructed using 'set_build'
     */
    ErrorCode closest_to_location( const double* point,
                                     EntityHandle tree_root,
                                     double* point_out,
                                     EntityHandle& facet_out,
                                     EntityHandle* set_out = 0, 
                                     TrvStats* accum = 0 );
                                     
    /**\brief Find closest facet(s) to input position.
     *
     * Find the closest location(s) in the tree to the specified location.
     *\param point Location to search from
     *\param facets_out Closest 2D elements to input position are appended to this list
     *\param sets_out If non-null, sets owning facets are appended to this list.
     */
    ErrorCode closest_to_location( const double* point,
                                     EntityHandle tree_root,
                                     double tolerance,
                                     std::vector<EntityHandle>& facets_out,
                                     std::vector<EntityHandle>* sets_out = 0, 
                                     TrvStats* accum = 0 );
    
    /**\brief Find facets intersected by a sphere 
     *
     * Find facets intersected by a spherical volume.
     *\param center     Sphere center
     *\param radius     Sphere radius
     *\param tree_root  Root of OBB tree
     *\param facets_out List of triangles intersecting sphere
     *\param sets_out   If not null, sets owning facets are appended to this
     *                  list in an order corresponding to the entries in 
     *                  facets_out.
     */
    ErrorCode sphere_intersect_triangles( const double* center,
                                        double radius,
                                        EntityHandle tree_root,
                                        std::vector<EntityHandle>& facets_out,
                                        std::vector<EntityHandle>* sets_out = 0, 
                                        TrvStats* accum = 0 );
    
    /**\brief Get oriented box at node in tree
     *
     * Get the oriented box for a node in an oriented bounding box tree.
     */
    ErrorCode box( EntityHandle node_set,
                     double center[3],
                     double axis1[3],
                     double axis2[3],
                     double axis3[3] );
                         
    ErrorCode delete_tree( EntityHandle root_set );

    /**\brief Print out tree
     *
     * Print the tree to an output stream in a human-readable form.
     *\param tree_root_set  Entity set representing tree root.
     *\param list_contents  If true, list entities in each tree node,
     *                      If false, just list number of entities.
     *\param id_tag_name    If specified, must be the name of an existing
     *                      integer tag containing an ID for the entities.
     *                      Not used if list_contents is false.
     */
    void print( EntityHandle tree_root_set, 
                std::ostream& stream,
                bool list_contents = false,
                const char* id_tag_name = 0 );
                
    /**\brief Print tree statistics
     *
     * Print misc. stats. describing tree
     */
    ErrorCode stats( EntityHandle tree_root_set, std::ostream& stream );
  
    /**\brief Get tree statistics
     *
     * Get summary stats. describing tree
     * \param set Root of tree for which data is requested
     * \param total_entities Entities in tree
     * \param root_volume Total volume of root box
     * \param tot_node_volume Total volume in all nodes
     * \param tot_to_root_volume Ratio of total / root volume
     * \param tree_height Maximum height of tree, from root to leaf
     * \param node_count Number of nodes in tree
     * \param num_leaves Number of leaf nodes in tree
     */
  ErrorCode stats( EntityHandle set, 
                     unsigned &entities_in_tree,
                     double &root_volume,
                     double &tot_node_volume,
                     double &tot_to_root_volume,
                     unsigned &tree_height,
                     unsigned &node_count,
                     unsigned &num_leaves);
  
    /** \brief Implement this and pass instance to preorder_traverse
     * 
     * This interface may be implemented and an instance passed to
     * preorder_traverse to define some operation to do when traversing
     * the tree.
     */
    class Op {
      public:

        /**\brief Visit a node in the tree during a traversal.
         *
         * This method is called for each node in the tree visited
         * during a pre-order traversal.  
         *\param node The EntityHandle for the entity set for the tree node.
         *\param depth The current depth in the tree.
         *\param descend Output: if false, traversal will skip children
         *             of the current node, or if the current node is a
         *             leaf, the 'leaf' method will not be called.
         */
        virtual ErrorCode visit( EntityHandle node,
                                 int depth,
				 bool& descend ) = 0;
       
        /**\brief Process a leaf node during tree traversal */
        virtual ErrorCode leaf( EntityHandle node ) = 0;

        virtual ~Op(); // probably isn't necessary in this case, and
                       // does nothing, but lots of compilers warn if
                       // virtual function but no virtual destructor.
    };
    
    /**\brief Visitor pattern - do operation for each tree node
     *
     * Do a preorder traversal of the tree, calling the method
     * in the passed operation instance for each node in the tree.
     * Parent node is visited before either child (pre-order traversal).
     * If operator method passes back the 'descend' argument as false,
     * traversal will not descend to the children of the current node.
     */
    ErrorCode preorder_traverse( EntityHandle root_set,
                                   Op& operation, 
                                   TrvStats* accum = 0 );
  
    Interface* get_moab_instance() const { return instance; }
  
    struct SetData;
    
    /**\brief Get oriented box at node in tree
     *
     * Get the oriented box for a node in an oriented bounding box tree.
     *
     * NOTE: This function is provided for internal MOAB use only.
     *       The definition of OrientedBox is not available as a part
     *       of the MOAB API
     */
    ErrorCode box( EntityHandle node_set,
                     OrientedBox& box );
  private:
  
    ErrorCode build_tree( const Range& entities, 
                            EntityHandle& set, 
                            int depth,
                            const Settings& settings );
  
    ErrorCode build_sets( std::list<SetData>& sets,
                            EntityHandle& node_set,
                            int depth,
                            const Settings& settings );

    ErrorCode recursive_stats( OrientedBoxTreeTool* tool,
                                    Interface* instance,
                                    EntityHandle set,
                                    int depth,
                                    StatData& data,
                                    unsigned& count_out,
                                    CartVect& dimensions_out );
  
    Interface* instance;
    Tag tagHandle;
 
    bool cleanUpTrees;
    std::vector<EntityHandle> createdTrees;
};

} // namespace moab 

#endif
