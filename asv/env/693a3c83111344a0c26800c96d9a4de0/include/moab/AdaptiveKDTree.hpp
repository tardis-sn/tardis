/**\file AdaptiveKDTree.hpp
 * \class moab::AdaptiveKDTree
 * \brief Adaptive KD tree, for sorting and searching entities spatially
 */

#ifndef MOAB_ADAPTIVE_KD_TREE_HPP
#define MOAB_ADAPTIVE_KD_TREE_HPP

#include "moab/Types.hpp"
#include "moab/Tree.hpp"

#include <string>
#include <vector>
#include <math.h>

namespace moab {

    class AdaptiveKDTreeIter;
    class Interface;
    class Range;

    class AdaptiveKDTree : public Tree
    {
  public:

      AdaptiveKDTree( Interface* iface);

        /** \brief Constructor (build the tree on construction)
         * Construct a tree object, and build the tree with entities input.  See comments
         * for build_tree() for detailed description of arguments.
         * \param iface MOAB instance 
         * \param entities Entities to build tree around
         * \param tree_root Root set for tree (see function description)
         * \param opts Options for tree (see function description)
         */
      AdaptiveKDTree(Interface* iface, const Range &entities, 
                     EntityHandle *tree_root_set = NULL, FileOptions *opts = NULL);

      ~AdaptiveKDTree();

        /** \brief Parse options for tree creation
         * \param options Options passed in by application
         * \return Failure is returned if any options were passed in and not interpreted; could mean
         * inappropriate options for a particular tree type
         */
      ErrorCode parse_options(FileOptions &options);

        /** Build the tree
         * Build a tree with the entities input.  If a non-NULL tree_root_set pointer is input, 
         * use the pointed-to set as the root of this tree (*tree_root_set!=0) otherwise construct 
         * a new root set and pass its handle back in *tree_root_set.  Options vary by tree type;
         * see Tree.hpp for common options; options specific to AdaptiveKDTree:
         * SPLITS_PER_DIR: number of candidate splits considered per direction; default = 3
         * PLANE_SET: method used to decide split planes; see CandidatePlaneSet enum (below)
         *          for possible values; default = 1 (SUBDIVISION_SNAP)
         * \param entities Entities with which to build the tree
         * \param tree_root Root set for tree (see function description)
         * \param opts Options for tree (see function description)
         * \return Error is returned only on build failure
         */
      virtual ErrorCode build_tree(const Range& entities,
                                   EntityHandle *tree_root_set = NULL,
                                   FileOptions *options = NULL);

        //! Reset the tree, optionally checking we have the right root
      virtual ErrorCode reset_tree();

        /** \brief Get leaf containing input position.
         *
         * Does not take into account global bounding box of tree.
         * - Therefore there is always one leaf containing the point.
         * - If caller wants to account for global bounding box, then
         * caller can test against that box and not call this method
         * at all if the point is outside the box, as there is no leaf
         * containing the point in that case.
         * \param point Point to be located in tree
         * \param leaf_out Leaf containing point
         * \param iter_tol Tolerance for convergence of point search
         * \param inside_tol Tolerance for inside element calculation
         * \param multiple_leaves Some tree types can have multiple leaves containing a point;
         *          if non-NULL, this parameter is returned true if multiple leaves contain
         *          the input point
         * \param start_node Start from this tree node (non-NULL) instead of tree root (NULL)
         * \return Non-success returned only in case of failure; not-found indicated by leaf_out=0
         */
      virtual ErrorCode point_search(const double *point,
                                     EntityHandle& leaf_out,
                                     const double iter_tol = 1.0e-10,
                                     const double inside_tol = 1.0e-6,
                                     bool *multiple_leaves = NULL,
                                     EntityHandle *start_node = NULL,
                                     CartVect *params = NULL);

        /** \brief Get leaf containing input position.
         *
         * Does not take into account global bounding box of tree.
         * - Therefore there is always one leaf containing the point.
         * - If caller wants to account for global bounding box, then
         * caller can test against that box and not call this method
         * at all if the point is outside the box, as there is no leaf
         * containing the point in that case.
         * \param point Point to be located in tree
         * \param leaf_it Iterator to leaf containing point
         * \param iter_tol Tolerance for convergence of point search
         * \param inside_tol Tolerance for inside element calculation
         * \param multiple_leaves Some tree types can have multiple leaves containing a point;
         *          if non-NULL, this parameter is returned true if multiple leaves contain
         *          the input point
         * \param start_node Start from this tree node (non-NULL) instead of tree root (NULL)
         * \return Non-success returned only in case of failure; not-found indicated by leaf_out=0
         */
      ErrorCode point_search(const double *point,
                             AdaptiveKDTreeIter& leaf_it,
                             const double iter_tol = 1.0e-10,
                             const double inside_tol = 1.0e-6,
                             bool *multiple_leaves = NULL,
                             EntityHandle *start_node = NULL);
      
        /** \brief Find all leaves within a given distance from point
         * If dists_out input non-NULL, also returns distances from each leaf; if
         * point i is inside leaf, 0 is given as dists_out[i].
         * If params_out is non-NULL and myEval is non-NULL, will evaluate individual entities
         * in tree nodes and return containing entities in leaves_out.  In those cases, if params_out
         * is also non-NULL, will return parameters in those elements in that vector.
         * \param point Point to be located in tree
         * \param distance Distance within which to query
         * \param leaves_out Leaves within distance or containing point
         * \param iter_tol Tolerance for convergence of point search
         * \param inside_tol Tolerance for inside element calculation
         * \param dists_out If non-NULL, will contain distsances to leaves
         * \param params_out If non-NULL, will contain parameters of the point in the ents in leaves_out
         * \param start_node Start from this tree node (non-NULL) instead of tree root (NULL)
         */
      virtual ErrorCode distance_search(const double *point,
                                        const double distance,
                                        std::vector<EntityHandle>& leaves_out,
                                        const double iter_tol = 1.0e-10,
                                        const double inside_tol = 1.0e-6,
                                        std::vector<double> *dists_out = NULL,
                                        std::vector<CartVect> *params_out = NULL,
                                        EntityHandle *start_node = NULL);
      
      ErrorCode get_info(EntityHandle root,
                         double min[3], double max[3], 
                         unsigned int &dep);
  
        //! Enumeriate split plane directions
      enum Axis { X = 0, Y = 1, Z = 2 };
  
        //! Split plane 
      struct Plane {
        double coord;  //!< Location of plane as coordinate on normal axis
        int norm; //!< The principal axis that is the normal of the plane;
    
          /** return true if point is below/to the left of the split plane */
        bool left_side( const double point[3] ) {
          return point[norm] < coord;
        }
          /** return true if point is above/to the right of the split plane */
        bool right_side( const double point[3] ) {
          return point[norm] > coord;
        }
          /** return distance from point to plane */
        double distance( const double point[3] ) const {
          return fabs(point[norm] - coord);
        }
      };
  
        //! Get split plane for tree node
      ErrorCode get_split_plane( EntityHandle node, Plane& plane );
  
        //! Set split plane for tree node
      ErrorCode set_split_plane( EntityHandle node, const Plane& plane );
  
        //! Get iterator for tree
      ErrorCode get_tree_iterator( EntityHandle tree_root,
                                   AdaptiveKDTreeIter& result );
  
        //! Get iterator at right-most ('last') leaf.
      ErrorCode get_last_iterator( EntityHandle tree_root,
                                   AdaptiveKDTreeIter& result );

        //! Get iterator for tree or subtree
      ErrorCode get_sub_tree_iterator( EntityHandle tree_root,
                                       const double box_min[3], 
                                       const double box_max[3],
                                       AdaptiveKDTreeIter& result );

        //! Split leaf of tree
        //! Updates iterator location to point to first new leaf node.
      ErrorCode split_leaf( AdaptiveKDTreeIter& leaf, Plane plane );

        //! Split leaf of tree
        //! Updates iterator location to point to first new leaf node.
      ErrorCode split_leaf( AdaptiveKDTreeIter& leaf, 
                            Plane plane,
                            EntityHandle& left_child,
                            EntityHandle& right_child );
        //! Split leaf of tree
        //! Updates iterator location to point to first new leaf node.
      ErrorCode split_leaf( AdaptiveKDTreeIter& leaf, 
                            Plane plane,
                            const Range& left_entities,
                            const Range& right_entities );

        //! Split leaf of tree
        //! Updates iterator location to point to first new leaf node.
      ErrorCode split_leaf( AdaptiveKDTreeIter& leaf, 
                            Plane plane,
                            const std::vector<EntityHandle>& left_entities,
                            const std::vector<EntityHandle>& right_entities );
  
        //! Merge the leaf pointed to by the current iterator with it's
        //! sibling.  If the sibling is not a leaf, multiple merges may
        //! be done.
      ErrorCode merge_leaf( AdaptiveKDTreeIter& iter );
  
        //! Find triangle closest to input position. 
        //!\param from_coords  The input position to test against
        //!\param closest_point_out  The closest point on the set of triangles in the tree
        //!\param triangle_out The triangle closest to the input position
      ErrorCode closest_triangle( EntityHandle tree_root,
                                  const double from_coords[3],
                                  double closest_point_out[3],
                                  EntityHandle& triangle_out );

      ErrorCode sphere_intersect_triangles( EntityHandle tree_root,
                                            const double center[3],
                                            double radius,
                                            std::vector<EntityHandle>& triangles );

      ErrorCode ray_intersect_triangles( EntityHandle tree_root,
                                         const double tolerance,
                                         const double ray_unit_dir[3],
                                         const double ray_base_pt[3],
                                         std::vector<EntityHandle>& triangles_out,
                                         std::vector<double>& distance_out,
                                         int result_count_limit = 0,
                                         double distance_limit = -1.0);

      ErrorCode compute_depth( EntityHandle root, 
                               unsigned int& min_depth,
                               unsigned int& max_depth );
  
        //! methods for selecting candidate split planes
      enum CandidatePlaneSet {
            //! Candidiate planes at evenly spaced intervals 
          SUBDIVISION=0,
            //! Like SUBDIVISION, except snap to closest vertex coordinate
          SUBDIVISION_SNAP, // = 1
            //! Median vertex coodinate values
          VERTEX_MEDIAN, // = 2
            //! Random sampling of vertex coordinate values
          VERTEX_SAMPLE // = 3
      };
  
        //! print various things about this tree
      virtual ErrorCode print();
      
  private:
      friend class AdaptiveKDTreeIter;

      ErrorCode init();
  
        /**\brief find a triangle near the input point */
      ErrorCode find_close_triangle( EntityHandle root,
                                     const double from_point[3],
                                     double pt[3],
                                     EntityHandle& triangle );

      ErrorCode make_tag( Interface* iface, std::string name, TagType storage, 
                          DataType type, int count, void* default_val, Tag& tag_handle,
                          std::vector<Tag>& created_tags );
  
      ErrorCode intersect_children_with_elems(
          const Range& elems,
          AdaptiveKDTree::Plane plane,
          double eps,
          CartVect box_min,
          CartVect box_max,
          Range& left_tris,
          Range& right_tris,
          Range& both_tris,
          double& metric_value );

      ErrorCode best_subdivision_snap_plane( int num_planes,
                                                    const AdaptiveKDTreeIter& iter,
                                                    Range& best_left,
                                                    Range& best_right,
                                                    Range& best_both,
                                                    AdaptiveKDTree::Plane& best_plane,
                                                    std::vector<double>& tmp_data,
                                                    double eps );
  
      ErrorCode best_subdivision_plane( int num_planes,
                                               const AdaptiveKDTreeIter& iter,
                                               Range& best_left,
                                               Range& best_right,
                                               Range& best_both,
                                               AdaptiveKDTree::Plane& best_plane,
                                               double eps );
  
      ErrorCode best_vertex_median_plane( int num_planes,
                                                 const AdaptiveKDTreeIter& iter,
                                                 Range& best_left,
                                                 Range& best_right,
                                                 Range& best_both,
                                                 AdaptiveKDTree::Plane& best_plane,
                                                 std::vector<double>& coords,
                                                 double eps);
  
      ErrorCode best_vertex_sample_plane( int num_planes,
                                                 const AdaptiveKDTreeIter& iter,
                                                 Range& best_left,
                                                 Range& best_right,
                                                 Range& best_both,
                                                 AdaptiveKDTree::Plane& best_plane,
                                                 std::vector<double>& coords,
                                                 std::vector<EntityHandle>& indices,
                                                 double eps );

      static const char *treeName;
      
      Tag planeTag, axisTag;

      unsigned splitsPerDir;
  
      CandidatePlaneSet planeSet;
    };
                    

//! Iterate over leaves of an adapative kD-tree
    class AdaptiveKDTreeIter
    {
  public:

      enum Direction { LEFT = 0, RIGHT = 1 };

  private:
  
      struct StackObj {
        StackObj( EntityHandle e, double c ) : entity(e), coord(c) {}
        StackObj() {}
        EntityHandle entity; //!< handle for tree node
        double coord;          //!< box coordinate of parent
      };
  
      enum { BMIN = 0, BMAX = 1 };  //!< indices into mBox and child list
  
      CartVect mBox[2];                //!< min and max corners of bounding box
      AdaptiveKDTree* treeTool;       //!< tool for tree
      std::vector<StackObj> mStack;     //!< stack storing path through tree
      mutable std::vector<EntityHandle> childVect; //!< tempory storage of child handles
  
        //! Descend tree to left most leaf from current position
        //! No-op if at leaf.
      ErrorCode step_to_first_leaf( Direction direction );

      friend class AdaptiveKDTree;
  public:

      AdaptiveKDTreeIter() : treeTool(0), childVect(2) {}
  
      ErrorCode initialize( AdaptiveKDTree* tool,
                            EntityHandle root,
                            const double box_min[3],
                            const double box_max[3],
                            Direction direction );

      AdaptiveKDTree* tool() const
          { return treeTool; }

        //! Get handle for current leaf
      EntityHandle handle() const
          { return mStack.back().entity; }
  
        //! Get min corner of axis-aligned box for current leaf
      const double* box_min() const 
          { return mBox[BMIN].array(); }
    
        //! Get max corner of axis-aligned box for current leaf
      const double* box_max() const 
          { return mBox[BMAX].array(); }
  
      double volume() const
          { return (mBox[BMAX][0] - mBox[BMIN][0]) * 
                (mBox[BMAX][1] - mBox[BMIN][1]) * 
                (mBox[BMAX][2] - mBox[BMIN][2]); }
  
        //! test if a plane intersects the leaf box
      bool intersects( const AdaptiveKDTree::Plane& plane ) const
          { return mBox[BMIN][plane.norm] <= plane.coord &&
                mBox[BMAX][plane.norm] >= plane.coord; }
  
        //! Get depth in tree. root is at depth of 1.
      unsigned depth() const
          { return mStack.size(); }
  
        //! Advance the iterator either left or right in the tree
        //! Note:  stepping past the end of the tree will invalidate
        //!        the iterator.  It will *not* be work step the
        //!        other direction.
      ErrorCode step( Direction direction );

        //! Advance to next leaf
        //! Returns MB_ENTITY_NOT_FOUND if at end.
        //! Note: steping past the end of the tree will invalidate
        //!       the iterator. Calling back() will not work.
      ErrorCode step() { return step(RIGHT); }

        //! Move back to previous leaf
        //! Returns MB_ENTITY_NOT_FOUND if at beginning.
        //! Note: steping past the start of the tree will invalidate
        //!       the iterator. Calling step() will not work.
      ErrorCode back() { return step(LEFT); }
  
  
        //! Return the side of the box bounding this tree node
        //! that is shared with the immediately adjacent sibling
        //! (the tree node that shares a common parent node with
        //! this node in the binary tree.)
        //!
        //!\param axis_out The principal axis orthogonal to the side of the box
        //!\param neg_out  true if the side of the box is toward the decreasing
        //!                direction of the principal axis indicated by axis_out,
        //!                false if it is toward the increasing direction.
        //!\return MB_ENTITY_NOT FOUND if root node.
        //!        MB_FAILURE if internal error.
        //!        MB_SUCCESS otherwise.
      ErrorCode sibling_side( AdaptiveKDTree::Axis& axis_out, bool& neg_out ) const;

        //! Get adjacent leaf nodes on side indicated by norm and neg.
        //!
        //! E.g. if norm == X and neg == true, then get neighbor(s)
        //! adjacent to the side of the box contained in the plane
        //! with normal to the X axis and with the x coordinate equal 
        //! to the minimum x of the bounding box.
        //!
        //! E.g. if norm == Y and neg == false, then get neighbor(s)
        //! adjacent to the side of the box with y = maximum y of bounding box.
        //!
        //!\param norm  Normal vector for box side (X, Y, or Z)
        //!\param neg   Which of two planes with norm (true->smaller coord, 
        //!             false->larget coord)
        //!\param results List to which to append results.  This function does
        //!             *not* clear existing values in list.
        //!\param epsilon Tolerance on overlap.  A positive value E will
        //!              result in nodes that are separated by as much as E
        //!              to be considered touching.  A negative value -E will
        //!              cause leaves that do not overlap by at least E to be
        //!              considered non-overlapping.  Amongst other things, 
        //!              this value can be used to control whether or not
        //!              leaves adjacent at only their edges or corners are
        //!              returned.
      ErrorCode get_neighbors( AdaptiveKDTree::Axis norm, bool neg,
                               std::vector<AdaptiveKDTreeIter>& results,
                               double epsilon = 0.0 ) const;
  
        //! Get split plane that separates this node from its immediate sibling.
      ErrorCode get_parent_split_plane( AdaptiveKDTree::Plane& plane ) const;
  
        //! Return true if thos node and the passed node share the
        //! same immediate parent.
      bool is_sibling( const AdaptiveKDTreeIter& other_leaf ) const;
  
        //! Return true if thos node and the passed node share the
        //! same immediate parent.
      bool is_sibling( EntityHandle other_leaf ) const;
  
        //! Returns true if calling step() will advance to the
        //! immediate sibling of the current node.  Returns false
        //! if current node is root or back() will move to the 
        //! immediate sibling.
      bool sibling_is_forward( ) const;
  
        //! Find range of overlap between ray and leaf.
        //!
        //!\param ray_point Coordinates of start point of ray
        //!\param ray_vect  Directionion vector for ray such that
        //!                 the ray is defined by r(t) = ray_point + t * ray_vect
        //!                 for t > 0.
        //!\param t_enter   Output: if return value is true, this value
        //!                 is the parameter location along the ray at which
        //!                 the ray entered the leaf.  If return value is false,
        //!                 then this value is undefined.
        //!\param t_exit    Output: if return value is true, this value
        //!                 is the parameter location along the ray at which
        //!                 the ray exited the leaf.  If return value is false,
        //!                 then this value is undefined.
        //!\return true if ray intersects leaf, false otherwise.
      bool intersect_ray( const double ray_point[3],
                          const double ray_vect[3],
                          double& t_enter, double& t_exit ) const;
    };

    inline ErrorCode AdaptiveKDTree::reset_tree()
    {
      return delete_tree_sets();
    }

} // namespace moab 

#endif
