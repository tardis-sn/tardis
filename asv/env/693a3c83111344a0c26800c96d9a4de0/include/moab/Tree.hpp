/**\file Tree.hpp
 * \class moab::Tree
 * \brief Parent class of various tree types in MOAB
 */

#ifndef MOAB_TREE_HPP
#define MOAB_TREE_HPP

#include "moab/Interface.hpp"
#include "moab/BoundBox.hpp"
#include "moab/CartVect.hpp"
#include "moab/FileOptions.hpp"
#include "moab/TreeStats.hpp"

#include <string>
#include <vector>
#include <math.h>
#include <assert.h>

namespace moab {

    class Interface;
    class Range;
    class ElemEvaluator;

    class Tree
    {
  public:
        /** \brief Constructor (bare)
         * \param iface MOAB instance 
         */
      Tree(Interface* iface);

        /** \brief Destructor
         */
      virtual ~Tree();

        /** \brief Destroy the tree maintained by this object, optionally checking we have the right root.
         * \param root If non-NULL, check that this is the root, return failure if not
         */
      virtual ErrorCode reset_tree() = 0;

        /** \brief Delete the entity sets associated with the tree, starting with the root and traversing children
         */
      ErrorCode delete_tree_sets();
      
        /** Build the tree
         * Build a tree with the entities input.  If a non-NULL tree_root_set pointer is input, 
         * use the pointed-to set as the root of this tree (*tree_root_set!=0) otherwise construct 
         * a new root set and pass its handle back in *tree_root_set.  Options vary by tree type, 
         * with a few common to all types of trees.  Common options:
         * MAX_PER_LEAF: max entities per leaf; default = 6
         * MAX_DEPTH: max depth of the tree; default = 30
         * MIN_WIDTH: minimum width of box, used like a tolerance; default = 1.0e-10
         * MESHSET_FLAGS: flags passed into meshset creation for tree nodes; should be a value from
         *          ENTITY_SET_PROPERTY (see Types.hpp); default = MESHSET_SET
         * CLEAN_UP: if false, do not delete tree sets upon tree class destruction; default = true
         * TAG_NAME: tag name to store tree information on tree nodes; default determined by tree type
         * \param entities Entities with which to build the tree
         * \param tree_root Root set for tree (see function description)
         * \param opts Options for tree (see function description)
         * \return Error is returned only on build failure
         */
      virtual ErrorCode build_tree(const Range& entities,
                                   EntityHandle *tree_root_set = NULL,
                                   FileOptions *options = NULL) = 0;

        /** \brief Get bounding box for tree below tree_node, or entire tree
         * If no tree has been built yet, returns +/- DBL_MAX for all dimensions.  Note for some tree types,
         * boxes are not available for non-root nodes, and this function will return failure if non-root
         * is passed in
         * \param box The box for this tree
         * \param tree_node If non-NULL, node for which box is requested, tree root if NULL
         * \return Only returns error on fatal condition
         */
      virtual ErrorCode get_bounding_box(BoundBox &box, EntityHandle *tree_node = NULL) const;
  
        /** \brief Return some basic information about the tree
         * Stats are returned for tree starting from input node or tree root (root = 0)
         * \param root If non-0, give stats below and including root
         * \param min Minimum corner of bounding box
         * \param max Maximum corner of bounding box
         * \param max_dep Maximum depth of tree below root
         */
      virtual ErrorCode get_info(EntityHandle root,
                                 double min[3], double max[3], 
                                 unsigned int &max_dep);
  
        /** \brief Find all trees, by bounding box tag
         */
      ErrorCode find_all_trees( Range& results );
      
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
                                     CartVect *params = NULL) = 0;

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
                                        EntityHandle *start_node = NULL) = 0;

        /** \brief Return the MOAB interface associated with this tree
         */
      Interface* moab() { return mbImpl; }

        /** \brief Return the MOAB interface associated with this tree
         */
      const Interface* moab() const { return mbImpl; }

        /** \brief Get max depth set on tree */
      double get_max_depth() {return maxDepth;}
      
        /** \brief Get max entities per leaf set on tree */
      double get_max_per_leaf() {return maxPerLeaf;}

        /** \brief Get tree traversal stats object */
      TreeStats &tree_stats() {return treeStats;}
      
        /** \brief Get tree traversal stats object */
      const TreeStats &tree_stats() const {return treeStats;}
      
        /** \brief Create tree root and tag with bounding box
         */
      ErrorCode create_root( const double box_min[3], const double box_max[3],
                             EntityHandle& root_handle);
      
        //! print various things about this tree
      virtual ErrorCode print() = 0;
      
        //! get/set the ElemEvaluator
      inline ElemEvaluator *get_eval() {return myEval;}
      
        //! get/set the ElemEvaluator
      inline void set_eval(ElemEvaluator *eval) {myEval = eval;}
      
        /** \brief Parse options for this tree, including common options for all trees
         * \param opts Options
         */
      virtual ErrorCode parse_options(FileOptions &opts) = 0;

  protected:

        /** \brief Parse options common to all trees
         * \param options Options for representing tree; see Tree::build_tree() and subclass build_tree()
         *          functions for allowed options
         * \return Non-success returned from base class function only under catastrophic circumstances;
         *          derived classes also can recognize subclass-specific options
         */
      ErrorCode parse_common_options(FileOptions &options);

        /** \brief Get the box tag, possibly constructing it first
         * \param create_if_missing If true and it has not been made yet, make it
         */
      Tag get_box_tag(bool create_if_missing = true);

        // moab instance
      Interface *mbImpl;

        // bounding box for entire tree
      BoundBox boundBox;
      
        // max entities per leaf
      int maxPerLeaf;
      
        // max depth of tree
      int maxDepth;
      
        // tree depth, set by build_tree
      int treeDepth;
      
        // min width of box, handled like tolerance
      double minWidth;

        // meshset creation flags
      unsigned int meshsetFlags;

        // clean up flag
      bool cleanUp;

        // tree root
      EntityHandle myRoot;

        // tag used to mark bounding box of nodes
      Tag boxTag;

        // tag name used for boxTag
      std::string boxTagName;

        // tree traversal stats
      TreeStats treeStats;

        // element evaluator
      ElemEvaluator *myEval;
    };

    inline Tree::Tree(Interface* iface) 
            : mbImpl(iface), maxPerLeaf(6), maxDepth(30), treeDepth(-1), minWidth(1.0e-10),
              meshsetFlags(0), cleanUp(true), myRoot(0), boxTag(0), myEval(0)
    {}

    inline Tree::~Tree() 
    {
    }
    
    inline ErrorCode Tree::get_bounding_box(BoundBox &box, EntityHandle *tree_node) const
    {
      if ((tree_node && *tree_node == myRoot) || !tree_node) {
        box = boundBox;
        return MB_SUCCESS;
      }
      else return MB_FAILURE;
    }
  
    inline ErrorCode Tree::get_info(EntityHandle /* root */,
                                    double * /*min[3]*/, double * /* max[3]*/, 
                                    unsigned int &/*dep*/) 
    {
      return MB_NOT_IMPLEMENTED;
    }

    inline Tag Tree::get_box_tag(bool create_if_missing) 
    {
      if (!boxTag && create_if_missing) {
        assert(boxTagName.length() > 0);
        ErrorCode rval = moab()->tag_get_handle(boxTagName.c_str(), 6, MB_TYPE_DOUBLE, boxTag, MB_TAG_CREAT | MB_TAG_SPARSE);
        if (MB_INVALID_SIZE == rval) {
            // delete the tag and get it again, legacy file...
          rval = moab()->tag_delete(boxTag);
          if (MB_SUCCESS != rval) return 0;
          boxTag = 0;
          return get_box_tag(true);
        }
        
        if (MB_SUCCESS != rval) return 0;
      }
      
      return boxTag;
    }
    
} // namespace moab 

#endif
