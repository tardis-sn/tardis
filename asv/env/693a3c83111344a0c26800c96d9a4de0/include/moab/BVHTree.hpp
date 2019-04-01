/**\file BVHTree.hpp
 * \class moab::BVHTree
 * \brief Bounding Volume Hierarchy (sorta like a tree), for sorting and searching entities spatially
 */

#ifndef BVH_TREE_HPP
#define BVH_TREE_HPP

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/BoundBox.hpp"
#include "moab/Tree.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"

#include <vector>
#include <cfloat>
#include <climits>
#include <map>
#include <set>
#include <iostream>

namespace moab {

    class ElemEvaluator;
    
    class BVHTree : public Tree {
  public:
      BVHTree(Interface *impl);

      ~BVHTree() {reset_tree();}
      
        /** \brief Destroy the tree maintained by this object, optionally checking we have the right root.
         * \param root If non-NULL, check that this is the root, return failure if not
         */
      virtual ErrorCode reset_tree();

      virtual ErrorCode parse_options(FileOptions &opts);
      
        /** \brief Get bounding box for tree below tree_node, or entire tree
         * If no tree has been built yet, returns +/- DBL_MAX for all dimensions.  Note for some tree types,
         * boxes are not available for non-root nodes, and this function will return failure if non-root
         * is passed in
         * \param box The box for this tree
         * \param tree_node If non-NULL, node for which box is requested, tree root if NULL
         * \return Only returns error on fatal condition
         */
      virtual ErrorCode get_bounding_box(BoundBox &box, EntityHandle *tree_node = NULL) const;
  
        /** Build the tree
         * Build a tree with the entities input.  If a non-NULL tree_root_set pointer is input, 
         * use the pointed-to set as the root of this tree (*tree_root_set!=0) otherwise construct 
         * a new root set and pass its handle back in *tree_root_set.  Options vary by tree type;
         * see Tree.hpp for common options; options specific to AdaptiveKDTree:
         * SPLITS_PER_DIR: number of candidate splits considered per direction; default = 3
         * CANDIDATE_PLANE_SET: method used to decide split planes; see CandidatePlaneSet enum (below)
         *          for possible values; default = 1 (SUBDIVISION_SNAP)
         * \param entities Entities with which to build the tree
         * \param tree_root Root set for tree (see function description)
         * \param opts Options for tree (see function description)
         * \return Error is returned only on build failure
         */
      virtual ErrorCode build_tree(const Range& entities,
                                   EntityHandle *tree_root_set = NULL,
                                   FileOptions *options = NULL);

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

        /** \brief Find all leaves within a given distance from point
         * If dists_out input non-NULL, also returns distances from each leaf; if
         * point i is inside leaf, 0 is given as dists_out[i].
         * If params_out is non-NULL and myEval is non-NULL, will evaluate individual entities
         * in tree nodes and return containing entities in leaves_out.  In those cases, if params_out
         * is also non-NULL, will return parameters in those elements in that vector.
         * \param from_point Point to be located in tree
         * \param distance Distance within which to query
         * \param result_list Leaves within distance or containing point
         * \param iter_tol Tolerance for convergence of point search
         * \param inside_tol Tolerance for inside element calculation
         * \param result_dists If non-NULL, will contain distsances to leaves
         * \param result_params If non-NULL, will contain parameters of the point in the ents in leaves_out
         * \param tree_root Start from this tree node (non-NULL) instead of tree root (NULL)
         */
      virtual ErrorCode distance_search(const double from_point[3],
                                        const double distance,
                                        std::vector<EntityHandle>& result_list,
                                        const double iter_tol = 1.0e-10,
                                        const double inside_tol = 1.0e-6,
                                        std::vector<double> *result_dists = NULL,
                                        std::vector<CartVect> *result_params = NULL,
                                        EntityHandle *tree_root = NULL);

        //! print various things about this tree
      virtual ErrorCode print();

      
  private:
        // don't allow copy constructor, too complicated
      BVHTree(const BVHTree &s);

      class HandleData 
      {
    public:
        HandleData(EntityHandle h, const BoundBox &bx, const double dp) : myHandle(h), myBox(bx), myDim(dp) {}
        HandleData() : myHandle(0), myDim(-1) {}
        EntityHandle myHandle;
        BoundBox myBox;
        double myDim;
      };
      typedef std::vector<HandleData> HandleDataVec;
    
      class SplitData {
    public:
        SplitData(): dim(UINT_MAX), nl(UINT_MAX), nr(UINT_MAX), split(DBL_MAX), 
                     Lmax(-DBL_MAX), Rmin(DBL_MAX) {}
        SplitData( const SplitData &f): 
        dim(f.dim), nl(f.nl), nr(f.nr), split(f.split), Lmax(f.Lmax), Rmin(f.Rmin), 
        boundingBox(f.boundingBox), leftBox(f.leftBox), rightBox(f.rightBox){}
        unsigned int dim, nl, nr;
        double split;
        double Lmax, Rmin;
        BoundBox boundingBox, leftBox, rightBox;
        SplitData& operator=(const SplitData & f){
          dim = f.dim; nl = f.nl; nr = f.nr;
          split = f.split; Lmax = f.Lmax; Rmin = f.Rmin;
          boundingBox = f.boundingBox; leftBox = f.leftBox; rightBox = f.rightBox;
          return *this;
        }
      }; //SplitData

      class Split_comparator : 
          public std::binary_function< SplitData, SplitData, bool> {
        inline double objective(const SplitData & a) const{
          return a.Lmax*a.nl - a.Rmin*a.nr;
        }
    public:
        bool operator()(const SplitData &a, const SplitData &b) const {
          return  objective( a) < objective( b);
        }
      }; //Split_comparator

      class HandleData_comparator : 
          public std::binary_function<HandleData, HandleData, bool> {
    public:
        bool operator()(const HandleData &a, const HandleData &b) {
          return a.myDim < b.myDim ||
              (a.myDim == b.myDim && 
               a.myHandle < b.myHandle);
        }
      }; //Iterator_comparator

      class Bucket {
    public:
        Bucket() : mySize(0) {}
        Bucket( const Bucket &f): mySize(f.mySize), boundingBox(f.boundingBox) {}
        Bucket(const unsigned int sz): mySize(sz) {}
        static unsigned int bucket_index(int num_splits, const BoundBox &box, const BoundBox &interval, const unsigned int dim);
        unsigned int mySize;
        BoundBox boundingBox;
        Bucket &operator=(const Bucket & f){
          boundingBox = f.boundingBox; mySize = f.mySize;
          return *this;
        }
      }; //Bucket

      class Node {
    public:
        HandleDataVec entities;
        unsigned int dim, child;
        double Lmax, Rmin;
        BoundBox box;
        Node() : dim(UINT_MAX), child(UINT_MAX), Lmax(-DBL_MAX), Rmin(DBL_MAX) {}
        Node &operator=(const Node& f) {
          dim = f.dim; child = f.child;
          Lmax = f.Lmax; Rmin = f.Rmin;
          entities = f.entities;
          return *this;
        }
      }; // Node

      class TreeNode {
    public:
        unsigned int dim, child;
        double Lmax, Rmin;
        BoundBox box;
        TreeNode(int dm, int chld, double lmx, double rmn, BoundBox &bx) 
                : dim(dm), child(chld), Lmax(lmx), Rmin(rmn), box(bx) {}
        TreeNode &operator=(const TreeNode& f) {
          dim = f.dim; child = f.child;
          Lmax = f.Lmax; Rmin = f.Rmin;
          return *this;
        }
      }; // TreeNode

      void establish_buckets(HandleDataVec::const_iterator begin, 
                             HandleDataVec::const_iterator end, 
                             const BoundBox &interval, std::vector<std::vector<Bucket> > &buckets) const;

      unsigned int set_interval(BoundBox & interval, 
                                std::vector<Bucket>::const_iterator begin, 
                                std::vector<Bucket>::const_iterator end) const;

      void initialize_splits(std::vector<std::vector<SplitData> > &splits, 
                             const std::vector<std::vector<Bucket> > &buckets, 
                             const SplitData &data) const;

      void order_elements(HandleDataVec::iterator &begin, 
                          HandleDataVec::iterator &end, 
                          SplitData &data) const;

      void median_order(HandleDataVec::iterator &begin, 
                        HandleDataVec::iterator &end, 
                        SplitData & data) const;

      void choose_best_split(const std::vector<std::vector<SplitData> > &splits,
                             SplitData & data) const;


      void find_split(HandleDataVec::iterator &begin, 
                      HandleDataVec::iterator &end, 
                      SplitData &data) const;

      ErrorCode find_point(const std::vector<double> &point, 
                           const unsigned int &index,
                           const double iter_tol,
                           const double inside_tol,
                           std::pair<EntityHandle, CartVect> &result);
      
      EntityHandle bruteforce_find(const double *point, 
                                   const double iter_tol = 1.0e-10,
                                   const double inside_tol = 1.0e-6);

      int local_build_tree(std::vector<Node> &tree_nodes,
                           HandleDataVec::iterator begin, 
                           HandleDataVec::iterator end,
                           const int index, const BoundBox &box, 
                           const int depth=0);

        // builds up vector of HandleData, which caches elements' bounding boxes
      ErrorCode construct_element_vec(std::vector<HandleData> &handle_data_vec,
                                      const Range &elements, 
                                      BoundBox & bounding_box);
      
        // convert the std::vector<Node> to myTree and a bunch of entity sets
      ErrorCode convert_tree(std::vector<Node> &tree_nodes);

        // print tree nodes
      ErrorCode print_nodes(std::vector<Node> &nodes);
      
      Range entityHandles;
      std::vector<TreeNode> myTree;
      int splitsPerDir;
      EntityHandle startSetHandle;
      static const char *treeName;
    }; //class Bvh_tree

    inline unsigned int BVHTree::Bucket::bucket_index(int num_splits, const BoundBox &box, const BoundBox & interval, const unsigned int dim)
    {
//see FastMemoryEfficientCellLocationinUnstructuredGridsForVisualization.pdf 
//around page 9
      
//Paper arithmetic is over-optimized.. this is safer.
      const double min = interval.bMin[dim];
      const double length = (interval.bMax[dim]-min)/(num_splits+1);
      const double center = ((box.bMax[dim] + box.bMin[dim])/2.0)-min;
#ifndef NDEBUG
#ifdef BVH_SHOW_INDEX
      std::cout << "[" << min << " , " 
                << interval.max[dim] << " ]" <<std::endl;
      std::cout << "[" 
                << box.bMin[dim] << " , " << box.bMax[dim] << " ]" <<std::endl;
      std::cout << "Length of bucket" << length << std::endl;
      std::cout << "Center: " 
                << (box.bMax[dim] + box.bMin[dim])/2.0 << std::endl;
      std::cout << "Distance of center from min:  " << center << std::endl;
      std::cout << "ratio: " << center/length << std::endl;
      std::cout << "index: " << std::ceil(center/length)-1 << std::endl;
#endif
#endif
      unsigned int cl = std::ceil(center/length);
      return (cl > 0 ? cl-1 : 0);
    }

    inline BVHTree::BVHTree(Interface *impl) : 
            Tree(impl), splitsPerDir(3), startSetHandle(0) {boxTagName = treeName;}

    inline unsigned int BVHTree::set_interval(BoundBox &interval, 
                                              std::vector<Bucket>::const_iterator begin, 
                                              std::vector<Bucket>::const_iterator end) const
    {
      bool first=true;
      unsigned int count = 0;
      for(std::vector<Bucket>::const_iterator b = begin; b != end; ++b) {
        const BoundBox &box = b->boundingBox;
        count += b->mySize;
        if(b->mySize != 0){
          if(first){
            interval = box;
            first=false;
          }
          else
            interval.update(box);
        }
      }
      return count;
    }

    inline void BVHTree::order_elements(HandleDataVec::iterator &begin, 
                                        HandleDataVec::iterator &end, 
                                        SplitData &data) const  
    {
      for(HandleDataVec::iterator i = begin; i != end; ++i) 
      {
        const int index = Bucket::bucket_index(splitsPerDir, i->myBox, data.boundingBox, data.dim);
        i->myDim = (index<=data.split)?0:1;
      }
      std::sort(begin, end, HandleData_comparator());
    }

    inline void BVHTree::choose_best_split(const std::vector<std::vector<SplitData> > &splits,
                                           SplitData &data) const
    {
      std::vector<SplitData> best_splits;
      for(std::vector<std::vector<SplitData> >::const_iterator i = splits.begin(); i != splits.end(); ++i){
        std::vector<SplitData>::const_iterator j = std::min_element((*i).begin(), (*i).end(), 
                                                                    Split_comparator());
        best_splits.push_back(*j);
      }
      data = *std::min_element(best_splits.begin(), best_splits.end(), Split_comparator());
    }

    inline ErrorCode BVHTree::construct_element_vec(std::vector<HandleData> &handle_data_vec,
                                                    const Range &elements, 
                                                    BoundBox & bounding_box) 
    {
      std::vector<double> coordinate(3*CN::MAX_NODES_PER_ELEMENT);
      int num_conn;
      ErrorCode rval = MB_SUCCESS;
      std::vector<EntityHandle> storage;
      for(Range::iterator i = elements.begin(); i != elements.end(); ++i) {
          //TODO: not generic enough. Why dim != 3
          //Commence un-necessary deep copying.
        const EntityHandle* connect;
        rval = mbImpl->get_connectivity(*i, connect, num_conn, false, &storage);
        if (MB_SUCCESS != rval) return rval;
        rval = mbImpl->get_coords(connect, num_conn, &coordinate[0]);
        if (MB_SUCCESS != rval) return rval;
        BoundBox box;
        for(int j = 0; j < num_conn; j++)
          box.update(&coordinate[3*j]);
        if(i == elements.begin())
          bounding_box = box;
        else bounding_box.update(box);
        handle_data_vec.push_back(HandleData(*i, box, 0.0));
      }

      return rval;
    }

    inline ErrorCode BVHTree::reset_tree()
    {
      return delete_tree_sets();
    }

} // namespace moab

#endif //BVH_TREE_HPP
