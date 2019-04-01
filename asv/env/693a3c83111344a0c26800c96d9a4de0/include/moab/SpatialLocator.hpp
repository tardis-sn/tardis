/**\file SpatialLocator.hpp
 * \class moab::SpatialLocator
 * \brief Tool to facilitate spatial location of a point in a mesh
 *
 * SpatialLocator facilitates searching for points in or performing ray traces on collections of mesh entities
 * in 2D or 3D.  This searching is facilitated by a tree-based decomposition of the mesh.  Various types
 * of trees are implemented in MOAB and can be used by this tool, but by default it uses AdaptiveKDTree
 * (see child classes of Tree for which others are available).  Parallel and serial searching are both 
 * supported.
 *
 * SpatialLocator can either cache the search results for points or pass back this information in arguments.  
 * Cached information is kept in locTable, indexed in the same order as search points passed in.  This information
 * consists of the entity handle containing the point and the parametric coordinates inside that element.
 * Information about the points searched, e.g. the entities from which those points are derived, can be stored
 * in the calling application if desired.
 *
 * In parallel, there is a separation between the proc deciding which points to search for (the "target" proc), 
 * and the proc locating the point in its local mesh (the "source" proc).  On the source proc, location 
 * information is cached in locTable, as in the serial case.  By default, this location information (handle and
 * parametric coords) is not returned to the target proc, since it would be of no use there.  Instead, the rank
 * of the source proc locating the point, and the index of that location info in the source proc's locTable, is
 * returned; this information is stored on the target proc in this class's parLocTable variable.  Again, 
 * information about the points searched should be stored in the calling application, if desired.
 *
 * This class uses the ElemEvaluator class for specification and evaluation of basis functions (used for computing
 * parametric coords within an entity).  See documentation and examples for that class for usage information.
 * 
 */

#ifndef MOAB_SPATIALLOCATOR_HPP
#define MOAB_SPATIALLOCATOR_HPP

#include "moab/Types.hpp"
#include "moab/Tree.hpp"
#include "moab/Range.hpp"
#include "moab/TupleList.hpp"
#include "moab/BoundBox.hpp"
#include "moab/SpatialLocatorTimes.hpp"
#include "moab/CpuTimer.hpp"

#include <string>
#include <vector>
#include <map>
#include <math.h>

namespace moab {

    class Interface;
    class ElemEvaluator;
    class ParallelComm;

    class SpatialLocator
    {
  public:
        /* constructor */
      SpatialLocator(Interface *impl, Range &elems, Tree *tree = NULL, ElemEvaluator *eval = NULL);

        /* destructor */
      virtual ~SpatialLocator();

        /* add elements to be searched */
      ErrorCode add_elems(Range &elems);
      
        /* get bounding box of this locator */
      BoundBox &local_box() {return localBox;}
      
        /* get bounding box of this locator */
      const BoundBox &local_box() const {return localBox;}
      
        /* locate a set of vertices, Range variant */
      ErrorCode locate_points(Range &vertices,
                              EntityHandle *ents, double *params, int *is_inside = NULL,
                              const double rel_iter_tol = 1.0e-10, const double abs_iter_tol = 1.0e-10,
                              const double inside_tol = 1.0e-6);
      
        /* locate a set of points */
      ErrorCode locate_points(const double *pos, int num_points,
                              EntityHandle *ents, double *params, int *is_inside = NULL,
                              const double rel_iter_tol = 1.0e-10, const double abs_iter_tol = 1.0e-10,
                              const double inside_tol = 1.0e-6);
      
        /* locate a set of vertices or entity centroids, storing results on TupleList in this class
         * Locate a set of vertices or entity centroids, storing the detailed results in member 
         * variable (TupleList) locTable (see comments on locTable for structure of that tuple).
         */
      ErrorCode locate_points(Range &ents,
                              const double rel_iter_tol = 1.0e-10, const double abs_iter_tol = 1.0e-10,
                              const double inside_tol = 1.0e-6);
      
        /* locate a set of points, storing results on TupleList in this class
         * Locate a set of points, storing the detailed results in member variable (TupleList) locTable
         * (see comments on locTable for structure of that tuple).
         */
      ErrorCode locate_points(const double *pos, int num_points,
                              const double rel_iter_tol = 1.0e-10, const double abs_iter_tol = 1.0e-10,
                              const double inside_tol = 1.0e-6);

        /* Count the number of located points in locTable
         * Return the number of entries in locTable that have non-zero entity handles, which
         * represents the number of points in targetEnts that were inside one element in sourceEnts
         *
         */
      int local_num_located();

        /* Count the number of located points in parLocTable
         * Return the number of entries in parLocTable that have a non-negative index in on a remote
         * proc in parLocTable, which gives the number of points located in at least one element in a
         * remote proc's sourceEnts.
         */
      int remote_num_located();

#ifdef MOAB_HAVE_MPI      
        /* locate a set of vertices or entity centroids, storing results on TupleList in this class
         * Locate a set of vertices or entity centroids, storing the detailed results in member 
         * variables (TupleList) locTable and parLocTable (see comments on locTable and parLocTable for 
         * structure of those tuples).
         */
      ErrorCode par_locate_points(ParallelComm *pc,
                                  Range &vertices,
                                  const double rel_iter_tol = 1.0e-10, const double abs_iter_tol = 1.0e-10,
                                  const double inside_tol = 1.0e-6);
      
        /* locate a set of points, storing results on TupleList in this class
         * Locate a set of points, storing the detailed results in member 
         * variables (TupleList) locTable and parLocTable (see comments on locTable and parLocTable for 
         * structure of those tuples).
         */
      ErrorCode par_locate_points(ParallelComm *pc,
                                  const double *pos, int num_points,
                                  const double rel_iter_tol = 1.0e-10, const double abs_iter_tol = 1.0e-10,
                                  const double inside_tol = 1.0e-6);
#endif

        /** \brief Return the MOAB interface associated with this locator
         */
      Interface* moab() { return mbImpl; }

        /* locate a point */
      ErrorCode locate_point(const double *pos, 
                             EntityHandle &ent, double *params, int *is_inside = NULL,
                              const double rel_iter_tol = 1.0e-10, const double abs_iter_tol = 1.0e-10,
                              const double inside_tol = 1.0e-6);

        /* return the tree */
      Tree *get_tree() {return myTree;}

        /* get the locTable
         */
      TupleList &loc_table() {return locTable;}
      
        /* get the locTable
         */
      const TupleList &loc_table() const {return locTable;}
      
        /* get the parLocTable
         */
      TupleList &par_loc_table() {return parLocTable;}
      
        /* get the parLocTable
         */
      const TupleList &par_loc_table() const {return parLocTable;}

        /* get elemEval */
      ElemEvaluator *elem_eval() {return elemEval;}
      
        /* get elemEval */
      const ElemEvaluator *elem_eval() const {return elemEval;}
      
        /* set elemEval */
      void elem_eval(ElemEvaluator *eval) {elemEval = eval; if (myTree) myTree->set_eval(eval);}

        /** \brief Get spatial locator times object */
      SpatialLocatorTimes &sl_times() {return myTimes;}
      
        /** \brief Get spatial locator times object */
      const SpatialLocatorTimes &sl_times() const {return myTimes;}
        
  private:

#ifdef MOAB_HAVE_MPI
        /* MPI_ReduceAll source mesh bounding boxes to get global source mesh bounding box
         */
      ErrorCode initialize_intermediate_partition(ParallelComm *pc);
      
        /* for a given point in space, compute its ijk location in the intermediate decomposition; tolerance is
         * used only to test proximity to global box extent, not for local box test
         */
      inline ErrorCode get_point_ijk(const CartVect &point, const double abs_iter_tol, int *ijk) const;

#if 0
        /* given an ijk location in the intermediate partition, return the proc rank for that location 
         */
      inline int proc_from_ijk(const int *ijk) const;
#endif

        /* given a point in space, return the proc responsible for that point from the intermediate decomp; no tolerances
         * applied here, so first proc in lexicographic ijk ordering is returned
         */
      inline int proc_from_point(const double *pos, const double abs_iter_tol) const;
      
        /* register my source mesh with intermediate decomposition procs
         */
      ErrorCode register_src_with_intermediate_procs(ParallelComm *pc, double abs_iter_tol, TupleList &TLreg_o);
      
#endif

        /** Create a tree
         * Tree type depends on what's in myElems: if empty or all vertices, creates a kdtree,
         * otherwise creates a BVHTree.
         */
      void create_tree();
      
        /* MOAB instance */
      Interface* mbImpl;

        /* elements being located */
      Range myElems;

        /* dimension of entities in locator */
      int myDim;
      
        /* tree used for location */
      Tree *myTree;
      
        /* element evaluator */
      ElemEvaluator *elemEval;

        /* whether I created the tree or not (determines whether to delete it or not on destruction) */
      bool iCreatedTree;

        /* \brief local locations table
         * This table stores detailed local location data results from locate_points, that is, location data
         * for points located on the local mesh.  Data is stored
         * in a TupleList, where each tuple consists of (p_i, hs_ul, r[3]_d), where
         *   p_i = (int) proc from which request for this point was made (0 if serial)
         *   hs_ul = (unsigned long) source entity containing the point
         *   r[3]_d = (double) parametric coordinates of the point in hs 
         */
      TupleList locTable;

        /* \brief parallel locations table
         * This table stores information about points located on a local or remote processor.  For 
         * points located on this processor's local mesh, detailed location data is stored in locTable.
         * For points located on remote processors, more communication is required to retrieve specific
         * location data (usually that information isn't useful on this processor).
         *
         * The tuple structure of this TupleList is (p_i, ri_i), where:
         *   p_i = (int) processor rank containing this point
         *   ri_i = (int) index into locTable on remote proc containing this point's location information
         * The indexing of parLocTable corresponds to that of the points/entities passed in.
         */
      TupleList parLocTable;

        /* \brief Local bounding box, duplicated from tree
         */
      BoundBox localBox;

        /* \brief Global bounding box, used in parallel spatial location
         */
      BoundBox globalBox;

        /* \brief Regional delta xyz, used in parallel spatial location
         */
      CartVect regDeltaXYZ;

        /* \brief Number of regions in each of 3 directions
         */
      int regNums[3];

        /* \brief Map from source processor to bounding box of that proc's source mesh
         *
         */
      std::map<int, BoundBox> srcProcBoxes;

        /* \brief Timing object for spatial location
         */
      SpatialLocatorTimes myTimes;

        /* \brief Timer object to manage overloaded search functions
         */
      CpuTimer myTimer;
      
        /* \brief Flag to manage initialization of timer for overloaded search functions
         */
      bool timerInitialized;
    };

    inline SpatialLocator::~SpatialLocator() 
    {
      if (iCreatedTree && myTree) delete myTree;
    }
    
    inline ErrorCode SpatialLocator::locate_point(const double *pos, 
                                                  EntityHandle &ent, double *params, int *is_inside, 
                                                  const double rel_iter_tol, const double abs_iter_tol,
                                                  const double inside_tol)
    {
      return locate_points(pos, 1, &ent, params, is_inside, rel_iter_tol, abs_iter_tol, inside_tol);
    }

#ifdef MOAB_HAVE_MPI
    inline ErrorCode SpatialLocator::get_point_ijk(const CartVect &point, const double abs_iter_tol, int *ijk) const
    {
      for (int i = 0; i < 3; i++) {
        if (point[i] < globalBox.bMin[i]-abs_iter_tol || point[i] > globalBox.bMax[i]+abs_iter_tol)
          ijk[i] = -1;
        else {
          ijk[i] = point[i] - globalBox.bMin[i] / regDeltaXYZ[i];
          if (ijk[i] >= regNums[i] && point[i] <= globalBox.bMax[i]+abs_iter_tol)
            ijk[i] = regNums[i]-1;
        }
      }
      
      return (ijk[0] >= 0 && ijk[1] >= 0 && ijk[2] >= 0 ? MB_SUCCESS : MB_FAILURE);;
    }

#if 0
    inline int SpatialLocator::proc_from_ijk(const int *ijk) const
    {
      return ijk[2] * regNums[0]*regNums[1] + ijk[1] * regNums[0] + ijk[0];
    }
#endif

    inline int SpatialLocator::proc_from_point(const double *pos, const double abs_iter_tol) const
    {
      int ijk[3];
      ErrorCode rval = get_point_ijk(CartVect(pos), abs_iter_tol, ijk);
      if (MB_SUCCESS != rval) return -1;
      
      return ijk[2] * regNums[0]*regNums[1] + ijk[1] * regNums[0] + ijk[0];
    }
#endif
    
} // namespace moab 

#endif
