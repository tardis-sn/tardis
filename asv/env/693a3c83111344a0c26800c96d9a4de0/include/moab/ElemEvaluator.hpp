#ifndef ELEM_EVALUATOR_HPP
#define ELEM_EVALUATOR_HPP

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/Matrix3.hpp"
#include "moab/CN.hpp"

#include <vector>

namespace moab {

    typedef ErrorCode (*EvalFcn)(const double *params, const double *field, const int ndim, const int num_tuples, 
                          double *work, double *result);

    typedef ErrorCode (*JacobianFcn)(const double *params, const double *verts, const int nverts, const int ndim, 
                                     double *work, double *result);
        
    typedef ErrorCode (*IntegrateFcn)(const double *field, const double *verts, const int nverts, const int ndim,
                                      const int num_tuples, double *work, double *result);

    typedef ErrorCode (*InitFcn)(const double *verts, const int nverts, double *&work);

    typedef int (*InsideFcn)(const double *verts, const int ndims, const double tol);

    typedef ErrorCode (*ReverseEvalFcn)(EvalFcn eval, JacobianFcn jacob, InsideFcn ins, 
                                        const double *posn, const double *verts, const int nverts, const int ndim,
                                        const double iter_tol, const double inside_tol, 
                                        double *work, double *params, int *is_inside);
        
  typedef ErrorCode (*NormalFcn)(const int ientDim, const int facet, const int nverts, const double *verts,  double normal[3]);

    class EvalSet
    {
  public:
        /** \brief Forward-evaluation of field at parametric coordinates */
      EvalFcn evalFcn;
        
        /** \brief Reverse-evaluation of parametric coordinates at physical space position */
      ReverseEvalFcn reverseEvalFcn;
        
      /** \brief Evaluate the normal at a local facet (edge/face for 2D/3D) */
      NormalFcn normalFcn;

        /** \brief Evaluate the jacobian at a specified parametric position */
      JacobianFcn jacobianFcn;
        
        /** \brief Forward-evaluation of field at parametric coordinates */
      IntegrateFcn integrateFcn;

        /** \brief Initialization function for an element */
      InitFcn initFcn;

        /** \brief Function that returns whether or not the parameters are inside the natural space of the element */
      InsideFcn insideFcn;

        /** \brief Bare constructor */
      EvalSet() : evalFcn(NULL), reverseEvalFcn(NULL), normalFcn(NULL), jacobianFcn(NULL), integrateFcn(NULL), initFcn(NULL), insideFcn(NULL) {}

        /** \brief Constructor */
      EvalSet(EvalFcn eval, ReverseEvalFcn rev, NormalFcn normal, JacobianFcn jacob, IntegrateFcn integ, InitFcn initf, InsideFcn insidef)
              : evalFcn(eval), reverseEvalFcn(rev), normalFcn(normal), jacobianFcn(jacob), integrateFcn(integ), initFcn(initf), insideFcn(insidef)
          {}

        /** \brief Given an entity handle, get an appropriate eval set, based on type & #vertices */
      static ErrorCode get_eval_set(Interface *mb, EntityHandle eh, EvalSet &eval_set);
      
        /** \brief Given type & #vertices, get an appropriate eval set */
      static ErrorCode get_eval_set(EntityType tp, unsigned int num_vertices, EvalSet &eval_set);
      
        /** \brief Operator= */
      EvalSet &operator=(const EvalSet &eval) {
        evalFcn = eval.evalFcn;
        reverseEvalFcn = eval.reverseEvalFcn;
        normalFcn = eval.normalFcn;
        jacobianFcn = eval.jacobianFcn;
        integrateFcn = eval.integrateFcn;
        initFcn = eval.initFcn;
        insideFcn = eval.insideFcn;
        return *this;
      }

        /** \brief Common function to do reverse evaluation based on evaluation and jacobian functions */
      static ErrorCode evaluate_reverse(EvalFcn eval, JacobianFcn jacob, InsideFcn inside_f,
                                        const double *posn, const double *verts, const int nverts, 
                                        const int ndim, const double iter_tol, const double inside_tol, 
                                        double *work, double *params, int *inside);
        /** \brief Common function that returns true if params is in [-1,1]^ndims */
      static int inside_function(const double *params, const int ndims, const double tol);
    };

        /** \brief Given an entity handle, get an appropriate eval set, based on type & #vertices */
    inline ErrorCode EvalSet::get_eval_set(Interface *mb, EntityHandle eh, EvalSet &eval_set) 
    {
      int nv;
      EntityType tp = mb->type_from_handle(eh);
      const EntityHandle *connect;
      std::vector<EntityHandle> dum_vec;
      ErrorCode rval = mb->get_connectivity(eh, connect, nv, false, &dum_vec);
      if (MB_SUCCESS != rval) return rval;
      
      return get_eval_set(tp, nv, eval_set);
    }
      
/**\brief Class facilitating local discretization-related functions
 * \class ElemEvaluator
 * This class implements discretization-related functionality operating
 * on data in MOAB.  A member of this class caches certain data about the element
 * it's currently operating on, but is not meant to be instantiated one-per-element,
 * but rather one-per-search (or other operation on a collection of elements).
 *
 * Actual discretization functionality is accessed through function pointers,
 * allowing applications to specialize the implementation of specific functions
 * while still using this class.
 *
 * This class depends on MOAB functionality for gathering entity-based data; the functions
 * it calls through function pointers depend only on POD (plain old data, or intrinsic data types).
 * This allows the use of other packages for serving these functions, without having to modify
 * them to get data through MOAB.  This should also promote efficiency, since in many cases they
 * will be able to read data from its native storage locations.
 */

    class ElemEvaluator {
  public:
        /** \brief Constructor 
         * \param impl MOAB instance
         * \param ent Entity handle to cache on the evaluator; if non-zero, calls set_ent_handle, which does some other stuff.
         * \param tag Tag to cache on the evaluator; if non-zero, calls set_tag_handle, which does some other stuff too.
         * \param tagged_ent_dim Dimension of entities to be tagged to cache on the evaluator
         */
      ElemEvaluator(Interface *impl, EntityHandle ent = 0, Tag tag = 0, int tagged_ent_dim = -1);
      ~ElemEvaluator();

        /** \brief Evaluate cached tag at a given parametric location within the cached entity 
         * If evaluating coordinates, call set_tag(0, 0), which indicates coords instead of a tag.
         * \param params Parameters at which to evaluate field
         * \param result Result of evaluation
         * \param num_tuples Size of evaluated field, in number of values
         */
      ErrorCode eval(const double *params, double *result, int num_tuples = -1) const;
        
        /** \brief Reverse-evaluate the cached entity at a given physical position
         * \param posn Position at which to evaluate parameters
         * \param iter_tol Tolerance of reverse evaluation non-linear iteration, usually 10^-10 or so
         * \param inside_tol Tolerance of is_inside evaluation, usually 10^-6 or so
         * \param params Result of evaluation
         * \param is_inside If non-NULL, returns true of resulting parameters place the point inside the element
         *                  (in most cases, within [-1]*(dim)
         */
      ErrorCode reverse_eval(const double *posn, double iter_tol, double inside_tol, double *params, 
                             int *is_inside = NULL) const;

      /**
       * \brief Evaluate the normal to a facet of an entity
       * \param ientDim Dimension of the facet. Should be (d-1) for d-dimensional entities
       * \param facet Local id of the facet w.r.t the entity
       * \param normal Returns the normal.
       */
      ErrorCode get_normal(const int ientDim, const int facet, double normal[3]) const;
        
        /** \brief Evaluate the jacobian of the cached entity at a given parametric location
         * \param params Parameters at which to evaluate jacobian
         * \param result Result of evaluation, in the form of a 3x3 matrix, stored in column-major order
         */
      ErrorCode jacobian(const double *params, double *result) const;
        
        /** \brief Integrate the cached tag over the cached entity
         * \param result Result of the integration
         */
      ErrorCode integrate(double *result) const;
      
        /** \brief Return whether a physical position is inside the cached entity to within a tolerance
         * \param params Parameters at which to query the element
         * \param tol Tolerance, usually 10^-6 or so
         */
      int inside(const double *params, const double tol) const;

        /** \brief Given a list of entities, return the entity the point is in, or none
         * This function reverse-evaluates the entities, returning the first entity containing the point.
         * If no entity contains the point, containing_ent is returned as 0 and params are unchanged.
         * This function returns something other than MB_SUCCESS only when queries go wrong for some reason.
         * num_evals, if non-NULL, is always incremented for each call to reverse_eval.
         * This function calls set_ent_handle for each entity before calling reverse_eval, so the ElemEvaluator
         * object is changed.
         * \param entities Entities tested
         * \param point Point tested, must have 3 dimensions, even for edge and face entities
         * \param iter_tol Tolerance for non-linear reverse evaluation
         * \param inside_tol Tolerance for is_inside test
         * \param containing_ent Entity containing the point, returned 0 if no entity
         * \param params Parameters of point in containing entity, unchanged if no containing entity
         * \param num_evals If non-NULL, incremented each time reverse_eval is called
         * \return Returns non-success only if evaulation failed for some reason (point not in element is NOT a
         * reason for failure)
         */
      ErrorCode find_containing_entity(Range &entities, const double *point, 
                                       const double iter_tol, const double inside_tol, 
                                       EntityHandle &containing_ent, double *params, 
                                       unsigned int *num_evals = NULL);
      
        /** \brief Given an entity set, return the contained entity the point is in, or none
         * This function reverse-evaluates the entities, returning the first entity containing the point.
         * If no entity contains the point, containing_ent is returned as 0 and params are unchanged.
         * This function returns something other than MB_SUCCESS only when queries go wrong for some reason.
         * num_evals, if non-NULL, is always incremented for each call to reverse_eval.
         * This function calls set_ent_handle for each entity before calling reverse_eval, so the ElemEvaluator
         * object is changed.
         * \param ent_set Entity set containing the entities to be tested
         * \param point Point tested, must have 3 dimensions, even for edge and face entities
         * \param iter_tol Tolerance for non-linear reverse evaluation
         * \param inside_tol Tolerance for is_inside test
         * \param containing_ent Entity containing the point, returned 0 if no entity
         * \param params Parameters of point in containing entity, unchanged if no containing entity
         * \param num_evals If non-NULL, incremented each time reverse_eval is called
         * \return Returns non-success only if evaulation failed for some reason (point not in element is NOT a
         * reason for failure)
         */
      ErrorCode find_containing_entity(EntityHandle ent_set, const double *point, 
                                       const double iter_tol, const double inside_tol,
                                       EntityHandle &containing_ent, double *params, 
                                       unsigned int *num_evals = NULL);
      
        /** \brief Set the eval set for a given type entity
         * \param tp Entity type for which to set the eval set
         * \param eval_set Eval set object to use
         */
      ErrorCode set_eval_set(EntityType tp, const EvalSet &eval_set);
      
        /** \brief Set the eval set using a given entity handle
        * This function queries the entity type and number of vertices on the entity to decide which type
        * of shape function to use.
        * NOTE: this function should only be called after setting a valid MOAB instance on the evaluator
        * \param eh Entity handle whose type and #vertices are queried
        */
      ErrorCode set_eval_set(const EntityHandle eh);
      
        /** \brief Get the eval set for a given type entity */
      inline EvalSet get_eval_set(EntityType tp) {return evalSets[tp];}

        /** \brief Set entity handle & cache connectivty & vertex positions */
      inline ErrorCode set_ent_handle(EntityHandle ent);

        /** \brief Get entity handle for this ElemEval */
      inline EntityHandle get_ent_handle() const {return entHandle;}

        /* \brief Get vertex positions cached on this EE
         */
      inline double *get_vert_pos() {return vertPos[0].array();}
      
        /* \brief Get the vertex handles cached here */
      inline const EntityHandle *get_vert_handles() const {return vertHandles;}

        /* \brief Get the number of vertices for the cached entity */
      inline int get_num_verts() const {return numVerts;}

        /* \brief Get the tag handle cached on this entity */
      inline Tag get_tag_handle() const {return tagHandle;};

        /* \brief Set the tag handle to cache on this evaluator
         * To designate that vertex coordinates are the desired tag, pass in a tag handle of 0
         * and a tag dimension of 0.
         * \param tag Tag handle to cache, or 0 to cache vertex positions
         * \param tagged_ent_dim Dimension of entities tagged with this tag
         */
      inline ErrorCode set_tag_handle(Tag tag, int tagged_ent_dim = -1);

        /* \brief Set the name of the tag to cache on this evaluator
         * To designate that vertex coordinates are the desired tag, pass in "COORDS" as the name
         * and a tag dimension of 0.
         * \param tag_name Tag handle to cache, or 0 to cache vertex positions
         * \param tagged_ent_dim Dimension of entities tagged with this tag
         */
      inline ErrorCode set_tag(const char *tag_name, int tagged_ent_dim = -1);
      
        /* \brief Get the dimension of the entities on which tag is set */
      inline int get_tagged_ent_dim() const {return taggedEntDim;};

        /* \brief Set the dimension of entities having the tag */
      inline ErrorCode set_tagged_ent_dim(int dim);

        /* \brief Get work space, sometimes this is useful for evaluating data you don't want to set as tags on entities
         * Can't be const because most of the functions (evaluate, integrate, etc.) take non-const work space *'s
         */
      inline double *get_work_space() {return workSpace;}
            
        /* \brief MOAB interface cached on this evaluator */
      inline Interface *get_moab() {return mbImpl;}
      
  private:

        /** \brief Interface */
      Interface *mbImpl;
      
        /** \brief Entity handle being evaluated */
      EntityHandle entHandle;

        /** \brief Entity type */
      EntityType entType;

        /** \brief Entity dimension */
      int entDim;
            
        /** \brief Number of vertices cached here */
      int numVerts;

        /** \brief Cached copy of vertex handle ptr */
      const EntityHandle *vertHandles;
      
        /** \brief Cached copy of vertex positions */
      CartVect vertPos[CN::MAX_NODES_PER_ELEMENT];
      
        /** \brief Tag being evaluated */
      Tag tagHandle;

        /** \brief Whether tag is coordinates or something else */
      bool tagCoords;
      
        /** \brief Number of values in this tag */
      int numTuples;

        /** \brief Dimension of entities from which to grab tag */
      int taggedEntDim;

        /** \brief Tag space */
      std::vector<unsigned char> tagSpace;
      
        /** \brief Evaluation methods for elements of various topologies */
      EvalSet evalSets[MBMAXTYPE];

        /** \brief Work space for element-specific data */
      double *workSpace;

    }; // class ElemEvaluator

    inline ElemEvaluator::ElemEvaluator(Interface *impl, EntityHandle ent, Tag tag, int tagged_ent_dim) 
            : mbImpl(impl), entHandle(0), entType(MBMAXTYPE), entDim(-1), numVerts(0), 
              vertHandles(NULL), tagHandle(0), tagCoords(false), numTuples(0), 
              taggedEntDim(0), workSpace(NULL)
    {
      if (ent) set_ent_handle(ent);
      if (tag) set_tag_handle(tag, tagged_ent_dim);
    }

    inline ElemEvaluator::~ElemEvaluator()
    {
      if (workSpace)
        delete [] workSpace;
    }

    inline ErrorCode ElemEvaluator::set_ent_handle(EntityHandle ent) 
    {
      entHandle = ent;
      if (workSpace) {
        delete [] workSpace;
        workSpace = NULL;
      }

      entType = mbImpl->type_from_handle(ent);
      entDim = mbImpl->dimension_from_handle(ent);

      std::vector<EntityHandle> dum_vec;
      ErrorCode rval = mbImpl->get_connectivity(ent, vertHandles, numVerts, false, &dum_vec);
      if (MB_SUCCESS != rval) return rval;

      if (!evalSets[entType].evalFcn)
        EvalSet::get_eval_set(entType, numVerts, evalSets[entType]);

      rval = mbImpl->get_coords(vertHandles, numVerts, vertPos[0].array());
      if (MB_SUCCESS != rval) return rval;

      if (tagHandle) {
        rval = set_tag_handle(tagHandle);
        if (MB_SUCCESS != rval) return rval;
      }
      if (evalSets[entType].initFcn) return (*evalSets[entType].initFcn)(vertPos[0].array(), numVerts, workSpace);
      return MB_SUCCESS;
    }
    
    inline ErrorCode ElemEvaluator::set_tag_handle(Tag tag, int tagged_ent_dim) 
    {
      ErrorCode rval = MB_SUCCESS;
      if (!tag && !tagged_ent_dim) {
        tagCoords = true;
        numTuples = 3;
        taggedEntDim = 0;
        tagHandle = 0;
        return rval;
      }
      else if (tagHandle != tag) {
        tagHandle = tag;
        rval = mbImpl->tag_get_length(tagHandle, numTuples);
        if (MB_SUCCESS != rval) return rval;
        int sz;
        rval = mbImpl->tag_get_bytes(tag, sz);
        if (MB_SUCCESS != rval) return rval;
        tagSpace.resize(CN::MAX_NODES_PER_ELEMENT*sz);
        tagCoords = false;
      }

      taggedEntDim = (-1 == tagged_ent_dim ? 0 : tagged_ent_dim);
      
      if (entHandle) {
        if (0 == taggedEntDim) {
          rval = mbImpl->tag_get_data(tagHandle, vertHandles, numVerts, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
        else if (taggedEntDim == entDim) {
          rval = mbImpl->tag_get_data(tagHandle, &entHandle, 1, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
      }

      return rval;
    }

    inline ErrorCode ElemEvaluator::set_tag(const char *tag_name, int tagged_ent_dim) 
    {
      ErrorCode rval = MB_SUCCESS;
      if (!tag_name || strlen(tag_name) == 0) return MB_FAILURE;
      Tag tag;
      if (!strcmp(tag_name, "COORDS")) {
        tagCoords = true;
        taggedEntDim = 0;
        numTuples = 3;
        tagHandle = 0;
          // can return here, because vertex coords already cached when entity handle set
        return rval;
      }
      else {
        rval = mbImpl->tag_get_handle(tag_name, tag);
        if (MB_SUCCESS != rval) return rval;
      
        if (tagHandle != tag) {
          tagHandle = tag;
          rval = mbImpl->tag_get_length(tagHandle, numTuples);
          if (MB_SUCCESS != rval) return rval;
          int sz;
          rval = mbImpl->tag_get_bytes(tag, sz);
          if (MB_SUCCESS != rval) return rval;
          tagSpace.reserve(CN::MAX_NODES_PER_ELEMENT*sz);
          tagCoords = false;
        }

        taggedEntDim = (-1 == tagged_ent_dim ? entDim : tagged_ent_dim);
      }
      
      if (entHandle) {
        if (0 == taggedEntDim) {
          rval = mbImpl->tag_get_data(tagHandle, vertHandles, numVerts, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
        else if (taggedEntDim == entDim) {
          rval = mbImpl->tag_get_data(tagHandle, &entHandle, 1, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
      }

      return rval;
    }

    inline ErrorCode ElemEvaluator::set_eval_set(EntityType tp, const EvalSet &eval_set) 
    {
      evalSets[tp] = eval_set;
      if (entHandle && evalSets[entType].initFcn) {
        ErrorCode rval = (*evalSets[entType].initFcn)(vertPos[0].array(), numVerts, workSpace);
        if (MB_SUCCESS != rval) return rval;
      }
      return MB_SUCCESS;
    }
    
    inline ErrorCode ElemEvaluator::eval(const double *params, double *result, int num_tuples) const 
    {
      assert(entHandle && MBMAXTYPE != entType);
      return (*evalSets[entType].evalFcn)(params, 
                                          (tagCoords ? (const double*) vertPos[0].array(): (const double*)&tagSpace[0]), 
                                          entDim, 
                                          (-1 == num_tuples ? numTuples : num_tuples), workSpace, result);
    }
        
    inline ErrorCode ElemEvaluator::reverse_eval(const double *posn, const double iter_tol, const double inside_tol,
                                                 double *params, int *ins) const
    {
      assert(entHandle && MBMAXTYPE != entType);
      return (*evalSets[entType].reverseEvalFcn)(evalSets[entType].evalFcn, evalSets[entType].jacobianFcn, evalSets[entType].insideFcn,
                                                 posn, vertPos[0].array(), numVerts, 
                                                 entDim, iter_tol, inside_tol, workSpace, params, ins);
    }

      /** \brief Evaluate the normal of the cached entity at a given facet */
    inline ErrorCode ElemEvaluator::get_normal(const int ientDim, const int facet, double normal[]) const
    {
      assert(entHandle && MBMAXTYPE != entType);
      return (*evalSets[entType].normalFcn)( ientDim, facet, numVerts, vertPos[0].array(), normal);
    }
        
      /** \brief Evaluate the jacobian of the cached entity at a given parametric location */
    inline ErrorCode ElemEvaluator::jacobian(const double *params, double *result) const
    {
      assert(entHandle && MBMAXTYPE != entType);
      return (*evalSets[entType].jacobianFcn)(params, vertPos[0].array(), numVerts, entDim, workSpace, result);
    }
        
      /** \brief Integrate the cached tag over the cached entity */
    inline ErrorCode ElemEvaluator::integrate(double *result) const
    {
      assert(entHandle && MBMAXTYPE != entType && (tagCoords || tagHandle));
      ErrorCode rval = MB_SUCCESS;
      if (!tagCoords) {
        if (0 == taggedEntDim) rval = mbImpl->tag_get_data(tagHandle, vertHandles, numVerts, (void*)&tagSpace[0]);
        else rval = mbImpl->tag_get_data(tagHandle, &entHandle, 1, (void*)&tagSpace[0]);
        if (MB_SUCCESS != rval) return rval;
      }
      return (*evalSets[entType].integrateFcn)((tagCoords ? vertPos[0].array() : (const double *)&tagSpace[0]), 
                                               vertPos[0].array(), numVerts, entDim, numTuples, 
                                               workSpace, result);
    }

    inline ErrorCode ElemEvaluator::find_containing_entity(EntityHandle ent_set, const double *point, 
                                                           const double iter_tol, const double inside_tol,
                                                           EntityHandle &containing_ent, double *params, 
                                                           unsigned int *num_evals) 
    {
      assert(mbImpl->type_from_handle(ent_set) == MBENTITYSET);
      Range entities;
      ErrorCode rval = mbImpl->get_entities_by_handle(ent_set, entities);
      if (MB_SUCCESS != rval) return rval;
      else return find_containing_entity(entities, point, iter_tol, inside_tol, containing_ent, params, num_evals);
    }
        
    inline int ElemEvaluator::inside(const double *params, const double tol) const 
    {
      return (*evalSets[entType].insideFcn)(params, entDim, tol);
    }

    inline ErrorCode ElemEvaluator::set_eval_set(const EntityHandle eh) 
    {
      EvalSet eset;
      ErrorCode rval = EvalSet::get_eval_set(mbImpl, eh, evalSets[mbImpl->type_from_handle(eh)]); 
      return rval;
    }

} // namespace moab

#endif /*MOAB_ELEM_EVALUATOR_HPP*/
