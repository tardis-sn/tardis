#ifndef LINEAR_TRI_HPP
#define LINEAR_TRI_HPP
  /**\brief Shape function space for trilinear tetrahedron, obtained by a pushforward of the canonical linear (affine) functions. */

#include "moab/ElemEvaluator.hpp"
#include "moab/CN.hpp"

namespace moab 
{
    
class LinearTri 
{
public:
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode evalFcn(const double *params, const double *field, const int ndim, const int num_tuples, 
                           double *work, double *result);
        
    /** \brief Reverse-evaluation of parametric coordinates at physical space position */
  static ErrorCode reverseEvalFcn(EvalFcn eval, JacobianFcn jacob, InsideFcn ins, 
                                  const double *posn, const double *verts, const int nverts, const int ndim,
                                  const double iter_tol, const double inside_tol, double *work, 
                                  double *params, int *is_inside);

    /** \brief Evaluate the normal at a specified facet*/
  static ErrorCode normalFcn(const int ientDim, const int facet, const int nverts, const double *verts,  double normal[]);
        
    /** \brief Evaluate the jacobian at a specified parametric position */
  static ErrorCode jacobianFcn(const double *params, const double *verts, const int nverts, const int ndim, 
                               double *work, double *result);
        
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode integrateFcn(const double *field, const double *verts, const int nverts, const int ndim, const int num_tuples, 
                                double *work, double *result);

    /** \brief Initialize this EvalSet */
  static ErrorCode initFcn(const double *verts, const int nverts, double *&work);
      
        /** \brief Function that returns whether or not the parameters are inside the natural space of the element */
  static int insideFcn(const double *params, const int ndim, const double tol);
  
  static ErrorCode evaluate_reverse(EvalFcn eval, JacobianFcn jacob, InsideFcn inside_f,
                                    const double *posn, const double *verts, const int nverts, 
                                    const int ndim, const double iter_tol, const double inside_tol, double *work, 
                                    double *params, int *inside);

  static EvalSet eval_set() 
      {
        return EvalSet(evalFcn, reverseEvalFcn, normalFcn, jacobianFcn, integrateFcn, initFcn, insideFcn);
      }
      
  static bool compatible(EntityType tp, int numv, EvalSet &eset) 
      {
        if (tp == MBTRI && numv >= 3) {
          eset = eval_set();
          return true;
        }
        else return false;
      }
  

protected:
      
  static const double corner[3][2];
};// class LinearTri

} // namespace moab

#endif
