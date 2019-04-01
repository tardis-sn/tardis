#ifndef MOAB_LINEAR_HEX_HPP
#define MOAB_LINEAR_HEX_HPP

#include "moab/Matrix3.hpp"
#include "moab/CartVect.hpp"
#include <sstream>
#include <iomanip>
#include <iostream>

namespace moab { 

namespace element_utility {

class Linear_hex_map {
  public: 
    //Constructor
    Linear_hex_map() {}
    //Copy constructor
    Linear_hex_map( const Linear_hex_map &) {}

 public:
    //Natural coordinates
    template< typename Moab, typename Entity_handle, 
	      typename Points, typename Point>
    std::pair< bool, CartVect> evaluate_reverse(const double *verts,
                                                const double *eval_point,
                                                const double tol=1.e-6) const{
      CartVect params(3, 0.0);
      solve_inverse(eval_point, params, verts);
      bool point_found = solve_inverse(eval_point, params, verts, tol) && 
          is_contained(params, tol);
      return std::make_pair(point_found, params);
    }

  private:
    //This is a hack to avoid a .cpp file and C++11
    //reference_points(i,j) will be a 1 or -1;
    //This should unroll..
  inline const CartVect &reference_points(const std::size_t& i) const{
    const CartVect rpts[8] = { CartVect( -1, -1, -1),
                              CartVect( 1, -1, -1),
                              CartVect( 1,  1, -1),
                              CartVect(-1,  1, -1),
                              CartVect(-1, -1,  1),
                              CartVect( 1, -1,  1),
                              CartVect( 1,  1,  1),
                              CartVect(-1,  1,  1)};
    return rpts[ i];
  }

  bool is_contained( const CartVect & params, const double tol) const{
     //just look at the box+tol here
    return ( params[0]>=-1.-tol) && (params[0]<=1.+tol) &&
        ( params[1]>=-1.-tol) && (params[1]<=1.+tol) &&
        ( params[2]>=-1.-tol) && (params[2]<=1.+tol);
  }

  bool solve_inverse(const CartVect &point, 
                     CartVect &params,
                     const CartVect *verts, 
                     const double tol=1.e-6) const {
    const double error_tol_sqr = tol*tol;
    CartVect delta(0.0, 0.0, 0.0);
    params = delta;
    evaluate_forward(params, verts, delta);
    delta -= point;
    std::size_t num_iterations=0;
#ifdef LINEAR_HEX_DEBUG
    std::stringstream ss;
    ss << "CartVect: "; 
    ss << point[0] << ", " << point[1] << ", " << point [2] << std::endl;
    ss << "Hex: ";
    for(int i = 0; i < 8; ++i)
      ss << points[i][0] << ", " << points[i][1] << ", " << points[i][2] << std::endl;
    ss << std::endl;
#endif
    while ( delta.length_squared() > error_tol_sqr) {
#ifdef LINEAR_HEX_DEBUG
      ss << "Iter #: "  << num_iterations 
         << " Err: " << delta.length() << " Iterate: ";
      ss << params[0] << ", " << params[1] << ", " << params[2] << std::endl;
#endif
      if( ++num_iterations >= 5){ return false; }
      Matrix3 J;
      jacobian( params, verts, J);
      double det = moab::Matrix3::determinant3(J);
      if (fabs(det) < 1.e-10){
#ifdef LINEAR_HEX_DEBUG
        std::cerr << ss.str();
#endif
#ifndef LINEAR_HEX_DEBUG
        std::cerr << x[0] << ", " << x[1] << ", " << x [2] << std::endl;
#endif
        std::cerr << "inverse solve failure: det: " << det << std::endl;
        exit(-1);
      }
      params -= moab::Matrix3::inverse(J, 1.0/det) * delta;
      evaluate_forward(params, points, delta);
      delta -= x;
    }
    return true;
  }

  void evaluate_forward(const CartVect &p, const CartVect *verts, CartVect &f) const{ 
    typedef typename Points::value_type Vector;
    f.set(0.0, 0.0, 0.0);
    for (unsigned i = 0; i < 8; ++i) {
      const double N_i = (1 + p[0]*reference_points(i)[0])
          * (1 + p[1]*reference_points(i)[1])
          * (1 + p[2]*reference_points(i)[2]);
      f += N_i * verts[ i];
    }
    f *= 0.125;
    return f;
  }

  double integrate_scalar_field(const CartVect *points, 
                                const double *field_values) const {
  }

  template< typename Point, typename Points>
  Matrix3& jacobian( const Point & p, const Points & points, Matrix3 & J) const{
  }
private:
}; //Class Linear_hex_map

}// namespace element_utility

}// namespace moab
#endif //MOAB_LINEAR_HEX_nPP
