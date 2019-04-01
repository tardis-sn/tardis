#ifndef MOAB_QUADRATIC_HEX_HPP
#define MOAB_QUADRATIC_HEX_HPP

#include "moab/Matrix3.hpp"
#include "moab/CartVect.hpp"
#include <sstream>
#include <iomanip>
#include <iostream>

namespace moab { 

namespace element_utility {

namespace {

double SH(const int i, const double xi)
{
  switch (i)
  {
  case -1: return (xi*xi-xi)/2;
  case 0: return 1-xi*xi;
  case 1: return (xi*xi+xi)/2;
  default: return 0.;
  }
}
double DSH(const int i, const double xi)
{
  switch (i)
  {
  case -1: return xi-0.5;
  case 0: return -2*xi;
  case 1: return xi+0.5;
  default: return 0.;
  }
}


} //non-exported functionality

template< typename _Matrix>
class Quadratic_hex_map {
  public:
	typedef _Matrix Matrix;
  private:
	typedef Quadratic_hex_map< Matrix> Self;
  public: 
    //Constructor
    Quadratic_hex_map() {}
    //Copy constructor
    Quadratic_hex_map( const Self & f ) {}

 public:
    //Natural coordinates
    template< typename Moab, typename Entity_handle, 
	      typename Points, typename Point>
    std::pair< bool, Point> operator()( const Moab & /* moab */,
					const Entity_handle & /* h */,
					const Points & v, 
					const Point & p, 
					const double tol = 1.e-6) const {
      Point result(3, 0.0);
      bool point_found = solve_inverse( p, result, v, tol) &&
                is_contained( result, tol);
      return std::make_pair( point_found, result);
    }

  private:
    //This is a hack to avoid a .cpp file and C++11
    //reference_points(i,j) will be a 1 or -1;
    //This should unroll..
    inline double reference_points( const std::size_t& i,
          				  const std::size_t& j) const{
    const double rpts[27][3] = {
    	{ -1, -1, -1 },
    	{  1, -1, -1 },
    	{  1,  1, -1 },  // reference_points nodes: 0-7
    	{ -1,  1, -1 },  // mid-edge nodes: 8-19
    	{ -1, -1,  1 },  // center-face nodes 20-25  center node  26
    	{  1, -1,  1 },  //
    	{  1,  1,  1 },
    	{ -1,  1,  1 }, //                    4   ----- 19   -----  7
    	{  0, -1, -1 }, //                .   |                 .   |
    	{  1,  0, -1 }, //            16         25         18      |
    	{  0,  1, -1 }, //         .          |          .          |
    	{ -1,  0, -1 }, //      5   ----- 17   -----  6             |
    	{ -1, -1,  0 }, //      |            12       | 23         15
    	{  1, -1,  0 }, //      |                     |             |
    	{  1,  1,  0 }, //      |     20      |  26   |     22      |
    	{ -1,  1,  0 }, //      |                     |             |
    	{  0, -1,  1 }, //     13         21  |      14             |
    	{  1,  0,  1 }, //      |             0   ----- 11   -----  3
    	{  0,  1,  1 }, //      |         .           |         .
    	{ -1,  0,  1 }, //      |      8         24   |     10
    	{  0, -1,  0 }, //      |  .                  |  .
    	{  1,  0,  0 }, //      1   -----  9   -----  2
    	{  0,  1,  0 }, //
    	{ -1,  0,  0 },
    	{  0,  0, -1 },
    	{  0,  0,  1 },
    	{  0,  0,  0 }
	};
	  return rpts[ i][ j];
    }

    template< typename Point>
    bool is_contained( const Point & p, const double tol) const{
     //just look at the box+tol here
     return ( p[0]>=-1.-tol) && (p[0]<=1.+tol) &&
            ( p[1]>=-1.-tol) && (p[1]<=1.+tol) &&
            ( p[2]>=-1.-tol) && (p[2]<=1.+tol);
    }

    template< typename Point, typename Points>
    bool solve_inverse( const Point & x, 
			Point & xi,
			const Points & points, 
			const double tol=1.e-6) const {
      const double error_tol_sqr = tol*tol;
      Point delta(3,0.0);
      xi = delta;
      evaluate( xi, points, delta);
      vec_subtract( delta, x);
      std::size_t num_iterations=0;
      #ifdef QUADRATIC_HEX_DEBUG
 	std::stringstream ss;
	ss << "Point: "; 
       ss << x[ 0 ] << ", " << x[ 1] 
          << ", " << x [ 2] << std::endl;
	ss << "Hex: ";
	for(int i = 0; i < 8; ++i){
 	      	ss << points[ i][ 0] << ", " << points[ i][ 1] << ", "
		   << points[ i][ 2] << std::endl;
	}
	ss << std::endl;
      #endif
      while ( normsq( delta) > error_tol_sqr) {
	#ifdef QUADRATIC_HEX_DEBUG
	ss << "Iter #: "  << num_iterations 
	   << " Err: " << sqrt( normsq( delta)) << " Iterate: ";
	ss << xi[ 0 ] << ", " << xi[ 1] 
		<< ", " << xi[ 2] << std::endl;
	#endif
	if( ++num_iterations >= 5){ return false; }
        Matrix J;
	jacobian( xi, points, J);
        double det = moab::Matrix::determinant3( J);
        if (fabs(det) < 1.e-10){
		#ifdef QUADRATIC_HEX_DEBUG
			std::cerr << ss.str();
		#endif
		#ifndef QUADRATIC_HEX_DEBUG
		std::cerr << x[ 0 ] << ", " << x[ 1] 
			  << ", " << x [ 2] << std::endl;
		#endif
		std::cerr << "inverse solve failure: det: " << det << std::endl;
		exit( -1);
	}
        vec_subtract( xi, moab::Matrix::inverse(J, 1.0/det) * delta);
        evaluate( xi, points, delta);
	vec_subtract( delta, x);
      }
       return true;
    }

    template< typename Point, typename Points>
    Point& evaluate( const Point & p, const Points & points, Point & f) const{ 
	typedef typename Points::value_type Vector;
	Vector result;
	for(int i = 0; i < 3; ++i){ result[ i] = 0; }
	for (unsigned i = 0; i < 27; ++i) {
	    const double sh= SH(reference_points(i,0), p[0])*
       		             SH(reference_points(i,1), p[1])*
                	     SH(reference_points(i,2), p[2]);
	    result += sh * points[ i];
	}
	for (int i = 0; i < 3; ++i){ f[ i] = result[ i]; }
	return f;
    }
    template< typename Point, typename Field>
    double   evaluate_scalar_field( const Point & p, 
				    const Field & field) const {
      double x=0.0;
      for (int i=0; i<27; i++){
        const double sh= SH(reference_points(i,0), p[0])*
                  	 SH(reference_points(i,1), p[1])*
                  	 SH(reference_points(i,2), p[2]);
        x+=sh* field[i];
      }
      return x;
    }
    template< typename Field, typename Points>
    double   integrate_scalar_field( const Points & p, 
				     const Field & field_values) const { 
        // TODO: gaussian integration , probably 2x2x2
	return 0.; 
    }

    template< typename Point, typename Points>
    Matrix& jacobian( const Point & p, const Points & /* points */, Matrix & J) const {
    J = Matrix(0.0);
    for (int i = 0; i < 27; i++) {
      const double sh[3] = { SH(reference_points(i,0), p[0]),
                       		           SH(reference_points(i,1), p[1]),
                       		           SH(reference_points(i,2), p[2]) };
 	    const double dsh[3] = { DSH(reference_points(i,0), p[0]),
                        		  DSH(reference_points(i,1), p[1]),
                        	   	  DSH(reference_points(i,2), p[2]) };
      for (int j = 0; j < 3; j++) {
        // dxj/dr first column
        J(j, 0) += dsh[0]*sh[1]*sh[2]*reference_points(i, j);
        J(j, 1) += sh[0]*dsh[1]*sh[2]*reference_points(i, j); // dxj/ds
        J(j, 2) += sh[0]*sh[1]*dsh[2]*reference_points(i, j); // dxj/dt
      }
    }
    return J;
  }
  private:
}; //Class Quadratic_hex_map

}// namespace element_utility

}// namespace moab
#endif //MOAB_QUADRATIC_HEX_nPP
