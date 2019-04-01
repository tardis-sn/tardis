#ifndef MOAB_SPECTRAL_HEX_HPP
#define MOAB_SPECTRAL_HEX_HPP

#include "moab/Matrix3.hpp"
#include "moab/CartVect.hpp"
#include "moab/FindPtFuncs.h"
#include <sstream>
#include <iomanip>
#include <iostream>

namespace moab { 

namespace element_utility {

namespace {} //non-exported functionality

template< typename _Matrix>
class Spectral_hex_map {
  public:
	typedef _Matrix Matrix;
  private:
	typedef Spectral_hex_map< Matrix> Self;
  public: 
    //Constructor
    Spectral_hex_map() {};
    Spectral_hex_map( int order){ initialize_spectral_hex( order); }
    //Copy constructor
    Spectral_hex_map( const Self & f ) {}
  private:
    void initialize_spectral_hex( int order){
	if (_init && _n==order){ return; }
	if( _init && _n != order){ free_data();}
	_init = true;
	_n = order;
	for( int d = 0; d < 3; d++){
		lobatto_nodes(_z[ d], _n);
		lagrange_setup(&_ld[ d], _z[ d], _n);
	}
	opt_alloc_3(&_data, _ld);
	std::size_t nf = _n*_n, ne = _n, nw = 2*_n*_n + 3*_n;
	_odwork = tmalloc(real, 6*nf + 9*ne + nw);
    }

    void free_data(){
       for(int d=0; d<3; d++){
         free(_z[d]);
         lagrange_free(&_ld[d]);
       }
       opt_free_3(&_data);
       free(_odwork);
     }

 public:
    //Natural coordinates
    template< typename Moab, typename Entity_handle, 
	      typename Points, typename Point>
    std::pair< bool, Point> operator()( const Moab & /* moab */,
					const Entity_handle & /* h */,
					const Points & v, 
					const Point & p, 
					const double tol = 1.e-6) {
        Point result(3, 0.0);
        /*
        moab.tag_get_by_ptr(_xm1Tag, &eh, 1,(const void **) &_xyz[ 0] );
        moab.tag_get_by_ptr(_ym1Tag, &eh, 1,(const void **) &_xyz[ 1] );
        moab.tag_get_by_ptr(_zm1Tag, &eh, 1,(const void **) &_xyz[ 2] );
        */
        bool point_found = solve_inverse( p, result, v, tol) &&
                  is_contained( result, tol);
        return std::make_pair( point_found, result);
    }

  private:
    void set_gl_points( double * x, double * y, double * z){
	    _xyz[ 0] = x; _xyz[ 1] = y; _xyz[ 2] = z;
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
			const double tol=1.e-6) {
      const double error_tol_sqr = tol*tol;
      Point delta(3,0.0);
      xi = delta;
      evaluate( xi, points, delta);
      vec_subtract( delta, x);
      std::size_t num_iterations=0;
      #ifdef SPECTRAL_HEX_DEBUG
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
	#ifdef SPECTRAL_HEX_DEBUG
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
		#ifdef SPECTRAL_HEX_DEBUG
			std::cerr << ss.str();
		#endif
		#ifndef SPECTRAL_HEX_DEBUG
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
    Point& evaluate( const Point & p, const Points & /* points */, Point & f) {
      for (int d = 0; d < 3; ++d) { lagrange_0(&_ld[ d], p[ 0]); }
      for (int d = 0; d < 3; ++d) {
        f[ d] = tensor_i3( _ld[ 0].J, _ld[ 0].n,
        _ld[1].J, _ld[1].n,
        _ld[2].J, _ld[2].n,
        _xyz[ d],
        _odwork);
      }
      return f;
    }

   template< typename Point, typename Field>
   double   evaluate_scalar_field(const Point & p, const Field & field) const {
     int d;
     for(d=0; d<3; d++){ lagrange_0(&_ld[d], p[d]); }
     return tensor_i3( _ld[0].J,_ld[0].n,
           	       _ld[1].J,_ld[1].n,
                       _ld[2].J,_ld[2].n,
           	       field, _odwork);
   }
   template< typename Points, typename Field>
   double   integrate_scalar_field(const Points & p, 
				   const Field & field) const {  
   // set the position of GL points
   // set the positions of GL nodes, before evaluations
   _data.elx[0]=_xyz[0];
   _data.elx[1]=_xyz[1];
   _data.elx[2]=_xyz[2];
   double xi[3];
   //triple loop; the most inner loop is in r direction, then s, then t
   double integral = 0.;
   //double volume = 0;
   int index=0; // used fr the inner loop
   for (int k=0; k<_n; k++ ) {
     xi[2]=_ld[2].z[k];
     //double wk= _ld[2].w[k];
     for (int j=0; j<_n; j++) {
       xi[1]=_ld[1].z[j];
       //double wj= _ld[1].w[j];
       for (int i=0; i<_n; i++) {
         xi[0]=_ld[0].z[i];
         //double wi= _ld[0].w[i];
         opt_vol_set_intp_3((opt_data_3 *)&_data,xi); // cast away const-ness
         double wk= _ld[2].J[k];
         double wj= _ld[1].J[j];
         double wi= _ld[0].J[i];
         Matrix3 J(0.);
         for (int n = 0; n < 8; n++)
           J(n/3, n%3) = _data.jac[n];
         double bm = wk*wj*wi* J.determinant();
         integral+= bm*field[index++];
         //volume +=bm;
       }
     }
   }
   //std::cout << "volume: " << volume << "\n";
   return integral;
 }

  template< typename Point, typename Points>
  Matrix& jacobian( const Point & /* p */, const Points & /* points */, Matrix & J) {
   	real x[ 3];
    for (int i = 0; i < 3; ++i) { _data.elx[ i] = _xyz[ i]; }
    opt_vol_set_intp_3(& _data, x);
    for (int i = 0; i < 9; ++i) { J(i%3, i/3) = _data.jac[ i]; }
    return J;
  }

  private:
  bool _init;
  int _n;
  real * _z[ 3];
  lagrange_data _ld[ 3];
  opt_data_3 _data;
  real * _odwork;
  real * _xyz[ 3];
}; //Class Spectral_hex_map

}// namespace element_utility

}// namespace moab
#endif //MOAB_SPECTRAL_HEX_nPP
