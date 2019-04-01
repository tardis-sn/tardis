#ifndef MOAB_LINEAR_TET_HPP
#define MOAB_LINEAR_TET_HPP

#include "moab/Matrix3.hpp"

namespace moab { 
namespace element_utility {

template<  typename Entity_handle, typename Matrix>
class Linear_tet_map {
  private:
	typedef Linear_tet_map< Entity_handle, Matrix> Self;
  public: 
    //Constructor
    Linear_tet_map() : Tinv(), eh() {}
    //Copy constructor
    Linear_tet_map( const Self & f ) : Tinv( f.Tinv), eh( f.eh){}
    //Natural coordinates
    template< typename Moab, typename Points, typename Point>
    std::pair< bool, Point> operator()( const Moab & moab,
					const Entity_handle _eh,
				        const Points & v, 
					const Point & p, 
					const double tol=1e-6) {
      // Remove the warning about unused parameter
      if (NULL != &moab) {}

      set_tet( _eh, v);
      //TODO: Make sure this is correct
      Point result = Tinv*p;
      return std::make_pair( is_contained( result, tol), result);
    }
  
    private:
    template< typename Point>
    bool is_contained( const Point & result, const double tol=1e-6){
	double sum=0.0;
	for( std::size_t i = 0; i < 3; ++i){ 
		sum += result[ i]; 
		if( result[ i] < -tol){ return false; } 
	}
	return sum < 1.0+tol;
    }
    template< typename Point, typename Field>
    double evaluate_scalar_field( const Point & p , 
				  const Field & field_values) const{
	double f0 = field_values[ 0];
	double f = f0;
	for(std::size_t i = 1; i < 5; ++i){f+=(field_values[ i] - f0)*p[ i -1];}
	return f;
    }
    template< typename Points, typename Field>
    double integrate_scalar_field(const Points & v, const Field & field_values) const {
      double I(0.0);
      for(unsigned int i = 0; i < 4; ++i) { I += field_values[i]; }
      double det = Matrix( v[1][0]-v[0][0], v[2][0]-v[0][0], 
        		   v[3][0]-v[0][0],
        		   v[1][1]-v[0][1], v[2][1]-v[0][1], 
        		   v[3][1]-v[0][1],
        		   v[1][2]-v[0][2], v[2][2]-v[0][2], 
        		   v[3][2]-v[0][2]).determinant();
      I *= det/24.0; 
      return I;
    }
    
    template< typename Points>
    void set_tet( const Entity_handle _eh, const Points & v){
		if (eh != _eh){
		   eh = _eh;
		   Tinv = moab::Matrix::inverse( 
		       Matrix( v[1][0]-v[0][0], v[2][0]-v[0][0], 
			       v[3][0]-v[0][0],
                    	       v[1][1]-v[0][1], v[2][1]-v[0][1], 
			       v[3][1]-v[0][1],
                    	       v[1][2]-v[0][2], v[2][2]-v[0][2], 
			       v[3][2]-v[0][2]) );
		}
    }
  private:
	Matrix Tinv;
	Entity_handle eh;
	/* We don't need this!
	static const double reference_points[ 4][ 3] = { {0,0,0},
				        	  	 {1,0,0},
				        	  	 {0,1,0},
				        	  	 {0,0,1} };*/
	
}; //Class Linear_tet_map

}// namespace element_utility

} // namespace moab
#endif //MOAB_LINEAR_TET_HPP
