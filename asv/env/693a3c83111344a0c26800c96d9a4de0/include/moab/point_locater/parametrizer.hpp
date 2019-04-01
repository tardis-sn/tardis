#ifndef MOAB_PARAMETRIZER_HPP
#define MOAB_PARAMETRIZER_HPP
#include "moab/Matrix3.hpp"
#include "moab/CartVect.hpp"
#include "moab/ElemUtil.hpp"
namespace moab { 

namespace element_utility {
//non-exported functionality
namespace { 

template< typename Moab, 
	  typename Entity_handle, typename Points>
void get_moab_points( Moab & moab, 
		      Entity_handle eh, 
		      Points & points){
	const Entity_handle* connectivity_begin;
	int num_vertices;
	moab.get_connectivity( eh, connectivity_begin, num_vertices);
	//TODO: This is hacky, it works correctly since
	//CartVect is only double d[ 3], with a default
	//constructor.. get_coords() should be
	//flexible enough to allow other types..
	points.resize( num_vertices);
	moab.get_coords( connectivity_begin, num_vertices, &(points[ 0][ 0]));
}

} // non-exported functionality

template< typename Element_map> 
class Element_parametrizer{
	public:
		//public types
		typedef Element_map Map;
	private: 
		typedef Element_parametrizer< Map> Self;
	public: //public functionality
	Element_parametrizer(): map(){}
 	Element_parametrizer( const Self & f): map( f.map) {}
	public:
		template< typename Moab, typename Entity_handle, typename Point>
		std::pair< bool, Point> operator()( Moab & moab,
						    const Entity_handle & eh, 
						    const Point & point, 
						    const double tol){
			typedef std::vector< moab::CartVect> Points;
			Points points;
			get_moab_points( moab, eh, points);
			return map( moab, eh, points, point, tol);
		}
	private: 
	Element_map map;
}; //class Element_parametrizer

class Parametrizer{
	private: 
		typedef Parametrizer Self;
		typedef moab::EntityHandle Entity_handle;
	public: //public functionality
	Parametrizer(): hex_map(), tet_map(){}
 	Parametrizer( const Self & f): hex_map( f.hex_map), 
				       tet_map( f.tet_map) {}
	public:
		template< typename Moab, typename Entity_handle, typename Point>
		std::pair< bool, Point> operator()( Moab & moab,
						    const Entity_handle & eh, 
						    const Point & point){
			//get entity
			typedef std::vector< moab::CartVect> Points;
		        Points points;
			get_moab_points( moab, eh, points);
			//get type
			switch( moab.type_from_handle( eh)){
 				case moab::MBHEX:
					return hex_map( moab, eh, 
							points, point);
				case moab::MBTET:
					return tet_map( moab, eh, 
							points, point);
				//TODO: not correct..
				//TODO: add quadratic hex, and a proper case for
				// spectral hex
				default:
					quadratic_hex_map( moab, eh, 
							   points, point);
					return spectral_hex_map( moab, eh, 
								 points, point);
				   std::cerr << "Element type not supported" 
					     << std::endl;
				  return make_pair( false, Point(3, 0.0));
			}
		}
		template< typename Moab, typename Entity_handle, typename Point>
		void interpolate( Moab & moab, const Entity_handle & eh, 
				  const Point & natural_coords){
			//get field values from moab tag, 
			//then 
			//evaluate_scalar_field(); 
		}
	private: 
	Linear_hex_map< moab::Matrix3> hex_map;
	Linear_tet_map< Entity_handle, moab::Matrix3> tet_map;
	Spectral_hex_map< moab::Matrix3> spectral_hex_map;
	Quadratic_hex_map< moab::Matrix3> quadratic_hex_map;
}; //class Parametrizer

}// namespace element_utility
} // namespace moab
#endif //MOAB_PARAMETRIZER_HPP
