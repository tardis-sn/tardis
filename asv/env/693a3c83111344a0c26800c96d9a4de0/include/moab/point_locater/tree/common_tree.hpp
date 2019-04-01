/** 
 * common_tree.hpp
 * Ryan H. Lewis
 * (C) 2012
 * Functionality common to all trees.
 */
#include <vector>
#include <set>
#include <iostream>
#include <map>
#include <algorithm>
#include <bitset>
#include <numeric>
#include <cmath>
#include <tr1/unordered_map>
#include <limits>

#ifndef COMMON_TREE_HPP
#define COMMON_TREE_HPP
#define NUM_DIM 3
#define TREE_DEBUG
namespace moab {
namespace common_tree {
 
template< typename T, typename Stream>
void print_vector( const T & v, Stream & out){
	typedef typename T::const_iterator Iterator;
	out << "[ ";
	for(Iterator i = v.begin(); i != v.end(); ++i){
		out << *i;
		if( i+1 != v.end()){
			out << ", ";
		}
	}
	out << " ]" << std::endl;
}

#ifdef TREE_DEBUG
template< typename T>
void print_vector( const T & begin, const T & end){
	std::cout << "[ ";
	for(T i = begin; i != end; ++i){
		std::cout << (*i)->second.second.to_ulong();
		if( i+1 != end){
			std::cout << ", ";
		}
	}
	std::cout << " ]" << std::endl;
}
#endif

template< typename _Box, typename _Point>
bool box_contains_point(  const _Box & box, const _Point & p, 
			  const double tol){
	for( std::size_t i = 0; i < box.min.size(); ++i){
	     if( p[ i] < (box.min[ i]-tol) || 
	         p[ i] > (box.max[ i])+tol){
	     	return false;
	     }
	 }
	return true;
}

template< typename _Box>
bool box_contains_box(  const _Box & a, const _Box & b, const double tol){
	for( std::size_t i = 0; i < a.min.size(); ++i){
	     if( b.min[ i] < (a.min[ i]-tol)){
		return false;
	      }  
	      if( b.max[ i] > (a.max[ i]+tol)){
	     	return false;
	     }
	 }
	return true;
}

namespace {
	template< typename T> 
	struct Compute_center: public std::binary_function< T, T, T> {
		T operator()( const T a, const T b) const{
			return (a+b)/2.0;
		}
	}; //Compute_center
} //non-exported center computation.

template< typename Vector>
inline void compute_box_center(Vector & max, Vector & min, Vector & center){
	typedef typename Vector::value_type Unit;
	center = min;
	std::transform( max.begin(), max.end(), center.begin(),
			center.begin(), Compute_center< Unit>() );
}


template< typename Box>
inline typename Box::value_type compute_box_center( const Box & box, 
						    const int dim){
	return (box.max[ dim] + box.min[ dim])/2.0;
}


template< typename T = float>
class Box{
	public:
	typedef T value_type;
	typedef std::vector< T> Vector;
	Box(): max(3,0.0), min(3,0.0) {}
	Box( const Box & from): max( from.max), 
				min( from.min){}
	template< typename Iterator>
	Box( const Iterator begin, const Iterator end):
	max( begin, end), min(begin, end){}
	Box& operator=( const Box & from){
		max = from.max;
		min = from.min;
		return *this;
	}
	Vector max;
	Vector min;
}; //Box

template< typename T>
std::ostream& operator<<( std::ostream& out, const Box<T> & box){
	out << "Max: ";
	print_vector( box.max, out);
	out << "Min: ";
	print_vector( box.min, out);	
	return out;
}

//essentially a pair, but with an added constructor.
template< typename T1, typename T2>
struct _Element_data  {
	typedef T1 first_type;
	typedef T2 second_type;
	T1 first;
	T2 second;
	_Element_data(): first(T1()), second(T2()){}
	_Element_data( const T1 & x): first(x), second(T2()) {}
	_Element_data( const T1 & x, T2 & y): first( x), second( y) {}
	template< typename U, typename V>
	_Element_data( const _Element_data<U,V> &p ): first(p.first),
						      second( p.second) {}
}; //Element_data

template< typename Entities, typename Iterator>
void assign_entities( Entities & entities, 
		      const Iterator & begin, const Iterator & end){
	entities.reserve( std::distance( begin, end)); 
	for( Iterator i = begin; i != end; ++i){
		entities.push_back( std::make_pair((*i)->second.first, 
						    (*i)->first));
	}
}

template<typename Coordinate, typename Coordinate_iterator>
void update_bounding_max( Coordinate & max, Coordinate_iterator j){
	typedef typename Coordinate::iterator Iterator;
	for( Iterator i = max.begin(); i != max.end(); ++i, ++j){
		*i = std::max( *i, *j);
	}
}

template<typename Coordinate, typename Coordinate_iterator>
void update_bounding_min( Coordinate & min, Coordinate_iterator j){
	typedef typename Coordinate::iterator Iterator;
	for( Iterator i = min.begin(); i != min.end(); ++i, ++j){
		*i = std::min( *i, *j);
	}
}

template< typename Box>
void update_bounding_box( Box & a, const Box & b){
	update_bounding_max( a.max, b.max.begin());
	update_bounding_min( a.min, b.min.begin());
	#ifdef COMMON_TREE_DEBUG
	if( !box_contains_box( a, b)){
		std::cout << a << b << std::endl;	
	}
	#endif
}


template< typename Entity_map, typename Ordering>
void construct_ordering( Entity_map & entity_map, Ordering & entity_ordering){
	entity_ordering.reserve( entity_map.size());
	typedef typename Entity_map::iterator Map_iterator;
	for(Map_iterator i = entity_map.begin(); 
			 i != entity_map.end(); 
			 ++i){
		entity_ordering.push_back( i); 
	}
}

//Input: A bunch of entity handles
//Output: A map from handle -> Data 
//Requirements: Data contains at least a bounding box.
//And a non-default constructor which takes only a Box&
template< typename Entity_handles, 
	  typename Element_map, 
	  typename Bounding_box,
	  typename Moab>
void construct_element_map( const Entity_handles & elements, 
			    Element_map & map, 
			    Bounding_box & bounding_box,
			    Moab & moab){
	typedef typename Element_map::mapped_type Box_data;
	typedef typename Entity_handles::value_type Entity_handle;
	typedef typename Entity_handles::iterator Entity_handles_iterator;
	typedef typename Box_data::first_type::value_type Unit;
	typedef typename std::vector< Unit> Coordinates;
	typedef typename Coordinates::iterator Coordinate_iterator;
		
	for( Entity_handles_iterator i = elements.begin(); 
				     i != elements.end(); ++i){
		//TODO: not generic enough. Why dim != 3
		const int DIM = 3;	
		int num_vertices=0;
		//Commence un-necessary deep copying.
		const Entity_handle* vertex_handle;
		moab.get_connectivity( *i, vertex_handle, num_vertices);
		Coordinates coordinate(DIM*num_vertices, 0.0);
		moab.get_coords( vertex_handle, num_vertices, &coordinate[ 0]);
		Bounding_box box( coordinate.begin(), coordinate.begin()+3);
		if( i == elements.begin() ){ bounding_box = box;}
		for( Coordinate_iterator j = coordinate.begin()+DIM; 
				         j != coordinate.end(); j+=DIM){
			update_bounding_max( box.max, j);
			update_bounding_min( box.min, j);
		}
		update_bounding_box( bounding_box, box);
		map.insert( std::make_pair( *i, Box_data( box)));
	}
}

} //namspace common_tree

} // namespace moab

#endif //COMMON_TREE_HPP
