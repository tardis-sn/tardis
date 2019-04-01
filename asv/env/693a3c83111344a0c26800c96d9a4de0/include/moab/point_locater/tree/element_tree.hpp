/** 
 * element_tree.hpp
 * Ryan H. Lewis
 * (C) 2012
 *
 * An element tree partitions a mesh composed of elements.
 * We subdivide the bounding box of a mesh, and each element is
 * either entirely on the left, entirely on the right, or crossing
 * the diving line. We build a tree on the mesh with this property.
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

#include "common_tree.hpp"

#ifndef ELEMENT_TREE_HPP
#define ELEMENT_TREE_HPP
namespace moab {
//forward declarations

template< typename _Entity_handles, 
	  typename _Box, 
	  typename _Moab,
	  typename _Parametrizer> class Element_tree;

//non-exported functionality
namespace  {
namespace _element_tree { 
template< typename Iterator>
struct Iterator_comparator{
typedef typename Iterator::value_type Value;
bool operator()(const Value & a, const Value & b){
	return a->second.second.to_ulong() < b->second.second.to_ulong();
}
}; //Iterator_comparator

template< typename Data>
struct Split_comparator {
  //we minimizes ||left| - |right|| + |middle|^2 
  double split_objective( const Data & a) const {
	if (a.second.sizes[ 2]==0 || a.second.sizes[ 0] == 0){
		return std::numeric_limits< std::size_t>::max();
	}
	const double total = a.second.sizes[ 0]+a.second.sizes[ 2];
	const int max = 2*(a.second.sizes[ 2]>a.second.sizes[ 0]);
	
	return (a.second.sizes[ max] -a.second.sizes[ 2*(1-(max==2))])/total;
  }
  bool operator()( const Data & a, const Data & b) const {
  	return split_objective( a) < split_objective( b);
  }
}; //Split_comparator

template< typename Partition_data, typename Box>
void correct_bounding_box( const Partition_data & data, Box & box, 
			   const int child){
	const int dim = data.dim;
	switch( child){
		case 0:
			box.max[ dim] = data.left_rightline;
			break;
		case 1:
			box.max[ dim] = data.right_line;
			box.min[ dim] = data.left_line;
			break;
		case 2:
			box.min[ dim] = data.right_leftline;
			break;
	}
	#ifdef ELEMENT_TREE_DEBUG
	print_vector( data.bounding_box.max);
	print_vector( data.bounding_box.min);
	print_vector( box.max);
	print_vector( box.min);
	#endif
}

template< typename Box>
struct _Partition_data{
	typedef _Partition_data< Box> Self;
	//default constructor
	_Partition_data():sizes(3,0),dim(0){}
	_Partition_data( const Self & f){ *this=f; }
	_Partition_data( const Box & _box, int _dim): sizes(3,0),
	bounding_box( _box), split((_box.max[ _dim] + _box.min[ _dim])/2.0), 
	left_line( split), right_line( split), dim( _dim){}
	_Partition_data& operator=( const Self & f){
		sizes = f.sizes;
		bounding_box = f.bounding_box;
		split = f.split;
		left_line = f.left_line;
		right_line = f.right_line;
		right_leftline=f.right_leftline;
		left_rightline=f.left_rightline;
		dim = f.dim;
		return *this;
	}
	std::vector< std::size_t> sizes;
	Box bounding_box;
	double split;
	double left_line;
	double right_line;
	double right_leftline;
	double left_rightline;
	int dim;
	std::size_t left()  const { return sizes[ 0]; }
	std::size_t middle()const { return sizes[ 1]; }
	std::size_t right() const { return sizes[ 2]; }
}; // Partition_data

template< typename _Entity_handles, typename _Entities>
class _Node{
	//public types:
	public:
	typedef _Entity_handles Entity_handles;
	typedef _Entities Entities;

	//private types:
	private:
	typedef _Node< _Entity_handles, _Entities> Self;

	//Constructors
	public:
	//Default constructor
	_Node(): children(3,-1), 
	left_( children[ 0]), middle_( children[ 1]),
	right_( children[ 2]),
	dim( -1), split( 0),
		 left_line( 0), right_line( 0),
		entities( 0) {}

	//Copy constructor
	_Node( const Self & from):
	children( from.children), 
	left_( children[ 0]), middle_( children[ 1]),
	right_( children[ 2]),
	dim( from.dim), split( from.split),
	left_line( from.left_line), right_line( from.right_line), 
	entities( from.entities) {}

	public:
	template< typename Iterator>
	void assign_entities(const Iterator & begin, const Iterator & end){
		entities.reserve( std::distance( begin, end)); 
		for( Iterator i = begin; i != end; ++i){
			entities.push_back( std::make_pair((*i)->second.first, 
							    (*i)->first));
		}
	}

	// Functionality
	public: 
	bool leaf() const { return children[ 0] == -1 && 
			           children[ 1] == -1 && 
			    	   children[ 2] == -1; }
	Self& operator=( const Self & from){
		children=from.children;
		dim=from.dim;
		left_ = from.left_;
		middle_ = from.middle_;
		right_ = from.right_;
		split=from.split;
		left_line=from.left_line;
		right_line=from.right_line;
		entities=from.entities;
		return *this;
	}
	template< typename Box>
	Self& operator=( const _Partition_data< Box> & from){
		dim=from.dim;
		split=from.split;
		left_line=from.left_line;
		right_line=from.right_line;
		return *this;
	}
	//private data members:
	private:
	//indices of children
	std::vector< int> children;
	int & left_, middle_, right_;
	int dim; //split dimension
	double split; //split position
	double left_line;
	double right_line;
	Entities entities;

	//Element_tree can touch my privates.
	template< Entity_handles, typename B> friend class moab::Element_tree;
}; //class Node

} //namespace _element_tree
} // anon namespace 

template< typename _Entity_handles,
	  typename _Box,
	  typename _Moab, 
	  typename _Parametrizer>
class Element_tree {

//public types
public:
	typedef  _Entity_handles Entity_handles;
	typedef  _Box Box;
	typedef  _Moab Moab;
	typedef _Parametrizer Parametrizer;
	typedef typename Entity_handles::value_type Entity_handle;
	
//private types
private: 
	typedef Element_tree< _Entity_handles, 
			      _Box, 
			      _Moab, 
			      _Parametrizer> Self; 
	typedef std::pair< Box, Entity_handle>  Leaf_element;
	typedef _element_tree::_Node< Entity_handles, std::vector< Leaf_element> > Node;
	//int is because we only need to store 	
	#define MAX_ITERATIONS 2
	typedef  common_tree::_Element_data< Box, std::bitset<NUM_DIM*MAX_ITERATIONS*2> > 
								Element_data;
	typedef std::vector< Node> Nodes;
	//TODO: we really want an unordered map here, make sure this is kosher..
	typedef std::tr1::unordered_map< Entity_handle, Element_data> 
								Element_map;
	typedef typename std::vector< typename Element_map::iterator> 
								Element_list;
	typedef _element_tree::_Partition_data< Box> Partition_data;
//public methods
public:
//Constructor
Element_tree( Entity_handles & _entities, Moab & _moab, Box & _bounding_box, 
	      Parametrizer & _entity_contains): 
	entity_handles_( _entities), tree_(), moab( _moab), 
	bounding_box( _bounding_box), entity_contains( _entity_contains){
	tree_.reserve( _entities.size());
	Element_map element_map( _entities.size());
	Partition_data _data;
	common_tree::construct_element_map( entity_handles_, element_map, 
						_data.bounding_box, moab);
	bounding_box = _data.bounding_box;
	_bounding_box = bounding_box;
	Element_list element_ordering( element_map.size());
	std::size_t index = 0;
	for(typename Element_map::iterator i = element_map.begin(); 
					  i != element_map.end(); ++i, ++index){
		element_ordering[ index] = i;
	}
	//We only build nonempty trees
	if( element_ordering.size()){ 
		//initially all bits are set
		std::bitset< 3> directions( 7);
		tree_.push_back( Node());
		int depth = 0;
		build_tree( element_ordering.begin(), 
			    element_ordering.end(),
			    0, directions, _data, depth);
		std::cout << "depth: " << depth << std::endl; 
	}
}

//Copy constructor
Element_tree( Self & s): entity_handles_( s.entity_handles_), 
			 tree_( s.tree_), moab( s.moab), 
			 bounding_box( s.bounding_box){}

//private functionality
private:

template< typename Iterator, typename Split_data>
void compute_split( Iterator & begin, Iterator & end, 
			 Split_data & split_data, bool iteration=false){
	typedef typename Iterator::value_type::value_type Map_value_type;
	typedef typename Map_value_type::second_type::second_type Bitset;
	//we will update the left/right line
	double & left_line = split_data.left_line;
	double & right_line = split_data.right_line;
	double & split = split_data.split;
	const int & dim = split_data.dim;
	#ifdef ELEMENT_TREE_DEBUG
 	std::cout << std::endl; 
	std::cout << "-------------------" << std::endl; 
	std::cout << "dim: " << dim << " split: " << split << std::endl;
	std::cout << "bounding_box min: "; 
	print_vector( split_data.bounding_box.min); 
	std::cout << "bounding_box max: "; 
	print_vector( split_data.bounding_box.max);
	#endif
	//for each elt determine if left/middle/right
	for(Iterator i = begin; i != end; ++i){
		const Box & box =  (*i)->second.first;
		Bitset & bits =  (*i)->second.second;
		//will be 0 if on left, will be 1 if in the middle
		//and 2 if on the right;
		const bool on_left = (box.max[ dim] < split);
		const bool on_right = (box.min[ dim] > split);
		const bool in_middle = !on_left && !on_right;
		//set the corresponding bits in the bit vector
		// looks like: [x_1 = 00 | x_2 = 00 | .. | z_1 = 00 | z_2 = 00]
		// two bits, left = 00, middle = 01, right = 10
		const int index = 4*dim + 2*iteration;
		if( on_left){ split_data.sizes[ 0]++; }
		else if(in_middle){
			split_data.sizes[ 1]++;
			bits.set( index, 1);
			left_line  = std::min( left_line,  box.min[ dim]);
			right_line = std::max( right_line, box.max[ dim]);
		}else if( on_right){
			bits.set( index+1, 1);
			split_data.sizes[ 2]++;
		}
	}
	#ifdef ELEMENT_TREE_DEBUG
	std::size_t _count = std::accumulate( split_data.sizes.begin(), 
					      split_data.sizes.end(), 0);
	std::size_t total = std::distance( begin, end);
	if( total != _count ){
		std::cout << total << "vs. " <<  _count << std::endl;

	}
	std::cout <<  " left_line: " << left_line;
	std::cout <<  " right_line: " << right_line << std::endl;
	std::cout << "co/mputed partition size: ";
	print_vector( split_data.sizes);
	std::cout << "-------------------" << std::endl; 
	#endif
}

template< typename Split_data>
bool update_split_line( Split_data & data) const{
	const int max = 2*(data.sizes[ 2]>data.sizes[ 0]);
	const int min = 2*(1-(max==2));
	bool one_side_empty = data.sizes[ max]==0 || data.sizes[ min]==0;
	double balance_ratio = data.sizes[ max] - data.sizes[ min];
	//if ( !one_side_empty && balance_ratio < .05*total){ return false; }
	if( !one_side_empty){
		//if we have some imbalance on left/right 
		//try to fix the situation 
		balance_ratio /= data.sizes[ max];
		data.split += (max-1)*balance_ratio*(data.split/2.0);
	}else{
		//if the (left) side is empty move the split line just past the 
		//extent of the (left) line of the middle box.
		//if middle encompasses everything then wiggle 
		//the split line a bit and hope for the best..
		const double left_distance = 
					std::abs(data.left_line-data.split);  	
        	const double right_distance = 
					std::abs(data.right_line-data.split);
		if( (data.sizes[ 0] == 0) && (data.sizes[ 2] != 0)){
			data.split += right_distance;
		}else if (data.sizes[ 2]==0 && data.sizes[ 0] != 0){
			data.split -= left_distance;
		}else{
			data.split *=1.05;
		}
	}
	data.left_line = data.right_line = data.split;
	data.sizes.assign( data.sizes.size(), 0);
	return true;
}

template< typename Iterator, 
	  typename Split_data, 
	  typename Directions>
void determine_split( Iterator & begin, 
		      Iterator & end, 
		      Split_data & data, 
		      const Directions & directions){ 
	typedef typename Iterator::value_type Pair;
	typedef typename Pair::value_type Map_value_type;
	typedef typename Map_value_type::second_type::second_type Bitset;
	typedef typename Map_value_type::second_type::first_type Box;
	typedef typename std::map< std::size_t, Split_data> Splits;
	typedef typename Splits::value_type Split;	
	typedef _element_tree::Split_comparator< Split> Comparator;
	Splits splits;
	for (std::size_t dir = 0; dir < directions.size(); ++dir){
		if( directions.test( dir)){
			Split_data split_data( data.bounding_box, dir);
			compute_split( begin, end, split_data);
			splits.insert( std::make_pair(2*dir, split_data));
			if( update_split_line( split_data)){
				compute_split( begin, end, split_data, true);
				splits.insert( std::make_pair( 2*dir+1,
							       split_data) );
			}
		}
	}
	Split best = *std::min_element( splits.begin(), splits.end(), 
					Comparator());
	#ifdef ELEMENT_TREE_DEBUG
	std::cout << "best: " << Comparator().split_objective( best) << " ";
	print_vector( best.second.sizes);
	#endif 
	const int dir = best.first/2;
	const int iter = best.first%2;
	double & left_rightline = best.second.left_rightline= 
	         			best.second.bounding_box.min[ dir];
	double & right_leftline = best.second.right_leftline = 
					best.second.bounding_box.max[ dir];
	Bitset mask( 0);
	mask.flip( 4*dir+2*iter).flip( 4*dir+2*iter+1);
	for(Iterator i = begin; i != end; ++i){
		Bitset & bits = (*i)->second.second;
		const Box & box =  (*i)->second.first;
		//replace 12 bits with just two.
		bits &= mask;
		bits >>= 4*dir+2*iter;
		//if box is labeled left/right but properly contained 
		//in the middle, move the element into the middle.
		//we can shrink the size of left/right
		switch( bits.to_ulong()){
			case 0:
				if ( box.max[ dir] > best.second.left_line){ 
					left_rightline = std::max( 
						left_rightline, box.max[ dir]);
				}
				break;
			case 2:
				if ( box.min[ dir] < best.second.right_line){ 
					right_leftline = std::min( 
						right_leftline, box.max[ dir]);
				}
				break;
		}
	}
	data = best.second;
}

//define here for now.
#define ELEMENTS_PER_LEAF 30
#define MAX_DEPTH 30
#define EPSILON 1e-1
template< typename Iterator, typename Node_index, 
	  typename Directions, typename Partition_data>
void build_tree( Iterator begin, Iterator end, 
		 const Node_index node, 
		 const Directions & directions, 
		 Partition_data & _data,
		 int & depth, 
		 const bool is_middle = false){
	std::size_t number_elements = std::distance(begin, end);
	if ( depth < MAX_DEPTH && 
		number_elements > ELEMENTS_PER_LEAF && 
		(!is_middle || directions.any())){
		determine_split( begin, end,  _data, directions);
		//count_sort( begin, end, _data);
		std::sort( begin, end,
			 _element_tree::Iterator_comparator< Iterator>());
		//update the tree
		tree_[ node] = _data;
		Iterator middle_begin( begin+_data.left());
		Iterator middle_end( middle_begin+_data.middle());
		std::vector< int> depths(3, depth);
		//left subtree
		if( _data.left()>0){
			Partition_data data( _data);
			tree_.push_back( Node());
			tree_[ node].children[ 0] = tree_.size()-1;
			correct_bounding_box( _data, data.bounding_box, 0);
			Directions new_directions( directions);
			const bool axis_is_very_small = 
			 (data.bounding_box.max[ _data.dim] - 
			    data.bounding_box.min[ _data.dim] < EPSILON);
			new_directions.set( _data.dim, axis_is_very_small); 
			build_tree( begin, middle_begin, 
				    tree_[ node].children[ 0], 
				    new_directions, data, ++depths[ 0],
				    is_middle);
		}
		//middle subtree
		if( _data.middle()>0){
			Partition_data data( _data);
			tree_.push_back( Node());
			tree_[ node].children[ 1] = tree_.size()-1;
			correct_bounding_box( _data, data.bounding_box, 1);
			//force the middle subtree to split
			//in a different direction from this one
			Directions new_directions( directions);
			new_directions.flip( tree_[node].dim);
			bool axis_is_very_small = 
			 (data.bounding_box.max[ _data.dim] - 
			    data.bounding_box.min[ _data.dim] < EPSILON);
			new_directions.set( _data.dim, axis_is_very_small); 
			build_tree( middle_begin, middle_end,
				    tree_[ node].children[ 1], 
				    new_directions, data, 
				    ++depths[ 1], true);
		}
		//right subtree
		if( _data.right()>0){
			Partition_data data( _data);
			tree_.push_back( Node());
			tree_[ node].children[ 2] = tree_.size()-1;
			correct_bounding_box( _data, data.bounding_box, 2);
			Directions new_directions( directions);
			const bool axis_is_very_small = 
			 (data.bounding_box.max[ _data.dim] - 
			    data.bounding_box.min[ _data.dim] < EPSILON);
			new_directions.set( _data.dim, axis_is_very_small); 

			build_tree( middle_end, end, tree_[ node].children[ 2],
				    directions, data, ++depths[ 2], is_middle);
		}
		depth = *std::max_element(depths.begin(), depths.end());
	}
	if( tree_[ node].leaf()){
		common_tree::assign_entities(tree_[ node].entities, begin, end);
	}
}

template< typename Vector, typename Node_index, typename Result>
Result& _find_point( const Vector & point, 
	             const Node_index & index,
		     Result & result) const{
	typedef typename Node::Entities::const_iterator Entity_iterator;
	typedef typename std::pair< bool, Vector> Return_type;
	const Node & node = tree_[ index];
	if( node.leaf()){
		//check each node
		for( Entity_iterator i = node.entities.begin(); 
				     i != node.entities.end(); ++i){
			if( common_tree::box_contains_point( i->first, point)){
				Return_type r = entity_contains( moab, 
							         i->second, 
								 point);
				if( r.first){ result = 
					std::make_pair( i->second, r.second);
				}
				return result;
			}
		}
		return Result(0, point);
	}
	if( point[ node.dim] < node.left_line){
		return _find_point( point, node.left_, result);
	}else if( point[ node.dim] > node.right_line){
		return _find_point( point, node.right_, result);
	} else {
		Entity_handle middle =  _find_point( point, 
							node.middle_, result);
		if( middle != 0){ return result; }
		if( point[ node.dim] < node.split){ 
			return _find_point( point, node.left_, result); 
		}
		return	_find_point( point, node.right_, result);
	}
}

//public functionality
public:
template< typename Vector, typename Result>
Result& find( const Vector & point, Result & result) const{
	typedef typename Vector::const_iterator Point_iterator;
	typedef typename Box::Pair Pair; 
	typedef typename Pair::first_type Box_iterator;
	return  _find_point( point, 0, result);
}
	

//public accessor methods
public:

//private data members  
private:
	const Entity_handles & entity_handles_;
	Nodes tree_;
	Moab & moab;
	Box bounding_box;
	Parametrizer entity_contains;

}; //class Element_tree

} // namespace moab

#endif //ELEMENT_TREE_HPP
