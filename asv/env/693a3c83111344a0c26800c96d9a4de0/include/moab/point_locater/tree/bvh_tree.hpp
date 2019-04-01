/** 
 * bvh_tree.hpp
 * Ryan H. Lewis
 * (C) 2012
 *
 * An element tree partitions a mesh composed of elements.
 * We subdivide the bounding box of a me, by putting boxes
 * on the left if there center is on the left of a split line
 * and vice versa.
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

//#define BVH_TREE_DEBUG
#ifndef BVH_TREE_HPP
#define BVH_TREE_HPP

namespace ct = moab::common_tree;

namespace moab {

//forward declarations
template< typename _Entity_handles, 
	  typename _Box, 
	  typename _Moab,
	  typename _Parametrizer> class Bvh_tree;

//non-exported functionality
namespace {
	namespace _bvh {
	template< typename Box, typename Entity_handle>
	struct _Node{
		typedef typename  std::vector< std::pair< Box, Entity_handle> > 
								Entities;
		std::size_t dim;
		std::size_t child;
		double Lmax, Rmin;
		Entities entities;
		_Node& operator=( const _Node& f){
			dim = f.dim;
			child=f.child;
			Lmax=f.Lmax;
			Rmin=f.Rmin;
			entities=f.entities;
			return *this;
		}
	}; // _Node


	template< typename Split>
	class Split_comparator : 
			public std::binary_function< Split, Split, bool> {
		inline double objective( const Split & a) const{
			return a.Lmax*a.nl - a.Rmin*a.nr;
		}
		public:
		bool operator()( const Split & a, const Split & b) const{
			return  objective( a) < objective( b);
		}
	}; //Split_comparator

	template< typename Iterator>
	class Iterator_comparator : 
			public std::binary_function< Iterator, Iterator, bool> {
		public:
		bool operator()( const Iterator & a, const Iterator & b) const{
			return a->second.second < b->second.second ||
				( !(b->second.second < a->second.second) 
					&& a->first < b->first);
		}
	}; //Split_comparator


	class _Split_data {
		public:
		typedef ct::Box< double> Box;
		_Split_data(): dim( 0), nl( 0), nr( 0), split( 0.0), 
				Lmax( 0.0), Rmin( 0.0),bounding_box(), 
				left_box(), right_box(){}
       		_Split_data( const _Split_data & f): 
			dim( f.dim), nl( f.nl), nr( f.nr), 
			split( f.split), Lmax( f.Lmax), Rmin( f.Rmin),
			bounding_box( f.bounding_box),
			left_box( f.left_box), right_box( f.right_box){}
		std::size_t dim;
		std::size_t nl;
		std::size_t nr;
		double split;
		double Lmax, Rmin;
		Box bounding_box;
		Box left_box;
		Box right_box;
		_Split_data& operator=( const _Split_data & f){
			dim  	     = f.dim;
			nl   	     = f.nl; 
			nr   	     = f.nr;
        		split	     = f.split;
			Lmax 	     = f.Lmax;
			Rmin 	     = f.Rmin;
        		bounding_box = f.bounding_box;
			left_box     = f.left_box;
			right_box    = f.right_box;
			return *this;
		}
	}; //_Split_data

	class _Bucket {
		public:
		_Bucket(): size( 0), bounding_box(){}
		_Bucket( const _Bucket & f): 
		size( f.size), bounding_box(f.bounding_box){}
		_Bucket( const std::size_t size_): 
		size( size_), bounding_box(){}
		std::size_t size;
		ct::Box< double> bounding_box;
		_Bucket& operator=( const _Bucket & f){
			bounding_box = f.bounding_box;
			size = f.size;
			return *this;
		}
	}; //_Split_data
	} // namespace _bvh
} //private namespace

template< typename _Entity_handles,
	  typename _Box,
	  typename _Moab,
	  typename _Parametrizer>
class Bvh_tree {
//public types
public:
	typedef  _Entity_handles Entity_handles;
	typedef  _Box Box;
	typedef  _Moab Moab;
	typedef  _Parametrizer Parametrizer;
	typedef typename Entity_handles::value_type Entity_handle;
	
//private types
private: 
	typedef Bvh_tree< _Entity_handles, 
			      _Box, 
			      _Moab,
			      _Parametrizer> Self;
	typedef typename std::pair< Box, Entity_handle> Leaf_element;
	typedef _bvh::_Node< Box, Entity_handle> Node;
	typedef typename std::vector< Node> Nodes;
//public methods
public:
//Constructor
Bvh_tree( Entity_handles & _entities, 
	  Moab & _moab, 
	  Box & _bounding_box, 
	  Parametrizer & _entity_contains): entity_handles_( _entities), 
				tree_(), moab( _moab), 
				bounding_box( _bounding_box),
				entity_contains( _entity_contains){
	typedef typename Entity_handles::iterator Entity_handle_iterator;
	typedef  ct::_Element_data< const _Box, double > Element_data;
	typedef typename std::tr1::unordered_map< Entity_handle, 
						  Element_data> Entity_map;
	typedef typename Entity_map::iterator Entity_map_iterator;
	typedef std::vector< Entity_map_iterator> Vector;
	//a fully balanced tree will have 2*_entities.size()
	//which is one doubling away..
	tree_.reserve( entity_handles_.size());
	Entity_map entity_map( entity_handles_.size());
	ct::construct_element_map( entity_handles_, 
					    entity_map, 
					    bounding_box, 
					    moab);
	#ifdef BVH_TREE_DEBUG
	for(Entity_map_iterator i = entity_map.begin(); 
				i != entity_map.end(); ++i){
		if( !box_contains_box( bounding_box, i->second.first, 0)){
			std::cerr << "BB:" << bounding_box << "EB:" <<
			i->second.first << std::endl;
			std::exit( -1);
		}
	}
	#endif
 	//_bounding_box = bounding_box;
	Vector entity_ordering;
	construct_ordering( entity_map, entity_ordering); 
	//We only build nonempty trees
	if( entity_ordering.size()){ 
	 //initially all bits are set
	 tree_.push_back( Node());
	 const int depth = build_tree( entity_ordering.begin(), 
				       entity_ordering.end(), 0, bounding_box);
	 #ifdef BVH_TREE_DEBUG
		 typedef typename Nodes::iterator Node_iterator;
		 typedef typename Node::Entities::iterator Entity_iterator;
		 std::set< Entity_handle> entity_handles;
		 for(Node_iterator i = tree_.begin(); i != tree_.end(); ++i){
			for(Entity_iterator j = i->entities.begin(); 
					    j != i->entities.end(); ++j){
				entity_handles.insert( j->second);
			}
				
		 }
		if( entity_handles.size() != entity_handles_.size()){
			std::cout << "Entity Handle Size Mismatch!" 
				  << std::endl;
		}
		typedef typename Entity_handles::iterator Entity_iterator_;
		for( Entity_iterator_ i  = entity_handles_.begin(); 
				     i != entity_handles_.end(); ++i){
			if ( entity_handles.find( *i) == entity_handles.end()){
				std::cout << "Tree is missing an entity! " 
					  << std::endl;
			}
		} 
					       
	 #endif
	 std::cout << "max tree depth: " << depth << std::endl; 
	}
}

//Copy constructor
Bvh_tree( Self & s): entity_handles_( s.entity_handles_), 
			 tree_( s.tree_), moab( s.moab), 
			 bounding_box( s.bounding_box),
			 entity_contains( s.entity_contains){}

//see FastMemoryEfficientCellLocationinUnstructuredGridsForVisualization.pdf 
//around page 9
#define NUM_SPLITS 4
#define NUM_BUCKETS (NUM_SPLITS + 1) //NUM_SPLITS+1
#define SMAX 5
//Paper arithmetic is over-optimized.. this is safer.
template < typename Box>
std::size_t bucket_index( const Box & box, const Box & interval, 
			  const std::size_t dim) const{
	const double min = interval.min[ dim];
	const double length = (interval.max[ dim]-min)/NUM_BUCKETS;
	const double center = ((box.max[ dim] + box.min[ dim])/2.0)-min;
	#ifdef BVH_TREE_DEBUG
	#ifdef BVH_SHOW_INDEX
	std::cout << "[ " << min << " , " 
		  << interval.max[ dim] << " ]" <<std::endl;
	std::cout << "[ " 
		<< box.min[ dim] << " , " << box.max[ dim] << " ]" <<std::endl;
	std::cout << "Length of bucket" << length << std::endl;
	std::cout << "Center: " 
			<< (box.max[ dim] + box.min[ dim])/2.0 << std::endl;
	std::cout << "Distance of center from min:  " << center << std::endl;
	std::cout << "ratio: " << center/length << std::endl;
	std::cout << "index: " << std::ceil(center/length)-1 << std::endl;
	#endif
	#endif
	return std::ceil(center/length)-1;
}

template< typename Iterator, typename Bounding_box, typename Buckets>
void establish_buckets( const Iterator begin, const Iterator end, 
			const Bounding_box & interval, 
			Buckets & buckets) const{
	//put each element into its bucket
	for(Iterator i = begin; i != end; ++i){
		const Bounding_box & box = (*i)->second.first;
		for (std::size_t dim = 0; dim < NUM_DIM; ++dim){
			const std::size_t index = bucket_index( box, 
								interval, dim);
			_bvh::_Bucket & bucket = buckets[ dim][ index];
			if(bucket.size > 0){
			     ct::update_bounding_box( bucket.bounding_box, box);
			}else{ 
				bucket.bounding_box = box; 
			}
			bucket.size++;
		}
	}
	#ifdef BVH_TREE_DEBUG
	Bounding_box elt_union = (*begin)->second.first;
	for(Iterator i = begin; i != end; ++i){
		const Bounding_box & box = (*i)->second.first;
		ct::update_bounding_box( elt_union, box);
		for (std::size_t dim = 0; dim < NUM_DIM; ++dim){
			const std::size_t index = bucket_index( box, 
								interval, dim);
			_bvh::_Bucket & bucket = buckets[ dim][ index];
			if(!box_contains_box( bucket.bounding_box, box)){
				std::cerr << "Buckets not covering elements!"
					  << std::endl;
			}
		}
	}
	if( !box_contains_box( elt_union, interval) ){
		std::cout << "element union: " << std::endl << elt_union; 
		std::cout << "intervals: " << std::endl << interval;
		std::cout << "union of elts does not contain original box!" 
			  << std::endl;
	}
	if ( !box_contains_box( interval, elt_union) ){
		std::cout << "original box does not contain union of elts" 
			  << std::endl;
		std::cout << interval << std::endl;
		std::cout << elt_union << std::endl;
	}
	typedef typename Buckets::value_type Bucket_list;
	typedef typename Bucket_list::const_iterator Bucket_iterator;
	for(std::size_t d = 0; d < NUM_DIM; ++d){
		std::vector< std::size_t> nonempty;
		const Bucket_list & buckets_ = buckets[ d];
		std::size_t j = 0;
		for(  Bucket_iterator i = buckets_.begin(); 
				      i != buckets_.end(); ++i, ++j){
		  if( i->size > 0){ nonempty.push_back( j); }
		}
		Bounding_box test_box = buckets_[ nonempty.front()].
							   bounding_box;
		for( std::size_t i = 0; i < nonempty.size(); ++i){
			ct::update_bounding_box( test_box, 
						 buckets_[ nonempty[ i]].
								bounding_box);
		}
		if( !box_contains_box( test_box, interval) ){
			std::cout << "union of buckets in dimension: " << d 
				  << "does not contain original box!" 
				  << std::endl;
		}
		if ( !box_contains_box( interval, test_box) ){
			std::cout << "original box does "
				  << "not contain union of buckets" 
				  << "in dimension: " << d 
				  << std::endl;
			std::cout << interval << std::endl;
			std::cout << test_box << std::endl;
		}
	}
	#endif
}

template< typename Box, typename Iterator>
std::size_t set_interval( Box & interval, const Iterator begin, 
			 	   const Iterator end) const{
	bool first=true;
	std::size_t count = 0;
	for( Iterator b = begin; b != end; ++b){
		const Box & box = b->bounding_box;
		count += b->size;
		if( b->size != 0){
			if( first){
			   interval = box;
			   first=false;
			}else{
			   ct::update_bounding_box( interval, box);
			}
		}
	}
	return count;
}

template< typename Splits, typename Buckets, typename Split_data>
void initialize_splits( Splits & splits, 
			const Buckets & buckets, 
			const Split_data & data) const{
	typedef typename Buckets::value_type Bucket_list;
	typedef typename Bucket_list::value_type Bucket;
	typedef typename Bucket_list::const_iterator Bucket_iterator;
	typedef typename Splits::value_type Split_list; 
	typedef typename Split_list::value_type Split; 
	typedef typename Split_list::iterator Split_iterator;
	for(std::size_t d = 0; d < NUM_DIM; ++d){
		const Split_iterator splits_begin = splits[ d].begin();
		const Split_iterator splits_end = splits[ d].end();
		const Bucket_iterator left_begin = buckets[ d].begin();
		const Bucket_iterator _end = buckets[ d].end();
		Bucket_iterator left_end = buckets[ d].begin()+1;
		for( Split_iterator s = splits_begin; 
				    s != splits_end; ++s, ++left_end){
			s->nl = set_interval( s->left_box, 
					      left_begin, left_end);
			if( s->nl == 0){ 
				s->left_box = data.bounding_box;
				s->left_box.max[ d] = s->left_box.min[ d];
			}
			s->nr = set_interval( s->right_box,
					      left_end,  _end);
			if( s->nr == 0){ 
				s->right_box = data.bounding_box;
				s->right_box.min[ d] = s->right_box.max[ d];
			}
			s->Lmax = s->left_box.max[ d];
			s->Rmin = s->right_box.min[ d];
			s->split = std::distance( splits_begin, s);
			s->dim = d;
		}
		#ifdef BVH_TREE_DEBUG
		for( Split_iterator s = splits_begin; 
				    s != splits_end; ++s, ++left_end){
			typename Split::Box test_box = s->left_box;
			ct::update_bounding_box( test_box, s->right_box);
			if( !box_contains_box( data.bounding_box, test_box)){
				std::cout << "nr: " << s->nr << std::endl;
				std::cout << "Test box: " << std::endl << 
					  test_box;
				std::cout << "Left box: " << std::endl << 
		  			  s->left_box;
				std::cout << "Right box: " << std::endl << 
	  				  s->right_box;
				std::cout << "Interval: " << std::endl << 
	  				  data.bounding_box;
				std::cout << "Split boxes larger than bb" 
					  << std::endl;
			}
			if( !box_contains_box( test_box, data.bounding_box)){
				std::cout << "bb larger than union "
					  << "of split boxes" 
					  << std::endl;
			}         	
			}
		#endif 
	}
}

template< typename Iterator, typename Split_data>
void order_elements( const Iterator & begin, const Iterator & end, 
		     const Split_data & data) const{
	typedef typename Iterator::value_type Map_iterator;
	for(Iterator i = begin; i != end; ++i){
		const int index = bucket_index( (*i)->second.first,
						data.bounding_box, data.dim);
		(*i)->second.second = (index<=data.split)?0:1;
	}
	std::sort( begin, end, _bvh::Iterator_comparator< Map_iterator>());
}

template< typename Iterator, typename Split_data>
void median_order( const Iterator & begin, const Iterator & end, 
		      Split_data & data) const{
	typedef typename Iterator::value_type Map_iterator;
	for(Iterator i = begin; i != end; ++i){
		const double center = 
		       compute_box_center((*i)->second.first, data.dim);
		(*i)->second.second = center; 
	}
	std::sort( begin, end, _bvh::Iterator_comparator< Map_iterator>());
	const std::size_t total = std::distance( begin, end);
	Iterator middle = begin+(total/2);
	double middle_center = (*middle)->second.second;
	       middle_center += (*(++middle))->second.second;
	       middle_center /=2.0;
	data.split = middle_center;
	data.nl = std::distance( begin, middle)+1;
	data.nr = total-data.nl;
	middle++;
	data.left_box  = (*begin)->second.first;
	data.right_box = (*middle)->second.first;
	for(Iterator i = begin; i != middle; ++i){
		(*i)->second.second = 0;
		update_bounding_box( data.left_box, (*i)->second.first);
	}
	for(Iterator i = middle; i != end; ++i){
		(*i)->second.second = 1;
		update_bounding_box( data.right_box, 
				     (*i)->second.first);
	}
	data.Rmin = data.right_box.min[ data.dim];
	data.Lmax = data.left_box.max[ data.dim];
	#ifdef BVH_TREE_DEBUG
	typename Split_data::Box test_box = data.left_box;
	ct::update_bounding_box( test_box, data.right_box);
	if( !box_contains_box( data.bounding_box, test_box) ){
		std::cerr << "MEDIAN: BB Does not contain splits" << std::endl;
	}
	if( !box_contains_box( test_box, data.bounding_box) ){
		std::cerr << "MEDIAN: splits do not contain BB" << std::endl;
	}
	#endif
}

template< typename Splits, typename Split_data>
void choose_best_split( const Splits & splits, Split_data & data) const{
	typedef typename Splits::const_iterator List_iterator;
	typedef typename List_iterator::value_type::const_iterator 
							Split_iterator;
	typedef typename Split_iterator::value_type Split;
	std::vector< Split> best_splits;
	typedef typename _bvh::Split_comparator< Split> Comparator;
	Comparator compare;
	for( List_iterator i = splits.begin(); i != splits.end(); ++i){
		Split_iterator j = std::min_element( i->begin(), i->end(), 
								  compare);
		best_splits.push_back( *j);
	}
	data = *std::min_element( best_splits.begin(), 
				  best_splits.end(), compare);
}


template< typename Iterator, typename Split_data>
void find_split(const Iterator & begin, 
		const Iterator & end, Split_data & data) const{
	typedef typename Iterator::value_type Map_iterator;
	typedef typename Map_iterator::value_type::second_type Box_data;
	typedef typename Box_data::first_type _Bounding_box; // Note, not global typedef moab::common_tree::Box< double> Bounding_box;
	typedef typename std::vector< Split_data> Split_list;
	typedef typename std::vector< Split_list> Splits;
	typedef typename Splits::iterator Split_iterator;
	typedef typename std::vector< _bvh::_Bucket> Bucket_list;
	typedef typename std::vector< Bucket_list > Buckets;
	Buckets buckets( NUM_DIM, Bucket_list( NUM_BUCKETS) );
	Splits splits( NUM_DIM, Split_list( NUM_SPLITS, data));
	
	const _Bounding_box interval = data.bounding_box;
	establish_buckets( begin, end, interval, buckets);
	initialize_splits( splits, buckets, data);
	choose_best_split( splits, data);
	const bool use_median = (0 == data.nl) || (data.nr == 0);
	if (!use_median){ order_elements( begin, end, data); } 
	else{ median_order( begin, end, data); }
	#ifdef BVH_TREE_DEBUG
	bool seen_one=false,issue=false;
	bool first_left=true,first_right=true;
	std::size_t count_left=0, count_right=0;
	typename Split_data::Box left_box, right_box;
	for( Iterator i = begin; i != end; ++i){
		double order = (*i)->second.second;
		if( order != 0 && order != 1){
			std::cerr << "Invalid order element !";
			std::cerr << order << std::endl;
			std::exit( -1);
		}
		if(order == 1){
			seen_one=1;
			count_right++;
			if( first_right){
				right_box = (*i)->second.first;
				first_right=false;
			}else{
				ct::update_bounding_box( right_box, 
							 (*i)->second.first);
			}
			if(!box_contains_box( data.right_box, 
					     (*i)->second.first)){
				if(!issue){
				std::cerr << "Bounding right box issue!" 
					  << std::endl;
				}
				issue=true;
			}
		}
		if(order==0){
			count_left++;
			if( first_left){
				left_box = (*i)->second.first;
				first_left=false;
			}else{
				ct::update_bounding_box( left_box, 
							 (*i)->second.first);
			}
			if(!box_contains_box( data.left_box, 
					     (*i)->second.first)){
				if(!issue){
				std::cerr << "Bounding left box issue!" 
					 << std::endl;
				}
				issue=true;
			}
			if(seen_one){
				std::cerr << "Invalid ordering!" << std::endl;
				std::cout << (*(i-1))->second.second 
					  << order << std::endl;
				exit( -1);
			}
		}
	}
	if( !box_contains_box( left_box, data.left_box)){
		std::cout << "left elts do not contain left box" << std::endl;
	}
	if( !box_contains_box( data.left_box, left_box)){
		std::cout << "left box does not contain left elts" << std::endl;
	}
	if( !box_contains_box( right_box, data.right_box)){
		std::cout << "right elts do not contain right box" << std::endl;
	}
	if( !box_contains_box( data.right_box, right_box)){
		std::cout << "right box do not contain right elts" << std::endl;
	}
	if( count_left != data.nl || count_right != data.nr) {
		std::cerr << "counts are off!" << std::endl;
		std::cerr << "total: " 
			  << std::distance( begin, end) << std::endl;
		std::cerr << "Dim: " << data.dim << std::endl;
		std::cerr << data.Lmax << " , " << data.Rmin << std::endl;
		std::cerr << "Right box: " << std::endl << data.right_box 
			  << "Left box: " << std::endl << data.left_box ;
		std::cerr << "supposed to be: " << 
					data.nl << " " << data.nr << std::endl;
		std::cerr << "accountant says: " << 
			count_left << " " << count_right << std::endl;
		std::exit( -1);
	}
	#endif
}

//private functionality
private:
template< typename Iterator>
int build_tree( const Iterator begin, const Iterator end, 
		const int index, const Box & box, 
		const int depth=0){
	#ifdef BVH_TREE_DEBUG
	for(Iterator i = begin; 
		     i != end; ++i){
		if( !box_contains_box( box, (*i)->second.first, 0)){
			std::cerr << "depth: " << depth << std::endl;
			std::cerr << "BB:" << box << "EB:" <<
			(*i)->second.first << std::endl;
			std::exit( -1);
		}
	}
	#endif

	const std::size_t total_num_elements = std::distance( begin, end);
	Node & node = tree_[ index];
	//logic for splitting conditions
	if( total_num_elements > SMAX){
		_bvh::_Split_data data;
		data.bounding_box = box;
		find_split( begin, end, data);
		//assign data to node
		node.Lmax = data.Lmax; node.Rmin = data.Rmin;
		node.dim = data.dim; node.child = tree_.size();
		//insert left, right children;
		tree_.push_back( Node()); tree_.push_back( Node());
		const int left_depth=
		build_tree( begin, begin+data.nl, node.child, 
			    data.left_box, depth+1);
		const int right_depth=
		build_tree( begin+data.nl, end, node.child+1, 
			    data.right_box, depth+1);
		return std::max( left_depth, right_depth);
	}
	node.dim = 3;
	ct::assign_entities( node.entities, begin, end);
	return depth;
}

template< typename Vector, typename Node_index, typename Result>
Result& _find_point( const Vector & point, 
			   const Node_index & index,
			   const double tol,
			   Result& result) const{
	typedef typename Node::Entities::const_iterator Entity_iterator;
	const Node & node = tree_[ index];
	if( node.dim == 3){
		for( Entity_iterator i = node.entities.begin(); 
				     i != node.entities.end(); ++i){
			if( ct::box_contains_point( i->first, point, tol)){
				const std::pair< bool, Vector> r = 
				entity_contains( moab, i->second, point);
				if (r.first){
				    return result = std::make_pair( i->second, 
								    r.second);
				}
			}
		}
		result = Result(0, point);
		return result;
	}
        //the extra tol here considers the case where
        //0 < Rmin - Lmax < 2tol
	if( (node.Lmax+tol) < (node.Rmin-tol)){
		if( point[ node.dim] <= (node.Lmax + tol)){            	
        		return _find_point( point, node.child, tol, result);
        	}else if( point[ node.dim] >= (node.Rmin - tol)){            	
			return _find_point( point, node.child+1, tol, result);
		}
		result = Result(0, point);
		return result; //point lies in empty space.
	}
	//Boxes overlap
	//left of Rmin, you must be on the left
	//we can't be sure about the boundaries since the boxes overlap
	//this was a typo in the paper which caused pain.
	if( point[ node.dim] < (node.Rmin - tol)){
		return _find_point( point, node.child, tol, result);
	//if you are on the right Lmax, you must be on the right
	}else if( point[ node.dim] > (node.Lmax+tol)){
		return _find_point( point, node.child+1, tol, result);
	}
	/* pg5 of paper
	 * However, instead of always traversing either subtree
	 * first (e.g. left always before right), we first traverse 
	 * the subtree whose bounding plane has the larger distance to the 
	 * sought point. This results in less overall traversal, and the correct
	 * cell is identified more quickly.
	 */
	//Sofar all testing confirms that this 'heuristic' is 
	//significantly slower.
	//I conjecture this is because it gets improperly
	//branch predicted..
	//bool dir = (point[ node.dim] - node.Rmin) <= 
	//				(node.Lmax - point[ node.dim]);
	bool dir=0;
	_find_point( point, node.child+dir, tol, result);
	if( result.first == 0 ){ 
		return _find_point( point, node.child+(!dir), tol, result);
	}
	return result;
}

//public functionality
public:
template< typename Vector, typename Result>
Result& find( const Vector & point, 
		    const double tol, Result & result) const{
	typedef typename Vector::const_iterator Point_iterator;
	return  _find_point( point, 0, tol, result);
}

//public functionality
public:
template< typename Vector>
Entity_handle bruteforce_find( const Vector & point, const double tol) const{
	typedef typename Vector::const_iterator Point_iterator;
	typedef typename Nodes::value_type Node;
	typedef typename Nodes::const_iterator Node_iterator;
	typedef typename Node::Entities::const_iterator Entity_iterator;
	for( Node_iterator i = tree_.begin(); i != tree_.end(); ++i){
		if( i->dim == 3){
			for( Entity_iterator j = i->entities.begin();
					     j != i->entities.end();
						++j){
				if( ct::box_contains_point( j->first, 
							    point, tol)){
				      const std::pair< bool, Vector> result = 
				      entity_contains( moab, j->second, point);
				      if (result.first){
				      	return j->second;
				      }
				}
			}
		}
	}
	return 0;
}


//public accessor methods
public:

//private data members  
private:
	const Entity_handles & entity_handles_;
	Nodes tree_;
	Moab & moab;
	Box bounding_box;
	Parametrizer & entity_contains;

}; //class Bvh_tree

} // namespace moab

#endif //BVH_TREE_HPP
