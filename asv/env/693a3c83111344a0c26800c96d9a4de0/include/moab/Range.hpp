/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/*!
  \brief Stores contiguous or partially contiguous values in an optimized
  fashion.  Partially contiguous accessing patterns is also optimized.
 
  \author Clinton Stimpson
 
  \date 15 April 2002

 *************  Range FAQ and tips ********************

 The purpose of this FAQ is to familiarize a user with
 the appropriate use of the Range template.

   ******* A few points about Range: *******
 1.  Range is not the be all of generic containers.
 2.  Range has its strengths and weaknesses as any other
     STL container has.
 3.  Strengths:
     a. For contiguous values, storage is extremely minimal.
     b. Searching through contiguous values, at best, is
        a constant time operation.
     b. Fairly compatible with most STL algorithms.
     c. Insertions of data from high value to low value
        is a linear operation (constant for each insertion).
     d. Removal of a value using an iterator is constant time.

 4.  Weaknesses:
     a. For non-contiguous values, storage is not minimal and is
        on the order of 4x the storage space as using a vector.
     b. Searching through non-contiguous values is linear
        time operation.
     c. Insertions of random data is VERY slow.

   Given the above characteristics of Ranges, you can now
   decide between Range and another STL container for your
   particular needs.


   ******* Tips *******
 1.  Be aware of what iterators are available.  Currently, there are
     three.  Range<T>::iterator, Range<T>::const_iterator,
     and Range<T>::RangeListIterator.
     iterator is derived from const_iterator.  const_iterator
     is derived from RangeListIterator.  RangeListIterator is a
     std::list<std::pair<T,T> >::const_iterator.
     If a particular algorithm could be more efficient by using
     RangeListIterator, do so.

     ie.
     
     Range<char> range1;
     ... put some stuff in range1
     Range<char> range2;
     
     // the SLOW way.
     std::copy(range1.begin(), range1.end(), range_inserter<...>(range2));

     // the FAST way.
     for(Range<char>::RangeListIterator iter = range1.begin(),
         iter != range1.end(); ++iter)
     {
       range2.insert(iter->first, iter->second);
     }

 2.  Prefer insert(val1, val2) over insert(val) where possible.
 
 3.  insert(val) and insert(val1, val2) have to perform searching
     to find out where to insert an item.  Searches are started
     from the beginning of the values.  Inserting larger values 
     before smaller values will increase efficiency.

     ie.
     std::set<int> my_set;
     Range<int> my_range;
     .. perform some operations which set does efficiently.

     // now copy the results from the set into the range.
     // copy from the end of the set to the beginning of the set
     std::copy(my_set.rbegin(), my_set.rend(), 
         range_inserter< Range<int> > ( my_range );

 4.  Use empty() instead of size() if you only need to find out
     if there is anything in the list.

 5.  Know what swap() does.  Sometimes it is useful.  It'll replace
     the contents of one list with another.

     void compute_and_get_some_set( 
         Range<char> range1, 
         Range<char> range2,
         Range<char>& results
         );
     {
       Range<char> tmp_results;
       .. perform some computation on range1 and range2
          and put results in tmp_results;
       .. filter tmp_results out for some special type.
       .. etc....
       // return results
       results.swap(tmp_results);
     }


   ******* FAQ *******
 1. Why won't this code compile?
    ------------------------
    class SomeClass
    {
    public:
      Range<int> range;
    };
    ....
    void some_function( const SomeClass& some_class )
    {
      Range<int>::iterator = some_class.range.begin();
    }
    ---------------------
    Solution:  you need to use
    Range<int>::const_iterator instead.

 2. Why doesn't this work right when I try to change the
    contents of an Range?

    // make a range that has the letters A,B,C in it.
    Range<char> my_chars('A', 'C');
    // make an iterator that points to 'A'.
    Range<char>::iterator iter = my_chars.begin();
    // change A to D
    *iter = 'D';
    // print contents of my_chars to stdout
    std::copy(my_chars.begin(), my_chars.end(),
              std::ostream_iterator(std::cout, " "));

    result is 
      A B C
    instead of
      B C D

    When one calls *iter, which returns 'A', the actual storage of the value 'A' 
    which is returned is in the iterator and not in the Range.  This allows
    for multiple iterators on a single Range and the Range does not have
    to keep track of which iterator is referencing which value.



*/


#ifndef MOAB_RANGE_HPP
#define MOAB_RANGE_HPP

#include <iterator>
#include <iosfwd>
#include <algorithm>
#include "moab/Types.hpp"

namespace moab {

struct range_iter_tag : public std::bidirectional_iterator_tag {};

struct range_base_iter
{
  typedef range_iter_tag iterator_category;
  typedef EntityID difference_type;
  typedef EntityHandle value_type;
  typedef EntityHandle* pointer;
  typedef EntityHandle& reference;
};


//! the class Range
class Range
{
public:

    // forward declare the iterators
  class const_iterator;
  class const_reverse_iterator;
  typedef const_iterator iterator;
  typedef const_reverse_iterator reverse_iterator;
 
  friend Range intersect( const Range&, const Range& );
  friend Range subtract( const Range&, const Range& );

    //! just like subtract, but as an operator
  Range &operator-=(const Range &rhs);
    
  //! for short hand notation, lets typedef the 
  //! container class that holds the ranges
  typedef EntityHandle value_type;

  //! default constructor
  Range();

    //! copy constructor
  Range(const Range& copy);

  //! another constructor that takes an initial range
  Range( EntityHandle val1, EntityHandle val2 );

    //! operator=
  Range& operator=(const Range& copy);
  
  //! destructor
  inline ~Range();

  //! return the beginning const iterator of this range
  inline const_iterator begin() const;
  
  //! return the beginning const reverse iterator of this range
  inline const_reverse_iterator rbegin() const;
 
  //! return the ending const iterator for this range
  inline const_iterator end() const;
  
  //! return the ending const reverse iterator for this range
  inline const_reverse_iterator rend() const;

  //! return the number of values this Ranges represents
  size_t size() const;

  //! return the number of range pairs in the list
  size_t psize() const;
  
  //! return whether empty or not 
  //! always use "if(!Ranges::empty())" instead of "if(Ranges::size())"
  inline bool empty() const;

  iterator insert( iterator hint, EntityHandle val );

  //! insert an item into the list and return the iterator for the inserted item
  iterator insert(EntityHandle val)
    { return insert( begin(), val ); }
  
  //! insert a range of items into this list and return the iterator for the first
  //! inserted item
  iterator insert( iterator hint, EntityHandle first, EntityHandle last );
  
  //! insert a range of items into this list and return the iterator for the first
  //! inserted item
  iterator insert(EntityHandle val1, EntityHandle val2)
    { return insert( begin(), val1, val2 ); }
    
  template <typename T> 
  iterator insert_list( T begin_iter, T end_iter );
    
  template <class T> 
  iterator insert( typename T::const_iterator begin_iter, typename T::const_iterator end_iter )
    { return insert_list( begin_iter, end_iter ); }
    
  template <typename T> 
  iterator insert( const T* begin_iter, const T* end_iter )
    { return insert_list( begin_iter, end_iter ); }
  
    //! remove an item from this list and return an iterator to the next item
  iterator erase(iterator iter);

  //! remove a range of items from the list
  iterator erase( iterator iter1, iterator iter2);

  //! erases a value from this container
  inline iterator erase(EntityHandle val);
  
  //! get first entity in range
  inline const EntityHandle& front() const;
  //! get last entity in range
  inline const EntityHandle& back() const;
  //! remove first entity from range
  EntityHandle pop_front();
  //! remove last entity from range
  EntityHandle pop_back();
  
  //! find an item int the list and return an iterator at that value
  const_iterator find(EntityHandle val) const;

  //! return an iterator to the first value >= val
  static const_iterator lower_bound(const_iterator first,
                                    const_iterator last,
                                    EntityHandle val);
  static const_iterator upper_bound(const_iterator first,
                                    const_iterator last,
                                    EntityHandle val);
  
  const_iterator lower_bound( EntityHandle val ) const 
    { return lower_bound( begin(), end(), val ); }
  const_iterator upper_bound( EntityHandle val ) const 
    { return upper_bound( begin(), end(), val ); }
  const_iterator lower_bound( EntityType type ) const;
  const_iterator upper_bound( EntityType type ) const;
  std::pair<const_iterator, const_iterator> equal_range( EntityType type ) const;
  const_iterator lower_bound( EntityType type, const_iterator first ) const;
  const_iterator upper_bound( EntityType type, const_iterator first ) const;
  
  //! True if all entities in range are of passed type 
  //! (also true if range is empty)
  bool all_of_type( EntityType type ) const;
  //! True if all entities in range are of passed dimension 
  //! (also true if range is empty)
  bool all_of_dimension( int dimension ) const;
  
  unsigned num_of_type( EntityType type ) const;
  unsigned num_of_dimension( int dim ) const;
  
  //! clears the contents of the list 
  void clear();
  
  //! for debugging
  void print(const char *indent_prefix = NULL) const;
  void print(std::ostream& s, const char *indent_prefix = NULL) const;
  
  unsigned long get_memory_use() const;

  double compactness() const;
  
  void insert( Range::const_iterator begin,
               Range::const_iterator end );

  //! merges this Range with another range
  void merge( const Range& range )
    { insert( range.begin(), range.end() ); }
  
  //! merge a subset of some other range
  void merge( Range::const_iterator beginr,
              Range::const_iterator endr )
    { insert( beginr, endr ); }

  //! swap the contents of this range with another one
  void swap( Range &range );

    //! check for internal consistency
  void sanity_check() const;
  
    //! Check if this range is a non-strict superset of some other range
  bool contains( const Range& other ) const;

    //! return a subset of this range, by type
  Range subset_by_type(EntityType t) const;
  
    //! return a subset of this range, by dimension
  Range subset_by_dimension(int dim) const;
  
  struct PairNode : public std::pair<EntityHandle,EntityHandle>
  {

    PairNode() : std::pair<EntityHandle,EntityHandle>(0, 0), mNext(NULL), mPrev(NULL) {}
    PairNode(PairNode* next, PairNode* prev, 
             EntityHandle _first, EntityHandle _second)
      : std::pair<EntityHandle,EntityHandle>(_first,_second), mNext(next), mPrev(prev) {}

    PairNode* mNext;
    PairNode* mPrev;
  };

  
  EntityHandle operator[](EntityID index) const;

  int index(EntityHandle handle) const;
  
protected:

  //! the head of the list that contains pairs that represent the ranges 
  //! this list is sorted and unique at all times
  PairNode mHead;
  
  //! if dead_node is not mHead, remove it from the list and free it's memory.
  void delete_pair_node( PairNode* dead_node );

public:

    //! used to iterate over sub-ranges of a range
  class pair_iterator : public range_base_iter
  {
    friend class Range;
  public:
    pair_iterator() : mNode(NULL) {}
    pair_iterator(PairNode *nodep) : mNode(nodep) {}
    pair_iterator(const pair_iterator& copy)
      : mNode(copy.mNode) {}
    pair_iterator(const const_iterator& copy)
      : mNode(copy.mNode) {}

    std::pair<EntityHandle,EntityHandle>* operator->() { return mNode; }
    
    pair_iterator& operator++()
    {
      mNode = mNode->mNext;
      return *this;
    }
    pair_iterator operator++(int)
    {
      pair_iterator tmp(*this);
      this->operator ++();
      return tmp;
    }

    pair_iterator& operator--()
    {
      mNode = mNode->mPrev;
      return *this;
    }
    pair_iterator operator--(int)
    {
      pair_iterator tmp(*this);
      this->operator--();
      return tmp;
    }
    bool operator==(const pair_iterator& other) const
    {
      return mNode == other.mNode;
    }

    bool operator!=(const pair_iterator& other) const
    {
      return mNode != other.mNode;
    }
  
    PairNode* node() { return mNode; }

  private:
    
    PairNode* mNode;
  };

  class const_pair_iterator;

  //! a const iterator which iterates over an Range
  class const_iterator : public range_base_iter
  {
    friend class Range;
    friend class pair_iterator;
    friend class const_pair_iterator;
    friend EntityID operator-( const const_iterator&, const const_iterator& );
  public:
    //! default constructor - intialize base default constructor
    const_iterator() : mNode(NULL), mValue(0) {}

    //! constructor used by Range
    const_iterator( const PairNode* iter, const EntityHandle val) 
      : mNode(const_cast<PairNode*>(iter)), mValue(val)  {} 

    //! dereference that value this iterator points to
    //! returns a const reference
    const EntityHandle& operator*() const { return  mValue; }

    //! prefix incrementer
    const_iterator& operator++()
    {
      // see if we need to increment the base iterator
      if(mValue == mNode->second)
      {
        mNode = mNode->mNext;
        mValue = mNode->first;
      }
      // if not, just increment the value in the range
      else
        ++mValue;
      return *this;
    }

    //! postfix incrementer
    const_iterator operator++(int)
    {
      // make a temporary copy
      const_iterator tmp(*this);
      // increment self
      this->operator ++();
      // return the copy
      return tmp;
    }

    //! prefix decrementer
    const_iterator& operator--()
    {
      // see if we need to decrement the base iterator
      if(mValue == mNode->first)
      {
        mNode = mNode->mPrev;;
        mValue = mNode->second;
      }
      // if not, just decrement the value
      else
        --mValue;
      return *this;
    }

    //! postfix decrementer
    const_iterator operator--(int)
    {
      // make a copy of this
      const_iterator tmp(*this);
      // decrement self
      this->operator --();
      // return the copy
      return tmp;
    }
    
    //! Advance iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator++ step times.
    const_iterator& operator+=( EntityID step );
    
    //! Regress iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator-- step times.
    const_iterator& operator-=( EntityID step );

    //! equals operator
    bool operator==( const const_iterator& other ) const
    {
      // see if the base iterator is the same and the
      // value of this iterator is the same
      return (mNode == other.mNode) && (mValue == other.mValue);
    }

    //! not equals operator
    bool operator!=( const const_iterator& other ) const
    {
      // call == operator and not it.
      return (mNode != other.mNode) || (mValue != other.mValue);
    }
    
    /**\brief get an iterator at the end of the block
     *
     * Get an iterator at the end of the block of consecutive
     * handles that this iterator is currently contained in.
     * That is, if the range contains blocks of consecutive 
     * handles of the form { [1,5], [7,100], ... } and this
     * iterator is at any handle in the range [7,100], return
     * an iterator at the '100' handle.
     *
     * Never returns begin() or end() unless this iterator is
     * at begin() or end().  May return the same location as
     * this iterator.
     */
    inline const_iterator end_of_block() const;
    
    /**\brief get an iterator at the start of the block
     *
     * Get an iterator at the start of the block of consecutive
     * handles that this iterator is currently contained in.
     * That is, if the range contains blocks of consecutive 
     * handles of the form { [1,5], [7,100], ... } and this
     * iterator is at any handle in the range [7,100], return
     * an iterator at the '7' handle.
     *
     * Never returns end() unless this iterator is
     * at end().  May return the same location as
     * this iterator.
     */
    inline const_iterator start_of_block() const;

  protected:

    //! the node we are pointing at
    PairNode* mNode;
    //! the value in the range
    EntityHandle mValue;
  };

  //! a const reverse iterator which iterates over an Range
  class const_reverse_iterator : public range_base_iter
  {
    friend class Range;
    friend class pair_iterator;
  public:
    //! default constructor - intialize base default constructor
    const_reverse_iterator() {}
    
    const_reverse_iterator( const_iterator fwd_iter ) : myIter(fwd_iter) {}

    //! constructor used by Range
    const_reverse_iterator( const PairNode* iter, const EntityHandle val) 
      : myIter(iter, val)  {} 

    //! dereference that value this iterator points to
    //! returns a const reference
    const EntityHandle& operator*() const { return  *myIter; }

    //! prefix incrementer
    const_reverse_iterator& operator++()
    {
      --myIter;
      return *this;
    }

    //! postfix incrementer
    const_reverse_iterator operator++(int)
    {
      return const_reverse_iterator( myIter-- );
    }

    //! prefix decrementer
    const_reverse_iterator& operator--()
    {
      ++myIter;
      return *this;
    }

    //! postfix decrementer
    const_reverse_iterator operator--(int)
    {
      return const_reverse_iterator( myIter++ );
    }
    
    //! Advance iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator++ step times.
    const_reverse_iterator& operator+=( EntityID step )
    {
      myIter -= step;
      return *this;
    }

    //! Regress iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator-- step times.
    const_reverse_iterator& operator-=( EntityID step )
    {
      myIter += step;
      return *this;
    }

    //! equals operator
    bool operator==( const const_reverse_iterator& other ) const
    {
      return myIter == other.myIter;
    }

    //! not equals operator
    bool operator!=( const const_reverse_iterator& other ) const
    {
      return myIter != other.myIter;
    }

  protected:

    //! the node we are pointing at
    const_iterator myIter;
  };

public:

  class const_pair_iterator {
    public:
      const_pair_iterator() : myNode(NULL) {}
      const_pair_iterator( const PairNode* node ) : myNode(node) {}
      const_pair_iterator( const const_iterator& i ) : myNode(i.mNode) {}
      
      const PairNode& operator*() const
        { return *myNode; }
      
      const PairNode* operator->() const
        { return myNode; }
      
      const_pair_iterator& operator--()
        { myNode = myNode->mPrev; return *this; }
      
      const_pair_iterator& operator++()
        { myNode = myNode->mNext; return *this; }
      
      const_pair_iterator operator--(int)
        { const_pair_iterator rval(*this); this->operator--(); return rval; }
      
      const_pair_iterator operator++(int)
        { const_pair_iterator rval(*this); this->operator++(); return rval; }
        
      bool operator==( const const_pair_iterator& other ) const
        { return other.myNode == myNode; }
    
      bool operator!=( const const_pair_iterator& other ) const
        { return other.myNode != myNode; }

    private:
      const PairNode* myNode;
  };
  
  pair_iterator pair_begin() { return pair_iterator(mHead.mNext); }
  pair_iterator pair_end() { return pair_iterator(&mHead); }

  const_pair_iterator const_pair_begin() const { return const_pair_iterator( mHead.mNext ); }
  const_pair_iterator const_pair_end() const { return const_pair_iterator( &mHead ); }
  const_pair_iterator pair_begin() const { return const_pair_iterator( mHead.mNext ); }
  const_pair_iterator pair_end() const { return const_pair_iterator( &mHead ); }

};

 
    //! intersect two ranges, placing the results in the return range
Range intersect( const Range&, const Range& );

    //! subtract range2 from this, placing the results in the return range
Range subtract( const Range& from, const Range& );

    //! unite two ranges, placing the results in the return range
inline Range unite( const Range& r1, const Range& r2 )
  { Range r(r1); r.insert(r2.begin(), r2.end()); return r; }


inline Range::const_iterator 
operator+( const Range::const_iterator& it, EntityID step )
  { Range::const_iterator tmp(it); return tmp += step; }
  
inline Range::const_iterator 
operator+( EntityID step, const Range::const_iterator& it )
  { Range::const_iterator tmp(it); return tmp += step; }
  
inline Range::const_iterator 
operator-( const Range::const_iterator& it, EntityID step )
  { Range::const_iterator tmp(it); return tmp -= step; }
  
inline Range::const_iterator 
operator-( EntityID step, const Range::const_iterator& it )
  { Range::const_iterator tmp(it); return tmp -= step; }
  
EntityID 
operator-(  const Range::const_iterator& it1, const Range::const_iterator& it2 );

//! Use as you would an STL back_inserter
/**
 *  e.g. std::copy(list.begin(), list.end(), range_inserter(my_range);
 * Also, see comments/instructions at the top of this class declaration
 */
class range_inserter 
{
  
protected:
  Range* container;
 
public:
  //constructor
  explicit range_inserter(Range& x) : container(&x) {}
  range_inserter&
  operator=(const Range::value_type& value) 
  {
    container->insert(value);
    return *this;
  }

  range_inserter& operator*() { return *this; }
  range_inserter& operator++() { return *this; }
  range_inserter& operator++(int) { return *this; }

  typedef EntityHandle            value_type;
  typedef EntityID                difference_type;
  typedef std::output_iterator_tag  iterator_category;
  typedef EntityHandle*           pointer;
  typedef EntityHandle&           reference;
};


inline Range::Range()
{
    // set the head node to point to itself
  mHead.mNext = mHead.mPrev = &mHead;
  mHead.first = mHead.second = 0;
}
  
  //! destructor
inline Range::~Range()
{
  clear();
}

  //! return the beginning const iterator of this range
inline Range::const_iterator Range::begin() const
{
  return const_iterator(mHead.mNext, mHead.mNext->first);
}
  
  //! return the beginning const reverse iterator of this range
inline Range::const_reverse_iterator Range::rbegin() const
{
  return const_reverse_iterator(mHead.mPrev, mHead.mPrev->second);
}
 
  //! return the ending const iterator for this range
inline Range::const_iterator Range::end() const
{
  return const_iterator(&mHead, mHead.first);
}
  
  //! return the ending const reverse iterator for this range
inline Range::const_reverse_iterator Range::rend() const
{
  return const_reverse_iterator(&mHead, mHead.second);
}

  //! return whether empty or not 
  //! always use "if(!Ranges::empty())" instead of "if(Ranges::size())"
inline bool Range::empty() const
{
  return (mHead.mNext == &mHead);
}

  //! erases a value from this container
inline Range::iterator Range::erase(EntityHandle val) 
{ 
  return erase(find(val)); 
}
  
inline Range::const_iterator Range::const_iterator::end_of_block() const
  { return Range::const_iterator( mNode, mNode->second ); }

inline Range::const_iterator Range::const_iterator::start_of_block() const
  { return Range::const_iterator( mNode, mNode->first ); }

  //! get first entity in range
inline const EntityHandle& Range::front() const
  { return mHead.mNext->first; }
  //! get last entity in range
inline const EntityHandle& Range::back() const
  { return mHead.mPrev->second; }

inline std::ostream& operator<<( std::ostream& s, const Range& r )
  { r.print(s); return s; }
  
bool operator==( const Range& r1, const Range& r2 );
inline bool operator!=( const Range& r1, const Range& r2 )
  { return !(r1 == r2); }

inline EntityHandle Range::operator[](EntityID indexp) const
{
  Range::const_iterator i = begin();
  i += indexp;
  return *i;
}

inline int Range::index(EntityHandle handle) const 
{
  if (handle < *begin() || handle > *rbegin()) return -1;
  
  unsigned int i = 0;
  Range::const_pair_iterator pit = const_pair_begin(); 
  while (handle > (*pit).second && pit != const_pair_end()) {
    i += (*pit).second - (*pit).first + 1;
    ++pit;
  }
  if (handle < (*pit).first || pit == const_pair_end()) return -1;
  
  return i + handle - (*pit).first;
}

inline double Range::compactness() const 
{
  unsigned int num_ents = size();
  return (num_ents ? ((double)get_memory_use() / (double) (num_ents*sizeof(EntityHandle))) : -1);
}
    
template <typename Iterator> 
Range::iterator Range::insert_list( Iterator begin_iter, Iterator end_iter )
{
  size_t n = std::distance(begin_iter, end_iter);
  EntityHandle* sorted = new EntityHandle[n];
  std::copy( begin_iter, end_iter, sorted );
  std::sort( sorted, sorted + n );
  iterator hint = begin();
  size_t i = 0;
  while (i < n) {
    size_t j = i + 1;
    while (j < n && sorted[j] == 1+sorted[j-1])
      ++j;
    hint = insert( hint, sorted[i], sorted[i] + ((j - i) - 1) );
    i = j;
  }
  delete [] sorted;
  return hint;
}

inline size_t Range::psize() const 
{
  size_t i = 0;
  Range::const_pair_iterator pit;
  for (pit = const_pair_begin(), i = 0; 
       pit != const_pair_end(); ++pit, i++);

  return i;
}

} // namespace moab 

#endif // MOAB_RANGE_HPP



