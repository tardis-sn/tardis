/*-----------------------------------------------------------------------------

  Tuple list definition and utilities

  Conceptually, a tuple list is a list of n records or tuples,
  each with mi integers, ml longs, mul Ulongs, and mr reals
  (these types are defined in "types.h" as sint, slong, moab::EntityHandle, real;
  it may be that sint==slong)

  There are four arrays, one for each type (vi,vl,vul,vr),
  with records layed out contiguously within each array

  -----------------------------------------------------------------------------*/

#ifndef TUPLE_LIST_HPP
#define TUPLE_LIST_HPP

#include <limits.h>
#include <stdlib.h>

#include "moab/Types.hpp"
#include <string>

/* Integral types defined here to ensure variable type sizes are consistent */
/* integer type to use for everything */
#if   defined(USE_LONG)
#  define INTEGER long
#elif defined(USE_LONG_LONG)
#  define INTEGER long long
#elif defined(USE_SHORT)
#  define INTEGER short
#else
#  define INTEGER int
#endif

/* when defined, use the given type for global indices instead of INTEGER */
#if   defined(USE_GLOBAL_LONG_LONG)
#  define GLOBAL_INT long long
#elif defined(USE_GLOBAL_LONG)
#  define GLOBAL_INT long
#else
#  define GLOBAL_INT long
#endif

/* floating point type to use for everything */
#if   defined(USE_FLOAT)
typedef float realType;
#elif defined(USE_LONG_DOUBLE)
typedef long double realType;
#else
typedef double realType;
#endif

/* apparently uint and ulong can be defined already in standard headers */
#define uint uint_
//#define ulong ulong_
#define sint sint_
#define slong slong_

typedef   signed INTEGER sint;
typedef unsigned INTEGER uint;
#undef INTEGER

#ifdef GLOBAL_INT
typedef   signed GLOBAL_INT slong;
//typedef unsigned GLOBAL_INT ulong;
#else
typedef sint slong;
// typedef uint ulong;
#endif

typedef moab::EntityHandle Ulong;


namespace moab 
{
  void fail(const char *fmt, ...);

  class TupleList
  {
  public:
    /*---------------------------------------------------------------------------

      buffer: a simple class which can be used to store data
  
      The ptr points to the chunk of memory allocated for the buffer's use.  The 
      size denotes the size of the allocated memory; the user must take care to 
      note how much of the buffer they are using.

      ---------------------------------------------------------------------------*/
    class buffer
    {
    public:
      size_t buffSize;
      char *ptr;
  
      /**Constructor which sets an initial capacity of the buffer
       */
      buffer(size_t sz);

      /**Default constructor (Note:  buffer must be initialized before use!)
       */
      buffer();
  
      ~buffer () { this->reset(); };
  
      /**Initializes the buffer to have a capacity of size
       */
      void buffer_init_(size_t sz, const char *file);

      /**Ensures that the buffer has at least a capacity of min
       */
      void buffer_reserve_(size_t min, const char *file);
    
      /**Frees any allocated memory used by the buffer
       */
      void reset ();
  
      //Aliases for using the buffer methods
#define buffer_init(sz) buffer_init_(sz,__FILE__)
#define buffer_reserve(min) buffer_reserve_(min,__FILE__)

    };

  public:

    /**Constructor that takes all parameters and initializes the TupleList
     */
    TupleList(uint mi, uint ml,
	      uint mul, uint mr, uint max);

    /**Default constructor (Note:  TupleList must be initialized before use!)
     */
    TupleList();

    ~TupleList() { reset(); };

    /**Initializes the starting memory to be used by the TupleList
     * Note:  TupleLists must be initialized before they can be used
     *
     * param mi   number of signed ints in each tuple
     * param ml   number of long ints in each tuple
     * param mul  number of unsigned long ints in each tuple
     * param mr   number of reals in each tuple
     * param max  starting capacity of max tuples in the TupleList
     */
    void initialize(uint mi, uint ml,
		    uint mul, uint mr, uint max);

    /**Resizes the TupleList to a new max
     *
     * param max   the new max capacity of the TupleList
     */
    ErrorCode resize(uint max);

    /**Sorts the TupleList by 'key'
     * if key<mi:  key represents the key'th index in vi
     * if mi<key<ml:  key represents the (key-mi)'th index in vl
     * if ml<key<mul:  key represents the (key-mi-ml)'th index in vul
     *
     * param key   index to be sorted by
     * param *buf  buffer space used for sorting
     */
    /*------------------------------------------------------------------------------
  
      Hybrid Stable Sort
  
      low-overhead merge sort when n is small,
      otherwise asymptotically superior radix sort

      result = O(n) sort with good performance for all n
  
      A, n, stride : specifices the input
  
      sort:
      uint out[n] : the sorted values (output)
      uint work[n]: scratch area
  
      index_sort:
      uint idx[n]  : the sorted indices (output)
      sort_data work[2*n]: scratch area

      ----------------------------------------------------------------------------*/
    ErrorCode sort(uint key, TupleList::buffer *buf);

    /**Frees all allocated memory in use by the TupleList
     */
    void reset();

    /**Adds one to the total number of in-use tuples and resizes if necessary
     */
    void reserve();

    /**Finds index of the tuple containing 'value' at the key_numth index of 
     * said tuple; return -1 if key_num is out of bounds or if 'value' not found
     * Uses binary search if TupleList is sorted by the key_numth field, seqential
     * otherwise (very slow for large TupleLists; please sort before you search)
     *
     * param key_num   index of the tuple where to search for value
     * param value     value to search for at the given key_num
     * return the index of the tuple that contains value
     */
    int find(unsigned int key_num, sint value);
    int find(unsigned int key_num, slong value);
    int find(unsigned int key_num, Ulong value);
    int find(unsigned int key_num, realType value);

    /**get the mth number of return type in the index'th tuple
     * returns 0 if m or index is out of bounds
     *
     * param index     index of the tuple within the TupleList
     * param m         index of the value within the tuple
     * return the value at the given position
     */
    sint get_sint(unsigned int index, unsigned int m);
    slong get_int(unsigned int index, unsigned int m);
    Ulong get_ulong(unsigned int index, unsigned int m);
    realType get_double(unsigned int index, unsigned int m);

    /**get pointers to the data for the index'th tuple; ptr is
     * NULL if that type is not part of this tuple
     *
     * param index     index of the tuple needed
     * param *&sp, *&ip, *&lp, *&dp   pointers to each piece of the tuple
     */
    ErrorCode get(unsigned int index, const sint *&sp, 
		  const slong *&ip, const Ulong *&lp, const realType *&dp);

    /**push back a new tuple on the TupleList;
     *
     * param *sp   pointer to a list of signed ints
     * param *ip   pointer to a list of signed longs
     * param *lp   pointer to a list of unsigned longs
     * param *dp   pointer to a list of reals
     * return index of that tuple
     */
    unsigned int push_back(sint *sp, slong *ip,
			   Ulong *lp, realType *dp);

    /*Enable or disable direct write access to arrays
      This is important so that we know whether or not
      the user of this class can directly modify data
      which can affect operations such as
      whether or not we know the tuple list is sorted
      (for a binary search)*/
    void enableWriteAccess();
    void disableWriteAccess();

    /*Get information on the Tuple Sizes
     * param &mi_out  Count of uints in a tuple
     * param &ml_out  Count of uints in a tuple
     * param &mul_out Count of uints in a tuple
     * param &mr_out  Count of uints in a tuple
     */
    void getTupleSize(uint &mi_out, uint &ml_out, 
		      uint &mul_out, uint &mr_out) const;

    /*Set the count of Tuples in the Tuple List
     * Warning, automatically calls enableWriteAccess()
     * param n_in     New count of Tuples
     */
    void set_n(uint n_in);
  
    /* Get the count of Tuples in the Tuple List */
    uint get_n() const;
    
    /*Get the maximum number of Tuples currently allocated for*/
    uint get_max() const;

    bool get_writeEnabled() const;

    /*Increment n by 1
     * Warning, automatically calls enableWriteAccess()
     * returns current TupleList.n after the increment */
    uint inc_n();
    
    void print(const char *) const;

    //Variables to allow for direct write access
    sint *vi_wr; slong *vl_wr; Ulong *vul_wr; realType *vr_wr;

    //Variables to allow for direct read access
    const sint *vi_rd; slong *vl_rd; Ulong *vul_rd; realType *vr_rd;
    
  private:
    /* storage layed out as: vi[max][mi], vl[max][ml], vul[max][mul], 
     * vr[max][mr] where "tuple" i is given by 
     * (vi[i][0:mi-1],vl[i][0:ml-1],vul[i][0:mul-1],vr[i][0:mr-1]).
     * only the first n tuples are in use */
    uint mi,ml,mul,mr;
    uint n, max;
    sint *vi; slong *vl; Ulong *vul; realType *vr;

    // Used by sort:  see .cpp for more details
    //void sort_bits(uint *work, uint key);
    void permute(uint *perm, void *work);
    
    /* last_sorted = the last sorted position in the tuple (if the
     * TupleList has not been sorted, or has become unsorted--i.e.
     * by adding a tuple--last_sorted = -1) */
    int last_sorted;
    //Whether or not the object is currently allowing direct
    //write access to the arrays
    bool writeEnabled;

    typedef uint Index;

    template <typename Value>
    struct SortData {Value v; Index i;};

#define DIGIT_BITS   8
#define DIGIT_VALUES (1<<DIGIT_BITS)
#define DIGIT_MASK   ((Value)(DIGIT_VALUES-1))
#define CEILDIV(a,b) (((a)+(b)-1)/(b))
#define DIGITS       CEILDIV(CHAR_BIT*sizeof(Value),DIGIT_BITS)
#define VALUE_BITS   (DIGIT_BITS*DIGITS)
#define COUNT_SIZE   (DIGITS*DIGIT_VALUES)

    template<class Value> 
    static Value radix_count(const Value *A, const Value *end, Index stride,
			     Index count[DIGITS][DIGIT_VALUES]);
    
    static void radix_offsets(Index *c);

    template<class Value>
    static unsigned radix_zeros(Value bitorkey, Index count[DIGITS][DIGIT_VALUES],
				unsigned *shift, Index **offsets);

    template<class Value> 
    static void radix_index_pass_b(const Value *A, Index n, Index stride,
				   unsigned sh, Index *off, SortData<Value> *out);

    template<class Value>     
    static void radix_index_pass_m(const SortData<Value> *src, const SortData<Value> *end,
				   unsigned sh, Index *off, SortData<Value> *out);

    template<class Value> 
    static void radix_index_pass_e(const SortData<Value> *src, const SortData<Value> *end,
				   unsigned sh, Index *off, Index *out);

    template<class Value>
    static void radix_index_pass_be(const Value *A, Index n, Index stride,
				    unsigned sh, Index *off, Index *out);
    

    /*------------------------------------------------------------------------------

  
      Radix Sort
  
      stable; O(n) time

      ----------------------------------------------------------------------------*/
    template<class Value>
    static void radix_index_sort(const Value *A, Index n, Index stride,
				 Index *idx, SortData<Value> *work);

    /*------------------------------------------------------------------------------

  
      Merge Sort
  
      stable; O(n log n) time

      ----------------------------------------------------------------------------*/
    template<class Value>
    static void merge_index_sort(const Value *A, const Index An, Index stride,
				 Index *idx, SortData<Value> *work);

    template<class Value>
    static void index_sort(const Value *A, Index n, Index stride,
			   Index *idx, SortData<Value> *work);


#undef DIGIT_BITS
#undef DIGIT_VALUES
#undef DIGIT_MASK
#undef CEILDIV
#undef DIGITS
#undef VALUE_BITS
#undef COUNT_SIZE
  };

  inline uint TupleList::get_max() const{ return max; }
  inline uint TupleList::get_n() const{ return n; }
  inline bool TupleList::get_writeEnabled() const{ return writeEnabled; }

} //namespace
#endif
#include <stdlib.h>
