#ifndef MB_SETITERATOR_HPP
#define MB_SETITERATOR_HPP

#include "moab/Interface.hpp"

namespace moab {

class Core;
    
/** \class Meshset iterator
 * \brief A class to iterator over MOAB Meshsets
 */
class SetIterator
{
public:

  friend class Core;
  
  //! destructor
  virtual ~SetIterator();

    //! get the ent set for this iterator
  inline EntityHandle ent_set() const {return entSet;};
  
    //! get the chunk size of this iterator
  inline unsigned int chunk_size() const {return chunkSize;};

    //! get the entity type for this iterator
  inline EntityType ent_type() const {return entType;};

    //! get the dimension for this iterator
  inline int ent_dimension() const {return entDimension;};

    /** \brief get the next chunkSize entities
     * Return the next chunkSize entities.  
     * \param arr Array of entities returned.
     * \param atend Returns true if iterator is at the end of iterable values, otherwise false
     */
  virtual ErrorCode get_next_arr(std::vector<EntityHandle> &arr,
                                 bool &atend) = 0;
  
    //! reset the iterator to the beginning of the set
  virtual ErrorCode reset() = 0;
  
protected:

    /** \brief Constructor
     * \param core MOAB Core instance
     * \param ent_set EntitySet to which this iterator corresponds
     * \param chunk_size Chunk size of this iterator
     * \param ent_type Entity type for this iterator
     * \param ent_dim Entity dimension for this iterator
     */
  inline SetIterator(Core *core, EntityHandle eset, unsigned int chunk_sz, 
                     EntityType ent_tp, int ent_dim, bool check_valid = false) 
          : myCore(core), entSet(eset), chunkSize(chunk_sz), 
            entType(ent_tp), entDimension(ent_dim), checkValid(check_valid) {};
  
    //! Core instance 
  Core *myCore;
  
    //! handle for entity set corresponding to this iterator
  EntityHandle entSet;
  
    //! chunk size of this iterator
  unsigned int chunkSize;
  
    //! entity type this iterator iterates over
  EntityType entType;
  
    //! dimension this iterator iterates over
  int entDimension;

    //! check for entity validity before returning handles
  bool checkValid;

};

/** \class Set-type set iterator
 * \brief A class to iterator over MOAB set-type meshsets
 */
class RangeSetIterator : public SetIterator
{
public:
  friend class Core;
  
    /** \brief Destructor
     */
  virtual ~RangeSetIterator();
  
    /** \brief get the next chunkSize entities
     * Return the next chunkSize entities.  
     * \param arr Array of entities returned.
     * \param atend Returns true if iterator is at the end of iterable values, otherwise false
     */
  virtual ErrorCode get_next_arr(std::vector<EntityHandle> &arr,
                                 bool &atend);
  
    //! reset the iterator to the beginning of the set
  virtual ErrorCode reset();
  
protected:
    /** \brief Constructor
     * \param core MOAB Core instance
     * \param ent_set EntitySet to which this iterator corresponds
     * \param chunk_size Chunk size of this iterator
     * \param ent_type Entity type for this iterator
     * \param ent_dim Entity dimension for this iterator
     */
  RangeSetIterator(Core *core, EntityHandle ent_set, int chunk_size, 
                   EntityType ent_type, int ent_dimension, bool check_valid = false);
  
private:
  ErrorCode get_next_by_type(const EntityHandle *&ptr, int count,
                             std::vector<EntityHandle> &arr, bool &atend);
  
  ErrorCode get_next_by_dimension(const EntityHandle *&ptr, int count,
                                  std::vector<EntityHandle> &arr, bool &atend);
  
    //! Build the special pair vector for the root set
  ErrorCode build_pair_vec();
  
    //! Current iterator position, 0 if at beginning
  EntityHandle iterPos;

    //! Special range pair ptr for root set
  EntityHandle *pairPtr;

    //! Number of range pairs
  int numPairs;
};
    
/** \class List-type set iterator
 * \brief A class to iterator over MOAB list-type meshsets
 */
class VectorSetIterator : public SetIterator
{
public:
  friend class Core;

    /** \brief get the next chunkSize entities
     * Return the next chunkSize entities.  
     * \param arr Array of entities returned.
     * \param atend Returns true if iterator is at the end of iterable values, otherwise false
     */
  virtual ErrorCode get_next_arr(std::vector<EntityHandle> &arr,
                                 bool &atend);
  
    //! reset the iterator to the beginning of the set
  virtual ErrorCode reset();
  
    //! decrement the position by the specified number; returns MB_FAILURE if resulting index is < 0
  inline ErrorCode decrement(int num);
  
protected:
    /** \brief Constructor
     * \param core MOAB Core instance
     * \param ent_set EntitySet to which this iterator corresponds
     * \param chunk_size Chunk size of this iterator
     * \param ent_type Entity type for this iterator
     * \param ent_dim Entity dimension for this iterator
     */
  inline VectorSetIterator(Core *core, EntityHandle eset, int chunk_sz, 
                           EntityType ent_tp, int ent_dim, bool check_valid = false)
          : SetIterator(core, eset, chunk_sz, ent_tp, ent_dim, check_valid),
            iterPos(0)
      {}
  
private:
    //! Current iterator position, 0 if at beginning
  int iterPos;
  
};

inline ErrorCode VectorSetIterator::decrement(int num) 
{
  iterPos -= num; 
  return (iterPos < 0 ? MB_FAILURE : MB_SUCCESS);
}
    
} // namespace moab

#endif
