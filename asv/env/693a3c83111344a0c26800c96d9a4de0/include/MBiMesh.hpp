#ifndef MBIMESH_HPP
#define MBIMESH_HPP

#include "moab/Core.hpp"
#include <vector>
#include <algorithm>
#include <cstring>

using namespace moab;

/* map from MOAB's ErrorCode to tstt's */
extern "C" const iBase_ErrorType iBase_ERROR_MAP[MB_FAILURE+1];

class MBiMesh
{
private:
  bool haveDeletedEntities;
  bool iCreatedInterface;
  std::vector<Tag> setHandleTags, entHandleTags;
public:
  MBiMesh(moab::Interface *mbImpl = NULL);

  virtual ~MBiMesh();
  bool have_deleted_ents( bool reset ) {
    bool result = haveDeletedEntities;
    if (reset)
      haveDeletedEntities = false;
    return result;
  }

  virtual ErrorCode delete_mesh();
  virtual ErrorCode delete_entities( const EntityHandle*, const int );
  virtual ErrorCode delete_entities( const Range& );
  iBase_AdjacencyCost AdjTable[16];
  moab::Interface *mbImpl;
  int lastErrorType;
  char lastErrorDescription[120];
  
  inline void note_set_handle_tag( Tag );
  inline void note_ent_handle_tag( Tag );
  inline void note_tag_destroyed( Tag );
  inline bool is_set_handle_tag( Tag ) const;
  inline bool is_ent_handle_tag( Tag ) const;

  inline int set_last_error( int, const char* );
  inline int set_last_error( ErrorCode, const char* );
};

static inline MBiMesh *mbimeshi_instance(iMesh_Instance instance) {return reinterpret_cast<MBiMesh*>(instance);}
#define MBIMESHI mbimeshi_instance(instance)
#define MOABI MBIMESHI->mbImpl

inline MBiMesh::MBiMesh(Interface *impl)
  : haveDeletedEntities(false), iCreatedInterface(false), mbImpl(impl),
    lastErrorType(iBase_SUCCESS)
{
  lastErrorDescription[0] = '\0';

  iBase_AdjacencyCost tmp_table[] = {
      iBase_ALL_ORDER_1, iBase_SOME_ORDER_1,    iBase_SOME_ORDER_1,    iBase_ALL_ORDER_1,
      iBase_ALL_ORDER_1, iBase_UNAVAILABLE,     iBase_SOME_ORDER_LOGN, iBase_SOME_ORDER_LOGN,
      iBase_ALL_ORDER_1, iBase_SOME_ORDER_LOGN, iBase_UNAVAILABLE,     iBase_SOME_ORDER_LOGN,
      iBase_ALL_ORDER_1, iBase_SOME_ORDER_LOGN, iBase_SOME_ORDER_LOGN, iBase_ALL_ORDER_1
  };
  memcpy(AdjTable, tmp_table, 16*sizeof(iBase_AdjacencyCost));

  if (!mbImpl) {
    mbImpl = new Core();
    iCreatedInterface = true;
  }
}

inline MBiMesh::~MBiMesh() 
{
  if (iCreatedInterface) delete mbImpl;
}

inline ErrorCode MBiMesh::delete_mesh() {
  haveDeletedEntities = true;
  return mbImpl->delete_mesh();
}

inline ErrorCode MBiMesh::delete_entities( const EntityHandle* a, const int n )
{
  if (n > 0)
    haveDeletedEntities = true;
  return mbImpl->delete_entities( a, n );
}

inline ErrorCode MBiMesh::delete_entities( const Range& r )
{
  if (!r.empty())
    haveDeletedEntities = true;
  return mbImpl->delete_entities( r );
}


void MBiMesh::note_set_handle_tag( Tag t )
{
  std::vector<Tag>::iterator i;
  i = std::lower_bound( entHandleTags.begin(), entHandleTags.end(), t );
  if (i != entHandleTags.end() && *i == t)
    entHandleTags.erase(i);
  i = std::lower_bound( setHandleTags.begin(), setHandleTags.end(), t );
  if (i == setHandleTags.end() || *i != t)
    setHandleTags.insert( i, t );
}

void MBiMesh::note_ent_handle_tag( Tag t )
{
  std::vector<Tag>::iterator i;
  i = std::lower_bound( setHandleTags.begin(), setHandleTags.end(), t );
  if (i != setHandleTags.end() && *i == t)
    setHandleTags.erase(i);
  i = std::lower_bound( entHandleTags.begin(), entHandleTags.end(), t );
  if (i == entHandleTags.end() || *i != t)
    entHandleTags.insert( i, t );
}

void MBiMesh::note_tag_destroyed( Tag t )
{
  std::vector<Tag>::iterator i;
  i = std::lower_bound( setHandleTags.begin(), setHandleTags.end(), t );
  if (i != setHandleTags.end() && *i == t)
    setHandleTags.erase(i);
  i = std::lower_bound( entHandleTags.begin(), entHandleTags.end(), t );
  if (i != entHandleTags.end() && *i == t)
    entHandleTags.erase(i);
}

bool MBiMesh::is_set_handle_tag( Tag t ) const
{
  return std::binary_search( setHandleTags.begin(), setHandleTags.end(), t );
}

bool MBiMesh::is_ent_handle_tag( Tag t ) const
{
  return std::binary_search( entHandleTags.begin(), entHandleTags.end(), t );
}

int MBiMesh::set_last_error( int code, const char* msg )
{
  std::strncpy( lastErrorDescription, msg, sizeof(lastErrorDescription) );
  lastErrorDescription[sizeof(lastErrorDescription)-1] = '\0';
  return (lastErrorType = static_cast<iBase_ErrorType>(code));
}

int MBiMesh::set_last_error( ErrorCode code, const char* msg )
{
  std::string message(msg);
  message += "  (MOAB Error Code: ";
  message += mbImpl->get_error_string(code);
  message += ")";
  return set_last_error( iBase_ERROR_MAP[code], message.c_str() );
}

#endif
