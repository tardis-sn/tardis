/* src/moab/EntityHandle.hpp.  Generated from EntityHandle.hpp.in by configure.  */
#ifndef MOAB_ENTITY_HANDLE_HPP
#define MOAB_ENTITY_HANDLE_HPP

#include "moab/MOABConfig.h"

#ifdef MOAB_HAVE_INTTYPES_H
# include <inttypes.h>
#elif defined (MOAB_HAVE_STDINT_H)
# include <stdint.h>
#elif defined (_MSC_VER)
  typedef __int8 int8_t;
  typedef __int16 int16_t;
  typedef __int32 int32_t;
  typedef __int64 int64_t;
  typedef unsigned __int8 uint8_t;
  typedef unsigned __int16 uint16_t;
  typedef unsigned __int32 uint32_t;
  typedef unsigned __int64 uint64_t;
#endif

#ifdef MOAB_HAVE_STDDEF_H
# include <stddef.h>
#elif defined (MOAB_HAVE_STDLIB_H)
# include <stdlib.h>
#elif defined (MOAB_HAVE_SYS_TYPES_H)
# include <sys/types.h>
#endif

namespace moab {

#ifdef MOAB_FORCE_64_BIT_HANDLES
  typedef uint64_t EntityHandle;
  typedef  int64_t EntityID;
#elif defined (MOAB_FORCE_32_BIT_HANDLES)
  typedef uint32_t EntityHandle;
  typedef  int32_t EntityID;
#else
# ifdef MOAB_HAVE_SIZE_T
    typedef size_t EntityHandle;
# else
    typedef unsigned long EntityHandle;
# endif
# ifdef MOAB_HAVE_PTRDIFF_T
    typedef ptrdiff_t EntityID;
# else
    typedef long EntityID;
# endif
#endif 

} // namespace moab

#endif
