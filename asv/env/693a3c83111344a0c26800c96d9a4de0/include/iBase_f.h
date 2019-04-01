#ifndef IBASE_F_H
#define IBASE_F_H

#ifdef POINTER_SIZE
!  user-defined override
#  define IBASE_HANDLE_T INTEGER*POINTER_SIZE
#elif defined(__SIZEOF_SIZE_T__)
!  gfortran
#  define IBASE_HANDLE_T integer*__SIZEOF_SIZE_T__
#elif defined(__GNUC__) && defined (_LANGUAGE_FORTRAN)
!  gfortran 4.2 and earlier
#  ifdef _LP64
#    define IBASE_HANDLE_T integer*8
#  else
!    NOTE: not integer*4 (see meaning of _LP64 in GNU docs)
#    define IBASE_HANDLE_T integer
#  endif
#elif defined(__INTEL_COMPILER)
#  ifdef __x86_64
#    define IBASE_HANDLE_T integer*8
#  else
#    define IBASE_HANDLE_T integer*4
#  endif
#elif defined(__XLCPP__) || defined(__bg__) || defined(__bgp__) || defined(__bgq__)
!  __XLCPP__ doesn't work for IBM.  AT least make things work for BGP and BGQ
!        USE, INTRINSIC :: ISO_C_BINDING
!#  define IBASE_HANDLE_T TYPE(C_PTR)
#  define IBASE_HANDLE_T integer*8
#else
#  WARNING "Unknown compiler; using SIZEOF_VOID_P."
#  define IBASE_HANDLE_T integer*8  
#endif

#define iBase_EntityHandle IBASE_HANDLE_T
#define iBase_EntitySetHandle IBASE_HANDLE_T
#define iBase_TagHandle IBASE_HANDLE_T
#define iBase_EntityIterator IBASE_HANDLE_T
#define iBase_EntityArrIterator IBASE_HANDLE_T

#endif

      integer iBase_VERTEX
      integer iBase_EDGE
      integer iBase_FACE
      integer iBase_REGION
      integer iBase_ALL_TYPES

      parameter (iBase_VERTEX = 0) 
      parameter (iBase_EDGE = 1) 
      parameter (iBase_FACE = 2) 
      parameter (iBase_REGION = 3)
      parameter (iBase_ALL_TYPES = 4)

      integer iBase_AdjacencyCost_MIN
      integer iBase_UNAVAILABLE
      integer iBase_ALL_ORDER_1
      integer iBase_ALL_ORDER_LOGN
      integer iBase_ALL_ORDER_N
      integer iBase_SOME_ORDER_1
      integer iBase_SOME_ORDER_LOGN
      integer iBase_SOME_ORDER_N
      integer iBase_AVAILABLE
      integer iBase_AdjacencyCost_MAX

      parameter (iBase_AdjacencyCost_MIN = 0)
      parameter (iBase_UNAVAILABLE = 0)
      parameter (iBase_ALL_ORDER_1 = 1)
      parameter (iBase_ALL_ORDER_LOGN = 2)
      parameter (iBase_ALL_ORDER_N = 3)
      parameter (iBase_SOME_ORDER_1 = 4)
      parameter (iBase_SOME_ORDER_LOGN = 5)
      parameter (iBase_SOME_ORDER_N = 6)
      parameter (iBase_AVAILABLE = 7)
      parameter (iBase_AdjacencyCost_MAX = 7)

      integer iBase_NEW
      integer iBase_ALREADY_EXISTED
      integer iBase_CREATED_DUPLICATE
      integer iBase_CREATION_FAILED

      parameter (iBase_NEW = 0)               
      parameter (iBase_ALREADY_EXISTED = 1)   
      parameter (iBase_CREATED_DUPLICATE = 2)
      parameter (iBase_CREATION_FAILED = 3)


      integer iBase_SILENT
      integer iBase_WARN_ONLY
      integer iBase_THROW_ERROR

      parameter (iBase_SILENT = 0) 
      parameter (iBase_WARN_ONLY = 1) 
      parameter (iBase_THROW_ERROR = 4)


      integer iBase_SUCCESS
      integer iBase_MESH_ALREADY_LOADED
      integer iBase_FILE_NOT_FOUND
      integer iBase_FILE_WRITE_ERROR
      integer iBase_NIL_ARRAY
      integer iBase_BAD_ARRAY_SIZE
      integer iBase_BAD_ARRAY_DIMENSION
      integer iBase_INVALID_ENTITY_HANDLE
      integer iBase_INVALID_ENTITY_COUNT
      integer iBase_INVALID_ENTITY_TYPE
      integer iBase_INVALID_ENTITY_TOPOLOGY
      integer iBase_BAD_TYPE_AND_TOPO
      integer iBase_ENTITY_CREATION_ERROR
      integer iBase_INVALID_TAG_HANDLE
      integer iBase_TAG_NOT_FOUND
      integer iBase_TAG_ALREADY_EXISTS
      integer iBase_TAG_IN_USE
      integer iBase_INVALID_ENTITYSET_HANDLE
      integer iBase_INVALID_ITERATOR_HANDLE
      integer iBase_INVALID_ARGUMENT
      integer iBase_MEMORY_ALLOCATION_FAILED
      integer iBase_NOT_SUPPORTED
      integer iBase_FAILURE

      parameter (iBase_SUCCESS = 0)
      parameter (iBase_MESH_ALREADY_LOADED = 1)
      parameter (iBase_FILE_NOT_FOUND = 2)
      parameter (iBase_FILE_WRITE_ERROR = 3)
      parameter (iBase_NIL_ARRAY = 4)
      parameter (iBase_BAD_ARRAY_SIZE = 5)
      parameter (iBase_BAD_ARRAY_DIMENSION = 6)
      parameter (iBase_INVALID_ENTITY_HANDLE = 7)
      parameter (iBase_INVALID_ENTITY_COUNT = 8)
      parameter (iBase_INVALID_ENTITY_TYPE = 9)
      parameter (iBase_INVALID_ENTITY_TOPOLOGY = 10)
      parameter (iBase_BAD_TYPE_AND_TOPO = 11)
      parameter (iBase_ENTITY_CREATION_ERROR = 12)
      parameter (iBase_INVALID_TAG_HANDLE = 13)
      parameter (iBase_TAG_NOT_FOUND = 14)
      parameter (iBase_TAG_ALREADY_EXISTS = 15)
      parameter (iBase_TAG_IN_USE = 16)
      parameter (iBase_INVALID_ENTITYSET_HANDLE = 17)
      parameter (iBase_INVALID_ITERATOR_HANDLE = 18)
      parameter (iBase_INVALID_ARGUMENT = 19)
      parameter (iBase_MEMORY_ALLOCATION_FAILED = 20)
      parameter (iBase_NOT_SUPPORTED = 21)
      parameter (iBase_FAILURE = 22)


      integer iBase_BLOCKED
      integer iBase_INTERLEAVED
      integer iBase_UNDETERMINED

      parameter (iBase_BLOCKED = 0)
      parameter (iBase_INTERLEAVED = 1)
      parameter (iBase_UNDETERMINED = 2)


      integer iBase_INTEGER
      integer iBase_DOUBLE
      integer iBase_ENTITY_HANDLE
      integer iBase_ENTITY_SET_HANDLE
      integer iBase_BYTES

      parameter (iBase_BYTES = 0)
      parameter (iBase_INTEGER = 1)
      parameter (iBase_DOUBLE = 2)
      parameter (iBase_ENTITY_HANDLE = 3)
      parameter (iBase_ENTITY_SET_HANDLE = 3)

