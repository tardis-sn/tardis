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

#ifndef MOAB_TYPES_HPP
#define MOAB_TYPES_HPP

#ifdef __cplusplus
#include "moab/EntityType.hpp"
#include "moab/EntityHandle.hpp"
#endif

/**\name Types and names
 * Types used in the MOAB interface
 */
/*@{*/

#ifdef __cplusplus
namespace moab {
#endif

/** MOAB error codes */
enum ErrorCode { MB_SUCCESS = 0,
                   MB_INDEX_OUT_OF_RANGE,
                   MB_TYPE_OUT_OF_RANGE,
                   MB_MEMORY_ALLOCATION_FAILED,
                   MB_ENTITY_NOT_FOUND,
                   MB_MULTIPLE_ENTITIES_FOUND,
                   MB_TAG_NOT_FOUND,
                   MB_FILE_DOES_NOT_EXIST,
                   MB_FILE_WRITE_ERROR,
                   MB_NOT_IMPLEMENTED,
                   MB_ALREADY_ALLOCATED,
                   MB_VARIABLE_DATA_LENGTH,
                   MB_INVALID_SIZE,
                   MB_UNSUPPORTED_OPERATION,
                   MB_UNHANDLED_OPTION,
                   MB_STRUCTURED_MESH,
                   MB_FAILURE};

#ifdef __cplusplus
extern const char* const ErrorCodeStr[];
#endif

/** Misc. integer constants, declared in enum for portability */
enum Constants {
  MB_VARIABLE_LENGTH = -1 /**< Length value for variable-length tags */ 
};

/** Specify storage type for tags.  See MOAB users guide for more information. */
enum TagType {
    MB_TAG_BIT   = 0,    /**< Tag size specified in bits, tag value is 8 bits or less */
  MB_TAG_SPARSE= 1<<0, /**< Storage optimized for tags on a few entities */
  MB_TAG_DENSE = 1<<1, /**< Storage optimized for tags on most entities of a type */
  MB_TAG_MESH  = 1<<2, /**< Storage for tags on no entities, only the root set/whole mesh. */
  MB_TAG_BYTES = 1<<3, /**< Size is in number of bytes rather than number of values of \c DataType */
  MB_TAG_VARLEN= 1<<4, /**< Create variable-length tag */
  MB_TAG_CREAT = 1<<5, /**< Create tag if it does not already exist */
  MB_TAG_EXCL  = 1<<6, /**< Fail if TAG_CREATE and tag already exists */
  MB_TAG_STORE = 1<<7, /**< Fail if tag exists and has different storage type */
  MB_TAG_ANY   = 1<<8, /**< Do not fail if size, type, or default value do not match. */
  MB_TAG_NOOPQ = 1<<9, /**< Do not accept MB_TYPE_OPAQUE as a match for any type. */
  MB_TAG_DFTOK = 1<<10 /**< Do not fail for mismatched default values */
/**<  MB_TAG_CNVRT = 1<<11,  Convert storage type if it does not match */
};

/** Specify data type for tags. */
enum DataType {
  MB_TYPE_OPAQUE  = 0, /**< byte array */
  MB_TYPE_INTEGER = 1, /**< native 'int' type */
  MB_TYPE_DOUBLE  = 2, /**< native 'double' type */
  MB_TYPE_BIT     = 3, /**< mandatory type for tags with MB_TAG_BIT storage */
  MB_TYPE_HANDLE  = 4, /**< EntityHandle */
  MB_MAX_DATA_TYPE = MB_TYPE_HANDLE
};

#ifdef __cplusplus
extern const char* const DataTypeStr[];
#endif

/** Used to reference tags; since they're so different from entities, we
 *  use void** instead of a uint to prevent them from being confused as 
 *  entity handles.
 */
#ifdef __cplusplus
class TagInfo;
typedef TagInfo* Tag;
#else
struct TagInfo;
typedef struct TagInfo* Tag;
#endif

/** Meshset options: properties for meshset creation.
 *  Values are bit flags that may be combined with a bitwise OR (|)
 */
enum EntitySetProperty {
  MESHSET_TRACK_OWNER = 0x1, /**< create entity to meshset adjacencies */
  MESHSET_SET         = 0x2, /**< set contents are unique */
  MESHSET_ORDERED     = 0x4  /**< order of set contents is preserved */
};

enum SenseType {
  SENSE_INVALID      = -2, /**< default, invalid, not defined */
  SENSE_REVERSE      = -1, /**< reversed */
  SENSE_BOTH         =  0, /**< both senses valid  */
  SENSE_FORWARD      =  1  /**< forward  */
};

#ifdef __cplusplus
extern const char* const* const SenseTypeStr;
#endif

#ifdef __cplusplus
} /* namespace moab */
#endif

/*@}*/

#endif
