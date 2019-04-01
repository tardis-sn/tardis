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

#ifndef MB_ENTITY_TYPE_H
#define MB_ENTITY_TYPE_H

/* This file can be used to define several different things.
 *
 * A) If included in C code (not C++), it defines:
 *    1) An enum named MBEntityType, guarded by the MB_ENTITY_TYPE_H
 *        include guards and 
 *    2) a typedef MBEntityType guarded by MOAB_ENTITY_TYPE_C include guards.
 *
 * B) If included in C++ code, it defines:
 *    1) An enum named EntiyType in the MOAB namespace, guarded
 *       by the MB_ENTITY_TYPE include guards
 *    2) Increment and decrement oeprators for the moab::EntityType enum,
 *       also guarded by the MB_ENTITY_TYPE include guards
 *    3) A typedef for moab::EntityType in the global namespace 
 *        named MBEntityType, guarded by the MOAB_ENTITY_TYPE_NS_ONLY
 *        include guards
 *
 * The C and C++ code should be entirely independent.  They are defined
 * in the same file only to avoid code duplication and inconsistent enum
 * values.  OTOH, the C++ definitions must be in the same file because
 * the compiler must treat both the namespaced and non-namespaced names
 * as the same type.
 *
 * The C++ code must be able to provide:
 *  a) An enum in the moab namespace
 *  b) An enum in the global namespace that is the *same type*
 *      as a) as far as the compiler is concerned.
 *  c) Nothing in the global namespace unless requested
 *  d) No breakage if both namespaced and non-namespaced headers
 *      are both included.
 * 
 * This is acheived with the somewhat complicated set of multiple
 * included guards described above, where moab/EntityType.hpp will
 * include this file with MOAB_ENTITY_TYPE_NS_OLNY temporarily defined
 * so as to pull in only the namespaced version at that time, without
 * prohibiting the non-namespaced version from being pulled in previously
 * or later.
 */
#ifdef __cplusplus
namespace moab { 
# define MOAB_ENTITY_TYPE_NAME EntityType
# else /* __cplusplus */
# define MOAB_ENTITY_TYPE_NAME MBEntityType
#endif /* __cplusplus */

/*! Entity types defined in MOAB and MBCN
 *  The ordering here must ensure that all element types are 
 *  grouped together and all elements of similar dimension are
 *  grouped together.
 */
enum MOAB_ENTITY_TYPE_NAME
{
  MBVERTEX = 0, /**< Mesh Vertex AKA node */
  MBEDGE,       /**< Mesh Edge */
  MBTRI,        /**< Triangular element (including shells) */
  MBQUAD,       /**< Quadrilateral element (including shells) */
  MBPOLYGON,    /**< Polygon */
  MBTET,        /**< Tetrahedral element */
  MBPYRAMID,    /**< Pyramid element (where are the face ids for this defined?) */
  MBPRISM,      /**< Wedge element (Exodus has one, Cubit doesn't. Does Mesh need it?) */
  MBKNIFE,      /**< Knife element */
  MBHEX,        /**< Hexahedral element */
  MBPOLYHEDRON, /**< Polyhedron */
  MBENTITYSET,    /**< MeshSet */
  MBMAXTYPE  /**< Just a place keeper - must be the # of entities, for array */
    /**< dimensioning purposes  */
};

#ifdef __cplusplus
/** prefix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME & operator++(MOAB_ENTITY_TYPE_NAME &type)
{
  return type = static_cast<MOAB_ENTITY_TYPE_NAME>(type+1);
}

/** postfix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME operator++(MOAB_ENTITY_TYPE_NAME &type, int)
{
  MOAB_ENTITY_TYPE_NAME oldval = type;
  ++type;
  return oldval;
}

/** prefix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME & operator--(MOAB_ENTITY_TYPE_NAME &type)
{
  return type = static_cast<MOAB_ENTITY_TYPE_NAME>(type-1);
}

/** postfix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME operator--(MOAB_ENTITY_TYPE_NAME &type, int)
{
  MOAB_ENTITY_TYPE_NAME oldval = type;
  --type;
  return oldval;
}

} /* namespace moab*/
#endif /* __cplusplus */

#undef MOAB_ENTITY_TYPE_NAME
#endif /* MB_ENTITY_TYPE_H */

#ifdef __cplusplus
#  ifndef MOAB_ENTITY_TYPE_NS_ONLY
#    define MOAB_ENTITY_TYPE_NS_ONLY
     typedef moab::EntityType MBEntityType;
     using moab::MBVERTEX;
     using moab::MBEDGE;
     using moab::MBTRI;
     using moab::MBQUAD;
     using moab::MBPOLYGON;
     using moab::MBTET;
     using moab::MBPYRAMID;
     using moab::MBPRISM;
     using moab::MBKNIFE;
     using moab::MBHEX;
     using moab::MBPOLYHEDRON;
     using moab::MBENTITYSET;
     using moab::MBMAXTYPE;
#  endif /* MOAB_ENTITY_TYPE_NS_ONLY */
#else /* __cplusplus */
#  ifndef MOAB_ENTITY_TYPE_C
#    define MOAB_ENTITY_TYPE_C
     typedef enum MBEntityType MBEntityType;
#  endif /* MOAB_ENTITY_TYPE_C */
#endif /* __cplusplus */
