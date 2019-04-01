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
 
#ifndef MB_TAG_CONVENTIONS_HPP
#define MB_TAG_CONVENTIONS_HPP

//! Conventional tag names used for some often-used sets

/* MATERIAL_SET_TAG_NAME tag:
 * Represents sets of elements having a common material (corresponds to
 * element blocks in ExodusII)
 * size = sizeof(int)
 * type = int
 * value = integer id for this set (block id from ExodusII)
 * default value = -1
 */
#define MATERIAL_SET_TAG_NAME  "MATERIAL_SET"

/* DIRICHLET_SET_TAG_NAME tag:
 * Represents dirichlet-type boundary condition, usually contains only mesh vertices
 * (corresponds to nodesets in ExodusII)
 * size = sizeof(int)
 * type = int
 * value = integer id for this set (nodeset id from ExodusII)
 * default value = -1
 */
#define DIRICHLET_SET_TAG_NAME "DIRICHLET_SET"

/* NEUMANN_SET_TAG_NAME  tag:
 * Represents neumann-type boundary condition, usually contains elements with dimension
 * one lower than those found in material sets (i.e. edges in FE quad/tri models, quads/tris
 * in FE hex/tet models) (corresponds to sidesets in ExodusII)
 * size = sizeof(int)
 * type = int
 * value = integer id for this set (sideset id from ExodusII)
 * default value = -1
 */
#define NEUMANN_SET_TAG_NAME   "NEUMANN_SET"

/* HAS_MID_NODES_TAG_NAM tag:
 * Flags telling whether elements in a given set have mid-(edge, face, region) vertices/nodes;
 * index 0 is a place holder, so this datum can be indexed by dimension, e.g. has_mid_nodes[dim]
 * indicates whether mesh entities of dimension dim have mid nodes
 * size = 4*sizeof(int)
 * type = int[4]
 * value = 1 (has mid nodes), 0 (does not have mid nodes)
 * default value = [-1, -1, -1, -1]
 */
#define HAS_MID_NODES_TAG_NAME "HAS_MID_NODES"

/* GEOM_DIMENSION tag: 
 * Represents entities "owned" by a given topological entity in a geometric model
 * size = sizeof(int)
 * type = int
 * value = dimension of geom entity 
 * default value = -1
 */
#define GEOM_DIMENSION_TAG_NAME "GEOM_DIMENSION"

/* MESH_TRANSFORM tag:
 * Represents homogeneous transform to be applied to mesh; used in ExodusII writer to apply
 * transform before writing nodal coordinates
 * size = 16*sizeof(double)
 * type = double[16]
 * value = 4x4 homogenous transform matrix
 */
#define MESH_TRANSFORM_TAG_NAME "MESH_TRANSFORM"

/* GLOBAL_ID tag:
 * Represents global id of entities (sets or mesh entities); this id is different than the id
 * embedded in the entity handle
 * size = sizeof(int)
 * type = int
 * value = global id
 * default value = 0 // not -1 to allow gids stored in unsigned data types
 */
#define GLOBAL_ID_TAG_NAME "GLOBAL_ID"

/* CATEGORY tag:
 * String name indicating generic "category" if the entity to which it is assigned (usually
 * sets); used e.g. to indicate a set represents geometric vertex/edge/face/region, 
 * dual surface/curve, etc.
 * size = CATEGORY_TAG_NAME_LENGTH (defined below)
 * type = char[CATEGORY_TAG_NAME_LENGTH]
 * value = NULL-terminated string denoting category name
 */
#define CATEGORY_TAG_NAME "CATEGORY"
#define CATEGORY_TAG_SIZE 32

/* NAME tag:
 * A fixed length NULL-padded string containing a name.
 * All values should be assumed to be of type char[NAME_TAG_SIZE].
 * The string need not be null terminated.  All values used for
 * storing or searching for a value must be padded with '\0' chars.
 */
#define NAME_TAG_NAME "NAME"
#define NAME_TAG_SIZE 32

/* BLOCK_HEADER: tag
 * A fixex lenght tag containg block header data
 * BlockColor, MaterialId and BlockDimension
 */
#define BLOCK_HEADER "BLOCK_HEADER"

/* BLOCK_ATTRIBUTES: tag
 * A varible lenght tag of doubles
 * Tag contains attributes set to BlockSet in cubit file
 */
#define BLOCK_ATTRIBUTES "BLOCK_ATTRIBUTES"

#ifndef MB_PARALLEL_CONVENTIONS_H
#define MB_PARALLEL_CONVENTIONS_H

/** Tag conventions for naming parallel things.  Note this header
 * file belongs in the main MOAB directory because even serial
 * applications (e.g. partitioners) may write tags for use in
 * parallel applications.
 */

/** \brief Global identifier for interface mesh
 *
 * An integer identifier common to the corresponding mesh entity
 * instances on each processor for a mesh entity on the interface.
 */
#define PARALLEL_GID_TAG_NAME "GLOBAL_ID"

/** \brief Tag on a meshset representing a parallel partition.
 *
 * When the mesh is partitioned for use in a parallel environment,
 * the each CPUs partiiton of the mesh is stored in a meshset with
 * this tag.  The value of the tag is an integer "part identifier".
 */
#define PARALLEL_PARTITION_TAG_NAME "PARALLEL_PARTITION"
#define PARALLEL_PART_TAG_NAME PARALLEL_PARTITION_TAG_NAME

/** \brief Tag that groups the set of parts/partitions that are
 *         a covering of the mesh.
 *
 * This tag labels an entity set for which the child sets are part(ition)s
 * that together are a single partitioning of the mesh.  I.e. There should
 * be no mesh entity that is contained in more than one child part(ition)
 * set, and typically every mesh entity of the dimenion used to partition
 * the mesh is contained in exactly one of the child sets.
 *
 * The data for this tag is a single integer value.  The value of
 * the tag is undefined.
 */
#define PARALLEL_PARITIONING_TAG_NAME "PARALLEL_MESH_PARITIONING"

/** \brief Tag storing which other processor a given entity is shared with
 *
 * This single-valued tag implies an entity is shared with one other proc
 */
#define PARALLEL_SHARED_PROC_TAG_NAME "__PARALLEL_SHARED_PROC"
 
/** \brief Tag storing which other processorS a given entity is shared with
 *
 * This multiple-valued tag implies an entity is shared with multiple
 * other processors.  Length of tag is application-dependent, and depends on
 * what the maximum number of processors is which share an entity
 */
#define PARALLEL_SHARED_PROCS_TAG_NAME "__PARALLEL_SHARED_PROCS"
 
/** \brief Tag storing the handle of a shared entity on the other proc
 *
 * This single-valued tag implies an entity is shared with one other proc
 */
#define PARALLEL_SHARED_HANDLE_TAG_NAME "__PARALLEL_SHARED_HANDLE"
 
/** \brief Tag storing handles of a shared entity on other processors
 *
 * This multiple-valued tag implies an entity is shared with multiple
 * other processors.  Length of tag is application-dependent, and depends on
 * what the maximum number of processors is which share an entity
 */
#define PARALLEL_SHARED_HANDLES_TAG_NAME "__PARALLEL_SHARED_HANDLES"
 
/** \brief Tag storing parallel status (as bits in this tag)
 *
 * This tag stores various aspects of parallel status in bits; see also 
 * #define's following, to be used in bit mask operations.  If an entity is
 * not shared with any other processors, the pstatus is 0, otherwise it's > 0
 *
 * bit 0: !owned (0=owned, 1=not owned)
 * bit 1: shared (0=not shared, 1=shared)
 * bit 2: multishared (shared by > 2 procs; 0=not shared, 1=shared)
 * bit 3: interface (0=not interface, 1=interface)
 * bit 4: ghost (0=not ghost, 1=ghost)
 * default value = 0
 */
#define PARALLEL_STATUS_TAG_NAME "__PARALLEL_STATUS"

#define PSTATUS_NOT_OWNED 0x1
#define PSTATUS_SHARED 0x2
#define PSTATUS_MULTISHARED 0x4
#define PSTATUS_INTERFACE 0x8
// note, these numbers are in hex, so 0x10 is the 4th bit, or 2^4.
#define PSTATUS_GHOST 0x10

#define PSTATUS_AND 0x1
#define PSTATUS_OR 0x2
#define PSTATUS_NOT 0x3
#endif

#endif
