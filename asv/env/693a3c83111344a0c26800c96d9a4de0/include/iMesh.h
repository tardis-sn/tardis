#ifndef _ITAPS_iMesh
#define _ITAPS_iMesh

#include "iBase.h"
#include "iMesh_protos.h"

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compile time version number digits
 *
 * iMesh maintains a major, minor and patch digit in its version number.
 * Technically speaking, there is not much practical value in patch digit
 * for an interface specification. A patch release is typically only used
 * for bug fix releases. Although it is rare, sometimes a bug fix
 * necessitates an API change. So, we define a patch digit for iMesh.
 ******************************************************************************/
#define IMESH_VERSION_MAJOR ITAPS_VERSION_MAJOR
#define IMESH_VERSION_MINOR ITAPS_VERSION_MINOR
#define IMESH_VERSION_PATCH ITAPS_VERSION_PATCH

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Maintain backward compatibility with old version symbol names
 ******************************************************************************/
#define IMESH_MAJOR_VERSION IMESH_VERSION_MAJOR
#define IMESH_MINOR_VERSION IMESH_VERSION_MINOR
#define IMESH_PATCH_VERSION IMESH_VERSION_PATCH

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Version Comparison
 *
 * Evaluates to true at CPP time if the version of iMesh currently being
 * compiled is greater than or equal to the version specified.
 ******************************************************************************/
#define IMESH_VERSION_GE(Maj,Min,Pat) ITAPS_VERSION_GE(Maj,Min,Pat)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose string represention of the iMesh version number
 ******************************************************************************/
#define IMESH_VERSION_STRING ITAPS_VERSION_STRING_(iMesh)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose a symbol name derived from the current iMesh version number.
 ******************************************************************************/
#define IMESH_VERSION_TAG ITAPS_VERSION_TAG_(iMesh)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Define iMesh_newMesh symbol such that it depends on version number.
 *
 * Note: We ran into problems with this as it influences or is influenced by
 * fortran name mangling and so breaks fortran compilation. So, this is
 * currently disabled.
 ******************************************************************************/
#define IMESH_NEW_MESH_NAME__(A,B,C) A##_##B##_##C
#define IMESH_NEW_MESH_NAME_(A,B,C) IMESH_NEW_MESH_NAME__(A,B,C)
#define IMESH_NEW_MESH_NAME(A) IMESH_NEW_MESH_NAME_(A,IMESH_VERSION_MAJOR,IMESH_VERSION_MINOR)
/*
#undef  iMesh_newMesh
#define iMesh_newMesh IMESH_NEW_MESH_NAME(iMesh_newMesh)
*/

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 * \ingroup Datatypes
 * \brief iMesh instance
 ******************************************************************************/
typedef struct iMesh_Instance_Private* iMesh_Instance;

/***************************************************************************//**
 * \ingroup Datatypes
 * \brief Entity Topology 
 ******************************************************************************/
enum iMesh_EntityTopology {
    iMesh_EntityTopology_MIN = 0,
        /**< MIN symbol used to facilitate iteration over topologies */
    iMesh_POINT = iMesh_EntityTopology_MIN,
        /**< a 0D entity (e.g. a vertex) */
    iMesh_LINE_SEGMENT,
        /**< a 1D entity (e.g. an edge) */
    iMesh_POLYGON,
        /**< a general 2D entity */
    iMesh_TRIANGLE,
        /**< a specific 2D entity bounded by 3 edge entities */
    iMesh_QUADRILATERAL,
        /**< a specific 2D entity bounded by 4 edge entities */
    iMesh_POLYHEDRON,
        /**< a general 3D entity */
    iMesh_TETRAHEDRON,
        /**< a specific 3D entity bounded by 4 triangle entities */
    iMesh_HEXAHEDRON,
        /**< a specific 3D entity bounded by 6 quadrilateral entities */
    iMesh_PRISM,
        /**< a specific 3D entity bounded by a combination of 3 quadrilateral
            entities and 2 triangle entities */
    iMesh_PYRAMID,
        /**< a specific 3D entity bounded by a combination of 1 quadrilateral
             entity and 4 triangle entities */
    iMesh_SEPTAHEDRON,
        /**< a hexahedral entity with one collapsed edge */
    iMesh_ALL_TOPOLOGIES,
        /**< used only in queries to request information about all topologies */
    iMesh_EntityTopology_MAX = iMesh_ALL_TOPOLOGIES
        /**< MAX symbol used to facilitate iteration over topologies */
};

/***************************************************************************//**
 * \ingroup ErrorHandling
 * \brief  Get the error type returned from the last iMesh function
 ******************************************************************************/
void iMesh_getErrorType(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int* error_type  
        /**< [out] Error type returned from last iMesh function
             (see iBase_ErrorType) */

);

/***************************************************************************//**
 * \ingroup ErrorHandling
 * \brief Get a description of the error returned from the last iMesh function
 ******************************************************************************/
void iMesh_getDescription(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    char* descr, 
        /**< [in,out] Pointer to a character string to be filled with a
           description of the error from the last iMesh function */
    int descr_len  
        /**< [in] Length of the character string pointed to by descr
             (\ref strlen) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Construct a new iMesh instance
 *
 ******************************************************************************/

void iMesh_newMesh(
    const char* options, 
        /**< [in] Pointer to implementation-specific options string
             (\ref options) */
    iMesh_Instance* instance, 
        /**< [in] iMesh instance handle */
    int* err, 
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int options_len  
        /**< [in] Length of the character string pointed to by options
             (\ref strlen) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Destroy an iMesh instance
 *
 ******************************************************************************/

void iMesh_dtor(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Load a mesh from a file
 *
 * Load a mesh from a file.  If entity set is specified, loaded mesh
 * is added to that set; specify root set if that is not desired.
 ******************************************************************************/

void iMesh_load(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Set to which loaded mesh will be added, otherwise root */
    const char* name, 
        /**< [in] File name from which mesh is to be loaded */
    const char* options, 
        /**< [in] Pointer to implementation-specific options string
             (\ref options) */
    int* err, 
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int name_len, 
        /**< [in] Length of the file name character string
             (\ref strlen) */
    int options_len  
        /**< [in] Length of the options character string
             (\ref strlen) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Save a mesh to a file
 *
 * Save a mesh to a file.  If entity set is specified, save only the
 * mesh contained in that set.
 ******************************************************************************/

void iMesh_save(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being saved */
    const char* name, 
        /**< [in] File name to which mesh is to be saved */
    const char* options, 
        /**< [in] Pointer to implementation-specific options string
             (\ref options) */
    int* err, 
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int name_len, 
        /**< [in] Length of the file name character string
             (\ref strlen) */
    int options_len  
        /**< [in] Length of the options character string
             (\ref strlen) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Get handle of the root set for this instance
 *
 * Get handle of the root set for this instance.  All mesh entities in
 * this instance can be accessed from this set.
 ******************************************************************************/

void iMesh_getRootSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle* root_set, 
        /**< [out] Pointer to set handle returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Get the geometric dimension of mesh represented in this instance
 *
 ******************************************************************************/

void iMesh_getGeometricDimension(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int* geom_dim, 
        /**< [out] Pointer to dimension returned from this function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Set geometric dimension of vertex coordinates
 *
 * Set the geometric dimension of the mesh.  Note:  An application
 * should not expect this function to succeed unless the mesh instance
 * is empty. An empty mesh instance is any mesh instance in which there are 
 * no vertex entities.
 ******************************************************************************/

void iMesh_setGeometricDimension(
     iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int geom_dim, 
        /**< [in] Requested geometric dimension. */
    int* err   
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Initialization
 * \brief  Get the default storage order used by this implementation
 *
 * Get the default storage order used by this implementation.  Value
 * returned is a member of the iBase_StorageOrder enumeration.
 ******************************************************************************/

void iMesh_getDfltStorage(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int* order, 
        /**< [out] Pointer to storage order returned from function
             (\sa iBase_StorageOrder) */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Adjacencies
 * \brief  Get the adjacency table for this implementation
 *
 * Get the adjacency table for this implementation.  This table 
 * is a 4x4 array whose entries characterize how the implementation
 * behaves when adjacent and intermediate entities are queried.
 * Entry [i,j] (i=row, j=col, 0-based indexing) will have one of
 * the values in the iBase_AdjacencyCost enumeration. Off-diagonal
 * entres, i!=j, represents the relative cost of retrieving 
 * entities of dimension i adjacent to entities of dimension j.
 * Diagonal entries, i==j, indicate whether or not handles to ALL
 * entities of that dimension are obtainable from calls that return
 * entity handles. This is always true by definition for i==j==0
 * as well as for i==j==2 in a 2D mesh and i==j==3 in a 3D mesh.
 * However, diagonal entries [1,1] for a 2D mesh and both [1,1]
 * and [2,2] for a 3D mesh indicate whether or not handles to ALL
 * intermediate dimensioned entities (ALL edges in a 2D mesh or
 * ALL edges and ALL faces in a 3D mesh) are obtainable from calls
 * that return entity handles. A value of iBase_AVAILABLE for a
 * diagonal entry indicates that handles are obtainable for ALL
 * entities of that dimension while a value of iBase_UNAVAILABLE
 * indicates that handles are not obtainable for ALL entities of that
 * dimension.
 ******************************************************************************/

void iMesh_getAdjTable (
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int** adjacency_table, 
        /**< [out] Pointer to array representing adjacency table
             \ref trio) */
    int* adjacency_table_allocated, 
        /**< [in,out] Pointer to allocated size of adjacency_table */
    int* adjacency_table_size, 
        /**< [out] Pointer to occupied size of adjacency_table */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Adjacencies
 * \brief  Set the adjacency table as requested by the application
 *
 * Set the adjacency table as requested by the application. 
 * See iMesh_getAdjTable for a description of the meaning of the entries
 * in this table. This call to set the adjacency table allows an application
 * to request how it would like an implementation to behave when adjacent
 * and/or intermediate entities are queried. If the implementation is not
 * able to accommodate the requested query behavior associated with ANY
 * entry in the table, the call will fail and return an error of
 * iBase_NOT_SUPPORTED. Otherwise, the implementation is able to accommodate
 * the requested query behavior associated with ALL table entries and the
 * call will succeed. In either case, however, the implementation will
 * over-write the entries in the adjacency_table argument with the same
 * values that would be obtained in a succeeding call to iMesh_getAdjTable.
 ******************************************************************************/

void iMesh_setAdjTable (
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int* adjacency_table, 
        /**< [in,out] Array representing adjacency table requested by
             application */
    int adjacency_table_size, 
        /**< [in] Size of adj table array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Get the number of entities of specified type in the instance or set
 *
 * Get the number of entities with the specified type in the instance 
 * or set.  If entity set handle is root set, return information for instance,
 * otherwise for set.  Value of entity type must be from the
 * iBase_EntityType enumeration.  If iBase_ALL_TYPES is specified,
 * total number of entities (excluding entity sets) is returned.
 ******************************************************************************/

void iMesh_getNumOfType(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being queried */
    const int entity_type, 
        /**< [in] Type of entity requested */
    int* num_type, 
        /**< [out] Pointer to number of entities, returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Get the number of entities of specified topology in instance or set
 *
 * Get the number of entities with the specified topology in the instance 
 * or set.  If entity set handle is root set, return information for instance,
 * otherwise for set.  Value of entity topology must be from the
 * iMesh_EntityTopology enumeration.  If iMesh_ALL_TOPOLOGIES is specified,
 * total number of entities (excluding entity sets) is returned.
 ******************************************************************************/

void iMesh_getNumOfTopo(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being queried */
    const int entity_topology, 
        /**< [in] Topology of entity requested */
    int* num_topo, 
        /**< [out] Pointer to number of entities, returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  MeshModification
 * \brief Permit implementation to 'optimize' the mesh instance 
 *
 * Its concievable that after a series of operations modifying the mesh
 * instance, the implementation's internal representation of the mesh may
 * include data and tie up memory resources that could be 'optimized' away.
 * For example, if the implementation only marks removed entities for deletion
 * but does not actually free up memory resources associated those entities,
 * a call to this function gives the implementation the 'ok' to go ahead and
 * free up such memory.
 *
 * Depending on the implementation as well as the amount and kind of changes
 * to the mesh that have been made prior to calling it, this call may be
 * very expensive in time complexity. On the other hand, it is also perfectly
 * acceptable that an implementation basically treat this operation as a no-op.
 *
 * In any event, any entity, entity set, iterator or tag handle obtained prior 
 * to calling this method will may be invalidated as a result of this call. In
 * that case, the caller must re-acquire handles after making this call. A
 * return argument indicates if handles have been invalidated.
 ******************************************************************************/

void iMesh_optimize(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    int* handles_invalidated,
        /**< [out] Returned flag indicating if any handles the caller held
             before calling this function have been invalidated as a result of
             this call. A value of non-zero returned here indicates that any
             handles the caller had prior to this call must not be trusted and
             must be re-acquired by caller. */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Get entities of specific type and/or topology in set or instance
 *
 * Get entities of specific type and/or topology in set or instance.  All 
 * entities of a given type or topology are requested by specifying
 * iBase_ALL_TOPOLOGIES or iBase_ALL_TYPES, respectively.  Specified type
 * or topology must be a value in the iBase_EntityType or iBase_EntityTopology
 * enumeration, respectively.
 ******************************************************************************/

void iMesh_getEntities(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being queried */
    const int entity_type, 
        /**< [in] Type of entities being requested */
    const int entity_topology, 
        /**< [in] Topology of entities being requested */
    iBase_EntityHandle** entity_handles, 
        /**< [in,out] Pointer to array of entity handles returned \ref trio) */
    int* entity_handles_allocated, 
        /**< [in,out] Pointer to allocated size of entity_handles array */
    int* entity_handles_size, 
        /**< [out] Pointer to occupied size of entity_handles array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  VertexEntities
 * \brief  Get coordinates of specified vertices
 *
 * Get coordinates of specified vertices. Coordinates are returned
 * in the storage order indicated by the storage_order argument.
 ******************************************************************************/

void iMesh_getVtxArrCoords(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* vertex_handles, 
        /**< [in] Array of mesh vertex handles whose coordinates are being
             requested */
    const int vertex_handles_size, 
        /**< [in] Number of vertices in vertex_handles array */
    const int storage_order, 
        /**< [in] Requested storage order of returned coordinates
             (see iBase_StorageOrder) */
    double** coords, 
        /**< [in,out] Pointer to array of coordinates returned from function
             \ref trio) */
    int* coords_allocated, 
        /**< [in,out] Pointer to allocated size of coords array */
    int* coords_size, 
        /**< [out] Pointer to occupied size of coords array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Initialize an array iterator over specified entity type, topology,
 *  and size
 *
 * Initialize an array iterator over specified entity type, topology, and 
 * size, for a specified set or instance.  Iterator returned can be used 
 * as input to functions returning entities for the iterator.  If all 
 * entities of a specified type and/or topology are to be iterated, 
 * specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  
 * Specified type or topology must be a value in the iBase_EntityType or 
 * iMesh_EntityTopology enumerations, respectively.
 ******************************************************************************/

void iMesh_initEntArrIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being iterated */
    const int requested_entity_type, 
        /**< [in] Type of entity to iterate */
    const int requested_entity_topology, 
        /**< [in] Topology of entity to iterate */
    const int requested_array_size, 
        /**< [in] Size of chunks of handles returned for each value of the
             iterator */
    const int resilient,
        /**< [in] If zero, return a non-resilient iterator.
                  Otherwise, a resilient iterator (\ref resilient) */
    iBase_EntityArrIterator* entArr_iterator, 
        /**< [out] Pointer to iterator returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Get entities contained in array iterator and increment iterator
 *
 * Get the entities corresponding to an array iterator (e.g. dereference
 * the array iterator), and increment the iterator. The dereferenced
 * value(s) are returned in entity_handles. If the iterator is at the
 * end of the iteration, the dereferenced value(s) are undefined and
 * has_data will be returned with a value of zero. Otherwise, has_data
 * will be returned with a non-zero value.
 ******************************************************************************/

void iMesh_getNextEntArrIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityArrIterator entArr_iterator, 
        /**< [in] Iterator being queried */
    iBase_EntityHandle** entity_handles, 
        /**< [in,out] Pointer to array of entity handles contained in current
           value of iterator \ref trio) */
    int* entity_handles_allocated, 
        /**< [in,out] Pointer to allocated size of entity_handles */
    int* entity_handles_size, 
        /**< [out] Pointer to occupied size of entity_handles  */
    int* has_data, 
        /**< [out] Pointer to a flag indicating if the value(s) returned
           in entity_handles are valid. A non-zero value indicates the value(s)
           are valid. A zero value indicates the value(s) are NOT valid. */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Reset the array iterator
 *
 ******************************************************************************/

void iMesh_resetEntArrIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityArrIterator entArr_iterator, 
        /**< [in] Iterator to reset */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Destroy the specified array iterator
 *
 ******************************************************************************/

void iMesh_endEntArrIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityArrIterator entArr_iterator, 
        /**< [in] Iterator which gets destroyed */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Get the entity topology for the specified entities
 *
 * Get the entity topology for the specified entities.  Topologies 
 * returned are values in the iMesh_EntityTopology enumeration.
 ******************************************************************************/

void iMesh_getEntArrTopo(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Array of entity handles being queried */
    const int entity_handles_size, 
        /**< [in] Number of entities in entity_handles array */
    int** topology, 
        /**< [in,out] Pointer to array of entity topologies returned
             \ref trio) */
    int* topology_allocated, 
        /**< [in,out] Pointer to allocated size of topology array */
    int* topology_size, 
        /**< [out] Pointer to occupied size of topology array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Get the entity type for the specified entities
 *
 * Get the entity type for the specified entities.  Types
 * returned are values in the iBase_EntityType enumeration.
 ******************************************************************************/

void iMesh_getEntArrType(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Array of entity handles being queried */
    const int entity_handles_size, 
        /**< [in] Number of entities in entity_handles array */
    int** type, 
        /**< [in,out] Pointer to array of types returned from function
             \ref trio) */
    int* type_allocated, 
        /**< [in,out] Pointer to allocated size of type array */
    int* type_size, 
        /**< [out] Pointer to occupied size of type array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Adjacencies
 * \brief  Get entities of specified type adjacent to entities
 *
 * Get entities of specified type adjacent to entities.  Specified type
 * must be value in the iBase_EntityType enumeration.  \em offset(i) is
 * index of first entity in adjacentEntityHandles array adjacent to 
 * entity_handles[i].  More precisely, the entities adjacent to the
 * ith entity in entity_handles are found in adjacentEntityHandles
 * running from offset[i] to offset[i+1] - 1.  This implies that the
 * offset_size will be entity_handles_size + 1.
 *
 * Note 1: Because 'adjacent' as defined by the iMesh data model refers
 *         to those entities that bound another, the entities being queried
 *         here (in entity_handles arg) are NEVER ALSO returned in
 *         adjacentEntityHandles even if the entity_type_requested
 *         matches the entity type(s) in entity_handles.
 ******************************************************************************/

void iMesh_getEntArrAdj(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Array of entity handles being queried */
    const int entity_handles_size, 
        /**< [in] Number of entities in entity_handles array */
    const int entity_type_requested, 
        /**< [in] Type of adjacent entities requested */
    iBase_EntityHandle** adjacentEntityHandles, 
        /**< [in,out] Pointer to array of adjacentEntityHandles \ref trio)
           returned from function. Note that the implicit INTERLEAVED storage
           order rule applies (see iBase_StorageOrder) */
    int* adjacentEntityHandles_allocated, 
        /**< [in,out] Pointer to allocated size of adjacentEntityHandles array */
    int* adj_entity_handles_size, 
        /**< [out] Pointer to occupied size of adjacentEntityHandles array */
    int** offset, 
        /**< [in,out] Pointer to array of offsets returned from function
             \ref trio) */
    int* offset_allocated, 
        /**< [in,out] Pointer to allocated size of offset array */
    int* offset_size, 
        /**< [out] Pointer to occupied size of offset array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Adjacencies
 * \brief  Get "2nd order" adjacencies to an array of entities
 *
 * Get "2nd order" adjacencies to an array of entities, that is, from each 
 * entity, through other entities of a specified "bridge" dimension, to other 
 * entities of another specified "to" dimension.
 * Note 1:  If the "bridge" dimension is the same as the "to" dimension,
 *    the output will be empty (and an error code of
 *    iBase_INVALID_ARGUMENT returned).  If the type of a particular
 *    entity matches the "bridge" dimension, there will be no entities
 *    returned for that input entity.  This is consistent with the
 *    definition of adjacencies and the behavior of iMesh first
 *    adjacency calls. 
 * Note 2:  An entity will never be returned as a second adjacency of
 *    itself, on the grounds that this is the most likely expectation of
 *    applications, and that it is easier for an application to add the
 *    original entity to the returned data than to find and remove it.
 * Note 3: The entities adjacent to the ith entity in entity_handles are
 *    found in adj_entity_handles running from offset[i] to
 *    offset[i+1] - 1.  This implies that the offset_size will be
 *    entity_handles_size + 1. 
 ******************************************************************************/

void iMesh_getEntArr2ndAdj(
     iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle const* entity_handles, 
        /**< [in] Entities from which adjacencies are requested */
    int entity_handles_size, 
        /**< [in] Number of entities whose adjacencies are requested */
    int bridge_entity_type, 
        /**< [in]  Type of bridge entity for 2nd order adjacencies */
    int requested_entity_type, 
        /**< [in] Type of adjacent entities returned */
    iBase_EntityHandle** adj_entity_handles, 
        /**< [in,out] Adjacent entities. Note that the implicit INTERLEAVED
             storage order rule applies (see iBase_StorageOrder)
             \ref trio) */
    int* adj_entity_handles_allocated, 
        /**< [in,out] Allocated size of returned array */
    int* adj_entity_handles_size, 
        /**< [out] Occupied size of returned array */
    int** offset, 
        /**< [in,out] Offset[i] is offset into adj_entity_handles of 2nd order
           adjacencies of ith entity in entity_handles \ref trio) */
    int* offset_allocated, 
        /**< [in,out] Allocated size of offset array */
    int* offset_size, 
        /**< [out] Occupied size of offset array */
    int* err   
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Adjacencies
 * \brief  Get indexed representation of mesh or subset of mesh
 *
 * Given an entity set and optionally a type or topology, return:
 * - The entities in the set of the specified type or topology
 * - The entities adjacent to those entities with a specified
 *    type, as a list of unique handles.
 * - For each entity in the first list, the adjacent entities,
 *    specified as indices into the second list.
 *
 * Note 1: Because 'adjacent' as defined by the iMesh data model refers
 *         to those entities that bound another, the entities being queried
 *         here (in entity_set_handle arg) are NEVER ALSO returned in
 *         adj_entity_handles even if the entity_type_requested
 *         matches the entity type(s) in entity_set_handle.
 * Note 2: The entities adjacent to the ith entity in entity_handles are
 *         found in adj_entity_handles running from offset[i] to
 *         offset[i+1] - 1.  This implies that the offset_size will
 *         be entity_handles_size + 1.
 * Note 3: This function will fail and return an error if the caller passes
 *         a combination of entity_type and entity_topology that are
 *         not consistent.
 ******************************************************************************/

void iMesh_getAdjEntIndices(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set_handle, 
        /**< [in] The set of entities from which to query */
    int entity_type_requester, 
        /**< [in] If not iBase_ALL_TYPES, act only on the subset of entities
           contained in entity_set_handle of the specified type. */
    int entity_topology_requester, 
        /**< [in] If not iMesh_ALL_TOPOLOGIES, act only on the subset of
           entities contained in entity_set_handle of specified topology. */
    int entity_type_requested, 
        /**< [in] The type of the adjacent entities to return */
    iBase_EntityHandle** entity_handles, 
        /**< [in,out] The handles of the (non-struct) subset of the entities
             contained in entity_set_handle indicated by the optional type and
             topology filtering arguments. \ref trio) */
    int* entity_handles_allocated, 
        /**< [in,out] Allocated size of entity_handles array */
    int* entity_handles_size, 
        /**< [out] Occupied size of entity_handles array */
    iBase_EntityHandle** adj_entity_handles, 
        /**< [in,out] The union of the unique entities of type
           requested_entity_type adjacent to each entity in entity_handles.
           Note that the implicit INTERLEAVED storage order rule
           applies (see iBase_StorageOrder) \ref trio) */
    int* adj_entity_handles_allocated, 
        /**< [in,out] Allocated size of adj_entity_handles array */
    int* adj_entity_handles_size, 
        /**< [out] Occupied size of adj_entity_handles array */
    int** adj_entity_indices, 
        /**< [in,out] For each entity in entity_handles, the adjacent
           entities of type entity_type_requested, specified as indices into
           adj_entity_handles.  The values are concatenated into a single
           array in the order of the entity handles in entity_handles.
           \ref trio) */
    int* adj_entity_indices_allocated, 
        /**< [in,out] Allocated size of adj_entity_indices array */
    int* adj_entity_indices_size, 
        /**< [out] Occupied size of adj_entity_indices array */
    int** offset, 
        /**< [in,out] For each entity in the corresponding position in
           entity_handles, the position in adj_entity_indices at which
           values for that entity are stored \ref trio) */
    int* offset_allocated, 
        /**< [in,out] Allocated size of offset array */
    int* offset_size, 
        /**< [out] Occipied size of offset array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Create an entity set
 *
 * Create an entity set, either ordered (isList=1) or unordered (isList=0).
 * Unordered entity sets can contain a given entity or set only once.
 *
 * An entity set in iMesh supports managing its contained members in one of
 * two modes. When an entitiy set is first created, the caller is required
 * to indicate which mode of membership the set will support. The two modes
 * that are supported are a) order-preserving (isList=1) or
 * b) duplicate-preventing (isList=0).
 *
 * For order-preserving membership, the implementation will permit duplicate
 * entities. However, the implementation will guarantee that the order in which
 * entities are added to the set will be the same as the order in which they are
 * queried by the various methods that return the entities of a set. This order
 * preserving guarantee holds across removals. However, it does not hold across
 * removals followed by re-additions of the previously removed entities. This
 * kind of an entity set behaves like an STL vector or STL list and is created
 * by setting isList=1 when creating an entity set.
 *
 * For duplicate-preventing membership, the implementation will guarantee that
 * duplicate entities are prevented. Any attempts to add duplicate entities to
 * such a set will be detected, prevented and silently ignored by the
 * implementation.  This kind of entity set behaves like an STL set and
 * is created by setting isList=0 when creating an entity set.
 *
 * Finally, all of the above comments apply only to entity members of an entity
 * set and do not apply to the entity set members. Order-preserving and
 * duplicate preventing behavior for entity set members is unspecified. Each
 * implementation may behave differently for entity set members.
 * This design was chosen because we could not identify any use cases where
 * order-preserving behavior for set members was necessary. However, if users
 * encounter situations where such behavior is desirable or necessary, then
 * the ITAPS development team can certainly consider adjusting the interface
 * specification to support it.
 ******************************************************************************/

void iMesh_createEntSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const int isList, 
        /**< [in] If non-zero, an ordered list is created, otherwise an
           unordered set is created. */
    iBase_EntitySetHandle* entity_set_created, 
        /**< [out] Entity set created by function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Destroy an entity set
 *
 ******************************************************************************/

void iMesh_destroyEntSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set to be destroyed */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Return whether a specified set is ordered or unordered
 *
 * Return whether a specified set is ordered (*is_list=1) or 
 * unordered (*is_list=0)
 ******************************************************************************/

void iMesh_isList(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set being queried */
    int* is_list, 
        /**< [out] Pointer to flag returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Get the number of entity sets contained in a set or interface
 *
 * Get the number of entity sets contained in a set or interface.
 ******************************************************************************/

void iMesh_getNumEntSets(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being queried */
    const int num_hops, 
        /**< [in] Maximum hops from entity_set_handle to contained set, not
             inclusive of the contained set.  \ref numhops) */
    int* num_sets, 
        /**< [out] Pointer to the number of sets returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Get the entity sets contained in a set or interface
 *
 * Get the entity sets contained in a set or interface.
 ******************************************************************************/

void iMesh_getEntSets(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being queried */
    const int num_hops, 
        /**< [in] Maximum hops from entity_set_handle to contained set, not
           inclusive of the contained set \ref numhops) */
    iBase_EntitySetHandle** contained_set_handles, 
        /**< [in,out] Pointer to array of set handles returned \ref trio) */
    int* contained_set_handles_allocated, 
        /**< [in,out] Pointer to allocated length of */
    int* contained_set_handles_size, 
        /**< [out] Pointer to occupied length of */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Add an entity to a set
 *
 ******************************************************************************/

void iMesh_addEntToSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] The entity being added */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Pointer to the set being added to */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Remove an entity from a set
 *
 ******************************************************************************/

void iMesh_rmvEntFromSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] The entity being removed */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Pointer to the set being removed from */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Add an array of entities to a set
 *
 ******************************************************************************/

void iMesh_addEntArrToSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Array of entities being added */
    int entity_handles_size, 
        /**< [in] Number of entities in entity_handles array */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Pointer to the set being added to */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Remove an array of entities from a set
 *
 ******************************************************************************/

void iMesh_rmvEntArrFromSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Array of entities being remove */
    int entity_handles_size, 
        /**< [in] Number of entities in entity_handles array */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Pointer to the set being removed from */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Add an entity set to a set
 *
 * Add an entity set to a set (\ref cycles)
 ******************************************************************************/

void iMesh_addEntSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set_to_add, 
        /**< [in] The entity set being added */
    iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Pointer to the set being added to */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Remove an entity set from a set
 *
 * Remove an entity set from a set
 ******************************************************************************/

void iMesh_rmvEntSet(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set_to_remove, 
        /**< [in] The entity set being removed */
    iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Pointer to the set being removed from */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Return whether an entity is contained in another set
 *
 * Return whether an entity is contained (*is_contained=1) or not 
 * contained (*is_contained=0) in another set
 ******************************************************************************/

void iMesh_isEntContained(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle containing_entity_set, 
        /**< [in] Entity set being queried */
    iBase_EntityHandle contained_entity, 
        /**< [in] Entity potentially contained in containing_entity_set */
    int* is_contained, 
        /**< [out] Pointer to flag returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Return whether entities are contained in a set
 *
 ******************************************************************************/

void iMesh_isEntArrContained(
     iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle containing_entity_set, 
        /**< [in] Entity set being queried */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] List of entities for which to check containment. */
    int num_entity_handles, 
        /**< [in] Size of entity_handles array of entities to be checked. */
    int** is_contained, 
        /**< [in,out] One value for each input entity, 1 if contained in set,
             zero otherwise.  \ref trio) */
    int* is_contained_allocated, 
        /**< [in,out] Allocated size of is_contained array */
    int* is_contained_size, 
        /**< [out] Occupied size of is_contained array */
    int* err   
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySets
 * \brief  Return whether an entity set is contained in another set
 *
 * Return whether a set is contained (*is_contained=1) or not contained
 * (*is_contained=0) in another set
 ******************************************************************************/

void iMesh_isEntSetContained(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle containing_entity_set, 
        /**< [in] Entity set being queried */
    const iBase_EntitySetHandle contained_entity_set, 
        /**< [in] Entity set potentially contained in containing_entity_set */
    int* is_contained, 
        /**< [out] Pointer to flag returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  ParentChildLinks
 * \brief  Add parent/child links between two sets
 *
 * Add parent/child links between two sets.  Makes parent point to child
 * and child point to parent. (\ref cycles)
 ******************************************************************************/

void iMesh_addPrntChld(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle parent_entity_set, 
        /**< [in] Pointer to parent set */
    iBase_EntitySetHandle child_entity_set, 
        /**< [in] Pointer to child set */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  ParentChildLinks 
 * \brief  Remove parent/child links between two sets
 *
 * Remove parent/child links between two sets.
 ******************************************************************************/

void iMesh_rmvPrntChld(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle parent_entity_set, 
        /**< [in] Pointer to parent set */
    iBase_EntitySetHandle child_entity_set, 
        /**< [in] Pointer to child set */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  ParentChildLinks
 * \brief  Return whether two sets are related by parent/child links
 *
 * Return whether two sets are related (*is_child=1) or not (*is_child=0)
 * by parent/child links
 ******************************************************************************/

void iMesh_isChildOf(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle parent_entity_set, 
        /**< [in] Pointer to parent set */
    const iBase_EntitySetHandle child_entity_set, 
        /**< [in] Pointer to child set */
    int* is_child, 
        /**< [out] Pointer to flag returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  ParentChildLinks
 * \brief  Get the number of child sets linked from a specified set
 *
 * Get the number of child sets linked from a specified set.
 ******************************************************************************/

void iMesh_getNumChld(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set being queried */
    const int num_hops, 
        /**< [in] Maximum hops from entity_set_handle to child set, not
           inclusive of the child set.  \ref numhops) */
    int* num_child, 
        /**< [out] Pointer to number of children returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  ParentChildLinks
 * \brief  Get the number of parent sets linked from a specified set
 *
 * Get the number of parent sets linked from a specified set.
 ******************************************************************************/

void iMesh_getNumPrnt(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set being queried */
    const int num_hops, 
        /**< [in] Maximum hops from entity_set_handle to parent set, not
           inclusive of the parent set.  \ref numhops) */
    int* num_parent, 
        /**< [out] Pointer to number of parents returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  ParentChildLinks
 * \brief  Get the child sets linked from a specified set
 *
 * Get the child sets linked from a specified set.
 ******************************************************************************/

void iMesh_getChldn(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle from_entity_set, 
        /**< [in] Entity set being queried */
    const int num_hops, 
        /**< [in] Maximum hops from entity_set_handle to child set,
           \ref numhops) */
    iBase_EntitySetHandle** entity_set_handles, 
        /**< [in,out] Pointer to array of child sets \ref trio) */
    int* entity_set_handles_allocated, 
        /**< [in,out] Pointer to allocated size of  */
    int* entity_set_handles_size, 
        /**< [out] Pointer to occupied size of  */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  ParentChildLinks
 * \brief  Get the parent sets linked from a specified set
 *
 * Get the parent sets linked from a specified set.
 ******************************************************************************/

void iMesh_getPrnts(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle from_entity_set, 
        /**< [in] Entity set being queried */
    const int num_hops, 
        /**< [in] Maximum hops from entity_set_handle to parent set,
           \ref numhops) */
    iBase_EntitySetHandle** entity_set_handles, 
        /**< [in,out] Pointer to array of parent sets \ref trio) */
    int* entity_set_handles_allocated, 
        /**< [in,out] Pointer to allocated size of  */
    int* entity_set_handles_size, 
        /**< [out] Pointer to occupied size of  */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  VertexEntities
 * \brief  Set coordinates for an array of vertices
 *
 * Set coordinates for an array of vertices.  Specified storage 
 * order must be either iBase_INTERLEAVED or iBase_BLOCKED, and 
 * indicates order of x, y, and z coordinates in coordinate array.
 ******************************************************************************/

void iMesh_setVtxArrCoords(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* vertex_handles, 
        /**< [in] Array of vertex handles */
    const int vertex_handles_size, 
        /**< [in] Number of vertex handles in array */
    const int storage_order, 
        /**< [in] Storage order of coordinates in coordinate array */
    const double* new_coords, 
        /**< [in] Coordinate array */
    const int new_coords_size, 
        /**< [in] Size of coordinate array; should be  */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  VertexEntities
 * \brief  Create an array of new vertices at specified coordinates
 *
 * Create an array of new vertices at specified coordinates.  Value of
 * storage_order must be either iBase_INTERLEAVED or iBase_BLOCKED.
 ******************************************************************************/

void iMesh_createVtxArr(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const int num_verts, 
        /**< [in] Number of new vertices to be created */
    const int storage_order, 
        /**< [in] Storage order of coordinates in new_coords array
             (see iBase_StorageOrder) */
    const double* new_coords, 
        /**< [in] Array of coordinates of new vertices. Should be G*num_verts
           in length where G is geometric dimension of the mesh. */
    const int new_coords_size, 
        /**< [in] Number of coordinates in new_coords array, should */
    iBase_EntityHandle** new_vertex_handles, 
        /**< [in,out] Pointer to array of new vertex handles \ref trio) */
    int* new_vertex_handles_allocated, 
        /**< [in,out] Pointer to allocated size of  */
    int* new_vertex_handles_size, 
        /**< [out] Pointer to occupied size of  */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Create an array of new entities with specified lower-order topology
 *
 * Create an array of new entities with specified lower-order topology.  
 * Specified new_entity_topology must be value in iMesh_EntityTopology
 * enumeration.  Values return in status array must be values in the
 * iBase_CreationStatus enumeration.
 ******************************************************************************/

void iMesh_createEntArr(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const int new_entity_topology, 
        /**< [in] Topology of created entity */
    const iBase_EntityHandle* lower_order_entity_handles, 
        /**< [in] Array of lower order entity handles to be used to construct
             new entities */
    const int lower_order_entity_handles_size, 
        /**< [in] Number of entities in lower_order_entity_handles array */
    iBase_EntityHandle** new_entity_handles, 
        /**< [in,out] Pointer to array of new_entity_handles
             \ref trio) */
    int* new_entity_handles_allocated, 
        /**< [in,out] Pointer to allocated size of  */
    int* new_entity_handles_size, 
        /**< [out] Pointer to occupied size of  */
    int** status, 
        /**< [in,out] Pointer to array of creation status returned from
           \ref trio) (see iBase_CreationStatus) */
    int* status_allocated, 
        /**< [in,out] Pointer to allocated size of status array */
    int* status_size, 
        /**< [out] Pointer to occupied size of status array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Delete specified entities
 *
 * Delete specified entities
 ******************************************************************************/

void iMesh_deleteEntArr(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Array of entity handles to be deleted */
    const int entity_handles_size, 
        /**< [in] Number of entities in array to be deleted */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Create a tag with specified name, size, and type
 *
 * Create a tag with specified name, size, and type.  Tag size is in
 * units of size of tag_type data types.  Value input for tag type must be 
 * value in iBase_TagType enumeration.
 ******************************************************************************/

void iMesh_createTag(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const char* tag_name, 
        /**< [in] Character string indicating tag name */
    const int tag_size, 
        /**< [in] Size of each tag value, in units of number of
             tag_type datums. */
    const int tag_type, 
        /**< [in] Data type for data stored in this tag */
    iBase_TagHandle* tag_handle, 
        /**< [out] Pointer to tag handle returned from function */
    int* err, 
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int tag_name_len  
        /**< [in] Length of tag name string (\ref strlen) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Destroy a tag
 *
 * Destroy a tag.  If forced is non-zero and entities still have values
 * set for this tag, tag is deleted anyway and those values disappear,
 * otherwise tag is not deleted.
 ******************************************************************************/

void iMesh_destroyTag(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_TagHandle tag_handle, 
        /**< [in] Handle of tag to be deleted */
    const int forced, 
        /**< [in] If non-zero, delete the tag even if entities have values set
           for the tag. */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Get the name for a given tag handle
 *
 * Get the name for a given tag handle
 ******************************************************************************/

void iMesh_getTagName(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag handle being queried */
    char* name, 
        /**< [in,out] Pointer to character string to store name returned from */
    int* err, 
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int name_len  
        /**< [in] Length of character string input to function
             (\ref strlen) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Get size of a tag in units of numbers of tag data type
 *
 * Get size of a tag in units of numbers of tag data type
 ******************************************************************************/

void iMesh_getTagSizeValues(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_TagHandle tag_handle, 
        /**< [in] Handle of tag being queried */
    int* tag_size, 
        /**< [out] Pointer to tag size returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Get size of a tag in units of bytes
 *
 * Get size of a tag in units of bytes
 ******************************************************************************/

void iMesh_getTagSizeBytes(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_TagHandle tag_handle, 
        /**< [in] Handle of tag being queried */
    int* tag_size, 
        /**< [out] Pointer to tag size returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Get a the handle of an existing tag with the specified name
 *
 * Get a the handle of an existing tag with the specified name
 ******************************************************************************/

void iMesh_getTagHandle(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const char* tag_name, 
        /**< [in] Name of tag being queried */
    iBase_TagHandle* tag_handle, 
        /**< [out] Pointer to tag handle returned from function */
    int* err, 
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int tag_name_len  
        /**< [in] Length of tag name string (\ref strlen) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Get the data type of the specified tag handle
 *
 * Get the data type of the specified tag handle.  Tag type is a value in
 * the iBase_TagType enumeration.
 ******************************************************************************/

void iMesh_getTagType(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_TagHandle tag_handle, 
        /**< [in] Handle for the tag being queried */
    int* tag_type, 
        /**< [out] Pointer to tag type returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Set a tag value of arbitrary type on an entity set
 *
 * Set a tag value of arbitrary type on an entity set. The tag data
 * is passed as void*. tag_value_size specifies the size of the memory
 * pointed to by tag_value in terms of bytes. Applications are free to
 * use this function to set data of any type, not just iBase_BYTES.
 * However, in all cases, the size specified by tag_value_size is
 * always in terms of bytes.
 ******************************************************************************/

void iMesh_setEntSetData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    const void* tag_value, 
        /**< [in] Pointer to tag data being set on entity set */
    const int tag_value_size, 
        /**< [in] Size in bytes of tag data */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Set a tag value of integer type on an entity set
 *
 * Set a tag value of integer type on an entity set.
 ******************************************************************************/

void iMesh_setEntSetIntData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    const int tag_value, 
        /**< [in] Tag value being set on entity set */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Set a tag value of double type on an entity set
 *
 * Set a tag value of double type on an entity set.
 ******************************************************************************/

void iMesh_setEntSetDblData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    const double tag_value, 
        /**< [in] Tag value being set on entity set */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Set a tag value of entity handle type on an entity set
 *
 * Set a tag value of entity handle type on an entity set.
 ******************************************************************************/

void iMesh_setEntSetEHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    const iBase_EntityHandle tag_value, 
        /**< [in] Tag value being set on entity set */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Set a tag value of entity set handle type on an
 *
 *         entity set
 * Set a tag value of entity set handle type on an entity set.
 ******************************************************************************/

void iMesh_setEntSetESHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    const iBase_EntitySetHandle tag_value, 
        /**< [in] Tag value being set on entity set */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Get the value of a tag of arbitrary type on an entity set
 *
 * Get the value of a tag of arbitrary type on an entity set.  Tag data 
 * is returned back as void*. tag_value_size specifies the size of the
 * memory pointed to by tag_value in terms of bytes. Applications may
 * use this function to get data of any type, not just iBase_BYTES.
 * However because this function supports data of arbitrary type,
 * in all cases the size specified by tag_value_size is always in terms
 * of bytes.
 ******************************************************************************/

void iMesh_getEntSetData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    void* tag_value, 
        /**< [in,out] Pointer to tag data array being queried
             \ref trio) */
    int* tag_value_allocated, 
        /**< [in,out] Pointer allocated size, in bytes, of tag_value array. */
    int* tag_value_size, 
        /**< [out] Pointer to occupied size, in bytes, of tag_value array. */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Get the value of a tag of integer type on an entity set
 *
 * Get the value of a tag of integer type on an entity set.
 ******************************************************************************/

void iMesh_getEntSetIntData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    int* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Get the value of a tag of double type on an entity set
 *
 * Get the value of a tag of double type on an entity set.
 ******************************************************************************/

void iMesh_getEntSetDblData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    double* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Get the value of a tag of entity handle type on an entity set
 *
 * Get the value of a tag of entity handle type on an entity set.
 ******************************************************************************/

void iMesh_getEntSetEHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    iBase_EntityHandle* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnSets
 * \brief  Get the value of a tag of entity set handle type on an
 *
 *         entity set
 * Get the value of a tag of entity set handle type on an entity set.
 ******************************************************************************/

void iMesh_getEntSetESHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set, 
        /**< [in] Entity set on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity set */
    iBase_EntitySetHandle* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags 
 * \brief  Get all the tags associated with a specified entity set
 *
 * Get all the tags associated with a specified entity set
 ******************************************************************************/

void iMesh_getAllEntSetTags(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity being queried */
    iBase_TagHandle** tag_handles, 
        /**< [in,out] Pointer to array of tag_handles returned from
             \ref trio) */
    int* tag_handles_allocated, 
        /**< [in,out] Pointer to allocated size of tag_handles  */
    int* tag_handles_size, 
        /**< [out] Pointer to occupied size of tag_handles array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Remove a tag value from an entity set
 *
 * Remove a tag value from an entity set
 ******************************************************************************/

void iMesh_rmvEntSetTag(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set from which tag is being removed */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag handle of tag being removed */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  VertexEntities
 * \brief  Set coordinates for a vertex
 *
 * Set coordinates for a vertex.
 ******************************************************************************/

void iMesh_setVtxCoord(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle vertex_handle, 
        /**< [in] vertex handle being set */
    const double x, 
        /**< [in] x coordinate being set */
    const double y, 
        /**< [in] y coordinate being set */
    const double z, 
        /**< [in] z coordinate being set */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  VertexEntities
 * \brief  Create a new vertex at specified coordinates
 *
 * Create a new vertex at specified coordinates.
 ******************************************************************************/

void iMesh_createVtx(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const double x, 
        /**< [in] x coordinate of new vertex */
    const double y, 
        /**< [in] y coordinate of new vertex */
    const double z, 
        /**< [in] z coordinate of new vertex */
    iBase_EntityHandle* new_vertex_handle, 
        /**< [out] Pointer to new vertex handles returned from  */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities 
 * \brief  Create a new entity with specified lower-order topology
 *
 * Create a new entity with specified lower-order topology.  
 * Specified new_entity_topology must be value in iMesh_EntityTopology
 * enumeration.  Value returned as status must be a value in the
 * iBase_CreationStatus enumeration.
 ******************************************************************************/

void iMesh_createEnt(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const int new_entity_topology, 
        /**< [in] Topology of created entity */
    const iBase_EntityHandle* lower_order_entity_handles, 
        /**< [in] Array of lower order entity handles to be used to construct
           new entity. */
    const int lower_order_entity_handles_size, 
        /**< [in] Number of entities in lower_order_entity_handles array */
    iBase_EntityHandle* new_entity_handle,
        /**< [out] Pointer to new entity handle returned from  */
    int* status, 
        /**< [out] Pointer to creation status returned from function
           (see iBase_CreationStatus) */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Delete specified entity
 *
 * Delete specified entity
 ******************************************************************************/

void iMesh_deleteEnt(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity to be deleted */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Get tag values of arbitrary type for an array of entities
 *
 * Get tag values of arbitrary type for an array of entities.  Tag data 
 * is returned as void*. tag_values_size specifies the size of the
 * memory pointed to by tag_values in terms of bytes. Applications may
 * use this function to get data of any type, not just iBase_BYTES.
 * However, because this function supports data of arbitrary type, in
 * all cases the size specified by tag_values_size always in terms of
 * bytes.
 ******************************************************************************/

void iMesh_getArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    void* tag_values, 
        /**< [in,out] Pointer to tag data array being returned. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) \ref trio) */
    int* tag_values_allocated, 
        /**< [in,out] Pointer to allocated size of tag data array */
    int* tag_values_size, 
        /**< [out] Pointer to occupied size in bytes of tag data */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Get tag values of integer type for an array of entities
 *
 * Get tag values of integer type for an array of entities.
 ******************************************************************************/

void iMesh_getIntArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    int** tag_values, 
        /**< [in,out] Pointer to tag data array being returned. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) \ref trio) */
    int* tag_values_allocated, 
        /**< [in,out] Pointer to allocated size of tag data array */
    int* tag_values_size, 
        /**< [out] Pointer to occupied size of tag data array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Get tag values of double type for an array of entities
 *
 * Get tag values of double type for an array of entities.
 ******************************************************************************/

void iMesh_getDblArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    double** tag_values, 
        /**< [in,out] Pointer to tag data array being returned. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) \ref trio) */
    int* tag_values_allocated, 
        /**< [in,out] Pointer to allocated size of tag data array */
    int* tag_values_size, 
        /**< [out] Pointer to occupied size of tag data array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Get tag values of entity handle type for an array of entities
 *
 * Get tag values of entity handle type for an array of entities.
 ******************************************************************************/

void iMesh_getEHArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    iBase_EntityHandle** tag_value, 
        /**< [in,out] Pointer to tag data array being returned. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) \ref trio) */
    int* tag_value_allocated, 
        /**< [in,out] Pointer to allocated size of tag data array */
    int* tag_value_size, 
        /**< [out] Pointer to occupied size of tag data array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Get tag values of entity set handle type for an array of
 *
 *         entities
 * Get tag values of entity set handle type for an array of entities.
 ******************************************************************************/

void iMesh_getESHArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    iBase_EntitySetHandle** tag_value, 
        /**< [in,out] Pointer to tag data array being returned. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) \ref trio) */
    int* tag_value_allocated, 
        /**< [in,out] Pointer to allocated size of tag data array */
    int* tag_value_size, 
        /**< [out] Pointer to occupied size of tag data array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Set tag values of arbitrary type on an array of entities
 *
 * Set tag values of arbitrary type on an array of entities.  Tag data
 * is passed as void*. tag_values_size specifies the size of the
 * memory pointed to by tag_values in terms of bytes. Applications may
 * use this function to set data of any type, not just iBase_BYTES.
 * However, because this function supports data of arbitrary type, in all
 * cases the size specified by tag_values_size is always in terms of
 * bytes.
 ******************************************************************************/

void iMesh_setArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const void* tag_values, 
        /**< [in] Pointer to tag data being set on entity. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) */
    const int tag_values_size, 
        /**< [in] Size in bytes of tag data */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Set tag values of integer type on an array of entities
 *
 * Set tag values of integer type on an array of entities.
 ******************************************************************************/

void iMesh_setIntArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const int* tag_values, 
        /**< [in] Pointer to tag data being set on entities. Note that the
             implicit INTERLEAVED storage order rule applies
             (see iBase_StorageOrder) */
    const int tag_values_size, 
        /**< [in] Size in total number of integers of tag data */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Set tag values of double type on an array of entities
 *
 * Set tag values of double type on an array of entities.
 ******************************************************************************/

void iMesh_setDblArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const double* tag_values, 
        /**< [in] Pointer to tag data being set on entities. Note that the
             implicit INTERLEAVED storage order rule applies
             (see iBase_StorageOrder) */
    const int tag_values_size, 
        /**< [in] Size in total number of doubles of tag data */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Set tag values of entity handle type on an array of entities
 *
 * Set tag values of entity handle type on an array of entities.
 ******************************************************************************/

void iMesh_setEHArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const iBase_EntityHandle* tag_values, 
        /**< [in] Pointer to tag data being set on entities. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) */
    const int tag_values_size, 
        /**< [in] Size in total number of entity handles of tag  */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnArr
 * \brief  Set tag values of entity set handle type on an array of
 *
 *         entities
 * Set tag values of entity set handle type on an array of entities.
 ******************************************************************************/

void iMesh_setESHArrData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity array on which tag is being set */
    const int entity_handles_size, 
        /**< [in] Number of entities in array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const iBase_EntitySetHandle* tag_values, 
        /**< [in] Pointer to tag data being set on entities. Note that the
           implicit INTERLEAVED storage order rule applies
           (see iBase_StorageOrder) */
    const int tag_values_size, 
        /**< [in] Size in total number of entity handles of tag  */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Remove a tag value from an array of entities
 *
 * Remove a tag value from an array of entities
 ******************************************************************************/

void iMesh_rmvArrTag(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle* entity_handles, 
        /**< [in] Entity from which tag is being removed */
    const int entity_handles_size, 
        /**< [in] Number of entities in entity array */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag handle of tag being removed */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Get the value of a tag of arbitrary type on an entity
 *
 * Get the value of a tag of arbitrary type on an entity.  Tag data 
 * is passed back as void*. tag_value_size specifies the size of the
 * memory pointed to by tag_value in terms of bytes. Applications may
 * use this function to get data of any type, not just iBase_BYTES.
 * However, because this function supports arbitrary type, in all
 * cases the size specified by tag_value_size is always in terms of
 * bytes.
 ******************************************************************************/

void iMesh_getData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    void* tag_value, 
        /**< [in,out] Pointer to tag data array being queried
             \ref trio) */
    int* tag_value_allocated, 
        /**< [in,out] Pointer to tag data array allocated size */
    int* tag_value_size, 
        /**< [out] Pointer to occupied size in bytes of tag data */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Get the value of a tag of integer type on an entity
 *
 * Get the value of a tag of integer type on an entity.
 ******************************************************************************/

void iMesh_getIntData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    int* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Get the value of a tag of double type on an entity
 *
 * Get the value of a tag of double type on an entity.
 ******************************************************************************/

void iMesh_getDblData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    double* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Get the value of a tag of entity handle type on an entity
 *
 * Get the value of a tag of entity handle type on an entity.
 ******************************************************************************/

void iMesh_getEHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    iBase_EntityHandle* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Get the value of a tag of entity set handle type on an
 *
 *         entity
 * Get the value of a tag of entity set handle type on an entity.
 ******************************************************************************/

void iMesh_getESHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    iBase_EntitySetHandle* out_data, 
        /**< [out] Pointer to tag value returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Set a tag value of arbitrary type on an entity
 *
 * Set a tag value of arbitrary type on an entity.  Tag data
 * is passed as void*. tag_value_size specifies the size of the
 * memory pointed to by tag_value in terms of bytes. Applications may
 * use this function to set data of any type, not just iBase_BYTES.
 * However, because this function supports data of arbitrary type, in
 * all cases the size specified by tag_value_size is always in terms
 * of bytes.
 ******************************************************************************/

void iMesh_setData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const void* tag_value, 
        /**< [in] Pointer to tag data being set on entity */
    const int tag_value_size, 
        /**< [in] Size in bytes of tag data */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Set a tag value of integer type on an entity
 *
 * Set a tag value of integer type on an entity.
 ******************************************************************************/

void iMesh_setIntData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const int tag_value, 
        /**< [in] Tag value being set on entity */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Set a tag value of double type on an entity
 *
 * Set a tag value of double type on an entity.
 ******************************************************************************/

void iMesh_setDblData(
    iMesh_Instance instance,
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const double tag_value, 
        /**< [in] Tag value being set on entity */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Set a tag value of entity handle type on an entity
 *
 * Set a tag value of entity handle type on an entity.
 ******************************************************************************/

void iMesh_setEHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const iBase_EntityHandle tag_value, 
        /**< [in] Tag value being set on entity */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  TagsOnEnts
 * \brief  Set a tag value of entity set handle type on an entity
 *
 * Set a tag value of entity set handle type on an entity.
 ******************************************************************************/

void iMesh_setESHData(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity on which tag is being set */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag being set on an entity */
    const iBase_EntitySetHandle tag_value, 
        /**< [in] Tag value being set on entity */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Get all the tags associated with a specified entity handle
 *
 * Get all the tags associated with a specified entity handle
 ******************************************************************************/

void iMesh_getAllTags(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity being queried */
    iBase_TagHandle** tag_handles, 
        /**< [in,out] Pointer to array of tag_handles returned from
             \ref trio) */
    int* tag_handles_allocated, 
        /**< [in,out] Pointer to allocated size of tag_handles  */
    int* tag_handles_size, 
        /**< [out] Pointer to occupied size of tag_handles array */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Remove a tag value from an entity
 *
 * Remove a tag value from an entity
 ******************************************************************************/

void iMesh_rmvTag(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity from which tag is being removed */
    const iBase_TagHandle tag_handle, 
        /**< [in] Tag handle of tag being removed */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators 
 * \brief  Initialize an iterator over specified entity type, topology, and size
 *
 * Initialize an iterator over specified entity type, topology, and size,
 * for a specified set or instance.  Iterator returned can be used as input
 * to functions returning the entity for the iterator.  If all entities of 
 * a specified type and/or topology are to be iterated, specify 
 * iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  Specified type 
 * or topology must be a value in the iBase_EntityType or 
 * iMesh_EntityTopology enumerations, respectively.
 *
 * Note: This function will fail and return an error if the caller passes a
 * combination of entity_type and entity_topology that are not consistent.
 ******************************************************************************/

void iMesh_initEntIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_handle, 
        /**< [in] Entity set being iterated */
    const int requested_entity_type, 
        /**< [in] Type of entity to iterate */
    const int requested_entity_topology, 
        /**< [in] Topology of entity to iterate */
    const int resilient,
        /**< [in] If zero, return a non-resilient iterator.
                  Otherwise, a resilient iterator (\ref resilient) */
    iBase_EntityIterator* entity_iterator, 
        /**< [out] Pointer to iterator returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators 
 * \brief  Get entity corresponding to an iterator and increment iterator
 *
 * Get the entity corresponding to an iterator (that is, dereference
 * the iterator), and increment the iterator. The dereferenced value
 * is returned in 'entity_handle'. If the iterator is at the end of the
 * iteration, the dereferenced value will be undefined and 'has_data'
 * will have a value of zero. Otherwise, 'has_data' will have a non-zero
 * value.
 ******************************************************************************/

void iMesh_getNextEntIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityIterator entity_iterator, 
        /**< [in] Iterator being queried */
    iBase_EntityHandle* entity_handle, 
        /**< [out] Pointer to an entity handle corresponding to the current
           value of iterator just prior to the call. */
    int* has_data, 
        /**< [out] Pointer to a flag indicating if the value returned in
           entity_handle is valid. A non-zero value indicates the value is
           valid. A zero value indicates the value is NOT valid. */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Reset the iterator
 *
 * Reset the iterator
 ******************************************************************************/

void iMesh_resetEntIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityIterator entity_iterator, 
        /**< [in] Iterator to reset */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Destroy the specified iterator
 *
 * Destroy the specified iterator
 ******************************************************************************/

void iMesh_endEntIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityIterator entity_iterator, 
        /**< [in] Iterator which gets destroyed */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Get the entity topology for the specified entity
 *
 * Get the entity topology for the specified entity.  Topology
 * returned is a value in the iMesh_EntityTopology enumeration.
 ******************************************************************************/

void iMesh_getEntTopo(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity handle being queried */
    int* out_topo, 
        /**< [out] Pointer to entity topology returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Entities
 * \brief  Get the entity type for the specified entity
 *
 * Get the entity type for the specified entity.  Type returned is a value
 * in the iBase_EntityType enumeration.
 ******************************************************************************/

void iMesh_getEntType(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity handle being queried */
    int* out_type, 
        /**< [out] Pointer to entity type returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  VertexEntities
 * \brief  Get coordinates of specified vertex
 *
 * Get coordinates of specified vertex.
 ******************************************************************************/

void iMesh_getVtxCoord(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle vertex_handle, 
        /**< [in] Mesh vertex being queried */
    double* x, 
        /**< [out] Pointer to x coordinate returned from function */
    double* y, 
        /**< [out] Pointer to y coordinate returned from function */
    double* z, 
        /**< [out] Pointer to z coordinate returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Adjacencies
 * \brief  Get entities of specified type adjacent to an entity
 *
 * Get entities of specified type adjacent to an entity.  Specified type
 * must be value in the iBase_EntityType enumeration.
 *
 * Note 1: Because 'adjacent' as defined by the iMesh data model refers
 *         to those entities that bound another, the entity being queried
 *         here (in entity_handle arg) is NEVER ALSO returned in
 *         adj_entity_handles even if the entity_type_requested
 *         matches the entity type in entity_handle.
 ******************************************************************************/

void iMesh_getEntAdj(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntityHandle entity_handle, 
        /**< [in] Entity handle being queried */
    const int entity_type_requested, 
        /**< [in] Type of adjacent entities requested */
    iBase_EntityHandle** adj_entity_handles, 
        /**< [in,out] Pointer to array of adjacent entities \ref trio) */
    int* adj_entity_handles_allocated, 
        /**< [in,out] Pointer to allocated size of adj_entity_handles */
    int* adj_entity_handles_size, 
        /**< [out] Pointer to occupied size of adj_entity_handles */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Adjacencies
 * \brief  Get "2nd order" adjacencies to an entity
 *
 * Get "2nd order" adjacencies to an entity, that is, from an entity, through
 * other entities of a specified "bridge" dimension, to other entities of
 * another specified "to" dimension.
 * Note 1: If the "bridge" dimension is the same as the "to" dimension
 *    or the dimension of the input entity, the output will be empty
 *    (and an error code of iBase_INVALID_ARGUMENT returned).  This is
 *    consistent with the definition of adjacencies and the behavior of
 *    iMesh first adjacency calls.
 * Note 2: An entity will never be returned as a second adjacency of
 *    itself, on the grounds that this is the most likely expectation of
 *    applications, and that it is easier for an application to add the
 *    original entity to the returned data than to find and remove it.
 ******************************************************************************/

void iMesh_getEnt2ndAdj(
     iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityHandle entity_handle, 
        /**< [in] Entity from which adjacencies are requested */
    int bridge_entity_type, 
        /**< [in]  Type of bridge entity for 2nd order adjacencies */
    int requested_entity_type, 
        /**< [in] Type of adjacent entities returned */
    iBase_EntityHandle** adjacent_entities, 
        /**< [in,out] Adjacent entities \ref trio) */
    int* adjacent_entities_allocated, 
        /**< [in,out] Allocated size of returned array */
    int* adjacent_entities_size, 
        /**< [out] Occupied size of returned array */
    int* err   
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySetOperators
 * \brief  Subtract contents of one entity set from another
 *
 * Subtract contents of one entity set from another
 ******************************************************************************/

void iMesh_subtract(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_1, 
        /**< [in] Entity set from which other set is being subtracted */
    const iBase_EntitySetHandle entity_set_2, 
        /**< [in] Entity set being subtracted from other set */
    iBase_EntitySetHandle* result_entity_set, 
        /**< [out] Pointer to entity set returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySetOperators
 * \brief  Intersect contents of one entity set with another
 *
 * Intersect contents of one entity set with another
 ******************************************************************************/

void iMesh_intersect(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_1, 
        /**< [in] Entity set being intersected with another */
    const iBase_EntitySetHandle entity_set_2, 
        /**< [in] Entity set being intersected with another */
    iBase_EntitySetHandle* result_entity_set, 
        /**< [out] Pointer to entity set returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntitySetOperators
 * \brief  Unite contents of one entity set with another
 *
 ******************************************************************************/

void iMesh_unite(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    const iBase_EntitySetHandle entity_set_1, 
        /**< [in] Entity set being united with another */
    const iBase_EntitySetHandle entity_set_2, 
        /**< [in] Entity set being united with another */
    iBase_EntitySetHandle* result_entity_set, 
        /**< [out] Pointer to entity set returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \page imesh iMesh: ITAPS Serial Mesh Interface
 *
 * The ITAPS Mesh Interface iMesh provides a common interface for
 * accessing mesh and data associated with a mesh.  Applications written
 * to use this interface can use a variety of implementations, choosing
 * the one that best meets its needs.  They can also use tools written
 * to this interface, for example mesh smoothing, adaptive mesh refinement,
 * and parallel mesh support.
 *
 * The ITAPS interfaces use a data model composed of four basic data types:
 * /Entity/: basic topological entities in a mesh, e.g. vertices, 
 * triangles, hexahedra.
 * /Entity Set/: arbitrary grouping of other entities and sets. 
 * Entity sets also support parent/child relations with other sets which
 * are distinct from entities contained in those sets.  Parent/child links
 * can be used to embed graph relationships between sets, e.g. to 
 * represent topological relationships between the sets.
 * /Interface/: the object with which mesh is associated and on which
 * functions in iMesh are called.
 * /Tag/: application data associated with objects of any of the other 
 * data types.  Each tag has a designated name, size, and data type.
 *
 * ITAPS Entity Type, Topology
 * Each entity has a specific Entity Type and Entity Topology.  The Entity 
 * Type is one of VERTEX, EDGE, FACE, and REGION, and is synonymous with
 * the topological dimension of the entity.  The Entity Topology denotes
 * the specific shape, for example TRIANGLE, QUADRILATERAL, TETRAHEDRON,
 * and HEXAHEDRON.  Entity Type and Entity Topology exist as enumerated
 * types, Entity Type in the iBase_EntityType enumeration, and
 * Entity Topology in the iMesh_EntityTopology enumeration.
 *
 * ITAPS Entity-, Array-, and Iterator-Based Access
 * The iMesh interface provides functions for accessing entities
 * individually, as arrays of entities, or using iterators.  These access
 * methods have different memory versus execution time tradeoffs, 
 * depending on the implementation.
 *
 * \image html example_pic.jpeg
 *
 * \subpage cycles
 *
 * \page cycles Cycles in Set-Inclusion and Parent-Child structures.
 *
 * There are two graph-like structures in the iMesh interface and data
 * model; the set-inclusion structure and the parent-child link structure.
 * Whether these structures support cycles is relevant to implementors.
 *
 * Over the evolution of the iMesh data model and API, both of these
 * structures have been viewed more or less like a tree and so cycles seem
 * incompatible with that notion.
 *
 * Allowing a cycle in the set inclusion structure implies all entity sets
 * in the cycle are all equal to each other. That is the only rational,
 * model-level view that would allow them all to be (improper) subsets of
 * each other. On the other hand if the iMesh specification excludes cycles
 * from the set inclusion structure, the time complexity (performance) as a
 * function of the number of entity sets may be prohibitive for
 * implementations to detect and prevent them.
 *
 * Allowing a cycle in the parent-child link structure seems likewise hard
 * to justify. However, when the parent-child structure is viewed more like
 * a general graph (a view that the current API itself supports even if the
 * function names themselves do not suggest that) than specifically a tree,
 * the admission of cycles there is potentially more natural and useful.
 *
 * Implementations are required to support cycles in the Parent-Child
 * structure. Implementations are neither required to support nor required
 * to explicitly prevent cycles in the Set-Inclusion structure. Portable
 * applications should NOT rely on implementations support for cycles
 * in the set-inclusion structure.
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iMesh iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup Initialization Initialization
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup Entities Entities 
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup VertexEntities Vertex Entities 
 * \ingroup Entities
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup EntitySets Entity Sets
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup EntitySetOperators Entity Set Operators
 * \ingroup EntitySets
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup Adjacencies Adjacencies
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup EntityIterators Entity Iterators
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup Tags Tags
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup TagData Tag Data
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup TagsOnEnts Tag Data On Entities
 * \ingroup TagData
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup TagsOnSets Tag Data On Entity Sets
 * \ingroup TagData
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup TagsOnArr Tag Data On Arrays of Entities
 * \ingroup TagData
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup MeshModification Mesh Modification
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup ParentChildLinks Parent Child Links
 * \ingroup iMesh
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup Datatypes Datatypes
 * \ingroup iMesh
 ******************************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ifndef _ITAPS_iMesh */

