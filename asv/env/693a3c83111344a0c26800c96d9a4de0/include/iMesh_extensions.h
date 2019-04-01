#ifndef IMESH_REC_CBIND_H__
#define IMESH_REC_CBIND_H__

#include "moab/MOABConfig.h"
#include "iMesh.h"
#include "iMesh_protos.h"
#ifdef MOAB_HAVE_MPI
#include "iMeshP.h"
#include "moab_mpi.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

    /**\brief  Get entities of specific type and/or topology in set or instance, recursive
     *
     * Get entities of specific type and/or topology in set or instance.  If recursive
     * is passed in non-zero, includes entities in owned sets.  All 
     * entities of a given type or topology are requested by specifying
     * iBase_ALL_TOPOLOGIES or iBase_ALL_TYPES, respectively.  Specified type
     * or topology must be a value in the iBase_EntityType or iMesh_EntityTopology
     * enumeration, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entities being requested
     * \param entity_topology Topology of entities being requested
     * \param recursive If non-zero, gets entities in owned sets too
     * \param *entity_handles Pointer to array of entity handles returned 
     *        from function
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntitiesRec(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const int entity_type,
                            /*in*/ const int entity_topology,
                            /*in*/ const int recursive,
                            /*out*/ iBase_EntityHandle** entity_handles,
                            /*out*/ int* entity_handles_allocated,
                            /*out*/ int* entity_handles_size,
                            /*out*/ int *err);

    /**\brief  Get the number of entities with the specified type in the instance or set, recursive
     *
     * Get the number of entities with the specified type in the instance 
     * or set.  If recursive is passed in non-zero, includes entities in owned sets.  
     * If entity set handle is zero, return information for instance, 
     * otherwise for set.  Value of entity type must be from the
     * iBase_EntityType enumeration.  If iBase_ALL_TYPES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entity requested
     * \param recursive If non-zero, includes entities in owned sets too
     * \param num_type Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfTypeRec(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set_handle,
                             /*in*/ const int entity_type,
                             /*in*/ const int recursive,
                             /*out*/ int *num_type, 
                             /*out*/ int *err);

    /**\brief  Get the number of entities with the specified topology in the instance or set
     *
     * Get the number of entities with the specified topology in the instance 
     * or set.  If recursive is passed in non-zero, includes entities in owned sets.  
     * If entity set handle is zero, return information for instance,
     * otherwise for set.  Value of entity topology must be from the
     * iMesh_EntityTopology enumeration.  If iMesh_ALL_TOPOLOGIES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_topology Topology of entity requested
     * \param recursive If non-zero, includes entities in owned sets too
     * \param num_topo Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfTopoRec(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set_handle,
                             /*in*/ const int entity_topology,
                             /*in*/ const int recursive,
                             /*out*/ int *num_topo, 
                             /*out*/ int *err);


    /**\brief  Get entities with specified type, topology, tag(s) and (optionally) tag value(s)
     *
     * Get entities with the specified type, topology, tag(s), and optionally tag value(s).
     * If tag values pointer is input as zero, entities with specified tag(s) are returned,
     * regardless of their value.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entities being requested
     * \param entity_topology Topology of entities being requested
     * \param tag_handles Array of tag handles
     * \param tag_vals Array of tag values (zero if values not requested)
     * \param num_tags_vals Number of tags and optionally values
     * \param recursive If non-zero, gets entities in owned sets too
     * \param *entity_handles Pointer to array of entity handles returned 
     *        from function
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntsByTagsRec(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set_handle,
                              /*in*/ const int entity_type,
                              /*in*/ const int entity_topology,
                              /*in*/ const iBase_TagHandle *tag_handles,
                              /*in*/ const char * const *tag_vals,
                              /*in*/ const int num_tags_vals,
                              /*in*/ const int recursive,
                              /*out*/ iBase_EntityHandle** entity_handles,
                              /*out*/ int* entity_handles_allocated,
                              /*out*/ int* entity_handles_size,
                              /*out*/ int *err);

    /**\brief  Get entity sets with specified tag(s) and (optionally) tag value(s)
     *
     * Get entity sets with the specified tag(s) and optionally tag value(s).
     * If tag values pointer is input as zero, entities with specified tag(s) are returned,
     * regardless of their value.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param tag_handles Array of tag handles
     * \param tag_vals Array of tag values (zero if values not requested)
     * \param num_tags_vals Number of tags and optionally values
     * \param recursive If non-zero, gets entities in owned sets too
     * \param *set_handles Pointer to array of entity handles returned 
     *        from function
     * \param *set_handles_allocated Pointer to allocated size of 
     *        set_handles array
     * \param *set_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetsByTagsRec(iMesh_Instance instance,
                                 /*in*/ const iBase_EntitySetHandle entity_set_handle,
                                 /*in*/ const iBase_TagHandle *tag_handles,
                                 /*in*/ const char * const *tag_vals,
                                 /*in*/ const int num_tags_vals,
                                 /*in*/ const int recursive,
                                 /*out*/ iBase_EntitySetHandle** set_handles,
                                 /*out*/ int* set_handles_allocated,
                                 /*out*/ int* set_handles_size,
                                 /*out*/ int *err);

    /**\brief Get MBCN type corresponding to iMesh topology value
     *
     * Get MBCN type corresponding to iMesh topology value.  Required for input
     * to MBCN canonical numbering functions, which are written in terms of 
     * MBCN entity types.  Returns -1 for type if entity topology is out of
     * bounds, or MBMAXTYPE if no corresponding MBCN type exists.
     * \param imesh_entity_topology iMesh_EntityTopology value
     * \param mbcn_type MBEntityType corresponding to entity topology
     */
  void iMesh_MBCNType(/*in*/ const int imesh_entity_topology,
                      /*out*/ int *mbcn_type);
    

  /**\brief Access tag data via direct pointer into contiguous blocks
   *
   * Iteratively obtain direct access to contiguous blocks of tag
   * storage.  This function cannot be used with bit tags because
   * of the compressed bit storage.  This function cannot be used
   * with variable length tags because it does not provide a mechanism
   * to determine the length of the value for each entity.  This
   * function may be used with sparse tags, but if it is used, it
   * will return data for a single entity at a time.  
   *
   *\Note If this function is called for entities for which no tag value
   *      has been set, but for which a default value exists, it will 
   *      force the allocation of explicit storage for each such entity
   *      even though MOAB would normally not explicitly store tag values
   *      for such entities.
   *
   *\Example:
   *\code
   *\endcode
   */
  void iMesh_tagIterate(iMesh_Instance instance,
                          /**< [in] iMesh instance */
                        const iBase_TagHandle tag_handle,
                          /**< [in] Tag being queried */
                        iBase_EntityArrIterator entArr_iterator, 
                          /**< [in] Iterator being queried */
                        void* tag_value, 
                          /**< [out] Pointer to pointer that will be set to tag data memory */
                        int* count,
                          /**< [out] Number of contiguous entities in this subrange */
                        int* err  
                          /**< [out] Returned Error status (see iBase_ErrorType) */
                        );

  /**\brief Access connectivity data via direct pointer into contiguous blocks
   *
   * Iteratively obtain direct access to contiguous blocks of connectivity
   * storage.
   *
   */
  void iMesh_connectIterate(iMesh_Instance instance,
                              /**< [in] iMesh instance */
                            iBase_EntityArrIterator entArr_iterator, 
                              /**< [in] Iterator being queried */
                            iBase_EntityHandle **connect,
                              /**< [out] Pointer to pointer that will be set to connectivity data memory */
                            int* verts_per_entity,
                              /**< [out] Number of vertices per entity in this subrange */
                            int* count,
                              /**< [out] Number of contiguous entities in this subrange */
                            int* err  
                              /**< [out] Returned Error status (see iBase_ErrorType) */
                            );

  /**\brief Access coordinates data via direct pointer into contiguous blocks
   *
   * Iteratively obtain direct access to contiguous blocks of coordinate
   * storage.
   *
   */
  void iMesh_coordsIterate(iMesh_Instance instance,
                             /**< [in] iMesh instance */
                           iBase_EntityArrIterator entArr_iterator, 
                             /**< [in] Iterator being queried */
                           double **coordsx,
                             /**< [out] Pointer to pointer x coordinates */
                           double **coordsy,
                             /**< [out] Pointer to pointer y coordinates */
                           double **coordsz,
                             /**< [out] Pointer to pointer z coordinates */
                           int* count,
                             /**< [out] Number of contiguous entities in this subrange */
                           int* err  
                             /**< [out] Returned Error status (see iBase_ErrorType) */
                           );

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Step the iterator a specified number of entities
 *
 * Step the iterator a specified number of entities.  If this number is greater
 * than the number of entities left in the iterator, the iterator is placed
 * at the end and at_end is returned non-zero; otherwise at_end is returned zero.
 ******************************************************************************/

void iMesh_stepEntIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityIterator ent_iterator, 
        /**< [in] Iterator being queried */
    int step_length, 
        /**< [in] Number of entities to step the iterator */
    int* at_end, 
        /**< [out] Non-zero if iterator is at the end of the iteration */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

void iMesh_stepEntArrIter(
    iMesh_Instance instance, 
        /**< [in] iMesh instance handle */
    iBase_EntityArrIterator entArr_iterator, 
        /**< [in] Iterator being queried */
    int step_length, 
        /**< [in] Number of entities to step the iterator */
    int* at_end, 
        /**< [out] Non-zero if iterator is at the end of the iteration */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  EntityIterators
 * \brief  Initialize an array iterator over specified entity type, topology,
 *  and size, with an optional recursive flag.
 *
 * Initialize an array iterator over specified entity type, topology, and 
 * size, for a specified set or instance.  Iterator returned can be used 
 * as input to functions returning entities for the iterator.  If all 
 * entities of a specified type and/or topology are to be iterated, 
 * specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  
 * Specified type or topology must be a value in the iBase_EntityType or 
 * iMesh_EntityTopology enumerations, respectively.  If recursive is true,
 * entities are retrieved recursively through contained (but not child) sets.
 ******************************************************************************/

void iMesh_initEntArrIterRec(
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
    const int recursive,
      /**< [in] If non-zero, entities retrieved recursively */
    iBase_EntityArrIterator* entArr_iterator, 
        /**< [out] Pointer to iterator returned from function */
    int* err  
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup  Tags 
 * \brief  Get all the tags associated with the entire interface
 *
 * Get all the tags associated with the entire interface
 ******************************************************************************/

void iMesh_getAllIfaceTags(iMesh_Instance instance,
                           /*inout*/ iBase_TagHandle **tag_handles,
                           /*inout*/ int *tag_handles_allocated,
                           /*out*/ int *tag_handles_size,
                           /*out*/ int *err
);

/***************************************************************************//**
 * \ingroup  Tags
 * \brief  Create a tag with options
 *
 * Create a tag with options; allows creation of Dense and Bit tags through iMesh
 * Allowable options are:
 * TAG_STORAGE_TYPE={DENSE | SPARSE | BIT | MESH}
 * TAG_DEFAULT_VALUE=<value> (data type of value should match tag data type)
 ******************************************************************************/

void iMesh_createTagWithOptions(iMesh_Instance instance,
                                  /**< [in] iMesh instance handle */
                                  /*in*/ const char* tag_name,
                                  /**< [in] tag name*/
                                  /*in*/ const char* tmp_tag_options,
                                  /**< [in] options string */
                                  /*in*/ const int tag_size,
                                  /**< [in] tag size, in number of values */
                                  /*in*/ const int tag_type,
                                  /**< [in] tag data type (int, double, etc.) */
                                  /*out*/ iBase_TagHandle* tag_handle, 
                                  /**< [out] handle of new tag */
                                  /*out*/ int *err,
                                  /**< [out] error */
                                  /*in*/ const int tag_name_len,
                                  /**< [in] length of tag name string */
                                  /*in*/ const int tag_options_len);
                                  /**< [in] length of options string */
    
/***************************************************************************//**
 * \ingroup  ScdMesh
 * \brief  Create a structured mesh
 *
 * Create a structured mesh, with local and (optionally) global ijk parameters and
 * optional physical positions.  If running in parallel, can request shared vertex resolution
 * and optional number and type of ghost layers of elements.  Global parameters are used to compute
 * global ids, which are used in shared vertex resolution.
 ******************************************************************************/

void iMesh_createStructuredMesh(
        iMesh_Instance instance,
          /**< [in] iMesh instance handle */
        int *local_dims,
          /**< [in] Min/max corners of local ijk parameters, -1 for unused dimensions; specified as
                    ilo, jlo, klo, ihi, jhi, khi. */
        int *global_dims,
          /**< [in] Min/max corners of global ijk parameters, -1 for unused dimensions; NULL if running in serial. 
                    Order similar to local_dims. */
        double *i_vals,
          /**< [in] Physical positions of i values, NULL if not placed in physical space. */
        double *j_vals,
          /**< [in] Physical positions of j values, NULL if not placed in physical space. */
        double *k_vals,
          /**< [in] Physical positions of k values, NULL if not placed in physical space. */
        int resolve_shared,
          /**< [in] Non-zero if running in parallel and resolution of shared vertices is desired, zero otherwise. */
        int ghost_dim,
          /**< [in] Dimension of entities to ghost, -1 if none desired. */
        int bridge_dim,
          /**< [in] Dimension of bridge entities used to compute ghosts, -1 if no ghosts desired. */
        int num_layers,
          /**< [in] Number of layers of ghosts desired, -1 if no ghosts desired. */
        int addl_ents,
          /**< [in] Dimension of addition entities adjacent to ghosts to exchange. */
        int vert_gids,
          /**< [in] If non-zero, assigns global ids to vertices, according to global parameterization. */
        int elem_gids,
          /**< [in] If non-zero, assigns global ids to elements, according to global parameterization. */
        iBase_EntitySetHandle* set_handle, 
          /**< [inout] A set to which the underlying ScdBox set will be added.  NULL if not desired. 
           *           If *NULL, will be set directly to the underlying ScdBox's set. */
        int *err
          /**< [out] Error flag. */
);
/***************************************************************************//**
 * \brief  Free memory allocated with malloc
 *
 ******************************************************************************/

void iMesh_freeMemory(
        iMesh_Instance instance,
          /**< [in] iMesh instance handle */
         void ** ptrToMem);

/***************************************************************************//**
 * \defgroup ScdMesh Structured Mesh
 * \ingroup iMeshExtensions
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iMeshExtensions iMesh Extensions
 * \ingroup iMesh
 ******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif
