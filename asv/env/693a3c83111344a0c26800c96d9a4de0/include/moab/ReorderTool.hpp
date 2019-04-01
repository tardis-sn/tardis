/** \file   ReorderTool.hpp
 *  \author Jason Kraftcheck 
 *  \date   2011-05-23
 */

#ifndef moab_REORDER_TOOL_HPP
#define moab_REORDER_TOOL_HPP

#include "moab/Types.hpp"
#include <vector>

namespace moab {

class Core;
class Range;

class ReorderTool 
{
  public:
  
    ReorderTool( Core* moab ) : mMB(moab) {}
    

    /**\brief Calculate new handle order by tag value.
     *
     * Given a tag containing integer values, calculate new order for entities 
     * in the database (except entity sets) such that all entities
     * with tag value A occur before all entities with tag value B
     * in the handle space where A < B.  Ordering will be stable for
     * entities with the same tag value.
     *
     *\param ordering_tag  Sinlge integer tag, where value on each entity
     *                     determines the new position in the handle ordering.
     *                     Entities may have duplicate values.  
     *\param ordering_tag_skip_value Do not reorder entities with this tag
     *                     value. This is typically the default value of
     *                     ordering_tag.  Specifying this limits the re-ordering
     *                     to only those entities for which ordering_tag has 
     *                     been set.
     *\param new_handle_tag_out  Passed back new tag handle containing the 
     *                     entity mapping.  The returned tag will be anonymous.
     *                     The caller is responsible for releasing the tag.  
     *                     The value of this tag on each handle is the new 
     *                     handle that the entity will will be moved to. The 
     *                     tag value will be zero for entities that were not 
     *                     re-ordered.
     */
    ErrorCode handle_order_from_int_tag( Tag ordering_tag,
                                         int ordering_tag_skip_value,
                                         Tag& new_handle_tag_out );
    
    /**\brief Calculate new handle order by tag value.
     *
     * Given a tag containing integer values, calculate new order for entities 
     * in the database (except entity sets) such that all entities
     * with tag value A occur before all entities with tag value B
     * in the handle space where A < B.  Ordering will be stable for
     * entities with the same tag value.
     *
     *\param type          Entity type for which to calculate re-ordering.
     *\param vals_per_ent  Zero for vertices.  Connectivity length for elements.
     *\param ordering_tag  Sinlge integer tag, where value on each entity
     *                     determines the new position in the handle ordering.
     *                     Entities may have duplicate values.  
     *\param ordering_tag_skip_value Do not reorder entities with this tag
     *                     value. This is typically the default value of
     *                     ordering_tag.  Specifying this limits the re-ordering
     *                     to only those entities for which ordering_tag has 
     *                     been set.
     *\param new_handle_tag Tag into which to store new handle for each
     *                     entity.  Tag must be defined to store a single
     *                     entity handle and must have a default value of 
     *                     zero.
     */
    ErrorCode handle_order_from_int_tag( EntityType type,
                                         int vals_per_ent,
                                         Tag ordering_tag,
                                         int ordering_tag_skip_value,
                                         Tag new_handle_tag );
    
    
    /**\brief Calculate new handle order by set containment
     *
     * Given a list of sets, re-order entities such that handles are
     * grouped contiguously by set.  Will also group all adjacent mesh
     * entities, such that entities that are are adjacent to members of
     * two or more of the input sets will ge grouped by the combined
     * list of sets (e.g. if the input sets contain elements then all
     * vertices that are adjacent to only elements in the first two sets
     * will be grouped together).
     *
     *\param sets          Entity sets by which to group entities.  
     *\param new_handle_tag_out  Passed back new tag handle containing the 
     *                     entity mapping.  The returned tag will be anonymous.
     *                     The caller is responsible for releasing the tag.  
     *                     The value of this tag on each handle is the new 
     *                     handle that the entity will will be moved to. The 
     *                     tag value will be zero for entities that were not 
     *                     re-ordered.
     */
    ErrorCode handle_order_from_sets_and_adj( const Range& sets,
                                              Tag& new_handle_tag_out );
    
    /**\brief Do the re-ordering indicated by the passed handle tag.
     *
     * The specified re-ordering must be a permutation.  Each existing
     * entity must be moved to a new, existing handle such that no
     * two entities are moved to the same new handle.
     *
     * Given a tag storing handles that define a permution, apply the
     * described re-ordering.  The passed tag must contain one entity
     * handle per entity.  The value of the tag must be zero for all
     * entities that are not to be re-ordered.  For entities to be 
     * re-ordered, the tag must contain the new handle that the entity
     * is to be moved to.  No two entities may have the same value for
     * this tag (other than a value of zero.)
     *
     *\param new_handle_tag Tag containing new handles for entities to
     *                      reorder.  Typically the output of 
     *                      handle_order_from_int_tag or similar.
     */
    ErrorCode reorder_entities( Tag new_handle_tag );

  private:
   
    /**\brief helper function for reorder_entities
     *
     * Reorder tag data for all entities of specified type.
     *
     * Also updates tag values for MB_TYPE_HANDLE tags.
     *\param type  Entity type to reorder
     *\param new_handles Tag containing old->new handle mapping
     *\param reorder_tag The tag data to reorder
     */
    ErrorCode reorder_tag_data( EntityType type, Tag new_handles, Tag reorder_tag );
    
    /**\brief helper function for reorder_entities
     *
     * Update set contents for changed handles.
     *\param new_handles Tag containing old->new handle mapping
     */
    ErrorCode update_set_contents( Tag new_handles );
    
    /**\brief Get all entities of specified type and size
     *\param t the type of entity to retreive
     *\param vals_per_ent entity size (connectivity length for elements,
     *         dimension for vertices)
     */
    void get_entities( EntityType t, int vals_per_ent, Range& result );
    
    /**\brief Get new handles corresponding to old handles
     *\param tag Tag containing old->new mapping
     */
    ErrorCode get_reordered_handles( Tag tag, 
                                     const Range& old_handles, 
                                     std::vector<EntityHandle>& new_handles );
    
    /**\brief Get new handles corresponding to old handles
     *\param tag Tag containing old->new mapping
     */
    ErrorCode get_reordered_handles( Tag tag, 
                                     const std::vector<EntityHandle>& old_handles,
                                     std::vector<EntityHandle>& new_handles );
    
    /**\brief Get new handles corresponding to old handles
     *\param tag Tag containing old->new mapping
     */
    ErrorCode get_reordered_handles( Tag tag, 
                                     const EntityHandle* old_handles,
                                     EntityHandle* new_handles,
                                     size_t num_handles );
    
    /**\brief Remove any non-ordered handles and return new handles for remaining
     *\param tag Tag containing old->new mapping
     */
    ErrorCode get_new_handles( Tag tag,
                               Range& old_handles,
                               std::vector<EntityHandle>& newhandles );
    
    /**\brief convert from input for \chandle_order_from_sets_and_adj to
     *        \c input for handle_order_from_int_tag
     */
    ErrorCode int_order_from_sets_and_adj( const Range& sets,
                                           Tag order_tag,
                                           int skip_val,
                                           std::vector<std::vector<EntityHandle>*>& data );
    
    Core* mMB;
};



} // namespace moab

#endif // moab_REORDER_TOOL_HPP
