
#include "MBEntityType.h"
#include "MBCN_protos.h"

#ifdef __cplusplus
extern "C" {
#endif    

    void MBCN_GetBasis(int *rval);
  
    void MBCN_SetBasis(const int in_basis);

    void MBCN_EntityTypeName(const int this_type, char *rval, int rval_len);
  
    void MBCN_EntityTypeFromName(const char *name, int *rval);
  
    void MBCN_Dimension(const int t, int *rval);

    void MBCN_VerticesPerEntity(const int t, int *rval);
  
    void MBCN_NumSubEntities(const int t, const int d, int *rval);

    void MBCN_SubEntityType(const int this_type,
                            const int sub_dimension,
                            const int index, int *rval);
  
    void MBCN_SubEntityVertexIndices(const int this_type, 
                                     const int sub_dimension,
                                     const int sub_index,
                                     int sub_entity_conn[]);

    void MBCN_AdjacentSubEntities(const int this_type,
                                  const int *source_indices,
                                  const int num_source_indices,
                                  const int source_dim,
                                  const int target_dim,
                                  int *index_list,
                                  int *num_indices,
                                  const int operation_type, int *rval);

    void MBCN_SideNumberInt(const int *parent_conn, const MBEntityType parent_type,
                            const int *child_conn, const int child_num_verts,
                            const int child_dim,
                            int *side_no, int *sense, int *offset);

    void MBCN_SideNumberUint(const unsigned int *parent_conn, const MBEntityType parent_type,
                             const unsigned int *child_conn, const int child_num_verts,
                             const int child_dim,
                             int *side_no, int *sense, int *offset);

    void MBCN_SideNumberLong(const long *parent_conn, const MBEntityType parent_type,
                             const long *child_conn, const int child_num_verts,
                             const int child_dim,
                             int *side_no, int *sense, int *offset);

    void MBCN_SideNumberUlong(const unsigned long *parent_conn, const MBEntityType parent_type,
                              const unsigned long *child_conn, const int child_num_verts,
                              const int child_dim,
                              int *side_no, int *sense, int *offset);

    void MBCN_SideNumberVoid(void * const *parent_conn, const MBEntityType parent_type,
                             void * const *child_conn, const int child_num_verts,
                             const int child_dim,
                             int *side_no, int *sense, int *offset);

    void MBCN_SideNumber(const int parent_type,
                         const int *child_conn_indices, const int child_num_verts,
                         const int child_dim,
                         int *side_no, int *sense, int *offset);

    void MBCN_OppositeSide(const int parent_type,
                           const int child_index,
                           const int child_dim,
                           int *opposite_index,
                           int *opposite_dim, int *rval);

    void MBCN_ConnectivityMatchInt(const int *conn1,
                                   const int *conn2,
                                   const int num_vertices,
                                   int *direct, int *offset, int *rval);
    void MBCN_ConnectivityMatchUint(const unsigned int *conn1,
                                    const unsigned int *conn2,
                                    const int num_vertices,
                                    int *direct, int *offset, 
                                    int *rval);
    void MBCN_ConnectivityMatchLong(const long* conn1,
                                    const long* conn2,
                                    const int num_vertices,
                                    int* direct, int* offset , int *rval);
    void MBCN_ConnectivityMatchUlong(const unsigned long* conn1,
                                     const unsigned long* conn2,
                                     const int num_vertices,
                                     int *direct, int* offset,
                                     int *rval);
    void MBCN_ConnectivityMatchVoid(void* const* conn1,
                                    void* const* conn2,
                                    const int num_vertices,
                                    int* direct, int* offset , int *rval);

    void MBCN_setPermutation(const MBEntityType t, const int dim, int *pvec, 
                             const int num_entries, const int is_reverse);
  
    void MBCN_resetPermutation(const MBEntityType t, const int dim);

    void MBCN_permuteThisInt(const MBEntityType t, const int dim, int *pvec, 
                             const int num_indices, const int num_entries, int *rval);

    void MBCN_permuteThisUint(const MBEntityType t, const int dim, unsigned int *pvec, 
                              const int num_indices, const int num_entries, int *rval);

    void MBCN_permuteThisLong(const MBEntityType t, const int dim, long *pvec, 
                              const int num_indices, const int num_entries, int *rval);

    void MBCN_permuteThisVoid(const MBEntityType t, const int dim, void **pvec, 
                              const int num_indices, const int num_entries, int *rval);
  
    void MBCN_revPermuteThisInt(const MBEntityType t, const int dim, int *pvec, 
                             const int num_indices, const int num_entries, int *rval);

    void MBCN_revPermuteThisUint(const MBEntityType t, const int dim, unsigned int *pvec, 
                              const int num_indices, const int num_entries, int *rval);

    void MBCN_revPermuteThisLong(const MBEntityType t, const int dim, long *pvec, 
                              const int num_indices, const int num_entries, int *rval);

    void MBCN_revPermuteThisVoid(const MBEntityType t, const int dim, void **pvec, 
                              const int num_indices, const int num_entries, int *rval);

    void MBCN_HasMidEdgeNodes(const int this_type, 
                              const int num_verts, int *rval);

    void MBCN_HasMidFaceNodes(const int this_type, 
                              const int num_verts, int *rval);

    void MBCN_HasMidRegionNodes(const int this_type, 
                                const int num_verts, int *rval);

    void MBCN_HasMidNodes(const int this_type, 
                          const int num_verts, 
                          int mid_nodes[4]);

    void MBCN_HONodeParent( int elem_type,
                            int num_nodes, 
                            int ho_node_index,
                            int *parent_dim, 
                            int *parent_index );

    void MBCN_HONodeIndex(const int this_type, const int num_verts,
                          const int subfacet_dim, const int subfacet_index, int *rval);

#ifdef __cplusplus
}
#endif
