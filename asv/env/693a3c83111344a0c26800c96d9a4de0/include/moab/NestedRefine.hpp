/*! \file NestedRefine.hpp
 * This class implements the generation of a hierarchy of meshes via uniform refinement from an input mesh.
 *  The internal upper bound on the number of levels is set to 20.
 */


#ifndef NESTED_REFINE_HPP
#define NESTED_REFINE_HPP

#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include <map>

namespace moab
{
  
#define MAX_DEGREE 3
#define MAX_VERTS 64
#define MAX_CHILDRENS 27
#define MAX_HE 12
#define MAX_HF 6
#define MAX_CONN 8
#define MAX_VHF 20
#define MAX_LEVELS 20


  class Core;
  class HalfFacetRep;
  class ParallelComm;
  class CpuTimer;

  class NestedRefine
  {
    
  public:

    NestedRefine(Core *impl,  ParallelComm *comm=0, EntityHandle rset=0);

    ~NestedRefine();
    
    ErrorCode initialize();

    //! \brief Generate a mesh hierarchy.
    /** Given a mesh, generate a sequence of meshes via uniform refinement. The inputs are: a) an array(level_degrees) storing the degrees which will be
       *  used to refine the previous level mesh to generate a new level and b) the total number of levels(should be same length as that of the
       *  array in a). Each mesh level in the hierarchy are stored in different meshsets whose handles are returned after the hierarchy generation.
       *  These handles can be used to work with a specific mesh level.
       * \param level_degrees Integer array storing the degrees used in each level.
       * \param num_level The total number of levels in the hierarchy.
       * \param hm_set EntityHandle STL vector that returns the handles of the sets created for each mesh level.
      */

    ErrorCode generate_mesh_hierarchy( int num_level, int *level_degrees,  std::vector<EntityHandle>& level_sets);

    //! Given an entity and its level, return its connectivity.
    /** Given an entity at a certain level, it finds the connectivity via direct access to a stored internal pointer to the memory to connectivity sequence for the given level.
       * \param ent EntityHandle of the entity
       * \param level Integer level of the entity for which connectivity is requested
       * \param conn std::vector returning the connectivity of the entity
      */

    ErrorCode get_connectivity(EntityHandle ent, int level, std::vector<EntityHandle> &conn);

    //! Given a vector of vertices and their level, return its coordinates.
    /** Given a vector of vertices at a certain level, it finds the coordinates via direct access to a stored internal pointer to the memory to coordinate sequence for the given level.
       * \param verts std::vector of the entity handles of the vertices
       * \param num_verts The number of vertices
       * \param level Integer level of the entity for which connectivity is requested
       * \param coords double pointer returning the coordinates of the vertices
      */

    ErrorCode get_coordinates(EntityHandle *verts, int num_verts,  int level, double *coords);

    //! Get the adjacencies associated with an entity.
    /** Given an entity of dimension <em>d</em>, gather all the adjacent <em>D</em> dimensional entities where <em>D >, = , < d </em>.
       *
       * \param source_entity EntityHandle to which adjacent entities have to be found.
       * \param target_dimension Int Dimension of the desired adjacent entities.
       * \param target_entities Vector in which the adjacent EntityHandle are returned.
       */

    ErrorCode get_adjacencies(const EntityHandle source_entity,
                              const unsigned int target_dimension,
                              std::vector<EntityHandle> &target_entities);


    //Interlevel parent-child or vice-versa queries
    /** Given an entity from a certain level, it returns a pointer to its parent at the requested parent level.
      *  NOTE: This query does not support vertices.
      *
      * \param child EntityHandle of the entity whose parent is requested
      * \param child_level Mesh level where the child exists
      * \param parent_level Mesh level from which parent is requested
      * \param parent Pointer to the parent in the requested parent_level
    */

    ErrorCode child_to_parent(EntityHandle child, int child_level, int parent_level, EntityHandle *parent);

    /** Given an entity from a certain level, it returns a std::vector of all its children from the requested child level.
      *  NOTE: This query does not support vertices.
      *
      * \param parent EntityHandle of the entity whose children in subsequent level is requested
      * \param parent_level Mesh level where the parent exists
      * \param child_level Mesh level from which its children are requested
      * \param children Vector containing all childrens from the requested child_level
    */

    ErrorCode parent_to_child(EntityHandle parent, int parent_level, int child_level,  std::vector<EntityHandle> &children);

    /** Given a vertex from a certain level, it returns a std::vector of all entities from the previous level that contains it.
      *
      * \param vertex EntityHandle of the vertex
      * \param level Mesh level of the vertex
      * \param incident_entities Vector containing entities from the previous level incident on the vertex
    */

    ErrorCode vertex_to_entities(EntityHandle vertex, int level, std::vector<EntityHandle> &incident_entities);

    /** Given an entity at a certain level, it returns a boolean value true if it lies on the domain boundary.
        * \param entity
      */

      bool is_entity_on_boundary(const EntityHandle &entity);

      ErrorCode exchange_ghosts(std::vector<EntityHandle> &lsets, int num_glayers);

      struct codeperf{
        double tm_total;
        double tm_refine;
        double tm_resolve;
      };

      codeperf timeall;

  protected:
    Core *mbImpl;
    ParallelComm *pcomm;
    HalfFacetRep *ahf;
    CpuTimer *tm;

    EntityHandle _rset;
    Range _inverts, _inedges, _infaces, _incells;

    EntityType elementype;
    int meshdim;
    int level_dsequence[MAX_LEVELS];
    std::map<int,int> deg_index;
    bool hasghost;

    /*! \struct refPatterns
     * Refinement patterns w.r.t the reference element. It consists of a locally indexed vertex list along with their natural coordinates, the connectivity of the subdivided entities with local indices, their local AHF maps along with other helper fields to aid in general book keeping such as avoiding vertex duplication during refinement. The entity and degree specific values are stored in the Templates.hpp.
     * \sa Templates.hpp
     */

    //! refPatterns
    struct refPatterns{
      //! Number of new vertices on edge
      short int nv_edge;
      //! Number of new vertices on face, does not include those on edge
      short int nv_face;
      //! Number of new vertices in cell
      short int nv_cell;
      //! Total number of new vertices per entity
      short int total_new_verts;
      //! Total number of new child entities
      short int total_new_ents;
      //! Lower and upper indices of the new vertices
      int vert_index_bnds[2];
      //! Natural coordinates of the new vertices w.r.t reference
      double vert_nat_coord[MAX_VERTS][3];
      //! Connectivity of the new entities
      int ents_conn[MAX_CHILDRENS][MAX_CONN];
      //! Vertex to half-facet map of the new vertices
      int v2hf[MAX_VERTS][2];
      //! Opposite half-facet map of the new entities
      int ents_opphfs[MAX_CHILDRENS][2*MAX_CONN];
      //! Helper: storing the local ids of vertices on each local edge
      int vert_on_edges[MAX_HE][MAX_VHF];
      //!  Helper: storing local ids of verts on each local face, doesnt include those on edges of the face
      int vert_on_faces[MAX_HF][MAX_VHF];
      //! Helper: stores child half-facets incident on parent half-facet. First column contain the number of such children
      int ents_on_pent[MAX_HF][MAX_CHILDRENS];
    };
    //! refPatterns

    static const refPatterns refTemplates[9][MAX_DEGREE];

    //! Helper
    struct intFEdge{
      //! Number of edges interior to a face
      short int nie;
      //! Local connectivity of the interior edges
      int ieconn[12][2];
    };
    static const intFEdge intFacEdg[2][2];

    int get_index_from_degree(int degree);

    // HM Storage Helper
    struct level_memory{
      int num_verts, num_edges, num_faces, num_cells;
      EntityHandle start_vertex, start_edge, start_face, start_cell;
      std::vector<double *> coordinates;
      EntityHandle *edge_conn, *face_conn, *cell_conn;
      Range verts, edges, faces, cells;
    };

    level_memory level_mesh[MAX_LEVELS];

    //Basic Functions

    //Estimate and create storage for the levels
    ErrorCode estimate_hm_storage(EntityHandle set, int level_degree, int cur_level, int hmest[4]);
    ErrorCode create_hm_storage_single_level(EntityHandle *set, int cur_level, int estL[4]);

    //Generate HM : Construct the hierarchical mesh: 1D, 2D, 3D
    ErrorCode generate_hm(int *level_degrees, int num_level, EntityHandle *hm_set);
    ErrorCode construct_hm_entities(int cur_level, int deg);
    ErrorCode construct_hm_1D(int cur_level, int deg);
    ErrorCode construct_hm_1D(int cur_level, int deg, EntityType type, std::vector<EntityHandle> &trackverts);
    ErrorCode construct_hm_2D(int cur_level, int deg);
    ErrorCode construct_hm_2D(int cur_level, int deg, EntityType type, std::vector<EntityHandle> &trackvertsC_edg, std::vector<EntityHandle> &trackvertsF);
    ErrorCode construct_hm_3D(int cur_level, int deg);

    ErrorCode subdivide_cells(EntityType type,int cur_level, int deg);
    ErrorCode subdivide_tets(int cur_level, int deg);

    // General helper functions
    ErrorCode copy_vertices_from_prev_level(int cur_level);
    ErrorCode count_subentities(EntityHandle set, int cur_level, int *nedges, int *nfaces);
    ErrorCode get_octahedron_corner_coords(int cur_level, int deg, EntityHandle *vbuffer, double * ocoords);
    int find_shortest_diagonal_octahedron(int cur_level, int deg, EntityHandle *vbuffer);
    int get_local_vid(EntityHandle vid, EntityHandle ent, int level);

    // Book-keeping functions
    ErrorCode update_tracking_verts(EntityHandle cid, int cur_level, int deg, std::vector<EntityHandle> &trackvertsC_edg, std::vector<EntityHandle> &trackvertsC_face, EntityHandle *vbuffer);
    ErrorCode reorder_indices(int cur_level, int deg, EntityHandle cell, int lfid, EntityHandle sib_cell, int sib_lfid, int index, int *id_sib);
    ErrorCode reorder_indices(int deg, EntityHandle *face1_conn, EntityHandle *face2_conn, int nvF, std::vector<int> &lemap, std::vector<int> &vidx, int *leorient=NULL);
    ErrorCode get_lid_inci_child(EntityType type, int deg, int lfid, int leid, std::vector<int> &child_ids, std::vector<int> &child_lvids);

    //Permutation matrices
    struct pmat{
      short int num_comb; // Number of combinations
      int comb[MAX_HE][MAX_HE]; //Combinations
      int lemap[MAX_HE][MAX_HE]; //Local edge map
      int orient[MAX_HE]; //Orientation
      int porder2[MAX_HE][MAX_HE]; // Permuted order degree 2
      int porder3[MAX_HE][MAX_HE]; // Permuted order degree 3
    };

    static const pmat permutation[2];

    // Print functions
    ErrorCode print_maps_1D(int level);
    ErrorCode print_maps_2D(int level, EntityType type);
    ErrorCode print_maps_3D(int level, EntityType type);

    // Coordinates
    ErrorCode compute_coordinates(int cur_level, int deg, EntityType type, EntityHandle *vbuffer, int vtotal, double *corner_coords, std::vector<int> &vflag, int nverts_prev);

    // Update the ahf maps

    ErrorCode update_local_ahf(int deg, EntityType type, EntityHandle *vbuffer, EntityHandle *ent_buffer, int etotal);

    ErrorCode update_local_ahf(int deg, EntityType type, int pat_id, EntityHandle *vbuffer, EntityHandle *ent_buffer, int etotal);

    ErrorCode update_global_ahf(EntityType type, int cur_level, int deg);

    ErrorCode update_global_ahf(int cur_level, int deg, std::vector<int> &pattern_ids);

    ErrorCode update_global_ahf_1D(int cur_level, int deg);

    ErrorCode update_global_ahf_1D_sub(int cur_level, int deg);

    ErrorCode update_ahf_1D(int cur_level);

    ErrorCode update_global_ahf_2D(int cur_level, int deg);

    ErrorCode update_global_ahf_2D_sub(int cur_level, int deg);

    ErrorCode update_global_ahf_3D(int cur_level, int deg);

    ErrorCode update_global_ahf_3D(int cur_level, int deg, std::vector<int> &pattern_ids);

    /** Boundary extraction functions
        * Given a vertex at a certain level, it returns a boolean value true if it lies on the domain boundary.
        * Note: This is a specialization of the NestedRefine::is_entity_on_boundary function and applies only to vertex queries.
        * \param entity
      */
    bool is_vertex_on_boundary(const EntityHandle& entity);
    bool is_edge_on_boundary(const EntityHandle& entity);
    bool is_face_on_boundary(const EntityHandle& entity);
    bool is_cell_on_boundary(const EntityHandle& entity);

    //ErrorCode find_skin_faces(EntityHandle set, int level, int nskinF);

  };
} //name space moab
#endif
