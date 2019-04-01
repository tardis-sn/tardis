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


#ifndef MOAB_HALF_FACET_REP_HPP
#define MOAB_HALF_FACET_REP_HPP

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"

namespace moab
{

/*! 
 *  \brief   HalfFacetRep class implements the Array-Based Half-Facet(AHF) Mesh data structure on top of MOAB.
 *  \        It is based on a generalized notion of a half-facet derived from that a half-edge/half-face data structure for 2D/3D.
 *  \        The core construct of AHF are its two maps: sibling-half-facets(SIBHFS) and vertex-to-half-facet(V2HF)
 *  \        1. SIBHFS: Maps every half-facet in the mesh to its sibling half-facets,
 *  \        2. V2HF: Maps each vertex to an incident half-facet
 *  \        Using these two maps, a range of adjacency queries is performed. The maps are stored in dense tags over entities and vertices.
 *  \
 *  \        Adjacency functions:
 *  \        1. upward-incidence queries: vertex -> edge, edge -> faces, edge -> cells, face ->cells
 *  \        2. neighborhood (same-dimensional) adjacency queries: edge -> edges, face -> faces, cell -> cells, etc.
 *  \        3. downward adjacency queries: face -> edges, cell -> edges, etc.
 *  \
 *  \        Mesh types supported: 
 *  \        1D(edges), 2D(triangles, quads), 3D(tet, pyramid, prism, hex), Mixed dimensional meshes
 *  \
 *  \        CURRENTLY NOT SUPPORTED:
 *  \        1. Meshes with mixed entity types of same dimension. Ex. a volume mesh with both tets and prisms.
 *  \        2. create_if_missing = true
 *  \        3. Modified meshes
 *  \
 */ 

 typedef EntityHandle HFacet;

const int MAXSIZE = 200;

//! ENUM for the type of input mesh.
enum MESHTYPE{
    CURVE = 0, //Homogeneous curve mesh
    SURFACE, // Homogeneous surface mesh
    SURFACE_MIXED, // Mixed surface with embedded curves
    VOLUME, // Homogeneous volume mesh
    VOLUME_MIXED_1, // Volume mesh with embedded curves
    VOLUME_MIXED_2, // Volume mesh with embedded surface
    VOLUME_MIXED //Volume mesh with embedded curves and surfaces
};

enum {
  MAX_VERTICES = 8,
  MAX_EDGES = 12,
  MAX_FACES = 6,
  MAX_VERTS_HF = 4,
  MAX_INCIDENT_HF = 4
};

class Core;
class ParallelComm;
class HalfFacetRep{

public:

  HalfFacetRep(Core *impl,  ParallelComm *comm=0, moab::EntityHandle rset=0);

  ~HalfFacetRep();

  bool check_mixed_entity_type();

  // User interface functions

  //! Constructs the sibling-half-facet and vertex-to-incident-half-facet maps for each dimension present in the input. This routine should be called before any calls for adjacency is made.
  ErrorCode initialize();

  //! Deinitialize
  ErrorCode deinitialize();

  //! Prints the tag values.
  ErrorCode print_tags(int dim);

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


  //! Get the upward incidences associated with an entity.
  /** Given an entity of dimension <em>d</em>, gather all the incident <em>D(>d)</em> dimensional entities.
     * :
     * \param ent EntityHandle to which incident entities have to be found.
     * \param out_dim Dimension of the desired incidence information.
     * \param adjents Vector in which the incident entities are returned.
     * \param local_id Set to false by default. If true, returns the local id's of the half-facets
     * \param lids Vector in which the local id's are returned.
     */

  ErrorCode get_up_adjacencies(EntityHandle ent,
                               int out_dim,
                               std::vector<EntityHandle> &adjents,
                               std::vector<int> * lids = NULL );

  //! Get the same-dimensional entities connected with an entity.
  /** Given an entity of dimension <em>d</em>, gather all the entities connected via <em>d-1</em> dimensional entities.
     *  Same as bridge_adjacencies in MOAB.
     *
     * \param ent EntityHandle to which neighbor entities have to be found.
     * \param adjents Vector in which the neighbor entities are returned.
     */

  ErrorCode get_neighbor_adjacencies(EntityHandle ent,
                                     std::vector<EntityHandle> &adjents);

  //! Get the downward adjacent entities connected with an entity.
  /** Given an entity of dimension <em>d</em>, gather all the <em>d-1</em> dimensional entities.
     *
     * \param ent EntityHandle to which neighbor entities have to be found.
     * \param out_dim Dimension of the desired downward adjacency.
     * \param adjents Vector in which the neighbor entities are returned.
     */

  ErrorCode get_down_adjacencies(EntityHandle ent, int out_dim, std::vector<EntityHandle> &adjents);


  // 1D Maps and queries

  //! Given a range of edges, determines the map for sibling half-verts and stores them into SIBHVS_EID, SIBHVS_LVID tags.
  /** Compute all sibling half-vertices for all half-vertices in the given curve. The sibling half-verts is
     *  defined in terms of the containing edge and the local id of the vertex w.r.t that edge.
     *  That is, the map consists of two pieces of information: <EntityHandle eid, int lvid>
     *
     * \param edges Range of edges.
    */

  ErrorCode determine_sibling_halfverts(Range &verts, Range &edges);

  //! Given a range of edges, determines the map for incident half-verts and stores them into V2HV_EID, V2HV_LVID tags.
  /** Compute a map between a vertex and an incident half-vertex. This map is not always required, but is
     * essential for local neighborhood searching as it acts like an anchor to start the search.
     *
     * \param edges Range of edges
    */

  ErrorCode determine_incident_halfverts(Range &edges);

  //! Given a vertex, finds the edges incident on it.
  /** Given a vertex handle, it starts by first finding an incident half-vert by using the incident
     * half-vert map, and then obtaining all the sibling half-verts of the corresponding half-vertex.
     *
     * \param vid EntityHandle of the query vertex
     * \param adjents Vector returning the incident edges
     * \param local_id False by default. If true, returns the local vertex id's corresponding to vid
     * \param lvids Vector returning the local vertex id's
    */

  ErrorCode get_up_adjacencies_1d(EntityHandle vid,
                                  std::vector<EntityHandle> &adjents,
                                  std::vector<int> * lvids = NULL);

  //! Given an edge, finds vertex-connected neighbor edges
  /** Given an edge, it gathers all the incident edges of each vertex of the edge.
     *
     * \param eid EntityHandle of the query edge
     * \param adjents Vector returning neighbor edges
    */

  ErrorCode get_neighbor_adjacencies_1d(EntityHandle eid,
                                        std::vector<EntityHandle> &adjents);


  // 2D Maps and queries

  //! Given a range of faces, determines the sibling half-edges and stores them into SIBHES_FID, SIBHES_LEID tags.
  /** Compute all sibling half-edges for all half-edges in the given surface.
     * The sibling half-edges is defined in terms of the containing face and the local id of the edge w.r.t that entity.
     * That is, the map consists of two pieces of information: <EntityHandle fid, int leid>
     *
     * \param faces Range of faces
    */

  ErrorCode determine_sibling_halfedges(Range &faces);

  //! Given a range of faces, determines the incident half-edges and stores them into V2HE_FID, V2HE_LEID tags.
  /** Compute a map between a vertex and an incident half-edge.
     * This map is not always required, but is essential for local neighborhood searching as it acts
     * like an anchor to start the search.
     *
     * \param faces Range of faces
    */

  ErrorCode determine_incident_halfedges(Range &faces);

  //! Given a vertex, finds the faces incident on it.
  /** Given a vertex, it first finds an incident half-edge via v2he map, and then
     * collects all the incident half-edges/faces via the sibhes map.
     *
     * \param vid EntityHandle of the query vertex
     * \param adjents Vector returning the incident faces
    */

  ErrorCode get_up_adjacencies_vert_2d(EntityHandle vid, std::vector<EntityHandle> &adjents);

  //! Given an edge, finds the faces incident on it.
  /** Given an edge, it first finds a matching half-edge corresponding to eid, and then
     * collects all the incident half-edges/faces via the sibhes map.
     *
     * \param eid EntityHandle of the query edge
     * \param adjents Vector returning the incident faces
     * \param local_id By default false. If true, returns the local edge id's corresponding to the input edge
     * \param leids Vector returning local edge ids
    */

  ErrorCode get_up_adjacencies_2d(EntityHandle eid,
                                  std::vector<EntityHandle> &adjents,
                                  std::vector<int> * leids = NULL);

  //! Given a half-edge <fid, leid>, finds the faces incident on it.
  /**
     *
     * \param fid EntityHandle of the containing face
     * \param leid local id of the edge w.r.t to the face
     * \param add_inent If true, adds the input fid into the returning vector of adjents.
     * \param adjents Vector returning the incident faces
     * \param local_id By default false. If true, returns the local edge id's as well.
     * \param leids Vector returning local edge ids
    */

  ErrorCode get_up_adjacencies_2d(EntityHandle fid,
                                  int leid,
                                  bool add_inent,
                                  std::vector<EntityHandle> &adj_ents,
                                  std::vector<int>  *adj_leids = NULL, std::vector<int> *adj_orients = NULL);

  //! Given an edge, finds edge-connected neighbor face
  /** Given an face, it gathers all the neighbor faces of each local edge of the face.
     *
     * \param fid EntityHandle of the query face
     * \param adjents Vector returning neighbor faces
    */

  ErrorCode get_neighbor_adjacencies_2d(EntityHandle fid,
                                        std::vector<EntityHandle> &adjents);

  //! Given a face, finds its edges.
  /** Given a face, it first finds incident edges on each vertex of the face, and then
     *  it performs a set intersection to gather all the edges of the given face.
     *
     * \param fid EntityHandle of the query face
     * \param adjents Vector returning its edges
    */

  ErrorCode get_down_adjacencies_2d(EntityHandle fid,
                                    std::vector<EntityHandle> &adjents);

  //! Given a range of faces, finds the total number of edges.

  int find_total_edges_2d(Range &faces);

  // 3D Maps and queries

  //! Given a range of cells, determines the sibling half-faces and stores them into SIBHFS_CID, SIBHFS_LFID tags.
  /** Compute all sibling half-faces for all half-faces in the given volume.
     * The sibling half-faces is defined in terms of the containing cell and the local id of the face w.r.t that cell.
     * That is, the map consists of two pieces of information: <EntityHandle cid, int lfid>
     *
     * \param faces Range of cells
    */

  ErrorCode determine_sibling_halffaces( Range &cells);

  //! Given a range of cells, determines the incident half-faces and stores them into V2HF_CID, V2HF_LFID tags.
  /** Compute a map between a vertex and an incident half-face.
     * This map is not always required, but is essential for local neighborhood searching as it acts
     * like an anchor to start the search.
     *
     * \param faces Range of cells
    */

  ErrorCode determine_incident_halffaces( Range &cells);

  //! Given a range of cells, tags all border vertices with a true value.
  /** Tag border vertices by using the sibhf_cid map. All vertices on half-faces with no sibling
     * half-faces are considered as border vertices.
     *
     * \param cells Range of cells
     * \param isborder: A dense tag over all vertices of size 1. Value is true for a border vertex, otherwise is false.
    */

  ErrorCode determine_border_vertices( Range &cells,
                                       Tag isborder);

  //! Given a vertex, finds the cells incident on it.
  /** Given a vertex, it first finds an incident half-face via v2hf map, and then
     * collects all the incident half-faces via the sibhfs map.
     *
     * \param vid EntityHandle of the query vertex
     * \param adjents Vector returning the incident cells
    */

  ErrorCode get_up_adjacencies_vert_3d(EntityHandle vid, std::vector<EntityHandle> &adjents);

  //! Given an edge, finds the cells incident on it.
  /** Given an edge, it first finds a matching local edge in a cell corresponding to eid, and then
     * collects all the incident cells via the sibhfs map.
     *
     * \param eid EntityHandle of the query edge
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local edge id's corresponding to the input edge
     * \param leids Vector returning local edge ids
    */

  ErrorCode get_up_adjacencies_edg_3d(EntityHandle eid,
                                      std::vector<EntityHandle> &adjents,
                                      std::vector<int> * leids = NULL);

  //! Given a local edge <cid, leid>, finds the cells incident on it.
  /** Given a local edge, it gathers all the incident cells via the sibhfs map.
     *
     * \param cid EntityHandle of the cell containing the local edge
     * \param leid local edge id w.r.t the cell
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local edge id's corresponding to the input edge
     * \param leids Vector returning local edge ids
    */

  ErrorCode get_up_adjacencies_edg_3d(EntityHandle cid,
                                      int leid, std::vector<EntityHandle> &adjents,
                                      std::vector<int> * leids = NULL, std::vector<int> *adj_orients = NULL);

  //! Given an face, finds the cells incident on it.
  /** Given an face, it first finds a matching half-face in a cell corresponding to face, and then
     * collects all the incident cells via the sibhfs map.
     *
     * \param fid EntityHandle of the query face
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local face id's corresponding to the input face
     * \param leids Vector returning local face ids
    */

  ErrorCode get_up_adjacencies_face_3d(EntityHandle fid,
                                       std::vector<EntityHandle> &adjents,
                                       std::vector<int>  * lfids = NULL);

  //! Given a local face <cid, lfid>, finds the cells incident on it.
  /** Given a local face, it gathers all the incident cells via the sibhfs map.
     *
     * \param cid EntityHandle of the cell containing the local edge
     * \param lfid local face id w.r.t the cell
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local face id's corresponding to the input face
     * \param lfids Vector returning local face ids
    */

  ErrorCode get_up_adjacencies_face_3d(EntityHandle cid,
                                       int lfid,
                                       std::vector<EntityHandle> &adjents,
                                       std::vector<int> * lfids = NULL);

  //! Given a cell, finds face-connected neighbor cells
  /** Given a cell, it gathers all the neighbor cells of each local face of the cell.
     *
     * \param cid EntityHandle of the query cell
     * \param adjents Vector returning neighbor cells
    */

  ErrorCode get_neighbor_adjacencies_3d(EntityHandle cid,
                                        std::vector<EntityHandle> &adjents);

  //! Given a cell, finds its edges.
  /** Given a cell, it first finds incident edges on each vertex of the cell, and then
     *  it performs a set intersection to gather all the edges of the given cell.
     *
     * \param cid EntityHandle of the query cell
     * \param adjents Vector returning its edges
    */

  ErrorCode get_down_adjacencies_edg_3d(EntityHandle cid, std::vector<EntityHandle> &adjents);

  //! Given a cell, finds its faces.
  /** Given a cell, it first finds incident faces on each vertex of the cell, and then
     *  performs a set intersection to gather all the faces of the given cell.
     *
     * \param cid EntityHandle of the query cell
     * \param adjents Vector returning its faces
    */

  ErrorCode get_down_adjacencies_face_3d(EntityHandle cid, std::vector<EntityHandle> &adjents);

  /* Find the number of edges and faces of given range of cells
     * */
  ErrorCode find_total_edges_faces_3d(Range cells, int *nedges, int *nfaces);

  ErrorCode count_subentities(Range &edges, Range &faces, Range &cells, int *nedges, int *nfaces);

  void get_memory_use(unsigned long long& entity_total, unsigned long long& memory_total);

  ErrorCode get_half_facet_in_comp(EntityHandle cid, int leid, std::vector<EntityHandle> &ents, std::vector<int> &lids,  std::vector<int> &lfids);


  /**************************
     *  Interface to AHF maps   *
     **************************/
  HFacet create_halffacet(EntityHandle handle, int lid);

  EntityHandle fid_from_halfacet(const HFacet facet, EntityType type);

  int lid_from_halffacet(const HFacet facet);

  ErrorCode update_entity_ranges(EntityHandle fileset);

  ErrorCode resize_hf_maps(EntityHandle start_vert, int nverts, EntityHandle start_edge, int nedges, EntityHandle start_face, int nfaces, EntityHandle start_cell, int ncells);

  ErrorCode get_sibling_map(EntityType type, EntityHandle ent, EntityHandle *sib_entids, int *sib_lids, int num_halffacets);

  ErrorCode get_sibling_map(EntityType type, EntityHandle ent, int lid, EntityHandle &sib_entid, int &sib_lid);

  ErrorCode set_sibling_map(EntityType type, EntityHandle ent,  EntityHandle *set_entids, int *set_lids, int num_halffacets);

  ErrorCode set_sibling_map(EntityType type, EntityHandle ent, int lid, EntityHandle &set_entid, int &set_lid);

  ErrorCode get_incident_map(EntityType type, EntityHandle vid, std::vector<EntityHandle> &inci_entid, std::vector<int> &inci_lid);

  ErrorCode set_incident_map(EntityType type, EntityHandle vid, std::vector<EntityHandle> &set_entid, std::vector<int> &set_lid);

  /**********************
   *         Local Maps
   * ********************/
  //! 2D local maps
  struct LocalMaps2D{
    //! Number of vertices in a face
    short int num_verts_in_face;
    //! Local ids of the next half-edge
    int next[MAX_INCIDENT_HF];
    //! Local ids of the previous half-edge
    int prev[MAX_INCIDENT_HF];
  };
  static const LocalMaps2D lConnMap2D[2];

//! 3D local maps
  struct LocalMaps3D{
    //! Number of vertices in cell
    short int num_verts_in_cell;
    //! Number of edges in cell
    short int num_edges_in_cell;
    // Number of faces in cell
    short int num_faces_in_cell;
    //! Number of vertices in each half-face
    int hf2v_num[MAX_FACES];
    //! Map: Half-face to local vertex ids
    int hf2v[MAX_FACES][MAX_VERTS_HF];
    //! Number of incident half-faces on each vertex
    int v2hf_num[MAX_VERTICES];
    //! Map: Local vertices to incident half-facets
    int v2hf[MAX_VERTICES][MAX_INCIDENT_HF];
    //!  Map: Local edges to local vertices
    int e2v[MAX_EDGES][2];
    //! Map: Local edges to incident half-faces
    int e2hf[MAX_EDGES][2];
    //! Map: Half-faces to local edges
    int f2leid[MAX_FACES][MAX_VERTS_HF];
    //! Map: local vertices to local edges
    int lookup_leids[MAX_VERTICES][MAX_VERTICES];
    //! Search edge verts:
    int search_everts[5];
    int search_fverts[2];
    int v2le[4][5];
  };
  static const LocalMaps3D lConnMap3D[4];
  MESHTYPE thismeshtype;

   std::map<EntityType, int> cell_index;

  int get_index_in_lmap(EntityHandle cid);
  ErrorCode get_entity_ranges(Range &verts, Range &edges, Range &faces, Range &cells);
  MESHTYPE get_mesh_type(int nverts, int nedges, int nfaces, int ncells);

  EntityHandle *get_rset() { return &_rset; }

protected:

  HalfFacetRep();

  Core *mb;
  ParallelComm *pcomm;
  bool mInitAHFmaps;
  EntityHandle _rset;
  Range _verts, _edges, _faces, _cells;

  //AHF map storage containers for 1D, 2D, 3D
  std::vector<HFacet> sibhvs, v2hv;
  std::vector<HFacet> sibhes, v2he;
  std::vector<HFacet> sibhfs, v2hf;

  //AHF maps for non-manifold vertices 2D, 3D
  std::multimap<EntityHandle, HFacet> v2hes, v2hfs;

  //Auxiliary storage for local searches.
  EntityHandle queue_fid[MAXSIZE], Stkcells[MAXSIZE], cellq[MAXSIZE];
  EntityHandle trackfaces[MAXSIZE], trackcells[MAXSIZE];
  int queue_lid[MAXSIZE];

  struct adj_matrix{
    int val[4][4];
  };

  static const adj_matrix adjMatrix[7];
  int get_index_for_meshtype(MESHTYPE mesh_type);

  // These two flags are for checking mixed entity type meshes
  bool is_mixed;
  bool chk_mixed;

  ErrorCode init_curve();
  ErrorCode init_surface();
  ErrorCode init_volume();

  //! Contains the local information for 2D entities
  /** Given a face, find the face type specific information
     *
     * \param face EntityHandle. Used to gather info about the type of face for which local info is required
     * \param nepf: Returns the number of vertices/edges for given face type.
    */

  //! Contains the local information for 2D entities
  /** Given number of edges, returns local indices of next and previous local edges.
     *
     * \param nepf: The number of vertices/edges for given face type.
     * \param next, prev: Local ids of next and previous edges w.r.t to the face
     *
     * Note: For 2D entities, the number of vertices and edges are same and each local edge is outgoing
     * corresponding to the local vertex, i.e,

            v2        v3 __e2__v2
            /\          |      |
        e2 /  \ e1    e3|      |e1
          /____\        |______|
        v0  e0  v1     v0  e0  v1
    */


    //! Given a half-edge as <he_fid,he_lid> , finds the half-edges incident on it and adds them
    //  to an input queue if it not already added.
    /** Given an half-edge, obtain all the incident half-edges via the sibhes map and add them to a given
     * queue of half-edges, if they do not already exist in the queue. This function is used to increment the
     * search space for finding a matching half-edge.
     *
     * \param he_fid EntityHandle of query half-edge
     * \param he_lid Local id of query half-edge
     */

    ErrorCode get_up_adjacencies_2d(EntityHandle he_fid,
                                    int he_lid,
                                    int *qsize, int *count);

    //! Given an edge, finds a matching half-edge in the surface.
    /** Given an edge eid, it first collects few half-edges belonging to one-ring neighborhood of
     * the starting vertex of the given edge, and then simultaneously searches and adds to the local list
     * of half-edges for searching, till it finds a matching half-edge.
     *
     * \param eid EntityHandle of the query edge
     * \param hefid, helid: Returns the matching half-edge corresponding to the query edge.
    */

    bool find_matching_halfedge( EntityHandle eid,
                                 EntityHandle *hefid,
                                 int *helid);

    //! Gather half-edges to a queue of half-edges.
    /** Given a vertex vid, and a half-edge <he_fid,he_lid>, add another half-edge in the same face sharing the vertex
    */

    ErrorCode gather_halfedges(EntityHandle vid,
                                EntityHandle he_fid,
                                int he_lid,
                                int *qsize, int *count);

    //! Obtains another half-edge belonging to the same face as the input half-edge
    /** It uses the local maps to find another half-edge that is either incident or outgoing depending
     * on vid and input half-edge
     */

    ErrorCode another_halfedge( EntityHandle vid,
                                EntityHandle he_fid,
                                int he_lid,
                                EntityHandle *he2_fid,
                                int *he2_lid);

    ErrorCode mark_halfedges(EntityHandle vid, EntityHandle he_fid, int he_lid, Range &faces, std::vector<char> &markHEdgs, HFacet &bnd_hf);

    //! Collect and compare to find a matching half-edge with the given edge connectivity.
    /** Given edge connectivity, compare to an input list of half-edges to find a matching half-edge
     * and add a list of half-edges belonging to the one-ring neighborhood to a queue till it finds a match.
    */

    bool collect_and_compare(const EntityHandle vid, const EntityHandle* edg_vert,
                             int *qsize, int *count,
                             EntityHandle *he_fid,
                             int *he_lid);


    ErrorCode add_cells_of_single_component(EntityHandle vid, EntityHandle curcid, int curlid, std::multimap<EntityHandle, EntityHandle> &comps, HFacet &hf);
    bool find_cell_in_component(EntityHandle vid, EntityHandle cell, std::multimap<EntityHandle, EntityHandle> &comps);

    //! Given an edge, finds a matching local edge in an incident cell.
    /** Find a local edge with the same connectivity as the input edge, belonging to an incident cell.
     *
     * \param eid EntityHandle of the edge
     * \param cid Returns EntityHandle of the incident cell
     * \param leid Returns the local id of the edge corresponding to the input edge w.r.t the incident cell.
    */

    bool find_matching_implicit_edge_in_cell(EntityHandle eid,
                                              std::vector<EntityHandle> &cid,
                                              std::vector<int> &leid);

    //! Given a face, finds a matching local face in an incident cell.
    /** Find a local face with the same connectivity as the input face, belonging to an incident cell.
     *
     * \param fid EntityHandle of the face
     * \param cid Returns EntityHandle of the incident cell
     * \param lfid Returns the local id of the face corresponding to the input face w.r.t the incident cell.
    */

    bool find_matching_halfface(EntityHandle fid,
                                EntityHandle *cid,
                                int *leid);

    bool find_match_in_array(EntityHandle ent,
                             EntityHandle *ent_list,
                             int count,
                             bool get_index = false,
                             int *index = NULL);

  };

} // namespace moab 

#endif

