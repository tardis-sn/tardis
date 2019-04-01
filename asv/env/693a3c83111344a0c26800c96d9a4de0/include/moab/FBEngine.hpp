#ifndef FBENGINE_HPP_
#define FBENGINE_HPP_
#include <stdlib.h>

#include <vector>
#include <map>

#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"

namespace moab {
class GeomTopoTool;

// some forward declarations
class SmoothFace;
class SmoothCurve;

/*
 *  Facet Based engine class for mesh-based geometry
 */
class FBEngine {
public:
  FBEngine(Interface *impl, GeomTopoTool* geomTopoTool = NULL,
      const bool smooth = false);

  ~FBEngine();

  ErrorCode Init();

  ErrorCode getRootSet(EntityHandle * root_set);

  ErrorCode getNumEntSets(EntityHandle set, int num_hops, int * all_sets);

  ErrorCode createEntSet(int isList, EntityHandle * pSet);

  ErrorCode addEntSet(EntityHandle entity_set_to_add,
      EntityHandle entity_set_handle);

  ErrorCode getEntities(EntityHandle root_set, int ent_type, Range & gentities);

  ErrorCode addEntArrToSet(Range entities, EntityHandle set);

  ErrorCode getNumOfType(EntityHandle set, int ent_type, int * pNum);

  ErrorCode getEntType(EntityHandle gent, int * type);

  ErrorCode getEntBoundBox(EntityHandle this_gent, double * x0, double * y0,
      double * z0, double * x1, double *y1, double * z1);
  ErrorCode getEntClosestPt(EntityHandle this_gent, double x, double y,
      double z, double * x1, double * y1, double *y3);

  ErrorCode getVtxCoord(EntityHandle this_gent, double * x0, double * y0, double * z0);

  ErrorCode gsubtract(EntityHandle entity_set_1, EntityHandle entity_set_2,
      EntityHandle result_entity_set);

  ErrorCode getEntNrmlXYZ( EntityHandle entity_handle, double x, double y, double z,
      double* nrml_i, double* nrml_j, double* nrml_k);

  ErrorCode getPntRayIntsct( double x, double y, double z,
      double dir_x, double dir_y, double dir_z,
      std::vector<EntityHandle> &intersect_entity_handles,
      /* int storage_order,*/
      std::vector<double> & intersect_coords,
      std::vector<double> & param_coords);

  // some new methods, that are needed

  ErrorCode createTag( const char* tag_name,
                    int tag_num_type_values,
                    int tag_type,
                    Tag & tag_handle_out );

  Interface * moab_instance () { return _mbImpl; }

  ErrorCode getArrData(const moab::EntityHandle* entity_handles,
      int entity_handles_size,
      Tag tag_handle,
      void* tag_values_out);

  ErrorCode setArrData(const EntityHandle* entity_handles,
        int entity_handles_size,
        Tag tag_handle,
        const void* tag_values);

  ErrorCode getEntAdj(EntityHandle handle,
        int type_requested, Range & adjEnts  );

  ErrorCode getEgFcSense(EntityHandle mbedge, EntityHandle mbface, int & sense );

  ErrorCode measure(const EntityHandle * moab_entities, int entities_size,
      double * measures);

  // to do
  ErrorCode getEntNrmlSense( EntityHandle face, EntityHandle region,
      int& sense );

  ErrorCode getEgEvalXYZ( EntityHandle edge,
                                 double x, double y, double z,
                                 double& on_x, double& on_y, double& on_z,
                                 double& tngt_i, double& tngt_j, double& tngt_k,
                                 double& cvtr_i, double& cvtr_j, double& cvtr_k );
  ErrorCode getFcEvalXYZ( EntityHandle face,
                                 double x, double y, double z,
                                 double& on_x, double& on_y, double& on_z,
                                 double& nrml_i, double& nrml_j, double& nrml_k,
                                 double& cvtr1_i, double& cvtr1_j, double& cvtr1_k,
                                 double& cvtr2_i, double& cvtr2_j, double& cvtr2_k );

  ErrorCode getEgVtxSense( EntityHandle edge, EntityHandle vtx1, EntityHandle vtx2, int& sense );

  ErrorCode getEntURange( EntityHandle edge,
                                 double& u_min, double& u_max );

  ErrorCode getEntUtoXYZ( EntityHandle edge, double u,
                                 double& x, double& y, double& z );

  ErrorCode getEntTgntU( EntityHandle edge,
                                    double u,
                                    double& i, double& j, double& k );

  ErrorCode isEntAdj( EntityHandle entity1, EntityHandle entity2,
      bool& adjacent_out );

  ErrorCode split_surface_with_direction(EntityHandle face, std::vector<double> & xyz, double * direction,
      int closed, double min_dot, EntityHandle & oNewFace );
  // these new points will be on edges or triangles, if in interior of triangles
  ErrorCode split_surface(EntityHandle face,
      std::vector<EntityHandle> & chainedEdges,
      std::vector<EntityHandle> & splittingNodes, EntityHandle & newFace);

  ErrorCode split_edge_at_point(EntityHandle edge, CartVect & point, EntityHandle & new_edge);

  ErrorCode split_edge_at_mesh_node(EntityHandle edge, EntityHandle node, EntityHandle & new_edge);

  ErrorCode split_bedge_at_new_mesh_node(EntityHandle b_edge, EntityHandle atNode, EntityHandle brokenEdge,
      EntityHandle & new_edge);
  // helper for cleaning the stuff
  // will be called if the topology is modified
  void clean();

  void delete_smooth_tags();

  // access to the geom topo tool
  // be careful what you do with it
  GeomTopoTool * get_gtt() { return this->_my_geomTopoTool ; }

  ErrorCode create_volume_with_direction(EntityHandle newFace1, EntityHandle newFace2, double * direction,
      EntityHandle & volume);

  // get nodes from edge in order
  ErrorCode get_nodes_from_edge(EntityHandle gedge, std::vector<EntityHandle> & nodes);

  ErrorCode  weave_lateral_face_from_edges(EntityHandle bEdge, EntityHandle tEdge,  double * direction,
      EntityHandle & newLatFace);

  // chain "chain"-able edges
  // this could be useful if there are too many points / edges in the splitting
  // polyline
  // 2 edges are "chain"-able if
  /*
   *  1. they have the same adjacent faces
   *  2. their orientation is such as the end of one is the start of the other
   *  3. at meeting point, their tangents make an angle smaller then something
   *   (cos(angle) > cos (max_angle) = min_dot)
   */
  ErrorCode chain_edges(double min_dot);

  // 2 edges will be chained, along with modification of the topology
  ErrorCode chain_two_edges(EntityHandle edge, EntityHandle next_edge);

  ErrorCode get_vert_edges(EntityHandle edge, EntityHandle & v1, EntityHandle & v2);

  void set_smooth() { _smooth = true;}
private:

  ErrorCode initializeSmoothing();

  ErrorCode getAdjacentEntities(const EntityHandle from, const int to_dim,
      Range &adj_ents);

  ErrorCode compute_intersection_points(EntityHandle & face,
      EntityHandle from, EntityHandle to, CartVect & Dir, std::vector<CartVect> & points,
      std::vector<EntityHandle> & entities, std::vector<EntityHandle> & triangles);

  ErrorCode BreakTriangle(EntityHandle tri, EntityHandle e1, EntityHandle e3, EntityHandle n1,
      EntityHandle n2, EntityHandle n3);// nodesAlongPolyline are on entities!

  ErrorCode BreakTriangle2(EntityHandle tri, EntityHandle e1, EntityHandle e2, EntityHandle n1,
        EntityHandle n2);// nodesAlongPolyline are on entities!

  void print_debug_triangle(EntityHandle triangle);

  ErrorCode  create_new_gedge(std::vector<EntityHandle> &nodesAlongPolyline,
      EntityHandle & new_geo_edge);

  // used for splitting surfaces
  ErrorCode separate (EntityHandle face,
      std::vector<EntityHandle> & chainedEdges, Range & first,
      Range & second);

  ErrorCode smooth_new_intx_points(EntityHandle face,
      std::vector<EntityHandle> & chainedEdges);

  // having a node, split boundary along that node
  ErrorCode  split_boundary(EntityHandle face, EntityHandle atNode);

  // see if the node is already part of a vertex set, do not create another one
  bool  find_vertex_set_for_node(EntityHandle iNode, EntityHandle & oVertexSet);

  // if the splitting edge is not a loop, the original boundary edges will belong to
  // either original face, or new face
  // only the new geo edge (splitting) will be part of both, with the
  // orientation already decided
  //
  ErrorCode redistribute_boundary_edges_to_faces(EntityHandle face, EntityHandle newFace,
      std::vector<EntityHandle> & chainedEdges);

  // used as a way to isolate the faces after splitting
  ErrorCode set_neumann_tags(EntityHandle face, EntityHandle newFace);

  // split the quads if needed; it will create a new gtt, which will
  // contain triangles instead of quads
  ErrorCode split_quads();

  ErrorCode boundary_nodes_on_face(EntityHandle face, std::vector<EntityHandle> & boundary_nodes);

  ErrorCode boundary_mesh_edges_on_face(EntityHandle face, Range & boundary_mesh_edges);

  // used for splitting an edge
  ErrorCode split_internal_edge(EntityHandle & edge, EntityHandle & newVertex);
  // triangle split
  ErrorCode divide_triangle(EntityHandle triangle, EntityHandle & newVertex);
  Interface * _mbImpl;

  // this will be used during volume creation
  ErrorCode set_default_neumann_tags();

  ErrorCode chain_able_edge(EntityHandle edge, double min_dot, EntityHandle & next_edge, bool & chainable);

  GeomTopoTool* _my_geomTopoTool;
  bool _t_created;
  bool _smooth;
  bool _initialized;
  // these are initial ranges, that should not change during geometry gimmicks
  // those that are added are changing gtt ranges, but they are not yet "smoothed"
  // when new geometry is created, these ranges are not yet updated
  Range _my_gsets[5];
  // these are used only for smooth evaluations
  // these smooth faces and edges will be initialized after reading the file
  // the maps keep the link between EH in moab (geom sets) and
  //   their corresponding smooth counterparts
  std::map<EntityHandle, SmoothFace*> _faces;
  std::map<EntityHandle, SmoothCurve*> _edges;
  SmoothFace ** _smthFace;
  SmoothCurve ** _smthCurve;

  Range _piercedTriangles; // triangles to delete
  Range _newTriangles;
  Range _piercedEdges;// edges to delete
  std::map<EntityHandle, EntityHandle> _brokenEdges; // this is a map between splitting nodes and edges that
  // are broken on the boundary, by the polyline; the edges will be deleted in the end
  // new edges?

};

} // namespace moab
#endif /* FBENGINE_HPP_ */
