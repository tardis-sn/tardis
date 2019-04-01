#ifndef MOAB_BSP_TREE_POLY_HPP
#define MOAB_BSP_TREE_POLY_HPP

#include "moab/Types.hpp"
#include <vector>

namespace moab {

class CartVect;

/**\brief Convex polyhedron
 *
 * This class is used to represent the convex polyhedron that bounds
 * a node in a general plane-based BSP-tree.
 */
class BSPTreePoly 
{
  public:
    struct Vertex;
    struct VertexUse;
    struct Edge;
    struct EdgeUse;
    struct Face;
  private:
    Face* faceList;

    void set_vertex_marks( int value );

    BSPTreePoly( const BSPTreePoly& copy ); // not implemented
    BSPTreePoly& operator=( const BSPTreePoly& copy ); // not implemented

  public:
  
      /**\brief Initialize as a planar-faced hexahedron
       *\param hex_corners Corner coordinates for a hexahedron, in Exodus/Patran order
       */
    BSPTreePoly( const CartVect hex_corners[8] ) : faceList(0) { set(hex_corners); }
    BSPTreePoly( ) : faceList(0) { }
    ~BSPTreePoly() { clear(); }
    
      /**\brief Initialize as a planar-faced hexahedron
       *\param hex_corners Corner coordinates for a hexahedron, in Exodus/Patran order
       */
    ErrorCode set( const CartVect hex_corners[8] );
    void clear();
    
      /**\brief Get handles for faces */
    void get_faces( std::vector<const Face*>& face_list ) const;
      /**\brief Get corner coordinates for a face */
    void get_vertices( const Face* face, 
                       std::vector<CartVect>& vertices ) const;
    
      /** Intersect a plane with a polyhedron, retaining
       * the portion of the polyhedron below the plane.
       * This will fail if polyhedron is not convex. 
       */
    bool cut_polyhedron( const CartVect& plane_normal,
                         double plane_coeff );
    
      /** Test if a point is contained in the polyhedron.
       *
       *\NOTE algorithm assumes *convex* polyhedron.
       */
    bool is_point_contained( const CartVect& point ) const;
    
      //! Assumes planar faces
    double volume() const;
    
      // Check that data structure is consistent and represents
      // a closed polyhedron
    bool is_valid() const;
    
      /** For debugging, does nothing unless debug feature is enabled */
    static void reset_debug_ids();
};

} // namespace moab 

#endif
