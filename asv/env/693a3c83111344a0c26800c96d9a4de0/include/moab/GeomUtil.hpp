/*
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

/**\file Geometry.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-27
 */

#ifndef MB_GEOM_UTIL_HPP
#define MB_GEOM_UTIL_HPP

#include "moab/CartVect.hpp"
#include <cmath>

namespace moab {

/** \class GeomUtil
 * \brief Functions for computational geometry on triangles, rays, and boxes
 */
namespace GeomUtil {

/** Given a line segment and an axis-aligned box, 
 *  return the sub-segment of the line segment that 
 *  itersects the box.
 *  
 *  Can be used to intersect ray with box by passing seg_end
 *  as HUGE_VAL or std::numeric_limits<double>::maximum().
 *
 *\param box_min   Minimum corner of axis-aligned box
 *\param box_max   Maximum corner of axis-aligned box
 *\param seg_pt    A point in the line containing the segement
 *\param seg_unit_dir A unit vector in the direction of the line 
 *                 containing the semgent.
 *\param seg_start The distance from seg_pt in the direction of
 *                 seg_unit_dir at which the segment begins.
 *                 As input, the start of the original segment, as output, the
 *                 start of the sub-segment intersecting the box.
 *                 Note:  seg_start must be less than seg_end
 *\param seg_end   The distance from seg_pt in the direction of 
 *                 seg_unit_dir at which the segment ends.
 *                 As input, the end of the original segment, as output, the
 *                 end of the sub-segment intersecting the box.
 *                 Note:  seg_start must be less than seg_end
 *\return true if line semgent intersects box, false otherwise.
 */
bool segment_box_intersect( CartVect box_min,
                            CartVect box_max,
                            const CartVect& seg_pt,
                            const CartVect& seg_unit_dir,
                            double& seg_start, double& seg_end );

/**\brief Test for intersection between a ray and a triangle.
 *\param ray_point  The start point of the ray.
 *\param ray_unit_direciton  The direction of the ray. Must be a unit vector.
 *\param tolerance  Absolute distance tolerance for point equality
 *\param t_out Output: The distance along the ray from ray_point in the
 *                  direction of ray_unit_direction at which the ray
 *                  itersected the triangle.
 *\param ray_length Optional:  If non-null, a pointer a maximum length
 *                  for the ray, converting this function to a segment-tri-
 *                  intersect test.
 *\return true if intersection, false otherwise.
 */
bool ray_tri_intersect( const CartVect vertices[3],
                        const CartVect& ray_point,
                        const CartVect& ray_unit_direction,
                        double tolerance,
                        double& t_out,
                        const double* ray_length = 0 );


/**\brief Plücker test for intersection between a ray and a triangle.
 *\param vertices            Nodes of the triangle.
 *\param ray_point           The start point of the ray.
 *\param ray_unit_direction  The direction of the ray. Must be a unit vector.
 *\param tolerance           Absolute distance tolerance for point equality
 *\param t_out               Output: The distance along the ray from ray_point in the
 *                           direction of ray_unit_direction at which the ray
 *                           intersected the triangle.
 *\param nonneg_ray_length   Optional: If non-null, a maximum length for the ray, 
 *                           converting this function to a segment-tri-intersect 
 *                           test.
 *\param neg_ray_length      Optional: If non-null, a maximum length for the ray 
 *                           behind the origin, converting this function to a 
 *                           segment-tri-intersect test.
 *\param orientation         Optional: Reject intersections without the specified
 *                           orientation of ray with respect to triangle normal
 *                           vector. Indicate desired orientation by passing
 *                           1 (forward), -1 (reverse), or 0 (no preference).
 *\param int_type            Optional Output: The type of intersection; used to
 *                           identify edge/node intersections.
 *\return true if intersection, false otherwise.
 */
enum intersection_type {NONE, INTERIOR, NODE0, NODE1, NODE2, EDGE0, EDGE1, EDGE2};
bool plucker_ray_tri_intersect( const CartVect vertices[3],
                                const CartVect& ray_point,
                                const CartVect& ray_unit_direction,
                                double tolerance, /* unused */
                                double& dist_out,
                                const double* nonneg_ray_length = 0,
                                const double* neg_ray_length = 0,
                                const int* orientation = 0,
                                intersection_type* int_type = 0);

    //! Find range of overlap between ray and axis-aligned box.
    //!
    //!\param box_min   Box corner with minimum coordinate values
    //!\param box_max   Box corner with minimum coordinate values
    //!\param ray_pt    Coordinates of start point of ray
    //!\param ray_dir   Directionion vector for ray such that
    //!                 the ray is defined by r(t) = ray_point + t * ray_vect
    //!                 for t > 0.
    //!\param t_enter   Output: if return value is true, this value
    //!                 is the parameter location along the ray at which
    //!                 the ray entered the leaf.  If return value is false,
    //!                 then this value is undefined.
    //!\param t_exit    Output: if return value is true, this value
    //!                 is the parameter location along the ray at which
    //!                 the ray exited the leaf.  If return value is false,
    //!                 then this value is undefined.
    //!\return true if ray intersects leaf, false otherwise.
bool ray_box_intersect( const CartVect& box_min,
                        const CartVect& box_max,
                        const CartVect& ray_pt,
                        const CartVect& ray_dir,
                        double& t_enter, double& t_exit );

/**\brief Test if plane intersects axis-aligned box
 *
 * Test for intersection between an unbounded plane and
 * an axis-aligned box.
 *\param plane_normal Vector in plane normal direction (need *not*
 *                    be a unit vector).  The N in 
 *                    the plane equation: N . X + D = 0
 *\param plane_coeff  The scalar 'D' term in the plane equation:
 *                    N . X + D = 0
 *\param box_min_corner The smallest coordinates of the box along each
 *                    axis.  The corner of the box for which all three
 *                    coordinate values are smaller than those of any
 *                    other corner.  The X, Y, Z values for the planes
 *                    normal to those axes and bounding the box on the
 *                    -X, -Y, and -Z sides respectively.
 *\param box_max_corner The largest coordinates of the box along each
 *                    axis.  The corner of the box for which all three
 *                    coordinate values are larger than those of any
 *                    other corner.  The X, Y, Z values for the planes
 *                    normal to those axes and bounding the box on the
 *                    +X, +Y, and +Z sides respectively.
 *\return true if overlap, false otherwise.
 */
bool box_plane_overlap( const CartVect& plane_normal, 
                        double            plane_coeff,
                        CartVect        box_min_corner, 
                        CartVect        box_max_corner );

/**\brief Test if triangle intersects axis-aligned box
 *
 * Test if a triangle intersects an axis-aligned box.
 *\param triangle_corners  The corners of the triangle.
 *\param box_min_corner The smallest coordinates of the box along each
 *                    axis.  The corner of the box for which all three
 *                    coordinate values are smaller than those of any
 *                    other corner.  The X, Y, Z values for the planes
 *                    normal to those axes and bounding the box on the
 *                    -X, -Y, and -Z sides respectively.
 *\param box_max_corner The largest coordinates of the box along each
 *                    axis.  The corner of the box for which all three
 *                    coordinate values are larger than those of any
 *                    other corner.  The X, Y, Z values for the planes
 *                    normal to those axes and bounding the box on the
 *                    +X, +Y, and +Z sides respectively.
 *\param tolerance    The tolerance used in the intersection test.  The box
 *                    size is increased by this amount before the intersection
 *                    test.
 *\return true if overlap, false otherwise.
 */
bool box_tri_overlap( const CartVect  triangle_corners[3],
                      const CartVect& box_min_corner,
                      const CartVect& box_max_corner,
                      double            tolerance );

/**\brief Test if triangle intersects axis-aligned box
 *
 * Test if a triangle intersects an axis-aligned box.
 *\param triangle_corners  The corners of the triangle.
 *\param box_center   The center of the box.
 *\param box_hanf_dims The distance along each axis, respectively, from the
 *                    box_center to the boundary of the box.
 *\return true if overlap, false otherwise.
 */
bool box_tri_overlap( const CartVect  triangle_corners[3],
                      const CartVect& box_center,
                      const CartVect& box_half_dims );

bool box_point_overlap( const CartVect& box_min_corner,
                        const CartVect& box_max_corner,
                        const CartVect& point,
                        double tolerance );

/**\brief Test if the specified element intersects an axis-aligned box.
 *
 * Test if element intersects axis-aligned box.  Use element-specific
 * optimization if available, otherwise call box_general_elem_overlap.
 *
 *\param elem_corners The coordinates of the element vertices
 *\param elem_type    The toplogy of the element.
 *\param box_center   The center of the axis-aligned box
 *\param box_half_dims Half of the width of the box in each axial
 *                     direction.
 */
bool box_elem_overlap( const CartVect *elem_corners,
                       EntityType elem_type,
                       const CartVect& box_center,
                       const CartVect& box_half_dims ); 

/**\brief Test if the specified element intersects an axis-aligned box.
 *
 * Uses MBCN and separating axis theorem for general algorithm that
 * works for all fixed-size elements (not poly*).
 *
 *\param elem_corners The coordinates of the element vertices
 *\param elem_type    The toplogy of the element.
 *\param box_center   The center of the axis-aligned box
 *\param box_half_dims Half of the width of the box in each axial
 *                     direction.
 */
bool box_linear_elem_overlap( const CartVect *elem_corners,
                              EntityType elem_type,
                              const CartVect& box_center,
                              const CartVect& box_half_dims ); 

/**\brief Test if the specified element intersects an axis-aligned box.
 *
 * Uses MBCN and separating axis theorem for general algorithm that
 * works for all fixed-size elements (not poly*).  Box and element
 * vertices must be translated such that box center is at origin.
 *
 *\param elem_corners The coordinates of the element vertices, in 
 *                    local coordinate system of box.
 *\param elem_type    The toplogy of the element.
 *\param box_half_dims Half of the width of the box in each axial
 *                     direction.
 */
bool box_linear_elem_overlap( const CartVect *elem_corners,
                              EntityType elem_type,
                              const CartVect& box_half_dims ); 

void closest_location_on_box( const CartVect& box_min_corner,
                              const CartVect& box_max_corner,
                              const CartVect& point,
                              CartVect& closest );

/**\brief find closest location on triangle
 *
 * Find closest location on linear triangle.
 *\param location  Input position to evaluate from
 *\param vertices  Array of three corner vertex coordinates.
 *\param closest_out Result position 
 */
void closest_location_on_tri( const CartVect& location,
                              const CartVect* vertices,
                              CartVect& closest_out );

/**\brief find closest location on polygon
 *
 * Find closest location on polygon
 *\param location  Input position to evaluate from
 *\param vertices  Array of corner vertex coordinates.
 *\param num_vertices Length of 'vertices' array.
 *\param closest_out Result position 
 */
void closest_location_on_polygon( const CartVect& location,
                                  const CartVect* vertices,
                                  int num_vertices,
                                  CartVect& closest_out );

/**\brief find closest topological location on triangle
 *
 * Find closest location on linear triangle.
 *\param location  Input position to evaluate from
 *\param vertices  Array of three corner vertex coordinates.
 *\param tolerance Tolerance to use when comparing to corners and edges
 *\param closest_out Result position 
 *\param closest_topo Closest topological entity
 *                     0-2 : vertex index
 *                     3-5 : edge beginning at closest_topo - 3
 *                       6 : triangle interior
 */
void closest_location_on_tri( const CartVect& location,
                              const CartVect* vertices,
                              double tolerance,
                              CartVect& closest_out,
                              int& closest_topo );

// Finds whether or not a box defined by the center and the half
// width intersects a trilinear hex defined by its eight vertices.
bool box_hex_overlap( const CartVect hexv[8],
                      const CartVect& box_center,
                      const CartVect& box_dims);

// Finds whether or not a box defined by the center and the half
// width intersects a linear tetrahedron defined by its four vertices.
bool box_tet_overlap( const CartVect tet_corners[4],
                      const CartVect& box_center,
                      const CartVect& box_dims);

// tests if 2 boxes overlap within a tolerance
// assume that data is valid, box_min1 << box_max1, and box_min2<< box_max2
// they are overlapping if they are overlapping in all directions
// for example in x direction:
//   box_max2[0] + tolerance < box_min1[0] -- false
bool boxes_overlap( const CartVect & box_min1, const CartVect & box_max1,
    const CartVect & box_min2, const CartVect & box_max2, double tolerance);

// see if boxes formed by 2 lists of "CartVect"s overlap
bool bounding_boxes_overlap (const CartVect * list1, int num1, const CartVect * list2, int num2,
      double tolerance);
//
// point_in_trilinear_hex
// Tests if a point in xyz space is within a hex element defined with
// its eight vertex points forming a trilinear basis function.  Computes
// the natural coordinates with respect to the hex of the xyz point 
// and checks if each are between +/-1.  If anyone is outside the range
// the function returns false, otherwise it returns true.
//
bool point_in_trilinear_hex(const CartVect *hex, 
                            const CartVect& xyz,
                            double etol);

bool nat_coords_trilinear_hex( const CartVect* corner_coords,
                               const CartVect& x,
                               CartVect& xi,
                               double tol );
} // namespace GeomUtil

} // namespace moab

#endif
