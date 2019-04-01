/*
 * VerdictWrapper.hpp
 *
 *  Created on: Nov 18, 2014
 *      Author: iulian
 */

#ifndef SRC_VERDICT_MOAB_VERDICTWRAPPER_HPP_
#define SRC_VERDICT_MOAB_VERDICTWRAPPER_HPP_

#include <map>

namespace moab
{

class Interface;

enum QualityType {
  // order exactly from HexMetricVals
  MB_UNDEFINED_QUALITY = -1,
  MB_EDGE_RATIO = 0,        // 0  MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_MAX_EDGE_RATIO ,       // 1  MBHEX,                           MBQUAD
  MB_SKEW,                  // 2  MBHEX,                           MBQUAD
  MB_TAPER,                 // 3  MBHEX,                           MBQUAD
  MB_VOLUME,                // 4  MBHEX, MBTET, MBPRISM, MBKNIFE
  MB_STRETCH,               // 5  MBHEX,                           MBQUAD
  MB_DIAGONAL,              // 6  MBHEX,
  MB_DIMENSION,             // 7  MBHEX,
  MB_ODDY,                  // 8  MBHEX,                           MBQUAD
  MB_MED_ASPECT_FROBENIUS,  // 9  MBHEX,                           MBQUAD
  MB_MAX_ASPECT_FROBENIUS,  // 10 MBHEX, MBTET (aspect_frobenius)  MBQUAD,  MBTRI (aspect_frobenius)
  MB_CONDITION,             // 11 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_JACOBIAN,              // 12 MBHEX, MBTET,                    MBQUAD
  MB_SCALED_JACOBIAN,       // 13 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_SHEAR,                 // 14 MBHEX,                           MBQUAD,  MBTRI
  MB_SHAPE,                 // 15 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_RELATIVE_SIZE_SQUARED, // 16 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_SHAPE_AND_SIZE,        // 17 MBHEX, MBTET,                    MBQUAD
  MB_SHEAR_AND_SIZE,        // 18 MBHEX,                           MBQUAD
  MB_DISTORTION,            // 19 MBHEX, MBTET,                    MBQUAD
  // length for edge:
  MB_LENGTH,                // 20 only for MBEDGE
  // specific to tets
  MB_RADIUS_RATIO,          // 21        MBTET,                    MBQUAD,  MBTRI
  MB_ASPECT_BETA,           // 22        MBTET
  MB_ASPECT_RATIO,          // 23        MBTET,                    MBQUAD,  MBTRI
  MB_ASPECT_GAMMA,          // 24        MBTET
  MB_MINIMUM_ANGLE,         // 25        MBTET,                    MBQUAD,  MBTRI
  MB_COLLAPSE_RATIO,        // 26        MBTET
  // specific to quads
  MB_WARPAGE,               // 27                                  MBQUAD
  MB_AREA,                  // 28                                  MBQUAD,  MBTRI
  MB_MAXIMUM_ANGLE,         // 29                                  MBQUAD,  MBTRI
  MB_QUALITY_COUNT // used to size the arrays

};

inline
std::string QualityType_ToString(QualityType qtype)
{
  switch(qtype)
  {
    case MB_UNDEFINED_QUALITY:
      return "MB_UNDEFINED_QUALITY";
    case MB_EDGE_RATIO:
      return "MB_EDGE_RATIO";
    case MB_MAX_EDGE_RATIO:
      return "MB_MAX_EDGE_RATIO";
    case MB_SKEW:
      return "MB_SKEW";
    case MB_TAPER:
      return "MB_TAPER";
    case MB_VOLUME:
      return "MB_VOLUME";
    case MB_STRETCH:
      return "MB_STRETCH";
    case MB_DIAGONAL:
      return "MB_DIAGONAL";
    case MB_DIMENSION:
      return "MB_DIMENSION";
    case MB_ODDY:
      return "MB_ODDY";
    case MB_MED_ASPECT_FROBENIUS:
      return "MB_MED_ASPECT_FROBENIUS";
    case MB_MAX_ASPECT_FROBENIUS:
      return "MB_MAX_ASPECT_FROBENIUS";
    case MB_CONDITION:
      return "MB_CONDITION";
    case MB_JACOBIAN:
      return "MB_JACOBIAN";
    case MB_SCALED_JACOBIAN:
      return "MB_SCALED_JACOBIAN";
    case MB_SHEAR:
      return "MB_SHEAR";
    case MB_SHAPE:
      return "MB_SHAPE";
    case MB_RELATIVE_SIZE_SQUARED:
      return "MB_RELATIVE_SIZE_SQUARED";
    case MB_SHAPE_AND_SIZE:
      return "MB_SHAPE_AND_SIZE";
    case MB_SHEAR_AND_SIZE:
      return "MB_SHEAR_AND_SIZE";
    case MB_DISTORTION:
      return "MB_DISTORTION";
    case MB_LENGTH:
      return "MB_LENGTH";
    case MB_RADIUS_RATIO:
      return "MB_RADIUS_RATIO";
    case MB_ASPECT_BETA:
      return "MB_ASPECT_BETA";
    case MB_ASPECT_RATIO:
      return "MB_ASPECT_RATIO";
    case MB_ASPECT_GAMMA:
      return "MB_ASPECT_GAMMA";
    case MB_MINIMUM_ANGLE:
      return "MB_MINIMUM_ANGLE";
    case MB_COLLAPSE_RATIO:
      return "MB_COLLAPSE_RATIO";
    case MB_WARPAGE:
      return "MB_WARPAGE";
    case MB_AREA:
      return "MB_AREA";
    case MB_MAXIMUM_ANGLE:
      return "MB_MAXIMUM_ANGLE";
    default:
      return "MB_QUALITY_COUNT";
  }
}

class VerdictWrapper {
public:
  VerdictWrapper(Interface * mb);
  virtual ~VerdictWrapper();
  //! return a quality for an entity
  /** compute the quality for an element; the coordinates and number of nodes can be passed
   * if available
  \param  eh element entity handle.
  \param q quality requested
  \param quality output
  \param num_nodes optional, number of vertices
  \param coords options, interleaved coordinates
  return MB_SUCCESS
  Example: \code
  EntityHandle hex;
  double jac;
  rval = quality_measure(hex, MB_JACOBIAN, jac); \endcode
  */
  ErrorCode quality_measure(EntityHandle eh, QualityType q, double & quality,
      int num_nodes=0, EntityType etype=MBMAXTYPE, double *coords=NULL);
  //! return a quality name
    /** return quality name (convert an enum QualityType to a string)
    \param  q quality type
    return string
    Example: \code

    const char * name = quality_name(MB_JACOBIAN); \endcode
    */
  const char * quality_name (QualityType q);
  //! return a string with entity type name
  const char * entity_type_name(EntityType etype);
  //! return an int with total available qualities for type
  int num_qualities(EntityType etype);
  //! return true if quality possible
  int possible_quality(EntityType et, QualityType q);
  // relative size needs a base size, that is set at global level, one for each major type (hex, tet, quad, tri)
  ErrorCode set_size(double size);
  //! return all qualities for an element
  /** compute all qualities for an element
    \param  eh element entity handle.
    \param qs list of QualityType
    \param qualities list of qualities
    return MB_SUCCESS
    Example: \code
    EntityHandle hex;
    std::vector<QualityType> qs;
    std::vector<double> qualities;
    all_quality_measures(hex, qs, qualities); \endcode
    */
  ErrorCode all_quality_measures(EntityHandle eh, std::map<QualityType, double> & qualities);
private:
  Interface * mbImpl;

};
} // namespace moab
#endif /* SRC_VERDICT_MOAB_VERDICTWRAPPER_HPP_ */
