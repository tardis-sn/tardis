#ifndef MOABMC_HPP
#define MOABMC_HPP

#include "moab/Interface.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <map>
#include <string>
#include <assert.h>

#include "moab/OrientedBoxTreeTool.hpp"

class RefEntity;

struct DagmcVolData {
  int mat_id;
  double density, importance;
  std::string comp_name;
};


namespace moab {

class CartVect;

#define DAGMC_VERSION 1.0
#define DAGMC_VERSION_STRING "1.0"
#define DAGMC_INTERFACE_REVISION "$Rev$"

/**\brief
 *
 * In sectction 1, the public interface you will find all the functions needed
 * for problem setup. For the typical MC code, the order of function calls
 * required to fully populate DAGMC ready to run are
 *
 *    1) DAG->load_file();
 *    2) DAG->init_OBBTree();
 *
 * Modifications were made to init_OBBTree which allows the functions of
 * init_OBBTree to be called without having used init_OBBTree. For example
 * if you would like access to be able to call DAG->point_in_volume() but without
 * having an implicit compliment you need only call
 *
 *   1) DAG->load_file();
 *   2) DAG->setup_obb();
 *
 * Similarly, if you need access to problem indices only then, one may call
 * load_file followed by setup_indices.
 *
 *   1) DAG->load_file();
 *   2) DAG->setup_indices();
 */

class DagMC
{
public:
  static DagMC *instance(Interface *mb_impl = NULL);
  static void destroy();


  ~DagMC() {}

  /** Return the version of this library */
  static float version(std::string *version_string = NULL);
  /** Get subversion revision of this file (DagMC.hpp) */
  static unsigned int interface_revision();


  /* SECTION I: Geometry Initialization */

  /**\brief Load a geometry description regardless of format
   *
   * This method will load the geometry file with name cfile.
   * In case this is a solid model geometry file, it will pass
   * the facet_tolerance option as guidance for the faceting engine.
   *\param cfile the file name to be loaded
   *\param facet_tolerance the faceting tolerance guidance for the faceting engine
   *\return - MB_SUCCESS if file loads correctly
   *        - other MB ErrorCodes returned from MOAB
   */
  ErrorCode load_file(const char* cfile,
                      const double facet_tolerance = 0);

  /*\brief Use pre-loaded geometry set
   *
   * Works like load_file, but using data that has been externally
   * loaded into DagMC's MOAB instance.
   * Only one of the two functions should be called.
   *
   * TODO: this function should accept a parameter, being the
   * entity set to use for DagMC's data.  Currently DagMC always
   * assumes that all the contents of its MOAB instance belong to it.
   */
  ErrorCode load_existing_contents();

  /**\brief initialize the OBB tree structure for ray firing acceleration
   *
   * This method generates an OBB tree from the faceted representation of
   * the geometry.  It also calls internal methods to generate the implicit
   * complement and to build the cross-referencing indices.
   */
  ErrorCode init_OBBTree();

  /**\brief sets up storage for the implicit complimennt
   *
   * This method generates the implicit compliment storage, in normal situations
   * this will be called just before init_OBBTree();
   */
  ErrorCode setup_impl_compl();

  /**\brief sets up ranges of the volume and surface entity sets
   *
   * Helper function for setup_obbs & setup_indices, sets ranges containing
   * all volumes and surfaces
   */
  ErrorCode setup_geometry(Range &surfs, Range &vols);

  /**\brief thin wrapper around build_obbs()
   *
   * Very thin wrapper around build_obbs, allows the setup_OBBTree() function to be
   * made into atomic chunks
   */
  ErrorCode setup_obbs();

  /**\brief thin wrapper around build_indices()
   *
   * Very thin wrapper around build_indices, allows the setup_OBBTree() function to be
   * made into atomic chunks
   */
  ErrorCode setup_indices();


private:
  /** loading code shared by load_file and load_existing_contents */
  ErrorCode finish_loading();

  /** test for existing OBB Tree */
  bool have_obb_tree();

  /** test for exisiting implicit compliment */
  bool have_impl_compl();

  /** test for pre-existing implicit complement definition, or return a new one */
  ErrorCode get_impl_compl();

  /** build obb structure for each surface and volume */
  ErrorCode build_obbs(Range &surfs, Range &vols);

  /** build obb structure for the implicit complement */
  ErrorCode build_obb_impl_compl(Range &surfs);


  /* SECTION II: Fundamental Geometry Operations/Queries */
public:

  /**\brief State object used in calls to ray_fire()
   *
   * Storage for the "history" of a ray.  This represents the surface facets
   * that the ray is known to have crossed, which cannot be crossed again
   * as long as the ray does not change direction.  It is intended to be used
   * with a series of consecutive calls to ray_fire(), in which a ray passes
   * over potentially many surfaces.
   */
  class RayHistory {

  public:
    /**
     * Clear this entire history-- logically equivalent to creating a new history,
     * but probably more efficient.
     */
    void reset();

    /**
     * Clear the history up to the most recent intersection.  This should be
     * called when a ray changes direction at the site of a surface crossing,
     * a situation that most commonly occurs at a reflecting boundary.
     */
    void reset_to_last_intersection();

    /**
     * Remove the most recent intersection.  This allows a subsequent call
     * along the same ray to return the same intersection.
     */
    void rollback_last_intersection();

    /**
     * @return the number of surface crossings currently represented by this ray history
     */
    int size() const { return prev_facets.size(); }

  private:
    std::vector<EntityHandle> prev_facets;

    friend class DagMC;

  };

  /**\brief find the next surface crossing from a given point in a given direction
   *
   * This is the primary method of DagMC, enabling ray tracing through a geometry.
   * Given a volume and a ray, it computes the surface ID and distance to the
   * nearest intersection on that volume.  The caller can compute the location of
   * the intersection by adding the distance to the ray.
   *
   * When a series of calls to this function are made along the same ray (e.g. for
   * the purpose of tracking a ray through several volumes), the optional history
   * argument should be given.  The history prevents previously intersected facets
   * from being intersected again.  A single history should be used as long as a
   * ray is proceeding forward without changing direction.  This situation is
   * sometimes referred to as "streaming."
   *
   * If a ray changes direction at an intersection site, the caller should call
   * reset_to_last_intersection() on the history object before the next ray fire.
   *
   * @param volume The volume to fire the ray at.
   * @param ray_start An array of x,y,z coordinates from which to start the ray.
   * @param ray_dir An array of x,y,z coordinates indicating the direction of the ray.
   *                Must be of unit length.
   * @param next_surf Output parameter indicating the next surface intersected by the ray.
   *                If no intersection is found, will be set to 0.
   * @param next_surf_dist Output parameter indicating distance to next_surf.  If next_surf is
   *                0, this value is undefined and should not be used.
   * @param history Optional RayHistory object.  If provided, the facets in the history are
   *                assumed to not intersect with the given ray.  The facet intersected
   *                by this query will also be added to the history.
   * @param dist_limit Optional distance limit.  If provided and > 0, no intersections at a
   *                distance further than this value will be returned.
   * @param ray_orientation Optional ray orientation. If provided determines intersections
   *                along the normal provided, e.g. if -1 allows intersections back along the
   *                the ray direction, Default is 1, i.e. exit intersections
   * @param stats Optional TrvStats object used to measure performance of underlying OBB
   *              ray-firing query.  See OrientedBoxTreeTool.hpp for details.
   *
   */
  ErrorCode ray_fire(const EntityHandle volume,
                     const double ray_start[3], const double ray_dir[3],
                     EntityHandle& next_surf, double& next_surf_dist,
                     RayHistory* history = NULL, double dist_limit = 0,
		     int ray_orientation = 1,
                     OrientedBoxTreeTool::TrvStats* stats = NULL );

  /**\brief Test if a point is inside or outside a volume
   *
   * This method finds the point on the boundary of the volume that is nearest
   * the test point (x,y,z).  If that point is "close" to a surface, a boundary test
   * is performed based on the normal of the surface at that point and the
   * optional ray direction (u,v,w).
   * @param volume The volume to test
   * @param xyz The location to test for volume containment
   * @param result Set to 0 if xyz it outside volume, 1 if inside, and -1 if on boundary.
   * @param Optional direction to use for underlying ray fire query.  Used to ensure
   *        consistent results when a ray direction is known.  If NULL or {0,0,0} is
   *        given, a random direction will be used.
   * @param history Optional RayHistory object to pass to underlying ray fire query.
   *        The history is not modified by this call.
   */
  ErrorCode point_in_volume(const EntityHandle volume,
                            const double xyz[3],
                            int& result,
                            const double* uvw = NULL,
                            const RayHistory* history = NULL );

  /**\brief Robust test if a point is inside or outside a volume using unit sphere area method
   *
   * This test may be more robust that the standard point_in_volume, but is much slower.
   * It does not detect 'on boundary' situations as point_in_volume does.
   * @param volume The volume to test
   * @param xyz The location to test for volume containment
   * @param result Set to 0 if xyz it outside volume, 1 if inside.
   */
  ErrorCode point_in_volume_slow( const EntityHandle volume, const double xyz[3], int& result );


  /** \brief Given a ray starting at a surface of a volume, check whether the ray enters or exits the volume
   *
   * This function is most useful for rays that change directions at a surface crossing.
   * It can be used to check whether a direction change redirects the ray back into the originating
   * volume.
   *
   * @param volume The volume to test
   * @param surface A surface on volume
   * @param xyz A point location on surface
   * @param uvw A (unit) direction vector
   * @param result Set to 1 if ray is entering volume, or 0 if it is leaving
   * @param history Optional ray history object from a previous call to ray_fire.  If present and non-empty,
   *        the history is used to look up the surface facet at which the ray begins.  Absent a
   *        history, the facet nearest to xyz will be looked up.  The history should always be
   *        provided if available, as it avoids the computational expense of a nearest-facet query.
   */
  ErrorCode test_volume_boundary( const EntityHandle volume, const EntityHandle surface,
                                  const double xyz[3], const double uvw[3], int& result,
                                  const RayHistory* history = NULL );

  /**\brief Find the distance to the point on the boundary of the volume closest to the test point
   *
   * @param volume Volume to query
   * @param point Coordinates of test point
   * @param result Set to the minimum distance from point to a surface in volume
   */
  ErrorCode closest_to_location( EntityHandle volume, const double point[3], double& result);

  /** Calculate the volume contained in a 'volume' */
  ErrorCode measure_volume( EntityHandle volume, double& result );

  /** Calculate sum of area of triangles */
  ErrorCode measure_area( EntityHandle surface, double& result );

  /** Get the sense of surfaces wrt a volume.  Sense values are:
   *  {-1 -> reversed, 0 -> both, 1 -> forward}
   */
  ErrorCode surface_sense( EntityHandle volume,
                             int num_surfaces,
                             const EntityHandle* surfaces,
                             int* senses_out );

  /** Get the sense of a single surface wrt a volume.  Sense values are:
   *  {-1 -> reversed, 0 -> both, 1 -> forward}
   */
  ErrorCode surface_sense( EntityHandle volume, EntityHandle surface, int& sense_out );

  /** Get the normal to a given surface at the point on the surface closest to a given point
   *
   * @param surf Surface on which to get normal
   * @param xyz Point on surf
   * @param angle Set to coordinates of surface normal nearest xyz
   * @param history Optional ray history from a previous call to ray_fire().
   *        If present and non-empty, return the normal
   *        of the most recently intersected facet, ignoring xyz.
   */
  ErrorCode get_angle(EntityHandle surf, const double xyz[3], double angle[3],
                      const RayHistory* history = NULL );

  /** Get the volume on the other side of a surface
   *
   * @param A surface to query
   * @param old_volume A volume on one side of surface
   * @param new_volume Output parameter for volume on the other side of surface
   * @return MB_SUCCESS if new_volume was set successfully, error if not.
   */
  ErrorCode next_vol( EntityHandle surface, EntityHandle old_volume,
                      EntityHandle& new_volume );

private:
  /**\brief pass the ray_intersection test to the solid modeling engine
   *
   * The user has the options to specify that ray tracing should ultimately occur on the
   * true CAD model rather than just on the faceted representation.  This is called from
   * within ray_fire if the user has selected that option
   */
  ErrorCode CAD_ray_intersect(const double *point,
                                const double *dir,
                                const double huge_val,
                                std::vector<double> &distances,
                                std::vector<EntityHandle> &surfaces,
                                double &len);

  /**\brief determine the point membership when the point is effectively on the boundary
   *
   * Called by point_in_volume when the point is with tolerance of the boundary. Compares the
   * ray direction with the surface normal to determine a volume membership.
   */
  ErrorCode boundary_case( EntityHandle volume, int& result,
                             double u, double v, double w,
                             EntityHandle facet,
                             EntityHandle surface);

  /** get the solid angle projected by a facet on a unit sphere around a point
   *  - used by point_in_volume_slow
   */
  ErrorCode poly_solid_angle( EntityHandle face, const CartVect& point, double& area );

  /* SECTION III: Indexing & Cross-referencing */
public:
  /* Most calling apps refer to geometric entities with a combination of
   * base-1/0 ordinal index (or rank) and global ID (or name).
   * DagMC also has an internal EntityHandle reference to each geometric entity.
   * These method provide ways to translate from one to the other.
   */

  /** map from dimension & base-1 ordinal index to EntityHandle */
  EntityHandle entity_by_index( int dimension, int index );
  /** map from dimension & base-1 ordinal index to global ID */
  int id_by_index( int dimension, int index );
  /** map from dimension & global ID to EntityHandle */
  EntityHandle entity_by_id( int dimension, int id );
  /** PPHW: Missing dim & global ID ==> base-1 ordinal index */
  /** map from EntityHandle to base-1 ordinal index */
  int index_by_handle( EntityHandle handle );
  /** map from EntityHandle to global ID */
  int get_entity_id(EntityHandle this_ent);

  /**\brief get number of geometric sets corresponding to geometry of specified dimension
   *
   * For a given dimension (e.g. dimension=3 for volumes, dimension=2 for surfaces)
   * return the number of entities of that dimension
   *\param dimension the dimensionality of the entities in question
   *\return integer number of entities of that dimension
   */
  int num_entities( int dimension );

private:
  /** build internal index vectors that speed up handle-by-id, etc. */
  ErrorCode build_indices(Range &surfs, Range &vols);


  /* SECTION IV: Handling DagMC settings */
public:

  /** retrieve overlap thickness */
  double overlap_thickness() {return overlapThickness;}
  /** retrieve numerical precision */
  double numerical_precision() {return numericalPrecision;}
  /** retrieve faceting tolerance */
  double faceting_tolerance() {return facetingTolerance;}
  /** retrieve use CAD toggle */
  bool use_CAD() {return useCAD;}

  /** Attempt to set a new overlap thickness tolerance, first checking for sanity */
  void set_overlap_thickness( double new_overlap_thickness );

  /** Attempt to set a new numerical precision , first checking for sanity
   *  Use of this function is discouraged; see top of DagMC.cpp
   */
  void set_numerical_precision( double new_precision );

  /** attempt to set useCAD, first checking for availability */
  void set_use_CAD( bool use_cad );

  /* SECTION V: Metadata handling */
  /** Detect all the property keywords that appear in the loaded geometry
   *
   * @param keywords_out The result list of keywords.  This list could be
   *        validly passed to parse_properties().
   */
  ErrorCode detect_available_props( std::vector<std::string>& keywords_out, const char *delimiters = "_" );

  /** Parse properties from group names per metadata syntax standard
   *
   * @param keywords A list of keywords to parse.  These are considered the canonical
   *                 names of the properties, and constitute the valid inputs to
   *                 has_prop() and prop_value().
   * @param delimiters An array of characters the routine will use to split the groupname
   *                   into properties.
   * @param synonyms An optional mapping of synonym keywords to canonical keywords.
   *                 This allows more than one group name keyword to take on the same
   *                 meaning
   *                 e.g. if synonyms["rest.of.world"] = "graveyard", then volumes
   *                 in the "rest.of.world" group will behave as if they were in a
   *                 group named "graveyard".
   */
  ErrorCode parse_properties( const std::vector<std::string>& keywords,
                              const std::map<std::string,std::string>& synonyms = no_synonyms,
			      const char* delimiters = "_" );

  /** Get the value of a property on a volume or surface
   *
   * @param eh The entity handle to get a property value on
   * @param prop The canonical property name
   * @param value Output parameter, the value of the property.  If no value was
   *              set on the handle, this will be the empty string.
   * @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *         MOAB, or MB_SUCCESS if successful
   */
  ErrorCode prop_value( EntityHandle eh, const std::string& prop, std::string& value );

  /** Get the value of a property on a volume or surface
   *
   * @param eh The entity handle to get a property value on
   * @param prop The canonical property name
   * @param values Output parameter, the values of the property will be appended to this list.  If no value was
   *               set on the handle, no entries will be added.
   * @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *         MOAB, or MB_SUCCESS if successful
   */
  ErrorCode prop_values( EntityHandle eh, const std::string& prop,
                         std::vector< std::string >& value );

  /** Return true if a volume or surface has the named property set upon it
   *
   * @param eh The entity handle to query
   * @param prop The canonical property name
   * @retrun True if the handle has the property set, or false if not.
   *         False is also returned if a MOAB error occurs.
   */
  bool has_prop( EntityHandle eh, const std::string& prop );

  /** Get a list of all unique values assigned to a named property on any entity
   *
   * @param prop The canonical property name
   * @param return_list Output param, a list of unique strings that are set as values for this property
   * @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *         MOAB, or MB_SUCCESS if succesful
   */
  ErrorCode get_all_prop_values( const std::string& prop, std::vector<std::string>& return_list );

  /** Get a list of all entities which have a given property
   *
   * @param prop The canonical property name
   * @param return_list Output param, a list of entity handles that have this property
   * @param dimension If nonzero, entities returned will be restricted to the given dimension,
   *                  i.e. 2 for surfaces and 3 for volumes
   * @parm value If non-NULL, only entities for which the property takes on this value will be returned.
   * @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *         MOAB, or MB_SUCCESS if succesful
   */
  ErrorCode entities_by_property( const std::string& prop, std::vector<EntityHandle>& return_list,
                                  int dimension = 0, const std::string* value = NULL );

  bool is_implicit_complement(EntityHandle volume);

  /** get the tag for the "name" of a surface == global ID */
  Tag name_tag() {return nameTag;}

    // Get the tag used to associate OBB trees with geometry in load_file(..).
  Tag obb_tag() {return obbTag;}
  Tag geom_tag() {return geomTag;}
  Tag id_tag() {return idTag;}
  Tag sense_tag() { return senseTag; }

private:
  /** tokenize the metadata stored in group names - basically borroed from ReadCGM.cpp */
  void tokenize( const std::string& str,
                 std::vector<std::string>& tokens,
                 const char* delimiters = "_" ) const;

  // a common type within the property and group name functions
  typedef std::map<std::string, std::string> prop_map;

  /** Store the name of a group in a string */
  ErrorCode get_group_name( EntityHandle group_set, std::string& name );
  /** Parse a group name into a set of key:value pairs */
  ErrorCode parse_group_name( EntityHandle group_set, prop_map& result, const char* delimiters = "_");
  /** Add a string value to a property tag for a given entity */
  ErrorCode append_packed_string( Tag, EntityHandle, std::string& );
  /** Convert a property tag's value on a handle to a list of strings */
  ErrorCode unpack_packed_string( Tag tag, EntityHandle eh,
                                  std::vector< std::string >& values );

  std::vector<EntityHandle>& surf_handles() {return entHandles[2];}
  std::vector<EntityHandle>& vol_handles() {return entHandles[3];}
  std::vector<EntityHandle>& group_handles() {return entHandles[4];}

  Tag get_tag( const char* name, int size, TagType store, DataType type,
                 const void* def_value = NULL, bool create_if_missing = true);

  /* SECTION VI: Other */
public:
  OrientedBoxTreeTool *obb_tree() {return &obbTree;}

  ErrorCode write_mesh(const char* ffile,
                       const int flen);

    // get the corners of the OBB for a given volume
  ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]);

    // get the center point and three vectors for the OBB of a given volume
  ErrorCode getobb(EntityHandle volume, double center[3],
                     double axis1[3], double axis2[3], double axis3[3]);

    // get the root of the obbtree for a given entity
  ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle &root);

    // Get the instance of MOAB used by functions in this file.
  Interface* moab_instance() {return mbImpl;}


private:

  DagMC(Interface *mb_impl);

  static void create_instance(Interface *mb_impl = NULL);

  /* PRIVATE MEMBER DATA */

  static DagMC *instance_;
  static Interface *moab_instance_created;
  Interface *mbImpl;

  OrientedBoxTreeTool obbTree;
  EntityHandle impl_compl_handle;
  Tag obbTag, geomTag, idTag, nameTag, senseTag, facetingTolTag;

  std::vector<EntityHandle> entHandles[5];
    // store some lists indexed by handle
    // this is the lowest-valued handle among entity sets representing
    // surfaces & volumes
  EntityHandle setOffset;
    // list of obbTree root sets for surfaces and volumes,
    // indexed by [surf_or_vol_handle - setOffset]
  std::vector<EntityHandle> rootSets;
    // entity index (contiguous 1-N indices) indexed like rootSets are
  std::vector<int> entIndices;

    // corresponding geometric entities indexed like rootSets are
  std::vector<RefEntity *> geomEntities;

  // metadata
  // an empty synonym map to provide as a default argument to parse_properties()
  static const std::map<std::string,std::string> no_synonyms;
  // a map from the canonical property names to the tags representing them
  std::map<std::string, Tag> property_tagmap;

  char implComplName[NAME_TAG_SIZE];

  double overlapThickness;
  double numericalPrecision;
  double facetingTolerance, defaultFacetingTolerance;
  bool useCAD;         /// true if user requested CAD-based ray firing
  bool have_cgm_geom;  /// true if CGM contains problem geometry; required for CAD-based ray firing.

  // temporary storage so functions don't have to reallocate vectors
  // for ray_fire:
  std::vector<double> distList;
  std::vector<EntityHandle> prevFacetList, surfList, facetList;
  // for point_in_volume:
  std::vector<double> disList;
  std::vector<int>    dirList;
  std::vector<EntityHandle> surList, facList;

  // for (optional) counting
  long long int n_pt_in_vol_calls, n_ray_fire_calls;

};

inline DagMC *DagMC::instance(Interface *mb_impl)
{
  if (NULL == instance_) create_instance(mb_impl);

  return instance_;
}

inline void DagMC::destroy()
{
  if (NULL != instance_) {
    delete instance_;
    instance_ = NULL;
  }

  if (NULL != moab_instance_created) {
    delete moab_instance_created;
    moab_instance_created = NULL;
  }
}

inline EntityHandle DagMC::entity_by_index( int dimension, int index )
{
  assert(2 <= dimension && 3 >= dimension && (unsigned) index < entHandles[dimension].size());
  return entHandles[dimension][index];
}

inline int DagMC::index_by_handle( EntityHandle handle )
{
  assert(handle-setOffset < entIndices.size());
  return entIndices[handle-setOffset];
}

inline int DagMC::num_entities( int dimension )
{
  assert(0 <= dimension && 3 >= dimension);

  return entHandles[dimension].size() - 1;
}

    // get the root of the obbtree for a given entity
inline ErrorCode DagMC::get_root(EntityHandle vol_or_surf, EntityHandle &root)
{
  unsigned int index = vol_or_surf - setOffset;
  root = (index < rootSets.size() ? rootSets[index] : 0);
  return (root ? MB_SUCCESS : MB_INDEX_OUT_OF_RANGE);
}

} // namespace moab

#endif
