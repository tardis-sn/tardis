/** \mainpage The Mesh-Oriented datABase (MOAB)
 *
 * MOAB is a component for representing and evaluating mesh data.  MOAB can store 
 * structured and unstructured mesh, consisting of elements in the finite element “zoo”, 
 * along with polygons and polyhedra.  The functional interface to MOAB is simple, consisting 
 * of only four fundamental data types.  This data is quite powerful, allowing the representation 
 * of most types of metadata commonly found on the mesh.  MOAB is optimized for efficiency in 
 * space and time, based on access to mesh in chunks rather than through individual entities, 
 * while also versatile enough to support individual entity access.
 *
 * The MOAB data model consists of the following four fundamental types: mesh interface instance, 
 * mesh entities (vertex, edge, tri, etc.), sets, and tags.  Entities are addressed through handles 
 * rather than pointers, to allow the underlying representation of an entity to change without 
 * changing the handle to that entity.  Sets are arbitrary groupings of mesh entities and other 
 * sets.  Sets also support parent/child relationships as a relation distinct from sets containing 
 * other sets.  The directed-graph provided by set parent/child relationships is useful for modeling 
 * topological relations from a geometric model and other metadata.  Tags are named data which can 
 * be assigned to the mesh as a whole, individual entities, or sets.  Tags are a mechanism for 
 * attaching data to individual entities and sets are a mechanism for describing relations between 
 * entities; the combination of these two mechanisms is a powerful yet simple interface for 
 * representing metadata or application-specific data.  For example, sets and tags can be used 
 * together to describe geometric topology, boundary condition, and inter-processor interface 
 * groupings in a mesh.
 *
 * MOAB's API is documented in the moab::Interface class.  Questions and comments should be sent to moab-dev
 * _at_ mcs.anl.gov.
 *
 * \ref userguide "User's Guide"
 *
 * \ref developerguide "Developer's Guide"
 *
 * \ref metadata "I/O and Meta-Data Storage Conventions in MOAB"
 *
 * <a href="pages.html">Full List of Documents</a>
 */

#ifndef MOAB_INTERFACE_HPP
#define MOAB_INTERFACE_HPP

#define MOAB_API_VERSION 1.01
#define MOAB_API_VERSION_STRING "1.01"

#include "moab/MOABConfig.h"
#include "moab/Forward.hpp"
#include "moab/Range.hpp"
#include "moab/Compiler.hpp"
#include "moab/ErrorHandler.hpp"

// include files
#include <string>
#include <functional>
#include <typeinfo>

//! component architecture definitions
#ifdef XPCOM_MB

#ifndef __gen_nsISupports_h__
#include "nsISupports.h"
#endif

#ifndef NS_NO_VTABLE
#define NS_NO_VTABLE
#endif

#define MBINTERFACE_IID_STR "f728830e-1dd1-11b2-9598-fb9f414f2465"

#define MBINTERFACE_IID \
  {0xf728830e, 0x1dd1, 0x11b2, \
    { 0x95, 0x98, 0xfb, 0x9f, 0x41, 0x4f, 0x24, 0x65 }}

#endif


#include "moab/UnknownInterface.hpp"
#define MB_INTERFACE_VERSION "2.0.0"
namespace moab {

static const MBuuid IDD_MBCore = MBuuid( 0x8956e0a, 0xc300, 0x4005,
                                         0xbd, 0xf6, 0xc3, 0x4e, 0xf7, 0x1f, 0x5a, 0x52 );



/**
 * \class Interface Interface.hpp "moab/Interface.hpp"
 * \brief Main interface class to MOAB
 * \nosubgrouping
 */
#if defined(XPCOM_MB)
class NS_NO_VTABLE Interface : public nsISupports {
#else
class Interface : public UnknownInterface {
#endif

public:

#ifdef XPCOM_MB
  NS_DEFINE_STATIC_IID_ACCESSOR(MBINTERFACE_IID)
#endif

        /** \name Interface */

        /**@{*/

      //! constructor
  Interface() {}

    //! destructor
  virtual ~Interface() {}

    //! return the entity set representing the whole mesh
  virtual EntityHandle get_root_set()=0;
  
    //! Get a pointer to an internal MOAB interface
    //!\return NULL if not found, iterface pointer otherwise
  virtual ErrorCode query_interface_type( const std::type_info& iface_type, void*& iface ) = 0;
  
    //! Get a pointer to an internal MOAB interface
    //!\return NULL if not found, iterface pointer otherwise
  template <class IFace> ErrorCode query_interface(IFace*& ptr)
    { 
      void* tmp_ptr;
      ErrorCode result = query_interface_type(typeid(IFace), tmp_ptr);
      ptr = reinterpret_cast<IFace*>(tmp_ptr);
      return result;
    }
 
    //! Release reference to MB interface
  virtual ErrorCode release_interface_type( const std::type_info& iface_type, void* iface ) = 0;
  
  template <class IFace> ErrorCode release_interface(IFace* interface)
    { return release_interface_type( typeid(IFace), interface ); }
 
    //! Release reference to MB interface

    //! Returns the major.minor version number of the interface
    /**
       \param version_string If non-NULL, will be filled in with a string, possibly 
       containing implementation-specific information
    */
  virtual float api_version(std::string *version_string = NULL);
    
    //! Returns the major.minor version number of the implementation
    /**
       \param version_string If non-NULL, will be filled in with a string, possibly 
       containing implementation-specific information
    */
  virtual float impl_version(std::string *version_string = NULL)=0;

    /**@}*/

    /** \name Type and id */

    /**@{*/

    //! Returns the entity type of an EntityHandle.
    /** Returns the EntityType (ie, MeshVertex, MeshQuad, MeshHex ) of <em>handle</em>.
        \param handle The EntityHandle you want to find the entity type of.
        \return type The entity type of <em>handle</em>. 

        Example: \code
        EntityType type = type_from_handle( handle); 
        if( type == MeshHex ) ...  \endcode 
    */
  virtual EntityType type_from_handle(const EntityHandle handle) const = 0;
 
    //! Returns the id from an EntityHandle.
    /** \param handle The EntityHandle you want to find the id of. 
        \return id Id of <em>handle</em>
     
        Example: \code
        int id = id_from_handle(handle); \endcode 
    */
  virtual EntityID id_from_handle(const EntityHandle handle) const =0;

    //! Returns the topological dimension of an entity
    /** Returns the topological dimension of an entity.
        \param handle The EntityHandle you want to find the dimension of.
        \return type The topological dimension of <em>handle</em>. 

        Example: \code
        int dim = dimension_from_handle( handle); 
        if( dim == 0 ) ...  \endcode 
    */
  virtual int dimension_from_handle(const EntityHandle handle) const = 0;

    //! Gets an entity handle from the data base, if it exists, according to type and id.
    /** Given an EntiyType and an id, this function gets the existent EntityHandle. 
        If no such EntityHandle exits, it returns MB_ENTITY_NOT_FOUND 
        and sets handle to zero.
        \param type The type of the EntityHandle to retrieve from the database.
        \param id The id of the EntityHandle to retrieve from the database.
        \param handle An EntityHandle of type <em>type</em> and <em>id</em>. 

        Example: \code
        EntityType handle;
        ErrorCode error_code = handle_from_id(MeshTri, 204, handle );
        if( error_code == MB_ENTITY_NOT_FOUND ) ... \endcode
    */
  virtual ErrorCode handle_from_id(const EntityType type, 
                                     const EntityID, 
                                     EntityHandle& handle) const =0;

    /**@}*/

    /** \name Mesh input/output */

    /**@{*/

    //! Loads a mesh file into the database.
    /** Loads the file 'file_name'; types of mesh which can be loaded 
        depend on modules available at MB compile time.  If 
        active_block_id_list is NULL, all material sets (blocks in the 
        ExodusII jargon) are loaded.  Individual material sets  can be 
        loaded by specifying their ids in 'active_block_id_list'.  All 
        nodes are loaded on first call for a given file.  Subsequent 
        calls for a file load any material sets not loaded in previous 
        calls.
        \param file_name Name of file to load into database.
        \param active_block_id_list Material set/block ids to load.  
                If NULL, ALL blocks of <em>file_name</em> are loaded.
        \param num_blocks Number of blocks in active_block_id_list

        Example: \code
        std::vector<int> active_block_id_list;
        int active_block_id_list[] = {1, 4, 10};
        load_mesh( "temp.gen", active_block_id_list, 3 );  //load blocks 1, 4, 10 
        \endcode 
    */
  virtual ErrorCode load_mesh(const char *file_name,
                                const int *active_block_id_list = NULL,
                                const int num_blocks = 0)=0;

  /**\brief Load or import a file.
   *
   * Load a MOAB-native file or import data from some other supported
   * file format.
   *
   *\param file_name The location of the file to read.
   *\param file_set  If non-null, this argument must be a pointer to
   *                 a valid entity set handle.  All entities read from
   *                 the file will be added to this set.  File metadata
   *                 will be added to tags on the set.
   *\param options A list of string options, separated by semicolons (;).
   *               See README.IO for more information.  Options are typically
   *               format-specific options or parallel options.  If an
   *               option value is unrecognized but the file read otherwise
   *               succeeded, MB_UNHANDLED_OPTION will be returned.
   *\param set_tag_name The name of a tag used to designate the subset
   *               of the file to read.  The name must correspond to 
   *               data in the file that will be instantiated in MOAB
   *               as a tag.  
   *\param set_tag_values If the name specified in 'set_tag_name'
   *               corresponds to a tag with a single integer value,
   *               the values in this tag can be used to further
   *               limit the subset of data written from the file to
   *               only those entities or sets that have a value for
   *               the tag that is one of the values in this array.
   *\param num_set_tag_values The length of set_tag_values.
   *
   *\Note file_set is passed by pointer rather than by value (where a 
   *      zero handle value would indicate no set) so as to intentionally
   *      break compatibility with the previous version of this function
   *      because the behavior with respect to the file set was changed.
   *      The file_set is now an input-only argument.  The previous 
   *      version of this function unconditionally created a set and
   *      passed it back to the caller via a non-const reference argument.
   */
  virtual ErrorCode load_file( const char* file_name,
                                 const EntityHandle* file_set = 0,
                                 const char* options = 0,
                                 const char* set_tag_name = 0,
                                 const int* set_tag_values = 0,
                                 int num_set_tag_values = 0 ) = 0;

    //! Writes mesh to a file.
    /** Write mesh to file 'file_name'; if output_list is non-NULL, only 
        material sets contained in that list will be written.
        \param file_name Name of file to write.
        \param output_list 1d array of material set handles to write; if 
                           NULL, all sets are written
        \param num_sets Number of sets in output_list array

        Example: \code
        EntityHandle output_list[] = {meshset1, meshset2, meshset3}; 
        write_mesh( "output_file.gen", output_list, 3 ); \endcode 
    */
  virtual ErrorCode write_mesh(const char *file_name,
                                 const EntityHandle *output_list = NULL,
                                 const int num_sets = 0) = 0;

  /**\brief Write or export a file.
   * 
   * Write a MOAB-native file or export data to some other supported
   * file format.
   *
   *\param file_name The location of the file to write.
   *\param file_type The type of the file.  If this value is NULL, 
   *                 then file type will be determined using the 
   *                 file name suffix.  
   *\param options   A semicolon-separated list of options.
   *                 See README.IO for more information.  Typical options
   *                 include the file type, parallel options, and options
   *                 specific to certain file formats.
   *\param output_sets A list of entity sets to write to the file.  If
   *                 no sets are sepcified, the default behavior is to
   *                 write all data that is supported by the target file
   *                 type.
   *\param num_output_sets The length of the output_sets array.
   *\param tag_list A list of tags for which to write the tag data.  The
   *                write may fail if a tag list is specified but the 
   *                target file type is not capable of representing the
   *                data.  If no tags are specified, the default is to
   *                write whatever data the target file format supports.
   *\param num_tags The length of tag_list.
   */
  virtual ErrorCode write_file( const char* file_name,
                                  const char* file_type = 0,
                                  const char* options = 0,
                                  const EntityHandle* output_sets = 0,
                                  int num_output_sets = 0,
                                  const Tag* tag_list = 0,
                                  int num_tags = 0 ) = 0;

  /**\brief Write or export a file.
   * 
   * Write a MOAB-native file or export data to some other supported
   * file format.
   *
   *\param file_name The location of the file to write.
   *\param file_type The type of the file.  If this value is NULL, 
   *                 then file type will be determined using the 
   *                 file name suffix.  
   *\param options   A semicolon-separated list of options.
   *                 See README.IO for more information.  Typical options
   *                 include the file type, parallel options, and options
   *                 specific to certain file formats.
   *\param output_sets A list of entity sets to write to the file.  If
   *                 no sets are sepcified, the default behavior is to
   *                 write all data that is supported by the target file
   *                 type.
   *\param tag_list A list of tags for which to write the tag data.  The
   *                write may fail if a tag list is specified but the 
   *                target file type is not capable of representing the
   *                data.  If no tags are specified, the default is to
   *                write whatever data the target file format supports.
   *\param num_tags The length of tag_list.
   */
  virtual ErrorCode write_file( const char* file_name,
                                  const char* file_type,
                                  const char* options,
                                  const Range& output_sets,
                                  const Tag* tag_list = 0,
                                  int num_tags = 0 ) = 0;

    //! Deletes all mesh entities from this MB instance
  virtual ErrorCode delete_mesh()=0;

    /**@}*/

    /** \name Coordinates and dimensions */

    /**@{*/

    //! Get blocked vertex coordinates for all vertices
    /** Blocked = all x, then all y, etc. 
          
    Example: \code
    std::vector<double> coords;
    get_vertex_coordinates(coords);
    double xavg = 0;
    for (int i = 0; i < coords.size()/3; i++) xavg += coords[i]; \endcode
    */
  virtual ErrorCode get_vertex_coordinates(std::vector<double> &coords) const =0;

    //! get pointers to coordinate data
    /** BEWARE, THIS GIVES ACCESS TO MOAB'S INTERNAL STORAGE, USE WITH CAUTION!
     * This function returns pointers to MOAB's internal storage for vertex coordinates.
     * Access is similar to tag_iterate, see documentation for that function for details
     * about arguments and a coding example.
     */
  virtual ErrorCode coords_iterate(Range::const_iterator iter,
                                     /**< Iterator to first entity you want coordinates for */
                                   Range::const_iterator end,
                                     /**< Iterator to last entity you want coordinates for */
                                   double*& xcoords_ptr,
                                     /**< Pointer to x coordinate storage for these entities */
                                   double*& ycoords_ptr,
                                     /**< Pointer to y coordinate storage for these entities */
                                   double*& zcoords_ptr,
                                     /**< Pointer to z coordinate storage for these entities */
                                   int& count
                                     /**< Number of entities for which returned pointers are valid/contiguous */
                                   ) = 0;

    //! Gets xyz coordinate information for range of vertices
    /** Length of 'coords' should be at least 3*<em>entity_handles.size()</em> before making call.
        \param entity_handles Range of vertex handles (error if not of type MeshVertex)
        \param coords Array used to return x, y, and z coordinates.
   
        Example: \code 
        double coords[3];
        get_coords( vertex_handle, coords ); 
        std::cout<<"x = "<<coords[0]<<std::endl;
        std::cout<<"y = "<<coords[1]<<std::endl;
        std::cout<<"z = "<<coords[2]<<std::endl; \endcode 
    */
  virtual ErrorCode  get_coords(const Range& entity_handles, 
                                  double *coords) const =0;
    
    //! Gets xyz coordinate information for vector of vertices
    /** Identical to range-based function, except entity handles are specified using a 1d vector
        and vector length.
    */
  virtual ErrorCode  get_coords(const EntityHandle* entity_handles, 
                                  const int num_entities, 
                                  double *coords) const =0;
  
  /**\brief Get vertex coordinates in blocks by dimension.
   *
   * Get the X, Y, and Z coordinates of a group of vertices.  
   * Coordinates are returned in separate arrays, one for each 
   * dimension.  Each coordinate array must be of sufficient
   * length to hold the coordinate value for each vertex.  Array
   * pointers may be NULL if coordinates in the the respective 
   * dimension are not desired.
   *\param entity_handles  The group of vertex handles for which to get the coordiantes.
   *\param x_coords        Output: the X coordinate of each vertex.  May be NULL.
   *\param y_coords        Output: the Y coordinate of each vertex.  May be NULL.
   *\param z_coords        Output: the Z coordinate of each vertex.  May be NULL.
   */
  virtual ErrorCode get_coords( const Range& entity_handles,
                                  double* x_coords,
                                  double* y_coords,
                                  double* z_coords ) const = 0;
  
  
    //! Sets the xyz coordinates for a vector of vertices
    /** An error is returned if any entities in the vector are not vertices.
        \param entity_handles EntityHandle's to set coordinates of. (Must be of type MeshVertex)
        \param num_entities Number of entities in entity_handles
        \param coords Array containing new xyz coordinates.
 
        Example: \code
        double coords[3] = {0.234, -2.52, 12.023};
        set_coords( entity_handle, 1, coords ); \endcode 
    */
  virtual ErrorCode  set_coords(const EntityHandle *entity_handles, 
                                  const int num_entities,
                                  const double *coords)=0;

    //! Sets the xyz coordinates for a vector of vertices
    /** An error is returned if any entities in the vector are not vertices.
        \param entity_handles EntityHandle's to set coordinates of. (Must be of type MeshVertex)
        \param num_entities Number of entities in entity_handles
        \param coords Array containing new xyz coordinates.
 
        Example: \code
        double coords[3] = {0.234, -2.52, 12.023};
        set_coords( entity_handle, 1, coords ); \endcode 
    */
  virtual ErrorCode  set_coords(Range entity_handles,
                                  const double *coords)=0;

    //! Get overall geometric dimension
  virtual ErrorCode get_dimension(int &dim) const =0;

    //! Set overall geometric dimension
    /** Returns error if setting to 3 dimensions, mesh has been created, and 
     *  there are only 2 dimensions on that mesh
     */
  virtual ErrorCode set_dimension(const int dim)=0;

    /**@}*/

    /** \name Connectivity */

    /**@{*/

    //! get pointers to connectivity data
    /** BEWARE, THIS GIVES ACCESS TO MOAB'S INTERNAL STORAGE, USE WITH CAUTION!
     * This function returns a pointer to MOAB's internal storage for entity connectivity.
     * For each contiguous sub-range of entities, those entities are guaranteed to have
     * the same number of vertices (since they're in the same ElementSequence).  Count
     * is given in terms of entities, not elements of the connectivity array.
     * Access is similar to tag_iterate, see documentation for that function for details
     * about arguments and a coding example.
     */
  virtual ErrorCode connect_iterate(Range::const_iterator iter,
                                      /**< Iterator to first entity you want coordinates for */
                                    Range::const_iterator end,
                                      /**< Iterator to last entity you want coordinates for */
                                    EntityHandle *&connect,
                                      /**< Pointer to connectivity storage for these entities */
                                    int &verts_per_entity,
                                      /**< Number of vertices per entity in this block of entities */
                                    int& count
                                      /**< Number of entities for which returned pointers are valid/contiguous */
                                    ) = 0;
  
    //! Get the connectivity array for all entities of the specified entity type
    /**  This function returns the connectivity of just the corner vertices, no higher order nodes
         \param type The entity type of elements whose connectivity is to be returned
         \param connect an STL vector used to return connectivity array (in the form of entity handles)
    */
  virtual ErrorCode get_connectivity_by_type(const EntityType type, 
                                               std::vector<EntityHandle> &connect) const =0;

    //! Gets the connectivity for a vector of elements
    /** Same as vector-based version except range is returned (unordered!)
    */
  virtual ErrorCode  get_connectivity(const EntityHandle *entity_handles, 
                                        const int num_handles,
                                        Range &connectivity, 
                                        bool corners_only = false) const =0;

    //! Gets the connectivity for elements
    /** Same as vector-based version except range is returned (unordered!)
    */
  virtual ErrorCode get_connectivity( const Range& entity_handles, 
                                        Range &connectivity, 
                                        bool corners_only = false) const =0;
 
    //! Gets the connectivity for a vector of elements
    /** Corner vertices or all vertices (including higher-order nodes, if any) are returned.
        For non-element handles (ie, MB_MeshSets), returns an error. Connectivity data is copied 
        from the database into the vector.  Connectivity of a vertex is the same vertex.
        The nodes in <em>connectivity</em> are properly ordered according to that element's 
        canonical ordering.
        \param entity_handles Vector of element handles to get connectivity of.
        \param num_handles Number of entity handles in <em>entity_handles</em>
        \param connectivity Vector in which connectivity of <em>entity_handles</em> is returned.  
        \param corners_only If true, returns only corner vertices, otherwise returns all of them (including any higher-order vertices)
        \param offsets If non-NULL, offsets->[i] stores the index of the start of entity i's connectivity,
                with the last value in offsets one beyond the last entry
    */
  virtual ErrorCode  get_connectivity(const EntityHandle *entity_handles, 
                                      const int num_handles,
                                      std::vector<EntityHandle> &connectivity, 
                                      bool corners_only = false,
                                      std::vector<int> *offsets = NULL) const =0;
 
    //! Gets a pointer to constant connectivity data of <em>entity_handle</em> 
    /** Sets <em>number_nodes</em> equal to the number of nodes of the <em> 
        entity_handle </em>.  Faster then the other <em>get_connectivity</em> function because no
        data is copied.  The nodes in 'connectivity' are properly ordered according to the 
        element's canonical ordering.
        

          Example: \code 
          const EntityHandle* conn;
          int number_nodes = 0;
          get_connectivity( entity_handle, conn, number_nodes ); \endcode 
          
          Example2: \code
          std::vector<EntityHandle> sm_storage;
          const EntityHandle* conn;
          int number_nodes;
          get_connectivity( handle, conn, number_nodes, false, &sm_storage );
          if (conn == &sm_storage[0])
            std::cout << "Structured mesh element" << std::endl;
          \endcode
        
        \param entity_handle EntityHandle to get connectivity of.
        \param connectivity Array in which connectivity of <em>entity_handle</em> is returned.
        \param num_nodes Number of MeshVertices in array <em>connectivity</em>. 
        \param corners_only If true, returns only corner vertices, otherwise returns all of them (including any higher-order vertices)
        \param storage Some elements (e.g. structured mesh) may not have an
                       explicit connectivity list.  This function will normally
                       return MB_NOT_IMPLEMENTED for such elements.  However,
                       if the caller passes in a non-null value for this 
                       argument, space will be allocated in this vector for
                       the connectivity data and the connectivity pointer will
                       be set to the data in this vector.
    */
  virtual ErrorCode  get_connectivity(const EntityHandle entity_handle, 
                                        const EntityHandle *&connectivity, 
                                        int &num_nodes, 
                                        bool corners_only = false,
                                        std::vector<EntityHandle>* storage = 0
                                        ) const =0;

    //! Sets the connectivity for an EntityHandle.  For non-element handles, return an error.
    /** Connectivity is stored exactly as it is ordered in vector <em>connectivity</em>. 
        \param entity_handle EntityHandle to set connectivity of.
        \param connect Vector containing new connectivity of <em>entity_handle</em>.
        \param num_connect Number of vertices in <em>connect</em>
   
        Example: \code 
        EntityHandle conn[] = {node1, node2, node3};
        set_connectivity( tri_element, conn, 3 ); \endcode 
    */
  virtual ErrorCode  set_connectivity(const EntityHandle entity_handle, 
                                        EntityHandle *connect,
                                        const int num_connect)=0;

    /**@}*/

    /** \name Adjacencies */

    /**@{*/

    //! Get the adjacencies associated with a vector of entities to entities of a specfied dimension.
    /** \param from_entities Vector of EntityHandle to get adjacencies of.
        \param num_entities Number of entities in <em>from_entities</em>
        \param to_dimension Dimension of desired adjacencies
        \param create_if_missing If true, MB will create any entities of the specfied dimension
        which have not yet been created (only useful when <em>to_dimension < dim(*from_entities)</em>)
        \param adj_entities STL vector to which adjacent entities are appended. 
        \param operation_type Enum of INTERSECT or UNION.  Defines whether to take
        the intersection or union of the set of adjacencies recovered for the from_entities.

        The adjacent entities in vector <em>adjacencies</em> are not in any particular 
        order. 

        Example: \code
        std::vector<EntityHandle> adjacencies, from_entities = {hex1, hex2};
          // generate all edges for these two hexes
          get_adjacencies( from_entities, 2, 1, true, adjacencies, Interface::UNION); 
          adjacencies.clear();
            // now find the edges common to both hexes
            get_adjacencies( from_entities, 2, 1, false, adjacencies, Interface::INTERSECT); 
            \endcode 
    */

  virtual ErrorCode get_adjacencies(const EntityHandle *from_entities,
                                      const int num_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      std::vector<EntityHandle>& adj_entities,
                                      const int operation_type = Interface::INTERSECT) = 0;

    //! Get the adjacencies associated with a vector of entities to entities of a specfied dimension.
    /** Identical to vector-based get_adjacencies function, except results are returned in a
        range instead of a vector.
    */
  virtual ErrorCode get_adjacencies(const EntityHandle *from_entities,
                                      const int num_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      Range &adj_entities,
                                      const int operation_type = Interface::INTERSECT) = 0;

    //! Get the adjacencies associated with a range of entities to entities of a specfied dimension.
    /** Identical to vector-based get_adjacencies function, except "from" entities specified in a
        range instead of a vector.
    */
  virtual ErrorCode get_adjacencies(const Range &from_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      Range &adj_entities,
                                      const int operation_type = Interface::INTERSECT) = 0;

    //! Adds adjacencies between "from" and "to" entities.
    /** \param from_handle Entities on which the adjacencies are placed
        \param to_handles Vector of entities referenced by new adjacencies added to <em>from_handle</em>
        \param num_handles Number of entities in <em>to_handles</em>
        \param both_ways If true, add the adjacency information in both directions; if false,
        adjacencies are added only to <em>from_handle</em>
    */
  virtual ErrorCode add_adjacencies(const EntityHandle from_handle, 
                                      const EntityHandle *to_handles,
                                      const int num_handles,
                                      bool both_ways) = 0;

    //! Adds adjacencies; same as vector-based, but with range instead
  virtual ErrorCode add_adjacencies(const EntityHandle from_handle, 
                                      Range &adjacencies,
                                      bool both_ways) = 0;

    //! Removes adjacencies between handles.
    /** Adjacencies in both directions are removed.
        \param from_handle Entity from which adjacencies are being removed.
        \param to_handles Entities to which adjacencies are being removed.
        \param num_handles Number of handles in <em>to_handles</em>
    */
  virtual ErrorCode remove_adjacencies(const EntityHandle from_handle, 
                                         const EntityHandle *to_handles,
                                         const int num_handles) = 0;

    /**\brief Get a ptr to adjacency lists
     * Get a pointer to adjacency lists.  These lists are std::vector<EntityHandle>, which are pointed
     * to by adjs[i].  Adjacencies are not guaranteed to be in order of increasing dimension.  Only a 
     * const version of this function is given, because adjacency data is managed more carefully in MOAB
     * and should be treated as read-only by applications.  If adjacencies have not yet been initialized,
     * adjs_ptr will be NULL (i.e. adjs_ptr == NULL).  There may also be NULL entries for individual entities,
     * i.e. adjs_ptr[i] == NULL.
     * \param iter Iterator to beginning of entity range desired
     * \param end End iterator for which adjacencies are requested
     * \param adjs_ptr Pointer to pointer to const std::vector<EntityHandle>; each member of that array is 
     *                 the vector of adjacencies for this entity
     * \param count Number of entities in the contiguous chunk starting from *iter
     */
  virtual ErrorCode adjacencies_iterate(Range::const_iterator iter,
                                        Range::const_iterator end,
                                        const std::vector<EntityHandle> **& adjs_ptr,
                                        int& count) = 0;
    /**@}*/

    //! Enumerated type used in get_adjacencies() and other functions
  enum {INTERSECT, UNION};

    /** \name Getting entities */

    /**@{*/

    //! Retrieves all entities of a given topological dimension in the database or meshset.
    /** Appends entities to list passed in.
        \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param dimension Topological dimension of entities desired.
        \param entities Range in which entities of dimension <em>dimension</em> are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.

        Example: \code
          // get 1d (edge) elements in the entire mesh
          Range edges;
          get_entities_by_dimension( 0, 1, edges );
          \endcode 
    */
  virtual ErrorCode get_entities_by_dimension(const EntityHandle meshset,
                                                const int dimension, 
                                                Range &entities,
                                                const bool recursive = false)  const = 0;

  virtual ErrorCode get_entities_by_dimension(const EntityHandle meshset,
                                                const int dimension, 
                                                std::vector<EntityHandle> &entities,
                                                const bool recursive = false)  const = 0;

    //! Retrieve all entities of a given type in the database or meshset.
    /** Appends entities to list passed in.
        \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param entities Range in which entities of type <em>type</em> are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBENTITYSET is an error, as it would always 
                         result in an empty list.

        Example: \code
          // get the quadrilateral elements in meshset
          Range quads;
          get_entities_by_type( meshset, MeshQuad, quads );
          \endcode 
    */
  virtual ErrorCode get_entities_by_type(const EntityHandle meshset,
                                           const EntityType type, 
                                           Range &entities,
                                           const bool recursive = false) const = 0;

  virtual ErrorCode get_entities_by_type(const EntityHandle meshset,
                                           const EntityType type, 
                                           std::vector<EntityHandle> &entities,
                                           const bool recursive = false) const = 0;

    //! Retrieve entities in the database or meshset which have any or all of the tag(s) and (optionally)
    //! value(s) specified.
    /** \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param tag_handles Vector of tag handles entities must have
        \param values Vector of pointers to values of tags in <em>tag_handles</em>
        \param num_tags Number of tags and values in <em>tag_handles</em> and <em>values</em>
        \param entities Range in which entities are returned.
        \param condition Boolean condition, either Interface::UNION or Interface::INTERSECT
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBENTITYSET is an error, as it would always 
                         result in an empty list.

        If Interface::UNION is specified as the condition, entities with <em>any</em> of the tags
        and values specified are returned.  If Interface::INTERSECT is specified, only entities with
        <em>all</em> of the tags/values are returned.

        If <em>values</em> is NULL, entities with the specified tags and any corresponding values are
        returned.  Note that if <em>values</em> is non-NULL, it is a vector of <em>pointers</em> to
        tag values.

        Example: \code
          // get the dirichlet sets in a mesh
          Range dir_sets;
          Tag dir_tag;
          tag_get_handle(DIRICHLET_SET_TAG_NAME, dir_tag, 1, MB_TYPE_INTEGER);
          get_entities_by_type_and_tag(0, MeshEntitySet, &dir_tag, NULL, 1, dir_sets, 
          Interface::UNION);
          \endcode 
    */
  virtual ErrorCode get_entities_by_type_and_tag(const EntityHandle meshset,
                                                   const EntityType type,
                                                   const Tag *tag_handles,
                                                   const void* const* values,
                                                   const int num_tags,
                                                   Range &entities,
                                                   const int condition = Interface::INTERSECT,
                                                   const bool recursive = false) const = 0;

    //! Returns all entities in the data base or meshset, in a range (order not preserved)
    /** Appends entities to list passed in.
        \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param entities Range in which entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.

        Example: \code
        Range entities;
          // get all non-meshset entities in meshset, including in contained meshsets
          get_entities_by_handle(meshset, entities, true);
          \endcode 
    */
  virtual ErrorCode get_entities_by_handle(const EntityHandle meshset,
                                             Range &entities,
                                             const bool recursive = false) const = 0;

    //! Returns all entities in the data base or meshset, in a vector (order preserved)
    /** Appends entities to list passed in.
        \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param entities STL vector in which entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.

        Example: \code
        std::vector<EntityHandle> entities;
          // get all non-meshset entities in meshset, including in contained meshsets
          get_entities_by_handle(meshset, entities, true);
          \endcode 
    */
  virtual ErrorCode get_entities_by_handle(const EntityHandle meshset,
                                             std::vector<EntityHandle> &entities,
                                             const bool recursive = false) const = 0;

    //! Return the number of entities of given dimension in the database or meshset
    /** \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param dimension Dimension of entities desired.
        \param num_entities Number of entities of the given dimension
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual ErrorCode get_number_entities_by_dimension(const EntityHandle meshset,
                                                       const int dimension, 
                                                       int &num_entities,
                                                       const bool recursive = false) const = 0;

    //! Retrieve the number of entities of a given type in the database or meshset.
    /** Identical to get_entities_by_dimension, except returns number instead of entities
        \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param num_entities Number of entities of type <em>type</em>
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBENTITYSET is an error, as it would always 
                         result in an empty list.
    */
  virtual ErrorCode get_number_entities_by_type(const EntityHandle meshset,
                                                  const EntityType type, 
                                                  int &num_entities,
                                                  const bool recursive = false) const = 0;

    //! Retrieve number of entities in the database or meshset which have any or all of the 
    //! tag(s) and (optionally) value(s) specified.
    /** Identical to get_entities_by_type_and_tag, except number instead of entities are returned
        \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param tag_handles Vector of tag handles entities must have
        \param values Vector of pointers to values of tags in <em>tag_handles</em>
        \param num_tags Number of tags and values in <em>tag_handles</em> and <em>values</em>
        \param num_entities Range in which number of entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBENTITYSET is an error, as it would always 
                         result in an empty list.
    */
  virtual ErrorCode get_number_entities_by_type_and_tag(const EntityHandle meshset,
                                                          const EntityType type,
                                                          const Tag *tag_handles,
                                                          const void* const* values,
                                                          const int num_tags,
                                                          int &num_entities,
                                                          const int condition = Interface::INTERSECT,
                                                          const bool recursive = false) const = 0;

    //! Returns number of entities in the data base or meshset
    /** Identical to get-entities_by_handle, except number instead of entities are returned
        \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param num_entities Range in which num_entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual ErrorCode get_number_entities_by_handle(const EntityHandle meshset,
                                                    int &num_entities,
                                                    const bool recursive = false) const = 0;

    /**@}*/

    /** \name Mesh modification */

    /**@{*/

    //! Create an element based on the type and connectivity. 
    /** Create a new element in the database.  Vertices composing this element must already exist,
        and connectivity must be specified in canonical order for the given element type.  If 
        connectivity vector is not correct for EntityType <em>type</em> (ie, a vector with 
        3 vertices is passed in to make an MeshQuad), the function returns MB_FAILURE. 
        \param type Type of element to create. (MeshTet, MeshTri, MeshKnife, etc.) 
        \param connectivity 1d vector containing connectivity of element to create.
        \param num_vertices Number of vertices in element
        \param element_handle Handle representing the newly created element in the database.

        Example: \code
        EntityHandle quad_conn[] = {vertex0, vertex1, vertex2, vertex3};
        EntityHandle quad_handle = 0;
        create_element( MeshQuad, quad_conn, 4, quad_handle ); \endcode 
    */
  virtual ErrorCode create_element(const EntityType type, 
                                     const EntityHandle *connectivity,
                                     const int num_vertices, 
                                     EntityHandle &element_handle) = 0;

    //! Creates a vertex with the specified coordinates.  
    /**
       \param coordinates Array that has 3 doubles in it.
       \param entity_handle EntityHandle representing the newly created vertex in the database.

       Example: \code
       double coordinates[] = {1.034, 23.23, -0.432};
       EntityHandle new_handle = 0;
       create_vertex( coordinates, entity_handle ); \endcode 
    */
  virtual ErrorCode create_vertex(const double coordinates[3], 
                                    EntityHandle &entity_handle ) = 0;

    //! Create a set of vertices with the specified coordinates
    /**
       \param coordinates Array that has 3*n doubles in it.
       \param nverts Number of vertices to create
       \param entity_handles Range passed back with new vertex handles
    */
  virtual ErrorCode create_vertices(const double *coordinates, 
                                      const int nverts,
                                      Range &entity_handles ) = 0;

    //! Merge two entities into a single entity
    /** Merge two entities into a single entities, with <em>entity_to_keep</em> receiving
        adjacencies that were on <em>entity_to_remove</em>.
        \param entity_to_keep Entity to be kept after merge
        \param entity_to_remove Entity to be merged into <em>entity_to_keep</em>
        \param auto_merge If false, <em>entity_to_keep</em> and <em>entity_to_remove</em> must share
        the same lower-dimensional entities; if true, MB tries to merge those entities automatically
        \param delete_removed_entity If true, <em>entity_to_remove</em> is deleted after merge is complete
    */
  virtual ErrorCode merge_entities(EntityHandle entity_to_keep, 
                                     EntityHandle entity_to_remove,
                                     bool auto_merge,
                                     bool delete_removed_entity) = 0;

    //! Removes entities in a vector from the data base.  
    /** If any of the entities are contained in any meshsets, it is removed from those meshsets 
        which were created with MESHSET_TRACK_OWNER option bit set.  Tags for <em>entity</em> are 
        removed as part of this function.
        \param entities 1d vector of entities to delete
        \param num_entities Number of entities in 1d vector
    */ 
  virtual ErrorCode delete_entities(const EntityHandle *entities,
                                      const int num_entities) = 0;

    //! Removes entities in a range from the data base.  
    /** If any of the entities are contained in any meshsets, it is removed from those meshsets 
        which were created with MESHSET_TRACK_OWNER option bit set.  Tags for <em>entity</em> are 
        removed as part of this function.
        \param entities Range of entities to delete
    */ 
  virtual ErrorCode delete_entities(const Range &entities) = 0;

    /**@}*/

    /** \name Information */

    /**@{*/

    //! List entities to standard output
    /** Lists all data pertaining to entities (i.e. vertex coordinates if vertices, connectivity if
        elements, set membership if set).  Useful for debugging, but output can become quite long
        for large databases.
    */
  virtual ErrorCode list_entities(const Range &entities) const = 0;
  
    //! List entities, or number of entities in database, to standard output
    /** Lists data pertaining to entities to standard output.  If <em>entities</em> is NULL and
        <em>num_entities</em> is zero, lists only the number of entities of each type in the 
        database.  If <em>entities</em> is NULL and <em>num_entities</em> is non-zero, lists all
        information for all entities in the database.
        \param entities 1d vector of entities to list
        \param num_entities Number of entities in 1d vector
    */
  virtual ErrorCode list_entities(const EntityHandle *entities,
                                    const int num_entities) const = 0;

    //! List a single entity; no header printed
    /** Lists a single entity, including its connectivity and its adjacencies.
     *  No header is printed, because calling function might print information between header
     *  and information printed by this function.
     *  \param entity The entity to be listed.
     */
  virtual ErrorCode list_entity(const EntityHandle entity) const = 0;

    //! Return information about the last error
    /** \param info std::string into which information on the last error is written.
     */
  virtual ErrorCode get_last_error(std::string& info) const = 0;

    //! Return string representation of given error code
    /** \param code Error code for which string is wanted
     */
  virtual std::string get_error_string(const ErrorCode code) const = 0;

  /**\brief Calculate amount of memory used to store MOAB data
   *
   * This function calculates the amount of memory used to store
   * MOAB data.  
   *
   * There are two possible values for each catagory of memory use.
   * The exact value and the amortized value.  The exact value is the
   * amount of memory used to store the data for the specified entities.
   * The amortized value includes the exact value and an amortized 
   * estimate of the memory consumed in overhead for storing the values
   * (indexing structures, access structures, etc.)  
   *
   * Note: If ent_array is NULL and total_amortized_storage is *not* NULL,
   *       the total memory used by MOAB for storing data all will be 
   *       returned in the address pointed to by total_amortized_storage.
   *
   *\param ent_array Array of entities for which to estimate the memory use.
   *                 If NULL, estimate is done for all entities.
   *\param num_ents The length of ent_array.  Not used if ent_rray is NULL.
   *\param total_(amortized_)storage The sum of the memory entity, adjacency,
   *                   and all tag storage.
   *\param (amortized_)entity_storage The storage for the entity definitions
   *                   (connectivity arrays for elements, coordinates for 
   *                   vertices, list storage within sets, etc.)
   *\param (amortized_)adjacency_storage The storage for adjacency data.
   *\param tag_array   An array of tags for which to calculate the memory use.
   *\param num_tags    The lenght of tag_array
   *\param (amortized_)tag_storage If tag_array is not NULL, then one value
   *                   for each tag specifying the memory used for storing 
   *                   that tag.  If tag_array is NULL and this value is not,
   *                   the location at which to store the total memory used
   *                   for all tags.
   */
  virtual void estimated_memory_use( const EntityHandle* ent_array = 0,
                             unsigned long  num_ents = 0,
                             unsigned long long* total_storage = 0,
                             unsigned long long* total_amortized_storage = 0,
                             unsigned long long* entity_storage = 0,
                             unsigned long long* amortized_entity_storage = 0,
                             unsigned long long* adjacency_storage = 0,
                             unsigned long long* amortized_adjacency_storage = 0,
                             const Tag*   tag_array = 0,
                             unsigned       num_tags = 0,
                             unsigned long long* tag_storage = 0,
                             unsigned long long* amortized_tag_storage = 0 ) = 0;

  /**\brief Calculate amount of memory used to store MOAB data
   *
   * This function calculates the amount of memory used to store
   * MOAB data.  
   *
   * There are two possible values for each catagory of memory use.
   * The exact value and the amortized value.  The exact value is the
   * amount of memory used to store the data for the specified entities.
   * The amortized value includes the exact value and an amortized 
   * estimate of the memory consumed in overhead for storing the values
   * (indexing structures, access structures, etc.)  
   *
   *\param ents        Entities for which to estimate the memory use.
   *\param total_(amortized_)storage The sum of the memory entity, adjacency,
   *                   and all tag storage.
   *\param (amortized_)entity_storage The storage for the entity definitions
   *                   (connectivity arrays for elements, coordinates for 
   *                   vertices, list storage within sets, etc.)
   *\param (amortized_)adjacency_storage The storage for adjacency data.
   *\param tag_array   An array of tags for which to calculate the memory use.
   *\param num_tags    The lenght of tag_array
   *\param (amortized_)tag_storage If tag_array is not NULL, then one value
   *                   for each tag specifying the memory used for storing 
   *                   that tag.  If tag_array is NULL and this value is not,
   *                   the location at which to store the total memory used
   *                   for all tags.
   */
  virtual void estimated_memory_use( const Range& ents,
                             unsigned long long* total_storage = 0,
                             unsigned long long* total_amortized_storage = 0,
                             unsigned long long* entity_storage = 0,
                             unsigned long long* amortized_entity_storage = 0,
                             unsigned long long* adjacency_storage = 0,
                             unsigned long long* amortized_adjacency_storage = 0,
                             const Tag*   tag_array = 0,
                             unsigned       num_tags = 0,
                             unsigned long long* tag_storage = 0,
                             unsigned long long* amortized_tag_storage = 0 ) = 0;
    /**@}*/

    /** \name Higher-order elements */

    /**@{*/

    //! function object for recieving events from MB of higher order nodes added to entities
  class HONodeAddedRemoved
  {
  public:
      //! Constructor
    HONodeAddedRemoved(){}
 
      //! Destructor
    virtual ~HONodeAddedRemoved(){}

      //! node_added called when a node was added to an element's connectivity array
      //! note: connectivity array of element may be incomplete (corner nodes will exist always)
      /** 
       * \param node Node being added
       * \param element Element node is being added to
       */
    virtual void node_added(EntityHandle node, EntityHandle element) = 0;

      //! node_added called when a node was added to an element's connectivity array
      //! note: connectivity array of element may be incomplete (corner nodes will exist always)
      /**
       * \param node Node being removed.
       */
    virtual void node_removed(EntityHandle node) = 0;
  };
  
    //! Convert entities to higher-order elements by adding mid nodes
    /** This function causes MB to create mid-nodes on all edges, faces, and element interiors 
        for all entities in <em>meshset</em>.  Higher order nodes appear in an element's connectivity
        array according to the algorithm described in the documentation for Mesh.  If 
        <em>HONodeAddedRemoved</em> function is input, this function is called to notify the application
        of nodes being added/removed from the mesh.
        \param meshset The set of entities being converted
        \param mid_edge If true, mid-edge nodes are created 
        \param mid_face If true, mid-face nodes are created 
        \param mid_region If true, mid-element nodes are created 
        \param function_object If non-NULL, the node_added or node_removed functions on this object 
        are called when nodes are added or removed from an entity, respectively
    */
  virtual ErrorCode convert_entities(const EntityHandle meshset, 
                                       const bool mid_edge,
                                       const bool mid_face, 
                                       const bool mid_region, 
                                       HONodeAddedRemoved* function_object = 0) = 0;

    //! Returns the side number, in canonical ordering, of <em>child</em> with respect to <em>parent</em>
    /** Given a parent and child entity, returns the canonical ordering information side number, sense, 
        and offset of <em>child</em> with respect to <em>parent</em>.  This function returns
        MB_FAILURE if <em>child</em> is not related to <em>parent</em>.  This function does *not* 
        create adjacencies between <em>parent</em> and <em>child</em>.
        \param parent Parent entity to be compared
        \param child Child entity to be compared
        \param side_number Side number in canonical ordering of <em>child</em> with respect to 
        <em>parent</em>
        \param sense Sense of <em>child</em> with respect to <em>parent</em>, assuming ordering of 
        <em>child</em> as given by get_connectivity called on <em>child</em>; sense is 1, -1
        for forward/reverse sense, resp.
        \param offset Offset between first vertex of <em>child</em> and first vertex of side 
        <em>side_number</em> on <em>parent</em>
    */
  virtual ErrorCode side_number(const EntityHandle parent,
                                  const EntityHandle child,
                                  int &side_number,
                                  int &sense,
                                  int &offset) const = 0;

    //! Find the higher-order node on a subfacet of an entity
    /** Given an entity and the connectivity and type of one of its subfacets, find the
        high order node on that subfacet, if any.  The number of vertices in <em>subfacet_conn</em>
        is derived from <em>subfacet_type</em> and the canonical numbering for that type.
        \param parent_handle The element whose subfacet is being queried
        \param subfacet_conn The connectivity of the subfacet being queried
        \param subfacet_type The type of subfacet being queried
        \param high_order_node If the subfacet has a high-order node defined on <em>parent_handle</em>,
        the handle for that node.
    */
  virtual ErrorCode high_order_node(const EntityHandle parent_handle,
                                      const EntityHandle *subfacet_conn,
                                      const EntityType subfacet_type,
                                      EntityHandle &high_order_node) const = 0;

    //! Return the handle of the side element of a given dimension and index
    /** Given a parent entity and a target dimension and side number, return the handle of
        the entity corresponding to that side.  If an entity has not been created to represent
        that side, one is not created by this function, and zero is returned in <em>target_entity</em>.
        \param source_entity The entity whose side is being queried.
        \param dim The topological dimension of the side being queried.
        \param side_number The canonical index of the side being queried.
        \param target_entity The handle of the entity representing this side, if any.
    */
  virtual ErrorCode side_element(const EntityHandle source_entity,
                                   const int dim, 
                                   const int side_number,
                                   EntityHandle &target_entity) const = 0;

    /**@}*/

    /** \name Tags */

    /**@{*/

    /**\brief Get a tag handle, possibly creating the tag
     *
     * Get a handle used to associate application-defined values
     * with MOAB entities.  If the tag does not already exist then
     * \c flags should contain exactly one of \c MB_TAG_SPARSE, 
     * \c MB_TAG_DENSE, \c MB_TAG_MESH unless \c type is MB_TYPE_BIT,
     * which implies \c MB_TAG_BIT storage.  
     * .
     *\param name          The tag name
     *\param size          Tag size as number of values of of data type per entity
     *                     (or number of bytes if \c MB_TAG_BYTES is passed in flags).  If \c MB_TAG_VARLEN
     *                     is specified, this value is taken to be the size of the
     *                     default value if one is specified and is otherwise ignored.
     *\param type          The type of the data (used for IO)
     *\param tag_handle    Output: the resulting tag handle.
     *\param flags         Bitwise OR of values from TagType
     *\param default_value Optional default value for tag.
     *\param created       Optional returned boolean indicating that the tag
     *                     was created.
     *\return - \c MB_ALREADY_ALLOCATED     if tag exists and \c MB_TAG_EXCL is specified, or default values
     *                                      do not match (and \c MB_TAG_ANY or \c MB_TAG_DFTOK not specified).
     *        - \c MB_TAG_NOT_FOUND         if tag does not exist and \c MB_TAG_CREAT is not specified
     *        - \c MB_INVALID_SIZE          if tag value size is not a multiple of the size of the data type
     *                                      (and \c MB_TAG_ANY not specified).
     *        - \c MB_TYPE_OUT_OF_RANGE     invalid or inconsistent parameter
     *        - \c MB_VARIABLE_DATA_LENGTH  if \c MB_TAG_VARLEN and \c default_value is non-null and
     *                                      \c default_value_size is not specified.
     *
     *\NOTE A call to tag_get_handle that includes a default value will fail
     * if the tag already exists with a different default value.  A call without
     * a default value will succeed if the tag already exists, regardless of 
     * whether or not the existing tag has a default value.
     *
     * Examples:
     *
     * Retrieve a handle for an existing tag, returning a non-success error
     * code if the tag does not exist or does not store 1 integer value per
     * entity:
     *\code
     * Tag git_tag;
     * mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag );
     * \endcode
     * Get the tag handle, or create it as a dense tag if it does not already 
     * exist:
     *\code
     * Tag gid_tag;
     * mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_CREAT|MB_TAG_BIT );
     * \endcode
     * Create the tag or *fail* if it already exists:
     *\code
     * Tag gid_tag;
     * mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_EXCL|MB_TAG_DENSE );
     * \endcode
     * Get an existing variable length tag, failing if it does not exist or
     * is not variable-length or does not contain double values.
     *\code
     * Tag vtag;
     * mb.tag_get_handle( tag_name, 0, MB_TYPE_DOUBLE, vtag );
     * \endcode
     * Get the same variable-length tag, but create it with a default value
     * if it doesn't exist.  Note that if the tag already exists this call
     * will return a non-success error code if the existing tag has a different
     * default value.
     *\code
     * Tag vtag;
     * const double default_val = M_PI;
     * const int def_val_len = 1;
     * mb.tag_get_handle( tag_name, def_val_len, MB_TYPE_DOUBLE, vtag,
     *                    MB_TAG_SPARSE|MB_TAG_VARLEN|MB_TAG_CREAT, &default_val );
     * \endcode
     */
  virtual ErrorCode tag_get_handle( const char* name,
                                    int size,
                                    DataType type,
                                    Tag& tag_handle,
                                    unsigned flags = 0,
                                    const void* default_value = 0,
                                    bool* created = 0 ) = 0;
  
    /**\brief same as non-const version, except that MB_TAG_CREAT flag is ignored. */
  virtual ErrorCode tag_get_handle( const char* name,
                                    int size,
                                    DataType type,
                                    Tag& tag_handle,
                                    unsigned flags = 0,
                                    const void* default_value = 0 ) const = 0;

    //! Get the name of a tag corresponding to a handle
    /** \param tag_handle Tag you want the name of.  
        \param tag_name Name string for <em>tag_handle</em>. 
    */
  virtual ErrorCode  tag_get_name(const Tag tag_handle, 
                                    std::string& tag_name) const = 0;

    /**\brief Gets the tag handle corresponding to a name
       
        If a tag of that name does not exist, returns MB_TAG_NOT_FOUND
        \param tag_name Name of the desired tag. 
        \param tag_handle Tag handle corresponding to <em>tag_name</em>
    */ 
  virtual ErrorCode tag_get_handle( const char *tag_name, 
                                    Tag &tag_handle ) const = 0;

    //! Get the size of the specified tag in bytes
    /** Get the size of the specified tag, in bytes (MB_TAG_SPARSE, MB_TAG_DENSE, MB_TAG_MESH)
        \note always returns 1 for bit tags.
        \param tag Handle of the desired tag. 
        \param bytes_per_tag Size of the specified tag
        \return - MB_TAG_NOT_FOUND for invalid tag handles
                - MB_VARIABLE_DATA_LENGTH for variable-length tags
                - MB_SUCCESS otherwise
    */ 
  virtual ErrorCode tag_get_bytes(const Tag tag, int& bytes_per_tag) const = 0;

    //! Get the array length of a tag
    /** Get the size of the specified tag, in as the number of values of
        the basic type (e.g. number of integer values for each tag value if
        the data type is MB_TYPE_INTEGER).  Gives number of bits for bit tags
        and is the same as \c tag_get_bytes for MB_TYPE_OPAQUE tags.
        \param tag Handle of the desired tag. 
        \param length Size of the specified tag
        \return - MB_TAG_NOT_FOUND for invalid tag handles
                - MB_VARIABLE_DATA_LENGTH for variable-length tags
                - MB_SUCCESS otherwise
    */ 
  virtual ErrorCode tag_get_length(const Tag tag, int &length) const = 0;

    //! Get the type of the specified tag
    /** Get the type of the specified tag
        \param tag Handle of the desired tag. 
        \param tag_type Type of the specified tag
    */ 
  virtual ErrorCode tag_get_type(const Tag tag, TagType &tag_type) const = 0;

    /** \brief Get data type of tag.
     *
     * Get the type of the tag data.  The tag is data is assumed to
     * be a vector of this type.  If the tag data vetcor contains 
     * more than one value, then the tag size must be a multiple of
     * the size of this type.
     * \param tag  The tag 
     * \param type The type of the specified tag (output).
     */
   virtual ErrorCode tag_get_data_type(const Tag tag, DataType& type) const = 0;

    //! Get the default value of the specified tag
    /** Get the default value of the specified tag
        \param tag Handle of the desired tag. 
        \param def_value Pointer to memory where default value of the specified tag is written
        \return - MB_ENTITY_NOT_FOUND If no default value is set for tag.
                - MB_SUCCESS          If success.
                - MB_FAILURE          If <code>def_val</code> is NULL.
                - MB_TAG_NOT_FOUND    If <code>tag_handle</code> is invalid.
    */ 
  virtual ErrorCode tag_get_default_value(const Tag tag, void *def_val) const = 0;
  virtual ErrorCode tag_get_default_value( Tag tag, const void*& def_val, int& size) const = 0;

    //! Get handles for all tags defined in the mesh instance
    /** Get handles for all tags defined on the mesh instance.
        \param tag_handles STL vector of all tags
    */
  virtual ErrorCode tag_get_tags(std::vector<Tag> &tag_handles) const = 0;

    //! Get handles for all tags defined on this entity
    /** Get handles for all tags defined on this entity; if zero, get all tags defined 
        on mesh instance
        \param entity Entity for which you want tags
        \param tag_handles STL vector of all tags defined on <em>entity</em>
    */
  virtual ErrorCode tag_get_tags_on_entity(const EntityHandle entity,
                                             std::vector<Tag> &tag_handles) const = 0;

    //! Get the value of the indicated tag on the specified entities in the specified vector
    /** Get the value of the indicated tag on the specified entities; <em>tag_data</em> must contain
        enough space (i.e. tag_size*num_entities bytes) to hold all tag data.  MOAB does <em>not</em>
        check whether this space is available before writing to it.
        \note For bit tags, tag_data must contain one byte per entity.  For each
              entity, the corresponding byte will contain the tag bits in the
              lower bit positions and zero bits in the higher.
        \param tag_handle Tag whose values are being queried.
        \param entity_handles 1d vector of entity handles whose tag values are being queried
        \param num_entities Number of entities in 1d vector of entity handles
        \param tag_data Pointer to memory into which tag data will be written
    */
  virtual ErrorCode  tag_get_data(const Tag tag_handle, 
                                  const EntityHandle* entity_handles, 
                                  int num_entities, 
                                  void *tag_data) const = 0;

    //! Get the value of the indicated tag on the specified entities in the specified range
    /** Identical to previous function, except entities are specified using a range instead of a 1d vector.
        \param tag_handle Tag whose values are being queried.
        \param entity_handles Range of entity handles whose tag values are being queried
        \param tag_data Pointer to memory into which tag data will be written
    */
  virtual ErrorCode  tag_get_data(const Tag tag_handle, 
                                    const Range& entity_handles, 
                                    void *tag_data) const = 0;

    //! Set the value of the indicated tag on the specified entities in the specified vector
    /** Set the value of the indicated tag on the specified entities; <em>tag_data</em> contains the
        values, <em>one value per entity in <em>entity_handles</em></em>.
        \note For bit tags, tag_data must contain one byte per entity.  For each
              entity, the tag bits will be read from the lower bits of the 
              corresponding byte.
        \param tag_handle Tag whose values are being set
        \param entity_handles 1d vector of entity handles whose tag values are being set
        \param num_entities Number of entities in 1d vector of entity handles
        \param tag_data Pointer to memory holding tag values to be set, <em>one entry per entity handle</em>
    */
  virtual ErrorCode  tag_set_data( Tag tag_handle, 
                                   const EntityHandle* entity_handles, 
                                   int num_entities,
                                   const void *tag_data ) = 0;
  
    //! Set the value of the indicated tag on the specified entities in the specified range
    /** Identical to previous function, except entities are specified using a range instead of a 1d vector.
        \param tag_handle Tag whose values are being set
        \param entity_handles Range of entity handles whose tag values are being set
        \param tag_data Pointer to memory holding tag values to be set, <em>one entry per entity handle</em>
    */
  virtual ErrorCode  tag_set_data( Tag tag_handle, 
                                    const Range& entity_handles,
                                    const void *tag_data ) = 0;

    /**\brief Get pointers to tag data
     *
     * For a tag, get the values for a list of passed entity handles.
     *\note  This function may not be used for bit tags.
     *\param tag_handle     The tag
     *\param entity_handles An array of entity handles for which to retreive tag values.
     *\param num_entities   The length of the 'entity_handles' array.
     *\param tag_data       An array of 'const void*'.  Array must be at least
     *                      'num_entitities' long.  Array is populated (output)
     *                      with pointers to the internal storage for the
     *                      tag value corresponding to each entity handle.
     *\param tag_sizes      The length of each tag value.  Optional for 
     *                      fixed-length tags.  Required for variable-length tags.
     */
  virtual ErrorCode  tag_get_by_ptr(const Tag tag_handle, 
                                  const EntityHandle* entity_handles, 
                                  int num_entities, 
                                  const void** tag_data,
                                  int* tag_sizes = 0 ) const = 0;

    /**\brief Get pointers to tag data
     *
     * For a tag, get the values for a list of passed entity handles.
     *\note  This function may not be used for bit tags.
     *\param tag_handle     The tag
     *\param entity_handles The entity handles for which to retreive tag values.
     *\param tag_data       An array of 'const void*'.  Array is populated (output)
     *                      with pointers to the internal storage for the
     *                      tag value corresponding to each entity handle.
     *\param tag_sizes      The length of each tag value.  Optional for 
     *                      fixed-length tags.  Required for variable-length tags.
     */
  virtual ErrorCode  tag_get_by_ptr(const Tag tag_handle, 
                                    const Range& entity_handles, 
                                    const void** tag_data,
                                    int* tag_sizes = 0 ) const = 0;

    /**\brief Set tag data given an array of pointers to tag values.
     *
     * For a tag, set the values for a list of passed entity handles.
     *\note  This function may not be used for bit tags.
     *\param tag_handle     The tag
     *\param entity_handles An array of entity handles for which to set tag values.
     *\param num_entities   The length of the 'entity_handles' array.
     *\param tag_data       An array of 'const void*'.  Array must be at least
     *                      'num_entitities' long.  Array is expected to
     *                      contain pointers to tag values for the corresponding
     *                      EntityHandle in 'entity_handles'.
     *\param tag_sizes      The length of each tag value.  Optional for 
     *                      fixed-length tags.  Required for variable-length tags.
     */
  virtual ErrorCode  tag_set_by_ptr( Tag tag_handle, 
                                   const EntityHandle* entity_handles, 
                                   int num_entities,
                                   void const* const* tag_data,
                                   const int* tag_sizes = 0 ) = 0;
  
    /**\brief Set tag data given an array of pointers to tag values.
     *
     * For a tag, set the values for a list of passed entity handles.
     *\note  This function may not be used for bit tags.
     *\param tag_handle     The tag
     *\param entity_handles The entity handles for which to set tag values.
     *\param tag_data       An array of 'const void*'.  Array is expected to
     *                      contain pointers to tag values for the corresponding
     *                      EntityHandle in 'entity_handles'.
     *\param tag_sizes      The length of each tag value.  Optional for 
     *                      fixed-length tags.  Required for variable-length tags.
     */
  virtual ErrorCode  tag_set_by_ptr( Tag tag_handle, 
                                    const Range& entity_handles,
                                    void const* const* tag_data,
                                    const int* tag_sizes = 0 ) = 0;

    /**\brief Set tag data given value.
     *
     * For a tag, set the values for a list of passed entity handles to
     * the same, specified value.
     *
     *\param tag_handle     The tag
     *\param entity_handles The entity handles for which to set tag values.
     *\param tag_data       A pointer to the tag value.
     *\param tag_sizes      For variable-length tags, the length of the
     *                      tag value.  This argument will be ignored for
     *                      fixed-length tags.
     */
  virtual ErrorCode tag_clear_data( Tag tag_handle,
                                    const Range& entity_handles,
                                    const void* value,
                                    int value_size = 0 ) = 0;

    /**\brief Set tag data given value.
     *
     * For a tag, set the values for a list of passed entity handles to
     * the same, specified value.
     *
     *\param tag_handle     The tag
     *\param entity_handles The entity handles for which to set tag values.
     *\param tag_data       A pointer to the tag value.
     *\param tag_sizes      For variable-length tags, the length of the
     *                      tag value.  This argument will be ignored for
     *                      fixed-length tags.
     */
  virtual ErrorCode tag_clear_data( Tag tag_handle,
                                    const EntityHandle* entity_handles,
                                    int num_entity_handles,
                                    const void* value,
                                    int value_size = 0 ) = 0;

    //! Delete the data of a vector of entity handles and sparse tag
    /** Delete the data of a tag on a vector of entity handles.  Only sparse tag data are deleted with this
        function; dense tags are deleted by deleting the tag itself using tag_delete.
        \param tag_handle Handle of the (sparse) tag being deleted from entity
        \param entity_handles 1d vector of entity handles from which the tag is being deleted
        \param num_handles Number of entity handles in 1d vector
    */
  virtual ErrorCode  tag_delete_data( Tag tag_handle, 
                                      const EntityHandle *entity_handles,
                                      int num_handles) = 0;

    //! Delete the data of a range of entity handles and sparse tag
    /** Delete the data of a tag on a range of entity handles.  Only sparse tag data are deleted with this
        function; dense tags are deleted by deleting the tag itself using tag_delete.
        \param tag_handle Handle of the (sparse) tag being deleted from entity
        \param entity_range Range of entities from which the tag is being deleted
    */
  virtual ErrorCode  tag_delete_data( Tag tag_handle, 
                                      const Range &entity_range) = 0;

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
   *\param tag_handle  The handle of the tag for which to access data
   *\param iter        The first entity for which to return data. 
   *\param end         One past the last entity for which data is desired.
   *\param count       The number of entities for which data was returned
   *\param data_ptr    Output: pointer to tag storage.
   *\param allocate    If true, space for this tag will be allocated, if not it wont
   *  
   *\Note If this function is called for entities for which no tag value
   *      has been set, but for which a default value exists, it will 
   *      force the allocation of explicit storage for each such entity
   *      even though MOAB would normally not explicitly store tag values
   *      for such entities.
   *
   *\Example:
   *\code
   * Range ents; // range to iterate over
   * Tag tag; // tag for which to access data
   * int bytes;
   * ErrorCode err = mb.tag_get_size( tag, bytes );
   * if (err) { ... }
   * 
   * ...
   * Range::iterator iter = ents.begin();
   * while (iter != ents.end()) {
   *   int count;
   *    // get contiguous block of tag dat
   *   void* ptr;
   *   err = mb.tag_iterate( tag, iter, ents.end(), count, ptr );
   *   if (err) { ... }
   *    // do something with tag data
   *   process_Data( ptr, count );
   *    // advance to next block of data
   *   iter += count;
   * }
   *\endcode
   */
  virtual ErrorCode tag_iterate( Tag tag_handle,
                                 Range::const_iterator begin,
                                 Range::const_iterator end,
                                 int& count,
                                 void*& data_ptr,
                                 bool allocate = true) = 0;

    //! Remove a tag from the database and delete all of its associated data
    /** Deletes a tag and all associated data.
     */
  virtual ErrorCode  tag_delete(Tag tag_handle) = 0;

    /**@}*/

    /** \name Sets */

    /**@{*/

    //! Create a new mesh set
    /** Create a new mesh set.  Meshsets can store entities ordered or unordered. A set can include entities
        at most once (MESHSET_SET) or more than once.  Meshsets can optionally track its members using
        adjacencies (MESHSET_TRACK_OWNER); if set, entities are deleted from tracking meshsets before
        being deleted.  This adds data to mesh entities, which can be expensive.
        \param options Options bitmask for the new meshset, possible values defined above
        \param ms_handle Handle for the meshset created
    */
  virtual ErrorCode create_meshset(const unsigned int options, 
                                     EntityHandle &ms_handle,
                                     int start_id = 0) = 0;

    //! Empty a vector of mesh set
    /** Empty a mesh set.
        \param ms_handles 1d vector of handles of sets being emptied
        \param num_meshsets Number of entities in 1d vector
    */
  virtual ErrorCode clear_meshset( const EntityHandle *ms_handles, 
                                     const int num_meshsets) = 0;

    //! Empty a range of mesh set
    /** Empty a mesh set.
        \param ms_handles Range of handles of sets being emptied
    */
  virtual ErrorCode clear_meshset( const Range &ms_handles) = 0;

    //! Get the options of a mesh set
    /** Get the options of a mesh set.
        \param ms_handle Handle for mesh set being queried
        \param options Bit mask in which mesh set options are returned
    */
  virtual ErrorCode get_meshset_options(const EntityHandle ms_handle, 
                                          unsigned int& options) const = 0;

    //! Set the options of a mesh set
    /** Set the options of a mesh set.
        \param ms_handle Handle for meshset whose options are being changed
        \param options Bit mask of options to be used
    */
  virtual ErrorCode set_meshset_options(const EntityHandle ms_handle, 
                                          const unsigned int options) = 0;

    //! Subtract meshsets
    /** Subtract <em>meshset2</em> from <em>meshset1</em>, placing the results in meshset1.
        \param meshset1 Mesh set being subtracted from, also used to pass back result
        \param meshset2 Mesh set being subtracted from <em>meshset1</em>
    */
  virtual ErrorCode subtract_meshset(EntityHandle meshset1, 
                                       const EntityHandle meshset2) = 0;

    //! Intersect meshsets
    /** Intersect <em>meshset1</em> with <em>meshset2</em>, placing the results in meshset1.
        \param meshset1 Mesh set being intersected, also used to pass back result
        \param meshset2 Mesh set being intersected with <em>meshset1</em>
    */
  virtual ErrorCode intersect_meshset(EntityHandle meshset1, 
                                        const EntityHandle meshset2) = 0;
    
    //! Unite meshsets
    /** Unite <em>meshset1</em> with <em>meshset2</em>, placing the results in meshset1.
        \param meshset1 Mesh set being united, also used to pass back result
        \param meshset2 Mesh set being united with <em>meshset1</em>
    */
  virtual ErrorCode unite_meshset(EntityHandle meshset1, 
                                    const EntityHandle meshset2) = 0;

    //! Add to a meshset entities in specified range
    /** Add to a meshset entities in specified range.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies are also added to entities in <em>entities</em>.
        \param meshset Mesh set being added to
        \param entities Range of entities being added to meshset
    */
  virtual ErrorCode add_entities(EntityHandle meshset, 
                                   const Range &entities) = 0;

    //! Add to a meshset entities in specified vector
    /** Add to a meshset entities in specified vector.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies are also added to entities in <em>entities</em>.
        \param meshset Mesh set being added to
        \param entities 1d vector of entities being added to meshset
        \param num_entities Number of entities in 1d vector
    */
  virtual ErrorCode add_entities(EntityHandle meshset, 
                                   const EntityHandle *entities,
                                   const int num_entities) = 0;
  
    //! Remove from a meshset entities in specified range
    /** Remove from a meshset entities in specified range.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies in entities in <em>entities</em> are updated.
        \param meshset Mesh set being removed from
        \param entities Range of entities being removed from meshset
    */
  virtual ErrorCode remove_entities(EntityHandle meshset, 
                                      const Range &entities) = 0;

    //! Remove from a meshset entities in specified vector
    /** Remove from a meshset entities in specified vector.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies in entities in <em>entities</em> are updated.
        \param meshset Mesh set being removed from
        \param entities 1d vector of entities being removed from meshset
        \param num_entities Number of entities in 1d vector
    */
  virtual ErrorCode remove_entities(EntityHandle meshset, 
                                      const EntityHandle *entities,
                                      const int num_entities) = 0;

    //! Return whether a set contains entities
    /** Return whether a set contains entities.  Returns true only
     * if ALL entities are contained
     * \param meshset Mesh set being queried
     * \param entities Entities being queried
     * \param num_entities Number of entities
     * \return bool If true, all entities are contained in set
    */
  virtual bool contains_entities(EntityHandle meshset, 
                                 const EntityHandle *entities,
                                 int num_entities,
                                 const int operation_type = Interface::INTERSECT) = 0;

    //! Replace entities in a set with other entities
    /** Replace entities in a set with other entities
     * 
     * \note  Behavior is undefined if an entity handle exists in both the
     *        old_entities and the new_entities arrays or old_entities
     *        contains multiple copies of an entity.
     * \note  If an entity occurs multiple times in an ordered set, all
     *        occurances will be replaced.
     * \note  For list-based sets, if not all handles in old_entities 
     *        occcur in the set, the corresponding new_entities will not 
     *        be added and  MB_ENTITY_NOT_FOUND will be returned.
     *        For set-based sets, all entities in new_entities wll be
     *        added and any contained entities in old_entities will be 
     *        removed, and the return value will be MB_SUCCESS.
     * \param meshset Mesh set being modified
     * \param old_entities Entities to replace
     * \param new_entities New entities to add
     * \param num_entities Number of entities in input arrays
     * \return - MB_SUCCESS : all entities in old_entities replaced
     *         - MB_ENTITY_NOT_FOUND : one or more entities in new_entities 
     *                        not added to set because corresponding entity
     *                        in old_entities did not occur in the ordered
     *                        set.
     *         - MB_FAILURE : internal error
    */
  virtual ErrorCode replace_entities(EntityHandle meshset, 
                                       const EntityHandle *old_entities,
                                       const EntityHandle *new_entities,
                                       int num_entities) = 0;
    /**@}*/

    /** \name Set parents/children */

    /**@{*/
  
    //! Get parent mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate parents are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose parents are being queried
        \param parents STL vector holding the parents returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual ErrorCode get_parent_meshsets(const EntityHandle meshset,
                                          std::vector<EntityHandle> &parents, 
                                          const int num_hops = 1) const = 0;

    //! Get parent mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate parents are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose parents are being queried
        \param parents Range holding the parents returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual ErrorCode get_parent_meshsets(const EntityHandle meshset,
                                          Range &parents,
                                          const int num_hops = 1) const = 0;

    //! Get child mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate children are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose children are being queried
        \param children STL vector holding the children returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual ErrorCode get_child_meshsets(const EntityHandle meshset, 
                                         std::vector<EntityHandle> &children, 
                                         const int num_hops = 1) const = 0;

    //! Get child mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate children are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose children are being queried
        \param children Range holding the children returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual ErrorCode get_child_meshsets(const EntityHandle meshset, 
                                         Range &children, 
                                         const int num_hops = 1) const = 0;

    /**\brief Get mesh sets contained in a mesh set
     *
     * If <em>num_hops</em> is 1, only immediate contents are returned. 
     * Otherwise a recursive query of all contained sets is performed, 
     * returning every visted set.  The value of num_hops limites the
     * depth of the search, with zero indicating no depth limit. 
     *
     *\param meshset The mesh set whose contents are being queried
     *\param contained The result list.
     *\param num_hops Number of generations to traverse (0 = all)
     */
  virtual ErrorCode get_contained_meshsets(const EntityHandle meshset, 
                                         std::vector<EntityHandle> &contained, 
                                         const int num_hops = 1) const = 0;

    /**\brief Get mesh sets contained in a mesh set
     *
     * If <em>num_hops</em> is 1, only immediate contents are returned. 
     * Otherwise a recursive query of all contained sets is performed, 
     * returning every visted set.  The value of num_hops limites the
     * depth of the search, with zero indicating no depth limit. 
     *
     *\param meshset The mesh set whose contents are being queried
     *\param contained The result list.
     *\param num_hops Number of generations to traverse (0 = all)
     */
  virtual ErrorCode get_contained_meshsets(const EntityHandle meshset, 
                                             Range &contained, 
                                             const int num_hops = 1) const = 0;

    //! Get the number of parent mesh sets of a mesh set
    /** Identical to get_parent_meshsets, only number is returned instead of actual parents.
        \param meshset The mesh set whose parents are being queried
        \param number Number of parents
    */
  virtual ErrorCode num_parent_meshsets(const EntityHandle meshset,  
                                          int *number,
                                          const int num_hops = 1) const = 0;

    //! Get the number of child mesh sets of a mesh set
    /** Identical to get_child_meshsets, only number is returned instead of actual children.
        \param meshset The mesh set whose children are being queried
        \param number Number of children
    */
  virtual ErrorCode num_child_meshsets(const EntityHandle meshset, 
                                         int *number,
                                         const int num_hops = 1) const = 0;

    /**\brief Get the number of mesh sets contained in a mesh set
     *
     * Return the number of sets that would be returned by get_contained_meshsets
     *
     *\param meshset  The initial set to begin the query from.
     *\param number   (Output) The result count.
     *\param num_hops Search depth (0 => unbounded).
     */
  virtual ErrorCode num_contained_meshsets(const EntityHandle meshset, 
                                             int *number,
                                             const int num_hops = 1) const = 0;

    //! Add a parent mesh set to a mesh set
    /** Make <em>parent_meshset</em> a new parent of <em>child_meshset</em>.  This function does 
        <em>not</em> add a corresponding child link to <em>parent_meshset</em>.
        \param child_meshset The child mesh set being given a new parent.
        \param parent_meshset The parent being added to <em>child_meshset</em>
    */
  virtual ErrorCode add_parent_meshset(EntityHandle child_meshset, 
                                         const EntityHandle parent_meshset) = 0;

    //! Add a parent mesh sets to a mesh set
    /** Make <em>parent_meshset</em> a new parent of <em>child_meshset</em>.  This function does 
        <em>not</em> add a corresponding child link to <em>parent_meshset</em>.
        \param child_meshset The child mesh set being given a new parent.
        \param parent_meshset The parent being added to <em>child_meshset</em>
    */
  virtual ErrorCode add_parent_meshsets(EntityHandle child_meshset, 
                                          const EntityHandle* parent_meshsets,
                                          int num_parent_meshsets ) = 0;

    //! Add a child mesh set to a mesh set
    /** Make <em>child_meshset</em> a new child of <em>parent_meshset</em>.  This function does 
        <em>not</em> add a corresponding parent link to <em>child_meshset</em>.
        \param parent_meshset The parent mesh set being given a new child.
        \param child_meshset The child being added to <em>parent_meshset</em>
    */
  virtual ErrorCode add_child_meshset(EntityHandle parent_meshset, 
                                        const EntityHandle child_meshset) = 0;

    //! Add a child mesh sets to a mesh set
    /** Make <em>child_meshset</em> a new child of <em>parent_meshset</em>.  This function does 
        <em>not</em> add a corresponding parent link to <em>child_meshset</em>.
        \param parent_meshset The parent mesh set being given a new child.
        \param child_meshset The child being added to <em>parent_meshset</em>
    */
  virtual ErrorCode add_child_meshsets(EntityHandle parent_meshset, 
                                         const EntityHandle* child_meshsets,
                                         int num_child_meshsets ) = 0;

    //! Add parent and child links between mesh sets
    /** Makes <em>child_meshset</em> a new child of <em>parent_meshset</em>, and vica versa.
        \param parent The parent mesh set being given a new child, and the new parent
        \param child The child being given a new parent, and the new child
    */
  virtual ErrorCode add_parent_child( EntityHandle parent, 
                                        EntityHandle child ) = 0;

    //! Remove parent and child links between mesh sets
    /** Removes parent/child links between <em>child_meshset</em> and <em>parent_meshset</em>.
        \param parent The parent mesh set being removed from <em>child</em>
        \param child The child mesh set being removed from <em>parent</em>
    */
  virtual ErrorCode remove_parent_child( EntityHandle parent, 
                                           EntityHandle child ) = 0;

    //! Remove a parent mesh set from a mesh set
    /** Removes <em>parent_meshset</em> from the parents of <em>child_meshset</em>.  This function does 
        <em>not</em> remove a corresponding child link from <em>parent_meshset</em>.
        \param child_meshset The child mesh whose parent is being removed
        \param parent_meshset The parent being removed from <em>meshset</em>
    */
  virtual ErrorCode remove_parent_meshset(EntityHandle child_meshset, 
                                            const EntityHandle parent_meshset) = 0;
  
    //! Remove a child mesh set from a mesh set
    /** Removes <em>child_meshset</em> from the children of <em>parent_meshset</em>.  This function does 
        <em>not</em> remove a corresponding parent link from <em>child_meshset</em>.
        \param parent_meshset The parent mesh set whose child is being removed
        \param child_meshset The child being removed from <em>parent_meshset</em>
    */
  virtual ErrorCode remove_child_meshset(EntityHandle parent_meshset, 
                                           const EntityHandle child_meshset) = 0;

    /**@}*/

    /** \name Set iterators */

    /**@{*/

    /** \brief Create an iterator over the set
     * Create a new iterator that iterates over entities with the specified type or dimension.  
     * Only one of ent_type or dim can be set; use dim=-1 or ent_type=MBMAXTYPE for the other.
     * Iterators for list-type (ordered) sets are stable over set modification, unless entity
     * removed or deleted is the one at the current position of the iterator.  If the check_valid
     * parameter is passed as true, entities are checked for validity before being passed back by
     * get_next_entities function (checking entity validity can have a non-negligible cost).
     *
     * Iterators returned by this function can be deleted using the normal C++ delete function.
     * After creating the iterator through this function, further interactions are through methods
     * on the SetIterator class.
     * \param meshset The entity set associated with this iterator (use 0 for whole instance)
     * \param ent_type Entity type associated with this iterator
     * \param ent_dim Dimension associated with this iterator
     * \param chunk_size Chunk size of the iterator
     * \param check_valid If true, entities are checked for validity before being returned 
     */
  virtual ErrorCode create_set_iterator(EntityHandle meshset,
                                        EntityType ent_type,
                                        int ent_dim,
                                        int chunk_size,
                                        bool check_valid,
                                        SetIterator *&set_iter) = 0;
    /**@}*/

};

//! predicate for STL algorithms.  Returns true if the entity handle is
//! of the specified type.  For example, to remove all the tris out of a list
//! of 2D entities retrieved using get_adjacencies you could do
//! std::remove_if(list.begin(), list.end(), type_equals(gMB, MeshTri));
class type_equals : public std::unary_function<EntityHandle, bool>
{
public:
    //! interface object
  Interface* meshDB;

    //! type corresponding to this predicate
  const EntityType test_type;

    //! Constructor
  type_equals(Interface* mdb, const EntityType type) : meshDB(mdb), test_type(type){}

    //! operator predicate
  bool operator()(EntityHandle handle) const
    { 
      return (meshDB->type_from_handle(handle) ==  test_type); 
    } 
};

//! predicate for STL algorithms.  Returns true if the entity handle is not
//! of the specified type.  For example, to remove all but the tris out of a list
//! of 2D entities retrieved using get_adjacencies you could do
//! std::remove_if(list.begin(), list.end(), type_not_equals(gMB, MeshTri));
class type_not_equals : public std::unary_function<EntityHandle, bool>
{
public:

    //! interface object
  Interface* meshDB;

    //! type corresponding to this predicate
  const EntityType test_type;

    //! Constructor
  type_not_equals(Interface* mdb, const EntityType type) : meshDB(mdb), test_type(type){}

    //! operator predicate
  bool operator()(EntityHandle handle) const
    { 
      return (meshDB->type_from_handle(handle) !=  test_type); 
    } 
};

inline float Interface::api_version(std::string *version_string) 
{
  if (NULL != version_string)
    *version_string = std::string("MOAB API version ") + std::string(MOAB_API_VERSION_STRING);
  return MOAB_API_VERSION;
}

} // namespace moab 

#endif   // MB_INTERFACE_HPP

  
