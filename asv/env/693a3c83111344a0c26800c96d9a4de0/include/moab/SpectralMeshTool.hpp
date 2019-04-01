#ifndef SPECTRALMESHTOOL_HPP
#define SPECTRALMESHTOOL_HPP

#include "moab/Interface.hpp" // needs to be here to support inline query_interface
#include "moab/Error.hpp" // needs to be here to support inline query_interface
#include <vector>

namespace moab {

/** \class SpectralMeshTool
 * \brief Class with convenience functions for handling spectral mesh
 * Class with convenience functions for handling spectral meshes.  See description of spectral
 * mesh handling in doc/metadata_info.doc and in the MOAB user's guide.
 *
 * There are two primary representations of spectral meshes:
 * a) coarse elements: with SPECTRAL_VERTICES lexicographically-ordered array of fine vertices 
 *    on each element, and tags on vertices or on elements (with _LEX suffix)
 * b) fine elements: as linear elements made from fine vertices, with tags on vertices
 * 
 */
class SpectralMeshTool
{
public:

    /** \brief Constructor
     * \param impl MOAB Interface instance
     * \param order Spectral order, defaults to 0
     */
  SpectralMeshTool(Interface *impl, int order = 0);

    /** \brief Destructor
     */
  ~SpectralMeshTool();

    /** \brief Return tag used to store lexicographically-ordered vertex array
     * NOTE: If creating this tag with this call, this SpectralMeshTool instance must already have
     * a non-zero spectral order value set on it; the size of the spectral vertices tag depends on this order.
     * \param sv_tag Spectral vertices tag
     * \param create_if_missing If true, will create this tag if it doesn't exist already
     */
  Tag spectral_vertices_tag(const bool create_if_missing = false);
  
    /** \brief Return tag used to store spectral order
     * \param so_tag Spectral order tag
     * \param create_if_missing If true, will create this tag if it doesn't exist already
     */
  Tag spectral_order_tag(const bool create_if_missing = false);
  
    /** \brief Convert representation from coarse to fine
     * Each element in set, or in interface if set is not input, is converted to fine elements, using
     * vertices in SPECTRAL_VERTICES tagged array
     * \param spectral_set Set containing spectral elements
     */
  ErrorCode convert_to_fine(EntityHandle spectral_set);
  
    /** \brief Convert representation from fine to coarse
     * Each element in set, or in interface if set is not input, is converted to coarse elements, with
     * fine vertices put into SPECTRAL_VERTICES tagged array.  NOTE: This function assumes that each
     * order^d (fine) elements comprise each coarse element, and are in order of fine elements in each
     * coarse element.  If order is input as 0, looks for a SPECTRAL_ORDER tag on the mesh.
     * \param order Order of the spectral mesh
     * \param spectral_set Set containing spectral elements
     */
  ErrorCode convert_to_coarse(int order = 0, EntityHandle spectral_set = 0);

    /** \brief Create coarse spectral elements from fine elements pointed to by conn
     * This function creates the coarse elements by taking conn (assumed to be in FE ordering)
     * and picking out the corner vertices to make coarse connectivity, and the other vertices
     * (along with corners) to make SPECTRAL_VERTICES array pointed to by each entity.
     * \param conn Connectivity of fine (linear) elements, in FE ordering
     * \param verts_per_e Vertices per entity
     * \param num_fine_elems Number of fine elements represented by conn
     * \param spectral_set Set to which coarse elements should be added, if any
     * \param start_idx Starting index in conn (for parallel support)
     * \param local_gids If non-null, will insert all fine vertices into this range
     */
  template <class T>
  ErrorCode create_spectral_elems(const T *conn, int num_fine_elems, int dim, 
                                  Range &output_range, int start_idx = 0, Range *local_gids = NULL);
  
    /** \brief Set spectral order for this instance
     * \param order Order set on this instance
     */
  void spectral_order(int order) {spectralOrder = order; spectralOrderp1 = order+1;}
  
    /** \brief Get spectral order for this instance
     * \return order Order set on this instance
     */
  int spectral_order() {return spectralOrder;}
/*
  struct ConnMap 
  {
    const short a[16];
  };
  */
  static const short int permute_array[];
  
  static const short int lin_permute_array[];
  
private:

  //! the MB instance that this works with
  Interface* mbImpl;

    //! error object for this tool
  Error *mError;
  
    //! SPECTRAL_VERTICES tag
  Tag svTag;
  
    //! SPECTRAL_ORDER tag
  Tag soTag;

    //! order of the spectral mesh being accessed
  int spectralOrder;

    //! order of the spectral mesh being accessed plus one
  int spectralOrderp1;
};

inline SpectralMeshTool::SpectralMeshTool(Interface *impl, int order) 
        : mbImpl(impl), svTag(0), soTag(0), spectralOrder(order), spectralOrderp1(order+1)
{
  impl->query_interface(mError);
}
    
inline SpectralMeshTool::~SpectralMeshTool() 
{}

} // namespace moab 

#endif

