/** \file ScdInterface.hpp
 */
#ifndef SCD_INTERFACE
#define SCD_INTERFACE

#include "moab/Interface.hpp"
#include "moab/HomXform.hpp"

#include <iostream>
#include <vector>    
#include "assert.h"

namespace moab {

class StructuredElementSeq;
class EntitySequence;
class ScdVertexData;
class EntitySequence;
class ScdBox;
class ParallelComm;

/** \class ScdInterface ScdInterface.hpp "moab/ScdInterface.hpp"
 * \brief A structured mesh interface for MOAB-based data
 *
 * Structured mesh in MOAB is created and accessed through the ScdInterface and ScdBox classes.
 * 
 * \section Construction Construction and Representation
 * Structured mesh can be constructed in one of two ways.  First, a rectangular block of mesh, 
 * both vertices and edges/quads/hexes, can be created in one shot, using the construct_box method.
 * In this case, there are single sequences of vertices/entities.  The second method for creating
 * structured mesh is to create the structured blocks of vertices and elements separately.  In
 * this case, different blocks of elements can share blocks of vertices, and each block of
 * elements has its own independent parametric space.  The algorithms behind this representation
 * are described in T. Tautges, "MOAB-SD: Integrated structured and unstructured mesh representation",
 * Eng. w Comp, vol 20 no. 3.
 * 
 * Structured mesh is represented in MOAB down at the element sequence level, which is something
 * applications don't see.  In addition, when structured mesh is created, entity sets are also
 * created and tagged with information about the parametric space.  In particular, the BOX_DIMS
 * tag is used to indicate the lower and upper corners in parametric space (this
 * tag is integer size 6).  Structured mesh blocks are also available through ScdBox class objects
 * returned by ScdInterface.  These class objects should be treated only as references to the 
 * structured mesh blocks; that is, the structured mesh referenced by these objects is not deleted
 * when the ScdBox instance is destroyed.  Functions for destroying the actual mesh are are available 
 * on this class, though.
 *
 * Structured mesh blocks are returned in the form of ScdBox class objects.  Each ScdBox instance
 * represents a rectangular block of vertices and possibly elements (edges, quads, or hexes).  The
 * edge/quad/hex entity handles for a ScdBox are guaranteed to be contiguous, starting at a starting value
 * which is also available through the ScdBox class.  However, vertex handles may or may not be
 * contiguous, depending on the construction method.  The start vertex handle is also available from
 * the ScdBox class.
 *
 * \section Parameters Parametric Space
 *
 * Each structured box has a parametric (ijk) space, which can be queried through the ScdBox interface.
 * For non-periodic boxes, the edge/quad/hex parameter bounds are one less in each dimension than that
 * of the vertices, otherwise they are the same as the vertex parameter bounds.  In a parallel representation,
 * boxes are locally non-periodic by default, but global ids are assigned such that the last set of vertices
 * in a periodic direction match those of the first set of vertices in that direction.
 *
 * Entity handles are allocated with the i parameter varying fastest, then j, then k.
 *
 * \section Per Periodic Meshes
 * Boxes can be periodic in i, or j, or both i and j.  If only i or j is periodic, the corresponding mesh
 * is a strip or an annular cylinder; if both i and j are periodic, the corresponding mesh is an annular
 * torus.  A box cannot be periodic in all three parameters.  If i and/or j is periodic, and assuming
 * IMIN/JMIN is zero, the parameter extents in the/each periodic direction (IMAX/JMAX) for vertices and 
 * edges/faces/hexes are the same, and the vertices on the "top" end in the periodic direction are at
 * parameter value IMIN/JMIN.
 * 
 * \section Par Parallel Representation
 *
 * For parallel structured meshes, each local mesh (the mesh on a given process) looks like a non-periodic
 * structured mesh, and there are both local and global parameters of the structured mesh.  If the mesh is 
 * periodic in a given direction, the last process in the periodic direction has local IMAX/JMAX that is 
 * one greater than the global IMAX/JMAX.
 *

 * directions, the parameter extent is that for vertices, with edge parameter extents one fewer.
 *
 * In parallel, the periodicity described in the previous paragraph is "local periodicity"; there is also the
 * notion of global periodicity.  For serial meshes, those concepts are the same.  In parallel, a mesh can be
 * locally non-periodic but globally periodic in a given direction.  In that case, the local mesh is still
 * non-periodic, i.e. the parametric extents for edges is one fewer than that of vertices in that direction.
 * However, vertices are given global ids such that they match those of the parametric minimum in that direction.
 * Geometric positions of the vertices at the high end should still be greater than the ones just below.
 *
 * \section Adjs Adjacent Entities
 * This interface supports parametric access to intermediate-dimension entities, e.g. adjacent faces
 * and edges in a 3d mesh.  In this case, a direction parameter is added, to identify the parametric
 * direction of the entities being requested.  For example, to retrieve the faces adjacent to a hex
 * with parameters ijk, in the i parametric direction, you would use the parameters ijk0.  These 
 * intermediate entities are not stored in a structured representation, but their parametric positions
 * can be evaluated based on their adjacencies to higher-dimensional entities.  Thanks to Milad Fatenejad
 * for the thinking behind this.
 *
 * \section Evaluation Evaluation
 * The ScdBox class provides functions for evaluating the mesh based on the ijk parameter space.
 * These functions are inlined where possible, for efficiency.
*/

      //! struct for keeping parallel data in one place
class ScdParData {
public:
  ScdParData() : partMethod(NOPART), pComm(NULL) {
    gDims[0] = gDims[1] = gDims[2] = gDims[3] = gDims[4] = gDims[5] = 0;
    gPeriodic[0] = gPeriodic[1] = gPeriodic[2] = 0;
    pDims[0] = pDims[1] = pDims[2] = 0;
  }

    //! Partition method enumeration; these strategies are described in comments for
    //! compute_partition_alljorkori, compute_partition_alljkbal, compute_partition_sqij,
    //! compute_partition_sqjk, and compute_partition_sqijk
  enum PartitionMethod {ALLJORKORI = 0, ALLJKBAL, SQIJ, SQJK, SQIJK, TRIVIAL, RCBZOLTAN, NOPART};

    //! Partition method names
  static const char *PartitionMethodNames[NOPART + 1];

    //! partition method used to partition global parametric space
  int partMethod;
  
    //! lower and upper corners of global box
  int gDims[6];

    //! is globally periodic in i or j or k
  int gPeriodic[3];

    //! number of procs in each direction
  int pDims[3];

    //! parallel communicator object for this par scd mesh
  ParallelComm *pComm;
};
  
class ScdInterface 
{
public:
  friend class ScdBox;

    //! Constructor
    /** Constructor; if find_boxes is true, this will search for entity sets marked as
     * structured blocks, based on the BOX_DIMS tag.  Structured mesh blocks will be stored
     * in this interface class for future retrieval.  Structured mesh blocks created through
     * this interface will also be stored here.
     * \param impl MOAB instance
     * \param find_boxes If true, search all the entity sets, caching the structured mesh blocks
     */
  ScdInterface(Interface *impl, bool find_boxes = false);
  
    // Destructor
  ~ScdInterface();

    //! Return the MOAB Interface instance *
  Interface *impl() const;
  
    //! Construct new structured mesh box, including both vertices and elements
    /** Parameter range
     * for vertex box is [low-high], for elements is [low-high).  Construct quads by passing
     * in low[2] == high[2], and edges by passing in low[1] == high[1] and low[2] == high[2].
     * The result is passed back in a ScdBox*, which is a *reference* to the box of structured mesh.
     * That is, the actual mesh is retained in MOAB when the ScdBox is destroyed.  To actually destroy
     * the mesh, call the destroy_mesh function on the ScdBox object first, before destroying it.
     * \param low Lower corner in parameter space
     * \param high Higher corner in parameter space
     * \param coords Coordinates of vertices, interleaved (xyzxyz...); if NULL, coords are set to parametric values
     * \param num_coords Number of coordinate values
     * \param new_box Reference to box of structured mesh
     * \param lperiodic[3] If lperiodic[s] != 0, direction s is locally periodic
     * \param par_data If non-NULL, this will get stored on the ScdBox once created, contains info
     *                 about global parallel nature of ScdBox across procs
     * \param assign_global_ids If true, assigns 1-based global ids to vertices using GLOBAL_ID_TAG_NAME
     * \param resolve_shared_ents If != -1, resolves shared entities up to and including dimension equal to value
     */
  ErrorCode construct_box(HomCoord low, HomCoord high, const double * const coords, unsigned int num_coords,
                          ScdBox *& new_box, int * const lperiodic = NULL, 
                          ScdParData * const par_data = NULL,
                          bool assign_global_ids = false, int resolve_shared_ents = -1);

    //! Create a structured sequence of vertices, quads, or hexes
    /** Starting handle for the sequence is available from the returned ScdBox.  
     * If creating a structured quad or hex box, subsequent calls must be made to ScdBox::add_vbox, 
     * until all the vertices have been filled in for the box.
     * \param low Lower corner of structured box
     * \param high Higher corner of structured box
     * \param type EntityType, one of MBVERTEX, MBEDGE, MBQUAD, MBHEX
     * \param starting_id Requested start id of entities
     * \param new_box Reference to the newly created box of entities
     * \param is_periodic[3] If is_periodic[s] is non-zero, mesh should be periodic in direction s (s=[0,1,2])
     */
  ErrorCode create_scd_sequence(const HomCoord &low, const HomCoord &high, EntityType type,
                                int starting_id, ScdBox *&new_box, 
                                int *is_periodic = NULL);

    //! Return all the structured mesh blocks in this MOAB instance, as ScdBox objects
    /** Return the structured blocks in this MOAB instance.  If these were not searched for
     * at instantiation time, then the search is done now.
     * \param boxes Vector of ScdBox objects representing structured mesh blocks
     */
  ErrorCode find_boxes(std::vector<ScdBox*> &boxes);

    //! Return all the structured mesh blocks in this MOAB instance, as entity set handles
    /** Return the structured blocks in this MOAB instance.  If these were not searched for
     * at instantiation time, then the search is done now.
     * \param boxes Range of entity set objects representing structured mesh blocks
     */
  ErrorCode find_boxes(Range &boxes);

    //! Return all the structured mesh blocks known by ScdInterface (does not search)
    /** Return the structured blocks in this ScdInterface instance.  Does not search for new boxes,
     * just returns the contents of the list.
     * \param boxes Structured boxes
     */
  ErrorCode get_boxes(std::vector<ScdBox*> &boxes);

    //! Return the tag marking the lower and upper corners of boxes
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_dims_tag(bool create_if_missing = true);

    //! Return the tag marking the global lower and upper corners of boxes
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag global_box_dims_tag(bool create_if_missing = true);

    //! Return the tag marking the partitioning method used to partition the box in parallel
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag part_method_tag(bool create_if_missing = true);

    //! Return the tag marking whether box is periodic in i and j
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_periodic_tag(bool create_if_missing = true);

    //! Return the tag marking the ScdBox for a set
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_set_tag(bool create_if_missing = true);

    //! Return the ScdBox corresponding to the entity set passed in
    /** If the entity isn't a structured box set, NULL is returned.
     * \param eh Entity whose box is being queried
     */
  ScdBox *get_scd_box(EntityHandle eh);

    //! Compute a partition of structured parameter space
    /** Compute a partition of structured parameter space, based on data in the ScdParData
     * passed in.  Results are passed back in arguments, which application can set back into
     * par_data argument if they so desire.
     * \param np Number of processors
     * \param nr Rank of this processor
     * \param par_data ScdParData object that contains input global parameter space, desired
     *           partitioning method, and information about global periodicity.
     * \param ldims Local parameters for grid
     * \param lperiodic Whether or not a given dimension is locally periodic
     * \param pdims Number of procs in i, j, k directions
     */
  static ErrorCode compute_partition(int np, int nr, const ScdParData &par_data,
                                     int *ldims, int *lperiodic = NULL, int *pdims = NULL);
  
    //! Get information about the neighbor in the dijk[] direction, where dijk can be -1 or 1 for all 3 params
    /** Get information about the neighbor in the dijk[] direction, where dijk can be -1 or 1 for all 3 params
     * \param np (in) Total # procs
     * \param nr Processor from which neighbor is requested
     * \param spd (in) ScdParData containing part method, gdims, gperiodic data
     * \param dijk(*) (in) Direction being queried, = +/-1 or 0
     * \param pto (out) Processor holding the destination part
     * \param rdims(6) (out) Parametric min/max of destination part
     * \param facedims(6) (out) Parametric min/max of interface between pfrom and pto; if at the max in a periodic
     *                          direction, set to global min of that direction
     * \param across_bdy(3) (out) If across_bdy[i] is -1(1), interface with pto is across periodic lower(upper) bdy
     *                            in parameter i, 0 otherwise
     */
  static ErrorCode get_neighbor(int np, int nr, const ScdParData &spd, const int * const dijk,
                                int &pto, int *rdims, int *facedims, int *across_bdy);
  
    //! Tag vertices with sharing data for parallel representations
    /** Given the ParallelComm object to use, tag the vertices shared with other processors
     */
  ErrorCode tag_shared_vertices(ParallelComm *pcomm, EntityHandle seth);
  
    //! Tag vertices with sharing data for parallel representations
    /** Given the ParallelComm object to use, tag the vertices shared with other processors
     */
  ErrorCode tag_shared_vertices(ParallelComm *pcomm, ScdBox *box);
  
protected:
    //! Remove the box from the list on ScdInterface
  ErrorCode remove_box(ScdBox *box);
  
    //! Add the box to the list on ScdInterface
  ErrorCode add_box(ScdBox *box);
  
private:
    //! Create an entity set for a box, and tag with the parameters
    /** \param low Lower corner parameters for this box
     * \param high Upper corner parameters for this box
     * \param scd_set Entity set created
     * \param is_periodic[3] If is_periodic[s] is non-zero, mesh should be periodic in direction s (s=[0,1,2])
     */
  ErrorCode create_box_set(const HomCoord &low, const HomCoord &high,
                           EntityHandle &scd_set,
                           int *is_periodic = NULL);
  
    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by partitioning j, k, or i only.
     * If j is greater than #procs, partition that, else k, else i.
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_alljorkori(int np, int nr,
                                                       const int gijk[6], const int * const gperiodic,
                                                       int *lijk, int *lperiodic, int *pijk);
  
    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by partitioning j, and possibly k,
     * seeking square regions of jk space
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_alljkbal(int np, int nr, const int gijk[6], const int * const gperiodic,
                                                     int *lijk, int *lperiodic, int *pijk);

    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by seeking square ij partitions
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_sqij(int np, int nr, const int gijk[6], const int * const gperiodic,
                                                 int *lijk, int *lperiodic, int *pijk);
  
    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by seeking square jk partitions
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_sqjk(int np, int nr, const int gijk[6], const int * const gperiodic,
                                                 int *lijk, int *lperiodic, int *pijk);

    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by seeking square ijk partitions
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_sqijk(int np, int nr, const int gijk[6], const int * const gperiodic,
                                                  int *lijk, int *lperiodic, int *pijk);

    //! Get vertices shared with other processors
    /** Shared vertices returned as indices into each proc's handle space
     * \param box Box used to get parametric space info
     * \param procs Procs this proc shares vertices with
     * \param offsets Offsets into indices list for each proc
     * \param shared_indices local/remote indices of shared vertices
     */
  static ErrorCode get_shared_vertices(ParallelComm *pcomm, ScdBox *box, std::vector<int> &procs,
                                       std::vector<int> &offsets, std::vector<int> &shared_indices);

  static ErrorCode get_indices(const int * const ldims, const int * const rdims, const int * const across_bdy, 
                               int *face_dims, std::vector<int> &shared_indices);
  
  static ErrorCode get_neighbor_alljorkori(int np, int pfrom,
                                           const int * const gdims, const int * const gperiodic, const int * const dijk, 
                                           int &pto, int *rdims, int *facedims, int *across_bdy);
  
  static ErrorCode get_neighbor_alljkbal(int np, int pfrom,
                                         const int * const gdims, const int * const gperiodic, const int * const dijk, 
                                         int &pto, int *rdims, int *facedims, int *across_bdy);
  
  static ErrorCode get_neighbor_sqij(int np, int pfrom,
                                     const int * const gdims, const int * const gperiodic, const int * const dijk, 
                                     int &pto, int *rdims, int *facedims, int *across_bdy);
  
  static ErrorCode get_neighbor_sqjk(int np, int pfrom,
                                     const int * const gdims, const int * const gperiodic, const int * const dijk, 
                                     int &pto, int *rdims, int *facedims, int *across_bdy);
  
  static ErrorCode get_neighbor_sqijk(int np, int pfrom,
                                      const int * const gdims, const int * const gperiodic, const int * const dijk, 
                                      int &pto, int *rdims, int *facedims, int *across_bdy);
  
  static int gtol(const int *gijk, int i, int j, int k);

    //! assign global ids to vertices in this box
  ErrorCode assign_global_ids(ScdBox *box);
  
  //! interface instance
  Interface *mbImpl;

    //! whether we've searched the database for boxes yet
  bool searchedBoxes;
  
    //! structured mesh blocks; stored as ScdBox objects, can get sets from those
  std::vector<ScdBox*> scdBoxes;

    //! tag representing whether box is periodic in i and j
  Tag boxPeriodicTag;

    //! tag representing box lower and upper corners
  Tag boxDimsTag;

    //! tag representing global lower and upper corners
  Tag globalBoxDimsTag;

    //! tag representing partition method
  Tag partMethodTag;

    //! tag pointing from set to ScdBox
  Tag boxSetTag;
  
};

class ScdBox 
{
  friend class ScdInterface;
  
public:

    //! Destructor
  ~ScdBox();

    //! Return the ScdInterface responsible for this box
  inline ScdInterface *sc_impl() const;

    //! Add a vertex box to this box
    /* Add a vertex box to the element sequence referenced by this box.  The passed in vbox must
     * be a vertex box, with parametric extents no larger than that of this box.  This vbox is
     * oriented to this box by matching parameters from1-from3 in vbox to to1-to3 in this box.
     * If bb_input is true, only the part of the vertex sequence between bb_min and bb_max is referenced
     * \param vbox The vertex box being added to this box
     * \param from1 1st reference point on vbox
     * \param to1 1st reference point on this box
     * \param from2 2nd reference point on vbox
     * \param to2 2nd reference point on this box
     * \param from3 3rd reference point on vbox
     * \param to3 3rd reference point on this box
     * \param bb_input If true, subsequent parameters list extents of vbox referenced
     * \param bb_min Lower corner of rectangle referenced
     * \param bb_max Upper corner of rectangle referenced
     */
  ErrorCode add_vbox(ScdBox *vbox,
                     HomCoord from1, HomCoord to1, 
                     HomCoord from2, HomCoord to2,
                     HomCoord from3, HomCoord to3,
                     bool bb_input = false,
                     const HomCoord &bb_min = HomCoord::unitv[0],
                     const HomCoord &bb_max = HomCoord::unitv[0]);

    //! Return whether this box has all its vertices defined
    /** Tests whether vertex boxs added with add_vbox have completely defined the vertex parametric 
     * space for this box.
     *
     */
  bool boundary_complete() const;

    //! Return highest topological dimension of box
  inline int box_dimension() const;
  
    //! Starting vertex handle for this box
  inline EntityHandle start_vertex() const;
  
    //! Starting entity handle for this box
    /** If this is a vertex box, the start vertex handle is returned.
     */
  inline EntityHandle start_element() const;

    //! Return the number of elements in the box
    /* Number of elements is (boxSize[0]-1)(boxSize[1]-1)(boxSize[2]-1)
     */
  inline int num_elements() const;
  
    //! Return the number of vertices in the box
    /* Number of vertices is boxSize[0] * boxSize[1] * boxSize[2]
     */
  inline int num_vertices() const;
  
    //! Return the parametric coordinates for this box
    /**
     * \return IJK parameters of lower and upper corners
     */
  inline const int *box_dims() const;
  
    //! Return the lower corner parametric coordinates for this box
  inline HomCoord box_min() const;
  
    //! Return the upper corner parametric coordinates for this box
  inline HomCoord box_max() const;
  
    //! Return the parameter extents for this box
  inline HomCoord box_size() const;
  
    //! Return the parametric extents for this box
    /**
     * \param ijk IJK extents of this box
     */
  inline void box_size(int *ijk) const;
  
    //! Return the parametric extents for this box
    /**
     * \param i I extent of this box
     * \param j J extent of this box
     * \param k K extent of this box
     */
  void box_size(int &i, int &j, int &k) const;
  
    //! Get the element at the specified coordinates
    /**
     * \param ijk Parametric coordinates being evaluated
     */
  EntityHandle get_element(const HomCoord &ijk) const;
  
    //! Get the element at the specified coordinates
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  EntityHandle get_element(int i, int j = 0, int k = 0) const;
  
    //! Get the vertex at the specified coordinates
    /**
     * \param ijk Parametric coordinates being evaluated
     */
  EntityHandle get_vertex(const HomCoord &ijk) const;
  
    //! Get the vertex at the specified coordinates
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  EntityHandle get_vertex(int i, int j = 0, int k = 0) const;
  
    //! Get parametric coordinates of the specified entity
    /** This function returns MB_ENTITY_NOT_FOUND if the entity is not
     * in this ScdBox.
     * \param ent Entity being queried
     * \param i Parametric coordinates returned
     * \param j Parametric coordinates returned
     * \param k Parametric coordinates returned
     * \param dir Parametric coordinate direction returned (in case of getting adjacent
     *            edges (2d, 3d) or faces (3d); not modified otherwise
     */
  ErrorCode get_params(EntityHandle ent, int &i, int &j, int &k, int &dir) const;
  
    //! Get parametric coordinates of the specified entity, intermediate entities not allowed (no dir parameter)
    /** This function returns MB_ENTITY_NOT_FOUND if the entity is not
     * in this ScdBox, or MB_FAILURE if the entity is an intermediate-dimension entity.
     * \param ent Entity being queried
     * \param i Parametric coordinates returned
     * \param j Parametric coordinates returned
     * \param k Parametric coordinates returned
     */
  ErrorCode get_params(EntityHandle ent, int &i, int &j, int &k) const;
  
    //! Get parametric coordinates of the specified entity
    /** This function returns MB_ENTITY_NOT_FOUND if the entity is not
     * in this ScdBox.
     * \param ent Entity being queried
     * \param ijkd Parametric coordinates returned (including direction, in case of 
     *            getting adjacent edges (2d, 3d) or faces (3d))
     */
  ErrorCode get_params(EntityHandle ent, HomCoord &ijkd) const;
  
    /** \brief Get the adjacent edge or face at a parametric location
     * This function gets the left (i=0), front (j=0), or bottom (k=0) edge or face for a parametric element.
     * Left, front, or bottom is indicated by dir = 0, 1, or 2, resp.  All edges and faces in a structured
     * mesh block can be accessed using these parameters.
     * \param dim Dimension of adjacent entity being requested
     * \param i Parametric coordinates of cell being evaluated
     * \param j Parametric coordinates of cell being evaluated
     * \param k Parametric coordinates of cell being evaluated
     * \param dir Direction (0, 1, or 2), for getting adjacent edges (2d, 3d) or faces (3d) 
     * \param ent Entity returned from this function
     * \param create_if_missing If true, creates the entity if it doesn't already exist
     */
  ErrorCode get_adj_edge_or_face(int dim, int i, int j, int k, int dir, EntityHandle &ent,
                                 bool create_if_missing = true) const;

    //! Return whether the box contains the parameters passed in
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  bool contains(int i, int j, int k) const;

    //! Return whether the box contains the parameters passed in
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  bool contains(const HomCoord &ijk) const;

    //! Set/Get the entity set representing the box
  void box_set(EntityHandle this_set);
  EntityHandle box_set();

    //! Get coordinate arrays for vertex coordinates for a structured block
    /** Returns error if there isn't a single vertex sequence associated with this structured block
     * \param xc X coordinate array pointer returned
     * \param yc Y coordinate array pointer returned
     * \param zc Z coordinate array pointer returned
     */
  ErrorCode get_coordinate_arrays(double *&xc, double *&yc, double *&zc);
  
    //! Get read-only coordinate arrays for vertex coordinates for a structured block
    /** Returns error if there isn't a single vertex sequence associated with this structured block
     * \param xc X coordinate array pointer returned
     * \param yc Y coordinate array pointer returned
     * \param zc Z coordinate array pointer returned
     */
  ErrorCode get_coordinate_arrays(const double *&xc, const double *&yc, const double *&zc) const;

    //! Return whether box is locally periodic in i
    /** Return whether box is locally periodic in i
     * \return True if box is locally periodic in i direction
     */
  bool locally_periodic_i() const;
  
    //! Return whether box is locally periodic in j
    /** Return whether box is locally periodic in j
     * \return True if box is locally periodic in j direction
     */
  bool locally_periodic_j() const;
  
    //! Return whether box is locally periodic in k
    /** Return whether box is locally periodic in k
     * \return True if box is locally periodic in k direction
     */
  bool locally_periodic_k() const;
  
    //! Set local periodicity
    /** 
     * \param lperiodic Vector of ijk periodicities to set this box to
      */
  void locally_periodic(bool lperiodic[3]);

    //! Get local periodicity
    /** 
     * \return Vector of ijk periodicities for this box
     */
  const int *locally_periodic() const;
 
    //! Return parallel data 
    /** Return parallel data, if there is any
     * \return par_data Parallel data set on this box 
     */
  ScdParData &par_data() {return parData;}
  
    //! Return parallel data 
    /** Return parallel data, if there is any
     * \return par_data Parallel data set on this box 
     */
  const ScdParData &par_data() const {return parData;}
  
    //! set parallel data 
    /** Set parallel data for this box
     * \param par_data Parallel data to be set on this box 
     */
  void par_data(const ScdParData &par_datap) {parData = par_datap;}
  
private:
    //! Constructor
    /** Create a structured box instance; this constructor is private because it should only be called
     * from ScdInterface, a friend class.  This constructor takes two sequences, one of which can be
     * NULL.  If both sequences come in non-NULL, the first should be a VertexSequence* corresponding to
     * a structured vertex sequence and the second should be a StructuredElementSeq*.  If the 2nd is NULL,
     * the first can be either of those types.  The other members of this class are taken from the sequences
     * (e.g. parametric space) or the box set argument.  Tags on the box set should be set from the caller.
     * \param sc_impl A ScdInterface instance
     * \param box_set Entity set representing this rectangle
     * \param seq1 An EntitySequence (see ScdBox description)
     * \param seq2 An EntitySequence (see ScdBox description), or NULL
     */
  ScdBox(ScdInterface *sc_impl, EntityHandle box_set,
         EntitySequence *seq1, EntitySequence *seq2 = NULL);

    //! function to get vertex handle directly from sequence
    /** \param i Parameter being queried
     * \param j Parameter being queried
     * \param k Parameter being queried
     */
  EntityHandle get_vertex_from_seq(int i, int j, int k) const;

    //! set the vertex sequence
  ErrorCode vert_dat(ScdVertexData *vert_dat);
  
    //! get the vertex sequence
  ScdVertexData *vert_dat() const;
  
    //! set the element sequence
  ErrorCode elem_seq(EntitySequence *elem_seq);
  
    //! get the element sequence
  StructuredElementSeq *elem_seq() const;
  
    //! Set the starting vertex handle for this box
  void start_vertex(EntityHandle startv);
  
    //! Set the starting entity handle for this box
  void start_element(EntityHandle starte);

    //! interface instance
  ScdInterface *scImpl;

    //! entity set representing this box
  EntityHandle boxSet;

    //! vertex sequence this box represents, if there's only one, otherwise they're
    //! retrieved from the element sequence
  ScdVertexData *vertDat;

    //! element sequence this box represents
  StructuredElementSeq *elemSeq;
  
    //! starting vertex handle for this box
  EntityHandle startVertex;
  
    //! starting element handle for this box
  EntityHandle startElem;

    //! lower and upper corners
  int boxDims[6];

    //! is locally periodic in i or j or k
  int locallyPeriodic[3];

    //! parallel data associated with this box, if any
  ScdParData parData;
  
    //! parameter extents
  HomCoord boxSize;

    //! convenience parameters, (boxSize[1]-1)*(boxSize[0]-1) and boxSize[0]-1
  int boxSizeIJ;
  int boxSizeIJM1;
  int boxSizeIM1;
  
};

inline ErrorCode ScdInterface::compute_partition(int np, int nr, const ScdParData &par_data,
                                     int *ldims, int *lperiodic, int *pdims)
{
  ErrorCode rval = MB_SUCCESS;
  switch (par_data.partMethod) {
    case ScdParData::ALLJORKORI:
    case -1:
        rval = compute_partition_alljorkori(np, nr, par_data.gDims, par_data.gPeriodic, 
                                            ldims, lperiodic, pdims);
        break;
    case ScdParData::ALLJKBAL:
        rval = compute_partition_alljkbal(np, nr, par_data.gDims, par_data.gPeriodic, 
                                          ldims, lperiodic, pdims);
        break;
    case ScdParData::SQIJ:
        rval = compute_partition_sqij(np, nr, par_data.gDims, par_data.gPeriodic, 
                                      ldims, lperiodic, pdims);
        break;
    case ScdParData::SQJK:
        rval = compute_partition_sqjk(np, nr, par_data.gDims, par_data.gPeriodic, 
                                      ldims, lperiodic, pdims);
        break;
    case ScdParData::SQIJK:
        rval = compute_partition_sqijk(np, nr, par_data.gDims, par_data.gPeriodic, 
                                      ldims, lperiodic, pdims);
        break;
    default:
        rval = MB_FAILURE;
        break;
  }

  return rval;
}

inline ErrorCode ScdInterface::compute_partition_alljorkori(int np, int nr,
                                                            const int gijk[6], const int * const gperiodic,
                                                            int *ldims, int *lperiodic, int *pijk)
{
    // partition *the elements* over the parametric space; 1d partition for now, in the j, k, or i
    // parameters
  int tmp_lp[3], tmp_pijk[3];
  if (!lperiodic) lperiodic = tmp_lp;
  if (!pijk) pijk = tmp_pijk;
  
  for (int i = 0; i < 3; i++) lperiodic[i] = gperiodic[i];
  
  if (np == 1) {
    if (ldims) {
      ldims[0] = gijk[0];
      ldims[3] = gijk[3];
      ldims[1] = gijk[1];
      ldims[4] = gijk[4];
      ldims[2] = gijk[2];
      ldims[5] = gijk[5];
    }
    pijk[0] = pijk[1] = pijk[2] = 1;
  }
  else {
    if (gijk[4] - gijk[1] > np) {
        // partition j over procs
      int dj = (gijk[4] - gijk[1]) / np;
      int extra = (gijk[4] - gijk[1]) % np;
      ldims[1] = gijk[1] + nr*dj + 
          std::min(nr, extra);
      ldims[4] = ldims[1] + dj + (nr < extra ? 1 : 0);

      if (gperiodic[1] && np > 1) {
        lperiodic[1] = 0;
        ldims[4]++;
      }
      
      ldims[2] = gijk[2]; ldims[5] = gijk[5];
      ldims[0] = gijk[0]; ldims[3] = gijk[3];
      pijk[0] = pijk[2] = 1;
      pijk[1] = np;
    }
    else if (gijk[5] - gijk[2] > np) {
        // partition k over procs
      int dk = (gijk[5] - gijk[2]) / np;
      int extra = (gijk[5] - gijk[2]) % np;
      ldims[2] = gijk[2] + nr*dk + 
          std::min(nr, extra);
      ldims[5] = ldims[2] + dk + (nr < extra ? 1 : 0);

      ldims[1] = gijk[1]; ldims[4] = gijk[4];
      ldims[0] = gijk[0]; ldims[3] = gijk[3];
      pijk[0] = pijk[1] = 1;
      pijk[2] = np;
    }
    else if (gijk[3] - gijk[0] > np) {
        // partition i over procs
      int di = (gijk[3] - gijk[0]) / np;
      int extra = (gijk[3] - gijk[0]) % np;
      ldims[0] = gijk[0] + nr*di + 
          std::min(nr, extra);
      ldims[3] = ldims[0] + di + (nr < extra ? 1 : 0);

      if (gperiodic[0] && np > 1) {
        lperiodic[0] = 0;
        ldims[3]++;
      }

      ldims[2] = gijk[2]; ldims[5] = gijk[5];
      ldims[1] = gijk[1]; ldims[4] = gijk[4];

      pijk[1] = pijk[2] = 1;
      pijk[0] = np;
    }
    else {
        // Couldn't find a suitable partition...
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

inline ErrorCode ScdInterface::compute_partition_alljkbal(int np, int nr,
                                                          const int gijk[6], const int * const gperiodic,
                                                          int *ldims, int *lperiodic, int *pijk)
{
  int tmp_lp[3], tmp_pijk[3];
  if (!lperiodic) lperiodic = tmp_lp;
  if (!pijk) pijk = tmp_pijk;

  for (int i = 0; i < 3; i++) lperiodic[i] = gperiodic[i];

  if (np == 1) {
    if (ldims) {
      ldims[0] = gijk[0];
      ldims[3] = gijk[3];
      ldims[1] = gijk[1];
      ldims[4] = gijk[4];
      ldims[2] = gijk[2];
      ldims[5] = gijk[5];
    }
    pijk[0] = pijk[1] = pijk[2] = 1;
  }
  else {
      // improved, possibly 2-d partition
    std::vector<double> kfactors;
    kfactors.push_back(1);
    int K = gijk[5] - gijk[2];
    for (int i = 2; i < K; i++) 
      if (!(K%i) && !(np%i)) kfactors.push_back(i);
    kfactors.push_back(K);
  
      // compute the ideal nj and nk
    int J = gijk[4] - gijk[1];
    double njideal = sqrt(((double)(np*J))/((double)K));
    double nkideal = (njideal*K)/J;
  
    int nk, nj;
    if (nkideal < 1.0) {
      nk = 1;
      nj = np;
    }
    else {
      std::vector<double>::iterator vit = std::lower_bound(kfactors.begin(), kfactors.end(), nkideal);
      if (vit == kfactors.begin()) nk = 1;
      else nk = (int)*(--vit);
      nj = np / nk;
    }

    int dk = K / nk;
    int dj = J / nj;
  
    ldims[2] = gijk[2] + (nr % nk) * dk;
    ldims[5] = ldims[2] + dk;
  
    int extra = J % nj;
  
    ldims[1] = gijk[1] + (nr / nk) * dj + std::min(nr / nk, extra);
    ldims[4] = ldims[1] + dj + (nr / nk < extra ? 1 : 0);

    ldims[0] = gijk[0];
    ldims[3] = gijk[3];

    if (gperiodic[1] && np > 1) {
      lperiodic[1] = 0;
      if (nr/nk == nj-1) {
        ldims[1]++;
      }
    }

    pijk[0] = 1; pijk[1] = nj; pijk[2] = nk;
  }
  
  return MB_SUCCESS;
}

inline ErrorCode ScdInterface::compute_partition_sqij(int np, int nr,
                                                      const int gijk[6], const int * const gperiodic,
                                                      int *ldims, int *lperiodic, int *pijk)
{
  int tmp_lp[3], tmp_pijk[3];
  if (!lperiodic) lperiodic = tmp_lp;
  if (!pijk) pijk = tmp_pijk;

    // square IxJ partition

  for (int i = 0; i < 3; i++) lperiodic[i] = gperiodic[i];
  
  if (np == 1) {
    if (ldims) {
      ldims[0] = gijk[0];
      ldims[3] = gijk[3];
      ldims[1] = gijk[1];
      ldims[4] = gijk[4];
      ldims[2] = gijk[2];
      ldims[5] = gijk[5];
    }
    pijk[0] = pijk[1] = pijk[2] = 1;
  }
  else {
    std::vector<double> pfactors, ppfactors;
    for (int i = 2; i <= np/2; i++) 
      if (!(np%i)) {
        pfactors.push_back(i);
        ppfactors.push_back(((double)(i*i))/np);
      }
    pfactors.push_back(np);
    ppfactors.push_back( (double)np);
    
      // ideally, Px/Py = I/J
    double ijratio = ((double)(gijk[3]-gijk[0]))/((double)(gijk[4]-gijk[1]));
    
    unsigned int ind = std::lower_bound(ppfactors.begin(), ppfactors.end(), ijratio) - ppfactors.begin();
    if (ind && fabs(ppfactors[ind-1]-ijratio) < fabs(ppfactors[ind]-ijratio)) ind--;
    

      // VARIABLES DESCRIBING THE MESH:
      // pi, pj = # procs in i and j directions
      // nri, nrj = my proc's position in i, j directions
      // I, J = # edges/elements in i, j directions
      // iextra, jextra = # procs having extra edge in i/j direction
      // top_i, top_j = if true, I'm the last proc in the i/j direction
      // i, j = # edges locally in i/j direction, *not* including one for iextra/jextra
    int pi = pfactors[ind];
    int pj = np / pi;
    
    int I = (gijk[3] - gijk[0]), J = (gijk[4] - gijk[1]);
    int iextra = I%pi, jextra = J%pj, i = I/pi, j = J/pj;
    int nri = nr % pi, nrj = nr / pi;

    if (ldims) {
      ldims[0] = gijk[0] + i*nri + std::min(iextra, nri);
      ldims[3] = ldims[0] + i + (nri < iextra ? 1 : 0);
      ldims[1] = gijk[1] + j*nrj + std::min(jextra, nrj);
      ldims[4] = ldims[1] + j + (nrj < jextra ? 1 : 0);
      
      ldims[2] = gijk[2];
      ldims[5] = gijk[5];

      if (gperiodic[0] && pi > 1) {
        lperiodic[0] = 0;
        if (nri == pi-1) ldims[3]++;
      }
      if (gperiodic[1] && pj > 1) {
        lperiodic[1] = 0;
        if (nrj == pj-1) ldims[4]++;
      }
      
    }

    pijk[0] = pi; pijk[1] = pj; pijk[2] = 1;
  }  

  return MB_SUCCESS;
}

inline ErrorCode ScdInterface::compute_partition_sqjk(int np, int nr,
                                                      const int gijk[6], const int * const gperiodic,
                                                      int *ldims, int *lperiodic, int *pijk)
{
  int tmp_lp[3], tmp_pijk[3];
  if (!lperiodic) lperiodic = tmp_lp;
  if (!pijk) pijk = tmp_pijk;

    // square JxK partition
  for (int i = 0; i < 3; i++) lperiodic[i] = gperiodic[i];

  if (np == 1) {
    if (ldims) {
      ldims[0] = gijk[0];
      ldims[3] = gijk[3];
      ldims[1] = gijk[1];
      ldims[4] = gijk[4];
      ldims[2] = gijk[2];
      ldims[5] = gijk[5];
    }
    pijk[0] = pijk[1] = pijk[2] = 1;
  }
  else {
    std::vector<double> pfactors, ppfactors;
    for (int p = 2; p <= np; p++) 
      if (!(np%p)) {
        pfactors.push_back(p);
        ppfactors.push_back(((double)(p*p))/np);
      }
  
      // ideally, Pj/Pk = J/K
    int pj, pk;
    if (gijk[5] == gijk[2]) {
      pk = 1;
      pj = np;
    }
    else {
      double jkratio = ((double)(gijk[4]-gijk[1]))/((double)(gijk[5]-gijk[2]));

      std::vector<double>::iterator vit  = std::lower_bound(ppfactors.begin(), ppfactors.end(), jkratio);
      if (vit == ppfactors.end()) --vit;
      else if (vit != ppfactors.begin() && fabs(*(vit-1)-jkratio) < fabs((*vit)-jkratio)) --vit;
      int ind = vit - ppfactors.begin();
  
      pj = 1;
      if (ind >= 0 && !pfactors.empty()) pfactors[ind];
      pk = np / pj;
    }

    int K = (gijk[5] - gijk[2]), J = (gijk[4] - gijk[1]);
    int jextra = J%pj, kextra = K%pk, j = J/pj, k = K/pk;
    int nrj = nr % pj, nrk = nr / pj;
    ldims[1] = gijk[1] + j*nrj + std::min(jextra, nrj);
    ldims[4] = ldims[1] + j + (nrj < jextra ? 1 : 0);
    ldims[2] = gijk[2] + k*nrk + std::min(kextra, nrk);
    ldims[5] = ldims[2] + k + (nrk < kextra ? 1 : 0);

    ldims[0] = gijk[0];
    ldims[3] = gijk[3];

    if (gperiodic[1] && pj > 1) {
      lperiodic[1] = 0;
      if (nrj == pj-1) ldims[4]++;
    }

    pijk[0] = 1; pijk[1] = pj; pijk[2] = pk;
  }
  
  return MB_SUCCESS;
}

inline ErrorCode ScdInterface::compute_partition_sqijk(int np, int nr,
                                                       const int * const gijk, const int * const gperiodic, 
                                                       int *ldims, int *lperiodic, int *pijk)
{
  if (gperiodic[0] || gperiodic[1] || gperiodic[2]) return MB_FAILURE;
  
  int tmp_lp[3], tmp_pijk[3];
  if (!lperiodic) lperiodic = tmp_lp;
  if (!pijk) pijk = tmp_pijk;

    // square IxJxK partition

  for (int i = 0; i < 3; i++) lperiodic[i] = gperiodic[i];
  
  if (np == 1) {
    if (ldims) 
      for (int i = 0; i < 6; i++) ldims[i] = gijk[i];
    pijk[0] = pijk[1] = pijk[2] = 1;
    return MB_SUCCESS;
  }

  std::vector<int> pfactors;
  pfactors.push_back(1);
  for (int i = 2; i <= np/2; i++) 
    if (!(np%i)) pfactors.push_back(i);
  pfactors.push_back(np);

    // test for IJ, JK, IK
  int IJK[3], dIJK[3];
  for (int i = 0; i < 3; i++)
    IJK[i] = std::max(gijk[3+i] - gijk[i], 1);
    // order IJK from lo to hi
  int lo = 0, hi = 0;
  for (int i = 1; i < 3; i++) {
    if (IJK[i] < IJK[lo]) lo = i;
    if (IJK[i] > IJK[hi]) hi = i;
  }
  if (lo == hi) hi = (lo+1)%3;
  int mid = 3 - lo - hi;
    // search for perfect subdivision of np that balances #cells
  int perfa_best = -1, perfb_best = -1;
  double ratio = 0.0;
  for (int po = 0; po < (int)pfactors.size(); po++) {
    for (int pi = po; pi < (int)pfactors.size() && np/(pfactors[po]*pfactors[pi]) >= pfactors[pi]; pi++) {
      int p3 = std::find(pfactors.begin(), pfactors.end(), np/(pfactors[po]*pfactors[pi])) - pfactors.begin();
      if (p3 == (int)pfactors.size() || pfactors[po]*pfactors[pi]*pfactors[p3] != np) continue; // po*pi should exactly factor np
      assert(po <= pi && pi <= p3);
        // by definition, po <= pi <= p3
      double minl = std::min(std::min(IJK[lo]/pfactors[po], IJK[mid]/pfactors[pi]), IJK[hi]/pfactors[p3]),
          maxl = std::max(std::max(IJK[lo]/pfactors[po], IJK[mid]/pfactors[pi]), IJK[hi]/pfactors[p3]);
      if (minl/maxl > ratio) {
        ratio = minl/maxl;
        perfa_best = po; perfb_best = pi;
      }
    }
  }
  if (perfa_best == -1 || perfb_best == -1) 
    return MB_FAILURE;

    // VARIABLES DESCRIBING THE MESH:
    // pijk[i] = # procs in direction i
    // numr[i] = my proc's position in direction i
    // dIJK[i] = # edges/elements in direction i
    // extra[i]= # procs having extra edge in direction i
    // top[i] = if true, I'm the last proc in direction i
    
  pijk[lo] = pfactors[perfa_best]; pijk[mid] = pfactors[perfb_best]; 
  pijk[hi] = (np/(pfactors[perfa_best]*pfactors[perfb_best]));
  int extra[3] = {0, 0, 0}, numr[3];
  for (int i = 0; i < 3; i++) {
    dIJK[i] = IJK[i] / pijk[i];
    extra[i] = IJK[i]%pijk[i];
  }
  numr[2] = nr / (pijk[0]*pijk[1]);
  int rem = nr % (pijk[0] * pijk[1]);
  numr[1] = rem / pijk[0];
  numr[0] = rem % pijk[0];
  for (int i = 0; i < 3; i++) {
    extra[i] = IJK[i] % dIJK[i];
    ldims[i] = gijk[i] + numr[i]*dIJK[i] + std::min(extra[i], numr[i]);
    ldims[3+i] = ldims[i] + dIJK[i] + (numr[i] < extra[i] ? 1 : 0);
  }

  return MB_SUCCESS;
}

inline int ScdInterface::gtol(const int *gijk, int i, int j, int k) 
{
  return ((k-gijk[2])*(gijk[3]-gijk[0]+1)*(gijk[4]-gijk[1]+1) + (j-gijk[1])*(gijk[3]-gijk[0]+1) + i-gijk[0]);
}

inline ErrorCode ScdInterface::get_indices(const int * const ldims, const int * const rdims, const int * const across_bdy, 
                                           int *face_dims, std::vector<int> &shared_indices) 
{
    // check for going across periodic bdy and face_dims not in my ldims (I'll always be on top in that case)...
  if (across_bdy[0] > 0 && face_dims[0] != ldims[3]) face_dims[0] = face_dims[3] = ldims[3];
  else if (across_bdy[0] < 0 && face_dims[0] != ldims[0]) face_dims[0] = face_dims[3] = ldims[0];
  if (across_bdy[1] > 0 && face_dims[1] != ldims[4]) face_dims[1] = face_dims[4] = ldims[4];
  else if (across_bdy[1] < 0 && face_dims[1] != ldims[1]) face_dims[0] = face_dims[3] = ldims[1];
  
  for (int k = face_dims[2]; k <= face_dims[5]; k++)
    for (int j = face_dims[1]; j <= face_dims[4]; j++)
      for (int i = face_dims[0]; i <= face_dims[3]; i++)
        shared_indices.push_back(gtol(ldims, i, j, k));

  if (across_bdy[0] > 0 && face_dims[0] != rdims[0]) face_dims[0] = face_dims[3] = rdims[0];
  else if (across_bdy[0] < 0 && face_dims[0] != rdims[3]) face_dims[0] = face_dims[3] = rdims[3];
  if (across_bdy[1] > 0 && face_dims[1] != rdims[1]) face_dims[1] = face_dims[4] = rdims[1];
  else if (across_bdy[1] < 0 && face_dims[1] != rdims[4]) face_dims[0] = face_dims[3] = rdims[4];
  
  for (int k = face_dims[2]; k <= face_dims[5]; k++)
    for (int j = face_dims[1]; j <= face_dims[4]; j++)
      for (int i = face_dims[0]; i <= face_dims[3]; i++)
        shared_indices.push_back(gtol(rdims, i, j, k));

  return MB_SUCCESS;
}
  
inline ErrorCode ScdInterface::get_neighbor(int np, int pfrom, const ScdParData &spd, const int * const dijk, 
                                            int &pto, int *rdims, int *facedims, int *across_bdy)
{
  if (!dijk[0] && !dijk[1] && !dijk[2]) {
    // not going anywhere, return
    pto = -1;
    return MB_SUCCESS;
  }
  
  switch (spd.partMethod) {
    case ScdParData::ALLJORKORI:
    case -1:
        return get_neighbor_alljorkori(np, pfrom, spd.gDims, spd.gPeriodic, dijk,
                                       pto, rdims, facedims, across_bdy);
    case ScdParData::ALLJKBAL:
        return get_neighbor_alljkbal(np, pfrom, spd.gDims, spd.gPeriodic, dijk,
                                     pto, rdims, facedims, across_bdy);
    case ScdParData::SQIJ:
        return get_neighbor_sqij(np, pfrom, spd.gDims, spd.gPeriodic, dijk,
                                 pto, rdims, facedims, across_bdy);
    case ScdParData::SQJK:
        return get_neighbor_sqjk(np, pfrom, spd.gDims, spd.gPeriodic, dijk,
                                 pto, rdims, facedims, across_bdy);
    case ScdParData::SQIJK:
        return get_neighbor_sqijk(np, pfrom, spd.gDims, spd.gPeriodic, dijk,
                                  pto, rdims, facedims, across_bdy);
    default:
        break;
  }

  return MB_FAILURE;
}

inline ErrorCode ScdInterface::tag_shared_vertices(ParallelComm *pcomm, EntityHandle seth) 
{
  ScdBox *box = get_scd_box(seth);
  if (!box) {
      // look for contained boxes
    Range tmp_range;
    ErrorCode rval = mbImpl->get_entities_by_type(seth, MBENTITYSET, tmp_range);
    if (MB_SUCCESS != rval) return rval;
    for (Range::iterator rit = tmp_range.begin(); rit != tmp_range.end(); ++rit) {
      box = get_scd_box(*rit);
      if (box) break;
    }
  }
  
  if (!box) return MB_FAILURE;

  return tag_shared_vertices(pcomm, box);
}

inline ScdInterface *ScdBox::sc_impl() const 
{
  return scImpl;
}

inline EntityHandle ScdBox::start_vertex() const
{
  return startVertex;
}
    
inline void ScdBox::start_vertex(EntityHandle startv) 
{
  startVertex = startv;
}
    
inline EntityHandle ScdBox::start_element() const
{
  return startElem;
}
    
inline void ScdBox::start_element(EntityHandle starte) 
{
  startElem = starte;
}
    
inline int ScdBox::num_elements() const
{
  return (!startElem ? 0 : 
          (boxSize[0]- (locallyPeriodic[0] ? 0 : 1)) * 
          (-1 == boxSize[1] ? 1 : (boxSize[1]-(locallyPeriodic[1] ? 0 : 1))) * 
          (boxSize[2] == -1 || boxSize[2] == 1 ? 1 : (boxSize[2]-(locallyPeriodic[2] ? 0 : 1))));
}
    
inline int ScdBox::num_vertices() const
{
  return boxSize[0] * (!boxSize[1] ? 1 : boxSize[1]) * 
      (!boxSize[2] ? 1 : boxSize[2]);
}
    
inline const int *ScdBox::box_dims() const 
{
  return boxDims;
}

inline HomCoord ScdBox::box_min() const 
{
  return HomCoord(boxDims, 3);
}

inline HomCoord ScdBox::box_max() const 
{
  return HomCoord(boxDims+3, 3);
}

inline HomCoord ScdBox::box_size() const 
{
  return boxSize;
}

inline void ScdBox::box_size(int *ijk) const 
{
  ijk[0] = boxSize[0];
  ijk[1] = boxSize[1];
  ijk[2] = boxSize[2];
}

inline void ScdBox::box_size(int &i, int &j, int &k) const 
{
  i = boxSize[0];
  j = boxSize[1];
  k = boxSize[2];
}
  
inline EntityHandle ScdBox::get_element(int i, int j, int k) const 
{
  return (!startElem ? 0 : 
          startElem + (k-boxDims[2])*boxSizeIJM1 + (j-boxDims[1])*boxSizeIM1 + i-boxDims[0]);
}

inline EntityHandle ScdBox::get_element(const HomCoord &ijk) const
{
  return get_element(ijk[0], ijk[1], ijk[2]);
}
  
inline EntityHandle ScdBox::get_vertex(int i, int j, int k) const 
{
  return (vertDat ? startVertex + (boxDims[2] == -1 && boxDims[5] == -1 ? 0 : (k-boxDims[2]))*boxSizeIJ +
	  (boxDims[1] == -1 && boxDims[4] == -1 ? 0 : (j-boxDims[1]))*boxSize[0] + i-boxDims[0] : get_vertex_from_seq(i, j, k));
}

inline EntityHandle ScdBox::get_vertex(const HomCoord &ijk) const
{
  return get_vertex(ijk[0], ijk[1], ijk[2]);
}
  
inline bool ScdBox::contains(const HomCoord &ijk) const
{
  return (ijk >= HomCoord(boxDims, 3) && 
          ijk <= HomCoord(boxDims+3, 3));
}

inline bool ScdBox::contains(int i, int j, int k) const
{
  return contains(HomCoord(i, j, k));
}

inline void ScdBox::box_set(EntityHandle this_set) 
{
  boxSet = this_set;
}

inline EntityHandle ScdBox::box_set() 
{
  return boxSet;
}
    
inline ScdVertexData *ScdBox::vert_dat() const
{
  return vertDat;
}
  
inline StructuredElementSeq *ScdBox::elem_seq() const
{
  return elemSeq;
}
  
inline ErrorCode ScdBox::get_params(EntityHandle ent, int &i, int &j, int &k, int &dir) const
{
  HomCoord hc;
  ErrorCode rval = get_params(ent, hc);
  if (MB_SUCCESS == rval) {
    i = hc[0];
    j = hc[1];
    k = hc[2];
    dir = hc[3];
  }
  
  return rval;
}

inline ErrorCode ScdBox::get_params(EntityHandle ent, int &i, int &j, int &k) const 
{
  HomCoord hc;
  ErrorCode rval = get_params(ent, hc);
  if (MB_SUCCESS == rval) {
    i = hc[0];
    j = hc[1];
    k = hc[2];
  }
  
  return rval;
}

inline bool ScdBox::locally_periodic_i() const 
{
  return locallyPeriodic[0];
}

inline bool ScdBox::locally_periodic_j() const 
{
  return locallyPeriodic[1];
}

inline bool ScdBox::locally_periodic_k() const 
{
  return locallyPeriodic[2];
}

inline void ScdBox::locally_periodic(bool lperiodic[3])
{
   for (int i = 0; i < 3; i++) 
    locallyPeriodic[i] = lperiodic[i];
}

inline const int *ScdBox::locally_periodic() const
{
  return locallyPeriodic;
}
 
inline std::ostream &operator<<(std::ostream &str, const ScdParData &pd) 
{
  str << "Partition method = " << ScdParData::PartitionMethodNames[pd.partMethod] << ", gDims = (" 
      << pd.gDims[0] << "," << pd.gDims[1] << "," << pd.gDims[2] << ")-(" 
      << pd.gDims[3] << "," << pd.gDims[4] << "," << pd.gDims[5] << "), gPeriodic = (" 
      << pd.gPeriodic[0] << "," << pd.gPeriodic[1] << "," << pd.gPeriodic[2] << "), pDims = ("
      << pd.pDims[0] << "," << pd.pDims[1] << "," << pd.pDims[2] << ")" << std::endl;
  return str;
}

} // namespace moab
#endif
