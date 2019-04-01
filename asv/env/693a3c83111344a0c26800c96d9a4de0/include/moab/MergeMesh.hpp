#ifndef MERGEMESH_HPP
#define MERGEMESH_HPP

#include "moab/Interface.hpp"
#include "moab/Range.hpp"

namespace moab {

class AdaptiveKDTree;

class MergeMesh
{
public:
  /* \brief Constructor
   */
  MergeMesh(Interface *mbImpl, bool printErrorIn = true);

  /* \brief Destructor
   */
  virtual ~MergeMesh();

  /* \brief Merge vertices in elements passed in
   */
  ErrorCode merge_entities(EntityHandle *elems, int elems_size,
      const double merge_tol, const int do_merge = true, const int update_sets =
          false, Tag merge_tag = 0, bool do_higher_dim = true);

  ErrorCode merge_entities(Range &elems, const double merge_tol,
      const int do_merge = true, const int update_sets = false,
      Tag merge_tag = 0, bool do_higher_dim = true);

  //Identify higher dimension to be merged
  ErrorCode merge_higher_dimensions(Range &elems);

  // merge vertices according to an input tag
  ErrorCode merge_using_integer_tag(Range & verts, Tag user_tag, Tag merge_tag=0);

  //- perform the actual merge
  ErrorCode perform_merge(Tag merged_to);

  // new method, for overlapped meshes
  // meshset could be the whole mesh, represented by root set 0;
  ErrorCode merge_all(EntityHandle meshset, const double merge_tol);
private:
  //iMesh_Instance imeshImpl;

  //- given a kdtree, set tag on vertices in leaf nodes with vertices
  //- to which they should be merged
  ErrorCode find_merged_to(EntityHandle &tree_root,
      AdaptiveKDTree &tree, Tag merged_to);

  Interface *mbImpl;

  //- the tag pointing to the entity to which an entity will be merged
  Tag mbMergeTag;

  double mergeTol, mergeTolSq;

  //- entities which will go away after the merge
  Range deadEnts;

  // vertices that were merged with other vertices, and were left in the database
  Range mergedToVertices;

  //Allow a warning to be suppressed when no merging is done
  bool printError;
};

}

#endif

