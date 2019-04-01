!> @example DirectAccessNoHolesF90.F90
!! \brief Use direct access to MOAB data to avoid calling through API, in Fortran90 \n
!!
!! This example creates a 1d row of quad elements, such that all quad and vertex handles
!! are contiguous in the handle space and in the database.  Then it shows how to get access
!! to pointers to MOAB-native data for vertex coordinates, quad connectivity, and tag storage
!! (vertex to quad adjacency lists aren't accessible from Fortran because they are std::vector's).
!! This allows applications to access this data directly
!! without going through MOAB's API.  In cases where the mesh is not changing (or only mesh
!! vertices are moving), this can save significant execution time in applications.
!!
!! Using direct access (or MOAB in general) from Fortran is complicated in two specific ways:
!! 1) There is no way to use arrays with specified dimension starting/ending values with ISO_C_BINDING;
!!    therefore, all arrays passed back from MOAB/iMesh must use 1-based indices; this makes it a bit
!!    more difficult to traverse the direct arrays.  In this example, I explicitly use indices that
!!    look like my_array(1+v_ind...) to remind users of this.
!! 2) Arithmetic on handles is complicated by the fact that Fortran has no unsigned integer type.  I get
!!    around this by assigning to a C-typed variable first, before the handle arithmetic.  This seems to
!!    work fine.  C-typed variables are part of the Fortran95 standard.
!!
!!  ----------------------
!!  |      |      |      |
!!  |      |      |      | ...
!!  |      |      |      |
!!  ----------------------
!!
!!    -#  Initialize MOAB \n
!!    -#  Create a quad mesh, as depicted above
!!    -#  Create 2 dense tags (tag1, tag2) for avg position to assign to quads, and # verts per quad (tag3)
!!    -#  Get connectivity, coordinate, tag1 iterators
!!    -#  Iterate through quads, computing midpoint based on vertex positions, set on quad-based tag1
!!    -#  Iterate through vertices, summing positions into tag2 on connected quads and incrementing vertex count
!!    -#  Iterate through quads, normalizing tag2 by vertex count and comparing values of tag1 and tag2
!!
!! <b>To compile</b>: \n
!!    make DirectAccessNoHolesF90 MOAB_DIR=<installdir>  \n
!! <b>To run</b>: ./DirectAccessNoHolesF90 [-nquads <# quads>]\n
!!
!! VSM: Note that IBM xlf compilers do not like dynamic sizing functions for
!!      C_F_POINTER calls. This leads to internal compiler failure.
!!      http://www-01.ibm.com/support/docview.wss?uid=swg1IV49428
!!
#define CHECK(a)   if (a .ne. iBase_SUCCESS) call exit(a)

program DirectAccessNoHolesF90
  use, intrinsic :: ISO_C_BINDING
  implicit none
#include "iMesh_f.h"

  ! declarations
  iMesh_Instance :: imesh
  iBase_EntitySetHandle root_set
  integer :: nverts, tmpflag
  integer ierr, nquads, nquads_tmp
  iBase_TagHandle tag1h, tag2h
  TYPE(C_PTR) verts_ptr, quads_ptr, conn_ptr, x_ptr, y_ptr, z_ptr, tag1_ptr, tag2_ptr
  TYPE(C_PTR) tmp_quads_ptr, offsets_ptr  ! pointers that are equivalence'd to arrays using c_f_pointer
  real(C_DOUBLE), pointer :: x(:), y(:), z(:), tag1(:), tag2(:) ! arrays equivalence'd to pointers
  integer, pointer :: offsets(:)                       ! arrays equivalence'd to pointers
  iBase_EntityHandle, pointer :: quads(:), verts(:), conn(:), tmp_quads(:)  ! arrays equivalence'd to pointers
  iBase_EntityArrIterator viter, qiter
  integer(C_SIZE_T) :: startv, startq, v_ind, e_ind  ! needed to do handle arithmetic
  integer count, vpere, e, i, j, v, offsets_size, tmp_quads_size, n_dis
  character*120 opt

  ! initialize interface and get root set
  call iMesh_newMesh("", %REF(imesh), ierr)
  CHECK(ierr)
  call iMesh_getRootSet(%val(imesh), root_set, ierr)
  CHECK(ierr)

  ! create mesh
  nquads_tmp = 1000
  call create_mesh_no_holes(imesh, nquads_tmp, ierr)
  CHECK(ierr)

  ! get all vertices and quads as arrays
  nverts = 0
  call iMesh_getEntities(%val(imesh), %val(root_set), %val(0), %val(iBase_VERTEX), &
       verts_ptr, nverts, nverts, ierr)
  CHECK(ierr)
  call c_f_pointer(verts_ptr, verts, [nverts])
  nquads = 0
  call iMesh_getEntities(%val(imesh), %val(root_set), %val(2), %val(iMesh_QUADRILATERAL), &
       quads_ptr, nquads, nquads, ierr)
  CHECK(ierr)
  call c_f_pointer(quads_ptr, quads, [nquads])

  ! First vertex/quad is at start of vertex/quad list, and is offset for vertex/quad index calculation
  startv = verts(1)
  startq = quads(1)
  call iMesh_freeMemory(%val(imesh), quads_ptr)

  ! create tag1 (element-based avg), tag2 (vertex-based avg)
  opt = 'moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_valUE=0'
  call iMesh_createTagWithOptions(%val(imesh), 'tag1', opt, %val(3), %val(iBase_DOUBLE), &
       tag1h, ierr, %val(4), %val(56))
  CHECK(ierr)
  call iMesh_createTagWithOptions(%val(imesh), 'tag2', opt, %val(3), %val(iBase_DOUBLE), &
       tag2h, ierr, %val(4), %val(56))
  CHECK(ierr)

  ! Get iterator to verts, then pointer to coordinates
  call iMesh_initEntArrIter(%val(imesh), %val(root_set), %val(iBase_VERTEX), %val(iMesh_ALL_TOPOLOGIES), &
       %val(nverts), %val(0), viter, ierr)
  CHECK(ierr)
  call iMesh_coordsIterate(%val(imesh), %val(viter), x_ptr, y_ptr, z_ptr, count, ierr)
  CHECK(ierr)
  CHECK(count-nverts) ! check that all verts are in one contiguous chunk
  call c_f_pointer(x_ptr, x, [nverts])
  call c_f_pointer(y_ptr, y, [nverts])
  call c_f_pointer(z_ptr, z, [nverts])

  ! Get iterator to quads, then pointers to connectivity and tags
  call iMesh_initEntArrIter(%val(imesh), %val(root_set), %val(iBase_FACE), %val(iMesh_ALL_TOPOLOGIES), &
       %val(nquads), %val(0), qiter, ierr)
  CHECK(ierr)
  call iMesh_connectIterate(%val(imesh), %val(qiter), conn_ptr, vpere, count, ierr)
  CHECK(ierr)
  call c_f_pointer(conn_ptr, conn, [vpere*nquads])
  call iMesh_tagIterate(%val(imesh), %val(tag1h), %val(qiter), tag1_ptr, count, ierr)
  CHECK(ierr)
  call c_f_pointer(tag1_ptr, tag1, [3*nquads])
  call iMesh_tagIterate(%val(imesh), %val(tag2h), %val(qiter), tag2_ptr, count, ierr)
  CHECK(ierr)
  call c_f_pointer(tag2_ptr, tag2, [3*nquads])
  
  ! iterate over elements, computing tag1 from coords positions;
  ! use (1+... in array indices for 1-based indexing
  do i = 0, nquads-1
     tag1(1+3*i+0) = 0.0
     tag1(1+3*i+1) = 0.0
     tag1(1+3*i+2) = 0.0 ! initialize position for this element
     do j = 0, vpere-1 ! loop over vertices in this quad
        v_ind = conn(1+vpere*i+j) ! assign to v_ind to allow handle arithmetic
        v_ind = v_ind - startv ! vert index is just the offset from start vertex
        tag1(1+3*i+0) = tag1(1+3*i+0) + x(1+v_ind)
        tag1(1+3*i+1) = tag1(1+3*i+1) + y(1+v_ind) ! sum vertex positions into tag1...
        tag1(1+3*i+2) = tag1(1+3*i+2) + z(1+v_ind)
     end do
     tag1(1+3*i+0) = tag1(1+3*i+0) / vpere;
     tag1(1+3*i+1) = tag1(1+3*i+1) / vpere  ! then normalize
     tag1(1+3*i+2) = tag1(1+3*i+2) / vpere;
  end do ! loop over elements in chunk
    
  ! Iterate through vertices, summing positions into tag2 on connected elements;
  ! get adjacencies all at once,
  ! could also get them per vertex (time/memory tradeoff)
  tmp_quads_size = 0
  offsets_size = 0
  tmpflag = iMesh_QUADRILATERAL
  call iMesh_getEntArrAdj(%val(imesh), verts, %val(nverts), %val(tmpflag), &
       tmp_quads_ptr, tmp_quads_size, tmp_quads_size, &
       offsets_ptr, offsets_size, offsets_size, ierr)
  CHECK(ierr)

  call c_f_pointer(tmp_quads_ptr, tmp_quads, [tmp_quads_size])
  call c_f_pointer(offsets_ptr, offsets, [offsets_size])
  call iMesh_freeMemory(%val(imesh), verts_ptr)
  do v = 0, nverts-1
     do e = 1+offsets(1+v), 1+offsets(1+v+1)-1  ! 1+ because offsets are in terms of 0-based arrays
        e_ind = tmp_quads(e) ! assign to e_ind to allow handle arithmetic
        e_ind = e_ind - startq ! e_ind is the quad handle minus the starting quad handle
        tag2(1+3*e_ind+0) = tag2(1+3*e_ind+0) + x(1+v)
        tag2(1+3*e_ind+1) = tag2(1+3*e_ind+1) + y(1+v)   ! tag on each element is 3 doubles, x/y/z
        tag2(1+3*e_ind+2) = tag2(1+3*e_ind+2) + z(1+v)
     end do
  end do
  call iMesh_freeMemory(%val(imesh), tmp_quads_ptr)
  call iMesh_freeMemory(%val(imesh), offsets_ptr)

  ! Normalize tag2 by vertex count (vpere); 
  ! loop over elements using same approach as before
  ! At the same time, compare values of tag1 and tag2
  n_dis = 0
  do e = 0, nquads-1
     do j = 0, 2
        tag2(1+3*e+j) = tag2(1+3*e+j) / vpere  ! normalize by # verts
     end do
     if (tag1(1+3*e) .ne. tag2(1+3*e) .or. tag1(1+3*e+1) .ne. tag2(1+3*e+1) .or. tag1(1+3*e+2) .ne. tag2(1+3*e+2)) then
        print *, "Tag1, tag2 disagree for element ", e
        print *, "tag1: ", tag1(1+3*e), tag1(1+3*e+1), tag1(1+3*e+2)
        print *, "tag2: ", tag2(1+3*e), tag2(1+3*e+1), tag2(1+3*e+2)
        n_dis = n_dis + 1
     endif
  end do
  if (n_dis .eq. 0) then
    print *, 'All tags agree, success!'
  endif

    ! Ok, we are done, shut down MOAB
  call iMesh_dtor(%val(imesh), ierr)
  CHECK(ierr)
end program DirectAccessNoHolesF90

subroutine create_mesh_no_holes(imesh, nquads, ierr) 
  use, intrinsic :: ISO_C_BINDING
  implicit none
#include "iMesh_f.h"

  iMesh_Instance imesh
  integer nquads, ierr

  real(C_DOUBLE), pointer :: coords(:,:)
  TYPE(C_PTR) tmp_ents_ptr, stat_ptr
  iBase_EntityHandle, pointer :: connect(:), tmp_ents(:)
  integer nverts, tmp_size, stat_size, i

  ! first make the mesh, a 1d array of quads with left hand side x = elem_num;
  ! vertices are numbered in layers
  ! create verts, num is 2(nquads+1) because they're in a 1d row
  nverts = 2*(nquads+1)
  allocate(coords(0:2, 0:nverts-1))
  do i = 0, nquads-1
     coords(0,2*i) = i
     coords(0,2*i+1) = i ! x values are all i
     coords(1,2*i) = 0.0
     coords(1,2*i+1) = 1.0 ! y coords
     coords(2,2*i) = 0.0
     coords(2,2*i+1) = 0.0 ! z values, all zero (2d mesh)
  end do
  ! last two vertices
  coords(0, nverts-2) = nquads
  coords(0, nverts-1) = nquads
  coords(1, nverts-2) = 0.0
  coords(1, nverts-1) = 1.0 ! y coords
  coords(2, nverts-2) = 0.0
  coords(2, nverts-1) = 0.0 ! z values, all zero (2d mesh)
  tmp_size = 0
  !!!!!! VSM: This is the culprit call that screws up IBM xlf compiler
  call iMesh_createVtxArr(%val(imesh), %val(nverts), %val(iBase_INTERLEAVED), &
       coords, %val(3*nverts), tmp_ents_ptr, tmp_size, tmp_size, ierr)
  CHECK(ierr)
  call c_f_pointer(tmp_ents_ptr, tmp_ents, [nverts])
  deallocate(coords)
  allocate(connect(0:4*nquads-1))
  do i = 0, nquads-1
     connect(4*i+0) = tmp_ents(1+2*i)
     connect(4*i+1) = tmp_ents(1+2*i+2)
     connect(4*i+2) = tmp_ents(1+2*i+3)
     connect(4*i+3) = tmp_ents(1+2*i+1)
  end do
  stat_size = 0
  stat_ptr = C_NULL_PTR
  ! re-use tmp_ents here;
  ! we know it'll always be larger than nquads so it's ok
  call iMesh_createEntArr(%val(imesh), %val(iMesh_QUADRILATERAL), connect, %val(4*nquads), &
       tmp_ents_ptr, tmp_size, tmp_size, stat_ptr, stat_size, stat_size, ierr)
  CHECK(ierr)

  ierr = iBase_SUCCESS

  call iMesh_freeMemory(%val(imesh), tmp_ents_ptr)
  call iMesh_freeMemory(%val(imesh), stat_ptr)
  deallocate(connect)

  return
end subroutine create_mesh_no_holes

