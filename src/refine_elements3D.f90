!---------------------------------------------------------------------!
!                                PHAML                                !
!                                                                     !
! The Parallel Hierarchical Adaptive MultiLevel code for solving      !
! linear elliptic partial differential equations of the form          !
! (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
! boundary conditions, and eigenvalue problems where F is lambda*U.   !
!                                                                     !
! PHAML is public domain software.  It was produced as part of work   !
! done by the U.S. Government, and is not subject to copyright in     !
! the United States.                                                  !
!                                                                     !
!     William F. Mitchell                                             !
!     Applied and Computational Mathematics Division                  !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!---------------------------------------------------------------------!

module refine_elements

!----------------------------------------------------
! This module contains routines with the details of element refinement.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use linsystype_mod
use boundary_util
use hash_mod
use grid_util
use error_estimators
use make_linsys
use message_passing, only: fatal, warning
use omp_lib, only: omp_get_num_threads
!----------------------------------------------------

implicit none
private
public p_coarsen_elem, &            ! not thread safe
       p_coarsen_elem_interior, &   ! conditionally thread safe
       h_coarsen_element, &         ! conditionally thread safe
       p_refine_elem, &             ! not thread safe
       p_refine_element_interior, & ! conditionally thread safe
       enforce_edge_rule, &         ! conditionally thread safe
       before_h_refine, &           ! not thread safe
       h_refine_element, &          ! conditionally thread safe
       after_h_refine, &            ! not thread safe
       remove_from_errind_list, &   ! not thread safe
       create_element, &            ! not thread safe
       init_guess_ic, &             ! not thread safe
       init_guess_p, &              ! not thread safe
       remove_hanging_nodes         ! not thread safe

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          --------------
subroutine p_coarsen_elem(grid,elem,errcode,refine_control,delay_errind)
!          --------------

!----------------------------------------------------
! This routine reduces the degree of element elem by 1.
! errcode is 0 if successful
!            1 if an error occurred
!           -1 if cannot coarsen
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
integer, intent(out) :: errcode
type(refine_options), intent(in) :: refine_control
logical, intent(in), optional :: delay_errind
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

errcode = 0

call fatal("p_coarsen_elem not yet written for 3D")
stop
end subroutine p_coarsen_elem

!          -----------------------
subroutine p_coarsen_elem_interior(grid,elem,errcode,refine_control, &
                                   delta_dof,delta_dof_own,delay_errind)
!          -----------------------

!----------------------------------------------------
! This routine reduces the degree of element elem by 1, but does not change
! the edge degrees.
! errcode is 0 if successful
!            1 if an error occurred
!           -1 if cannot coarsen
!
! Thread safe if delta_dof and delta_dof_own are present, and the error
! estimator is not equilibrated residual if delay_errind is absent or false.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
integer, intent(out) :: errcode
type(refine_options), intent(in) :: refine_control
integer, intent(out), optional :: delta_dof, delta_dof_own
logical, intent(in), optional :: delay_errind
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

errcode = 0

call fatal("p_coarsen_elem_interior not yet written for 3D")
stop
end subroutine p_coarsen_elem_interior

!          ----------------------
subroutine h_coarsen_element(grid,parent,errcode,refcont, &
                                  delta_nelem, delta_nelem_leaf, delta_nedge, &
                                  delta_nvert, delta_nelem_leaf_own, &
                                  delta_nvert_own, delta_dof, delta_dof_own, &
                                  delay_errind)
!          ----------------------

!----------------------------------------------------
! This routine removes the bisection refinement of element parent and
! its mate.
! errcode is 0 if successful
!            1 if an error occurred
!           -1 if cannot unrefine because of grandchildren
!           -2 if nothing to derefine (no children)
!           -3 if cannot unrefine because children have different owners
!
! Thread safe if:
!  -parent and it's mate are not both among the elements to be done in parallel
!  -the deltas are present
!
! If the deltas are present, the corresponding scalars are not updated in
! grid, the edge degrees and solutions are not updated (enforce edge rule
! later) and the initial guess for the solution is not set.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent
integer, intent(out) :: errcode
type(refine_options), intent(in) :: refcont
integer, intent(out), optional :: delta_nelem, delta_nelem_leaf, &
                                  delta_nedge, delta_nvert, &
                                  delta_nelem_leaf_own, delta_nvert_own, &
                                  delta_dof, delta_dof_own 
logical, intent(in), optional :: delay_errind

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

errcode = 0

call fatal("h_coarsen_element not yet written for 3D")
stop
end subroutine h_coarsen_element

!          -------------
subroutine p_refine_elem(grid,elem,refine_control,elist,numpref, &
                         return_to_elist,errcode)
!          -------------

!----------------------------------------------------
! This routine performs p refinement of element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
type(refine_options), intent(in) :: refine_control
type(errind_list), optional, intent(inout) :: elist
integer, optional, intent(inout) :: numpref(:)
logical, optional, intent(in) :: return_to_elist
integer, optional, intent(out) :: errcode
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! if the error indicator list is provided, remove the element from it

errcode = 0

call fatal("p_refine_elem not yet written for 3D")
stop
end subroutine p_refine_elem

!          -------------------------
subroutine p_refine_element_interior(grid,refine_control,elem,errcode, &
                                     delta_dof,delta_dof_own)
!          -------------------------

!----------------------------------------------------
! This routine performs p refinement of the interior of element elem,
! but not the edges.
!
! Thread safe if delta_dof and delta_dof_own are present.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: elem
integer, intent(out) :: errcode
integer, intent(out), optional :: delta_dof, delta_dof_own
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

errcode = 0

call fatal("p_refine_element_interior not yet written for 3D")
stop
end subroutine p_refine_element_interior

!          -----------------
subroutine enforce_edge_rule(grid,refine_control,edge,delta_dof,delta_dof_own)
!          -----------------

!----------------------------------------------------
! This routine enforces either the minimum or maximum rule for edge
!
! Thread safe if:
!   delta_dof and delta_dof_own are present
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: edge
integer, intent(out), optional :: delta_dof, delta_dof_own
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call fatal("enforce_edge_rule not yet written for 3D")
stop
end subroutine enforce_edge_rule

!          ---------------
subroutine before_h_refine(grid,element_list,nelem,elist,reftype,new_p,numpref,&
                           numhref,vert_lid,edge_lid,elem_lid,desired_level, &
                           desired_degree,elem_list)
!          ---------------

!----------------------------------------------------
! This routine performs parts of h refinement that must be done by a single
! OpenMP thread before the OpenMP-parallel refinement.
! All elements in element_list must be on the same level, and there must only
! be one from a compatibly divisible pair.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: element_list(:), nelem
type(errind_list), optional, intent(inout) :: elist
character(len=*), optional, pointer :: reftype(:)
integer, optional, pointer :: numhref(:), numpref(:), new_p(:,:)
integer, intent(out) :: vert_lid(:,:), edge_lid(:,:), elem_lid(:,:)
! second dim >= nelem; first dim 2              6              4
integer, optional, pointer :: desired_level(:), desired_degree(:), elem_list(:)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

vert_lid=0; edge_lid=0; elem_lid=0

! TEMP110606 don't stop
return

call fatal("before_h_refine not yet written for 3D")
stop
end subroutine before_h_refine

!          ----------------
subroutine h_refine_element(grid,parent,errcode,did_refinement_edge,refcont, &
                     solver_control,elist,reftype,new_p,numhref,numpref, &
                     return_to_elist,reduce_p,reduce_p_max,vert_child_lid, &
                     edge_child_lid,elem_child_lid,delta_dof,delta_dof_own, &
                     delta_nelem,delta_nelem_leaf,delta_nelem_leaf_own, &
                     delta_nedge,delta_nedge_own,delta_nvert,delta_nvert_own, &
                     max_nlev,delay_errind,desired_level,desired_degree, &
                     elem_list)
!          ----------------

!----------------------------------------------------
! This routine refines triangle parent by bisection.
! errcode = 0 for success, 1 if the grid is full, 2 if too many levels,
! 3 if the element is not a leaf, 4 if it cannot be refined for another reason
!
! NOT THREAD SAFE, but ultimately
! thread safe if:
!   all elements being refined in parallel have the same h-level and are
!     compatibly divisible
!   an element and its mate are not both in the list of elements being refined
!   elist is not present
!   vert_child_lid, edge_child_lid, elem_child_lid are present
!     (in this case, the tasks performed in before_h_refine and after_h_refine
!     are not performed here)
!   delta_* and max_nlev are present
!     (in this case, the grid scalars are not updated)
!   the solution at preexisting vertices does not change
!   if numhref is present, it is at most 1 for parent and mate
!   if numpref is present, it is 0 for parent and mate
!   if return_to_elist is present, it is false
!   TEMP reduce_p is not present
!   TEMP no periodic boundary conditions
!   new_p is 0 or not present, i.e. refinement is not steepest slope hp-adaptive
!   some other conditions that I think hold; search for OPENMP
!
! Will need to remove desired_level, desired_degree and elem_list when making
! it thread safe.  They are only here to be reallocated if grid%element is
! reallocated, and that will eventually happen in before_h_refine.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent
integer, intent(inout) :: errcode
integer, intent(out) :: did_refinement_edge
type(refine_options), intent(in) :: refcont
type(solver_options), intent(in) :: solver_control
type(errind_list), optional, intent(inout) :: elist
character(len=*), optional, pointer :: reftype(:)
integer, optional, pointer :: new_p(:,:), numhref(:), numpref(:)
logical, optional, intent(in) :: return_to_elist
integer, optional, intent(in) :: reduce_p
integer, optional, intent(inout) :: reduce_p_max
integer, optional, intent(in) :: vert_child_lid(:), edge_child_lid(:), &
                                 elem_child_lid(:)
integer, optional, intent(out) :: delta_dof, delta_dof_own, delta_nelem, &
                                  delta_nelem_leaf, delta_nelem_leaf_own, &
                                  delta_nedge, delta_nedge_own, delta_nvert, &
                                  delta_nvert_own, max_nlev
logical, optional, intent(in) :: delay_errind
integer, optional, pointer :: desired_level(:), desired_degree(:), elem_list(:)

!----------------------------------------------------
! Local variables:

integer :: refinement_edge, refinement_face, child1, child2, face1, face2, &
           face3, face4, face5, edge1, edge2, edge3, edge4, vert1, i, j, astat,&
           nfound, parent_face(2), vert, edge, face, &
           child_vert(2), common_vert(2), child_face(2), child_edge(4), &
           common_edge, jerr, mate_face_1, mate_face_2, mate_face_3, &
           mate_face_4, mate_edge_1, mate_edge_2, mate_edge_3, mate_edge_4, &
           mate_vert_1, itemp, last, lev
integer :: bctype(grid%system_size)
real(my_real) :: bccoef(grid%system_size,grid%system_size), &
                 bcrhs(grid%system_size)
type(hash_key) :: parent_neighbors(NEIGHBORS_PER_ELEMENT)
logical :: new_ref_edge, new_face1, new_face2, do_face1, do_face2, &
           mate_edge_1_in_cycle, edge1_in_cycle, mate_vert_1_in_cycle, &
           vert1_in_cycle
!----------------------------------------------------
! Begin executable code

! check for things not yet supported

if (omp_get_num_threads() /= 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: OpenMP not yet supported in 3D")
   stop
endif

! some preliminary stuff

did_refinement_edge = -1
errcode = 0
if (present(delta_dof)) then
   delta_dof = 0
   delta_dof_own = 0
   delta_nelem = 0
   delta_nelem_leaf = 0
   delta_nelem_leaf_own = 0
   delta_nedge = 0
   delta_nedge_own = 0
   delta_nvert = 0
   delta_nvert_own = 0
   max_nlev = 0
endif
if (present(elist)) then
   print *,"h_refine_element didn't expect elist to be present"
   stop
endif

! can't refine a non-leaf element

child1 = get_child_lid(grid%element(parent)%gid,1,grid%elem_hash)
if (child1 /= NO_CHILD) then
   errcode = 3
   return
endif

! get the refinement edge

refinement_edge = grid%element(parent)%refinement_edge

! see if refinement of this element requires extending the number of levels.

!!!!if (.not. present(elem_child_lid)) then
   if (grid%element(parent)%level >= size(grid%head_level_elem)) then
      call extend_nlev(grid)
      if (ierr /= NO_ERROR) return
   endif
!!!!endif

! get the neighbors of the parent in the order of the face they share

parent_neighbors = get_ordered_neighbors(parent,jerr)

! If jerr is 1, we didn't find all matching faces; probably a neighbor has not
! been refined yet.  Those not found were assigned BOUNDARY.  parent_neighbors
! is only used to assign neighbor_hint, so use the parent's hint as the element
! that shares the faces not matched.  If the face really is BOUNDARY, then the
! parent's hint will be BOUNDARY, so that is not lost.

if (jerr == 1) then
   do i=1,FACES_PER_ELEMENT
      if (parent_neighbors(i) == BOUNDARY) then
         parent_neighbors(i) = grid%element(parent)%neighbor_hint(i)
      endif
   end do
endif

! get the lids for the first two children

child1 = get_next_free_elem(grid,errcode,elist,reftype,new_p,numhref, &
                            numpref,desired_level,desired_degree, &
                            elem_list)
if (errcode /= 0) return
child2 = get_next_free_elem(grid,errcode,elist,reftype,new_p,numhref, &
                            numpref,desired_level,desired_degree, &
                            elem_list)
if (errcode /= 0) return

! if the refinement edge has already been refined then locate the lids of
! the entities that were created by that bisection, otherwise get new lids

new_ref_edge = grid%edge(refinement_edge)%child(1) == NO_CHILD
if (new_ref_edge) then

! refinement edge is not refined, get new lids for the new vertex and edges
! that replace the refinement edge

   edge1 = get_next_free_edge(grid,errcode)
   if (errcode /= 0) return
   edge2 = get_next_free_edge(grid,errcode)
   if (errcode /= 0) return
   vert1 = get_next_free_vert(grid,errcode)
   if (errcode /= 0) return
   did_refinement_edge = refinement_edge

! refinement edge is refined, identify the lids

else

   if (grid%edge(grid%edge(refinement_edge)%child(1))%vertex(1) == &
       grid%edge(refinement_edge)%vertex(1)) then
      edge1 = grid%edge(refinement_edge)%child(1)
      edge2 = grid%edge(refinement_edge)%child(2)
      vert1 = grid%edge(edge1)%vertex(2)
   elseif (grid%edge(grid%edge(refinement_edge)%child(1))%vertex(2) == &
       grid%edge(refinement_edge)%vertex(1)) then
      edge1 = grid%edge(refinement_edge)%child(1)
      edge2 = grid%edge(refinement_edge)%child(2)
      vert1 = grid%edge(edge1)%vertex(1)
   elseif (grid%edge(grid%edge(refinement_edge)%child(2))%vertex(1) == &
       grid%edge(refinement_edge)%vertex(1)) then
      edge1 = grid%edge(refinement_edge)%child(2)
      edge2 = grid%edge(refinement_edge)%child(1)
      vert1 = grid%edge(edge1)%vertex(2)
   elseif (grid%edge(grid%edge(refinement_edge)%child(2))%vertex(2) == &
       grid%edge(refinement_edge)%vertex(1)) then
      edge1 = grid%edge(refinement_edge)%child(2)
      edge2 = grid%edge(refinement_edge)%child(1)
      vert1 = grid%edge(edge1)%vertex(1)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("h_refine_element: couldn't identify lids in refined refinement edge")
      stop
   endif

endif

! determine the two faces of the parent that contain the refinement edge

nfound = 0
do i=1,FACES_PER_ELEMENT
   if (any(grid%face(grid%element(parent)%face(i))%edge==refinement_edge)) then
      nfound = nfound + 1
      parent_face(nfound) = grid%element(parent)%face(i)
   endif
end do
if (nfound /= 2) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: didn't find two faces with refinement edge")
   stop
endif

! for each of those two faces, if the face has already been refined then
! locate the lids of the entities that were created by that bisection,
! otherwise get new lids

new_face1 = grid%face(parent_face(1))%child(1) == NO_CHILD
if (new_face1) then


   face1 = get_next_free_face(grid,errcode)
   if (errcode /= 0) return
   face2 = get_next_free_face(grid,errcode)
   if (errcode /= 0) return
   edge3 = get_next_free_edge(grid,errcode)
   if (errcode /= 0) return

else

   if (any(grid%face(grid%face(parent_face(1))%child(1))%edge==edge1)) then
      face1 = grid%face(parent_face(1))%child(1)
      face2 = grid%face(parent_face(1))%child(2)
   elseif (any(grid%face(grid%face(parent_face(1))%child(1))%edge==edge2)) then
      face1 = grid%face(parent_face(1))%child(2)
      face2 = grid%face(parent_face(1))%child(1)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("h_refine_element: failed to find refined edge in parent face 1")
      stop
   endif
   if (grid%face(face1)%edge(1) /= edge1 .and. &
       grid%face(face1)%edge(1) /= grid%face(parent_face(1))%edge(1) .and. &
       grid%face(face1)%edge(1) /= grid%face(parent_face(1))%edge(2) .and. &
       grid%face(face1)%edge(1) /= grid%face(parent_face(1))%edge(3)) then
      edge3 = grid%face(face1)%edge(1)
   elseif (grid%face(face1)%edge(2) /= edge1 .and. &
           grid%face(face1)%edge(2) /= grid%face(parent_face(1))%edge(1) .and. &
           grid%face(face1)%edge(2) /= grid%face(parent_face(1))%edge(2) .and. &
           grid%face(face1)%edge(2) /= grid%face(parent_face(1))%edge(3)) then
      edge3 = grid%face(face1)%edge(2)
   elseif (grid%face(face1)%edge(3) /= edge1 .and. &
           grid%face(face1)%edge(3) /= grid%face(parent_face(1))%edge(1) .and. &
           grid%face(face1)%edge(3) /= grid%face(parent_face(1))%edge(2) .and. &
           grid%face(face1)%edge(3) /= grid%face(parent_face(1))%edge(3)) then
      edge3 = grid%face(face1)%edge(3)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("h_refine_element: failed to find edge1 on face1")
      stop
   endif

endif

new_face2 = grid%face(parent_face(2))%child(1) == NO_CHILD
if (new_face2) then

   face3 = get_next_free_face(grid,errcode)
   if (errcode /= 0) return
   face4 = get_next_free_face(grid,errcode)
   if (errcode /= 0) return
   edge4 = get_next_free_edge(grid,errcode)
   if (errcode /= 0) return

else

   if (any(grid%face(grid%face(parent_face(2))%child(1))%edge==edge1)) then
      face3 = grid%face(parent_face(2))%child(1)
      face4 = grid%face(parent_face(2))%child(2)
   elseif (any(grid%face(grid%face(parent_face(2))%child(1))%edge==edge2)) then
      face3 = grid%face(parent_face(2))%child(2)
      face4 = grid%face(parent_face(2))%child(1)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("h_refine_element: failed to find refined edge in parent face 2")
      stop
   endif
   if (grid%face(face3)%edge(1) /= edge1 .and. &
       grid%face(face3)%edge(1) /= grid%face(parent_face(2))%edge(1) .and. &
       grid%face(face3)%edge(1) /= grid%face(parent_face(2))%edge(2) .and. &
       grid%face(face3)%edge(1) /= grid%face(parent_face(2))%edge(3)) then
      edge4 = grid%face(face3)%edge(1)
   elseif (grid%face(face3)%edge(2) /= edge1 .and. &
           grid%face(face3)%edge(2) /= grid%face(parent_face(2))%edge(1) .and. &
           grid%face(face3)%edge(2) /= grid%face(parent_face(2))%edge(2) .and. &
           grid%face(face3)%edge(2) /= grid%face(parent_face(2))%edge(3)) then
      edge4 = grid%face(face3)%edge(2)
   elseif (grid%face(face3)%edge(3) /= edge1 .and. &
           grid%face(face3)%edge(3) /= grid%face(parent_face(2))%edge(1) .and. &
           grid%face(face3)%edge(3) /= grid%face(parent_face(2))%edge(2) .and. &
           grid%face(face3)%edge(3) /= grid%face(parent_face(2))%edge(3)) then
      edge4 = grid%face(face3)%edge(3)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("h_refine_element: failed to find edge1 on face1")
      stop
   endif

endif

! identify which vertices from the parent go in each child or are common

child_vert(1) = grid%edge(refinement_edge)%vertex(1)
child_vert(2) = grid%edge(refinement_edge)%vertex(2)
nfound = 0
do i=1,VERTICES_PER_ELEMENT
   vert = grid%element(parent)%vertex(i)
   if (vert /= child_vert(1) .and. vert /= child_vert(2)) then
      nfound = nfound + 1
      common_vert(nfound) = vert
   endif
end do
if (nfound /= 2) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: didn't find 2 common vertices")
   stop
endif

! make sure the first common vertex is on the first face

if (any(grid%face(parent_face(2))%vertex == common_vert(1))) then
   i = common_vert(1)
   common_vert(1) = common_vert(2)
   common_vert(2) = i
endif

! determine which faces from the parent go in each child

nfound = 0
do i=1,FACES_PER_ELEMENT
   face = grid%element(parent)%face(i)
   if (face /= parent_face(1) .and. face /= parent_face(2)) then
      if (any(grid%face(face)%vertex == child_vert(1))) then
         nfound = nfound + 1
         child_face(1) = face
      endif
      if (any(grid%face(face)%vertex == child_vert(2))) then
         nfound = nfound + 1
         child_face(2) = face
      endif
   endif
end do
if (nfound /= 2) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: didn't find 2 child faces")
   stop
endif

! determine which edges from the parent go in each child and which is common

common_edge = -1
nfound = 0
do i=1,EDGES_PER_ELEMENT
   edge = grid%element(parent)%edge(i)
   if (edge == refinement_edge) cycle
   if (any(grid%edge(edge)%vertex == child_vert(1))) then
      nfound = nfound + 1
      child_edge(nfound) = edge
   elseif (any(grid%edge(edge)%vertex == common_vert(1)) .and. &
           any(grid%edge(edge)%vertex == common_vert(2))) then
      common_edge = edge
   endif
end do
do i=1,EDGES_PER_ELEMENT
   edge = grid%element(parent)%edge(i)
   if (edge == refinement_edge) cycle
   if (any(grid%edge(edge)%vertex == child_vert(2))) then
      nfound = nfound + 1
      child_edge(nfound) = edge
   endif
end do
if (nfound /= 4) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: didn't identify which edges go in which child")
   stop
endif
if (common_edge == -1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: didn't identify common edge for children")
   stop
endif

! make sure child edges 1 and 3 are on the first face

if (any(grid%face(parent_face(2))%edge == child_edge(1))) then
   i = child_edge(1)
   child_edge(1) = child_edge(2)
   child_edge(2) = i
endif
if (any(grid%face(parent_face(2))%edge == child_edge(3))) then
   i = child_edge(3)
   child_edge(3) = child_edge(4)
   child_edge(4) = i
endif

! define the new vertex and edges from bisecting the refinement edge, if they
! didn't already exist

if (new_ref_edge) then

   grid%vertex(vert1)%gid = 2*grid%element(parent)%gid
   grid%vertex(vert1)%assoc_elem = child1
   nullify(grid%vertex(vert1)%boundary_vertex)
   if (grid%edge(refinement_edge)%bmark == 0) then ! vert1 is interior
      grid%vertex(vert1)%coord%x = &
               (grid%vertex(grid%edge(refinement_edge)%vertex(1))%coord%x + &
                grid%vertex(grid%edge(refinement_edge)%vertex(2))%coord%x)/2
      grid%vertex(vert1)%coord%y = &
               (grid%vertex(grid%edge(refinement_edge)%vertex(1))%coord%y + &
                grid%vertex(grid%edge(refinement_edge)%vertex(2))%coord%y)/2
      grid%vertex(vert1)%coord%z = &
               (grid%vertex(grid%edge(refinement_edge)%vertex(1))%coord%z + &
                grid%vertex(grid%edge(refinement_edge)%vertex(2))%coord%z)/2
   else ! vert1 is a new boundary vertex
      call define_boundary_vertex(grid,refinement_edge,vert1)
   endif
   grid%vertex(vert1)%bmark = grid%edge(refinement_edge)%bmark
   grid%vertex(vert1)%tags = grid%edge(refinement_edge)%tags
   grid%vertex_type(vert1,:) = grid%edge_type(refinement_edge,:)
   if (any(grid%vertex_type(vert1,:) == DIRICHLET) .or. &
       any(grid%vertex_type(vert1,:) == PERIODIC_MASTER_DIR) .or. &
       any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE_DIR)) then
      call my_bconds(grid%vertex(vert1)%coord%x,grid%vertex(vert1)%coord%y, &
                     grid%vertex(vert1)%coord%z,grid%vertex(vert1)%bmark, &
                     bctype,bccoef,bcrhs)
   endif
   do i=1,grid%system_size
      if (grid%vertex_type(vert1,i) == DIRICHLET .or. &
          grid%vertex_type(vert1,i) == PERIODIC_MASTER_DIR .or. &
          grid%vertex_type(vert1,i) == PERIODIC_SLAVE_DIR) then
         grid%vertex_solution(vert1,i,:) = bcrhs(i)
      else
         grid%vertex_solution(vert1,i,:) = &
            (grid%vertex_solution(grid%edge(refinement_edge)%vertex(1),i,:) + &
             grid%vertex_solution(grid%edge(refinement_edge)%vertex(2),i,:))/2
      endif
   end do
   if (grid%have_true) then
      do j=1,max(1,grid%num_eval)
         do i=1,grid%system_size
            grid%vertex_exact(vert1,i,j)=my_trues(grid%vertex(vert1)%coord%x, &
                                                  grid%vertex(vert1)%coord%y, &
                                                  grid%vertex(vert1)%coord%z, &
                                                  i,j)
         end do
      end do
   endif
   grid%vertex(vert1)%next = vert1

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%vertex(vert1)%gid,vert1,grid%vert_hash)
!!!!!endif

! edge 1

   grid%edge(edge1)%gid = 4*grid%element(parent)%gid
   grid%edge(edge1)%vertex = (/ child_vert(1),vert1 /)
   grid%edge(edge1)%bmark = grid%edge(refinement_edge)%bmark
   grid%edge(edge1)%tags = grid%edge(refinement_edge)%tags
   grid%edge(edge1)%degree = grid%edge(refinement_edge)%degree
   grid%edge_type(edge1,:) = grid%edge_type(refinement_edge,:)
   grid%edge(edge1)%assoc_elem = child1
   grid%edge(edge1)%child = NO_CHILD
   if (grid%edge(edge1)%degree > 1) then
      allocate(grid%edge(edge1)%solution(grid%edge(edge1)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%edge(edge1)%exact(grid%edge(edge1)%degree-1, &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%edge(edge1)%exact = 0.0_my_real
      else
         nullify(grid%edge(edge1)%exact)
      endif
      do i=1,grid%system_size
         if (grid%edge_type(edge1,i) == DIRICHLET .or. &
             grid%edge_type(edge1,i) == PERIODIC_MASTER_DIR .or. &
             grid%edge_type(edge1,i) == PERIODIC_SLAVE_DIR) then
            call edge_exact(grid,edge1,i,"d")
         else
            grid%edge(edge1)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call edge_exact(grid,edge1,i,"t")
      end do
   else
      nullify(grid%edge(edge1)%solution)
      nullify(grid%edge(edge1)%exact)
   endif
   nullify(grid%edge(edge1)%oldsoln)
   grid%edge(edge1)%next = edge1

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%edge(edge1)%gid,edge1,grid%edge_hash)
!!!!!endif

! edge 2

   grid%edge(edge2)%gid = 4*grid%element(parent)%gid + 1
   grid%edge(edge2)%vertex = (/ child_vert(2),vert1 /)
   grid%edge(edge2)%bmark = grid%edge(refinement_edge)%bmark
   grid%edge(edge2)%tags = grid%edge(refinement_edge)%tags
   grid%edge(edge2)%degree = grid%edge(refinement_edge)%degree
   grid%edge_type(edge2,:) = grid%edge_type(refinement_edge,:)
   grid%edge(edge2)%assoc_elem = child2
   grid%edge(edge2)%child = NO_CHILD
   if (grid%edge(edge2)%degree > 1) then
      allocate(grid%edge(edge2)%solution(grid%edge(edge2)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%edge(edge2)%exact(grid%edge(edge2)%degree-1, &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%edge(edge2)%exact = 0.0_my_real
      else
         nullify(grid%edge(edge2)%exact)
      endif
      do i=1,grid%system_size
         if (grid%edge_type(edge2,i) == DIRICHLET .or. &
             grid%edge_type(edge2,i) == PERIODIC_MASTER_DIR .or. &
             grid%edge_type(edge2,i) == PERIODIC_SLAVE_DIR) then
            call edge_exact(grid,edge2,i,"d")
         else
            grid%edge(edge2)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call edge_exact(grid,edge2,i,"t")
      end do
   else
      nullify(grid%edge(edge2)%solution)
      nullify(grid%edge(edge2)%exact)
   endif
   nullify(grid%edge(edge2)%oldsoln)
   grid%edge(edge2)%next = edge2

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%edge(edge2)%gid,edge2,grid%edge_hash)
!!!!!endif

endif ! refinement edge is new

! define the new edge and faces from bisecting parent face 1, if they
! didn't already exist

if (new_face1) then

   refinement_face = parent_face(1)
   grid%edge(edge3)%gid = 4*grid%element(parent)%gid + 2
   grid%edge(edge3)%vertex = (/ common_vert(1),vert1 /)
   grid%edge(edge3)%bmark = grid%face(refinement_face)%bmark
   grid%edge(edge3)%tags = grid%face(refinement_face)%tags
   grid%edge(edge3)%degree = grid%face(refinement_face)%degree
   grid%edge_type(edge3,:) = grid%face_type(refinement_face,:)
   grid%edge(edge3)%assoc_elem = child1
   grid%edge(edge3)%child = NO_CHILD
   if (grid%edge(edge3)%degree > 1) then
      allocate(grid%edge(edge3)%solution(grid%edge(edge3)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%edge(edge3)%exact(grid%edge(edge3)%degree-1, &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%edge(edge3)%exact = 0.0_my_real
      else
         nullify(grid%edge(edge3)%exact)
      endif
      do i=1,grid%system_size
         if (grid%edge_type(edge3,i) == DIRICHLET) then
            call edge_exact(grid,edge3,i,"d")
         else
            grid%edge(edge3)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call edge_exact(grid,edge3,i,"t")
      end do
   else
      nullify(grid%edge(edge3)%solution)
      nullify(grid%edge(edge3)%exact)
   endif
   nullify(grid%edge(edge3)%oldsoln)
   grid%edge(edge3)%next = edge3

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%edge(edge3)%gid,edge3,grid%edge_hash)
!!!!!endif

! face 1

   grid%face(face1)%gid = 4*grid%face(parent_face(1))%gid
   grid%face(face1)%vertex = (/ child_vert(1), &
                                common_vert(1), &
                                vert1 /)
   grid%face(face1)%edge = (/ child_edge(1), &
                              edge1, edge3 /)
   grid%face(face1)%bmark = grid%face(refinement_face)%bmark
   grid%face(face1)%tags = grid%face(refinement_face)%tags
   grid%face(face1)%degree = grid%face(refinement_face)%degree
   grid%face_type(face1,:) = grid%face_type(refinement_face,:)
   grid%face(face1)%assoc_elem = child1
   grid%face(face1)%marked_edge = grid%face(face1)%edge(1)
   grid%face(face1)%child = NO_CHILD
   grid%face(face1)%next = 0
   if (grid%face(face1)%degree > 2) then
      allocate(grid%face(face1)%solution(face_dof(grid%face(face1)%degree), &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%face(face1)%exact(face_dof(grid%face(face1)%degree), &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%face(face1)%exact = 0.0_my_real
      else
         nullify(grid%face(face1)%exact)
      endif
      do i=1,grid%system_size
         if (grid%face_type(face1,i) == DIRICHLET) then
            call face_exact(grid,face1,i,"d")
         else
            grid%face(face1)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call face_exact(grid,face1,i,"t")
      end do
   else
      nullify(grid%face(face1)%solution)
      nullify(grid%face(face1)%exact)
   endif
   nullify(grid%face(face1)%oldsoln)
   grid%face(face1)%next = face1

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%face(face1)%gid,face1,grid%face_hash)
!!!!!endif

! face 2

   grid%face(face2)%gid = 4*grid%face(parent_face(1))%gid + 1
   grid%face(face2)%vertex = (/ child_vert(2), &
                                common_vert(1), &
                                vert1 /)
   grid%face(face2)%edge = (/ child_edge(3), &
                              edge2, edge3 /)
   grid%face(face2)%bmark = grid%face(refinement_face)%bmark
   grid%face(face2)%tags = grid%face(refinement_face)%tags
   grid%face(face2)%degree = grid%face(refinement_face)%degree
   grid%face_type(face2,:) = grid%face_type(refinement_face,:)
   grid%face(face2)%assoc_elem = child2
   grid%face(face2)%marked_edge = grid%face(face2)%edge(1)
   grid%face(face2)%child = NO_CHILD
   grid%face(face2)%next = 0
   if (grid%face(face2)%degree > 2) then
      allocate(grid%face(face2)%solution(face_dof(grid%face(face2)%degree), &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%face(face2)%exact(face_dof(grid%face(face2)%degree), &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%face(face2)%exact = 0.0_my_real
      else
         nullify(grid%face(face2)%exact)
      endif
      do i=1,grid%system_size
         if (grid%face_type(face2,i) == DIRICHLET) then
            call face_exact(grid,face2,i,"d")
         else
            grid%face(face2)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call face_exact(grid,face2,i,"t")
      end do
   else
      nullify(grid%face(face2)%solution)
      nullify(grid%face(face2)%exact)
   endif
   nullify(grid%face(face2)%oldsoln)
   grid%face(face2)%next = face2

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%face(face2)%gid,face2,grid%face_hash)
!!!!!endif

endif ! children of parent face 1 didn't exist

! define the new edge and faces from bisecting parent face 2, if they
! didn't already exist

if (new_face2) then

   refinement_face = parent_face(2)
   grid%edge(edge4)%gid = 4*grid%element(parent)%gid + 3
   grid%edge(edge4)%vertex = (/ common_vert(2),vert1 /)
   grid%edge(edge4)%bmark = grid%face(refinement_face)%bmark
   grid%edge(edge4)%tags = grid%face(refinement_face)%tags
   grid%edge(edge4)%degree = grid%face(refinement_face)%degree
   grid%edge_type(edge4,:) = grid%face_type(refinement_face,:)
   grid%edge(edge4)%assoc_elem = child1
   grid%edge(edge4)%child = NO_CHILD
   if (grid%edge(edge4)%degree > 1) then
      allocate(grid%edge(edge4)%solution(grid%edge(edge4)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%edge(edge4)%exact(grid%edge(edge4)%degree-1, &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%edge(edge4)%exact = 0.0_my_real
      else
         nullify(grid%edge(edge4)%exact)
      endif
      do i=1,grid%system_size
         if (grid%edge_type(edge4,i) == DIRICHLET) then
            call edge_exact(grid,edge4,i,"d")
         else
            grid%edge(edge4)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call edge_exact(grid,edge4,i,"t")
      end do
   else
      nullify(grid%edge(edge4)%solution)
      nullify(grid%edge(edge4)%exact)
   endif
   nullify(grid%edge(edge4)%oldsoln)
   grid%edge(edge4)%next = edge4

!!!!!if (.not. present(elem_child_lid)) then
   call hash_insert(grid%edge(edge4)%gid,edge4,grid%edge_hash)
!!!!!endif

! face 3

   grid%face(face3)%gid = 4*grid%face(parent_face(2))%gid
   grid%face(face3)%vertex = (/ child_vert(1), &
                                common_vert(2), &
                                vert1 /)
   grid%face(face3)%edge = (/ child_edge(2), &
                              edge1, edge4 /)
   grid%face(face3)%bmark = grid%face(refinement_face)%bmark
   grid%face(face3)%tags = grid%face(refinement_face)%tags
   grid%face(face3)%degree = grid%face(refinement_face)%degree
   grid%face_type(face3,:) = grid%face_type(refinement_face,:)
   grid%face(face3)%assoc_elem = child1
   grid%face(face3)%marked_edge = grid%face(face3)%edge(1)
   grid%face(face3)%child = NO_CHILD
   grid%face(face3)%next = 0
   if (grid%face(face3)%degree > 2) then
      allocate(grid%face(face3)%solution(face_dof(grid%face(face3)%degree), &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%face(face3)%exact(face_dof(grid%face(face3)%degree), &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%face(face3)%exact = 0.0_my_real
      else
         nullify(grid%face(face3)%exact)
      endif
      do i=1,grid%system_size
         if (grid%face_type(face3,i) == DIRICHLET) then
            call face_exact(grid,face3,i,"d")
         else
            grid%face(face3)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call face_exact(grid,face3,i,"t")
      end do
   else
      nullify(grid%face(face3)%solution)
      nullify(grid%face(face3)%exact)
   endif
   nullify(grid%face(face3)%oldsoln)
   grid%face(face3)%next = face3

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%face(face3)%gid,face3,grid%face_hash)
!!!!!endif

! face 4

   grid%face(face4)%gid = 4*grid%face(parent_face(2))%gid + 1
   grid%face(face4)%vertex = (/ child_vert(2), &
                                common_vert(2), &
                                vert1 /)
   grid%face(face4)%edge = (/ child_edge(4), &
                              edge2, edge4 /)
   grid%face(face4)%bmark = grid%face(refinement_face)%bmark
   grid%face(face4)%tags = grid%face(refinement_face)%tags
   grid%face(face4)%degree = grid%face(refinement_face)%degree
   grid%face_type(face4,:) = grid%face_type(refinement_face,:)
   grid%face(face4)%assoc_elem = child2
   grid%face(face4)%marked_edge = grid%face(face4)%edge(1)
   grid%face(face4)%child = NO_CHILD
   grid%face(face4)%next = 0
   if (grid%face(face4)%degree > 2) then
      allocate(grid%face(face4)%solution(face_dof(grid%face(face4)%degree), &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         stop
      endif
      if (grid%have_true) then
         allocate(grid%face(face4)%exact(face_dof(grid%face(face4)%degree), &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in h_refine_element")
            return
         endif
         grid%face(face4)%exact = 0.0_my_real
      else
         nullify(grid%face(face4)%exact)
      endif
      do i=1,grid%system_size
         if (grid%face_type(face4,i) == DIRICHLET) then
            call face_exact(grid,face4,i,"d")
         else
            grid%face(face4)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call face_exact(grid,face4,i,"t")
      end do
   else
      nullify(grid%face(face4)%solution)
      nullify(grid%face(face4)%exact)
   endif
   nullify(grid%face(face4)%oldsoln)
   grid%face(face4)%next = face4

!!!!!if (.not. present(elem_child_lid)) then
      call hash_insert(grid%face(face4)%gid,face4,grid%face_hash)
!!!!!endif

endif ! children of parent face 2 didn't exist

! define the new face interior to the parent

face5 = get_next_free_face(grid,errcode)
if (errcode /= 0) return

if (new_face2) then
   grid%face(face5)%gid = 4*grid%face(parent_face(2))%gid + 2
elseif (new_face1) then
   grid%face(face5)%gid = 4*grid%face(parent_face(1))%gid + 2
else
   grid%face(face5)%gid = 4*grid%face(parent_face(2))%gid + 3
endif
grid%face(face5)%vertex = (/ common_vert(1), &
                             common_vert(2), &
                             vert1 /)
grid%face(face5)%edge = (/ common_edge, &
                           edge3, edge4 /)
grid%face(face5)%bmark = 0
grid%face(face5)%tags = grid%element(parent)%tags
grid%face(face5)%degree = grid%element(parent)%degree
grid%face_type(face5,:) = INTERIOR
grid%face(face5)%assoc_elem = child1
grid%face(face5)%child = NO_CHILD
grid%face(face5)%next = 0
select case (grid%element(parent)%type)
case ("Pf")
   grid%face(face5)%marked_edge = grid%face(face5)%edge(2)
case default
   grid%face(face5)%marked_edge = grid%face(face5)%edge(1)
end select
if (grid%face(face5)%degree > 2) then
   allocate(grid%face(face5)%solution(face_dof(grid%face(face5)%degree), &
            grid%system_size,max(1,grid%num_eval)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in h_refine_element")
      stop
   endif
   if (grid%have_true) then
      allocate(grid%face(face5)%exact(face_dof(grid%face(face5)%degree), &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         return
      endif
      grid%face(face5)%exact = 0.0_my_real
   else
      nullify(grid%face(face5)%exact)
   endif
   do i=1,grid%system_size
      if (grid%face_type(face5,i) == DIRICHLET) then
         call face_exact(grid,face5,i,"d")
      else
         grid%face(face5)%solution(:,i,:) = 0.0_my_real
      endif
      if (grid%have_true) call face_exact(grid,face5,i,"t")
   end do
else
   nullify(grid%face(face5)%solution)
   nullify(grid%face(face5)%exact)
endif
nullify(grid%face(face5)%oldsoln)
grid%face(face5)%next = face5

!!!!!if (.not. present(elem_child_lid)) then
   call hash_insert(grid%face(face5)%gid,face5,grid%face_hash)
!!!!!endif

! define the two child elements

! first child

grid%element(child1)%gid = MAX_CHILD*grid%element(parent)%gid
grid%element(child1)%face = (/child_face(1), &
                              face1, face3, face5/)
grid%element(child1)%edge = (/child_edge(1), &
                              child_edge(2), &
                              common_edge, &
                              edge1, edge3, edge4/)
grid%element(child1)%vertex = (/child_vert(1), &
                                common_vert(1), &
                                common_vert(2), &
                                vert1/)
do i=1,FACES_PER_ELEMENT
   if (grid%element(parent)%face(i) == child_face(1)) then
      grid%element(child1)%neighbor_hint(1) = parent_neighbors(i)
   elseif (grid%element(parent)%face(i) == parent_face(1)) then
      grid%element(child1)%neighbor_hint(2) = parent_neighbors(i)
   elseif (grid%element(parent)%face(i) == parent_face(2)) then
      grid%element(child1)%neighbor_hint(3) = parent_neighbors(i)
   elseif (grid%element(parent)%face(i) == child_face(2)) then
      grid%element(child1)%neighbor_hint(4) = &
                     MAX_CHILD*grid%element(parent)%gid+1
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("h_refine_element: didn't match parent face when setting neighbor_hint for child1")
      stop
   endif
end do
grid%element(child1)%level = grid%element(parent)%level+1
! set mate to BOUNDARY so that processing of the mate is always skipped
grid%element(child1)%mate = BOUNDARY
grid%element(child1)%iown = grid%element(parent)%iown
grid%element(child1)%hrefined_unowned = .false.
grid%element(child1)%prefined_unowned = .false.
! TEMP3D doesn't matter until I do REFTREE partitioning
grid%element(child1)%in  = grid%element(child1)%vertex(1)
grid%element(child1)%out = grid%element(child1)%vertex(2)
!!!!!if (.not. present(elem_child_lid)) then
   grid%element(child1)%next = grid%head_level_elem(grid%element(child1)%level)
   grid%head_level_elem(grid%element(child1)%level) = child1
   if (grid%element(child1)%next /= END_OF_LIST) then
      grid%element(grid%element(child1)%next)%previous = child1
   endif
   grid%element(child1)%previous = END_OF_LIST
!!!!!endif
grid%element(child1)%isleaf = .true.
grid%element(child1)%oldleaf = .false.
grid%element(child1)%tags = grid%element(parent)%tags
grid%element(child1)%degree = grid%element(parent)%degree
if (refcont%reftype == HP_ADAPTIVE .and. refcont%hp_strategy == HP_SMOOTH_PRED) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: have not determined eta for SMOOTH_PRED yet")
   stop
endif
grid%element(child1)%sp_eta_pred = 0
if (associated(grid%element(parent)%solution)) then
   allocate(grid%element(child1)%solution(size(grid%element(parent)%solution, &
            dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in h_refine_element")
      stop
   endif
   grid%element(child1)%solution = 0.0_my_real
   if (grid%have_true) then
      allocate(grid%element(child1)%exact(size(grid%element(parent)%exact, &
               dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         return
      endif
      grid%element(child1)%exact = 0.0_my_real
   else
      nullify(grid%element(child1)%exact)
   endif
else
   nullify(grid%element(child1)%solution)
   nullify(grid%element(child1)%exact)
endif
nullify(grid%element(child1)%oldsoln)
if (grid%element(child1)%degree >= 4) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + grid%system_size * &
!!!!!                              element_dof(grid%element(child1)%degree)
!!!!!   else
      grid%dof = grid%dof + grid%system_size * &
                            element_dof(grid%element(child1)%degree)
!!!!!   endif
endif
if (grid%element(child1)%iown .and. grid%element(child1)%degree >= 4) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof_own = delta_dof_own + grid%system_size * &
!!!!!                                      element_dof(grid%element(child1)%degree)
!!!!!   else
      grid%dof_own = grid%dof_own + grid%system_size * &
                                    element_dof(grid%element(child1)%degree)
!!!!!   endif
endif
!!!!!if (.not. present(elem_child_lid)) then
   call hash_insert(grid%element(child1)%gid,child1,grid%elem_hash)
!!!!!endif
select case (grid%element(parent)%type)
case ("Pu")
   grid%element(child1)%type = "Pf"
case ("A")
   grid%element(child1)%type = "Pu"
case ("M")
   grid%element(child1)%type = "Pu"
case ("O")
   grid%element(child1)%type = "Pu"
case ("Pf")
   grid%element(child1)%type = "A"
end select
grid%element(child1)%refinement_edge = &
            grid%face(grid%element(child1)%face(1))%marked_edge
! TEMP120716 for longest edge bisection
!grid%element(child1)%refinement_edge = longest_edge(grid,child1)

! second child

grid%element(child2)%gid = MAX_CHILD*grid%element(parent)%gid + 1
grid%element(child2)%face = (/child_face(2), &
                              face2, face4, face5/)
grid%element(child2)%edge = (/child_edge(3), &
                              child_edge(4), &
                              common_edge, &
                              edge2, edge3, edge4/)
grid%element(child2)%vertex = (/child_vert(2), &
                                common_vert(1), &
                                common_vert(2), &
                                vert1/)
do i=1,FACES_PER_ELEMENT
   if (grid%element(parent)%face(i) == child_face(2)) then
      grid%element(child2)%neighbor_hint(1) = parent_neighbors(i)
   elseif (grid%element(parent)%face(i) == parent_face(1)) then
      grid%element(child2)%neighbor_hint(2) = parent_neighbors(i)
   elseif (grid%element(parent)%face(i) == parent_face(2)) then
      grid%element(child2)%neighbor_hint(3) = parent_neighbors(i)
   elseif (grid%element(parent)%face(i) == child_face(1)) then
      grid%element(child2)%neighbor_hint(4) = &
                     MAX_CHILD*grid%element(parent)%gid
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("h_refine_element: didn't match parent face when setting neighbor_hint for child2")
      stop
   endif
end do
grid%element(child2)%level = grid%element(parent)%level+1
grid%element(child2)%mate = BOUNDARY
grid%element(child2)%iown = grid%element(parent)%iown
grid%element(child2)%hrefined_unowned = .false.
grid%element(child2)%prefined_unowned = .false.
! TEMP3D doesn't matter until I do REFTREE partitioning
grid%element(child2)%in  = grid%element(child2)%vertex(1)
grid%element(child2)%out = grid%element(child2)%vertex(2)
!!!!!if (.not. present(elem_child_lid)) then
   grid%element(child2)%next = grid%head_level_elem(grid%element(child2)%level)
   grid%head_level_elem(grid%element(child2)%level) = child2
   if (grid%element(child2)%next /= END_OF_LIST) then
      grid%element(grid%element(child2)%next)%previous = child2
   endif
   grid%element(child2)%previous = END_OF_LIST
!!!!!endif
grid%element(child2)%isleaf = .true.
grid%element(child2)%oldleaf = .false.
grid%element(child2)%tags = grid%element(parent)%tags
grid%element(child2)%degree = grid%element(parent)%degree
if (refcont%reftype == HP_ADAPTIVE .and. refcont%hp_strategy == HP_SMOOTH_PRED) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("h_refine_element: have not determined eta for SMOOTH_PRED yet")
   stop
endif
grid%element(child2)%sp_eta_pred = 0
if (associated(grid%element(parent)%solution)) then
   allocate(grid%element(child2)%solution(size(grid%element(parent)%solution, &
            dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in h_refine_element")
      stop
   endif
   grid%element(child2)%solution = 0.0_my_real
   if (grid%have_true) then
      allocate(grid%element(child2)%exact(size(grid%element(parent)%exact, &
               dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_element")
         return
      endif
      grid%element(child2)%exact = 0.0_my_real
   else
      nullify(grid%element(child2)%exact)
   endif
else
   nullify(grid%element(child2)%solution)
   nullify(grid%element(child2)%exact)
endif
nullify(grid%element(child2)%oldsoln)
if (grid%element(child2)%degree >= 4) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + grid%system_size * &
!!!!!                              element_dof(grid%element(child2)%degree)
!!!!!   else
      grid%dof = grid%dof + grid%system_size * &
                            element_dof(grid%element(child2)%degree)
!!!!!   endif
endif
if (grid%element(child2)%iown .and. grid%element(child2)%degree >= 4) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof_own = delta_dof_own + grid%system_size * &
!!!!!                                      element_dof(grid%element(child2)%degree)
!!!!!   else
      grid%dof_own = grid%dof_own + grid%system_size * &
                                    element_dof(grid%element(child2)%degree)
!!!!!   endif
endif
!!!!!if (.not. present(elem_child_lid)) then
   call hash_insert(grid%element(child2)%gid,child2,grid%elem_hash)
!!!!!endif
select case (grid%element(parent)%type)
case ("Pu")
   grid%element(child2)%type = "Pf"
case ("A")
   grid%element(child2)%type = "Pu"
case ("M")
   grid%element(child2)%type = "Pu"
case ("O")
   grid%element(child2)%type = "Pu"
case ("Pf")
   grid%element(child2)%type = "A"
end select
grid%element(child2)%refinement_edge = &
            grid%face(grid%element(child2)%face(1))%marked_edge
! TEMP120716 for longest edge bisection
!grid%element(child2)%refinement_edge = longest_edge(grid,child2)

! the parent is no longer a leaf

grid%element(parent)%isleaf = .false.

! TEMP3D order doesn't matter until I do REFTREE partitioning

grid%element(parent)%order = (/1,2/)

! if parent face 1 is periodic, then set the "next" pointers if the matching
! face has already been refined, and likewise for parent face 2

do_face1 = is_periodic_face(parent_face(1),grid)
if (do_face1) then
   do_face1 = grid%face(grid%face(parent_face(1))%next)%child(1) /= NO_CHILD
endif
do_face2 = is_periodic_face(parent_face(2),grid)
if (do_face2) then
   do_face2 = grid%face(grid%face(parent_face(2))%next)%child(1) /= NO_CHILD
endif

if (do_face1 .or. do_face2) then

! The matching faces are the children of the periodic mate of the face that was
! refined.  The child that matches face1 is the one that contains the matching
! edge for face1's first edge, etc.  Note that the first edges of face1 and
! face2 can't both be doubly periodic because they share a common point, so
! one of them should have the correct mate as the first "next".

   if (do_face1) then
      mate_face_1 = grid%face(grid%face(parent_face(1))%next)%child(1)
      mate_face_2 = grid%face(grid%face(parent_face(1))%next)%child(2)
      if (any(grid%face(mate_face_2)%edge == &
              grid%edge(grid%face(face1)%edge(1))%next) .or. &
          any(grid%face(mate_face_1)%edge == &
              grid%edge(grid%face(face2)%edge(1))%next)) then
         itemp = mate_face_1
         mate_face_1 = mate_face_2
         mate_face_2 = itemp
      elseif (all(grid%face(mate_face_1)%edge /= &
                  grid%edge(grid%face(face1)%edge(1))%next) .and. &
              all(grid%face(mate_face_2)%edge /= &
                  grid%edge(grid%face(face2)%edge(1))%next)) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("couldn't find periodic matching edge in matching face")
         stop
      endif
   endif
   if (do_face2) then
      mate_face_3 = grid%face(grid%face(parent_face(2))%next)%child(1)
      mate_face_4 = grid%face(grid%face(parent_face(2))%next)%child(2)
      if (any(grid%face(mate_face_4)%edge == &
              grid%edge(grid%face(face3)%edge(1))%next) .or. &
          any(grid%face(mate_face_3)%edge == &
              grid%edge(grid%face(face4)%edge(1))%next)) then
         itemp = mate_face_3
         mate_face_3 = mate_face_4
         mate_face_4 = itemp
      elseif (all(grid%face(mate_face_3)%edge /= &
                  grid%edge(grid%face(face3)%edge(1))%next) .and. &
              all(grid%face(mate_face_4)%edge /= &
                  grid%edge(grid%face(face4)%edge(1))%next)) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("couldn't find periodic matching edge in matching face")
         stop
      endif
   endif

! The matching vertex of the new vertex is the third vertex of any matching
! face.

   if (do_face1) then
      mate_vert_1 = grid%face(mate_face_1)%vertex(3)
   else
      mate_vert_1 = grid%face(mate_face_3)%vertex(3)
   endif

! The matching edge of edge1 and edge2 is the second edge of the matching face.

   if (do_face1) then
      mate_edge_1 = grid%face(mate_face_1)%edge(2)
      mate_edge_2 = grid%face(mate_face_2)%edge(2)
   else
      mate_edge_1 = grid%face(mate_face_3)%edge(2)
      mate_edge_2 = grid%face(mate_face_4)%edge(2)
   endif

! For the edge that cuts the face the matching edge is the 3rd edge
! of the matching face.

   if (do_face1) then
      mate_edge_3 = grid%face(mate_face_1)%edge(3)
   endif
   if (do_face2) then
      mate_edge_4 = grid%face(mate_face_3)%edge(3)
   endif

! The faces simply point at each other; it doesn't matter which is master.
! The edges that split faces cannot be doubly periodic, so they just point at 
! at each other.

   if (do_face1) then
      grid%face(face1)%next = mate_face_1
      grid%face(mate_face_1)%next = face1
      grid%face(face2)%next = mate_face_2
      grid%face(mate_face_2)%next = face2
      grid%edge(edge3)%next = mate_edge_3
      grid%edge(mate_edge_3)%next = edge3
   endif
   if (do_face2) then
      grid%face(face3)%next = mate_face_3
      grid%face(mate_face_3)%next = face3
      grid%face(face4)%next = mate_face_4
      grid%face(mate_face_4)%next = face4
      grid%edge(edge4)%next = mate_edge_4
      grid%edge(mate_edge_4)%next = edge4
   endif

! three possibilities for the edges:
! 1) mate_edge_1 is already in a cycle of periodic edges which does not
!    contain edge1 (indicating it is doubly periodic)
! 2) mate_edge_1 and edge1 are already in a cycle together (previously assigned)
! 3) mate_edge_1 is not in a cycle; if we start at mate_edge_1 and don't come
!    back to it within 5 nexts, it is not in a cycle
! mate_edge_2 must be in the same state as mate_edge_1

   mate_edge_1_in_cycle = .false.
   edge1_in_cycle = .false.
   itemp = mate_edge_1
   do i=1,5
      itemp = grid%edge(itemp)%next
      if (itemp == END_OF_LIST) exit
      if (itemp == edge1) edge1_in_cycle = .true.
      if (itemp == mate_edge_1) then
         mate_edge_1_in_cycle = .true.
         exit
      endif
   end do
   edge1_in_cycle = edge1_in_cycle .and. mate_edge_1_in_cycle

! it doesn't matter where the master is in the cycle.  If mate_edge_1
! is already in a cycle and edge1 is not, insert edge1 after mate_edge_1.  If
! it is not, make them point to each other.

   if (.not. edge1_in_cycle) then
      if (mate_edge_1_in_cycle) then
         grid%edge(edge1)%next = grid%edge(mate_edge_1)%next
         grid%edge(mate_edge_1)%next = edge1
         grid%edge(edge2)%next = grid%edge(mate_edge_2)%next
         grid%edge(mate_edge_2)%next = edge2
      else
         grid%edge(edge1)%next = mate_edge_1
         grid%edge(mate_edge_1)%next = edge1
         grid%edge(edge2)%next = mate_edge_2
         grid%edge(mate_edge_2)%next = edge2
      endif
   endif

! if both face1 and face2 are periodic, we have only done the edges for the
! first one; do it again with the matching edge of the second face, but now
! the insertion is a little different

   if (do_face1 .and. do_face2) then

      mate_edge_1 = grid%face(mate_face_3)%edge(2)
      mate_edge_2 = grid%face(mate_face_4)%edge(2)

      mate_edge_1_in_cycle = .false.
      edge1_in_cycle = .false.
      itemp = mate_edge_1
      do i=1,5
         itemp = grid%edge(itemp)%next
         if (itemp == END_OF_LIST) exit
         if (itemp == edge1) edge1_in_cycle = .true.
         if (itemp == mate_edge_1) then
            mate_edge_1_in_cycle = .true.
            exit
         endif
      end do
      edge1_in_cycle = edge1_in_cycle .and. mate_edge_1_in_cycle

      if (.not. edge1_in_cycle) then
         if (mate_edge_1_in_cycle) then
            itemp = grid%edge(edge1)%next
            grid%edge(edge1)%next = grid%edge(mate_edge_1)%next
            grid%edge(mate_edge_1)%next = itemp
            itemp = grid%edge(edge2)%next
            grid%edge(edge2)%next = grid%edge(mate_edge_2)%next
            grid%edge(mate_edge_2)%next = itemp
         else
            grid%edge(mate_edge_1)%next = grid%edge(mate_edge_1)%next
            grid%edge(edge1)%next = mate_edge_1
            grid%edge(mate_edge_2)%next = grid%edge(mate_edge_2)%next
            grid%edge(edge2)%next = mate_edge_2
         endif
      endif

   endif

! same for vertices, except check 9 for a cycle

   mate_vert_1_in_cycle = .false.
   vert1_in_cycle = .false.
   itemp = mate_vert_1
   do i=1,9
      itemp = grid%vertex(itemp)%next
      if (itemp == END_OF_LIST) exit
      if (itemp == vert1) vert1_in_cycle = .true.
      if (itemp == mate_vert_1) then
         mate_vert_1_in_cycle = .true.
         exit
      endif
   end do
   vert1_in_cycle = vert1_in_cycle .and. mate_vert_1_in_cycle

   if (.not. vert1_in_cycle) then
      if (mate_vert_1_in_cycle) then
         grid%vertex(vert1)%next = grid%vertex(mate_vert_1)%next
         grid%vertex(mate_vert_1)%next = vert1
      else
         grid%vertex(vert1)%next = mate_vert_1
         grid%vertex(mate_vert_1)%next = vert1
      endif
   endif

   if (do_face1 .and. do_face2) then

      mate_vert_1 = grid%face(mate_face_3)%vertex(3)

      mate_vert_1_in_cycle = .false.
      vert1_in_cycle = .false.
      itemp = mate_vert_1
      do i=1,9
         itemp = grid%vertex(itemp)%next
         if (itemp == END_OF_LIST) exit
         if (itemp == vert1) vert1_in_cycle = .true.
         if (itemp == mate_vert_1) then
            mate_vert_1_in_cycle = .true.
            exit
         endif
      end do
      vert1_in_cycle = vert1_in_cycle .and. mate_vert_1_in_cycle

      if (.not. vert1_in_cycle) then
         if (mate_vert_1_in_cycle) then
            itemp = grid%vertex(vert1)%next
            grid%vertex(vert1)%next = grid%vertex(mate_vert_1)%next
            grid%vertex(mate_vert_1)%next = itemp
         else
            grid%vertex(mate_vert_1)%next = grid%vertex(mate_vert_1)%next
            grid%vertex(vert1)%next = mate_vert_1
         endif
      endif

   endif

endif ! do_face1 or do_face2

! associated element for the vertices of the parent

! TEMP3D when we make this MPI parallel, I think we need to make the
!        associated element for periodic entities be on the periodic master
!        side, so we don't have the master and slave assigned to different
!        processors.  It looks like in 2D I assign the associated element of
!        a periodic slave vertex (and edge?) to be the associated element of
!        the master, but check that.

! make this OpenMP critical just to be sure another thread doesn't try to
! read it while this thread is writing it

!$omp critical (reset_assoc_elem_critical)

if (grid%vertex(child_vert(1))%assoc_elem == parent) then
   grid%vertex(child_vert(1))%assoc_elem = child1
endif
if (grid%vertex(child_vert(2))%assoc_elem == parent) then
   grid%vertex(child_vert(2))%assoc_elem = child2
endif
if (grid%vertex(common_vert(1))%assoc_elem == parent) then
   grid%vertex(common_vert(1))%assoc_elem = child1
endif
if (grid%vertex(common_vert(2))%assoc_elem == parent) then
   grid%vertex(common_vert(2))%assoc_elem = child1
endif

! associated elements of the edges of the parent

if (grid%edge(child_edge(1))%assoc_elem == parent) then
   grid%edge(child_edge(1))%assoc_elem = child1
endif
if (grid%edge(child_edge(2))%assoc_elem == parent) then
   grid%edge(child_edge(2))%assoc_elem = child1
endif
if (grid%edge(child_edge(3))%assoc_elem == parent) then
   grid%edge(child_edge(3))%assoc_elem = child2
endif
if (grid%edge(child_edge(4))%assoc_elem == parent) then
   grid%edge(child_edge(4))%assoc_elem = child2
endif
if (grid%edge(common_edge)%assoc_elem == parent) then
   grid%edge(common_edge)%assoc_elem = child1
endif

! associated elements of the faces of the parent

if (grid%face(child_face(1))%assoc_elem == parent) then
   grid%face(child_face(1))%assoc_elem = child1
endif
if (grid%face(child_face(2))%assoc_elem == parent) then
   grid%face(child_face(2))%assoc_elem = child2
endif

!$omp end critical (reset_assoc_elem_critical)

! children of the bisected edge and faces

if (new_ref_edge) then
   grid%edge(refinement_edge)%child = (/edge1,edge2/)
endif
if (new_face1) then
   grid%face(parent_face(1))%child = (/face1,face2/)
endif
if (new_face2) then
   grid%face(parent_face(2))%child = (/face3,face4/)
endif

! grid scalars

if (new_ref_edge) then

   grid%nvert = grid%nvert + 1
   if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
      grid%nvert_own = grid%nvert_own + 1
   endif
   grid%nedge = grid%nedge + 2
! TEMP also need to change nedge_own

!!!!!if (present(delta_dof)) then
!!!!!   delta_dof = delta_dof + grid%system_size
!!!!!   if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) &
!!!!!      delta_dof_own = delta_dof_own + grid%system_size
!!!!!else
      grid%dof = grid%dof + grid%system_size
      if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) &
         grid%dof_own = grid%dof_own + grid%system_size
!!!!!endif

   if (grid%edge(edge1)%degree >= 2) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + grid%system_size*(grid%edge(edge1)%degree - 1)
!!!!!      if (grid%element(grid%edge(edge1)%assoc_elem)%iown) &
!!!!!         delta_dof_own = delta_dof_own + &
!!!!!                         grid%system_size*(grid%edge(edge1)%degree - 1)
!!!!!   else
         grid%dof = grid%dof + grid%system_size*(grid%edge(edge1)%degree - 1)
         if (grid%element(grid%edge(edge1)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*(grid%edge(edge1)%degree - 1)
!!!!!   endif
   endif

   if (grid%edge(edge2)%degree >= 2) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + grid%system_size*(grid%edge(edge2)%degree - 1)
!!!!!      if (grid%element(grid%edge(edge2)%assoc_elem)%iown) &
!!!!!         delta_dof_own = delta_dof_own + &
!!!!!                         grid%system_size*(grid%edge(edge2)%degree - 1)
!!!!!   else
         grid%dof = grid%dof + grid%system_size*(grid%edge(edge2)%degree - 1)
         if (grid%element(grid%edge(edge2)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*(grid%edge(edge2)%degree - 1)
!!!!!   endif
   endif

endif

if (new_face1) then

   grid%nedge = grid%nedge + 1
   grid%nface = grid%nface + 2

   if (grid%edge(edge3)%degree >= 2) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + grid%system_size*(grid%edge(edge3)%degree - 1)
!!!!!      if (grid%element(grid%edge(edge3)%assoc_elem)%iown) &
!!!!!         delta_dof_own = delta_dof_own + &
!!!!!                         grid%system_size*(grid%edge(edge3)%degree - 1)
!!!!!   else
         grid%dof = grid%dof + grid%system_size*(grid%edge(edge3)%degree - 1)
         if (grid%element(grid%edge(edge3)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*(grid%edge(edge3)%degree - 1)
!!!!!   endif
   endif

   if (grid%face(face1)%degree >= 3) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + &
!!!!!                  grid%system_size*face_dof(grid%face(face1)%degree)
!!!!!      if (grid%element(grid%face(face1)%assoc_elem)%iown) &
!!!!!         delta_dof_own = delta_dof_own + &
!!!!!                         grid%system_size*face_dof(grid%face(face1)%degree)
!!!!!   else
         grid%dof = grid%dof + grid%system_size*face_dof(grid%face(face1)%degree)
         if (grid%element(grid%face(face1)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*face_dof(grid%face(face1)%degree)
!!!!!   endif
   endif

   if (grid%face(face2)%degree >= 3) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + &
!!!!!                  grid%system_size*face_dof(grid%face(face2)%degree)
!!!!!      if (grid%element(grid%face(face2)%assoc_elem)%iown) &
!!!!!         delta_dof_own = delta_dof_own + &
!!!!!                         grid%system_size*face_dof(grid%face(face2)%degree)
!!!!!   else
         grid%dof = grid%dof + grid%system_size*face_dof(grid%face(face2)%degree)
         if (grid%element(grid%face(face2)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*face_dof(grid%face(face2)%degree)
!!!!!   endif
   endif

endif

if (new_face2) then

   grid%nface = grid%nface + 2
   grid%nedge = grid%nedge + 1

   if (grid%edge(edge4)%degree >= 2) then
!!!!!      if (present(delta_dof)) then
!!!!!         delta_dof = delta_dof + grid%system_size*(grid%edge(edge4)%degree - 1)
!!!!!         if (grid%element(grid%edge(edge4)%assoc_elem)%iown) &
!!!!!            delta_dof_own = delta_dof_own + &
!!!!!                            grid%system_size*(grid%edge(edge4)%degree - 1)
!!!!!      else
         grid%dof = grid%dof + grid%system_size*(grid%edge(edge4)%degree - 1)
         if (grid%element(grid%edge(edge4)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*(grid%edge(edge4)%degree - 1)
!!!!!      endif
   endif

   if (grid%face(face3)%degree >= 3) then
!!!!!      if (present(delta_dof)) then
!!!!!         delta_dof = delta_dof + &
!!!!!                     grid%system_size*face_dof(grid%face(face3)%degree)
!!!!!         if (grid%element(grid%face(face3)%assoc_elem)%iown) &
!!!!!            delta_dof_own = delta_dof_own + &
!!!!!                            grid%system_size*face_dof(grid%face(face3)%degree)
!!!!!      else
         grid%dof = grid%dof + grid%system_size*face_dof(grid%face(face3)%degree)
         if (grid%element(grid%face(face3)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*face_dof(grid%face(face3)%degree)
!!!!!      endif
   endif

   if (grid%face(face4)%degree >= 3) then
!!!!!      if (present(delta_dof)) then
!!!!!         delta_dof = delta_dof + &
!!!!!                     grid%system_size*face_dof(grid%face(face4)%degree)
!!!!!         if (grid%element(grid%face(face4)%assoc_elem)%iown) &
!!!!!            delta_dof_own = delta_dof_own + &
!!!!!                            grid%system_size*face_dof(grid%face(face4)%degree)
!!!!!      else
         grid%dof = grid%dof + grid%system_size*face_dof(grid%face(face4)%degree)
         if (grid%element(grid%face(face4)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*face_dof(grid%face(face4)%degree)
!!!!!      endif
   endif

endif

if (grid%face(face5)%degree >= 3) then
!!!!!   if (present(delta_dof)) then
!!!!!      delta_dof = delta_dof + &
!!!!!                  grid%system_size*face_dof(grid%face(face5)%degree)
!!!!!      if (grid%element(grid%face(face5)%assoc_elem)%iown) &
!!!!!         delta_dof_own = delta_dof_own + &
!!!!!                         grid%system_size*face_dof(grid%face(face5)%degree)
!!!!!   else
      grid%dof = grid%dof + grid%system_size*face_dof(grid%face(face5)%degree)
      if (grid%element(grid%face(face5)%assoc_elem)%iown) &
         grid%dof_own = grid%dof_own + &
                        grid%system_size*face_dof(grid%face(face5)%degree)
!!!!!   endif
endif

grid%nelem = grid%nelem + 2
grid%nelem_leaf = grid%nelem_leaf + 1
if (grid%element(parent)%iown) then
   grid%nelem_leaf_own = grid%nelem_leaf_own + 1
endif
grid%nface = grid%nface + 1
grid%nlev = max(grid%nlev,grid%element(child1)%level)

return

contains

function get_ordered_neighbors(elem,jerr)
!        ---------------------

!----------------------------------------------------
! This routine returns the neighbors of element elem in the same order as
! the face they share.  If that neighbor has already been refined, it is
! the parent that is returned.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
integer, intent(out) :: jerr
type(hash_key) :: get_ordered_neighbors(NEIGHBORS_PER_ELEMENT)
!----------------------------------------------------
! Local variables:

integer :: neigh(NEIGHBORS_PER_ELEMENT), num_boundary, f, n, &
           neigh_parent(NEIGHBORS_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

jerr = 0

! get the neighbors

neigh = get_neighbors(elem,grid,missing_OK=.true.)

! get the parents of the neighbors

do n=1,NEIGHBORS_PER_ELEMENT
   if (neigh(n) == BOUNDARY) then
      neigh_parent(n) = BOUNDARY
   elseif (grid%element(neigh(n))%level == 1) then
      neigh_parent(n) = BOUNDARY ! not really, but it keeps it from looking
   elseif (grid%element(neigh(n))%gid/MAX_CHILD == &
           grid%element(elem)%gid/MAX_CHILD) then
      neigh_parent(n) = BOUNDARY ! avoid looking at my own parent
   else
      neigh_parent(n) = hash_decode_key(grid%element(neigh(n))%gid/MAX_CHILD, &
                                        grid%elem_hash)
      if (neigh_parent(n) == HASH_NOT_FOUND) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("get_ordered_neighbors: parent of neighbor not found")
         stop
      endif
   endif
end do

! for each face, find the neighbor that shares that face

num_boundary = 0
do f=1,FACES_PER_ELEMENT
   get_ordered_neighbors(f) = BOUNDARY
   do n=1,NEIGHBORS_PER_ELEMENT
      if (neigh(n) /= BOUNDARY) then
         if (is_periodic_face(grid%element(elem)%face(f),grid)) then
            if (any(grid%element(neigh(n))%face == &
                grid%face(grid%element(elem)%face(f))%next)) then
               get_ordered_neighbors(f) = grid%element(neigh(n))%gid
               exit
            endif
         else
            if (any(grid%element(neigh(n))%face == &
                grid%element(elem)%face(f))) then
               get_ordered_neighbors(f) = grid%element(neigh(n))%gid
               exit
            endif
         endif
      endif
      if (neigh_parent(n) /= BOUNDARY) then
         if (is_periodic_face(grid%element(elem)%face(f),grid)) then
            if (any(grid%element(neigh_parent(n))%face == &
                grid%face(grid%element(elem)%face(f))%next)) then
               get_ordered_neighbors(f) = grid%element(neigh_parent(n))%gid
               exit
            endif
         else
            if (any(grid%element(neigh_parent(n))%face == &
                grid%element(elem)%face(f))) then
               get_ordered_neighbors(f) = grid%element(neigh_parent(n))%gid
               exit
            endif
         endif
      endif
   end do
   if (get_ordered_neighbors(f) == BOUNDARY) then
      num_boundary = num_boundary + 1
   endif
end do

if (count(neigh == BOUNDARY) /= num_boundary) then
   jerr = 1
endif

end function get_ordered_neighbors

end subroutine h_refine_element

! TEMP120716 for longest edge bisection
function longest_edge(grid,elem)
type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
integer :: longest_edge
real(my_real) :: longest, this_length, small=1000*epsilon(0.0_my_real)
integer :: j
longest_edge = grid%element(elem)%edge(1)
longest = edge_length(grid,grid%element(elem)%edge(1))
! TEMP change the tie breaker to that used in init_grid3D
do j=2,EDGES_PER_ELEMENT
   this_length = edge_length(grid,grid%element(elem)%edge(j))
   if (this_length > longest + small*this_length .or. &
       (abs(this_length-longest) <= small*this_length .and. &
        grid%element(elem)%edge(j) < longest_edge)) then
      longest_edge = grid%element(elem)%edge(j)
      longest = this_length
   endif
end do
end function longest_edge
function edge_length(grid,edge)
type(grid_type), intent(in) :: grid
integer, intent(in) :: edge
real(my_real) :: edge_length
edge_length = sqrt( (grid%vertex(grid%edge(edge)%vertex(1))%coord%x - &
                     grid%vertex(grid%edge(edge)%vertex(2))%coord%x)**2 + &
                    (grid%vertex(grid%edge(edge)%vertex(1))%coord%y - &
                     grid%vertex(grid%edge(edge)%vertex(2))%coord%y)**2 + &
                    (grid%vertex(grid%edge(edge)%vertex(1))%coord%z - &
                     grid%vertex(grid%edge(edge)%vertex(2))%coord%z)**2 )
end function edge_length
! end TEMP120716

!          --------------
subroutine after_h_refine(grid,element_list,nelem,vert_lid,edge_lid,elem_lid)
!          --------------

!----------------------------------------------------
! This routine performs parts of h refinement that must be done by a single
! OpenMP thread after the OpenMP-parallel refinement.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: element_list(:), nelem, vert_lid(:,:), edge_lid(:,:), &
                       elem_lid(:,:)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! TEMP110606 don't stop
return

call fatal("after_h_refine not yet written for 3D")
stop
end subroutine after_h_refine

!                    --------------
recursive subroutine create_element(grid,elemgid,refine_control, &
                                    solver_control,err)
!                    --------------

!----------------------------------------------------
! This routine performs the refinements needed to create the element
! with global id elemgid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(hash_key), intent(in) :: elemgid
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
integer, intent(inout) :: err
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call fatal("create_element not yet written for 3D")
stop
end subroutine create_element

!          -----------------------
subroutine remove_from_errind_list(elem,elist)
!          -----------------------

!----------------------------------------------------
! This routine removes elem from the error indicator lists
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
type(errind_list) :: elist

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call fatal("remove_from_errind_list not yet written for 3D")
stop
end subroutine remove_from_errind_list

!          -------------
subroutine init_guess_ic(grid,elem,children)
!          -------------

!----------------------------------------------------
! This routine sets the solution in either element elem or the descendants
! of elem and its mate from the function in iconds.  The descendants are
! set if and only if children is present.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
integer, optional, intent(in) :: children(4)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call fatal("init_guess_ic not yet written for 3D")
stop
end subroutine init_guess_ic

!          ------------
subroutine init_guess_p(grid,elem,old_edge_deg)
!          ------------

!----------------------------------------------------
! This routine sets the initial guess for the solution components associated
! with bases that have been added by the p refinement of element elem
! by solving a local Dirichlet residual problem for the red p-hierarchical bases
! over element elem using a domain consisting of elem and its three neighbors.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout), target :: grid
integer, intent(in) :: elem, old_edge_deg(:)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call fatal("init_guess_p not yet written for 3D")
stop
end subroutine init_guess_p

!          --------------------
subroutine remove_hanging_nodes(grid,refinement_edge_list,nrefedge, &
                                refine_control,solver_control, &
                                desired_level,desired_degree,elem_list)
!          --------------------

!----------------------------------------------------
! This routine removes any hanging nodes in the grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(inout) :: refinement_edge_list(:), nrefedge
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
integer, optional, pointer :: desired_level(:), desired_degree(:), elem_list(:)
!----------------------------------------------------
! Local variables:

logical :: any_refined, put_edge_back
integer :: elem, i, j, k, errcode, mate, new_nrefedge, refinement_edge, &
           did_refinement_edge, ref_edge
integer :: nlev, lev
integer :: new_refinement_edge_list(grid%biggest_edge), child(MAX_CHILD)
integer, pointer :: elements(:)
!----------------------------------------------------
! Begin executable code

! For each edge that was a refinement edge, check all elements that have
! that edge to make sure that edge has been refined.  In the process, make
! a new list of refinement edges and repeat until there are no refinements.

do while (nrefedge /= 0)

   new_nrefedge = 0
   do i=1,nrefedge
      put_edge_back = .false.
      refinement_edge = refinement_edge_list(i)
      call get_edge_elements(grid,refinement_edge,elements,.true.)
      do j=1,size(elements)
         elem = elements(j)
         if (grid%element(elem)%isleaf) then
            ref_edge = refinement_edge
            do
               if (any(grid%element(elem)%edge == ref_edge)) then
                  call h_refine_element(grid,elem,errcode, &
                                        did_refinement_edge, &
                                        refine_control,solver_control, &
                                        desired_level=desired_level, &
                                        desired_degree=desired_degree, &
                                        elem_list=elem_list)
                  if (did_refinement_edge /= -1) then
                     new_nrefedge = new_nrefedge + 1
                     new_refinement_edge_list(new_nrefedge)=did_refinement_edge
                  endif

! Also have to check the children for refinement edge in case a different
! edge was refined.

                  if (.not. put_edge_back) then
                     child = get_child_lid(grid%element(elem)%gid, &
                                           ALL_CHILDREN,grid%elem_hash)
                     do k=1,MAX_CHILD
                        if (any(grid%element(child(k))%edge==ref_edge)) then
                           new_nrefedge = new_nrefedge + 1
                           new_refinement_edge_list(new_nrefedge) = ref_edge
                           put_edge_back = .true.
                        exit
                        endif
                     end do
                  endif
               endif
               ref_edge = grid%edge(ref_edge)%next
               if (ref_edge == refinement_edge) exit
            end do
         endif
      end do
      deallocate(elements)
   end do

   nrefedge = new_nrefedge
   refinement_edge_list(1:nrefedge) = new_refinement_edge_list(1:nrefedge)

end do

end subroutine remove_hanging_nodes

end module refine_elements
