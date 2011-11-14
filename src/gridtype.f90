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

module gridtype_mod

!----------------------------------------------------
! This module contains data structures for the grid.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use hash_mod
use message_passing
!----------------------------------------------------

implicit none
private
public VERTICES_PER_ELEMENT, EDGES_PER_ELEMENT, MAX_CHILD, ALL_CHILDREN, &
       END_OF_LIST, NOT_ON_LIST, NO_CHILD, BOUNDARY, point, element_t, &
       edge_t, vertex_t, grid_type, refine_options, errind_list, &
       triangle_data, binw

!----------------------------------------------------
! The following parameters are defined:

! values that depend on the kind of elements, polynomial degree, etc.
! RESTRICTION This version for bisected triangles.

integer, parameter :: VERTICES_PER_ELEMENT = 3, &
                      EDGES_PER_ELEMENT    = 3, &
                      MAX_CHILD            = 2

! argument for get_child when all children are wanted

integer, private :: i
integer, parameter :: ALL_CHILDREN(MAX_CHILD) = (/ (i,i=1,MAX_CHILD) /)

! flags

integer, parameter :: END_OF_LIST = -10, & ! end of linked list
                      NOT_ON_LIST = -11, & ! link value when not on any list
                      NO_CHILD    = -12    ! child does not exist

! neighbor or mate of an element when that side is on the boundary

integer, parameter :: BOUNDARY = -6

!----------------------------------------------------
! The following types are defined:

type point
   real(my_real) :: x,y ! RESTRICTION 2D, add z for 3D
end type point

type element_t
   type(hash_key) :: gid
   type(hash_key) :: mate ! RESTRICTION 2D, multiple mates in 3D
   real(my_real) :: weight
! subscripts are basis rank, system rank, eigenvalue
   real(my_real), pointer :: solution(:,:,:), exact(:,:,:), oldsoln(:,:,:)
   real(my_real) :: work
   real(my_real) :: sp_eta_pred
   integer :: vertex(VERTICES_PER_ELEMENT)
   integer :: edge(EDGES_PER_ELEMENT)
   integer :: degree
   integer :: level
   integer :: in, out
   integer :: order(MAX_CHILD)
   integer :: next, previous ! links for available memory list when not in use
                             ! and elements of one level when in use
   logical(small_logical) :: isleaf, oldleaf
   logical(small_logical) :: iown ! true if I own all leaves below this element
   logical(small_logical) :: hrefined_unowned, prefined_unowned
end type element_t

! edge owner is the same as it's vertex 2
! there is no linked list of edges in use.  next is used for linked list of
! free memory and bidirectional pointers between PERIODIC_SLAVEs and
! PERIODIC_MASTERs

type edge_t
   type(hash_key) :: gid
   integer :: vertex(2)
   integer :: bmark, degree, assoc_elem
! subscripts are basis rank, system rank, eigenvalue
   real(my_real), pointer :: solution(:,:,:), exact(:,:,:), oldsoln(:,:,:)
   integer :: next
end type edge_t

type vertex_t
   type(hash_key) :: gid
   type(point) :: coord
   real(my_real) :: bparam
   integer :: bmark
   integer :: assoc_elem
   integer :: next, previous
end type vertex_t

! if not solving an eigenproblem, num_eval is 0 and nsoln is system_size.

type grid_type
   type(element_t), pointer :: element(:)
   type(edge_t), pointer :: edge(:)
   type(vertex_t), pointer :: vertex(:)
   type(hash_table) :: elem_hash, edge_hash, vert_hash
   type(point) :: boundbox_min, boundbox_max
! subscripts are lid, system rank, eigenvalue
   real(my_real), pointer :: vertex_solution(:,:,:), vertex_exact(:,:,:), &
                             vertex_oldsoln(:,:,:)
! subscripts are element, eigenvalue
   real(my_real), pointer :: element_errind(:,:)
   real(my_real), pointer :: eigenvalue(:)
   real(my_real) :: eigen_linsys_max_l2_resid, eigen_linsys_ave_l2_resid
   real(my_real), pointer :: eigenprob_l2_resid(:), eigenprob_variance(:)
   real(my_real), pointer :: errest_energy(:), errest_Linf(:), errest_L2(:), &
                             errest_eigenvalue(:)
   real(my_real) :: max_blen
   real(my_real), pointer :: bp_start(:), bp_finish(:)
   integer, pointer :: edge_type(:,:), vertex_type(:,:)
   integer, pointer :: initial_neighbor(:,:) ! (EDGES_PER_ELEMENT,nelem_init)
   integer, pointer :: head_level_elem(:), head_level_vert(:)
   integer :: next_free_elem, next_free_edge, next_free_vert
   integer :: partition
   integer :: system_size, num_eval, nsoln
   integer :: nelem, nelem_leaf, nelem_leaf_own, nedge, nedge_own, &
              nvert, nvert_own, nlev, dof, dof_own
   integer :: arpack_iter, arpack_nconv, arpack_numop, arpack_numopb, &
              arpack_numreo, arpack_info
   integer :: errtype ! really belongs in io_options but that doesn't get
                      ! passed where it is needed
   logical :: errind_up2date, oldsoln_exists, have_true
   character(len=FN_LEN) :: triangle_files
end type grid_type

! data structure for data from triangle files, read and derived

type triangle_data
   type(point), pointer :: vert_coord(:)
   real(my_real), pointer :: vert_bparam(:)
   integer, pointer :: tri_edge(:,:), tri_vert(:,:), tri_neigh(:,:), &
                       edge_tri(:,:), edge_vert(:,:), edge_bmark(:), &
                       vert_tri(:,:), vert_edge(:,:), vert_bmark(:), &
                       vert_master(:), vert_master2(:)
   integer :: ntri, nedge, nvert
end type triangle_data

type refine_options
   real(my_real) :: inc_factor, reftol, term_energy_err, term_Linf_err, &
                    term_L2_err, t3s_gamma, t3s_eta, t3s_h_target, &
                    t3s_p_target, tp_gamma, sp_gamma_h, sp_gamma_p, &
                    refsoln_pbias
   integer :: error_estimator, reftype, refterm, edge_rule, hp_strategy, &
              max_vert, max_elem, max_dof, max_lev, max_deg, t3s_nunif, &
              t3s_maxref, t3s_maxdeginc, t3s_reftype, nlp_max_h_dec, &
              nlp_max_h_inc, nlp_max_p_dec, nlp_max_p_inc
   logical :: derefine, stop_on_maxlev, stop_on_maxdeg
end type refine_options

! error indicator lists, to determine which elements get refined next during
! adaptive refinement.  List k contains leaf elements in this partition for
! which m/binw**(k-1) > e > m/binw**k where m is max_errind, binw is the
! bin width (2 for cutoffs at 1, 1/2, 1/4, 1/8, ...)  and e is the error
! indicator for the element, with the last set going down to 0.0
! Currently using 16 lists of width 4th root of 2.

type errind_list
   integer :: current_list
   integer :: head_errind(16), tail_errind(16)
   integer, pointer :: next_errind(:), prev_errind(:)
   real(my_real) :: max_errind
end type

!real(my_real), parameter :: binw = sqrt(sqrt(2.0_my_real))
real(my_real), parameter :: binw = 1.1892071150027_my_real

end module gridtype_mod
