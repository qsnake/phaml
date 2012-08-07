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

! This file contains two modules, the grid type module that is very much
! like the 2D grid type module, and a boundary type module that contains
! additional types for handling boundary surfaces

module boundtype_mod

!----------------------------------------------------
! This module contains data structures for the boundary.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
!----------------------------------------------------

implicit none
private
public bvert_type, line_type, lineloop_type, surf_type, &
       SURF_PLANAR, SURF_RULED_3, SURF_RULED_4, &
       LINE_LINE, LINE_CIRCLE, LINE_ELLIPSE

!----------------------------------------------------
! The following parameters are defined:

integer, parameter :: SURF_PLANAR  = 1, &
                      SURF_RULED_3 = 2, &
                      SURF_RULED_4 = 3

integer, parameter :: LINE_LINE    = 1, &
                      LINE_CIRCLE  = 2, &
                      LINE_ELLIPSE = 3

!----------------------------------------------------
! The following types are defined:

type bvert_type
   integer :: surface      ! ID of associated surface
   integer :: lines(2)     ! ID of lines containing point, or 0
   real(my_real) :: p1, p2 ! PLANAR: fraction of way from lines(1,2)%point(1)
                           !         to the other end point
                           ! RULED_3: (b1,b2) note b3=1-b1-b2
                           ! RULED_4: (xi,eta)
   type(bvert_type), pointer :: next
end type bvert_type

type line_type
   integer :: type
   integer, pointer :: point(:)
end type line_type

type lineloop_type
   integer, pointer :: line(:)
end type lineloop_type

type surf_type
   integer :: type
   integer, pointer :: lineloop(:)
end type surf_type

end module boundtype_mod

module gridtype_mod

!----------------------------------------------------
! This module contains data structures for the grid.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use hash_mod
use message_passing
use boundtype_mod
!----------------------------------------------------

implicit none
private
public VERTICES_PER_ELEMENT, EDGES_PER_ELEMENT, FACES_PER_ELEMENT, &
       NEIGHBORS_PER_ELEMENT, VERTICES_PER_FACE, EDGES_PER_FACE, &
       MAX_CHILD, ALL_CHILDREN, END_OF_LIST, NOT_ON_LIST, NO_CHILD, BOUNDARY, &
       point, element_t, face_t, edge_t, vertex_t, grid_type, refine_options, &
       eigen_monitor,errind_list, binw, element_kind, MAX_TAG, point_type, &
       operator(+), operator(-), operator(*), operator(/), dot_point, &
       operator(.eq.)

!----------------------------------------------------
! The following interface operators are defined:

interface operator(+)
   module procedure point_plus_point
end interface

interface operator(-)
   module procedure point_minus_point
end interface

interface operator(*)
   module procedure scalar_times_point
end interface

interface operator(/)
   module procedure point_div_scalar
end interface

interface operator(.eq.)
   module procedure point_equals_point
end interface

!----------------------------------------------------
! The following parameters are defined:

! values that depend on the kind of elements, polynomial degree, etc.

integer, parameter :: VERTICES_PER_ELEMENT  = 4, &
                      EDGES_PER_ELEMENT     = 6, &
                      FACES_PER_ELEMENT     = 4, &
                      NEIGHBORS_PER_ELEMENT = 4, &
                      VERTICES_PER_FACE     = 3, &
                      EDGES_PER_FACE        = 3, &
                      MAX_CHILD             = 2

! argument for get_child when all children are wanted

integer, private :: i
integer, parameter :: ALL_CHILDREN(MAX_CHILD) = (/ (i,i=1,MAX_CHILD) /)

! flags

integer, parameter :: END_OF_LIST = -10, & ! end of linked list
                      NOT_ON_LIST = -11, & ! link value when not on any list
                      NO_CHILD    = -12    ! child does not exist

! neighbor or mate of an element when that side is on the boundary

integer, parameter :: BOUNDARY = -6

! elements are triangles

integer, parameter :: element_kind = TETRAHEDRAL_ELEMENT

! maximum number of tags for grid entities

integer, parameter :: MAX_TAG = 2

!----------------------------------------------------
! The following types are defined:

type eigen_monitor
   real(my_real), pointer :: eigenvalue(:)
   real(my_real), pointer :: eigensolver_l2_resid(:), eigensolver_errbound(:)
   integer :: ncv, maxit, niter, nconv
end type eigen_monitor

type point
   real(my_real) :: x,y,z
end type point

! This really belongs in boundtype_mod, but then I couldn't use type(point)
type point_type
   type(point) :: coord
end type point_type

type element_t
   type(hash_key) :: gid
! mate is not used in the 3D h-adaptive refinement; instead it is set to
! BOUNDARY so that processing of the mate is always skipped
   type(hash_key) :: mate
! neighbor_hint is the neighbor, an ancestor of the neighbor or a descendant
! of the neighbor that shares the corresponding face
   type(hash_key) :: neighbor_hint(NEIGHBORS_PER_ELEMENT)
   real(my_real) :: weight
! subscripts are basis rank, system rank, eigenvalue
   real(my_real), pointer :: solution(:,:,:), exact(:,:,:), oldsoln(:,:,:)
   real(my_real) :: work
   real(my_real) :: sp_eta_pred
   integer :: vertex(VERTICES_PER_ELEMENT)
   integer :: edge(EDGES_PER_ELEMENT)
   integer :: face(FACES_PER_ELEMENT)
   integer :: degree
   integer :: level
   integer :: in, out
   integer :: order(MAX_CHILD)
   integer :: tags(MAX_TAG)
   integer :: next, previous ! links for available memory list when not in use
                             ! and elements of one level when in use
   integer :: refinement_edge
   logical(small_logical) :: isleaf, oldleaf
   logical(small_logical) :: iown ! true if I own all leaves below this element
   logical(small_logical) :: hrefined_unowned, prefined_unowned
   character(len=2) :: type
end type element_t

! there is no linked list of vertices, edges or faces in use.  next is used for
! linked list of free memory and bidirectional pointers between PERIODIC_SLAVEs
! and PERIODIC_MASTERs

type face_t
   type(hash_key) :: gid
   integer :: edge(EDGES_PER_FACE)
   integer :: vertex(VERTICES_PER_FACE)
   integer :: child(2)
   integer :: bmark, degree, assoc_elem, marked_edge, tags(MAX_TAG)
! subscripts are basis rank, system rank, eigenvalue
   real(my_real), pointer :: solution(:,:,:), exact(:,:,:), oldsoln(:,:,:)
   integer :: next
end type face_t

type edge_t
   type(hash_key) :: gid
   integer :: vertex(2)
   integer :: child(2)
   integer :: bmark, degree, assoc_elem, tags(MAX_TAG)
! subscripts are basis rank, system rank, eigenvalue
   real(my_real), pointer :: solution(:,:,:), exact(:,:,:), oldsoln(:,:,:)
   integer :: next
end type edge_t

type vertex_t
   type(hash_key) :: gid
   type(point) :: coord
   real(my_real) :: bparam
   integer :: bmark, assoc_elem, tags(MAX_TAG)
   integer :: next
   type(bvert_type), pointer :: boundary_vertex
end type vertex_t

! if not solving an eigenproblem, num_eval is 0 and nsoln is system_size.

type grid_type
   type(element_t), pointer :: element(:) => NULL()
   type(face_t), pointer :: face(:) => NULL()
   type(edge_t), pointer :: edge(:) => NULL()
   type(vertex_t), pointer :: vertex(:) => NULL()
   type(hash_table) :: elem_hash, face_hash, edge_hash, vert_hash
   type(point) :: boundbox_min, boundbox_max
! subscripts are lid, system rank, eigenvalue
   real(my_real), pointer :: vertex_solution(:,:,:), vertex_exact(:,:,:), &
                             vertex_oldsoln(:,:,:)
! subscripts are element, eigenvalue
   real(my_real), pointer :: element_errind(:,:)
   real(my_real), pointer :: eigenprob_l2_resid(:), eigenprob_variance(:)
   real(my_real), pointer :: errest_energy(:), errest_Linf(:), errest_L2(:), &
                             errest_eigenvalue(:)
   real(my_real) :: max_blen
   real(my_real), pointer :: bp_start(:), bp_finish(:)
   integer, pointer :: face_type(:,:), edge_type(:,:), vertex_type(:,:)
   integer, pointer :: initial_neighbor(:,:) ! (NEIGHBORS_PER_ELEMENT,nelem_init)
   integer, pointer :: head_level_elem(:), head_level_vert(:)
   integer :: next_free_elem, next_free_face, next_free_edge, next_free_vert
   integer :: partition
   integer :: system_size, num_eval, nsoln
   integer :: nelem, nelem_leaf, nelem_leaf_own, nface, nface_own, nedge, &
              nedge_own, nvert, nvert_own, nlev, dof, dof_own
   integer :: biggest_vert, biggest_edge, biggest_face, biggest_elem
   integer :: errtype ! really belongs in io_options but that doesn't get
                      ! passed where it is needed
   logical :: errind_up2date, oldsoln_exists, have_true, any_periodic
   character(len=FN_LEN) :: triangle_files
   type(eigen_monitor) :: eigen_results
   type(point_type), pointer :: point(:)
   type(line_type), pointer :: line(:)
   type(lineloop_type), pointer :: lineloop(:)
   type(surf_type), pointer :: surface(:)
   integer :: npoint, nline, nlineloop, nsurface
end type grid_type

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

contains

! some operators on points

function point_plus_point(p1,p2)
type(point), intent(in) :: p1, p2
type(point) :: point_plus_point
point_plus_point%x = p1%x + p2%x
point_plus_point%y = p1%y + p2%y
point_plus_point%z = p1%z + p2%z
end function point_plus_point

function point_minus_point(p1,p2)
type(point), intent(in) :: p1, p2
type(point) :: point_minus_point
point_minus_point%x = p1%x - p2%x
point_minus_point%y = p1%y - p2%y
point_minus_point%z = p1%z - p2%z
end function point_minus_point

function scalar_times_point(a,p1)
real(my_real), intent(in) :: a
type(point), intent(in) :: p1
type(point) :: scalar_times_point
scalar_times_point%x = a*p1%x
scalar_times_point%y = a*p1%y
scalar_times_point%z = a*p1%z
end function scalar_times_point

function point_div_scalar(p1,a)
type(point), intent(in) :: p1
real(my_real), intent(in) :: a
type(point) :: point_div_scalar
point_div_scalar%x = p1%x/a
point_div_scalar%y = p1%y/a
point_div_scalar%z = p1%z/a
end function point_div_scalar

function point_equals_point(p1,p2)
type(point), intent(in) :: p1, p2
logical :: point_equals_point
point_equals_point = p1%x == p2%x .and. p1%y == p2%y .and. p1%z == p2%z
end function point_equals_point

function dot_point(p1,p2)
! This routine computes the dot product of two points (meaning the dot
! product of the vectors from the origin to the points)
type(point), intent(in) :: p1, p2
real(my_real) :: dot_point
dot_point = p1%x*p2%x + p1%y*p2%y + p1%z*p2%z
end function dot_point

end module gridtype_mod
