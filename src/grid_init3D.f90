!---------------------------------------------------------------------t
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

module grid_init_mod

!----------------------------------------------------
! This module contains routines to create, initialize and destroy the grid.
!
! communication tags in this module are of the form 31xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use boundtype_mod
use boundary_util
use grid_util
use message_passing
use hash_mod
use sort_mod
use stack_mod
use sysdep
!----------------------------------------------------

implicit none
private
public init_grid

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

integer :: SUFFIX_LENGTH=6
integer, parameter :: MAX_CONTAINS = 200        ! # faces containing a vertex
integer, parameter :: MAX_INPUT_LINE_LEN = 1000 ! line length in .geo
integer, parameter :: MAX_EXPRESSION = 1000     ! expression length in .geo
integer, parameter :: MAX_EXP_LIST = 100        ! expressions in list in .geo
integer, parameter :: MAX_VAR_LEN = 100         ! variable name length in .geo
integer, parameter :: MAX_NUM_LEN = 100         ! literal number length in .geo

! States for the expression parser

integer, parameter :: START_EXPRESSION          = 1, &
                      NUMBER                    = 2, &
                      SIGNIFICAND               = 3, &
                      EXPONENT                  = 4, &
                      EXPONENT_RIGHT_AFTER_SIGN = 5, &
                      EXPONENT_AFTER_SIGN       = 6, &
                      VARIABLE                  = 7, &
                      END_PAREN_EXPRESSION      = 8, &
                      FINISH                    = 9

! Single character tokens for the expression parser

integer, parameter :: DIGIT                          =  1, &
                      LETTER_NOT_EXPONENT            =  2, &
                      EXPONENT_LETTER                =  3, &
                      OPERATOR_BINARY_NOT_PLUS_MINUS =  4, &
                      PLUS                           =  5, &
                      MINUS                          =  6, &
                      DOT                            =  7, &
                      UNDERSCORE                     =  8, &
                      OPEN_PAREN                     =  9, &
                      CLOSE_PAREN                    = 10, &
                      END_OF_EXPRESSION              = 11

! Operators for the expression parser

integer, parameter :: EXPONENTIATION = 1, &
                      UNARY_MINUS    = 2, &
                      TIMES          = 3, &
                      DIVIDE         = 4, &
                      MODULUS        = 5, &
                      ADD            = 6, &
                      SUBTRACT       = 7, &
                      BEGIN_PAREN    = 8, &
                      EMPTY_STACK    = 9

! Precidence of the operators.  precidence(a,b) means
!   +1 a has higher precidence than b, e.g. precidence(*,+) = +1
!    0 a and b have equal precidence
!   -1 a has lower precidence than b

integer, parameter :: precidence(9,9) = reshape( (/ &
                                        0,  1,  1,  1,  1,  1,  1,  1,  1, &
                                       -1,  0,  1,  1,  1,  1,  1,  1,  1, &
                                       -1, -1,  0,  0,  0,  1,  1,  1,  1, &
                                       -1, -1,  0,  0,  0,  1,  1,  1,  1, &
                                       -1, -1,  0,  0,  0,  1,  1,  1,  1, &
                                       -1, -1, -1, -1, -1,  0,  0,  1,  1, &
                                       -1, -1, -1, -1, -1,  0,  0,  1,  1, &
                                       -1, -1, -1, -1, -1, -1, -1,  0,  1, &
                                       -1, -1, -1, -1, -1, -1, -1, -1,  0 &
                                       /), (/9,9/))

! for comparison of real numbers that should be equal

real(my_real), parameter :: small = 1000*epsilon(0.0_my_real), &
                            really_small = 1000*tiny(0.0_my_real)
!----------------------------------------------------
! The following types are defined:

type msh_point
   integer :: tags(MAX_TAG)
   integer :: vertex
end type msh_point

type msh_line
   integer :: tags(MAX_TAG)
   integer :: vertex(2)
end type msh_line

type msh_tri
   integer :: tags(MAX_TAG)
   integer :: vertex(3)
end type msh_tri

type input_line
   character(len=MAX_INPUT_LINE_LEN) :: line
   integer :: ptr, len
   type(stack_integer) :: stack
end type input_line

type var_name_tree
   character(len=MAX_VAR_LEN) :: string
   real(my_real) :: value
   type(var_name_tree), pointer :: left
   type(var_name_tree), pointer :: right
end type var_name_tree

type point_uv
   real(my_real) :: u, v
end type point_uv

!----------------------------------------------------
! The following variables are defined:

! TEMP120327
real(my_real) :: holeHeight, holeBottomRadius, holeTopRadius, cylinderRadius
integer :: CYLINDERSIDE = -1, HOLESIDE
! end TEMP120327
!----------------------------------------------------

contains

!          ---------
subroutine init_grid(grid,procs,degree,partition,set_iown)
!          ---------

!----------------------------------------------------
! This routine initializes the grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
integer, intent(in) :: degree,partition
logical, intent(in) :: set_iown

!----------------------------------------------------
! Local variables:

character(len=SUFFIX_LENGTH) :: suffix
character(len=len(grid%triangle_files)+6) :: geo_filename, msh_filename
!----------------------------------------------------
! Begin executable code

suffix = get_suffix(grid%triangle_files)
select case (suffix)

case ("msh")
   grid%nsurface = 0
   call init_grid_msh(grid,grid%triangle_files,procs,degree,partition,set_iown)

case ("geo", "geomsh")
   call make_geo_msh_filenames(grid%triangle_files,geo_filename,msh_filename)
   if (suffix == "geo") msh_filename = "phaml_"//msh_filename
   call process_geo(grid,procs,degree,partition,set_iown,geo_filename)
   if (suffix == "geo" .and. &
       (my_proc(procs) == MASTER .or. PARALLEL==SEQUENTIAL)) then
      call my_system("gmsh -3 -v 0 -optimize -o "//trim(msh_filename)//" "//&
                     trim(geo_filename))
      call my_system("gmsh --version")
   endif
   call phaml_barrier_master(procs)
   call init_grid_msh(grid,msh_filename,procs,degree,partition,set_iown)

case default
   call fatal("unrecognized file type for initial grid")
   stop

end select

end subroutine init_grid

!        ----------
function get_suffix(filename)
!        ----------

!----------------------------------------------------
! This routine returns the suffix of filename
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: filename
character(len=SUFFIX_LENGTH) :: get_suffix
!----------------------------------------------------
! Local variables:

integer :: dot
!----------------------------------------------------
! Begin executable code

dot = index(filename,".",back=.true.)
if (dot == 0) then
   get_suffix = ""
else
   get_suffix = filename(dot+1:)
endif

end function get_suffix

!          ----------------------
subroutine make_geo_msh_filenames(base_filename,geo_filename,msh_filename)
!          ----------------------

!----------------------------------------------------
! This routine returns the root base_filename with .geo and .msh suffixes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: base_filename
character(len=*), intent(out) :: geo_filename, msh_filename
!----------------------------------------------------
! Local variables:

character(len=SUFFIX_LENGTH) :: suffix
integer :: dot
!----------------------------------------------------
! Begin executable code

suffix = get_suffix(base_filename)

dot = index(base_filename,".",back=.true.)

geo_filename = base_filename
geo_filename(dot+1:dot+3) = "geo"
geo_filename(dot+4:) = " "

msh_filename = base_filename
msh_filename(dot+1:dot+3) = "msh"
msh_filename(dot+4:) = " "

end subroutine make_geo_msh_filenames

!          -------------
subroutine init_grid_msh(grid,msh_filename,procs,degree,partition,set_iown)
!          -------------

!----------------------------------------------------
! This routine creates the initial grid from a Gmsh .msh file.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
character(len=*), intent(in) :: msh_filename
type(proc_info), intent(in) :: procs
integer, intent(in) :: degree,partition
logical, intent(in) :: set_iown
!----------------------------------------------------
! Local variables:

integer :: i, j, k, m, astat, vert, edge, face, elem, loc_vert, loc_edge, &
           loc_face, existing_edge, existing_face, neigh_elem, first_gid, &
           elem_nedge, errcode, nmsh_point, nmsh_line, nmsh_tri, ivert
integer :: temp_edge(2), temp_face(VERTICES_PER_FACE)
integer, allocatable :: contains_vert(:,:), contains_face(:), ncontain(:)
real(my_real) :: x, y, z, longest, this_length, xvert(2), yvert(2), zvert(2), &
                 xvert_longest(2), yvert_longest(2), zvert_longest(2)
real(my_real) :: bccoef(grid%system_size,grid%system_size), &
                 bcrhs(grid%system_size)
type(msh_point), pointer :: msh_points(:)
type(msh_line), pointer :: msh_lines(:)
type(msh_tri), pointer :: msh_tris(:)
logical, allocatable :: visited_vert(:)
logical :: change_to_this
!----------------------------------------------------
! Begin executable code

! read the mesh file.  this sets vertex%coord, nvert, element%vertex, nelem

call read_msh(grid,msh_filename,msh_points,nmsh_point,msh_lines,nmsh_line, &
              msh_tris,nmsh_tri)

! all elements are leaves and owned by everyone

grid%nelem_leaf = grid%nelem
grid%biggest_elem = grid%nelem
if (set_iown) then
   grid%nelem_leaf_own = grid%nelem
   grid%nvert_own = grid%nvert
else
   grid%nelem_leaf_own = 0
   grid%nvert_own = 0
endif
grid%nlev = 1
grid%next_free_vert = grid%nvert+1

! set the bounding box

x = minval(grid%vertex(1:grid%nvert)%coord%x)
y = minval(grid%vertex(1:grid%nvert)%coord%y)
z = minval(grid%vertex(1:grid%nvert)%coord%z)
grid%boundbox_min = point(x,y,z)
x = maxval(grid%vertex(1:grid%nvert)%coord%x)
y = maxval(grid%vertex(1:grid%nvert)%coord%y)
z = maxval(grid%vertex(1:grid%nvert)%coord%z)
grid%boundbox_max = point(x,y,z)

! check for true solution known

if (my_trues((grid%boundbox_max%x+grid%boundbox_min%x)/2, &
             (grid%boundbox_max%y+grid%boundbox_min%y)/2, &
             (grid%boundbox_max%z+grid%boundbox_min%z)/2, &
             1, 1) == huge(0.0_my_real)) then
   grid%have_true = .false.
else
   grid%have_true = .true.
   allocate(grid%vertex_exact(size(grid%vertex),grid%system_size, &
            max(1,grid%num_eval)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in init_grid")
      stop
   endif
endif

! Set the element linked list

grid%head_level_elem(1) = 1
grid%element(1)%next = 2
grid%element(1)%previous = END_OF_LIST
do i=2,grid%nelem-1
   grid%element(i)%next = i+1
   grid%element(i)%previous = i-1
end do
grid%element(grid%nelem)%next = END_OF_LIST
grid%element(grid%nelem)%previous = grid%nelem-1
grid%next_free_elem = grid%nelem+1

! finish defining the vertices

grid%biggest_vert = grid%nvert
allocate(visited_vert(grid%biggest_vert),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in init_grid_msh")
   stop
endif

visited_vert = .false.
do i=1,grid%nelem
   do ivert=1,VERTICES_PER_ELEMENT
      vert = grid%element(i)%vertex(ivert)
      if (visited_vert(vert)) cycle
      visited_vert(vert) = .true.
      grid%vertex(vert)%gid = vert
      call hash_insert(grid%vertex(vert)%gid,vert,grid%vert_hash)
      if (grid%have_true) then
         do j=1,grid%system_size
            do k=1,max(1,grid%num_eval)
               grid%vertex_exact(vert,j,k) = &
                  my_trues(grid%vertex(vert)%coord%x,grid%vertex(vert)%coord%y,&
                           grid%vertex(vert)%coord%z,j,k)
            end do
         end do
      endif
      nullify(grid%vertex(vert)%boundary_vertex)
      grid%vertex(vert)%assoc_elem = i
      grid%vertex(vert)%next = vert
   end do
end do

! remove any vertices that are not used

do i=1,grid%biggest_vert
   if (.not. visited_vert(i)) then
      call put_next_free_vert(grid,i)
      grid%nvert = grid%nvert - 1
      if (set_iown) grid%nvert_own = grid%nvert_own - 1
   endif
end do

! Create the faces by going through the elements and creating any of the
! four faces that don't already exist.  For each vertex, keep a list of the
! containing faces to check whether a face has already been created. Assign the
! faces to elements with first face opposite first vertex, etc.  This is
! a good time to set initial neighbors since we know which two elements
! share a face, and set the face associated element as the first element that
! created the face.

if (size(grid%face) < 8*size(grid%vertex)) then
   call more_faces(grid,errcode,size_wanted=8*size(grid%vertex))
endif

allocate(contains_vert(MAX_CONTAINS,grid%biggest_vert), &
         ncontain(grid%biggest_vert), contains_face(size(grid%face)), &
         grid%initial_neighbor(NEIGHBORS_PER_ELEMENT,grid%nelem),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in init_grid_msh")
   stop
endif

ncontain = 0
grid%initial_neighbor = BOUNDARY

face = 0
do elem=1,grid%nelem
   do loc_face=1,FACES_PER_ELEMENT
      i = 0
! vertices of this face
      do loc_vert=1,VERTICES_PER_ELEMENT
         if (loc_vert == loc_face) cycle
         i = i+1
         temp_face(i) = grid%element(elem)%vertex(loc_vert)
      end do
! does the face already exist?
      existing_face = -1
      outer1:do i=1,ncontain(temp_face(1))
         do j=1,VERTICES_PER_FACE
            if (all(grid%face(contains_vert(i,temp_face(1)))%vertex /= &
                temp_face(j))) then
               cycle outer1
            endif
         end do
         existing_face = contains_vert(i,temp_face(1))
         exit outer1
      end do outer1
! exists
      if (existing_face /= -1) then
         grid%element(elem)%face(loc_face) = existing_face
         neigh_elem = contains_face(existing_face)
         grid%initial_neighbor(loc_face,elem) = neigh_elem
         do i=1,FACES_PER_ELEMENT
            if (grid%element(neigh_elem)%face(i) == existing_face) exit
         end do
         if (i > FACES_PER_ELEMENT) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("didn't find existing face on neighbor in init_grid_msh")
            stop
         endif
         grid%initial_neighbor(i,neigh_elem) = elem
! does not exist
      else
         face = face + 1
         grid%face(face)%vertex = temp_face
         grid%element(elem)%face(loc_face) = face
         grid%face(face)%assoc_elem = elem
         contains_face(face) = elem
         do i=1,VERTICES_PER_FACE
            ncontain(temp_face(i)) = ncontain(temp_face(i)) + 1
            if (ncontain(temp_face(i)) > MAX_CONTAINS) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("more faces contain a vertex than allowed in init_grid_msh")
               stop
            endif
            contains_vert(ncontain(temp_face(i)),temp_face(i)) = face
         end do
      endif
   end do ! next face
end do ! next element
deallocate(contains_face)

grid%nface = face
grid%biggest_face = grid%nface
if (set_iown) then
   grid%nface_own = grid%nface
else
   grid%nface_own = 0
endif
grid%next_free_face = grid%nface+1

! Create the edges by going through the faces and creating any of the
! three edges that don't already exist.  For each vertex, keep a list of the
! containing edges to check whether an edge has already been created.  Assign
! the edges to faces with the first edge opposite the first vertex, etc.
! Set the edge associated element to be the same element as the associated
! element of the first face that created the edge.

if (size(grid%edge) < 6*size(grid%vertex)) then
   call more_edges(grid,errcode,size_wanted=6*size(grid%vertex))
endif

ncontain = 0

edge = 0
do face=1,grid%nface
   do loc_edge=1,EDGES_PER_FACE
      i = 0
! vertices of this edge
      do loc_vert=1,VERTICES_PER_FACE
         if (loc_vert == loc_edge) cycle
         i = i+1
         temp_edge(i) = grid%face(face)%vertex(loc_vert)
      end do
! does the edge already exist?
      existing_edge = -1
      do i=1,ncontain(temp_edge(1))
         if ((grid%edge(contains_vert(i,temp_edge(1)))%vertex(1) == temp_edge(1) .and. &
              grid%edge(contains_vert(i,temp_edge(1)))%vertex(2) == temp_edge(2)) .or. &
             (grid%edge(contains_vert(i,temp_edge(1)))%vertex(1) == temp_edge(2) .and. &
              grid%edge(contains_vert(i,temp_edge(1)))%vertex(2) == temp_edge(1))) then
            existing_edge = contains_vert(i,temp_edge(1))
            exit
         endif
      end do
! exists
      if (existing_edge /= -1) then
         grid%face(face)%edge(loc_edge) = existing_edge
! does not exist
      else
         edge = edge + 1
         grid%edge(edge)%vertex = temp_edge
         grid%face(face)%edge(loc_edge) = edge
         grid%edge(edge)%assoc_elem = grid%face(face)%assoc_elem
         do i=1,2
            ncontain(temp_edge(i)) = ncontain(temp_edge(i)) + 1
            if (ncontain(temp_edge(i)) > MAX_CONTAINS) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("more edges contain a vertex than allowed in init_grid_msh")
               stop
            endif
            contains_vert(ncontain(temp_edge(i)),temp_edge(i)) = edge
         end do
      endif
   end do ! next edge
end do ! next face
deallocate(contains_vert,ncontain)

grid%nedge = edge
grid%biggest_edge = grid%nedge
if (set_iown) then
   grid%nedge_own = grid%nedge
else
   grid%nedge_own = 0
endif
grid%next_free_edge = grid%nedge+1

! finish defining the edges

do i=1,grid%nedge
   grid%edge(i)%gid = i
   call hash_insert(grid%edge(i)%gid,i,grid%edge_hash)
   grid%edge(i)%degree = degree
   grid%edge(i)%child = NO_CHILD
   nullify(grid%edge(i)%solution)
   nullify(grid%edge(i)%exact)
   nullify(grid%edge(i)%oldsoln)
   grid%edge(i)%next = i
end do

! Finish defining the faces.  The marked edge is the longest edge.  As a tie
! breaker, use the one with the smallest x coordinate.  As a second tie
! breaker, use the smallest y coordinate, then z, and finally the smallest
! index.

first_gid = grid%nface/3

do i=1,grid%nface
   grid%face(i)%gid = first_gid + i
   call hash_insert(grid%face(i)%gid,i,grid%face_hash)
   grid%face(i)%degree = degree
   grid%face(i)%child = NO_CHILD
   grid%face(i)%marked_edge = grid%face(i)%edge(1)
   longest = edge_length(grid,grid%face(i)%edge(1))
   xvert_longest = grid%vertex(grid%edge(grid%face(i)%edge(1))%vertex)%coord%x
   yvert_longest = grid%vertex(grid%edge(grid%face(i)%edge(1))%vertex)%coord%y
   zvert_longest = grid%vertex(grid%edge(grid%face(i)%edge(1))%vertex)%coord%z
   do j=2,EDGES_PER_FACE
      this_length = edge_length(grid,grid%face(i)%edge(j))
      change_to_this = .false.
      if (this_length > longest + small*this_length) then        ! longer length
         change_to_this = .true.
      elseif (abs(this_length-longest) <= small*this_length) then ! equal length
         xvert = grid%vertex(grid%edge(grid%face(i)%edge(j))%vertex)%coord%x
         if (minval(xvert) < &
             minval(xvert_longest)-small*this_length) then       ! smaller x
            change_to_this = .true.
         elseif (abs(minval(xvert)-minval(xvert_longest)) < &
                 small*this_length) then                         ! equal x
            yvert = grid%vertex(grid%edge(grid%face(i)%edge(j))%vertex)%coord%y
            if (minval(yvert) < &
                minval(yvert_longest)-small*this_length) then    ! smaller y
               change_to_this = .true.
            elseif (abs(minval(yvert)-minval(yvert_longest)) < &
                    small*this_length) then                      ! equal y
               zvert=grid%vertex(grid%edge(grid%face(i)%edge(j))%vertex)%coord%z
               if (minval(zvert) < &
                   minval(zvert_longest)-small*this_length) then ! smaller z
                  change_to_this = .true.
               elseif (abs(minval(zvert)-minval(zvert_longest)) < &
                       small*this_length) then                   ! equal z
                  if (grid%face(i)%edge(j) < &
                      grid%face(i)%marked_edge) then             ! smaller index
                     change_to_this = .true.
                  endif
               endif
            endif
         endif
      endif
      if (change_to_this) then
         grid%face(i)%marked_edge = grid%face(i)%edge(j)
         longest = this_length
         xvert_longest = &
            grid%vertex(grid%edge(grid%face(i)%edge(j))%vertex)%coord%x
         yvert_longest = &
            grid%vertex(grid%edge(grid%face(i)%edge(j))%vertex)%coord%y
         zvert_longest = &
            grid%vertex(grid%edge(grid%face(i)%edge(j))%vertex)%coord%z
      endif
   end do
   nullify(grid%face(i)%solution)
   nullify(grid%face(i)%exact)
   nullify(grid%face(i)%oldsoln)
   grid%face(i)%next = i
end do

! Finish defining the elements.  Assign edges in arbitrary order.
! The refinement edge is the longest edge, with the same ties breakers as
! the marked edges (see above).

first_gid = max(max(grid%biggest_vert/2,grid%nedge/4),grid%nelem)

do i=1,grid%nelem
   elem_nedge = 0
   do j=1,FACES_PER_ELEMENT
      face = grid%element(i)%face(j)
      do k=1,EDGES_PER_FACE
         do m=1,elem_nedge
            if (grid%element(i)%edge(m) == grid%face(face)%edge(k)) exit
         end do
         if (m > elem_nedge) then
            elem_nedge = elem_nedge+1
            grid%element(i)%edge(elem_nedge) = grid%face(face)%edge(k)
         endif
      end do
   end do
   grid%element(i)%gid = first_gid+i
   call hash_insert(grid%element(i)%gid,i,grid%elem_hash)
   grid%element(i)%degree = degree
   grid%element(i)%level = 1
   grid%element(i)%mate = BOUNDARY
   grid%element(i)%in = 1
   grid%element(i)%out = 1
   grid%element(i)%order = (/1,2/)
   grid%element(i)%weight = 0.0_my_real
   grid%element(i)%isleaf = .true.
   grid%element(i)%oldleaf = .false.
   grid%element(i)%iown = set_iown
   grid%element(i)%hrefined_unowned = .false.
   grid%element(i)%prefined_unowned = .false.
   grid%element(i)%sp_eta_pred = 0.0_my_real
   do j=1,NEIGHBORS_PER_ELEMENT
      if (grid%initial_neighbor(j,i) == BOUNDARY) then
         grid%element(i)%neighbor_hint(j) = BOUNDARY
      else
         grid%element(i)%neighbor_hint(j) = &
                       grid%element(grid%initial_neighbor(j,i))%gid
      endif
   end do
   grid%element(i)%refinement_edge = grid%element(i)%edge(1)
   longest = edge_length(grid,grid%element(i)%edge(1))
   xvert_longest = &
      grid%vertex(grid%edge(grid%element(i)%edge(1))%vertex)%coord%x
   yvert_longest = &
      grid%vertex(grid%edge(grid%element(i)%edge(1))%vertex)%coord%y
   zvert_longest = &
      grid%vertex(grid%edge(grid%element(i)%edge(1))%vertex)%coord%z
   do j=2,EDGES_PER_ELEMENT
      this_length = edge_length(grid,grid%element(i)%edge(j))
      change_to_this = .false.
      if (this_length > longest + small*this_length) then        ! longer length
         change_to_this = .true.
      elseif (abs(this_length-longest) <= small*this_length) then ! equal length
         xvert = grid%vertex(grid%edge(grid%element(i)%edge(j))%vertex)%coord%x
         if (minval(xvert) < &
             minval(xvert_longest)-small*this_length) then       ! smaller x
            change_to_this = .true.
         elseif (abs(minval(xvert)-minval(xvert_longest)) < &
                 small*this_length) then                         ! equal x
            yvert=grid%vertex(grid%edge(grid%element(i)%edge(j))%vertex)%coord%y
            if (minval(yvert) < &
                minval(yvert_longest)-small*this_length) then    ! smaller y
               change_to_this = .true.
            elseif (abs(minval(yvert)-minval(yvert_longest)) < &
                    small*this_length) then                      ! equal y
               zvert = &
                  grid%vertex(grid%edge(grid%element(i)%edge(j))%vertex)%coord%z
               if (minval(zvert) < &
                   minval(zvert_longest)-small*this_length) then ! smaller z
                  change_to_this = .true.
               elseif (abs(minval(zvert)-minval(zvert_longest)) < &
                       small*this_length) then                   ! equal z
                  if (grid%element(i)%edge(j) < &
                      grid%element(i)%refinement_edge) then      ! smaller index
                     change_to_this = .true.
                  endif
               endif
            endif
         endif
      endif
      if (change_to_this) then
         grid%element(i)%refinement_edge = grid%element(i)%edge(j)
         longest = this_length
         xvert_longest = &
            grid%vertex(grid%edge(grid%element(i)%edge(j))%vertex)%coord%x
         yvert_longest = &
            grid%vertex(grid%edge(grid%element(i)%edge(j))%vertex)%coord%y
         zvert_longest = &
            grid%vertex(grid%edge(grid%element(i)%edge(j))%vertex)%coord%z
      endif
   end do
   call determine_element_type(grid,i)
   nullify(grid%element(i)%solution)
   nullify(grid%element(i)%exact)
   nullify(grid%element(i)%oldsoln)
end do

grid%dof = grid%nvert
if (set_iown) then
   grid%dof_own = grid%dof
else
   grid%dof_own = 0
endif

! set vertex, line and face tags

call set_vertex_tags(grid,msh_points,nmsh_point)
deallocate(msh_points)
call set_edge_tags(grid,msh_lines,nmsh_line)
deallocate(msh_lines)
call set_face_tags(grid,msh_tris,nmsh_tri)
deallocate(msh_tris)

! set the boundary markers

call set_bmark(grid)

! set the b.c. type and Dirichlet solution for vertices, edges and faces

visited_vert = .false.
do i=1,grid%nelem
   do ivert=1,VERTICES_PER_ELEMENT
      vert = grid%element(i)%vertex(ivert)
      if (visited_vert(vert)) cycle
      visited_vert(vert) = .true.
      if (grid%vertex(vert)%bmark /= 0) then
         call my_bconds(grid%vertex(vert)%coord%x,grid%vertex(vert)%coord%y, &
                        grid%vertex(vert)%coord%z,grid%vertex(vert)%bmark, &
                        grid%vertex_type(vert,:),bccoef,bcrhs)
         grid%vertex_solution(vert,:,:) = 0
         do j=1,max(1,grid%num_eval)
            where (grid%vertex_type(vert,:) == DIRICHLET) &
               grid%vertex_solution(vert,:,j) = bcrhs
         end do
      else
         grid%vertex_type(vert,:) = INTERIOR
         grid%vertex_solution(vert,:,:) = 0
      endif
   end do
end do
deallocate(visited_vert)

do i=1,grid%nedge
   if (grid%edge(i)%bmark /= 0) then
      call my_bconds((grid%vertex(grid%edge(i)%vertex(1))%coord%x + &
                      grid%vertex(grid%edge(i)%vertex(2))%coord%x)/2, &
                     (grid%vertex(grid%edge(i)%vertex(1))%coord%y + &
                      grid%vertex(grid%edge(i)%vertex(2))%coord%y)/2, &
                     (grid%vertex(grid%edge(i)%vertex(1))%coord%z + &
                      grid%vertex(grid%edge(i)%vertex(2))%coord%z)/2, &
                     grid%edge(i)%bmark,grid%edge_type(i,:),bccoef,bcrhs)
! TEMP3D will need to set Dirichlet solution for high order elements
   else
      grid%edge_type(i,:) = INTERIOR
! TEMP3D will need to 0 solution for high order elements
   endif
end do

do i=1,grid%nface
   if (grid%face(i)%bmark /= 0) then
      call my_bconds((grid%vertex(grid%face(i)%vertex(1))%coord%x + &
                      grid%vertex(grid%face(i)%vertex(2))%coord%x + &
                      grid%vertex(grid%face(i)%vertex(3))%coord%x)/3, &
                     (grid%vertex(grid%face(i)%vertex(1))%coord%y + &
                      grid%vertex(grid%face(i)%vertex(2))%coord%y + &
                      grid%vertex(grid%face(i)%vertex(3))%coord%y)/3, &
                     (grid%vertex(grid%face(i)%vertex(1))%coord%z + &
                      grid%vertex(grid%face(i)%vertex(2))%coord%z + &
                      grid%vertex(grid%face(i)%vertex(3))%coord%z)/3, &
                     grid%face(i)%bmark,grid%face_type(i,:),bccoef,bcrhs)
! TEMP3D will need to set Dirichlet solution for high order elements
   else
      grid%face_type(i,:) = INTERIOR
! TEMP3D will need to 0 solution for high order elements
   endif
end do

! set up periodic boundary conditions

call setup_periodic(grid)

! and all the rest

grid%partition = partition
grid%errind_up2date = .false.

! define the boundary points, lines, lineloops and surfaces of the grid (if
! not already set by process_geo) and the boundary_vertex list of each vertex

call fix_gmsh_circle(grid) ! TEMP120323
if (grid%nsurface == 0) then
   call setup_polygonal_boundary(grid)
endif
call set_boundary_vertex(grid)

end subroutine init_grid_msh

! TEMP120323
subroutine fix_gmsh_circle(grid)
type(grid_type), intent(inout) :: grid
integer :: i,j
real(my_real) :: goodr,r
if (CYLINDERSIDE == -1) return ! indicates not cylinderHole
! for cylinderHole Gmsh 2.5.0 is placing some boundary points on a line
! between two neighboring points instead of on the circular surfaces.  This
! routine identifies faces with bmark corresponding to circular sides
! and makes sure all three vertices are on the circular side
! This may have been fixed in 2.6.0, but I'm not positive, so keep the kludge.
do i=1,grid%nface
   if (grid%face(i)%bmark == CYLINDERSIDE .or. &
       grid%face(i)%bmark == HOLESIDE) then
      do j=1,VERTICES_PER_FACE
         if (grid%face(i)%bmark == CYLINDERSIDE) then
            goodr = cylinderRadius
         else ! HOLESIDE
            goodr = holeBottomRadius + (holeTopRadius-holeBottomRadius) * &
                    grid%vertex(grid%face(i)%vertex(j))%coord%z/holeHeight
         endif
         r = sqrt(grid%vertex(grid%face(i)%vertex(j))%coord%x**2 + &
                  grid%vertex(grid%face(i)%vertex(j))%coord%y**2)
         if (abs(r-goodr) > small) then
            grid%vertex(grid%face(i)%vertex(j))%coord%x = &
               grid%vertex(grid%face(i)%vertex(j))%coord%x*goodr/r
            grid%vertex(grid%face(i)%vertex(j))%coord%y = &
               grid%vertex(grid%face(i)%vertex(j))%coord%y*goodr/r
         endif
      end do
   endif
end do
end subroutine fix_gmsh_circle
! end TEMP120323

!          --------
subroutine read_msh(grid,msh_filename,msh_points,nmsh_point,msh_lines, &
                    nmsh_line,msh_tris,nmsh_tri)
!          --------

!----------------------------------------------------
! This routine reads the Gmsh .msh file named msh_filename
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
character(len=*), intent(in) :: msh_filename
type(msh_point), pointer :: msh_points(:)
type(msh_line), pointer :: msh_lines(:)
type(msh_tri), pointer :: msh_tris(:)
integer, intent(out) :: nmsh_point,nmsh_line,nmsh_tri
!----------------------------------------------------
! Local variables:

character(len=80) :: inline
integer :: iunit, i, ind, msh_nelem, elem, point, line, tri, etype, ntag, &
           mtag, iostat, errcode, astat
real(my_real) :: x, y, z
real :: msh_version
!----------------------------------------------------
! Begin executable code

! open the file

iunit = get_lun()
open(unit=iunit,file=msh_filename,action="READ",status="OLD",iostat=iostat)
if (iostat /= 0) then
   ierr = USER_INPUT_ERROR
   call fatal("failed to open .msh file in read_msh")
   stop
endif

! format for reading an input line into the character string "inline"

100 format(A80)

! check the mesh format

read(iunit,100) inline
if (inline /= "$MeshFormat") then
   ierr = USER_INPUT_ERROR
   call fatal("initial grid file does not appear to be in .msh format")
   stop
endif

! get the version number, ignore file-type and data-size

read (iunit,*) msh_version
if (abs(msh_version - 2.2) > .0001) then
   call warning("msh file format is not 2.2; reading initial grid may fail")
endif

! check the next two lines are as expected

read(iunit,100) inline
if (inline /= "$EndMeshFormat") then
   ierr = USER_INPUT_ERROR
   call fatal("missing $EndMeshFormat statement in msh file")
   stop
endif
read(iunit,100) inline
if (inline /= "$Nodes") then
   ierr = USER_INPUT_ERROR
   call fatal("expecting $Nodes statement in msh file")
   stop
endif

! get the number of nodes and make sure grid has sufficient space for them

read(iunit,*) grid%nvert
if (size(grid%vertex) < grid%nvert) then
   call more_verts(grid,errcode,size_wanted=grid%nvert)
endif

! read the vertices.  require that no indices are skipped

do i=1,grid%nvert
   read(iunit,*) ind,x,y,z
   if (ind > grid%nvert) then
      ierr = USER_INPUT_ERROR
      call fatal("found a node-number greater than number of nodes in msh file")
      stop
   endif
   grid%vertex(ind)%coord%x = x
   grid%vertex(ind)%coord%y = y
   grid%vertex(ind)%coord%z = z
end do

! check the next two lines are as expected

read(iunit,100) inline
if (inline /= "$EndNodes") then
   ierr = USER_INPUT_ERROR
   call fatal("last node was not followed by $EndNodes statement in msh file")
   stop
endif
read(iunit,100) inline
if (inline /= "$Elements") then
   ierr = USER_INPUT_ERROR
   call fatal("expecting $Elements statement in msh file")
   stop
endif

! get the number of elements and make sure grid has sufficient space for them

read(iunit,*) msh_nelem
if (size(grid%element) < msh_nelem) then
   call more_elements(grid,errcode,size_wanted=msh_nelem)
endif
allocate(msh_points(msh_nelem),msh_lines(msh_nelem),msh_tris(msh_nelem), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in read_msh")
   stop
endif

! read the elements.  elements that are not tetrahedra are not an element
! in my grid; consequently the element IDs will be different than in the msh
! file.

point = 0
line = 0
tri  = 0
elem = 0
do i=1,msh_nelem
   read(iunit,100) inline
   read(inline,*) ind,etype,ntag
   if (ntag > MAX_TAG) then
      ierr = USER_INPUT_ERROR
      call fatal("Number of element tags in msh file is greater than MAX_TAG",&
                 "Increase MAX_TAG in gridtype.f90 or gridtype3D.f90 and rebuild the PHAML library")
      stop
   endif
   select case (etype)
! point
   case (15)
      point = point + 1
      read(inline,*) ind,etype,mtag,msh_points(point)%tags(1:ntag), &
                   msh_points(point)%vertex
      msh_points(point)%tags(ntag+1:MAX_TAG) = 0
! line
   case (1)
      line = line + 1
      read(inline,*) ind,etype,mtag,msh_lines(line)%tags(1:ntag), &
                   msh_lines(line)%vertex
      msh_lines(line)%tags(ntag+1:MAX_TAG) = 0
! triangle
   case (2)
      tri = tri + 1
      read(inline,*) ind,etype,mtag,msh_tris(tri)%tags(1:ntag), &
                   msh_tris(tri)%vertex
      msh_tris(tri)%tags(ntag+1:MAX_TAG) = 0
! tetrahedron
   case (4)
      elem = elem + 1
      read(inline,*) ind,etype,mtag,grid%element(elem)%tags(1:ntag), &
                   grid%element(elem)%vertex
      grid%element(elem)%tags(ntag+1:MAX_TAG) = 0
! ignore all others
   case default
   end select
end do
nmsh_point = point
nmsh_line = line
nmsh_tri = tri
grid%nelem = elem

! check the next line is as expected and ignore any additional sections
! TEMP add reading NodeData section

read(iunit,100) inline
if (inline /= "$EndElements") then
   call warning("expecting $EndElements statement in msh file")
endif

close(iunit)

end subroutine read_msh

!          ----------------------
subroutine determine_element_type(grid,elem)
!          ----------------------

!----------------------------------------------------
! This routine determines the type of element elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
!----------------------------------------------------
! Local variables:

integer :: i, marked_edge(FACES_PER_ELEMENT), refinement_edge, vert, &
           num_intersect
!----------------------------------------------------
! Begin executable code

! flag to indicate the type has not yet been determined

grid%element(elem)%type = "?"

! get the marked edges and refinement edge for convenience

do i=1,FACES_PER_ELEMENT
   marked_edge(i) = grid%face(grid%element(elem)%face(i))%marked_edge
end do
refinement_edge = grid%element(elem)%refinement_edge

! it is type P (planar) if one of the element vertices is not included in
! any of the marked edges.  If it is type P, it is unflagged, i.e. Pu

do i=1,VERTICES_PER_ELEMENT
   vert = grid%element(elem)%vertex(i)
   if (.not. (grid%edge(marked_edge(1))%vertex(1) == vert .or. &
              grid%edge(marked_edge(1))%vertex(2) == vert .or. &
              grid%edge(marked_edge(2))%vertex(1) == vert .or. &
              grid%edge(marked_edge(2))%vertex(2) == vert .or. &
              grid%edge(marked_edge(3))%vertex(1) == vert .or. &
              grid%edge(marked_edge(3))%vertex(2) == vert .or. &
              grid%edge(marked_edge(4))%vertex(1) == vert .or. &
              grid%edge(marked_edge(4))%vertex(2) == vert)) then
      grid%element(elem)%type = "Pu"
      return
   endif
end do

! determine the number of marked edges that intersect the refinement edge

num_intersect = 0
do i=1,FACES_PER_ELEMENT
   if (grid%edge(marked_edge(i))%vertex(1) == grid%edge(refinement_edge)%vertex(1) .or. &
       grid%edge(marked_edge(i))%vertex(1) == grid%edge(refinement_edge)%vertex(2) .or. &
       grid%edge(marked_edge(i))%vertex(2) == grid%edge(refinement_edge)%vertex(1) .or. &
       grid%edge(marked_edge(i))%vertex(2) == grid%edge(refinement_edge)%vertex(2)) then
      num_intersect = num_intersect + 1
   endif
end do

select case (num_intersect)

case(1)

   ierr = PHAML_INTERNAL_ERROR
   call fatal("determine_element_type: only one marked edge intersects the refinement edge")
   stop

case(2)

! it is type O (opposite) if two edges are opposite the refinement edge

   grid%element(elem)%type = "O"

case(3)

! it is type M (mixed) if just one of the non-refinement marked edges intersects
! the refinement edge

   grid%element(elem)%type = "M"

case(4)

! it is type A (adjacent) if all marked edges intersect the refinement edge

   grid%element(elem)%type = "A"

end select

end subroutine determine_element_type

!        -----------
function edge_length(grid,edge)
!        -----------

!----------------------------------------------------
! This routine computes the length of edge
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: edge
real(my_real) :: edge_length
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

edge_length = sqrt( (grid%vertex(grid%edge(edge)%vertex(1))%coord%x - &
                     grid%vertex(grid%edge(edge)%vertex(2))%coord%x)**2 + &
                    (grid%vertex(grid%edge(edge)%vertex(1))%coord%y - &
                     grid%vertex(grid%edge(edge)%vertex(2))%coord%y)**2 + &
                    (grid%vertex(grid%edge(edge)%vertex(1))%coord%z - &
                     grid%vertex(grid%edge(edge)%vertex(2))%coord%z)**2 )

end function edge_length

!          ---------------
subroutine set_vertex_tags(grid,msh_points,nmsh_point)
!          ---------------

!----------------------------------------------------
! This routine copies the tags in msh_points to the corresponding vertices in
! grid.
! The vertex IDs in grid must be the same as those in msh_points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(msh_point), intent(in) :: msh_points(:)
integer, intent(in) :: nmsh_point
!----------------------------------------------------
! Local variables:

integer :: point, i
!----------------------------------------------------
! Begin executable code

! Missing tags are set to 0

do i=1,grid%biggest_vert
   grid%vertex(i)%tags = 0
end do

! if there aren't any mesh points, nothing to do here

if (nmsh_point == 0) return

! Go through the points of msh_points and copy the tags.

do point=1,nmsh_point
   grid%vertex(msh_points(point)%vertex)%tags = msh_points(point)%tags
end do

end subroutine set_vertex_tags

!          -------------
subroutine set_edge_tags(grid,msh_lines,nmsh_line)
!          -------------

!----------------------------------------------------
! This routine copies the tags in msh_lines to the corresponding edges in grid.
! The vertex IDs in grid must be the same as those in msh_tris.
! The edges in grid must fully occupy grid%edge(1:grid%nedge)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(msh_line), intent(in) :: msh_lines(:)
integer, intent(in) :: nmsh_line
!----------------------------------------------------
! Local variables:

integer :: i, errcode, msh_index, edge, match
integer :: msh_key(nmsh_line), msh_perm(nmsh_line), edge_key(grid%nedge), &
           edge_perm(grid%nedge)
!----------------------------------------------------
! Begin executable code

! tags that aren't in msh_lines are 0

do i=1,grid%nedge
   grid%edge(i)%tags = 0
end do

! if there aren't any mesh lines, there is nothing to do here

if (nmsh_line == 0) return

! sort the msh_lines by the smallest vertex index

do i=1,nmsh_line
   msh_key(i) = minval(msh_lines(i)%vertex)
end do
call sort(msh_key,nmsh_line,msh_perm,1,errcode)

! likewise sort the edges

do i=1,grid%nedge
   edge_key(i) = minval(grid%edge(i)%vertex)
end do
call sort(edge_key,grid%nedge,edge_perm,1,errcode)

! Go through the edges of grid in sorted order to set the tags.  Follow	
! along with msh_lines looking for a matching line.

msh_index = 1
do edge=1,grid%nedge

! move up in the msh_lines until the sorted index agrees or exceeds edge's

   do while (msh_key(msh_perm(msh_index)) < edge_key(edge_perm(edge)))
      msh_index = msh_index + 1
      if (msh_index > nmsh_line) exit
   end do

! if msh_index overflows, there are no more tags to copy

   if (msh_index > nmsh_line) exit

! if there is a matching msh_lines for this edge, copy the tags

   match = 0
   i = msh_index
   do while (msh_key(msh_perm(i)) == msh_key(msh_perm(msh_index)))
      if (lines_match(grid%edge(edge_perm(edge))%vertex, &
                      msh_lines(msh_perm(i))%vertex)) then
         match = msh_perm(i)
         exit
      endif
      i = i+1
      if (i > nmsh_line) exit
   end do

   if (match /= 0) grid%edge(edge_perm(edge))%tags = msh_lines(match)%tags

end do ! next edge

end subroutine set_edge_tags

!        -----------
function lines_match(l1,l2)
!        -----------

!----------------------------------------------------
! This routine returns true if the vertices of line l1 are the
! same as the vertices of line l2
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: l1(2), l2(2)
logical :: lines_match
!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

lines_match = (l1(1)==l2(1) .and. l1(2)==l2(2)) .or. &
              (l1(1)==l2(2) .and. l1(2)==l2(1))

end function lines_match

!          -------------
subroutine set_face_tags(grid,msh_tris,nmsh_tri)
!          -------------

!----------------------------------------------------
! This routine copies the tags in msh_tris to the corresponding faces in grid.
! The vertex IDs in grid must be the same as those in msh_tris.
! The faces in grid must fully occupy grid%face(1:grid%nface)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(msh_tri), intent(in) :: msh_tris(:)
integer, intent(in) :: nmsh_tri
!----------------------------------------------------
! Local variables:

integer :: i, errcode, msh_index, face, match
integer :: msh_key(nmsh_tri), msh_perm(nmsh_tri), face_key(grid%nface), &
           face_perm(grid%nface)
!----------------------------------------------------
! Begin executable code

! tags that aren't in msh_tris are 0

do i=1,grid%nface
   grid%face(i)%tags = 0
end do

! if there aren't any mesh faces, nothing to do here

if (nmsh_tri == 0) return

! sort the msh_tris by the smallest vertex index

do i=1,nmsh_tri
   msh_key(i) = minval(msh_tris(i)%vertex)
end do
call sort(msh_key,nmsh_tri,msh_perm,1,errcode)

! likewise sort the faces

do i=1,grid%nface
   face_key(i) = minval(grid%face(i)%vertex)
end do
call sort(face_key,grid%nface,face_perm,1,errcode)

! Go through the faces of grid in sorted order to set the tags.  Follow	
! along with msh_tris looking for a matching triangle.

msh_index = 1
do face=1,grid%nface

! move up in the msh_tris until the sorted index agrees or exceeds face's

   do while (msh_key(msh_perm(msh_index)) < face_key(face_perm(face)))
      msh_index = msh_index + 1
      if (msh_index > nmsh_tri) exit
   end do

! if msh_index overflows, there are no more tags to copy

   if (msh_index > nmsh_tri) exit

! if there is a matching msh_tris for this face, copy the tags

   match = 0
   i = msh_index
   do while (msh_key(msh_perm(i)) == msh_key(msh_perm(msh_index)))
      if (triangles_match(grid%face(face_perm(face))%vertex, &
                          msh_tris(msh_perm(i))%vertex)) then
         match = msh_perm(i)
         exit
      endif
      i = i+1
      if (i > nmsh_tri) exit
   end do

   if (match /= 0) grid%face(face_perm(face))%tags = msh_tris(match)%tags

end do ! next face

end subroutine set_face_tags

!        ---------------
function triangles_match(t1,t2)
!        ---------------

!----------------------------------------------------
! This routine returns true if the vertices of triangle t1 are the
! same as the vertices of triangle t2
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: t1(3), t2(3)
logical :: triangles_match
!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

triangles_match = .true.

do i=1,3
   if (all(t2 /= t1(i))) then
      triangles_match = .false.
      exit
   endif
end do

end function triangles_match

!          ---------
subroutine set_bmark(grid)
!          ---------

!----------------------------------------------------
! This routine sets bmark.
!
! bmark is set to 0 for interior vertices, edges and faces.
!
! If a boundary vertex, edge or face is in the .msh file with tags, the first
! tag is used as bmark.  By default, this is the number of the physical region
! containing the entity.  It is possible to assign all entities to physical
! regions, giving complete control over bmarks.
! 
! If a boundary face is not in the .msh file or has no tags, it is assigned
! bmark=the smallest positive integer that is not a first tag.
! 
! If a boundary edge is not in the .msh file or has no tags, it is assigned the
! bmark of one of the two faces that share the edge.  If one of the faces has
! bmark=the smallest positive integer that is not a first tag, the edge gets the
! bmark of the other face (which could be the smallest positive integer that is
! not a first tag).
! Otherwise, if one face has Dirichlet boundary conditions and the other does
! not, the edge gets the bmark of the Dirichlet face.  Otherwise, it gets the
! smaller of the two bmarks.
! 
! If a boundary vertex is not in the .msh file or has no tags, it is assigned
! the bmark of one of the edges that contains it.  If any edges have Dirichlet
! boundary conditions, it is assigned the smallest bmark of those edges.
! Otherwise, it is assigned the smallest bmark of all the edges.
!
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: i, default_bmark, elem, iface, face, iedge, edge, ivert, vert
integer :: containing_face(2,grid%nedge), face_type(2,grid%system_size), &
           containing_edge(MAX_CONTAINS,grid%biggest_vert), &
           ncontain(grid%biggest_vert), &
           edge_type(MAX_CONTAINS,grid%system_size)
real(my_real) :: bccoef(grid%system_size,grid%system_size), &
                 bcrhs(grid%system_size)
logical :: visited_vert(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

! find the smallest positive integer that is not a first tag

if (MAX_TAG < 1) then
   default_bmark = 1
else
   default_bmark = 1
   do
      if (any(grid%vertex(1:grid%biggest_vert)%tags(1)==default_bmark)) then
         default_bmark = default_bmark+1
         cycle
      endif
      if (any(grid%edge(1:grid%nedge)%tags(1)==default_bmark)) then
         default_bmark = default_bmark+1
         cycle
      endif
      if (any(grid%face(1:grid%nface)%tags(1)==default_bmark)) then
         default_bmark = default_bmark+1
         cycle
      endif
      exit
   end do
endif

! set all bmarks to 0, indicating an interior entity, and then identify the
! boundary entities

grid%vertex(1:grid%biggest_vert)%bmark = 0
grid%edge(1:grid%nedge)%bmark = 0
grid%face(1:grid%nface)%bmark = 0

! a face is boundary if and only if there is an element for which the
! neighbor across this face is BOUNDARY.  Assign boundary faces the first
! tag or, if the first tag is 0, the default bmark.


do elem=1,grid%nelem
   do iface=1,FACES_PER_ELEMENT
! TEMP120313 do all faces, to allow for interior interfaces
!      if (.TRUE. .OR. grid%initial_neighbor(iface,elem) == BOUNDARY) then
      if (grid%initial_neighbor(iface,elem) == BOUNDARY) then
         face = grid%element(elem)%face(iface)
         if (MAX_TAG < 1) then
            grid%face(face)%bmark = default_bmark
! TEMP120313 even more TEMP, don't change a 0 bmark
!         elseif (.FALSE. .AND. grid%face(face)%tags(1) == 0) then
         elseif (grid%face(face)%tags(1) == 0) then
            grid%face(face)%bmark = default_bmark
         else
            grid%face(face)%bmark = grid%face(face)%tags(1)
         endif
      endif
   end do
end do

! an edge is boundary if and only if it is an edge of a boundary face.  If
! the first tag is 0, keep track of the faces and in a second pass set the
! bmark from one of the faces.

containing_face = 0
do face=1,grid%nface
   if (grid%face(face)%bmark == 0) cycle
   do iedge=1,EDGES_PER_FACE
      edge = grid%face(face)%edge(iedge)
      if (MAX_TAG < 1) then
         grid%edge(edge)%bmark = default_bmark
         if (containing_face(1,edge) == 0) then
            containing_face(1,edge) = face
         else
            containing_face(2,edge) = face
         endif
      elseif (grid%edge(edge)%tags(1) == 0) then
         grid%edge(edge)%bmark = default_bmark
         if (containing_face(1,edge) == 0) then
            containing_face(1,edge) = face
         else
            containing_face(2,edge) = face
         endif
      else
         grid%edge(edge)%bmark = grid%edge(edge)%tags(1)
      endif
   end do
end do

! apply the rule contained in the opening comments

do edge=1,grid%nedge
   if (grid%edge(edge)%bmark == default_bmark) then
      if (grid%face(containing_face(1,edge))%bmark == default_bmark) then
         grid%edge(edge)%bmark = grid%face(containing_face(2,edge))%bmark
      elseif (grid%face(containing_face(2,edge))%bmark == default_bmark) then
         grid%edge(edge)%bmark = grid%face(containing_face(1,edge))%bmark
      else
         do i=1,2
            face = containing_face(i,edge)
            call my_bconds((grid%vertex(grid%face(face)%vertex(1))%coord%x + &
                            grid%vertex(grid%face(face)%vertex(2))%coord%x + &
                            grid%vertex(grid%face(face)%vertex(3))%coord%x)/3, &
                           (grid%vertex(grid%face(face)%vertex(1))%coord%y + &
                            grid%vertex(grid%face(face)%vertex(2))%coord%y + &
                            grid%vertex(grid%face(face)%vertex(3))%coord%y)/3, &
                           (grid%vertex(grid%face(face)%vertex(1))%coord%z + &
                            grid%vertex(grid%face(face)%vertex(2))%coord%z + &
                            grid%vertex(grid%face(face)%vertex(3))%coord%z)/3, &
                           grid%face(face)%bmark,face_type(i,:),bccoef,bcrhs)
         end do
         if (any(face_type(1,:) == DIRICHLET) .and. &
             all(face_type(2,:) /= DIRICHLET)) then
            grid%edge(edge)%bmark = grid%face(containing_face(1,edge))%bmark
         elseif (any(face_type(2,:) == DIRICHLET) .and. &
                 all(face_type(1,:) /= DIRICHLET)) then
            grid%edge(edge)%bmark = grid%face(containing_face(2,edge))%bmark
         else
            grid%edge(edge)%bmark = min( &
                                     grid%face(containing_face(1,edge))%bmark, &
                                     grid%face(containing_face(2,edge))%bmark)
         endif
      endif
   endif
end do

! a vertex is boundary if and only if it is a vertex of a boundary edge.  If
! the first tag is 0, keep track of the edges and in a second pass set the
! bmark from one of the edges.

ncontain = 0
do edge=1,grid%nedge
   if (grid%edge(edge)%bmark == 0) cycle
   do ivert=1,2
      vert = grid%edge(edge)%vertex(ivert)
      if (MAX_TAG < 1) then
         grid%vertex(vert)%bmark = default_bmark
         ncontain(vert) = ncontain(vert) + 1
         if (ncontain(vert) > MAX_CONTAINS) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("more edges contain a vertex than allowed in set_bmark")
            stop
         endif
         containing_edge(ncontain(vert),vert) = edge
      elseif (grid%vertex(vert)%tags(1) == 0) then
         grid%vertex(vert)%bmark = default_bmark
         ncontain(vert) = ncontain(vert) + 1
         if (ncontain(vert) > MAX_CONTAINS) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("more edges contain a vertex than allowed in set_bmark")
            stop
         endif
         containing_edge(ncontain(vert),vert) = edge
      else
         grid%vertex(vert)%bmark = grid%vertex(vert)%tags(1)
      endif
   end do
end do

! apply the rule contained in the opening comments

visited_vert = .false.
do elem=1,grid%nelem
   do ivert=1,VERTICES_PER_ELEMENT
      vert = grid%element(elem)%vertex(ivert)
      if (visited_vert(vert)) cycle
      visited_vert(vert) = .true.
      if (grid%vertex(vert)%bmark == default_bmark) then
         do i=1,ncontain(vert)
            edge = containing_edge(i,vert)
            call my_bconds((grid%vertex(grid%edge(edge)%vertex(1))%coord%x + &
                            grid%vertex(grid%edge(edge)%vertex(2))%coord%x)/2, &
                           (grid%vertex(grid%edge(edge)%vertex(1))%coord%y + &
                            grid%vertex(grid%edge(edge)%vertex(2))%coord%y)/2, &
                           (grid%vertex(grid%edge(edge)%vertex(1))%coord%z + &
                            grid%vertex(grid%edge(edge)%vertex(2))%coord%z)/2, &
                           grid%edge(edge)%bmark,edge_type(i,:),bccoef,bcrhs)
         end do
         grid%vertex(vert)%bmark = huge(0)
         do i=1,ncontain(vert)
            if (any(edge_type(i,:) == DIRICHLET)) then
               grid%vertex(vert)%bmark = min(grid%vertex(vert)%bmark, &
                                       grid%edge(containing_edge(i,vert))%bmark)
            endif
         end do
         if (grid%vertex(vert)%bmark == huge(0)) then
            do i=1,ncontain(vert)
               grid%vertex(vert)%bmark = min(grid%vertex(vert)%bmark, &
                                       grid%edge(containing_edge(i,vert))%bmark)
            end do
         endif
      endif
   end do
end do

end subroutine set_bmark

!          --------------
subroutine setup_periodic(grid)
!          --------------

!----------------------------------------------------
! This routine sets up periodic boundary conditions.  PHAML only supports
! periodic boundary conditions on a parallelepiped with edges parallel to
! the axes.  The periodic faces (i.e. the faces on the two sides with
! x=constant, y=constant and/or z=constant) must have a bmark for which
! bconds returns PERIODIC.  The edges and vertices on the boundary of those
! sides may be periodic or of the same type as another face that shares
! that edge or vertex.
!
! On input, the type of all vertices, edges and faces have already been set
! from calling bconds, but those with periodic boundary conditions are just
! labeled PERIODIC.  On output, PERIODIC will be replaced by PERIODIC_MASTER or
! PERIODIC_SLAVE.  Edges and vertices that are not PERIODIC (e.g. DIRICHLET)
! but are on a PERIODIC face will be relabled as partially periodic, e.g.
! PERIODIC_MASTER_DIR.  The linked list of vertices will be rearranged so that
! periodic masters are preceded immediately by their periodic slaves.  next in
! edges and faces is set such that periodic masters and slaves point to each
! other, forming a loop if there is more than one slave (doubly or triply
! periodic boundary conditions).  initial_neighbor for elements on a periodic
! boundary will be changed from BOUNDARY to the matching element.  The marked
! edges and refinement edge will be changed on the slave side to agree with
! the master.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer, parameter :: XPERIODIC=1, YPERIODIC=2, ZPERIODIC=3
integer :: first_face, last_face, next_face(grid%biggest_face), &
           prev_face(grid%biggest_face), face, edges(EDGES_PER_FACE), &
           verts(VERTICES_PER_FACE), periodicity, mate, v, vmate, &
           edges_mate(EDGES_PER_FACE), verts_mate(VERTICES_PER_FACE), e, f, &
           edge, v1, v2, emate, comp, last, lev, edge_mate
logical :: found_one, found_all, doubly_periodic1, doubly_periodic2, &
           did_edge(grid%biggest_edge,3), did_vert(grid%biggest_vert,3)
!----------------------------------------------------
! Begin executable code

! the "did"s indicate we already did an edge or vertex in a particular
! coordinate direction, to avoid duplication

did_edge = .false.
did_vert = .false.

! create a doubly linked list of all periodic faces

first_face = -1
last_face = -1

grid%any_periodic = .false.
do face=1,grid%nface
   if (any(grid%face_type(face,:) == PERIODIC)) then
      call add_to_linked_list(face,first_face,last_face,next_face,prev_face)
      grid%any_periodic = .true.
   endif
end do

! if no periodic faces were found, we're done

if (.not. grid%any_periodic) return

! For each periodic face, find the mate and set them up with the master as
! the one with the minumum value of the periodic coordinate.
! For each edge of a master periodic face, find the mate on the matching slave
! face, setup those two edges, and remove them from the list of periodic edges.
! Note that if an edge has already been set up, then we have a doubly periodic
! situation.  Likewise for the vertices.

do while (first_face /= -1)

! next periodic face

   face = first_face

! edges and vertices of the periodic face

   edges = grid%face(face)%edge
   verts = grid%face(face)%vertex

! determine the coordinate of periodicity

   if (grid%vertex(verts(1))%coord%x == grid%vertex(verts(2))%coord%x .and. &
       grid%vertex(verts(1))%coord%x == grid%vertex(verts(3))%coord%x) then
      periodicity = XPERIODIC
   elseif (grid%vertex(verts(1))%coord%y == grid%vertex(verts(2))%coord%y .and.&
           grid%vertex(verts(1))%coord%y == grid%vertex(verts(3))%coord%y) then
      periodicity = YPERIODIC
   elseif (grid%vertex(verts(1))%coord%z == grid%vertex(verts(2))%coord%z .and.&
           grid%vertex(verts(1))%coord%z == grid%vertex(verts(3))%coord%z) then
      periodicity = ZPERIODIC
   else
      ierr = USER_INPUT_ERROR
      call fatal("a face marked as periodic is not parallel to one of the coordinate planes")
      stop
   endif

! find the matching face, and set the matching vertices

   mate = next_face(face)
   do while (mate /= -1)
      found_all = .true.
      do v=1,VERTICES_PER_FACE
         found_one = .false.
         do vmate=1,VERTICES_PER_FACE
            select case (periodicity)
            case (XPERIODIC)
               if (myeq(grid%vertex(grid%face(face)%vertex(v))%coord%y, &
                   grid%vertex(grid%face(mate)%vertex(vmate))%coord%y) .and. &
                   myeq(grid%vertex(grid%face(face)%vertex(v))%coord%z, &
                   grid%vertex(grid%face(mate)%vertex(vmate))%coord%z)) then
                  found_one = .true.
                  verts_mate(v) = grid%face(mate)%vertex(vmate)
                  exit
               endif
            case (YPERIODIC)
               if (myeq(grid%vertex(grid%face(face)%vertex(v))%coord%x, &
                   grid%vertex(grid%face(mate)%vertex(vmate))%coord%x) .and. &
                   myeq(grid%vertex(grid%face(face)%vertex(v))%coord%z, &
                   grid%vertex(grid%face(mate)%vertex(vmate))%coord%z)) then
                  found_one = .true.
                  verts_mate(v) = grid%face(mate)%vertex(vmate)
                  exit
               endif
            case (ZPERIODIC)
               if (myeq(grid%vertex(grid%face(face)%vertex(v))%coord%x, &
                   grid%vertex(grid%face(mate)%vertex(vmate))%coord%x) .and. &
                   myeq(grid%vertex(grid%face(face)%vertex(v))%coord%y, &
                   grid%vertex(grid%face(mate)%vertex(vmate))%coord%y)) then
                  found_one = .true.
                  verts_mate(v) = grid%face(mate)%vertex(vmate)
                  exit
               endif
            end select
         end do
         if (.not. found_one) then
            found_all = .false.
            exit
         endif
      end do
      if (found_all) exit
      mate = next_face(mate)
   end do

   if (mate == -1) then
      ierr = USER_INPUT_ERROR
      call fatal("could not find a matching face for a periodic face")
      stop
   endif

! matching edges

   do e=1,EDGES_PER_FACE
      edge = grid%face(face)%edge(e)
      if (grid%edge(edge)%vertex(1) == verts(1)) then
         v1 = 1
      elseif (grid%edge(edge)%vertex(1) == verts(2)) then
         v1 = 2
      elseif (grid%edge(edge)%vertex(1) == verts(3)) then
         v1 = 3
      else
         ierr = PHAML_INTERNAL_ERROR
         call fatal("setup_periodic: failed to find vertex on face")
         stop
      endif
      if (grid%edge(edge)%vertex(2) == verts(1)) then
         v2 = 1
      elseif (grid%edge(edge)%vertex(2) == verts(2)) then
         v2 = 2
      elseif (grid%edge(edge)%vertex(2) == verts(3)) then
         v2 = 3
      else
         ierr = PHAML_INTERNAL_ERROR
         call fatal("setup_periodic: failed to find vertex on face")
         stop
      endif
      edges_mate(e) = -1
      do emate=1,EDGES_PER_FACE
         edge_mate = grid%face(mate)%edge(emate)
         if ((grid%edge(edge_mate)%vertex(1) == verts_mate(v1) .and. &
              grid%edge(edge_mate)%vertex(2) == verts_mate(v2)) .or. &
             (grid%edge(edge_mate)%vertex(2) == verts_mate(v1) .and. &
              grid%edge(edge_mate)%vertex(1) == verts_mate(v2))) then
            edges_mate(e) = edge_mate
            exit
         endif
      end do
      if (edges_mate(e) == -1) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("setup_periodic: could not find a matching edge")
         stop
      endif
   end do

! set up the faces

   select case (periodicity)
   case (XPERIODIC)
      if (grid%vertex(verts(1))%coord%x < &
          grid%vertex(verts_mate(1))%coord%x) then
         where (grid%face_type(face,:) == PERIODIC) &
            grid%face_type(face,:) = PERIODIC_MASTER
         where (grid%face_type(mate,:) == PERIODIC) &
            grid%face_type(mate,:) = PERIODIC_SLAVE
      else
         where (grid%face_type(face,:) == PERIODIC) &
            grid%face_type(face,:) = PERIODIC_SLAVE
         where (grid%face_type(mate,:) == PERIODIC) &
            grid%face_type(mate,:) = PERIODIC_MASTER
      endif
   case (YPERIODIC)
      if (grid%vertex(verts(1))%coord%y < &
          grid%vertex(verts_mate(1))%coord%y) then
         where (grid%face_type(face,:) == PERIODIC) &
            grid%face_type(face,:) = PERIODIC_MASTER
         where (grid%face_type(mate,:) == PERIODIC) &
            grid%face_type(mate,:) = PERIODIC_SLAVE
      else
         where (grid%face_type(face,:) == PERIODIC) &
            grid%face_type(face,:) = PERIODIC_SLAVE
         where (grid%face_type(mate,:) == PERIODIC) &
            grid%face_type(mate,:) = PERIODIC_MASTER
      endif
   case (ZPERIODIC)
      if (grid%vertex(verts(1))%coord%z < &
          grid%vertex(verts_mate(1))%coord%z) then
         where (grid%face_type(face,:) == PERIODIC) &
            grid%face_type(face,:) = PERIODIC_MASTER
         where (grid%face_type(mate,:) == PERIODIC) &
            grid%face_type(mate,:) = PERIODIC_SLAVE
      else
         where (grid%face_type(face,:) == PERIODIC) &
            grid%face_type(face,:) = PERIODIC_SLAVE
         where (grid%face_type(mate,:) == PERIODIC) &
            grid%face_type(mate,:) = PERIODIC_MASTER
      endif
   end select
   grid%face(face)%next = mate
   grid%face(mate)%next = face

! set up the edges

   do e=1,EDGES_PER_FACE
      if (did_edge(edges(e),periodicity)) cycle
      do comp=1,grid%system_size

! if this component of the face is not periodic, leave the edge alone

         if (.not. is_periodic_face(face,grid,comp)) cycle

! determine if the edge is doubly periodic and has already been assigned

         select case (grid%edge_type(edges(e),comp))
         case(DIRICHLET,NATURAL,MIXED,PERIODIC)
            doubly_periodic1 = .false.
         case(PERIODIC_MASTER,PERIODIC_SLAVE,PERIODIC_MASTER_DIR, &
              PERIODIC_SLAVE_DIR,PERIODIC_MASTER_NAT,PERIODIC_SLAVE_NAT, &
              PERIODIC_MASTER_MIX,PERIODIC_SLAVE_MIX )
            doubly_periodic1 = .true.
         case default
            ierr = PHAML_INTERNAL_ERROR
            call fatal("setup_periodic: unrecognized edge type")
            stop
         end select

! determine if the mate is doubly periodic and has already been assigned

         select case (grid%edge_type(edges_mate(e),comp))
         case(DIRICHLET,NATURAL,MIXED,PERIODIC)
            doubly_periodic2 = .false.
         case(PERIODIC_MASTER,PERIODIC_SLAVE,PERIODIC_MASTER_DIR, &
              PERIODIC_SLAVE_DIR,PERIODIC_MASTER_NAT,PERIODIC_SLAVE_NAT, &
              PERIODIC_MASTER_MIX,PERIODIC_SLAVE_MIX )
            doubly_periodic2 = .true.
         case default
            ierr = PHAML_INTERNAL_ERROR
            call fatal("setup_periodic: unrecognized edge type")
            stop
         end select

! if neither is doubly periodic, set the periodic type based on the type
! of the edge and whether the face is the master

         if (.not. doubly_periodic1 .and. .not. doubly_periodic2) then

            if (grid%face_type(face,comp) == PERIODIC_MASTER) then
               select case (grid%edge_type(edges(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
               select case (grid%edge_type(edges_mate(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
            else
               select case (grid%edge_type(edges(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
               select case (grid%edge_type(edges_mate(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
            endif

! if edge is doubly periodic but mate is not, set mate's type and set the
! edge to slave only if the mate is master

         elseif (doubly_periodic1 .and. .not. doubly_periodic2) then

            if (grid%face_type(face,comp) == PERIODIC_MASTER) then
               select case (grid%edge_type(edges_mate(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
            else
               select case (grid%edge_type(edges_mate(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
               select case (grid%edge_type(edges(e),comp))
               case(PERIODIC_MASTER_DIR, PERIODIC_SLAVE_DIR)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_DIR
               case(PERIODIC_MASTER_NAT, PERIODIC_SLAVE_NAT)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_NAT
               case(PERIODIC_MASTER_MIX, PERIODIC_SLAVE_MIX)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC_MASTER, PERIODIC_SLAVE)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
            endif

! if mate is doubly periodic but edge is not, leave the type for mate alone
! and set edge's type

         elseif (.not. doubly_periodic1 .and. doubly_periodic2) then

            if (grid%face_type(face,comp) == PERIODIC_MASTER) then
               select case (grid%edge_type(edges(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%edge_type(edges(e),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
               select case (grid%edge_type(edges_mate(e),comp))
               case(PERIODIC_MASTER_DIR, PERIODIC_SLAVE_DIR)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_DIR
               case(PERIODIC_MASTER_NAT, PERIODIC_SLAVE_NAT)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_NAT
               case(PERIODIC_MASTER_MIX, PERIODIC_SLAVE_MIX)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC_MASTER, PERIODIC_SLAVE)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
            else
               select case (grid%edge_type(edges(e),comp))
               case(DIRICHLET)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized edge type")
                  stop
               end select
            endif

! if both are doubly periodic ...

         else

! if they are both masters, remove the masterness of one of them depending on
! the masterness of face

            if (is_periodic_edge_master(edges(e),grid,comp) .and. &
                is_periodic_edge_master(edges_mate(e),grid,comp)) then

               if (grid%face_type(face,comp) == PERIODIC_MASTER) then
                  select case (grid%edge_type(edges_mate(e),comp))
                  case (PERIODIC_MASTER_DIR)
                     grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_DIR
                  case (PERIODIC_MASTER_NAT)
                     grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_NAT
                  case (PERIODIC_MASTER_MIX)
                     grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_MIX
                  case (PERIODIC_MASTER)
                     grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE
                  end select
               else
                  select case (grid%edge_type(edges(e),comp))
                  case (PERIODIC_MASTER_DIR)
                     grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_DIR
                  case (PERIODIC_MASTER_NAT)
                     grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_NAT
                  case (PERIODIC_MASTER_MIX)
                     grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_MIX
                  case (PERIODIC_MASTER)
                     grid%edge_type(edges(e),comp) = PERIODIC_SLAVE
                  end select
               endif

! if edge is a master and mate is not, then remove the masterness of edge and
! let mate's master be the master of the combined group

            elseif (is_periodic_edge_master(edges(e),grid,comp)) then

               select case (grid%edge_type(edges(e),comp))
               case (PERIODIC_MASTER_DIR)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_DIR
               case (PERIODIC_MASTER_NAT)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_NAT
               case (PERIODIC_MASTER_MIX)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE_MIX
               case (PERIODIC_MASTER)
                  grid%edge_type(edges(e),comp) = PERIODIC_SLAVE
               end select

! likewise if mate is a master and edge is not

            elseif (is_periodic_edge_master(edges_mate(e),grid,comp)) then

               select case (grid%edge_type(edges_mate(e),comp))
               case (PERIODIC_MASTER_DIR)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_DIR
               case (PERIODIC_MASTER_NAT)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_NAT
               case (PERIODIC_MASTER_MIX)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE_MIX
               case (PERIODIC_MASTER)
                  grid%edge_type(edges_mate(e),comp) = PERIODIC_SLAVE
               end select

! if edge and mate are both slaves, then leave the masters alone.  They should
! also show up as doubly periodic and get fixed then

            endif
         endif

      end do ! next component

! TEMP the following assumes double periodicity is the same in every
!      component that is periodic

! if neither is doubly periodic, point to each other

      if (.not. doubly_periodic1 .and. .not. doubly_periodic2) then
         grid%edge(edges(e))%next = edges_mate(e)
         grid%edge(edges_mate(e))%next = edges(e)

! if edge is doubly periodic and edge_mate is not, insert mate into edge's cycle

      elseif (doubly_periodic1 .and. .not. doubly_periodic2) then
         grid%edge(edges_mate(e))%next = grid%edge(edges(e))%next
         grid%edge(edges(e))%next = edges_mate(e)

! if edge_mate is doubly periodic and edge is not, insert edge into mate's cycle

      elseif (.not. doubly_periodic1 .and. doubly_periodic2) then
         grid%edge(edges(e))%next = grid%edge(edges_mate(e))%next
         grid%edge(edges_mate(e))%next = edges(e)

! If both are doubly periodic, merge the two cycles if they have not already
! been merged.  To do this, find the edge E that points to the mate.  If we
! see edge along the way, they have already been merged.  If we don't see edge,
! let E point to what edge points to and edge point to the mate.

      else
         last = grid%edge(edges_mate(e))%next
         do while (grid%edge(last)%next /= edges_mate(e))
            if (last == edges(e)) exit
            last = grid%edge(last)%next
         end do
         if (last /= edges(e)) then
            grid%edge(last)%next = grid%edge(edges(e))%next
            grid%edge(edges(e))%next = edges_mate(e)
         endif
      endif

      did_edge(edges(e),periodicity) = .true.
      did_edge(edges_mate(e),periodicity) = .true.
   end do ! next edge

! set up the vertices

   do v=1,VERTICES_PER_FACE
      if (did_vert(verts(v),periodicity)) cycle
      do comp=1,grid%system_size

! if this component of the face is not periodic, leave the vertex alone

         if (.not. is_periodic_face(face,grid,comp)) cycle

! determine if the vertex is doubly periodic and has already been assigned

         select case (grid%vertex_type(verts(v),comp))
         case(DIRICHLET,NATURAL,MIXED,PERIODIC)
            doubly_periodic1 = .false.
         case(PERIODIC_MASTER,PERIODIC_SLAVE,PERIODIC_MASTER_DIR, &
              PERIODIC_SLAVE_DIR,PERIODIC_MASTER_NAT,PERIODIC_SLAVE_NAT, &
              PERIODIC_MASTER_MIX,PERIODIC_SLAVE_MIX )
            doubly_periodic1 = .true.
         case default
            ierr = PHAML_INTERNAL_ERROR
            call fatal("setup_periodic: unrecognized vertex type")
            stop
         end select

! determine if the mate is doubly periodic and has already been assigned

         select case (grid%vertex_type(verts_mate(v),comp))
         case(DIRICHLET,NATURAL,MIXED,PERIODIC)
            doubly_periodic2 = .false.
         case(PERIODIC_MASTER,PERIODIC_SLAVE,PERIODIC_MASTER_DIR, &
              PERIODIC_SLAVE_DIR,PERIODIC_MASTER_NAT,PERIODIC_SLAVE_NAT, &
              PERIODIC_MASTER_MIX,PERIODIC_SLAVE_MIX )
            doubly_periodic2 = .true.
         case default
            ierr = PHAML_INTERNAL_ERROR
            call fatal("setup_periodic: unrecognized vertex type")
            stop
         end select

! if neither is doubly periodic, set the periodic type based on the type
! of the vertex and whether the face is the master

         if (.not. doubly_periodic1 .and. .not. doubly_periodic2) then

            if (grid%face_type(face,comp) == PERIODIC_MASTER) then
               select case (grid%vertex_type(verts(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
               select case (grid%vertex_type(verts_mate(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
            else
               select case (grid%vertex_type(verts(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
               select case (grid%vertex_type(verts_mate(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
            endif

! if vert is doubly periodic but mate is not, set mate's type and set the
! vert to a slave only if the mate is the master

         elseif (doubly_periodic1 .and. .not. doubly_periodic2) then

            if (grid%face_type(face,comp) == PERIODIC_MASTER) then
               select case (grid%vertex_type(verts_mate(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
            else
               select case (grid%vertex_type(verts_mate(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
               select case (grid%vertex_type(verts(v),comp))
               case(PERIODIC_MASTER_DIR, PERIODIC_SLAVE_DIR)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_DIR
               case(PERIODIC_MASTER_NAT, PERIODIC_SLAVE_NAT)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_NAT
               case(PERIODIC_MASTER_MIX, PERIODIC_SLAVE_MIX)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC_MASTER, PERIODIC_SLAVE)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
            endif

! if mate is doubly periodic but vert is not, set vert's type and set the
! type for mate to slave only if vert is the master

         elseif (.not. doubly_periodic1 .and. doubly_periodic2) then

            if (grid%face_type(face,comp) == PERIODIC_MASTER) then
               select case (grid%vertex_type(verts(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER_DIR
               case(NATURAL)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER_NAT
               case(MIXED)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER_MIX
               case(PERIODIC)
                  grid%vertex_type(verts(v),comp) = PERIODIC_MASTER
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
               select case (grid%vertex_type(verts_mate(v),comp))
               case(PERIODIC_MASTER_DIR, PERIODIC_SLAVE_DIR)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_DIR
               case(PERIODIC_MASTER_NAT, PERIODIC_SLAVE_NAT)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_NAT
               case(PERIODIC_MASTER_MIX, PERIODIC_SLAVE_MIX)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC_MASTER, PERIODIC_SLAVE)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
            else
               select case (grid%vertex_type(verts(v),comp))
               case(DIRICHLET)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_DIR
               case(NATURAL)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_NAT
               case(MIXED)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_MIX
               case(PERIODIC)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE
               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("setup_periodic: unrecognized vertex type")
                  stop
               end select
            endif

! if both are doubly periodic ...

         else

! if they are both masters, remove the masterness of one of them depending on
! the masterness of face

            if (is_periodic_vert_master(verts(v),grid,comp) .and. &
                is_periodic_vert_master(verts_mate(v),grid,comp)) then

               if (grid%face_type(face,comp) == PERIODIC_MASTER) then
                  select case (grid%vertex_type(verts_mate(v),comp))
                  case (PERIODIC_MASTER_DIR)
                     grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_DIR
                  case (PERIODIC_MASTER_NAT)
                     grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_NAT
                  case (PERIODIC_MASTER_MIX)
                     grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_MIX
                  case (PERIODIC_MASTER)
                     grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE
                  end select
               else
                  select case (grid%vertex_type(verts(v),comp))
                  case (PERIODIC_MASTER_DIR)
                     grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_DIR
                  case (PERIODIC_MASTER_NAT)
                     grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_NAT
                  case (PERIODIC_MASTER_MIX)
                     grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_MIX
                  case (PERIODIC_MASTER)
                     grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE
                  end select
               endif

! if vert is a master and mate is not, then remove the masterness of vert and
! let mate's master be the master of the combined group

            elseif (is_periodic_vert_master(verts(v),grid,comp)) then

               select case (grid%vertex_type(verts(v),comp))
               case (PERIODIC_MASTER_DIR)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_DIR
               case (PERIODIC_MASTER_NAT)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_NAT
               case (PERIODIC_MASTER_MIX)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE_MIX
               case (PERIODIC_MASTER)
                  grid%vertex_type(verts(v),comp) = PERIODIC_SLAVE
               end select

! likewise if mate is a master and vert is not

            elseif (is_periodic_vert_master(verts_mate(v),grid,comp)) then

               select case (grid%vertex_type(verts_mate(v),comp))
               case (PERIODIC_MASTER_DIR)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_DIR
               case (PERIODIC_MASTER_NAT)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_NAT
               case (PERIODIC_MASTER_MIX)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE_MIX
               case (PERIODIC_MASTER)
                  grid%vertex_type(verts_mate(v),comp) = PERIODIC_SLAVE
               end select

! if vert and mate are both slaves, then leave the masters alone.  They should
! also show up as doubly periodic and get fixed then

            endif
         endif

      end do ! next component

! TEMP the following assumes double periodicity is the same in every
!      component that is periodic

! if neither is doubly periodic, point to each other

      if (.not. doubly_periodic1 .and. .not. doubly_periodic2) then
         grid%vertex(verts(v))%next = verts_mate(v)
         grid%vertex(verts_mate(v))%next = verts(v)

! if vert is doubly periodic and vert_mate is not, insert mate into vert's cycle

      elseif (doubly_periodic1 .and. .not. doubly_periodic2) then
         grid%vertex(verts_mate(v))%next = grid%vertex(verts(v))%next
         grid%vertex(verts(v))%next = verts_mate(v)

! if vert_mate is doubly periodic and vert is not, insert vert into mate's cycle

      elseif (.not. doubly_periodic1 .and. doubly_periodic2) then
         grid%vertex(verts(v))%next = grid%vertex(verts_mate(v))%next
         grid%vertex(verts_mate(v))%next = verts(v)

! If both are doubly periodic, merge the two cycles if they have not already
! been merged.  To do this, find the vertex V that points to the mate.  If we
! see vert along the way, they have already been merged.  If we don't see vert,
! let V point to what vert points to and vert point to the mate.

      else
         last = grid%vertex(verts_mate(v))%next
         do while (grid%vertex(last)%next /= verts_mate(v))
            if (last == verts(v)) exit
            last = grid%vertex(last)%next
         end do
         if (last /= verts(v)) then
            grid%vertex(last)%next = grid%vertex(verts(v))%next
            grid%vertex(verts(v))%next = verts_mate(v)
         endif
      endif

      did_vert(verts(v),periodicity) = .true.
      did_vert(verts_mate(v),periodicity) = .true.
   end do ! next vertex

! remove the face and mate from the linked list

   call remove_from_linked_list(face,first_face,last_face,next_face,prev_face)
   call remove_from_linked_list(mate,first_face,last_face,next_face,prev_face)

end do ! next periodic face

! Set initial_neighbor for elements on periodic faces.  The neighbor element
! must be the associated element of the matching periodic face.

do e=1,grid%nelem
   do f=1,FACES_PER_ELEMENT
      face = grid%element(e)%face(f)
      if (is_periodic_face(face,grid)) then
         grid%initial_neighbor(f,e) = grid%face(grid%face(face)%next)%assoc_elem
      endif
   end do
end do

contains

function myeq(x,y)
real(my_real), intent(in) :: x,y
logical :: myeq
! TEMP120620 use single precision epsilon as small
!real(my_real), parameter :: small = 100*epsilon(0.0_my_real)
real(my_real), parameter :: small = epsilon(0.0)
real(my_real) :: normalize
normalize = max(abs(x),abs(y))
if (normalize == 0.0_my_real) then
   myeq = .true.
else
   myeq = abs(x-y)/normalize < small
endif
end function myeq

end subroutine setup_periodic

!          ------------------
subroutine add_to_linked_list(addme,first,last,next,prev)
!          ------------------

!----------------------------------------------------
! This routine adds addme to the linked list defined by the other arguments
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: addme
integer, intent(inout) :: first, last, next(:), prev(:)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (first == -1) then
   first = addme
   last = addme
   next(addme) = -1
   prev(addme) = -1
else
   next(last) = addme
   prev(addme) = last
   next(addme) = -1
   last = addme
endif

end subroutine add_to_linked_list

!          -----------------------
subroutine remove_from_linked_list(removeme,first,last,next,prev)
!          -----------------------

!----------------------------------------------------
! This routine removes removeme from the linked list defined by the other args
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: removeme
integer, intent(inout) :: first, last, next(:), prev(:)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (prev(removeme) == -1) then
   first = next(removeme)
else
   next(prev(removeme)) = next(removeme)
endif
if (next(removeme) == -1) then
   last = prev(removeme)
else
   prev(next(removeme)) = prev(removeme)
endif

end subroutine remove_from_linked_list

!          ------------------------
subroutine setup_polygonal_boundary(grid)
!          ------------------------

!----------------------------------------------------
! This routine sets the boundary data structures for a polygonal boundary.
! Each boundary face is defined as a surface.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: face, ivert, vert, iedge, edge, n, astat
integer :: vert_to_point(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

! initial allocation of data structures

allocate(grid%point(8), grid%line(8), grid%lineloop(8), grid%surface(8), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in setup_polygonal_boundary")
   stop
endif
grid%npoint = 0
grid%nline = 0
grid%nlineloop = 0
grid%nsurface = 0

! Make each boundary face be a surface.  Create a point and line for each
! vertex and edge of the surface; don't worry about identifying common
! points and lines, I'm pretty sure I don't look to see if surfaces have
! common points or lines

do face=1,grid%nface
   if (grid%face(face)%bmark == 0) cycle ! interior face

! create a point for each vertex

   do ivert=1,VERTICES_PER_FACE
      vert = grid%face(face)%vertex(ivert)
      grid%npoint = grid%npoint + 1
      if (grid%npoint > size(grid%point)) call more_points(grid)
      grid%point(grid%npoint)%coord = grid%vertex(vert)%coord
      vert_to_point(vert) = grid%npoint
   end do

! create a line for each edge

   do iedge = 1,EDGES_PER_FACE
      edge = grid%face(face)%edge(iedge)
      grid%nline = grid%nline + 1
      if (grid%nline > size(grid%line)) call more_lines(grid)
      grid%line(grid%nline)%type = LINE_LINE
      allocate(grid%line(grid%nline)%point(2),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in setup_polygonal_boundary")
         stop
      endif
      grid%line(grid%nline)%point(1) = vert_to_point(grid%edge(edge)%vertex(1))
      grid%line(grid%nline)%point(2) = vert_to_point(grid%edge(edge)%vertex(2))
   end do

! create a line loop; begin with grid%nline-2 and add grid%nline-1 and
! grid%nline in correct order and orientation to match endpoints

   grid%nlineloop = grid%nlineloop + 1
   if (grid%nlineloop > size(grid%lineloop)) call more_lineloops(grid)
   allocate(grid%lineloop(grid%nlineloop)%line(3),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in setup_polygonal_boundary")
      stop
   endif
   n = grid%nline
   grid%lineloop(grid%nlineloop)%line(1) = n-2
   if (grid%line(n-1)%point(1) == grid%line(n-2)%point(2)) then
      grid%lineloop(grid%nlineloop)%line(2) = n-1
      if (grid%line(n)%point(1) == grid%line(n-1)%point(2)) then
         grid%lineloop(grid%nlineloop)%line(3) = n
      else
         grid%lineloop(grid%nlineloop)%line(3) = -n
      endif
   elseif (grid%line(n-1)%point(2) == grid%line(n-2)%point(2)) then
      grid%lineloop(grid%nlineloop)%line(2) = -(n-1)
      if (grid%line(n)%point(1) == grid%line(n-1)%point(1)) then
         grid%lineloop(grid%nlineloop)%line(3) = n
      else
         grid%lineloop(grid%nlineloop)%line(3) = -n
      endif
   elseif (grid%line(n)%point(1) == grid%line(n-2)%point(2)) then
      grid%lineloop(grid%nlineloop)%line(2) = n
      if (grid%line(n-1)%point(1) == grid%line(n)%point(2)) then
         grid%lineloop(grid%nlineloop)%line(3) = n-1
      else
         grid%lineloop(grid%nlineloop)%line(3) = -(n-1)
      endif
   else
      grid%lineloop(grid%nlineloop)%line(2) = -n
      if (grid%line(n-1)%point(1) == grid%line(n)%point(1)) then
         grid%lineloop(grid%nlineloop)%line(3) = n-1
      else
         grid%lineloop(grid%nlineloop)%line(3) = -(n-1)
      endif
   endif

! create the surface

   grid%nsurface = grid%nsurface + 1
   if (grid%nsurface > size(grid%surface)) call more_surfaces(grid)
   allocate(grid%surface(grid%nsurface)%lineloop(1),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in setup_polygonal_boundary")
      stop
   endif
   grid%surface(grid%nsurface)%type = SURF_PLANAR
   grid%surface(grid%nsurface)%lineloop(1) = grid%nlineloop

end do

end subroutine setup_polygonal_boundary

!          -------------------
subroutine set_boundary_vertex(grid)
!          -------------------

!----------------------------------------------------
! This routine creates the linked list of bverts for each boundary vertex
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: vert, ivert, elem, surf
logical :: success, visited_vert(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

! for each boundary vertex ...

visited_vert = .false.
do elem=1,grid%nelem
   do ivert=1,VERTICES_PER_ELEMENT
      vert = grid%element(elem)%vertex(ivert)
      if (visited_vert(vert)) cycle
      visited_vert(vert) = .true.
      if (grid%vertex(vert)%bmark == 0) cycle ! interior vertex

! check each surface to see if this vertex is on it

      do surf=1,grid%nsurface

! See if it is one of the points (corners) of the surface.  If so, the
! bvert is defined and move on to the next surface.

         call try_vert_as_point(grid,vert,surf,success)
         if (success) cycle

! See if it is on one of the lines.  If so, the bvert is defined and move on
! to the next surface.

         call try_vert_on_line(grid,vert,surf,success)
         if (success) cycle

! See if it is on the interior of the surface.  If so, bvert is defined.

         call try_vert_on_surface(grid,vert,surf)

      end do
   end do
end do

end subroutine set_boundary_vertex

!          -----------------
subroutine try_vert_as_point(grid,vert,surf,success)
!          -----------------

!----------------------------------------------------
! This routine looks to see if vert is one of the points (corners) of surf.
! If so, it defines the next boundary_vertex of vert and returns success=true.
! Otherwise it returns success=false.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: vert, surf
logical, intent(out) :: success
!----------------------------------------------------
! Local variables:

integer :: ilineloop, iline, line, point, astat
type(bvert_type), pointer :: bvert
!----------------------------------------------------
! Begin executable code

success = .false.

! for each point (as the first endpoint of a line)

do ilineloop=1,size(grid%surface(surf)%lineloop)
   do iline=1,size(grid%lineloop(grid%surface(surf)%lineloop(ilineloop))%line)
      line = grid%lineloop(grid%surface(surf)%lineloop(ilineloop))%line(iline)
      if (line > 0) then
         point = grid%line(line)%point(1)
      else
         point = grid%line(-line)%point(size(grid%line(-line)%point))
      endif

! if this point is the vertex, define a new bvert

      if (vert_is_point(grid,vert,point)) then

! bvert points to a new bvert

         if (associated(grid%vertex(vert)%boundary_vertex)) then
            bvert => grid%vertex(vert)%boundary_vertex
            do while (associated(bvert%next))
               bvert => bvert%next
            end do
            allocate(bvert%next,stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("allocation failed in try_vert_as_point")
               stop
            endif
            bvert => bvert%next
         else
            allocate(grid%vertex(vert)%boundary_vertex,stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("allocation failed in try_vert_as_point")
               stop
            endif
            bvert => grid%vertex(vert)%boundary_vertex
         endif

! define the bvert

         bvert%surface = surf
         nullify(bvert%next)

! the point is the beginning of the line, or the end if the line ID is negative

         if (line > 0) then
            bvert%lines(1) = line
            bvert%p1 = 0.0_my_real
         else
            bvert%lines(1) = -line
            bvert%p1 = 1.0_my_real
         endif

! the point should also be the end of the preceding line in the line loop

         if (iline == 1) then
            line = grid%lineloop(grid%surface(surf)%lineloop(ilineloop))%line( &
               size(grid%lineloop(grid%surface(surf)%lineloop(ilineloop))%line))
         else
            line = &
             grid%lineloop(grid%surface(surf)%lineloop(ilineloop))%line(iline-1)
         endif

         if (line > 0) then
            if (grid%line(line)%point(size(grid%line(line)%point)) /= point) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("try_vert_as_point: mismatch of lineloop endpoints")
               stop
            endif
            bvert%lines(2) = line
            bvert%p2 = 1.0_my_real
         else
            if (grid%line(-line)%point(1) /= point) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("try_vert_as_point: mismatch of lineloop endpoints")
               stop
            endif
            bvert%lines(2) = -line
            bvert%p2 = 0.0_my_real
         endif

! above, (p1,p2) was set for a planar surface.  reset for ruled surfaces

         if (grid%surface(surf)%type == SURF_RULED_3) then

! barycentric coordinates are (p1,p2,1-p1-p2)

            select case(iline)
            case (1)
               bvert%p1 = 1.0_my_real
               bvert%p2 = 0.0_my_real
            case (2)
               bvert%p1 = 0.0_my_real
               bvert%p2 = 1.0_my_real
            case (3)
               bvert%p1 = 0.0_my_real
               bvert%p2 = 0.0_my_real
            end select

         elseif (grid%surface(surf)%type == SURF_RULED_4) then

! (xi,eta) for a corner

            select case(iline)
            case (1)
               bvert%p1 = 0.0_my_real
               bvert%p2 = 0.0_my_real
            case (2)
               bvert%p1 = 1.0_my_real
               bvert%p2 = 0.0_my_real
            case (3)
               bvert%p1 = 1.0_my_real
               bvert%p2 = 1.0_my_real
            case (4)
               bvert%p1 = 0.0_my_real
               bvert%p2 = 1.0_my_real
            end select

         endif

! no need to look at any more lines

         success = .true.
         exit
      end if
   end do
end do

end subroutine try_vert_as_point

!        -------------
function vert_is_point(grid,vert,point)
!        -------------

!----------------------------------------------------
! This routine returns true if the vertex and point are close enough to
! believe they are the same
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: vert, point
logical :: vert_is_point
!----------------------------------------------------
! Local variables:

real(my_real) :: norm_vert, norm_point, norm_diff, maxnorm
!----------------------------------------------------
! Begin executable code

norm_vert = sqrt(dot_point(grid%vertex(vert)%coord,grid%vertex(vert)%coord))
norm_point = sqrt(dot_point(grid%point(point)%coord,grid%point(point)%coord))
norm_diff = sqrt(dot_point(grid%vertex(vert)%coord-grid%point(point)%coord, &
                           grid%vertex(vert)%coord-grid%point(point)%coord))
maxnorm = max(norm_vert,norm_point)
if (maxnorm < really_small) then
   vert_is_point = .true.
elseif (norm_diff/maxnorm < small) then
   vert_is_point = .true.
else
   vert_is_point = .false.
endif

end function vert_is_point

!          ----------------
subroutine try_vert_on_line(grid,vert,surf,success)
!          ----------------

!----------------------------------------------------
! This routine looks to see if vert is on one of the lines of surf.
! If so, it defines the next boundary_vertex of vert and returns success=true.
! Otherwise it returns success=false.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: vert, surf
logical, intent(out) :: success
!----------------------------------------------------
! Local variables:

integer :: ilineloop, iline, line, sline, astat
real(my_real) :: param
type(bvert_type), pointer :: bvert
!----------------------------------------------------
! Begin executable code

success = .false.

! for each line

do ilineloop=1,size(grid%surface(surf)%lineloop)
   do iline=1,size(grid%lineloop(grid%surface(surf)%lineloop(ilineloop))%line)
      sline = grid%lineloop(grid%surface(surf)%lineloop(ilineloop))%line(iline)
      line = abs(sline)

! if this point is the vertex, define a new bvert

      if (vert_is_on_line(grid,vert,line,param)) then

! bvert points to a new bvert

         if (associated(grid%vertex(vert)%boundary_vertex)) then
            bvert => grid%vertex(vert)%boundary_vertex
            do while (associated(bvert%next))
               bvert => bvert%next
            end do
            allocate(bvert%next,stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("allocation failed in try_vert_as_point")
               stop
            endif
            bvert => bvert%next
         else
            allocate(grid%vertex(vert)%boundary_vertex,stat=astat)
            bvert => grid%vertex(vert)%boundary_vertex
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("allocation failed in try_vert_as_point")
               stop
            endif
         endif

! define the bvert

         bvert%surface = surf
         nullify(bvert%next)
         bvert%lines(1) = line
         bvert%lines(2) = 0
         select case (grid%surface(surf)%type)
         case (SURF_PLANAR)
            bvert%p1 = param
            bvert%p2 = 0.0_my_real
         case (SURF_RULED_3)
            if (sline < 0) param = 1.0_my_real - param
            select case (iline)
            case (1)
               bvert%p1 = 1.0_my_real - param
               bvert%p2 = param
            case (2)
               bvert%p1 = 0.0_my_real
               bvert%p2 = 1.0_my_real - param
            case (3)
               bvert%p1 = param
               bvert%p2 = 0.0_my_real
            end select
         case (SURF_RULED_4)
            if (sline < 0) param = 1.0_my_real - param
            select case (iline)
            case (1)
               bvert%p1 = param
               bvert%p2 = 0.0_my_real
            case (2)
               bvert%p1 = 1.0_my_real
               bvert%p2 = param
            case (3)
               bvert%p1 = 1.0_my_real - param
               bvert%p2 = 1.0_my_real
            case (4)
               bvert%p1 = 0.0_my_real
               bvert%p2 = 1.0_my_real - param
            end select
         end select

! no need to look at any more lines

         success = .true.
         exit
      end if
   end do
end do

end subroutine try_vert_on_line

!        ---------------
function vert_is_on_line(grid,vert,line,param)
!        ---------------

!----------------------------------------------------
! This routine returns true if vert is close enough to line to believe
! it is on it.  If so, it returns the parameter for how far along the
! line it is.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: vert, line
real(my_real), intent(out) :: param
logical :: vert_is_on_line
!----------------------------------------------------
! Local variables:

real(my_real) :: xdiff_line, ydiff_line, zdiff_line, xdiff_vert, ydiff_vert, &
                 zdiff_vert, size, frac, r1, rv
!----------------------------------------------------
! Begin executable code

vert_is_on_line = .false.
param = 0.0_my_real

select case (grid%line(line)%type)

! for a straight line, the point should be a linear combination of the
! endpoints, i.e., the fraction for x, y and z should be the same.
! Also, it should be between 0 and 1.

case (LINE_LINE)

   xdiff_line = grid%point(grid%line(line)%point(2))%coord%x - &
                grid%point(grid%line(line)%point(1))%coord%x
   ydiff_line = grid%point(grid%line(line)%point(2))%coord%y - &
                grid%point(grid%line(line)%point(1))%coord%y
   zdiff_line = grid%point(grid%line(line)%point(2))%coord%z - &
                grid%point(grid%line(line)%point(1))%coord%z
   xdiff_vert = grid%vertex(vert)%coord%x - &
                grid%point(grid%line(line)%point(1))%coord%x
   ydiff_vert = grid%vertex(vert)%coord%y - &
                grid%point(grid%line(line)%point(1))%coord%y
   zdiff_vert = grid%vertex(vert)%coord%z - &
                grid%point(grid%line(line)%point(1))%coord%z

   size = max(abs(xdiff_line),max(abs(ydiff_line),abs(zdiff_line)))

   if (abs(xdiff_line) == size) then
      frac = xdiff_vert/xdiff_line
   elseif (abs(ydiff_line) == size) then
      frac = ydiff_vert/ydiff_line
   else
      frac = zdiff_vert/zdiff_line
   endif

   vert_is_on_line = .true.

   if (frac < -small*size .or. frac > 1.0_my_real + small*size) then
      vert_is_on_line = .false.
   endif

   if (vert_is_on_line) then
      if (abs(xdiff_line) < small*size) then
         if (abs(xdiff_vert) > small*size) then
            vert_is_on_line = .false.
         endif
      elseif (abs(xdiff_vert/xdiff_line - frac) > small*size) then
         vert_is_on_line = .false.
      endif
   endif

   if (vert_is_on_line) then
      if (abs(ydiff_line) < small*size) then
         if (abs(ydiff_vert) > small*size) then
            vert_is_on_line = .false.
         endif
      elseif (abs(ydiff_vert/ydiff_line - frac) > small*size) then
         vert_is_on_line = .false.
      endif
   endif

   if (vert_is_on_line) then
      if (abs(zdiff_line) < small*size) then
         if (abs(zdiff_vert) > small*size) then
            vert_is_on_line = .false.
         endif
      elseif (abs(zdiff_vert/zdiff_line - frac) > small*size) then
         vert_is_on_line = .false.
      endif
   endif

   if (vert_is_on_line) param = frac

! for a circular line, the point must be in the plane of the circle, must be
! the right distance from the center, and must have an angle between the
! angles of the two endpoints

case (LINE_CIRCLE)

   vert_is_on_line = coord_is_in_plane(grid%vertex(vert)%coord, &
                                  grid%point(grid%line(line)%point(1))%coord, &
                                  grid%point(grid%line(line)%point(2))%coord, &
                                  grid%point(grid%line(line)%point(3))%coord)

   if (vert_is_on_line) then
      r1 = sqrt(dot_point(grid%point(grid%line(line)%point(1))%coord - &
                          grid%point(grid%line(line)%point(2))%coord, &
                          grid%point(grid%line(line)%point(1))%coord - &
                          grid%point(grid%line(line)%point(2))%coord))
      rv = sqrt(dot_point(grid%vertex(vert)%coord - &
                          grid%point(grid%line(line)%point(2))%coord, &
                          grid%vertex(vert)%coord - &
                          grid%point(grid%line(line)%point(2))%coord))
      if (abs(r1-rv) > small*max(abs(r1),abs(rv))) vert_is_on_line = .false.
   endif

   if (vert_is_on_line) then
      frac = circle_relative_angle(grid%vertex(vert)%coord,grid,line)
      if (frac < -small .or. frac > 1.0_my_real + small) then
         vert_is_on_line = .false.
      endif
   endif

   if (vert_is_on_line) param = frac

! for an elliptical line, the point must be in the plane of the ellipse,
! TEMP and I don't know what else

case (LINE_ELLIPSE)

   vert_is_on_line = coord_is_in_plane(grid%vertex(vert)%coord, &
                                  grid%point(grid%line(line)%point(1))%coord, &
                                  grid%point(grid%line(line)%point(2))%coord, &
                                  grid%point(grid%line(line)%point(4))%coord)

! TEMP120309
   if (vert_is_on_line) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("haven't written code to determine if a vertex is on an ellipse")
      stop
   endif

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("vert_is_on_line: unrecognized line type")
   stop

end select

end function vert_is_on_line

!        -----------------
function coord_is_in_plane(p,p1,p2,p3)
!        -----------------

!----------------------------------------------------
! This routine returns true if the point p is in the plane defined by p1,p2,p3
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: p, p1, p2, p3
logical :: coord_is_in_plane
!----------------------------------------------------
! Local variables:

real(my_real) :: a, b, c, d, check, scale
!----------------------------------------------------
! Begin executable code

! equation of plane a*x + b*y + c*z + d = 0

a = p1%y*(p2%z-p3%z) + p2%y*(p3%z-p1%z) + p3%y*(p1%z-p2%z)
b = p1%z*(p2%x-p3%x) + p2%z*(p3%x-p1%x) + p3%z*(p1%x-p2%x)
c = p1%x*(p2%y-p3%y) + p2%x*(p3%y-p1%y) + p3%x*(p1%y-p2%y)
d = -(p1%x*(p2%y*p3%z-p3%y*p2%z) + p2%x*(p3%y*p1%z-p1%y*p3%z) + &
      p3%x*(p1%y*p2%z-p2%y*p1%z))

scale = max(abs(a),max(abs(b),max(abs(c),abs(d))))

! check p for satisfying the equation

check = (a*p%x + b*p%y + c*p%z + d)/scale

coord_is_in_plane = (abs(check) < small)

end function coord_is_in_plane

!        ---------------------
function circle_relative_angle(coord,grid,line)
!        ---------------------

!----------------------------------------------------
! This routine returns the ratio of the angle between coord and the first
! endpoint of the circular line to the angle between the second endpoint and
! the first.  coord must lie on the circle defined by line (this is not
! checked), but does not have to be on line.  The result is between 0 and 1 iff
! it is on the line.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: coord
type(grid_type), intent(in) :: grid
integer, intent(in) :: line
real(my_real) :: circle_relative_angle
!----------------------------------------------------
! Local variables:

type(point) :: P1, P2, C
real(my_real) :: r, theta1, theta2, theta3, frac
!----------------------------------------------------
! Begin executable code

if (grid%line(line)%type /= LINE_CIRCLE) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("called circle_relative_angle with a line that is not a circle")
   stop
endif

P1 = grid%point(grid%line(line)%point(1))%coord
P2 = grid%point(grid%line(line)%point(3))%coord
C  = grid%point(grid%line(line)%point(2))%coord

! radius

r = sqrt(dot_point(P1-C,P1-C))
if (r == 0.0_my_real) then
   ierr = USER_INPUT_ERROR
   call fatal("circle_relative_angle: looks like a point on a circle is the center of the circle")
   stop
endif

! angle between P2-C and P1-C

theta1 = acos(dot_point(P1-C,P2-C)/r**2)
if (theta1 == 0.0_my_real) then
   ierr = USER_INPUT_ERROR
   call fatal("looks like the two endpoints of a circular arc are the same point")
   stop
endif

! angle between coord-C and P1-C

theta2 = acos(dot_point(P1-C,coord-C)/r**2)

! ratio of the angle

frac = theta2/theta1

! If frac is greater than one, coord is not between P1 and P2, and we can quit.
! frac is positive, because acos is between 0 and pi, but is should be
! negative if coord lies outside P1 compared to P2.  This would imply that the
! angle between coord-C and P2-C is bigger than theta1.

if (frac <= 1.0_my_real) then
   theta3 = acos(dot_point(P2-C,coord-C)/r**2)
   if (theta3 > theta1) frac = -frac
endif

circle_relative_angle = frac

end function circle_relative_angle

!          -------------------
subroutine try_vert_on_surface(grid,vert,surf)
!          -------------------

!----------------------------------------------------
! This routine looks to see if vert is on the surface.  If so, it defines the
! next boundary_vertex of vert.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: vert, surf
!----------------------------------------------------
! Local variables:

integer :: l1, l2, ip1, ip2, ip3, i, astat, ilineloop, lineloop, iline, line, &
           crossings
type(point) :: p1, p2, p3
type(bvert_type), pointer :: bvert
logical :: success, v_same_sign
real(my_real) :: param1, param2, resid
!----------------------------------------------------
! Begin executable code

select case(grid%surface(surf)%type)

! for a planar surface, see if the point is in the plane defined by
! three vertices of the line loop

case (SURF_PLANAR)

! find points to define the plane of the surface

   l1 = grid%lineloop(grid%surface(surf)%lineloop(1))%line(1)
   select case (grid%line(abs(l1))%type)

   case (LINE_LINE)

! for a line, the points are p1) the first endpoint of the first line segment,
! p2) the second endpoint of the first line segment, p3) the second endpoint
! of the second line segment

! TEMP this will fail if the three points are colinear

      if (l1 > 0) then
         ip1 = grid%line(abs(l1))%point(1)
         ip2 = grid%line(abs(l1))%point(2)
      else
         ip1 = grid%line(abs(l1))%point(2)
         ip2 = grid%line(abs(l1))%point(1)
      endif
      l2 = grid%lineloop(grid%surface(surf)%lineloop(1))%line(2)
      if (l2 > 0) then
         select case (grid%line(abs(l2))%type)
         case (LINE_LINE)
            ip3 = grid%line(abs(l2))%point(2)
         case (LINE_CIRCLE)
            ip3 = grid%line(abs(l2))%point(3)
         case (LINE_ELLIPSE)
            ip3 = grid%line(abs(l2))%point(4)
         end select
      else
         ip3 = grid%line(abs(l2))%point(1)
      endif
      p1 = grid%point(ip1)%coord
      p2 = grid%point(ip2)%coord
      p3 = grid%point(ip3)%coord

   case (LINE_CIRCLE, LINE_ELLIPSE)

! for a circle or ellipse, the points are p1) the first endpoint of the arc,
! p2) the center of the circle/ellipse, p3) the second endpoint of the arc

      p1 = grid%point(grid%line(abs(grid%lineloop(grid%surface(surf)%lineloop(1))%line(1)))%point(1))%coord
      p2 = grid%point(grid%line(abs(grid%lineloop(grid%surface(surf)%lineloop(1))%line(1)))%point(2))%coord
      if (grid%line(abs(l1))%type == LINE_CIRCLE) then
         p3 = grid%point(grid%line(abs(grid%lineloop(grid%surface(surf)%lineloop(1))%line(1)))%point(3))%coord
      else
         p3 = grid%point(grid%line(abs(grid%lineloop(grid%surface(surf)%lineloop(1))%line(1)))%point(4))%coord
      endif
   end select

! see if this point is in the plane of the surface

   success = coord_is_in_plane(grid%vertex(vert)%coord,p1,p2,p3)
   if (success) then

! In the plane, see if it is on the surface.

! Let vert be the origin of a coordinate system for the surface with
! p1-p2 and p3-p2 orthogonalized to p1-p2 as the axes, refered to as the u and
! v directions.  We are looking to see if vert is inside a polygon with edges
! that are line segments, circular arcs and elliptical arcs.
! The crossing test determines inside/outside: count the number of
! intersections between a ray from vert with the polygon.  If it is odd, the
! point is inside; if it is even, the point is outside.  Take the direction of
! the ray to be +u.

      crossings = 0

! for each edge of the surface ...

      do ilineloop=1,size(grid%surface(surf)%lineloop)
         lineloop = grid%surface(surf)%lineloop(ilineloop)
         do iline=1,size(grid%lineloop(lineloop)%line)
            line = abs(grid%lineloop(lineloop)%line(iline))

! increment crossings if the line intersects the positive u axis

            select case(grid%line(line)%type)
            case (LINE_LINE)
               crossings = crossings + line_crossing(grid,line,vert,p1,p2,p3)
            case (LINE_CIRCLE)
               crossings = crossings + circle_crossing(grid,line,vert,p1,p2,p3)
            case (LINE_ELLIPSE)
               crossings = crossings + ellipse_crossing(grid,line,vert,p1,p2,p3)
            end select

         end do
      end do

! check evenness of crossings

      success = crossings /= 2*(crossings/2)

      param1 = 0.0_my_real
      param2 = 0.0_my_real

   endif ! point is in the plane of the surface

case (SURF_RULED_3, SURF_RULED_4)

   success = coord_is_on_surface(grid%vertex(vert)%coord,surf,grid, &
                                 param1,param2,resid)

   if (success) then
      if (param1 < -small .or. param1 > 1.0_my_real+small .or. &
          param2 < -small .or. param2 > 1.0_my_real+small) success = .false.
      if (grid%surface(surf)%type == SURF_RULED_3 .and. &
          param1 + param2 > 1.0_my_real + small) success = .false.
   endif

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("try_vert_on_surface: unrecognized surface type")
   stop

end select

! if the point is on the surface ...

if (success) then

! bvert points to a new bvert

   if (associated(grid%vertex(vert)%boundary_vertex)) then
      bvert => grid%vertex(vert)%boundary_vertex
      do while (associated(bvert%next))
         bvert => bvert%next
      end do
      allocate(bvert%next,stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in try_vert_on_surface")
         stop
      endif
      bvert => bvert%next
   else
      allocate(grid%vertex(vert)%boundary_vertex,stat=astat)
      bvert => grid%vertex(vert)%boundary_vertex
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in try_vert_on_surface")
         stop
      endif
   endif

! define the bvert

   bvert%surface = surf
   nullify(bvert%next)
   bvert%lines(1) = 0
   bvert%lines(2) = 0
   bvert%p1 = param1
   bvert%p2 = param2

endif

end subroutine try_vert_on_surface

!        -------------
function line_crossing(grid,line,vert,p1,p2,p3)
!        -------------

!----------------------------------------------------
! This routine returns 1 if an edge that is a line segment crosses the positive
! u axis.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: line, vert
type(point), intent(in) :: p1, p2, p3
integer :: line_crossing
!----------------------------------------------------
! Local variables:

type(point_uv) :: end1, end2
!----------------------------------------------------
! Begin executable code

! If the v coordinates have the same sign, there is no intersection.  If they
! don't have the same sign, compute the intersection of the edge and u axis.
! There is an itersection iff the u coordinate is positive.
! If the v coordinate of a vertex is (nearly) exactly 0, then the ray is
! passing through the vertex.  vert should not be a vertex or on a line since
! those were checked before calling try_vert_on_surface.  A 0 v coordinate is
! increased a tiny amount so that the ray only intersects one of the edges.
! If both v coordinates of a line edge are 0, then the ray travels along the
! edge.  The tiny increase will cause the ray to miss the edge, but since there
! are two vertices not counted, it will not change the even-ness of the count.
! If the v coordinates have different signs and the u coordinate is 0, then the
! vert is on the edge, which should not happen.

! get the u-v coordinates of the endpoints of the line

end1 = to_uv(grid%point(grid%line(line)%point(1))%coord, &
             grid%vertex(vert)%coord,p1,p3,p2)
end2 = to_uv(grid%point(grid%line(line)%point(2))%coord, &
             grid%vertex(vert)%coord,p1,p3,p2)

! fudge 0 v coordinates

if (end1%v < small .and. end1%v > -small) end1%v = end1%v + 2*small
if (end2%v < small .and. end2%v > -small) end2%v = end2%v + 2*small

line_crossing = 0

! if the v coordinates have different signs

if ((end1%v > 0.0_my_real .and. end2%v < 0.0_my_real) .or. &
    (end1%v < 0.0_my_real .and. end2%v > 0.0_my_real)) then

! and the line intersects the positive u axis

   if (abs(end1%v-end2%v) > small) then
      if (end1%u - ((end2%u-end1%u)/(end2%v-end1%v))*end1%v > 0.0_my_real) then

! then we have a crossing

         line_crossing = 1

      endif
   endif
endif

end function line_crossing

!        ---------------
function circle_crossing(grid,line,vert,p1,p2,p3)
!        ---------------

!----------------------------------------------------
! This routine returns 1 if an edge that is a circle arc crosses the positive
! u axis once, and 0 if it doesn't cross or crosses twice.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: line, vert
type(point), intent(in) :: p1, p2, p3
integer :: circle_crossing
!----------------------------------------------------
! Local variables:

type(point_uv) :: end1, end2, center, uplus_c, uminus_c, end1_c, end2_c
real(my_real) :: r, discrim, uplus, uminus, norm_uplus, norm_uminus, &
                 norm_end1, norm_end2, angle_uplus_end1, angle_uplus_end2, &
                 angle_uminus_end1, angle_uminus_end2, angle_end1_end2
logical :: uplus_in_arc, uminus_in_arc
!----------------------------------------------------
! Begin executable code

! get the u-v coordinates of the endpoints of the arc and center of the circle

end1 = to_uv(grid%point(grid%line(line)%point(1))%coord, &
             grid%vertex(vert)%coord,p1,p3,p2)
end2 = to_uv(grid%point(grid%line(line)%point(3))%coord, &
             grid%vertex(vert)%coord,p1,p3,p2)
center = to_uv(grid%point(grid%line(line)%point(2))%coord, &
             grid%vertex(vert)%coord,p1,p3,p2)

! r is the radius of the circle

r = sqrt((end1%u-center%u)**2 + (end1%v-center%v)**2)

! fudge 0 v coordinates, and keep the point on the circle

if (end1%v < small*r .and. end1%v > -small*r) then
   end1%v = end1%v + 2*sqrt(small*r)
   end1%u = sign(sqrt(r**2-(end1%v-center%v)**2),end1%u-center%u)+center%u
endif
if (end2%v < small*r .and. end2%v > -small*r) then
   end2%v = end2%v + 2*sqrt(small*r)
   end2%u = sign(sqrt(r**2-(end2%v-center%v)**2),end2%u-center%u)+center%u
endif

circle_crossing = 0

! compute the u coordinates of the intersections of the circle with the u axis

! if the absolute value of the v component of the center is bigger than the
! radius, then there are no intersections, and hence no crossings.  If it is
! exactly 0 then the circle is tangent to the u axis and we don't count that
! as a crossing.

if (abs(center%v) < r) then

! compute the intersections (u+,0) and (u-,0)

   discrim = sqrt(r**2 - center%v**2)
   uplus  = center%u + discrim
   uminus = center%u - discrim

! if u+ is negative, the arc cannot cross the positive u axis

   if (uplus > 0.0_my_real) then

! Consider the line segments between points and the center of the circle and
! angles between them, with all angles between 0 and pi, since GMSH uses
! the arc between end1 and end2 that is less than pi.  (u+,0) is on the arc iff
! the angle between it and end1 plus the angle between it and end2 equals the
! angle between end1 and end2.  Same for (u-,0).

      uplus_c  = point_uv(uplus -center%u,       -center%v)
      uminus_c = point_uv(uminus-center%u,       -center%v)
      end1_c   = point_uv(end1%u-center%u, end1%v-center%v)
      end2_c   = point_uv(end2%u-center%u, end2%v-center%v)

      norm_uplus   = sqrt( uplus_c%u**2 +  uplus_c%v**2)
      norm_uminus  = sqrt(uminus_c%u**2 + uminus_c%v**2)
      norm_end1    = sqrt(  end1_c%u**2 +   end1_c%v**2)
      norm_end2    = sqrt(  end2_c%u**2 +   end2_c%v**2)

      angle_uplus_end1  = acos((uplus_c%u*end1_c%u+uplus_c%v*end1_c%v) / &
                              (norm_uplus*norm_end1))
      angle_uplus_end2  = acos((uplus_c%u*end2_c%u+uplus_c%v*end2_c%v) / &
                              (norm_uplus*norm_end2))
      angle_uminus_end1 = acos((uminus_c%u*end1_c%u+uminus_c%v*end1_c%v) / &
                              (norm_uminus*norm_end1))
      angle_uminus_end2 = acos((uminus_c%u*end2_c%u+uminus_c%v*end2_c%v) / &
                              (norm_uminus*norm_end2))
      angle_end1_end2   = acos((end1_c%u*end2_c%u+end1_c%v*end2_c%v) / &
                              (norm_end1*norm_end2))

      uplus_in_arc  = abs(angle_uplus_end1+angle_uplus_end2-angle_end1_end2) &
                      < 10*small
      uminus_in_arc  = abs(angle_uminus_end1+angle_uminus_end2-angle_end1_end2)&
                      < 10*small

! We already verified that u+ is positive.  If u- is positive, then the circle
! lies to the right of the origin and we have one crossing iff exactly one
! of the intersections lies on the arc.  If u- is negative then (u+,0) is
! the only candidate for a crossing.

      if (uminus < 0.0_my_real) then
         if (uplus_in_arc) circle_crossing = 1
      else
         if (uplus_in_arc .neqv. uminus_in_arc) circle_crossing = 1
      endif

   endif ! uplus > 0
endif ! abs(center%v) < r

end function circle_crossing

!        ----------------
function ellipse_crossing(grid,line,vert,p1,p2,p3)
!        ----------------

!----------------------------------------------------
! This routine returns 1 if an edge that is an ellipse arc crosses the positive
! u axis once, and 0 if it doesn't cross or crosses twice.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: line, vert
type(point), intent(in) :: p1, p2, p3
integer :: ellipse_crossing
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

ierr = PHAML_INTERNAL_ERROR
call fatal("have not written ellipse_crossing yet")
stop

end function ellipse_crossing

!        -----
function to_uv(p,origin,p1,p2,p3)
!        -----

!----------------------------------------------------
! This routine returns the u-v coordinates of p in the coordinate system
! where origin is the origin and the axes are in the directions of p1-p3
! and p2-p3 orthogonalized to p1-p3.  p, origin, p1, p2 and p3 are assumed to
! lie in the same plane.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: origin, p, p1, p2, p3
type(point_uv) :: to_uv
!----------------------------------------------------
! Local variables:

type(point) :: d1, d2
!----------------------------------------------------
! Begin executable code

d1 = p1-p3
d1 = d1/sqrt(dot_point(d1,d1))
d2 = p2-p3 - dot_point(p2-p3,d1)*d1
d2 = d2/sqrt(dot_point(d2,d2))
to_uv%u = dot_point(d1,p-origin)
to_uv%v = dot_point(d2,p-origin)

end function to_uv

!=====================================================================
! Gmsh geometry file
!
! This section parses a .geo file to define the boundary surface.
!=====================================================================

!          -----------
subroutine process_geo(grid,procs,degree,partition,set_iown,geo_filename)
!          -----------

!----------------------------------------------------
! This routine reads a Gmsh .geo file to set boundary information in the grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
integer, intent(in) :: degree,partition
logical, intent(in) :: set_iown
character(len=*), intent(in) :: geo_filename
!----------------------------------------------------
! Local variables:

integer, parameter :: MAX_WORD_LEN = 80

integer :: iunit, iostat, astat
type(input_line) :: inline
character(len=MAX_WORD_LEN) :: word
type(var_name_tree) :: variables
type(hash_table) :: hash_point, hash_line, hash_lineloop
!----------------------------------------------------
! Begin executable code

! initialize the variable name tree to empty

call create_var_name_tree(variables)

! initialize the 3D boundary in the grid

allocate(grid%point(8), grid%line(8), grid%lineloop(8), grid%surface(8), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in process_geo")
   stop
endif
grid%npoint = 0
grid%nline = 0
grid%nlineloop = 0
grid%nsurface = 0

! create hash tables to relate my indexing to IDs in the .geo file

call hash_table_init(hash_point,1000)
call hash_table_init(hash_line, 1000)
call hash_table_init(hash_lineloop,1000)

! open the file

iunit = get_lun()
open(unit=iunit,file=geo_filename,action="READ",status="OLD",iostat=iostat)
if (iostat /= 0) then
   ierr = USER_INPUT_ERROR
   call fatal("failed to open .geo file in process_geo")
   stop
endif

! create the input unit stack

call create_stack(inline%stack)

! read the first line

call read_input_line(iunit,inline)

! repeat until the end of the input file

do while (inline%len /= -1)

! if current input is the beginning of a comment, move past it and any
! contiguous comments, and check for end of data

   call check_and_delete_comments(iunit,inline)
   if (inline%len == -1) exit

! get the next word

   call get_word(iunit,inline,word)

! process recognized words

   select case (word)

   case ("Include")
      call process_include(iunit,inline)
   case ("Point")
      call process_point(iunit,inline,variables,grid,hash_point)
   case ("Circle")
      call process_circle(iunit,inline,variables,grid,hash_point,hash_line)
   case ("Ellipse")
      call process_ellipse(iunit,inline,variables,grid,hash_point,hash_line)
   case ("Line")
      call process_line(iunit,inline,variables,grid,hash_point,hash_line)
   case ("LineLoop")
      call process_line_loop(iunit,inline,variables,grid,hash_line, &
                             hash_lineloop)
   case ("PlaneSurface")
      call process_plane_surface(iunit,inline,variables,grid,hash_lineloop)
   case ("RuledSurface")
      call process_ruled_surface(iunit,inline,variables,grid,hash_lineloop)
   case default

! if the next character is =, it is "string = expression"
! otherwise it is an unrecognized command to be ignored

      if (is_next_char(iunit,inline,"=")) then
         call process_assign(grid,iunit,inline,word,variables)
      else
         call ignore_rest_of_statement(iunit,inline)
      endif

   end select

end do ! main loop for parsing .geo file

! free memory

call destroy_stack(inline%stack)
call destroy_var_name_tree(variables)
call hash_table_destroy(hash_point)
call hash_table_destroy(hash_line)
call hash_table_destroy(hash_lineloop)


end subroutine process_geo

!          ---------------
subroutine read_input_line(iunit,inline)
!          ---------------

!----------------------------------------------------
! This routine reads the next input line.  inline%len==-1 indicates end of data
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! The stack keeps track of i/o unit numbers for Included files.
! END_OF_STACK indicates the end of the original file

do while (iunit /= END_OF_STACK)

! Read until we find a nonblank line

   inline%line = ""
   do while (len_trim(inline%line) == 0)
      read(iunit,fmt="(A)",end=100) inline%line
   end do

! Got the next line

   inline%ptr = 1
   inline%len = len_trim(inline%line)
   exit

! end of file, pop the Include stack and continue with the including file

100 continue
   close(iunit)
   iunit = pop_stack(inline%stack)

end do

! check for end of the original file

if (iunit == END_OF_STACK) then
   inline%line = ""
   inline%ptr = 1
   inline%len = -1
endif

call remove_tabs(inline%line)

end subroutine read_input_line

!          -----------
subroutine remove_tabs(string)
!          -----------

!----------------------------------------------------
! This routine replaces all tabs in string with a blank
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(inout) :: string
!----------------------------------------------------
! Local variables:

character(len=1) :: tab = achar(9)
integer :: tabloc
!----------------------------------------------------
! Begin executable code

tabloc = index(string,tab)
do while (tabloc /= 0)
   string(tabloc:tabloc) = " "
   tabloc = index(string,tab)
end do

end subroutine remove_tabs

!          --------
subroutine get_word(iunit,inline,word)
!          --------

!----------------------------------------------------
! This routine returns the next word, i.e. the character string consisting of
! all alphabetic characters up to the next non-alphabetic character except
! blank, with blanks removed.  It does not wrap the character string across
! lines.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
character(len=*), intent(out) :: word
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: word_ptr
!----------------------------------------------------
! Begin executable code

! initialize the word

word = ""
word_ptr = 0

! repeat until a non-alphabetic non-blank character or end of line is found

do
   if (inline%ptr <= inline%len) then

! next character

      char = inline%line(inline%ptr:inline%ptr)
      inline%ptr = inline%ptr+1

! check for alphabet and add to word if it is

      if ((char >= "a" .and. char <= "z") .or. &
          (char >= "A" .and. char <= "Z")) then
         word_ptr = word_ptr + 1
         word(word_ptr:word_ptr) = char

! check for non-blank and exit if it is but don't consume it.
! blank continues to next character

      elseif (char /= " ") then
         inline%ptr = inline%ptr - 1
         exit
      endif

! end of line; read next line and return

   else
      call read_input_line(iunit,inline)
      exit
   endif

end do

end subroutine get_word

!          ------------------------
subroutine ignore_rest_of_statement(iunit,inline)
!          ------------------------

!----------------------------------------------------
! This routine moves the input to the first character after the next semicolon
! that is not in a comment or character string.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
!----------------------------------------------------
! Local variables:

integer :: semicolon, c_comment, cpp_comment, quote
!----------------------------------------------------
! Begin executable code

! repeat until we find the end of statement

do

! check for end of data input

   if (inline%len == -1) exit

! location of the semicolon, if it is in this line

   semicolon = index(inline%line(inline%ptr:inline%len),";") + inline%ptr-1

! see if a C comment starts before the semicolon, and get rid of it if it does

   c_comment = index(inline%line(inline%ptr:inline%len),"/*") + inline%ptr-1
   if (c_comment /= inline%ptr-1 .and. c_comment < semicolon) then
      call ignore_c_comment(iunit,inline)
      cycle
   endif

! likewise for C++ comment

   cpp_comment = index(inline%line(inline%ptr:inline%len),"//") + inline%ptr-1
   if (cpp_comment /= inline%ptr-1 .and. cpp_comment < semicolon) then
      call ignore_cpp_comment(iunit,inline)
      cycle
   endif

! and also a string in double quotes

   quote = index(inline%line(inline%ptr:inline%len),'"') + inline%ptr-1
   if (quote /= inline%ptr-1 .and. quote < semicolon) then
      call ignore_quote(iunit,inline)
      cycle
   endif

! if there is no semicolon, read the next line and repeat

   if (semicolon == inline%ptr-1) then
      call read_input_line(iunit,inline)
      cycle
   endif

! if the semicolon is the end of the line, read the next line and exit

   if (semicolon == inline%len) then
      call read_input_line(iunit,inline)
      exit
   endif

! otherwise, move to the first character after the semicolon

   inline%ptr = semicolon+1
   exit

end do

end subroutine ignore_rest_of_statement

!          -------------------------
subroutine check_and_delete_comments(iunit,inline)
!          -------------------------

!----------------------------------------------------
! This routine check to see if the input is currently at the beginning
! of a C or C++ style comment, and if so, moves to the first character
! after the comment and repeats until it is not starting a comment.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
!----------------------------------------------------
! Local variables:

integer :: c_comment, cpp_comment
!----------------------------------------------------
! Begin executable code

! repeat until we don't find a comment

do

! look for the start of a comment

   c_comment   = index(inline%line(inline%ptr:inline%len),"/*") + inline%ptr-1
   cpp_comment = index(inline%line(inline%ptr:inline%len),"//") + inline%ptr-1

! if no comment, then we are done

   if (c_comment == inline%ptr-1 .and. cpp_comment == inline%ptr-1) exit

! C comment starts first

   if (c_comment /= inline%ptr-1 .and. &
       (cpp_comment == inline%ptr-1 .or. c_comment < cpp_comment)) then

! if there are nonblank characters before the comment, then the comment is not
! the next thing and we are done

      if (inline%line(inline%ptr:c_comment-1) /= " ") exit

! otherwise, get rid of the comment and repeat

      call ignore_c_comment(iunit,inline)
! note we don't need to continue if we reach end of data
      if (inline%len == -1) exit

! C++ comment starts first

   elseif (cpp_comment /= inline%ptr-1 .and. &
       (c_comment == inline%ptr-1 .or. cpp_comment < c_comment)) then
      if (inline%line(inline%ptr:cpp_comment-1) /= " ") exit
      call ignore_cpp_comment(iunit,inline)
      if (inline%len == -1) exit

   endif
end do

end subroutine check_and_delete_comments

!          ----------------
subroutine ignore_c_comment(iunit,inline)
!          ----------------

!----------------------------------------------------
! This routine moves the input point to the first character after */
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
!----------------------------------------------------
! Local variables:

integer :: symbol
!----------------------------------------------------
! Begin executable code

! repeat until end is found

do

! look for the symbol

   symbol = index(inline%line(inline%ptr:inline%len),"*/") + inline%ptr-1
   if (symbol /= inline%ptr-1) exit

! not there, try next line

   call read_input_line(iunit,inline)
   if (inline%len == -1) return

end do

! if the symbol is the end of the line, read the next line
! otherwise, move to the first character after the symbol

if (symbol == inline%len-1) then
   call read_input_line(iunit,inline)
else
   inline%ptr = symbol+2
endif

end subroutine ignore_c_comment

!          ------------------
subroutine ignore_cpp_comment(iunit,inline)
!          ------------------

!----------------------------------------------------
! This routine moves to the beginning of the next line
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call read_input_line(iunit,inline)

end subroutine ignore_cpp_comment

!          ------------
subroutine ignore_quote(iunit,inline)
!          ------------

!----------------------------------------------------
! This routine moves the input point to the first character after "
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
!----------------------------------------------------
! Local variables:

integer :: symbol
!----------------------------------------------------
! Begin executable code

! repeat until end is found

do

! look for the symbol

   symbol = index(inline%line(inline%ptr:inline%len),'"') + inline%ptr-1
   if (symbol /= inline%ptr-1) exit

! not there, try next line

   call read_input_line(iunit,inline)
   if (inline%len == -1) return

end do

! if the symbol is the end of the line, read the next line
! otherwise, move to the first character after the symbol

if (symbol == inline%len) then
   call read_input_line(iunit,inline)
else
   inline%ptr = symbol+1
endif

end subroutine ignore_quote

!        ------------
function is_next_char(iunit,inline,char)
!        ------------

!----------------------------------------------------
! This routine returns true if the next nonblank character is char
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
character(len=1), intent(in) :: char
logical :: is_next_char
!----------------------------------------------------
! Local variables:

integer :: char_loc
!----------------------------------------------------
! Begin executable code

call check_and_delete_comments(iunit,inline)

char_loc = index(inline%line(inline%ptr:inline%len),char) + inline%ptr-1
if (char_loc == inline%ptr-1) then ! not there
   is_next_char = .false.
elseif (inline%line(inline%ptr:char_loc-1) == " ") then ! all blanks before it
   is_next_char = .true.
else ! some non-blank before it
   is_next_char = .false.
endif

end function is_next_char

!          ---------------------
subroutine consume_expected_char(iunit,inline,char)
!          ---------------------

!----------------------------------------------------
! This routine is called when char is expected to be the next character.
! If it is, it moves the pointer to the first character after char.
! If not, it generates an error.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
character(len=1), intent(in) :: char
!----------------------------------------------------
! Local variables:

integer :: char_loc
!----------------------------------------------------
! Begin executable code

! find char in the input

! repeat until we find it

char_loc = inline%ptr-1
do while (char_loc == inline%ptr-1)

! remove any comments that might precede the expected character

   call check_and_delete_comments(iunit,inline)

! look for char in this input line

   char_loc = index(inline%line(inline%ptr:inline%len),char) + inline%ptr-1
   if (char_loc /= inline%ptr-1) exit

! it is not in the current line of input

! make sure nothing else is

   if (inline%line(inline%ptr:inline%len) /= " ") then
      ierr = USER_INPUT_ERROR
      call fatal("in .geo file, expecting "//char,inline%line(1:inline%len), &
                 intlist=(/inline%ptr/))
      stop
   endif

! and get the next input line

   call read_input_line(iunit,inline)

end do

! move the pointer to the character after char

inline%ptr = char_loc + 1
if (inline%ptr > inline%len) then
   call read_input_line(iunit,inline)
endif

end subroutine consume_expected_char

!          ---------------------
subroutine get_expression_string(iunit,inline,expression_str)
!          ---------------------

!----------------------------------------------------
! This routine returns a character string with the expression starting at the
! current location of the input.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
character(len=*) :: expression_str
!----------------------------------------------------
! Local variables:

integer :: paren, comma, brace, semicolon, delimiter, paren_balance, i
!----------------------------------------------------
! Begin executable code

! TEMP I don't know how to do this yet.  For now, return everything before
! the next ) , } ; and assume it does not contain comments and is not split
! over more than one line.

! remove any comments that might precede the expression

call check_and_delete_comments(iunit,inline)

! find the location of the delimiters

comma = index(inline%line(inline%ptr:inline%len),",") + inline%ptr-1
brace = index(inline%line(inline%ptr:inline%len),"}") + inline%ptr-1
semicolon = index(inline%line(inline%ptr:inline%len),";") + inline%ptr-1

! the parenthesis is harder because it can be in an expression.  Look for the
! first closing parenthesis that is not balanced by an opening parenthesis

paren_balance = 0
paren = inline%ptr-1
do i=inline%ptr,inline%len
   if (inline%line(i:i) == "(") then
      paren_balance = paren_balance + 1
   elseif (inline%line(i:i) == ")") then
      if (paren_balance == 0) then
         paren = i
         exit
      else
         paren_balance = paren_balance - 1
      endif
   endif
end do

! determine which is the active delimiter

delimiter = inline%ptr - 1
if (paren /= inline%ptr-1) delimiter = paren
if (comma /= inline%ptr-1) then
   if (delimiter == inline%ptr-1) then
      delimiter = comma
   else
      if (comma < delimiter) delimiter = comma
   endif
endif
if (brace /= inline%ptr-1) then
   if (delimiter == inline%ptr-1) then
      delimiter = brace
   else
      if (brace < delimiter) delimiter = brace
   endif
endif
if (semicolon /= inline%ptr-1) then
   if (delimiter == inline%ptr-1) then
      delimiter = semicolon
   else
      if (semicolon < delimiter) delimiter = semicolon
   endif
endif

if (delimiter == inline%ptr-1) then
   ierr = USER_INPUT_ERROR
   call fatal("failed to find end of expression in get_expression_string")
   stop
endif

if (delimiter - inline%ptr > len(expression_str)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("get_expression_string: expression_str is not long enough", &
              "increase MAX_EXPRESSION in init_grid3D.f90 and recompile")
   stop
endif

expression_str = inline%line(inline%ptr:delimiter-1)
inline%ptr = delimiter

end subroutine get_expression_string

!        --------------
function get_expression(iunit,inline,variables)
!        --------------

!----------------------------------------------------
! This routine parses and evaluates the expression currently next in the
! input, and returns a real number.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
real(my_real) :: get_expression
type(var_name_tree), intent(in) :: variables
!----------------------------------------------------
! Local variables:

character(len=MAX_EXPRESSION) :: expression
!----------------------------------------------------
! Begin executable code

! get and parse the expression

call get_expression_string(iunit,inline,expression)
get_expression = parse_expression(expression,variables)

end function get_expression

!          -------------------
subroutine get_expression_list(iunit,inline,variables,explist,nexp)
!          -------------------

!----------------------------------------------------
! This routine evaluates a comma-separated list of expressions.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
real(my_real), intent(out) :: explist(:)
integer,intent(out) :: nexp
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! LIMITATION: an expression-list-item must be an expression.

! process expressions as long as the character following the expression
! is a comma

nexp = 1
do
   explist(nexp) = get_expression(iunit,inline,variables)
   if (.not. is_next_char(iunit,inline,",")) exit
   call consume_expected_char(iunit,inline,",")
   nexp = nexp + 1
   if (nexp > size(explist)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Expression list has more items than size of explist.", &
                 "Increase MAX_EXP_LIST in grid_init3D.f90 and recompile.")
      stop
   endif
end do

end subroutine get_expression_list

!          ---------------
subroutine process_include(iunit,inline)
!          ---------------

!----------------------------------------------------
! This routine processes an include statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
!----------------------------------------------------
! Local variables:

integer :: quote1, quote2, iostat
!----------------------------------------------------
! Begin executable code

! LIMITATION: the char-expression must be a string enclosed in
! double quotes.  The string must be a file name that can be passed to the
! Fortran OPEN statement, which usually means a file in the current directory
! or a file with the full path to it.

! find the first " and move to the first character after it

quote1 = index(inline%line(inline%ptr:inline%len),'"') + inline%ptr-1
if (quote1 == inline%ptr-1) then
   ierr = USER_INPUT_ERROR
   call fatal("Include statement in .geo file does not have two double quotes")
   stop
endif
inline%ptr = quote1 + 1

! find the second quote

quote2 = index(inline%line(inline%ptr:inline%len),'"') + inline%ptr-1
if (quote2 == inline%ptr-1) then
   ierr = USER_INPUT_ERROR
   call fatal("Include statement in .geo file does not have two double quotes")
   stop
endif

! push the current logical unit number on the unit number stack, get a new
! LUN, open the file and read the first line

call push_stack(inline%stack,iunit)
iunit = get_lun()
open(unit=iunit,file=inline%line(quote1+1:quote2-1),action="READ", &
     status="OLD",iostat=iostat)
if (iostat /= 0) then
   ierr = USER_INPUT_ERROR
   call fatal("failed to open Include file when processing .geo file")
   stop
endif
call read_input_line(iunit,inline)

end subroutine process_include

!          -------------
subroutine process_point(iunit,inline,variables,grid,hash_point)
!          -------------

!----------------------------------------------------
! This routine processes a Point statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
type(grid_type), intent(inout) :: grid
type(hash_table), intent(inout) :: hash_point
!----------------------------------------------------
! Local variables:

real(my_real) :: exp
type(hash_key) :: id
!----------------------------------------------------
! Begin executable code

! Point(expression) = {expression,expression,expression[,expression]};

! get ID of this point

call consume_expected_char(iunit,inline,"(")
exp = get_expression(iunit,inline,variables)
call consume_expected_char(iunit,inline,")")

! add this point to hash table

grid%npoint = grid%npoint + 1
if (grid%npoint > size(grid%point)) then
   call more_points(grid)
endif
id = nint(exp)
call hash_insert(id,grid%npoint,hash_point)

! get and store the three coordinates of the point

call consume_expected_char(iunit,inline,"=")
call consume_expected_char(iunit,inline,"{")
exp = get_expression(iunit,inline,variables)
grid%point(grid%npoint)%coord%x = exp
call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
grid%point(grid%npoint)%coord%y = exp
call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
grid%point(grid%npoint)%coord%z = exp

! we ignore the fourth expression if it is present, so at this point we can
! just ignore the rest of the statement

call ignore_rest_of_statement(iunit,inline)

end subroutine process_point

!          --------------
subroutine process_circle(iunit,inline,variables,grid,hash_point,hash_line)
!          --------------

!----------------------------------------------------
! This routine processes a Circle statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
type(grid_type), intent(inout) :: grid
type(hash_table), intent(in) :: hash_point
type(hash_table), intent(inout) :: hash_line
!----------------------------------------------------
! Local variables:

real(my_real) :: exp
integer :: astat, point
type(hash_key) :: id
!----------------------------------------------------
! Begin executable code

! Circle(expression) = {expression,expression,expression};

! get ID of this line

call consume_expected_char(iunit,inline,"(")
exp = get_expression(iunit,inline,variables)
call consume_expected_char(iunit,inline,")")

! add this line to hash table

grid%nline = grid%nline + 1
if (grid%nline > size(grid%line)) then
   call more_lines(grid)
endif
id = nint(exp)
call hash_insert(id,grid%nline,hash_line)

! set the line type and allocate the list of points

grid%line(grid%nline)%type = LINE_CIRCLE
allocate(grid%line(grid%nline)%point(3),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in process_circle")
   stop
endif

! get the ID of each point that defines this line, find its index
! in grid%point, and store it

call consume_expected_char(iunit,inline,"=")
call consume_expected_char(iunit,inline,"{")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(1) = point

call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(2) = point

call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(3) = point

call consume_expected_char(iunit,inline,"}")
call consume_expected_char(iunit,inline,";")

end subroutine process_circle

!          ---------------
subroutine process_ellipse(iunit,inline,variables,grid,hash_point,hash_line)
!          ---------------

!----------------------------------------------------
! This routine processes a Ellipse statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
type(grid_type), intent(inout) :: grid
type(hash_table), intent(in) :: hash_point
type(hash_table), intent(inout) :: hash_line
!----------------------------------------------------
! Local variables:

real(my_real) :: exp
integer :: astat, point
type(hash_key) :: id
!----------------------------------------------------
! Begin executable code

! Ellipse(expression) = {expression,expression,expression,expression};

! get ID of this line

call consume_expected_char(iunit,inline,"(")
exp = get_expression(iunit,inline,variables)
call consume_expected_char(iunit,inline,")")

! add this line to hash table

grid%nline = grid%nline + 1
if (grid%nline > size(grid%line)) then
   call more_lines(grid)
endif
id = nint(exp)
call hash_insert(id,grid%nline,hash_line)

! set the line type and allocate the list of points

grid%line(grid%nline)%type = LINE_ELLIPSE
allocate(grid%line(grid%nline)%point(4),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in process_ellipse")
   stop
endif

! get the ID of each point that defines this line, find its index
! in grid%point, and store it

call consume_expected_char(iunit,inline,"=")
call consume_expected_char(iunit,inline,"{")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(1) = point

call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(2) = point

call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(3) = point

call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(4) = point

call consume_expected_char(iunit,inline,"}")
call consume_expected_char(iunit,inline,";")

end subroutine process_ellipse

!          ------------
subroutine process_line(iunit,inline,variables,grid,hash_point,hash_line)
!          ------------

!----------------------------------------------------
! This routine processes a Line statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
type(grid_type), intent(inout) :: grid
type(hash_table), intent(in) :: hash_point
type(hash_table), intent(inout) :: hash_line
!----------------------------------------------------
! Local variables:

real(my_real) :: exp
integer :: astat, point
type(hash_key) :: id
!----------------------------------------------------
! Begin executable code

! Line(expression) = {expression,expression};

! get ID of this line

call consume_expected_char(iunit,inline,"(")
exp = get_expression(iunit,inline,variables)
call consume_expected_char(iunit,inline,")")

! add this line to hash table

grid%nline = grid%nline + 1
if (grid%nline > size(grid%line)) then
   call more_lines(grid)
endif
id = nint(exp)
call hash_insert(id,grid%nline,hash_line)

! set the line type and allocate the list of points

grid%line(grid%nline)%type = LINE_LINE
allocate(grid%line(grid%nline)%point(2),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in process_line")
   stop
endif

! get the ID of each point that defines this line, find its index
! in grid%point, and store it

call consume_expected_char(iunit,inline,"=")
call consume_expected_char(iunit,inline,"{")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(1) = point

call consume_expected_char(iunit,inline,",")
exp = get_expression(iunit,inline,variables)
id = nint(exp)
point = hash_decode_key(id,hash_point)
if (point == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a point must be defined before it is used", &
              intlist=(/nint(exp)/))
   stop
endif
grid%line(grid%nline)%point(2) = point

call consume_expected_char(iunit,inline,"}")
call consume_expected_char(iunit,inline,";")

end subroutine process_line

!          -----------------
subroutine process_line_loop(iunit,inline,variables,grid,hash_line, &
                             hash_lineloop)
!          -----------------

!----------------------------------------------------
! This routine processes a Line Loop statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
type(grid_type), intent(inout) :: grid
type(hash_table), intent(in) :: hash_line
type(hash_table), intent(inout) :: hash_lineloop
!----------------------------------------------------
! Local variables:

real(my_real) :: exp, explist(MAX_EXP_LIST)
integer :: nexp
integer :: i, line, astat
type(hash_key) :: id
!----------------------------------------------------
! Begin executable code

! Line Loop(expression) = {expression-list};

! get ID of this line loop

call consume_expected_char(iunit,inline,"(")
exp = get_expression(iunit,inline,variables)
call consume_expected_char(iunit,inline,")")

! add this line loop to hash table

grid%nlineloop = grid%nlineloop + 1
if (grid%nlineloop > size(grid%lineloop)) then
   call more_lineloops(grid)
endif
id = nint(exp)
call hash_insert(id,grid%nlineloop,hash_lineloop)

! get the list of lines that make up this line loop

call consume_expected_char(iunit,inline,"=")
call consume_expected_char(iunit,inline,"{")
call get_expression_list(iunit,inline,variables,explist,nexp)

! allocate the list of lines

allocate(grid%lineloop(grid%nlineloop)%line(nexp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in process_line")
   stop
endif

! for each line that defines this line loop, find its index
! in grid%line and store it

do i=1,nexp
   id = abs(nint(explist(i)))
   line = hash_decode_key(id,hash_line)
   if (line == HASH_NOT_FOUND) then
      ierr = USER_INPUT_ERROR
      call fatal("in a .geo file, a line must be defined before it is used", &
                 intlist=(/abs(nint(explist(i)))/))
      stop
   endif
   grid%lineloop(grid%nlineloop)%line(i) = sign(line,nint(explist(i)))
end do

call consume_expected_char(iunit,inline,"}")
call consume_expected_char(iunit,inline,";")

end subroutine process_line_loop

!          ---------------------
subroutine process_plane_surface(iunit,inline,variables,grid,hash_lineloop)
!          ---------------------

!----------------------------------------------------
! This routine processes a Plane Surface statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
type(grid_type), intent(inout) :: grid
type(hash_table), intent(in) :: hash_lineloop
!----------------------------------------------------
! Local variables:

real(my_real) :: exp, explist(MAX_EXP_LIST)
integer :: nexp, lineloop, ilineloop, astat
type(hash_key) :: id
!----------------------------------------------------
! Begin executable code

! Plane Surface(expression) = {expression-list};

! get the ID of this surface

call consume_expected_char(iunit,inline,"(")
exp = get_expression(iunit,inline,variables)
call consume_expected_char(iunit,inline,")")

! get the index of this surface in grid

grid%nsurface = grid%nsurface + 1
if (grid%nsurface > size(grid%surface)) then
   call more_surfaces(grid)
endif

! get the list of line loops that define this surface

call consume_expected_char(iunit,inline,"=")
call consume_expected_char(iunit,inline,"{")
call get_expression_list(iunit,inline,variables,explist,nexp)

! set the type of surface

grid%surface(grid%nsurface)%type = SURF_PLANAR

! get the IDs of the line loops that defines this surface, find their index
! in grid%lineloop, and store them

allocate(grid%surface(grid%nsurface)%lineloop(nexp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in process_plane_surface")
   stop
endif

do ilineloop=1,nexp
   id = nint(explist(ilineloop))
   lineloop = hash_decode_key(id,hash_lineloop)
   if (lineloop == HASH_NOT_FOUND) then
      ierr = USER_INPUT_ERROR
      call fatal("in a .geo file, a line loop must be defined before it is used", &
                 intlist=(/nint(explist(ilineloop))/))
      stop
   endif
   grid%surface(grid%nsurface)%lineloop(ilineloop) = lineloop
end do

call consume_expected_char(iunit,inline,"}")
call consume_expected_char(iunit,inline,";")

end subroutine process_plane_surface

!          ---------------------
subroutine process_ruled_surface(iunit,inline,variables,grid,hash_lineloop)
!          ---------------------

!----------------------------------------------------
! This routine processes a Ruled Surface statement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
type(var_name_tree), intent(in) :: variables
type(grid_type), intent(inout) :: grid
type(hash_table), intent(in) :: hash_lineloop
!----------------------------------------------------
! Local variables:

real(my_real) :: exp, explist(MAX_EXP_LIST)
integer :: nexp, lineloop, astat
type(hash_key) :: id
!----------------------------------------------------
! Begin executable code

! LIMITATION: the Ruled Surface statement must not have an In Sphere clause.
! If it does, the clause is ignored.
! LIMITATION: no holes in ruled surfaces

! Ruled Surface(expression) = {expression-list};

! get the ID of this surface

call consume_expected_char(iunit,inline,"(")
exp = get_expression(iunit,inline,variables)
call consume_expected_char(iunit,inline,")")

! get the index of this surface in grid

grid%nsurface = grid%nsurface + 1
if (grid%nsurface > size(grid%surface)) then
   call more_surfaces(grid)
endif

! get the list of line loops that define this surface

call consume_expected_char(iunit,inline,"=")
call consume_expected_char(iunit,inline,"{")
call get_expression_list(iunit,inline,variables,explist,nexp)
call consume_expected_char(iunit,inline,"}")

! check for holes

if (nexp > 1) then
   ierr = USER_INPUT_ERROR
   call fatal("holes in ruled surfaces are not supported")
   stop
endif

! get the ID of the line loop that defines this surface, find its index
! in grid%lineloop, and store it

allocate(grid%surface(grid%nsurface)%lineloop(1),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in process_ruled_surface")
   stop
endif

id = nint(explist(1))
lineloop = hash_decode_key(id,hash_lineloop)
if (lineloop == HASH_NOT_FOUND) then
   ierr = USER_INPUT_ERROR
   call fatal("in a .geo file, a line loop must be defined before it is used", &
              intlist=(/nint(explist(1))/))
   stop
endif
grid%surface(grid%nsurface)%lineloop(1) = lineloop

! set the type of surface

if (size(grid%lineloop(lineloop)%line) == 3) then
   grid%surface(grid%nsurface)%type = SURF_RULED_3
elseif (size(grid%lineloop(lineloop)%line) == 4) then
   grid%surface(grid%nsurface)%type = SURF_RULED_4
else
   ierr = USER_INPUT_ERROR
   call fatal("the line loop for a ruled surface must have either 3 or 4 lines")
   stop
endif

! if the next character is "I", then it is probably the beginning of an
! In Sphere clause, which is not allowed.

if (is_next_char(iunit,inline,"I")) then
   call warning("The In Sphere clause in a .geo file is not supported, and will be ignored.")
endif

! whether or not there is an In Sphere clause, ignore the rest of the statement

call ignore_rest_of_statement(iunit,inline)

end subroutine process_ruled_surface

!          --------------
subroutine process_assign(grid,iunit,inline,word,variables)
!          --------------

!----------------------------------------------------
! This routine processes an assignment statement.  word has already been
! read as the variable name.  It has been established that the next
! character is = but it has not been read.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(inout) :: iunit
type(input_line), intent(inout) :: inline
character(len=*) :: word
type(var_name_tree), intent(inout) :: variables
!----------------------------------------------------
! Local variables:

real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! consume the =

call consume_expected_char(iunit,inline,"=")

! evaluate the expression

value = get_expression(iunit,inline,variables)

! add the new variable and value to the variable name tree (or update
! the value of an existing variable)

call insert_var_name_tree(trim(word),variables,value)

! TEMP120327 keep certain values from the contactHole domain to be used
!            for fixing points that are supposed to be on a circular boundary

select case (word)
case ("holeHeight"); holeHeight = value
case ("holeBottomRadius"); holeBottomRadius = value
case ("holeTopRadius"); holeTopRadius = value
case ("cylinderRadius"); cylinderRadius = value
case ("CYLINDERSIDE"); CYLINDERSIDE = value
case ("HOLESIDE"); HOLESIDE = value
end select
! end TEMP120327

! there should be just a semicolon remaining

call consume_expected_char(iunit,inline,";")

end subroutine process_assign

!          -----------
subroutine more_points(grid)
!          -----------

!----------------------------------------------------
! This routine doubles the allocation of grid%point
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

type(point_type), pointer :: point(:)
integer :: astat
!----------------------------------------------------
! Begin executable code

allocate(point(2*size(grid%point)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in more_points")
   stop
endif
point(1:size(grid%point)) = grid%point
deallocate(grid%point)
grid%point => point

end subroutine more_points

!          ----------
subroutine more_lines(grid)
!          ----------

!----------------------------------------------------
! This routine doubles the allocation of grid%line
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

type(line_type), pointer :: line(:)
integer :: astat
!----------------------------------------------------
! Begin executable code

allocate(line(2*size(grid%line)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in more_lines")
   stop
endif
line(1:size(grid%line)) = grid%line
deallocate(grid%line)
grid%line => line

end subroutine more_lines

!          --------------
subroutine more_lineloops(grid)
!          --------------

!----------------------------------------------------
! This routine doubles the allocation of grid%lineloop
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

type(lineloop_type), pointer :: lineloop(:)
integer :: astat
!----------------------------------------------------
! Begin executable code

allocate(lineloop(2*size(grid%lineloop)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in more_lineloops")
   stop
endif
lineloop(1:size(grid%lineloop)) = grid%lineloop
deallocate(grid%lineloop)
grid%lineloop => lineloop

end subroutine more_lineloops

!          -------------
subroutine more_surfaces(grid)
!          -------------

!----------------------------------------------------
! This routine doubles the allocation of grid%surface
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

type(surf_type), pointer :: surface(:)
integer :: astat
!----------------------------------------------------
! Begin executable code

allocate(surface(2*size(grid%surface)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in more_surfaces")
   stop
endif
surface(1:size(grid%surface)) = grid%surface
deallocate(grid%surface)
grid%surface => surface

end subroutine more_surfaces

!          --------------------
subroutine create_var_name_tree(tree)
!          --------------------

!----------------------------------------------------
! This routine initializes a variable name tree to empty
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(var_name_tree), intent(out) :: tree
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

tree%string = ""
tree%value = 0
nullify(tree%left)
nullify(tree%right)

end subroutine create_var_name_tree

!                    ---------------------
recursive subroutine destroy_var_name_tree(tree)
!                    ---------------------

!----------------------------------------------------
! This routine frees the memory in a variable name tree
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(var_name_tree), intent(inout) :: tree
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (associated(tree%left)) then
   call destroy_var_name_tree(tree%left)
   deallocate(tree%left)
endif

if (associated(tree%right)) then
   call destroy_var_name_tree(tree%right)
   deallocate(tree%right)
endif

end subroutine destroy_var_name_tree

!                    --------------------
recursive subroutine insert_var_name_tree(string,tree,value)
!                    --------------------

!----------------------------------------------------
! This routine inserts (string,value) into a variable name tree.  If string
! is already in the tree, value is replaced by the given value.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: string
type(var_name_tree), intent(inout) :: tree
real(my_real), intent(in) :: value

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! make sure string isn't too long

if (len(string) > MAX_VAR_LEN) then
   ierr = USER_INPUT_ERROR
   call fatal("variable name too long in .geo file")
   stop
endif

if (tree%string == "") then

! this is the first entry in the tree; store it here

   tree%string = string
   tree%value  = value

elseif (tree%string == string) then

! found the string already in the tree; just assign value

   tree%value = value

elseif (string < tree%string) then

! string goes in left subtree.
! if left subtree exists, descend it; otherwise insert string there

   if (associated(tree%left)) then
      call insert_var_name_tree(string,tree%left,value)
   else
      allocate(tree%left)
      tree%left%string = string
      tree%left%value  = value
      nullify(tree%left%left)
      nullify(tree%left%right)
   endif

else

! string goes in right subtree.
! if right subtree exists, descend it; otherwise insert string there

   if (associated(tree%right)) then
      call insert_var_name_tree(string,tree%right,value)
   else
      allocate(tree%right)
      tree%right%string = string
      tree%right%value  = value
      nullify(tree%right%left)
      nullify(tree%right%right)
   endif

endif

end subroutine insert_var_name_tree

!                  ------------------
recursive function find_var_name_tree(string,tree) result (value)
!                  ------------------

!----------------------------------------------------
! This routine returns the value associated with string in a variable name tree
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: string
type(var_name_tree), intent(in) :: tree
real(my_real) :: value
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (tree%string == string) then

! found the string, return value

   value = tree%value

elseif (string < tree%string) then

! string should be in left subtree; error if it doesn't exist

   if (associated(tree%left)) then
      value = find_var_name_tree(string,tree%left)
   else
      ierr = USER_INPUT_ERROR
      call fatal("expression contains undefined variable in .geo file")
      stop
   endif

else

! string should be in right subtree; error if it doesn't exist

   if (associated(tree%right)) then
      value = find_var_name_tree(string,tree%right)
   else
      ierr = USER_INPUT_ERROR
      call fatal("expression contains undefined variable in .geo file")
      stop
   endif

endif

end function find_var_name_tree

!=====================================================================
! Parse expression
!
! This section parses an expression.
! Expressions are limited to
!   expression: number |
!               variable |
!               ( expression ) |
!               operator-unary-left expression |
!               expression operator-binary expression
!=====================================================================

!        ----------------
function parse_expression(expression,variables)
!        ----------------

!----------------------------------------------------
! This is the top level routine to parse expression and return a real number
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*) :: expression
type(var_name_tree), intent(in) :: variables
real(my_real) :: parse_expression
!----------------------------------------------------
! Local variables:

integer :: state, exp_ptr, num_ptr, var_ptr
type(stack_integer) :: operator_stack
type(stack_my_real) :: value_stack
character(len=MAX_NUM_LEN) :: number_string
character(len=MAX_VAR_LEN) :: variable_string
!----------------------------------------------------
! Begin executable code

! initialize the operator and value stacks to empty

call create_stack(operator_stack)
call create_stack(value_stack)

! begin in the start_expression state, at the beginning of the expression
! and with empty number and variable strings

state = START_EXPRESSION
exp_ptr = 0
number_string = ""
num_ptr = 0
variable_string = ""
var_ptr = 0

! repeat until end of expression; each routine processes the next character for
! the current state and returns the next state

do
   select case(state)
   case(START_EXPRESSION)
      state = state_start_expression(expression,exp_ptr,operator_stack, &
                                     value_stack,variables,number_string, &
                                     num_ptr,variable_string,var_ptr)
   case(NUMBER)
      state = state_number(expression,exp_ptr,operator_stack,value_stack, &
                           variables,number_string,num_ptr,variable_string, &
                           var_ptr)
   case(SIGNIFICAND)
      state = state_significand(expression,exp_ptr,operator_stack,value_stack, &
                                variables,number_string,num_ptr, &
                                variable_string,var_ptr)
   case(EXPONENT)
      state = state_exponent(expression,exp_ptr,operator_stack,value_stack, &
                             variables,number_string,num_ptr,variable_string, &
                             var_ptr)
   case(EXPONENT_RIGHT_AFTER_SIGN)
      state = state_exponent_right_after_sign(expression,exp_ptr, &
                                              operator_stack,value_stack, &
                                              variables,number_string,num_ptr, &
                                              variable_string,var_ptr)
   case(EXPONENT_AFTER_SIGN)
      state = state_exponent_after_sign(expression,exp_ptr,operator_stack, &
                                        value_stack,variables,number_string, &
                                        num_ptr,variable_string,var_ptr)
   case(VARIABLE)
      state = state_variable(expression,exp_ptr,operator_stack,value_stack, &
                             variables,number_string,num_ptr,variable_string, &
                             var_ptr)
   case(END_PAREN_EXPRESSION)
      state = state_end_paren_expression(expression,exp_ptr,operator_stack, &
                                         value_stack,variables,number_string, &
                                         num_ptr,variable_string,var_ptr)
   case(FINISH)
      parse_expression = state_finish(expression,exp_ptr,operator_stack, &
                                      value_stack,variables,number_string, &
                                      num_ptr,variable_string,var_ptr)
      exit
   end select
end do

! clean up memory

call destroy_stack(operator_stack)
call destroy_stack(value_stack)

end function parse_expression

!          --------
subroutine get_char(expression,exp_ptr,char,char_type)
!          --------

!----------------------------------------------------
! This routine gets the next nonblank character in the expression and its type
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
character(len=1), intent(out) :: char
integer, intent(out) :: char_type
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (exp_ptr >= len(trim(expression))) then
   char = " "
   char_type = END_OF_EXPRESSION
else
   char = " "
   do while (char == " ")
      exp_ptr = exp_ptr + 1
      char = expression(exp_ptr:exp_ptr)
   end do
   if (char >= "0" .and. char <= "9") then
      char_type = DIGIT
   elseif (char == "d" .or. char == "e" .or. char == "D" .or. char == "E") then
      char_type = EXPONENT_LETTER
   elseif ((char >= "a" .and. char <= "z") .or. &
           (char >= "A" .and. char <= "Z")) then
      char_type = LETTER_NOT_EXPONENT
   elseif (char == "^" .or. char == "*" .or. char == "/" .or. char == "%") then
      char_type = OPERATOR_BINARY_NOT_PLUS_MINUS
   elseif (char == "+") then
      char_type = PLUS
   elseif (char == "-") then
      char_type = MINUS
   elseif (char == ".") then
      char_type = DOT
   elseif (char == "_") then
      char_type = UNDERSCORE
   elseif (char == "(") then
      char_type = OPEN_PAREN
   elseif (char == ")") then
      char_type = CLOSE_PAREN
   else
      ierr = USER_INPUT_ERROR
      call fatal("unrecognized character in expression in .geo file", &
                 trim(expression),(/exp_ptr/))
      stop
   endif
endif

end subroutine get_char

!        ----------------------
function state_start_expression(expression,exp_ptr,operator_stack,value_stack, &
                                variables,number_string,num_ptr, &
                                variable_string,var_ptr) result(state)
!        ----------------------

!----------------------------------------------------
! This routine processes the start-expression state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "start-expression"
expecting = "digit, letter, minus, dot, open-paren"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = NUMBER

case(LETTER_NOT_EXPONENT, EXPONENT_LETTER)
   var_ptr = var_ptr + 1
   variable_string(var_ptr:var_ptr) = char
   state = VARIABLE

case(MINUS)
   call push_stack(operator_stack,UNARY_MINUS)
   state = START_EXPRESSION

case(OPERATOR_BINARY_NOT_PLUS_MINUS, PLUS)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(DOT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = SIGNIFICAND

case(UNDERSCORE)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPEN_PAREN)
   call push_stack(operator_stack,BEGIN_PAREN)
   state = START_EXPRESSION

case(CLOSE_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(END_OF_EXPRESSION)
   call expression_error(state_string,expression,exp_ptr,"end of expression", &
                         expecting)

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_start_expression

!        ------------
function state_number(expression,exp_ptr,operator_stack,value_stack, &
                      variables,number_string,num_ptr, &
                      variable_string,var_ptr) result(state)
!        ------------

!----------------------------------------------------
! This routine processes the number state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "number"
expecting = "digit, exponent-letter, operator-binary, dot, close-paren, end-of-expression"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = NUMBER

case(LETTER_NOT_EXPONENT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(EXPONENT_LETTER)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = EXPONENT

case(OPERATOR_BINARY_NOT_PLUS_MINUS, PLUS, MINUS)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_no_paren(operator_stack,value_stack, &
                                char_to_operator(char))
   state = START_EXPRESSION

case(DOT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = SIGNIFICAND

case(UNDERSCORE)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPEN_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(CLOSE_PAREN)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_with_paren(operator_stack,value_stack)
   state = END_PAREN_EXPRESSION

case(END_OF_EXPRESSION)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   state = FINISH

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_number

!        -----------------
function state_significand(expression,exp_ptr,operator_stack,value_stack, &
                           variables,number_string,num_ptr, &
                           variable_string,var_ptr) result(state)
!        -----------------

!----------------------------------------------------
! This routine processes the significand state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "significand"
expecting = "digit, exponent-letter, operator-binary, close-paren, end-of-expression"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = SIGNIFICAND

case(LETTER_NOT_EXPONENT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(EXPONENT_LETTER)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = EXPONENT

case(OPERATOR_BINARY_NOT_PLUS_MINUS, PLUS, MINUS)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_no_paren(operator_stack,value_stack, &
                                char_to_operator(char))
   state = START_EXPRESSION

case(DOT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(UNDERSCORE)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPEN_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(CLOSE_PAREN)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_with_paren(operator_stack,value_stack)
   state = END_PAREN_EXPRESSION

case(END_OF_EXPRESSION)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   state = FINISH

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_significand

!        --------------
function state_exponent(expression,exp_ptr,operator_stack,value_stack, &
                        variables,number_string,num_ptr, &
                        variable_string,var_ptr) result(state)
!        --------------

!----------------------------------------------------
! This routine processes the exponent state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "exponent"
expecting = "digit, plus, minus"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = EXPONENT_AFTER_SIGN

case(LETTER_NOT_EXPONENT, EXPONENT_LETTER)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPERATOR_BINARY_NOT_PLUS_MINUS)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(PLUS, MINUS)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = EXPONENT_AFTER_SIGN

case(DOT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(UNDERSCORE)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPEN_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(CLOSE_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(END_OF_EXPRESSION)
   call expression_error(state_string,expression,exp_ptr,"end of expression", &
                         expecting)

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_exponent

!        -------------------------------
function state_exponent_right_after_sign(expression,exp_ptr,operator_stack, &
                                         value_stack,variables,number_string, &
                                         num_ptr,variable_string,var_ptr)  &
                                         result(state)
!        -------------------------------

!----------------------------------------------------
! This routine processes the exponent-right-after-sign state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "exponent-right-after-sign"
expecting = "digit"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = EXPONENT_AFTER_SIGN

case(LETTER_NOT_EXPONENT, EXPONENT_LETTER)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPERATOR_BINARY_NOT_PLUS_MINUS, PLUS, MINUS)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(DOT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(UNDERSCORE)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPEN_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(CLOSE_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(END_OF_EXPRESSION)
   call expression_error(state_string,expression,exp_ptr,"end of expression", &
                         expecting)

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_exponent_right_after_sign

!        -------------------------
function state_exponent_after_sign(expression,exp_ptr,operator_stack, &
                                   value_stack,variables,number_string, &
                                   num_ptr,variable_string,var_ptr) &
                                   result(state)
!        -------------------------

!----------------------------------------------------
! This routine processes the exponent_after_sign state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "exponent_after_sign"
expecting = "digit, operator-binary, close-paren, end-of-expression"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   num_ptr = num_ptr + 1
   number_string(num_ptr:num_ptr) = char
   state = EXPONENT_AFTER_SIGN

case(LETTER_NOT_EXPONENT, EXPONENT_LETTER)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPERATOR_BINARY_NOT_PLUS_MINUS, PLUS, MINUS)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_no_paren(operator_stack,value_stack, &
                                char_to_operator(char))
   state = START_EXPRESSION

case(DOT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(UNDERSCORE)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPEN_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(CLOSE_PAREN)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_with_paren(operator_stack,value_stack)
   state = START_EXPRESSION

case(END_OF_EXPRESSION)
   read(number_string,*) value
   number_string = ""
   num_ptr = 0
   call push_stack(value_stack,value)
   state = FINISH

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_exponent_after_sign

!        --------------
function state_variable(expression,exp_ptr,operator_stack,value_stack, &
                        variables,number_string,num_ptr, &
                        variable_string,var_ptr) result(state)
!        --------------

!----------------------------------------------------
! This routine processes the variable state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "variable"
expecting = "digit, letter, operator-binary, underscore, close-paren, end-of-expression"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   var_ptr = var_ptr + 1
   variable_string(var_ptr:var_ptr) = char
   state = VARIABLE

case(LETTER_NOT_EXPONENT, EXPONENT_LETTER)
   var_ptr = var_ptr + 1
   variable_string(var_ptr:var_ptr) = char
   state = VARIABLE

case(OPERATOR_BINARY_NOT_PLUS_MINUS, PLUS, MINUS)
   value = find_var_name_tree(variable_string(1:var_ptr),variables)
   variable_string = ""
   var_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_no_paren(operator_stack,value_stack, &
                                char_to_operator(char))
   state = START_EXPRESSION

case(DOT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(UNDERSCORE)
   var_ptr = var_ptr + 1
   variable_string(var_ptr:var_ptr) = char
   state = VARIABLE

case(OPEN_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(CLOSE_PAREN)
   value = find_var_name_tree(variable_string(1:var_ptr),variables)
   variable_string = ""
   var_ptr = 0
   call push_stack(value_stack,value)
   call evaluate_stack_with_paren(operator_stack,value_stack)
   state = END_PAREN_EXPRESSION

case(END_OF_EXPRESSION)
   value = find_var_name_tree(variable_string(1:var_ptr),variables)
   variable_string = ""
   var_ptr = 0
   call push_stack(value_stack,value)
   state = FINISH

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_variable

!        ----------------------
function state_end_paren_expression(expression,exp_ptr,operator_stack, &
                                    value_stack,variables,number_string, &
                                    num_ptr,variable_string,var_ptr)  &
                                    result(state)
!        ----------------------

!----------------------------------------------------
! This routine processes the end-paren-expression state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
integer :: state
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: char_type
character(len=80) :: state_string, expecting
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! for error reporting

state_string = "end-paren-expression"
expecting = "operator-binary, close-paren, end-of-expression"

! get next character

call get_char(expression,exp_ptr,char,char_type)

! transition to next state

select case (char_type)

case(DIGIT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(LETTER_NOT_EXPONENT, EXPONENT_LETTER)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPERATOR_BINARY_NOT_PLUS_MINUS, PLUS, MINUS)
   call evaluate_stack_no_paren(operator_stack,value_stack, &
                                char_to_operator(char))
   state = START_EXPRESSION

case(DOT)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(UNDERSCORE)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(OPEN_PAREN)
   call expression_error(state_string,expression,exp_ptr-1,char,expecting)

case(CLOSE_PAREN)
   call evaluate_stack_with_paren(operator_stack,value_stack)
   state = END_PAREN_EXPRESSION

case(END_OF_EXPRESSION)
   state = FINISH

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized state in expression parser")
   stop

end select

end function state_end_paren_expression

!        ------------
function state_finish(expression,exp_ptr,operator_stack, &
                      value_stack,variables,number_string, &
                      num_ptr,variable_string,var_ptr)  &
                      result(value)
!        -----------

!----------------------------------------------------
! This routine processes the end-paren-expression state
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: expression
integer, intent(inout) :: exp_ptr
type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
type(var_name_tree), intent(in) :: variables
character(len=*), intent(inout) :: number_string
integer, intent(inout) :: num_ptr
character(len=*), intent(inout) :: variable_string
integer, intent(inout) :: var_ptr
real(my_real) :: value
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: old_operator
real(my_real) :: value1, value2, new_value
!----------------------------------------------------
! Begin executable code

old_operator = pop_stack(operator_stack)
do while (old_operator /= END_OF_STACK)
   if (old_operator == BEGIN_PAREN) then
      ierr = USER_INPUT_ERROR
      call fatal("unmatched parenthesis in expression in .geo file", &
                 expression)
      stop
   endif
   value1 = pop_stack(value_stack)
   if (old_operator /= UNARY_MINUS) then
      value2 = pop_stack(value_stack)
   endif
   select case (old_operator)
   case(EXPONENTIATION)
      new_value = value2**value1
   case(UNARY_MINUS)
      new_value = -value1
   case(TIMES)
      new_value = value2*value1
   case(DIVIDE)
      new_value = value2/value1
   case(MODULUS)
      new_value = mod(value2,value1)
   case(ADD)
      new_value = value2+value1
   case(SUBTRACT)
      new_value = value2-value1
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("unrecognized operator in state_finish")
      stop
   end select
   call push_stack(value_stack,new_value)
   old_operator = pop_stack(operator_stack)
end do

value = pop_stack(value_stack)
if (value == END_OF_STACK) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("premature end of value_stack in state_finish")
   stop
endif
if (peek_stack(value_stack) /= END_OF_STACK) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("value_stack not empty at end of parsing expression")
   stop
endif

end function state_finish

!          -----------------------
subroutine evaluate_stack_no_paren(operator_stack,value_stack,new_operator)
!          ----------------------

!----------------------------------------------------
! This routine evaluates exp op exp on the stack until the operator on top
! of the stack has lower precidence than new_operator, and puts new_operator
! on top of the stack
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
integer, intent(in) :: new_operator
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: old_operator
real(my_real) :: value1, value2, new_value
!----------------------------------------------------
! Begin executable code

old_operator = peek_stack(operator_stack)
if (old_operator == END_OF_STACK) old_operator = EMPTY_STACK

do while (precidence(old_operator,new_operator) >= 0)
   if (old_operator == EMPTY_STACK .or. &
       old_operator == BEGIN_PAREN) exit
   old_operator = pop_stack(operator_stack)
   value1 = pop_stack(value_stack)
   if (old_operator /= UNARY_MINUS) then
      value2 = pop_stack(value_stack)
   endif
   select case (old_operator)
   case(EXPONENTIATION)
      new_value = value2**value1
   case(UNARY_MINUS)
      new_value = -value1
   case(TIMES)
      new_value = value2*value1
   case(DIVIDE)
      new_value = value2/value1
   case(MODULUS)
      new_value = mod(value2,value1)
   case(ADD)
      new_value = value2+value1
   case(SUBTRACT)
      new_value = value2-value1
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("unrecognized operator in state_finish")
      stop
   end select
   call push_stack(value_stack,new_value)
   old_operator = peek_stack(operator_stack)
   if (old_operator == END_OF_STACK) old_operator = EMPTY_STACK
end do

call push_stack(operator_stack,new_operator)

end subroutine evaluate_stack_no_paren

!          -----------------------
subroutine evaluate_stack_with_paren(operator_stack,value_stack)
!          ----------------------

!----------------------------------------------------
! This routine evaluates exp op exp on the stack until an opening
! parenthesis is found.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(stack_integer), intent(inout) :: operator_stack
type(stack_my_real), intent(inout) :: value_stack
!----------------------------------------------------
! Local variables:

character(len=1) :: char
integer :: old_operator
real(my_real) :: value1, value2, new_value
!----------------------------------------------------
! Begin executable code

old_operator = pop_stack(operator_stack)
if (old_operator == END_OF_STACK) old_operator = EMPTY_STACK

do
   if (old_operator == EMPTY_STACK) then
      ierr = USER_INPUT_ERROR
      call fatal("didn't find matching open parenthesis in expression in .geo file")
      stop
   endif
   if (old_operator == BEGIN_PAREN) exit
   value1 = pop_stack(value_stack)
   if (old_operator /= UNARY_MINUS) then
      value2 = pop_stack(value_stack)
   endif
   select case (old_operator)
   case(EXPONENTIATION)
      new_value = value2**value1
   case(UNARY_MINUS)
      new_value = -value1
   case(TIMES)
      new_value = value2*value1
   case(DIVIDE)
      new_value = value2/value1
   case(MODULUS)
      new_value = mod(value2,value1)
   case(ADD)
      new_value = value2+value1
   case(SUBTRACT)
      new_value = value2-value1
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("unrecognized operator in state_finish")
      stop
   end select
   call push_stack(value_stack,new_value)
   old_operator = pop_stack(operator_stack)
   if (old_operator == END_OF_STACK) old_operator = EMPTY_STACK
end do

end subroutine evaluate_stack_with_paren

!        ----------------
function char_to_operator(char)
!        ----------------

!----------------------------------------------------
! This routine returns the operator corresponding to the given character
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=1), intent(in) :: char
integer :: char_to_operator
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

select case(char)
case ("^")
   char_to_operator = EXPONENTIATION
case ("*")
   char_to_operator = TIMES
case ("/")
   char_to_operator = DIVIDE
case ("%")
   char_to_operator = MODULUS
case ("+")
   char_to_operator = ADD
case ("-")
   char_to_operator = SUBTRACT
case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized character in char_to_operator")
end select

end function char_to_operator

!          ----------------
subroutine expression_error(state_string,expression,ptr,char,expecting)
!          ----------------

!----------------------------------------------------
! This routine reports an error in the expression
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: state_string, expression
integer, intent(in) :: ptr
character(len=*), intent(in) :: char, expecting
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

write(errunit,"(A)")
write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(3A)") "          PHAML Version ",version_number," ERROR"

write(errunit,"(A)")
write(errunit,"(A)") " Error while parsing expression from .geo file."
write(errunit,"(A)")
write(errunit,"(2A)") "expression: ",trim(expression)
write(errunit,"(2A)") "state: ",trim(state_string)
write(errunit,"(A,I4)") "position: ",ptr
write(errunit,"(2A)") "character: ",char
write(errunit,"(2A)") "expecting: ",trim(expecting)
write(errunit,"(A)")

write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(A)")

stop

end subroutine expression_error

end module grid_init_mod
