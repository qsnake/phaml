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

!----------------------------------------------------
! This file contains the program for graphics processes.  It consists of
! a module with the routines that start/terminate graphics windows and
! routines that draw the graphics, followed by a short main program, which
! is actually a subroutine.
!
! RESTRICTION 3D tetrahedra
!----------------------------------------------------

module graphics_mod

!----------------------------------------------------
! Other modules used:

use global
use message_passing
use opengl_gl
use opengl_glu
use opengl_glut
use hash_mod
use view_modifier
use gridtype_mod
use grid_util
use evaluate
use phaml_type_mod
use sort_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! The following types are defined:

!----------------------------------------------------
! The following parameters are defined:

! values for what function to plot on the grid

integer(glcint), parameter :: DRAW_NO_FUNCTION    = 20, &
                              DRAW_SOLUTION       = 21, &
                              DRAW_TRUE           = 22, &
                              DRAW_ERROR          = 23

! values to indicate we want the intersection of a plane with an element

integer(glcint), parameter :: YZ_PLANE = 30, &
                              XZ_PLANE = 31, &
                              XY_PLANE = 32

! values for preprocessing a function to draw

integer(glcint), parameter :: PREPROC_NONE = 0, &
                              PREPROC_ABS  = 1, &
                              PREPROC_LOG  = 2, &
                              PREPROC_SQ   = 3, &
                              PREPROC_NEG  = 4

! values for what color to use

integer(glcint), parameter :: COLOR_TRANSPARENT = 0, &
                              COLOR_WHITE       = 1, &
                              COLOR_BLACK       = 2, &
                              COLOR_SOLUTION    = 3, &
                              COLOR_TRUE        = 4, &
                              COLOR_ERROR       = 5, &
                              COLOR_OWNER       = 6, &
                              COLOR_VERT_OWNER  = 7, &
                              COLOR_SIZE        = 8, &
                              COLOR_PART_BOUND  = 9, &
                              COLOR_DEGREE      = 10

! values for color schemes

integer(glcint), parameter :: SCHEME_RAINBOW        = 0, &
                              SCHEME_GRAY           = 1, &
                              SCHEME_DOUBLE_RAINBOW = 3, &
                              SCHEME_STEP_SEQ       = 4

! values for what label to place on elements and vertices

integer(glcint), parameter :: LABEL_NOLABEL   = 0, &
                              LABEL_LID       = 1, &
                              LABEL_GID       = 2, &
                              LABEL_PARTITION = 3

! values for what to use for color key

integer, parameter :: NO_KEY =      1_glcint, &
                      EDGE_KEY =    2_glcint, &
                      FACE_KEY =    3_glcint, &
                      CUTTING_KEY = 4_glcint

! other parameters

integer, parameter :: REFTREE_ROOT = 0, &
                      NO_OWNER = 0
real (my_real), parameter :: myzero  = 0.0_my_real
real (glfloat), parameter :: glfzero  = 0.0_glfloat, &
                             glfhalf  = 0.5_glfloat, &
                             glfone   = 1.0_glfloat
real(gldouble), parameter :: gldzero = 0.0_gldouble, &
                             gldhalf = 0.5_gldouble, &
                             gldone  = 1.0_gldouble
real(gldouble), parameter :: MAX_DEGREE = 23
type(point), parameter :: NULL_POINT = point(huge(0.0_my_real), &
                                       huge(0.0_my_real), huge(0.0_my_real))

!----------------------------------------------------
! The following variables are defined:

! the grid data structure

type (grid_type), save :: grid

! convenience pointers into the grid data structure

type (element_t), pointer :: element(:)
type (face_t), pointer :: face(:)
type (edge_t), pointer ::    edge(:)
type (vertex_t), pointer ::  vertex(:)
real(my_real), pointer :: vertex_solution(:,:,:), vertex_exact(:,:,:)

integer, allocatable ::         elem_owner(:)
!type(hash_key), allocatable ::  neighbors(:,:)

! graphics choices

integer :: label_elements = LABEL_NOLABEL, &
           label_verts    = LABEL_NOLABEL, &
           label_edges    = LABEL_NOLABEL, &
           label_faces    = LABEL_NOLABEL, &
           color_faces    = COLOR_TRANSPARENT, &
           color_lines    = COLOR_BLACK, &
           color_cutting_plane = COLOR_SOLUTION, &
           color_cutting_plane_lines = COLOR_TRANSPARENT, &
           draw_isosurface= DRAW_NO_FUNCTION, &
           preproc_func   = PREPROC_NONE, &
           color_scheme   = SCHEME_RAINBOW, &
           step_scheme_steps = 4, &
           step_scheme_hues  = 6, &
           color_key_color = NO_KEY, &
           subelement_resolution = 0
logical :: vert_associated_element=.false., edge_associated_element=.false., &
           face_associated_element=.false., cutting_plane_contour=.false., &
           convert_cylindrical_to_cartesian=.false.

! cropping region

real(my_real) :: xcrop1 = -huge(myzero), &
                 xcrop2 =  huge(myzero), &
                 ycrop1 = -huge(myzero), &
                 ycrop2 =  huge(myzero), &
                 zcrop1 = -huge(myzero), &
                 zcrop2 =  huge(myzero)

! contour plot parameters

logical :: contour_values_given  = .false.
integer :: num_contour = 12
real(my_real), allocatable :: actual_contours(:)

! isosurface parameters

integer :: num_isosurface = 12

! parameters for exploding the grid into partitions

type(point), allocatable :: part_cent(:), explode_shift(:),explode_elem_shift(:)
real(gldouble) :: old_explode_factor = 1.0_gldouble, &
                  old_explode_elem_factor = 1.0_gldouble

! which lights are used

logical :: leftlighton=.false., rightlighton=.true., toplighton=.false., &
           bottomlighton=.false., movelighton=.false.

! useful min and max values

real(my_real) :: xmin, xmax, ymin, ymax, zmin, zmax, &
                 minsolut, maxsolut, maxabssolut, mintrue,  maxtrue, &
                 maxabstrue, minerror, maxerror, &
                 minabserr, maxabserr, minsize, maxsize, &
                 maxdomain
real(my_real), allocatable :: all_minsolut(:,:),all_maxsolut(:,:), &
                              all_mintrue (:,:),all_maxtrue (:,:), &
                              all_minerror(:,:),all_maxerror(:,:), &
                              all_maxabserr(:,:)

! components are scaled individually or all the same

logical :: indiv_compnt_scale = .true.

! other parameters for control of drawing the grid

logical :: draw_axes_flag  = .false., &
           draw_sfc_flag   = .false., &
           draw_sfcio_flag = .false., &
           show_key_flag   = .false., &
           grid_newview    = .true.
integer :: eigen          = 1, &
           compnt         = 1
real(glfloat) :: offset   = 1.0_glfloat
real(gldouble) :: edge_transparency = 1.0_gldouble, &
                  face_transparency = 1.0_gldouble, &
                  isosurface_transparency = 0.5_gldouble

! GLUT identifiers

integer(glcint) :: grid_win, grid_menu, eigen_menu=-1, component_menu=-1, &
                   scheme_menu, eigen_submenu, eigen_10s_menus(99)
integer(gluint) :: grid_list=1

! misc variables

integer :: my_processor, nproc, neigen, ncompnt, menu_neigen, menu_ncompnt
character (len=HOSTLEN) :: host
type (proc_info), pointer :: procs
real(gldouble), allocatable :: owner_color(:,:)
real(gldouble) :: lookat_x, lookat_y, lookat_z
logical :: window_initialized = .false.
integer(glcint) :: intzero=0

!----------------------------------------------------
! Non-module procedures used are:

contains

!          ---------------
subroutine process_message
!          ---------------

!----------------------------------------------------
! This routine checks for a message and processes it
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: proc, ni, nr, allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
type(phaml_solution_type) :: phaml_solution
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! see if modview has requested a redraw

if (modview_requests_redraw) then
   call draw_grid
   modview_requests_redraw = .false.
endif

! receive a message if there is one

if (PARALLEL == SEQUENTIAL) then
   call sequential_recv(imess,ni,rmess,nr)
else
   call phaml_recv(procs,proc,imess,ni,rmess,nr,101,noblock=.true.)
endif

! first entry in imess gives the job

if (ni > 0) then
   select case (imess(1))

! initialize graphics

   case (GRAPHICS_INIT)
      global_element_kind = element_kind
      if (global_element_kind /= TETRAHEDRAL_ELEMENT) then
         call fatal("using graphics3D.f90 for tetrahedra with a different kind of element")
         stop
      endif
      nullify(secret_grid)
      allocate(grid%head_level_elem(1),stat=allocstat)
      if (allocstat /= 0) then
         call fatal("allocation failed in process_message",intlist=(/allocstat/))
         stop
      endif
      call hash_table_init(grid%elem_hash,imess(2))
      call hash_table_init(grid%face_hash,imess(5))
      call hash_table_init(grid%edge_hash,imess(3))
      call hash_table_init(grid%vert_hash,imess(4))
      xmin = rmess(1)
      xmax = rmess(2)
      ymin = rmess(3)
      ymax = rmess(4)
      zmin = rmess(5)
      zmax = rmess(6)
      maxdomain = max(max(abs(xmax-xmin),abs(ymax-ymin)),abs(zmax-zmin))
      nullify(grid%element,grid%face,grid%edge,grid%vertex, &
              grid%vertex_solution,grid%vertex_exact)
      deallocate(imess,rmess,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed",intlist=(/allocstat/))
      endif
      call init_windows

! draw grid

   case (GRAPHICS_GRID)
      if (ni > 1) then ! new grid data sent
         call reallocate_grid(imess(2),imess(9),imess(3),imess(4),imess(5), &
                              imess(6))
         call unpack_grid(imess,rmess)
         call set_owner_color
      endif
      call draw_grid
      deallocate(imess,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed",intlist=(/allocstat/))
      endif
      if (nr > 0) then
         deallocate(rmess,stat=allocstat)
         if (allocstat /= 0) then
            call warning("deallocation failed",intlist=(/allocstat/))
         endif
      endif

! terminate

   case (GRAPHICS_TERMINATE)

      call reallocate_grid(0,0,0,0,0,0)

      if (allocated(part_cent)) deallocate(part_cent,stat=allocstat)
      if (allocated(explode_shift)) deallocate(explode_shift,stat=allocstat)
      if (allocated(explode_elem_shift)) deallocate(explode_elem_shift, &
                                                    stat=allocstat)
      call glutdestroywindow(grid_win)
      call terminate_comm(procs,.true.,.true.)
      deallocate(procs,stat=allocstat)
      deallocate(imess,stat=allocstat)
      stop

   case (9) ! update usermod

      phaml_solution%procs = procs
      call update_usermod(phaml_solution)
      if (ni > 0) deallocate(imess)

   case (10) ! increase process universe

      call increase_universe(procs,imess)
      deallocate(imess,stat=allocstat)

   case (11) ! decrease process universe

      call decrease_universe(procs)
      deallocate(imess,stat=allocstat)

   case (12) ! save as postscript
      call write_postscript(intzero)
      deallocate(imess,stat=allocstat)

   end select
endif

return
end subroutine process_message
 
!          ---------------
subroutine reallocate_grid(nelem,nface,nedge,nvert,soln2,soln3)
!          ---------------

!----------------------------------------------------
! This routine reallocates the grid arrays
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: nelem, nface, nedge, nvert, soln2, soln3
!----------------------------------------------------
! Local variables:

integer :: i, allocstat
!----------------------------------------------------
! Begin executable code

! Nullify convenience pointers

nullify(element,face,edge,vertex,vertex_solution,vertex_exact)

! Deallocate memory that is already allocated

if (associated(grid%element)) then
   do i=1,grid%biggest_elem
      if (associated(grid%element(i)%solution)) &
         deallocate(grid%element(i)%solution,stat=allocstat)
      if (associated(grid%element(i)%exact)) &
         deallocate(grid%element(i)%exact,stat=allocstat)
   end do
   deallocate(grid%element,stat=allocstat)
endif
if (associated(grid%face)) then
   do i=1,grid%biggest_face
      if (associated(grid%face(i)%solution)) &
         deallocate(grid%face(i)%solution,stat=allocstat)
      if (associated(grid%face(i)%exact)) &
         deallocate(grid%face(i)%exact,stat=allocstat)
   end do
   deallocate(grid%face,stat=allocstat)
endif
if (associated(grid%edge)) then
   do i=1,grid%biggest_edge
      if (associated(grid%edge(i)%solution)) &
         deallocate(grid%edge(i)%solution,stat=allocstat)
      if (associated(grid%edge(i)%exact)) &
         deallocate(grid%edge(i)%exact,stat=allocstat)
   end do
   deallocate(grid%edge,stat=allocstat)
endif
if (associated(grid%vertex)) deallocate(grid%vertex,stat=allocstat)
if (associated(grid%vertex_solution)) deallocate(grid%vertex_solution)
if (associated(grid%vertex_exact)) deallocate(grid%vertex_exact)
if (allocated(elem_owner)) deallocate(elem_owner,stat=allocstat)
!if (allocated(neighbors)) deallocate(neighbors,stat=allocstat)

! allocate to new sizes, if the new size is not 0

if (nelem > 0) then
   allocate(grid%element(nelem),elem_owner(nelem), &
!            neighbors(NEIGHBORS_PER_ELEMENT,nelem), &
            stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
endif

if (nface > 0) then
   allocate(grid%face(nface),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
endif

if (nedge > 0) then
   allocate(grid%edge(nedge),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
endif

if (nvert > 0) then
   allocate(grid%vertex(nvert),grid%vertex_solution(nvert,soln2,soln3), &
            grid%vertex_exact(nvert,soln2,soln3),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
endif

! nullify components to be allocated later

do i=1,nelem
   nullify(grid%element(i)%solution,grid%element(i)%exact)
end do
do i=1,nface
   nullify(grid%face(i)%solution,grid%face(i)%exact)
end do
do i=1,nedge
   nullify(grid%edge(i)%solution,grid%edge(i)%exact)
end do

! replace hash tables

call hash_table_destroy(grid%elem_hash)
call hash_table_destroy(grid%face_hash)
call hash_table_destroy(grid%edge_hash)
call hash_table_destroy(grid%vert_hash)
if (nelem > 0) call hash_table_init(grid%elem_hash,nelem)
if (nface > 0) call hash_table_init(grid%face_hash,nface)
if (nedge > 0) call hash_table_init(grid%edge_hash,nedge)
if (nvert > 0) call hash_table_init(grid%vert_hash,nvert)

! set convenience pointers into grid data structrue

if (nelem > 0) element => grid%element
if (nface > 0) face => grid%face
if (nedge > 0) edge => grid%edge
if (nvert > 0) then
   vertex => grid%vertex
   vertex_solution => grid%vertex_solution
   vertex_exact => grid%vertex_exact
endif

end subroutine reallocate_grid

!          -----------
subroutine unpack_grid(imess,rmess)
!          -----------

!----------------------------------------------------
! This subroutine unpacks the message containing the grid and flags 
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: imess(:)
real (my_real), intent(in) :: rmess(:)
!----------------------------------------------------
! Local variables:

integer :: ind,rind,elem,vert,i,allocstat,other_lid,part,oldvert,newvert, &
           num_elem,iedge,ssize,iface
real (my_real) :: perim,xcent,ycent,zcent,valcent(1),elem_size
logical :: duplicate
type(hash_key), allocatable :: vert_gid(:,:), edge_gid(:,:), &
!                               invert(:), outvert(:), &
                               face_gid(:,:)
real(my_real), allocatable :: sum_perim(:)
real(my_real) :: rho, phi, small=.005_my_real
integer, pointer :: head_temp(:)
logical, allocatable :: visited_vert(:)

!----------------------------------------------------
! Begin executable code

! initialize min and max values

xmin      =  huge(myzero)
xmax      = -huge(myzero)
ymin      =  huge(myzero)
ymax      = -huge(myzero)
zmin      =  huge(myzero)
zmax      = -huge(myzero)
minsolut  =  huge(myzero)
maxsolut  = -huge(myzero)
mintrue   =  huge(myzero)
maxtrue   = -huge(myzero)
minerror  =  huge(myzero)
maxerror  = -huge(myzero)
minsize   =  huge(myzero)
maxsize   = -huge(myzero)

allocate(vert_gid(VERTICES_PER_ELEMENT,size(element)), &
         edge_gid(EDGES_PER_ELEMENT,size(element)), &
         face_gid(FACES_PER_ELEMENT,size(element)), &
!         invert(size(element)), outvert(size(element)), &
         stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for unpacking grid",intlist=(/allocstat/))
   stop
endif

! unpack the processor info

nproc = imess(7)
my_processor  = imess(8)
ind = 9
host = " "
do i=1,HOSTLEN
   host(i:i) = char(imess(ind+i))
end do
ind = ind + HOSTLEN

! unpack number of eigenfunctions and components (system_size)

neigen = imess(ind+1)
if (eigen > neigen) eigen = neigen
ncompnt = imess(ind+2)
grid%system_size = ncompnt
ind = ind+2

! determine the eigenfunction and component to use

if (compnt > ncompnt) compnt = ncompnt
if (allocated(all_minsolut)) deallocate(all_minsolut,stat=allocstat)
if (allocated(all_maxsolut)) deallocate(all_maxsolut,stat=allocstat)
if (allocated(all_mintrue)) deallocate(all_mintrue,stat=allocstat)
if (allocated(all_maxtrue)) deallocate(all_maxtrue,stat=allocstat)
if (allocated(all_minerror)) deallocate(all_minerror,stat=allocstat)
if (allocated(all_maxerror)) deallocate(all_maxerror,stat=allocstat)
if (allocated(all_maxabserr)) deallocate(all_maxabserr,stat=allocstat)
allocate(all_minsolut(ncompnt,neigen), all_maxsolut(ncompnt,neigen), &
         all_mintrue(ncompnt,neigen), all_maxtrue(ncompnt,neigen), &
         all_minerror(ncompnt,neigen), all_maxerror(ncompnt,neigen), &
         all_maxabserr(ncompnt,neigen), stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed in unpack_grid",intlist=(/allocstat/))
   stop
endif
all_minsolut  =  huge(myzero)
all_mintrue   =  huge(myzero)
all_minerror  =  huge(myzero)
all_maxsolut  = -huge(myzero)
all_maxtrue   = -huge(myzero)
all_maxerror  = -huge(myzero)
all_maxabserr = -huge(myzero)

! unpack the elements

element%level = -999
num_elem = 0
grid%head_level_elem = END_OF_LIST
grid%nlev = 0
rind = 0
grid%biggest_elem = 0
elem = imess(ind+1)
do while (elem /= END_OF_ELEMENTS)
   grid%biggest_elem = max(grid%biggest_elem,elem)
   ind = ind + 1
   ssize = imess(ind+1)
   ind = ind + 1
   element(elem)%gid = hash_unpack_key(imess,ind+1)
   ind = ind + KEY_SIZE
! avoid duplicates from different processors by not adding it to the linked
! list (at end of loop), and hash table, but continue to process the input
   other_lid = hash_decode_key(element(elem)%gid,grid%elem_hash)
   duplicate = (other_lid /= HASH_NOT_FOUND)
   if (.not. duplicate) call hash_insert(element(elem)%gid,elem,grid%elem_hash)
   allocate(element(elem)%solution(ssize/(neigen*ncompnt),ncompnt,neigen), &
            element(elem)%exact(ssize/(neigen*ncompnt),ncompnt,neigen), &
            stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed in unpack_grid",intlist=(/allocstat/))
      stop
   endif
   element(elem)%solution=reshape(rmess(rind+1:rind+ssize), &
                                  (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   element(elem)%exact=reshape(rmess(rind+1:rind+ssize), &
                               (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   do i=1,VERTICES_PER_ELEMENT
      vert_gid(i,elem) = hash_unpack_key(imess,ind+1+(i-1)*KEY_SIZE)
   end do
   ind = ind + VERTICES_PER_ELEMENT*KEY_SIZE
   do i=1,EDGES_PER_ELEMENT
      edge_gid(i,elem) = hash_unpack_key(imess,ind+1+(i-1)*KEY_SIZE)
   end do
   ind = ind + EDGES_PER_ELEMENT*KEY_SIZE
   do i=1,FACES_PER_ELEMENT
      face_gid(i,elem) = hash_unpack_key(imess,ind+1+(i-1)*KEY_SIZE)
   end do
   ind = ind + FACES_PER_ELEMENT*KEY_SIZE
   element(elem)%degree = imess(ind+1)
   ind = ind + 1
   element(elem)%level = imess(ind+1)
   ind = ind + 1
!   invert(elem) = hash_unpack_key(imess,ind+1)
!   ind = ind + KEY_SIZE
!   outvert(elem) = hash_unpack_key(imess,ind+1)
!   ind = ind + KEY_SIZE
!   element(elem)%order = imess(ind+1:ind+MAX_CHILD)
!   ind = ind + MAX_CHILD
   element(elem)%isleaf = (imess(ind+1) == 1)
   ind = ind + 1
   elem_owner(elem) = imess(ind+1)
   ind = ind + 1
!   do i=1,NEIGHBORS_PER_ELEMENT
!      neighbors(i,elem) = hash_unpack_key(imess,ind+1+(i-1)*KEY_SIZE)
!   end do
!   ind = ind + NEIGHBORS_PER_ELEMENT*KEY_SIZE
   grid%nlev = max(grid%nlev,element(elem)%level)
   num_elem = max(num_elem,elem)
   if (element(elem)%level > size(grid%head_level_elem)) then
      head_temp => grid%head_level_elem
      allocate(grid%head_level_elem(element(elem)%level),stat=allocstat)
      if (allocstat /= 0) then
         call fatal("allocation failed in unpack_grid",intlist=(/allocstat/))
         stop
      endif
      grid%head_level_elem = END_OF_LIST
      grid%head_level_elem(1:size(head_temp)) = head_temp
      deallocate(head_temp)
   endif
   if (.not. duplicate) then
      element(elem)%next = grid%head_level_elem(element(elem)%level)
      grid%head_level_elem(element(elem)%level) = elem
   endif
   elem = imess(ind+1)
end do
ind = ind + 1

! unpack the vertices

grid%biggest_vert = 0
vert = imess(ind+1)
do while (vert /= END_OF_VERTICES)
   grid%biggest_vert = max(grid%biggest_vert,vert)
   ind = ind + 1
   vertex(vert)%assoc_elem = imess(ind+1)
   ind = ind + 1
   vertex(vert)%gid = hash_unpack_key(imess,ind+1)
   other_lid = hash_decode_key(vertex(vert)%gid,grid%vert_hash)
! insert in hash table if this is the first instance or I am the owner
   if (other_lid == HASH_NOT_FOUND) then
      call hash_insert(vertex(vert)%gid,vert,grid%vert_hash)
   elseif (elem_owner(vertex(vert)%assoc_elem) == my_processor) then
      call hash_remove(vertex(vert)%gid,grid%vert_hash)
      call hash_insert(vertex(vert)%gid,vert,grid%vert_hash)
   endif
   ind = ind + KEY_SIZE
   vertex(vert)%coord%x = rmess(rind+1)
   vertex(vert)%coord%y = rmess(rind+2)
   vertex(vert)%coord%z = rmess(rind+3)
   rind = rind + 3
   if (convert_cylindrical_to_cartesian) then
      rho = vertex(vert)%coord%x + small
      phi = min(vertex(vert)%coord%y,8*atan(1.0_my_real)-small)
      vertex(vert)%coord%x = rho*cos(phi)
      vertex(vert)%coord%y = rho*sin(phi)
   endif
   vertex_solution(vert,:,:) = reshape(rmess(rind+1:rind+neigen*ncompnt), &
                                      (/ncompnt,neigen/))
   rind = rind + neigen*ncompnt
   vertex_exact(vert,:,:) = reshape(rmess(rind+1:rind+neigen*ncompnt), &
                                    (/ncompnt,neigen/))
   rind = rind + neigen*ncompnt
   xmin     = min(xmin,vertex(vert)%coord%x)
   xmax     = max(xmax,vertex(vert)%coord%x)
   ymin     = min(ymin,vertex(vert)%coord%y)
   ymax     = max(ymax,vertex(vert)%coord%y)
   zmin     = min(zmin,vertex(vert)%coord%z)
   zmax     = max(zmax,vertex(vert)%coord%z)
   all_minsolut = min(all_minsolut,vertex_solution(vert,:,:))
   all_maxsolut = max(all_maxsolut,vertex_solution(vert,:,:))
   all_mintrue = min(all_mintrue,vertex_exact(vert,:,:))
   all_maxtrue = max(all_maxtrue,vertex_exact(vert,:,:))
   all_minerror = min(all_minerror,vertex_solution(vert,:,:)-vertex_exact(vert,:,:))
   all_maxerror = max(all_maxerror,vertex_solution(vert,:,:)-vertex_exact(vert,:,:))
   vert  = imess(ind+1)
end do
ind = ind+1

! unpack the faces

grid%biggest_face = 0
iface = imess(ind+1)
do while (iface /= END_OF_FACES)
   grid%biggest_face = max(grid%biggest_face,iface)
   ind = ind + 1
   ssize = imess(ind+1)
   ind = ind + 1
   face(iface)%assoc_elem = imess(ind+1)
   ind = ind + 1
   face(iface)%gid = hash_unpack_key(imess,ind+1)
   other_lid = hash_decode_key(face(iface)%gid,grid%face_hash)
! insert in hash table if this is the first instance or I am the owner
   if (other_lid == HASH_NOT_FOUND) then
      call hash_insert(face(iface)%gid,iface,grid%face_hash)
   elseif (elem_owner(face(iface)%assoc_elem) == my_processor) then
      call hash_remove(face(iface)%gid,grid%face_hash)
      call hash_insert(face(iface)%gid,iface,grid%face_hash)
   endif
   ind = ind + KEY_SIZE
   face(iface)%vertex(1) = hash_decode_key(hash_unpack_key(imess,ind+1), &
                                           grid%vert_hash)
   if (face(iface)%vertex(1) == HASH_NOT_FOUND) then
      call fatal("didn't find face vertex 1 in vertex hash table")
      stop
   endif
   ind = ind + KEY_SIZE
   face(iface)%vertex(2) = hash_decode_key(hash_unpack_key(imess,ind+1), &
                                           grid%vert_hash)
   if (face(iface)%vertex(2) == HASH_NOT_FOUND) then
      call fatal("didn't find face vertex 2 in vertex hash table")
      stop
   endif
   ind = ind + KEY_SIZE
   face(iface)%vertex(3) = hash_decode_key(hash_unpack_key(imess,ind+1), &
                                           grid%vert_hash)
   if (face(iface)%vertex(3) == HASH_NOT_FOUND) then
      call fatal("didn't find face vertex 3 in vertex hash table")
      stop
   endif
   ind = ind + KEY_SIZE
   face(iface)%degree = imess(ind+1)
   ind = ind + 1
   allocate(face(iface)%solution(ssize/(neigen*ncompnt),ncompnt,neigen), &
           face(iface)%exact(ssize/(neigen*ncompnt),ncompnt,neigen),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed in unpack_grid",intlist=(/allocstat/))
      stop
   endif
   face(iface)%solution=reshape(rmess(rind+1:rind+ssize), &
                                (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   face(iface)%exact=reshape(rmess(rind+1:rind+ssize), &
                             (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   iface = imess(ind+1)
end do
ind = ind+1

! unpack the edges

grid%biggest_edge = 0
iedge = imess(ind+1)
do while (iedge /= END_OF_EDGES)
   grid%biggest_edge = max(grid%biggest_edge,iedge)
   ind = ind + 1
   ssize = imess(ind+1)
   ind = ind + 1
   edge(iedge)%assoc_elem = imess(ind+1)
   ind = ind + 1
   edge(iedge)%gid = hash_unpack_key(imess,ind+1)
   other_lid = hash_decode_key(edge(iedge)%gid,grid%edge_hash)
! insert in hash table if this is the first instance or I am the owner
   if (other_lid == HASH_NOT_FOUND) then
      call hash_insert(edge(iedge)%gid,iedge,grid%edge_hash)
   elseif (elem_owner(edge(iedge)%assoc_elem) == my_processor) then
      call hash_remove(edge(iedge)%gid,grid%edge_hash)
      call hash_insert(edge(iedge)%gid,iedge,grid%edge_hash)
   endif
   ind = ind + KEY_SIZE
   edge(iedge)%vertex(1) = hash_decode_key(hash_unpack_key(imess,ind+1), &
                                           grid%vert_hash)
   if (edge(iedge)%vertex(1) == HASH_NOT_FOUND) then
      call fatal("didn't find edge vertex 1 in vertex hash table")
      stop
   endif
   ind = ind + KEY_SIZE
   edge(iedge)%vertex(2) = hash_decode_key(hash_unpack_key(imess,ind+1), &
                                           grid%vert_hash)
   if (edge(iedge)%vertex(2) == HASH_NOT_FOUND) then
      call fatal("didn't find edge vertex 2 in vertex hash table")
      stop
   endif
   ind = ind + KEY_SIZE
   edge(iedge)%degree = imess(ind+1)
   ind = ind + 1
   allocate(edge(iedge)%solution(ssize/(neigen*ncompnt),ncompnt,neigen), &
           edge(iedge)%exact(ssize/(neigen*ncompnt),ncompnt,neigen),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed in unpack_grid",intlist=(/allocstat/))
      stop
   endif
   edge(iedge)%solution=reshape(rmess(rind+1:rind+ssize), &
                                (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   edge(iedge)%exact=reshape(rmess(rind+1:rind+ssize), &
                             (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   iedge = imess(ind+1)
end do

maxdomain = max(max(abs(xmax-xmin),abs(ymax-ymin)),abs(zmax-zmin))
modview_xmin = xmin
modview_xmax = xmax
modview_ymin = ymin
modview_ymax = ymax
modview_zmin = zmin
modview_zmax = zmax

xcrop1 = max(xcrop1,xmin-.01_my_real)
xcrop2 = min(xcrop2,xmax+.01_my_real)
ycrop1 = max(ycrop1,ymin-.01_my_real)
ycrop2 = min(ycrop2,ymax+.01_my_real)
zcrop1 = max(zcrop1,zmin-.01_my_real)
zcrop2 = min(zcrop2,zmax+.01_my_real)

! set the faces of the elements from the face gids,
! set the edges of the elements from the edge gids,
! set the vertices of the elements from the vertex gids,
! decode the in and out vertices
! find min and max error indicator from owned leaf elements
! find min and max error at element midpoints
! find the min and max element edge size

do elem=1,num_elem
   if (element(elem)%level == -999) cycle
   do i=1,FACES_PER_ELEMENT
      iface = hash_decode_key(face_gid(i,elem),grid%face_hash)
      if (iface == HASH_NOT_FOUND) then
         call warning("face hash code not found; using face 1; element", &
                      intlist=(/elem/))
         call hash_print_key(face_gid(i,elem),errunit)
         element(elem)%face(i) = 1
      else
         element(elem)%face(i) = iface
      endif
   end do
   do i=1,EDGES_PER_ELEMENT
      iedge = hash_decode_key(edge_gid(i,elem),grid%edge_hash)
      if (iedge == HASH_NOT_FOUND) then
         call warning("edge hash code not found; using edge 1; element", &
                      intlist=(/elem/))
         call hash_print_key(edge_gid(i,elem),errunit)
         element(elem)%edge(i) = 1
      else
         element(elem)%edge(i) = iedge
      endif
   end do
   do i=1,VERTICES_PER_ELEMENT
      vert = hash_decode_key(vert_gid(i,elem),grid%vert_hash)
      if (vert == HASH_NOT_FOUND) then
         call warning("vertex hash code not found; using vertex 1; element", &
                      intlist=(/elem/))
         call hash_print_key(vert_gid(i,elem),errunit)
         element(elem)%vertex(i) = 1
      else
         element(elem)%vertex(i) = vert
      endif
   end do
!   element(elem)%in      = hash_decode_key(invert(elem),grid%vert_hash)
!   element(elem)%out     = hash_decode_key(outvert(elem),grid%vert_hash)
   if (element(elem)%isleaf) then
      elem_size = maxval( &
         (/ sqrt((vertex(element(elem)%vertex(1))%coord%x-vertex(element(elem)%vertex(2))%coord%x)**2 + &
                 (vertex(element(elem)%vertex(1))%coord%y-vertex(element(elem)%vertex(2))%coord%y)**2), &
            sqrt((vertex(element(elem)%vertex(2))%coord%x-vertex(element(elem)%vertex(3))%coord%x)**2 + &
                 (vertex(element(elem)%vertex(2))%coord%y-vertex(element(elem)%vertex(3))%coord%y)**2), &
            sqrt((vertex(element(elem)%vertex(3))%coord%x-vertex(element(elem)%vertex(1))%coord%x)**2 + &
                 (vertex(element(elem)%vertex(3))%coord%y-vertex(element(elem)%vertex(1))%coord%y)**2) /) )
      minsize = min(minsize,elem_size)
      maxsize = max(maxsize,elem_size)
   endif
end do
minsize = log(minsize)
maxsize = log(maxsize)

! Set max and min solution, true and error.
! Evaluate errors at element midpoints to avoid 0 error for initial condition.

do elem=1,num_elem
   if (element(elem)%level == -999) cycle
   if (.not. element(elem)%isleaf) cycle
   if (elem_owner(elem) /= my_processor .and.  my_processor/=MASTER) cycle
   xcent = sum(vertex(element(elem)%vertex)%coord%x)/VERTICES_PER_ELEMENT
   ycent = sum(vertex(element(elem)%vertex)%coord%y)/VERTICES_PER_ELEMENT
   zcent = sum(vertex(element(elem)%vertex)%coord%z)/VERTICES_PER_ELEMENT
   call graphics_evaluate_soln((/xcent/),(/ycent/),(/zcent/),elem,valcent)
! TEMP only makes the adjustment for the first eigen,compnt
   valcent(1) = valcent(1) - graphics_trues(xcent,ycent,zcent)
   if (valcent(1) < all_minerror(1,1)) all_minerror(1,1) = valcent(1)
   if (valcent(1) > all_maxerror(1,1)) all_maxerror(1,1) = valcent(1)
end do

all_maxabserr = max(abs(all_minerror),abs(all_maxerror))
call set_max_min

deallocate(vert_gid,edge_gid, &
!           invert,outvert, &
           stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed for vert_gid",intlist=(/allocstat/))
endif

! determine a center of each partition as a weighted average of the
! centers of the elements, using the square of the perimeter of the
! elements as the weights

if (allocated(part_cent)) deallocate(part_cent,stat=allocstat)
if (allocated(explode_shift)) deallocate(explode_shift,stat=allocstat)
if (allocated(explode_elem_shift)) deallocate(explode_elem_shift,stat=allocstat)
allocate(part_cent(nproc),explode_shift(NO_OWNER:nproc), &
         explode_elem_shift(size(element)),stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for explode data structures", &
              intlist=(/allocstat/))
   stop
endif
part_cent = point(myzero,myzero,myzero)
explode_shift = point(myzero,myzero,myzero)
allocate(sum_perim(nproc),stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for sum_perim",intlist=(/allocstat/))
   stop
endif
sum_perim = myzero
do elem=1,num_elem
   if (element(elem)%level == -999) cycle
   if (.not. element(elem)%isleaf) cycle
   part = elem_owner(elem)
   if (part == NO_OWNER) cycle
   oldvert = element(elem)%vertex(VERTICES_PER_ELEMENT)
   perim = myzero; xcent = myzero; ycent = myzero; zcent = myzero
   do i=1,VERTICES_PER_ELEMENT
      newvert = element(elem)%vertex(i)
      perim = perim + sqrt( &
               (vertex(newvert)%coord%x-vertex(oldvert)%coord%x)**2 + &
               (vertex(newvert)%coord%y-vertex(oldvert)%coord%y)**2)
      xcent = xcent + vertex(newvert)%coord%x
      ycent = ycent + vertex(newvert)%coord%y
      zcent = zcent + vertex(newvert)%coord%z
      oldvert = newvert
   end do
   part_cent(part)%x = part_cent(part)%x+perim*perim*xcent/VERTICES_PER_ELEMENT
   part_cent(part)%y = part_cent(part)%y+perim*perim*ycent/VERTICES_PER_ELEMENT
   part_cent(part)%z = part_cent(part)%z+perim*perim*zcent/VERTICES_PER_ELEMENT
   sum_perim(part) = sum_perim(part) + perim*perim
end do
where(sum_perim /= myzero) part_cent%x = part_cent%x / sum_perim
where(sum_perim /= myzero) part_cent%y = part_cent%y / sum_perim
where(sum_perim /= myzero) part_cent%z = part_cent%z / sum_perim
deallocate(sum_perim,stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed for sum_perim",intlist=(/allocstat/))
endif

! change the menu if the number of eigenfunctions has changed

if (menu_neigen /= neigen) then
   call update_eigen_menu
endif
if (menu_ncompnt /= ncompnt+2) then
   call update_compnt_menu
endif

grid_newview = .true.

return
end subroutine unpack_grid

!          -----------
subroutine set_max_min()
!          -----------

!----------------------------------------------------
! This routine set max and min values for the current eigen and compnt
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (allocated(all_minsolut)) then

   select case(compnt)

   case (-1)
      minsolut = 0.0_my_real
      maxsolut = sum(abs(all_maxsolut(:,eigen)))
      maxabssolut = maxsolut
      if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
      mintrue = 0.0_my_real
      maxtrue = sum(abs(all_maxtrue(:,eigen)))
      maxabstrue = maxtrue
      if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
      minerror = 0.0_my_real
      maxerror = sum(abs(all_maxerror(:,eigen)))
      maxabserr = sum(abs(all_maxabserr(:,eigen)))
      if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real

   case (-2)
      minsolut = 0.0_my_real
      maxsolut = sum(all_maxsolut(:,eigen)**2)
      maxabssolut = maxsolut
      if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
      mintrue = 0.0_my_real
      maxtrue = sum(all_maxtrue(:,eigen)**2)
      maxabstrue = maxtrue
      if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
      minerror = 0.0_my_real
      maxerror = sum(all_maxerror(:,eigen)**2)
      maxabserr = sum(all_maxabserr(:,eigen)**2)
      if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real

   case default
      if (indiv_compnt_scale) then
         minsolut = all_minsolut(compnt,eigen)
         maxsolut = all_maxsolut(compnt,eigen)
         maxabssolut = max(abs(minsolut),abs(maxsolut))
         if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
         mintrue = all_mintrue(compnt,eigen)
         maxtrue = all_maxtrue(compnt,eigen)
         maxabstrue = max(abs(mintrue),abs(maxtrue))
         if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
         minerror = all_minerror(compnt,eigen)
         maxerror = all_maxerror(compnt,eigen)
         maxabserr = all_maxabserr(compnt,eigen)
         if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real
      else
         minsolut = minval(all_minsolut(:,eigen))
         maxsolut = maxval(all_maxsolut(:,eigen))
         maxabssolut = max(abs(minsolut),abs(maxsolut))
         if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
         mintrue = minval(all_mintrue(:,eigen))
         maxtrue = maxval(all_maxtrue(:,eigen))
         maxabstrue = max(abs(mintrue),abs(maxtrue))
         if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
         minerror = minval(all_minerror(:,eigen))
         maxerror = maxval(all_maxerror(:,eigen))
         maxabserr = maxval(all_maxabserr(:,eigen))
         if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real
      endif

   end select
endif

end subroutine set_max_min

!          ----------------
subroutine update_eigen_menu
!          ----------------

!----------------------------------------------------
! This routine changes the number of entries in the eigenfunction selection
! menu when the number of eigenfunctions changes
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer(glcint) :: hold
!----------------------------------------------------
! Begin executable code

if (eigen_menu == -1) return ! no eigenfunction menu

hold = glutgetmenu()
call glutsetmenu(eigen_menu)

if (menu_neigen > neigen) then
   call decrease_eigen_menu(menu_neigen,neigen)
elseif (menu_neigen < neigen) then
   call increase_eigen_menu(menu_neigen,neigen)
endif

menu_neigen = neigen
call glutsetmenu(hold)
end subroutine update_eigen_menu

!          ------------------
subroutine increase_eigen_menu(oldn,newn)
!          ------------------

!----------------------------------------------------
! This routine increases the number of entries in the eigenfunction selection menu
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments:

integer, intent(in) :: oldn, newn

!----------------------------------------------------
! Local variables:

integer(glcint) :: i, j, oldone, newone
character(len=17) :: str
!----------------------------------------------------
! Begin executable code

! add new entries in primary menu

do i=oldn+1,min(newn,9)
   write(str,"(a14,i3)") "eigenfunction ",i
   call glutaddmenuentry(trim(str),i)
end do

! create submenu for selection of 10's if needed and not already there

if (oldn < 10 .and. newn >= 10) then
   eigen_submenu = glutcreatemenu(select_eigen)
   call glutsetmenu(eigen_menu)
   call glutaddsubmenu("more",eigen_submenu)
endif

! create menu for each of the 10's that are needed

do i=1+oldn/10,newn/10
   eigen_10s_menus(i) = glutcreatemenu(select_eigen)
end do
do i=1+oldn/10,newn/10
   call glutsetmenu(eigen_submenu)
   str = "  0's"
   write(str(1:2),"(i2)") i
   call glutaddsubmenu(trim(str),eigen_10s_menus(i))
end do

! create the menu entries in the 10's menus

oldone = oldn - 10*(oldn/10)
newone = newn - 10*(newn/10)

! the old last 10's menu

if (oldn >= 10) then
   i = oldn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j = oldone+1, min(9,newn-10*i)
      write(str,"(a14,i3)") "eigenfunction ",10*i+j
      call glutaddmenuentry(trim(str),10*i+j)
   end do
endif

! menus before the new last menu

do i=1+oldn/10,newn/10-1
   call glutsetmenu(eigen_10s_menus(i))
   do j=0,9
      write(str,"(a14,i3)") "eigenfunction ",10*i+j
      call glutaddmenuentry(trim(str),10*i+j)
   end do
end do

! new last menu

if (oldn/10 /= newn/10) then
   i = newn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j=0,newone
      write(str,"(a14,i3)") "eigenfunction ",10*i+j
      call glutaddmenuentry(trim(str),10*i+j)
   end do
endif

end subroutine increase_eigen_menu

!          ------------------
subroutine decrease_eigen_menu(oldn,newn)
!          ------------------

!----------------------------------------------------
! This routine decreases the number of entries in the eigenfunction selection menu
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments:

integer, intent(in) :: oldn, newn

!----------------------------------------------------
! Local variables:

integer(glcint) :: i, j, oldone, newone
!----------------------------------------------------
! Begin executable code

oldone = oldn - 10*(oldn/10)
newone = newn - 10*(newn/10)

! remove last menu if not needed

if (oldn >= 10 .and. oldn/10 /= newn/10) then
   i = oldn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j=oldone,0,-1
      call glutremovemenuitem(j+1)
   end do
   call glutsetmenu(eigen_submenu)
   call glutremovemenuitem(i)
   call glutsetmenu(eigen_submenu)
endif

! remove menus between the old and new last menu

do i=oldn/10-1,newn/10+1,-1
   call glutsetmenu(eigen_10s_menus(i))
   do j=9,0,-1
      call glutremovemenuitem(j+1)
   end do
   call glutsetmenu(eigen_submenu)
   call glutremovemenuitem(i)
   call glutdestroymenu(eigen_10s_menus(i))
end do

! remove menu items in the new last menu

if (newn >= 10) then
   i = newn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j = min(9,oldn-10*i), newone+1, -1
      call glutremovemenuitem(j+1)
   end do
endif

! remove submenu for selection of 10's if there and no longer needed

if (oldn >= 10 .and. newn < 10) then
   call glutsetmenu(eigen_menu)
   call glutremovemenuitem(10)
   call glutdestroymenu(eigen_submenu)
endif

! remove entries in primary menu

do i=min(9,oldn),newn+1,-1
   call glutremovemenuitem(i)
end do

end subroutine decrease_eigen_menu

!          ----------------
subroutine update_compnt_menu
!          ----------------

!----------------------------------------------------
! This routine changes the number of entries in the component selection
! menu when the number of components changes
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer(glcint) :: i, hold
character(len=15) :: str
!----------------------------------------------------
! Begin executable code

if (component_menu == -1) return ! no component menu

hold = glutgetmenu()
call glutsetmenu(component_menu)

str = ""
if (menu_ncompnt > ncompnt+2) then
   do i=menu_ncompnt,ncompnt+3,-1
      call glutremovemenuitem(i-2)
   end do

elseif (menu_ncompnt < ncompnt+2) then
   do i=menu_ncompnt+1,ncompnt+2
      write(str,"(a10,i3)") "component ",i-2
      call glutaddmenuentry(trim(str),i-2)
   end do
endif

menu_ncompnt = ncompnt+2
call glutsetmenu(hold)
return
end subroutine update_compnt_menu

!          ---------------
subroutine set_owner_color
!          ---------------

!----------------------------------------------------
! This routine defines the color to use for coloring by owner
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

integer :: proc, allocstat
real(gldouble) :: max_val
!----------------------------------------------------
! Begin executable code

if (allocated(owner_color)) then
   deallocate(owner_color,stat=allocstat)
   if (allocstat /= 0) then
      call warning("deallocation failed",intlist=(/allocstat/))
   endif
endif
allocate(owner_color(4,max(2,nproc)),stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for owner colors", &
              intlist=(/allocstat/))
   stop
endif

if (color_scheme == SCHEME_STEP_SEQ) then
   max_val = step_scheme_steps*step_scheme_hues
else
   max_val = max(2,nproc)
endif
do proc=1,max(2,nproc)
   call get_rainbow(real(proc,gldouble),1.0_gldouble,max_val,owner_color(:,proc))
end do

return
end subroutine set_owner_color

!          -----------------
subroutine preproc_and_scale(selector,val,ppval,ppssval,ppmin,ppssmin,ppmax, &
                             ppssmax)
!          -----------------

!----------------------------------------------------
! This routine converts val to the preprocessed value ppval and the
! preprocessed, scaled and shifted value ppssval
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: selector
real(my_real), intent(in) :: val
real(my_real), intent(out) :: ppval, ppssval, ppmin, ppssmin, ppmax, ppssmax
!----------------------------------------------------
! Local variables:

real(my_real) :: minval, maxval
!----------------------------------------------------
! Begin executable code

select case (preproc_func)

case (PREPROC_NONE)
   ppval = val
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      ppmin = minsolut
      ppmax = maxsolut
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmax
      else
         ppssval = (ppval-minsolut)/(maxsolut-minsolut)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      ppmin = mintrue
      ppmax = maxtrue
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmax
      else
         ppssval = (ppval-mintrue)/(maxtrue-mintrue)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      ppmin = minerror
      ppmax = maxerror
      if ((maxerror-minerror) == 0.0_my_real) then
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmax
      else
         ppssval = (ppval-minerror)/(maxerror-minerror)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   end select

case (PREPROC_ABS)
   ppval = abs(val)
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppmin = abs(minsolut)
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minsolut >= 0.0_my_real) then
            ppmin = minsolut
            ppmax = maxsolut
            ppssval = (ppval-minsolut)/(maxsolut-minsolut)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxsolut <= 0.0_my_real) then
            ppmin = -maxsolut
            ppmax = -minsolut
            ppssval = (ppval+maxsolut)/(maxsolut-minsolut)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(abs(minsolut),abs(maxsolut))
            ppssval = ppval/max(abs(minsolut),abs(maxsolut))
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppmin = abs(mintrue)
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (mintrue >= 0.0_my_real) then
            ppmin = mintrue
            ppmax = maxtrue
            ppssval = (ppval-mintrue)/(maxtrue-mintrue)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxsolut <= 0.0_my_real) then
            ppmin = -maxtrue
            ppmax = -mintrue
            ppssval = (ppval+maxtrue)/(maxtrue-mintrue)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(abs(mintrue),abs(maxtrue))
            ppssval = ppval/max(abs(mintrue),abs(maxtrue))
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      if ((maxerror-minerror) == 0.0_my_real) then
         ppmin = abs(minerror)
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minerror >= 0.0_my_real) then
            ppmin = minerror
            ppmax = maxerror
            ppssval = (ppval-minerror)/(maxerror-minerror)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxerror <= 0.0_my_real) then
            ppmin = -maxerror
            ppmax = -minerror
            ppssval = (ppval+maxerror)/(maxerror-minerror)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(abs(minerror),abs(maxerror))
            ppssval = ppval/max(abs(minerror),abs(maxerror))
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   end select

case (PREPROC_NEG)
   ppval = -val
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppmin = -minsolut
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         ppmin = -maxsolut
         ppmax = -minsolut
         ppssval = (ppval+maxsolut)/(maxsolut-minsolut)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppmin = -mintrue
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         ppmin = -maxtrue
         ppmax = -mintrue
         ppssval = (ppval+maxtrue)/(maxtrue-mintrue)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      if ((maxerror-minerror) == 0.0_my_real) then
         ppmin = -minerror
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         ppmin = -maxerror
         ppmax = -minerror
         ppssval = (ppval+maxerror)/(maxerror-minerror)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   end select

case (PREPROC_SQ)
   ppval = val**2
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppmin = minsolut**2
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minsolut > 0.0_my_real) then
            ppmin = minsolut**2
            ppmax = maxsolut**2
            ppssval = (ppval-minsolut**2)/(maxsolut**2-minsolut**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxsolut < 0.0_my_real) then
            ppmin = maxsolut**2
            ppmax = minsolut**2
            ppssval = (ppval-maxsolut**2)/(maxsolut**2-minsolut**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(maxsolut**2,minsolut**2)
            ppssval = ppval/max(maxsolut**2,minsolut**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppmin = mintrue**2
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (mintrue > 0.0_my_real) then
            ppmin = mintrue**2
            ppmax = maxtrue**2
            ppssval = (ppval-mintrue**2)/(maxtrue**2-mintrue**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxtrue < 0.0_my_real) then
            ppmin = maxtrue**2
            ppmax = mintrue**2
            ppssval = (ppval-maxtrue**2)/(maxtrue**2-mintrue**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(maxtrue**2,mintrue**2)
            ppssval = ppval/max(maxtrue**2,mintrue**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      if ((maxerror-minerror) == 0.0_my_real) then
         ppmin = minerror**2
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minerror > 0.0_my_real) then
            ppmin = minerror**2
            ppmax = maxerror**2
            ppssval = (ppval-minerror**2)/(maxerror**2-minerror**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxerror < 0.0_my_real) then
            ppmin = maxerror**2
            ppmax = minerror**2
            ppssval = (ppval-maxerror**2)/(maxerror**2-minerror**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(maxerror**2,minerror**2)
            ppssval = ppval/max(maxerror**2,minerror**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   end select

case (PREPROC_LOG)
   if (abs(val) < 1.0e-14_my_real) then
      ppval = -14.0_my_real
   else
      ppval = log10(abs(val))
   endif
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppval = 0.0_my_real
      minval = 1.0_my_real
      maxval = 1.0_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      minval = minsolut
      maxval = maxsolut
   case (DRAW_TRUE, COLOR_TRUE)
      minval = mintrue
      maxval = maxtrue
   case (DRAW_ERROR, COLOR_ERROR)
      minval = minerror
      maxval = maxerror
   end select
   if (minval < 0.0_my_real) then
      if (maxval > 0.0_my_real) then
         ppmin = -14
         ppmax = max(log10(abs(minval)),log10(abs(maxval)))
      else
         if (maxval > -1.0e-14_my_real) then
            ppmin = -14
            if (minval > -1.0e-14_my_real) then
               ppmax = -14
            else
               ppmax = log10(abs(minval))
            endif
         else
            ppmin = log10(abs(maxval))
            ppmax = log10(abs(minval))
         endif
      endif
   elseif (minval < 1.0e-14_my_real) then
      ppmin = -14
      if (maxval < 1.0e-14_my_real) then
         ppmax = -14
      else
         ppmax = log10(abs(maxval))
      endif
   else
      ppmin = log10(abs(minval))
      ppmax = log10(abs(maxval))
   endif
   if (ppmax == ppmin) then
      ppssval = ppval
      ppssmin = ppmin
      ppssmax = ppmax
   else
      ppssval = (ppval-ppmin)/(ppmax-ppmin)
      ppssmin = 0.0_my_real
      ppssmax = 1.0_my_real
   endif

end select

ppssval = 2*maxdomain*(ppssval-0.5_my_real)
ppssmin = 2*maxdomain*(ppssmin-0.5_my_real)
ppssmax = 2*maxdomain*(ppssmax-0.5_my_real)

end subroutine preproc_and_scale

!          ---------
subroutine draw_grid
!          ---------

!----------------------------------------------------
! This subroutine draws the grid.  Optionally, the elements, faces, and/or
! vertices can be labeled, the associated element can be marked, and the
! faces, vertices, and/or element edges can be colored.
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(glfloat) :: fblack(4) = (/glfzero,glfzero,glfzero,glfone/)
logical :: nolight
character(len=14+HOSTLEN) :: title
integer :: i, allocstat
logical(small_logical), allocatable :: labeled(:)

!----------------------------------------------------
! Begin executable code

! reset the graphics window

call reset_view
if (grid_newview) call recalcview

call glutsetwindow(grid_win)
write(title(12:13),'(i2)') my_processor
title(1:11)='PHAML grid '
title(14:14)=' '
title(15:14+HOSTLEN)=host
call glutSetWindowTitle(title)
call gldeletelists(grid_list, 1_glsizei)
call glnewlist(grid_list, gl_compile_and_execute)

! set lighting

call glDisable(GL_BLEND)
call glEnable(GL_DEPTH_TEST)
call gldisable(gl_lighting)
nolight = .true.

! amount to shift each partition for exploding by partition

do i=1,nproc
   explode_shift(i)%x = (part_cent(i)%x-lookat_x)*(explode_factor-gldone)
   explode_shift(i)%y = (part_cent(i)%y-lookat_y)*(explode_factor-gldone)
   explode_shift(i)%z = (part_cent(i)%z-lookat_z)*(explode_factor-gldone)
end do

! amount to shift each element for exploding by element

do i=1,size(element)
   if (element(i)%level == -999) cycle
   if (.not. element(i)%isleaf) cycle
   explode_elem_shift(i)%x = &
      (sum(vertex(element(i)%vertex)%coord%x)/VERTICES_PER_ELEMENT-lookat_x) * &
      (explode_elem_factor-gldone)
   explode_elem_shift(i)%y = &
      (sum(vertex(element(i)%vertex)%coord%y)/VERTICES_PER_ELEMENT-lookat_y) * &
      (explode_elem_factor-gldone)
   explode_elem_shift(i)%z = &
      (sum(vertex(element(i)%vertex)%coord%z)/VERTICES_PER_ELEMENT-lookat_z) * &
      (explode_elem_factor-gldone)
end do

! recursively draw the cutting planes

if (yz_cutting_plane_value /= CUTTING_PLANE_OFF .or. &
    xz_cutting_plane_value /= CUTTING_PLANE_OFF .or. &
    xy_cutting_plane_value /= CUTTING_PLANE_OFF) then

   call glenable(GL_POLYGON_OFFSET_FILL)
   call glpolygonmode(gl_front_and_back, gl_fill)
   call glbegin(gl_triangles)

   call draw_cutting_planes(grid,nolight)

   call glend
   call gldisable(GL_POLYGON_OFFSET_FILL)

   if (color_cutting_plane_lines /= COLOR_TRANSPARENT) then

      call glpolygonmode(gl_front_and_back, gl_line)
      call glbegin(gl_lines)

      call draw_cutting_plane_lines(grid,nolight)

      call glend

   endif

   if (cutting_plane_contour) then
      call glpolygonmode(gl_front_and_back, gl_line)
      call glbegin(gl_lines)

      call draw_cutting_plane_contours(nolight)

      call glend

   endif

endif

! recursively label elements

if (label_elements /= LABEL_NOLABEL .or. label_verts /= LABEL_NOLABEL .or. &
    label_faces /= LABEL_NOLABEL .or. &
    label_edges /= LABEL_NOLABEL .or. vert_associated_element .or. &
    edge_associated_element .or. face_associated_element) then

   allocate(labeled(size(vertex)),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for labeled",intlist=(/allocstat/))
      stop
   endif

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   labeled = .false.

   call draw_subgrid_label(REFTREE_ROOT,labeled,nolight)

   deallocate(labeled,stat=allocstat)
   if (allocstat /= 0) then
      call warning("deallocation failed for labeled",intlist=(/allocstat/))
   endif

endif

! recursively draw space filling curve

if (draw_sfc_flag) then

   call glpolygonmode(gl_front_and_back, gl_line)
   call glbegin(gl_line_strip)

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   call draw_sfc(REFTREE_ROOT)

   call glend
endif

! label in and out vertices

if (draw_sfcio_flag) then

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   call draw_sfcio(REFTREE_ROOT,1.0_gldouble)

endif

! draw axes

if (draw_axes_flag) then

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   call draw_axes

endif

! for edges, faces and isosurfaces, draw opaque objects before transparent ones

! recursively draw the element edges, if opaque

if (color_lines /= COLOR_TRANSPARENT .and. edge_transparency==1.0_my_real) then

   call glpolygonmode(gl_front_and_back, gl_line)
   call glbegin(gl_lines)

   call draw_subgrid_lines(REFTREE_ROOT,nolight)

   call glend

endif

! recursively draw the element faces, if opaque

if (color_faces /= COLOR_TRANSPARENT .and. face_transparency==1.0_my_real) then

   call glenable(GL_POLYGON_OFFSET_FILL)
   call glpolygonmode(gl_front_and_back, gl_fill)
   call glbegin(gl_triangles)

   call draw_faces(nolight)

   call glend
   call gldisable(GL_POLYGON_OFFSET_FILL)

endif

! recursively draw isosurfaces, if opaque.

if (draw_isosurface /= DRAW_NO_FUNCTION .and. &
    isosurface_transparency == 1.0_my_real) then

   call glpolygonmode(gl_front_and_back, gl_fill)
   call glbegin(gl_triangles)

   call draw_isosurfaces(nolight)

   call glend

endif

! recursively draw the element edges, if transparent

if (color_lines /= COLOR_TRANSPARENT .and. edge_transparency/=1.0_my_real) then

   call glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
   call glEnable(GL_BLEND)
   call glpolygonmode(gl_front_and_back, gl_line)
   call glbegin(gl_lines)

   call draw_subgrid_lines(REFTREE_ROOT,nolight)

   call glend
   call glDisable(GL_BLEND)

endif

! recursively draw the element faces, if transparent

if (color_faces /= COLOR_TRANSPARENT .and. face_transparency/=1.0_my_real) then

   call glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
   call glEnable(GL_BLEND)
   call glenable(GL_POLYGON_OFFSET_FILL)
   call glpolygonmode(gl_front_and_back, gl_fill)
   call glbegin(gl_triangles)

   call draw_faces(nolight)

   call glend
   call gldisable(GL_POLYGON_OFFSET_FILL)
   call glDisable(GL_BLEND)

endif

! recursively draw isosurfaces, if transparent.

if (draw_isosurface /= DRAW_NO_FUNCTION .and. &
    isosurface_transparency /= 1.0_my_real) then

   call glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
   call glEnable(GL_BLEND)
   call glpolygonmode(gl_front_and_back, gl_fill)
   call glbegin(gl_triangles)

   call draw_isosurfaces(nolight)

   call glend

   call glDisable(GL_BLEND)

endif

! draw the color key

if (show_key_flag) then
   call draw_key
endif

call glendlist
call glutpostredisplay

return
end subroutine draw_grid

!                    ------------------
recursive subroutine draw_subgrid_lines(elem,nolight)
!                    ------------------

!----------------------------------------------------
! This routine draws the element edges for the subtree below element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
logical, intent(in) :: nolight
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

logical :: part_bound(NEIGHBORS_PER_ELEMENT)
integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(my_real) :: x(VERTICES_PER_ELEMENT), y(VERTICES_PER_ELEMENT), &
                 z(VERTICES_PER_ELEMENT)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call draw_subgrid_lines(i,nolight)
      i = element(i)%next
   end do
   return
endif

! draw this element only if it is a leaf or we are drawing all levels

if (element(elem)%isleaf) then

! draw this element only if one of the vertices is inside the croping range

 if (inside_crop(element(elem)%vertex)) then

! draw the element edges

   if (color_lines /= COLOR_TRANSPARENT) then

! set the vertices into x,y,z arrays

      x = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x
      y = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y
      z = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%z

! determine if each edge is on a partition boundary

!      do i=1,NEIGHBORS_PER_ELEMENT
!         if (neighbors(i,elem) == BOUNDARY) then
!            part_bound(i) = .true.
!         elseif (elem_owner(hash_decode_key(neighbors(i,elem),grid%elem_hash)) &
!                 /= elem_owner(elem)) then
!            part_bound(i) = .true.
!         else
!            part_bound(i) = .false.
!         endif
!      end do

! draw the element

      call draw_element_lines(grid,x,y,z,elem, &
  elem_owner(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%assoc_elem), &
  elem_owner(edge(element(elem)%edge(1:EDGES_PER_ELEMENT))%assoc_elem), &
                              part_bound,nolight)

   endif
 endif
endif

! draw the elements in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call draw_subgrid_lines(children(i),nolight)
   endif
end do

return
end subroutine draw_subgrid_lines

!                    ------------------
recursive subroutine draw_element_lines(grid,xmr,ymr,zmr,elem,vown,eown, &
                                        part_bound,nolight)
!                    ------------------

!----------------------------------------------------
! This routine draws part of the edges of element elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: xmr(:),ymr(:),zmr(:)
integer, intent(in) :: elem,vown(:),eown(:)
logical, intent(in) :: part_bound(:),nolight
!----------------------------------------------------
! Local variables:

logical :: cpart_bound(NEIGHBORS_PER_ELEMENT)
integer :: i, j, k, cvown(VERTICES_PER_ELEMENT), ceown(EDGES_PER_ELEMENT)
real(my_real) :: xmid(EDGES_PER_ELEMENT),ymid(EDGES_PER_ELEMENT), &
                 zmid(EDGES_PER_ELEMENT), &
                 val1(1),ppval,ppssval,ppmin,ppssmin,ppmax,ppssmax
real(gldouble) :: x(VERTICES_PER_ELEMENT),y(VERTICES_PER_ELEMENT), &
                  z(VERTICES_PER_ELEMENT), normal(3), &
                  val2,elem_size,max_val
real(gldouble) :: color(VERTICES_PER_ELEMENT,4), ecolor(EDGES_PER_ELEMENT,4), &
                 grey(4)  = (/gldhalf,gldhalf,gldhalf,gldone/), &
                 white(4) = (/gldone, gldone, gldone, gldone/), &
                 black(4) = (/gldzero,gldzero,gldzero,gldone/)
real(glfloat) :: fcolor(EDGES_PER_ELEMENT,4)
integer :: edge_verts(EDGES_PER_ELEMENT,2)

!----------------------------------------------------
! Begin executable code

   call get_edge_verts(grid,elem,edge_verts)

! (xmr,ymr,zmr) are used for evaluating solution; (x,y,z) are used as vertices

   x = xmr + explode_shift(elem_owner(elem))%x + explode_elem_shift(elem)%x
   y = ymr + explode_shift(elem_owner(elem))%y + explode_elem_shift(elem)%y
   z = zmr + explode_shift(elem_owner(elem))%z + explode_elem_shift(elem)%z

! if coloring lines by size, determine the length of the longest side

   if (color_lines == COLOR_SIZE) then
      elem_size = maxval( &
         (/ sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2 + (z(1)-z(2))**2), &
            sqrt((x(2)-x(3))**2 + (y(2)-y(3))**2 + (z(2)-z(3))**2), &
            sqrt((x(3)-x(1))**2 + (y(3)-y(1))**2 + (z(3)-z(1))**2), &
            sqrt((x(4)-x(1))**2 + (y(4)-y(1))**2 + (z(4)-z(1))**2), &
            sqrt((x(4)-x(2))**2 + (y(4)-y(2))**2 + (z(4)-z(2))**2), &
            sqrt((x(4)-x(3))**2 + (y(4)-y(3))**2 + (z(4)-z(3))**2) /) )
      elem_size = log(elem_size)
   endif

! find the color of each vertex (in some cases, color of edge)

   do i=1,EDGES_PER_ELEMENT

      select case (color_lines)
      case (COLOR_VERT_OWNER)
         if (i > VERTICES_PER_ELEMENT) exit
         if (vown(i) /= NO_OWNER) then
            color(i,:) = owner_color(:,vown(i))
         else
            color(i,:) = grey
         endif
         color(i,4) = edge_transparency

      case (COLOR_BLACK, COLOR_PART_BOUND)
         ecolor(i,:) = black
         ecolor(i,4) = edge_transparency

      case (COLOR_SOLUTION)
         if (i > VERTICES_PER_ELEMENT) exit
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),zmr(i:i),elem,val1)
         call preproc_and_scale(color_lines,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble), &
                          color(i,:))
         color(i,4) = edge_transparency

      case (COLOR_TRUE)
         if (i > VERTICES_PER_ELEMENT) exit
         val1(1) = graphics_trues(xmr(i),ymr(i),zmr(i))
         call preproc_and_scale(color_lines,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble), &
                          color(i,:))
         color(i,4) = edge_transparency

      case (COLOR_ERROR)
         if (i > VERTICES_PER_ELEMENT) exit
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),zmr(i:i),elem,val1)
         val1(1) = val1(1) - graphics_trues(xmr(i),ymr(i),zmr(i))
         call preproc_and_scale(color_lines,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble), &
                          color(i,:))
         color(i,4) = edge_transparency

      case (COLOR_SIZE)
         call get_rainbow(elem_size,real(minsize,gldouble), &
                          real(maxsize,gldouble),ecolor(i,:))
         ecolor(i,4) = edge_transparency

      case (COLOR_OWNER)
         if (eown(i) /= NO_OWNER) then
            ecolor(i,:) = owner_color(:,eown(i))
         else
            ecolor(i,:) = grey
         endif
         ecolor(i,4) = edge_transparency

      end select

   end do

   if (color_lines == COLOR_DEGREE) then
      if (color_scheme == SCHEME_STEP_SEQ) then
         max_val = step_scheme_steps*step_scheme_hues
      else
         max_val = MAX_DEGREE
      endif
      do i=1,EDGES_PER_ELEMENT
         call get_rainbow(real(edge(element(elem)%edge(i))%degree,gldouble), &
                          1.0_gldouble,max_val,ecolor(i,:))
         ecolor(i,4) = edge_transparency
      end do
   endif

! draw each edge of the tetrahedron

   do i=1,EDGES_PER_ELEMENT

! add the vertex to the display list

      select case (color_lines)

! when coloring by edge owner, edge degree, size or black, draw edge i
! with color i

      case (COLOR_OWNER, COLOR_DEGREE, COLOR_BLACK, COLOR_SIZE)

         if (nolight) then
            call glcolor4dv(ecolor(i,:))
         else
            fcolor = ecolor
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fcolor(i,:))
         endif

         do k=1,2
            j = edge_verts(i,k)
            call glvertex3d(x(j),y(j),z(j))
         end do

! for drawing partition boundary, draw opposite line if opposite element is
! boundary or owned by a different processor

      case (COLOR_PART_BOUND)

         if (part_bound(i)) then

            if (nolight) then
               call glcolor4dv(ecolor(i,:))
            else
               fcolor = ecolor
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fcolor(i,:))
            endif

            do k=1,2
               j = edge_verts(i,k)
               call glvertex3d(x(j),y(j),z(j))
            end do

         endif

! for most cases, draw the edge with colors from the vertices

      case default

         do k=1,2
            j = edge_verts(i,k)
            if (nolight) then
               call glcolor4dv(color(j,:))
            else
               fcolor(j,:) = color(j,:)
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fcolor(j,:))
            endif
! TEMP do I need a normal here? If so, what are the vertices?  And there are
!      other places above where I would also need it.
!            normal = normcrossprod(xmr,ymr,zmr)
!            call glnormal3dv(normal)
            call glvertex3d(x(j),y(j),z(j))
         end do

      end select

   end do

end subroutine draw_element_lines

!                    ------------------
recursive subroutine draw_subgrid_label(elem,labeled,nolight)
!                    ------------------

!----------------------------------------------------
! This routine labels elements, vertices, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
logical(small_logical), intent(inout) :: labeled(:)
logical, intent(in) :: nolight
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: j,i,int_key(KEY_SIZE)
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(gldouble) :: x(VERTICES_PER_ELEMENT), y(VERTICES_PER_ELEMENT), &
                  z(VERTICES_PER_ELEMENT)
real(gldouble) :: xx,yy,zz,h,xmid,ymid,zmid
real, parameter :: delta=.010
real(gldouble) :: elem_lc(4) = (/gldzero,gldzero,gldone,gldone/), &
                 edge_lc(4) = (/gldzero,0.5_gldouble,gldzero,gldone/), &
                 vert_lc(4) = (/gldone,gldzero,gldzero,gldone/), &
                 face_lc(4) = (/gldzero,0.75_gldouble,0.75_gldouble,gldone/), &
                 black(4)   = (/gldzero,gldzero,gldzero,gldone/)
real(glfloat) :: felem_lc(4) = (/glfzero,glfzero,glfone,glfone/), &
                 fedge_lc(4) = (/glfzero,glfone,glfzero,glfone/), &
                 fvert_lc(4) = (/glfone,glfzero,glfzero,glfone/), &
                 fface_lc(4) = (/glfone,glfzero,glfzero,glfone/), &
                 fblack(4)   = (/glfzero,glfzero,glfzero,glfone/)
integer :: edge_verts(EDGES_PER_ELEMENT,2), &
           face_verts(FACES_PER_ELEMENT,VERTICES_PER_FACE)

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call draw_subgrid_label(i,labeled,nolight)
      i = element(i)%next
   end do
   return
endif

! label this element only if it is a leaf or we are drawing all levels

if (element(elem)%isleaf) then

! label this element only if one of the vertices is inside the croping range

 if (inside_crop(element(elem)%vertex)) then

   call get_edge_verts(grid,elem,edge_verts)
   call get_face_verts(grid,elem,face_verts)

! number the element and vertices

   if (label_elements /= LABEL_NOLABEL .or. label_verts /= LABEL_NOLABEL .or. &
       label_faces /= LABEL_NOLABEL .or. &
       label_edges /= LABEL_NOLABEL .or. vert_associated_element .or. &
       edge_associated_element .or. face_associated_element) then

! determine the location of the label

      x = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x + &
          explode_shift(elem_owner(elem))%x + explode_elem_shift(elem)%x
      y = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y + &
          explode_shift(elem_owner(elem))%y + explode_elem_shift(elem)%y
      z = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%z + &
          explode_shift(elem_owner(elem))%z + explode_elem_shift(elem)%z
      h = sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2) + &
          sqrt((x(2)-x(3))**2 + (y(2)-y(3))**2) + &
          sqrt((x(3)-x(1))**2 + (y(3)-y(1))**2)
      xx = sum(x)/size(x)
      yy = sum(y)/size(y)
      zz = sum(z)/size(z)

! element label

      select case (label_elements)
      case (LABEL_NOLABEL)
      case (LABEL_LID)
         if (nolight) then
            call glcolor3d(elem_lc(1),elem_lc(2),elem_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,felem_lc)
         endif
         call number(elem,xx,yy,zz,h)
      case (LABEL_GID)
         if (nolight) then
            call glcolor3d(elem_lc(1),elem_lc(2),elem_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,felem_lc)
         endif
         call hash_pack_key(element(elem)%gid,int_key,1)
         call number(int_key(1),xx,yy,zz,h)
      end select

! vertex label

      select case (label_verts)
      case (LABEL_NOLABEL)
      case (LABEL_LID)
         if (nolight) then
            call glcolor3d(vert_lc(1),vert_lc(2),vert_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fvert_lc)
         endif
         do j=1,VERTICES_PER_ELEMENT
           if (.not. labeled(element(elem)%vertex(j)) .or. &
               explode_elem_factor /= 1.0_glfloat) then
            xx = x(j)
            yy = y(j)
            zz = z(j)
            call number(element(elem)%vertex(j),xx,yy,zz,h)
            labeled(element(elem)%vertex(j)) = .true.
           endif
         end do
      case (LABEL_GID)
         if (nolight) then
            call glcolor3d(vert_lc(1),vert_lc(2),vert_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fvert_lc)
         endif
         do j=1,VERTICES_PER_ELEMENT
           if (.not. labeled(element(elem)%vertex(j)) .or. &
               explode_elem_factor /= 1.0_glfloat) then
            xx = x(j)
            yy = y(j)
            zz = z(j)
            call hash_pack_key(vertex(element(elem)%vertex(j))%gid,int_key,1)
            call number(int_key(1),xx,yy,zz,h)
            labeled(element(elem)%vertex(j)) = .true.
           endif
         end do
      end select

! draw a line segment from the vertices into the associated element

      if (vert_associated_element) then
         if (nolight) then
            call glcolor3d(gldzero,gldzero,gldzero)
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
         endif
         xx = sum(x)/size(x)
         yy = sum(y)/size(y)
         zz = sum(z)/size(z)
         do j=1,VERTICES_PER_ELEMENT
            if (vertex(element(elem)%vertex(j))%assoc_elem == elem) then
               call gllinewidth(2.0_glfloat)
               call glbegin(gl_lines)
               call glVertex3d(.25_gldouble*xx+.75_gldouble*x(j), &
                               .25_gldouble*yy+.75_gldouble*y(j), &
                               .25_gldouble*zz+.75_gldouble*z(j))
               call glVertex3d(x(j),y(j),z(j))
               call glend
               call gllinewidth(glfone)
            endif
         end do
      endif

! draw a line segment from the edges into the associated element

      if (edge_associated_element) then
         if (nolight) then
            call glcolor3d(gldzero,gldzero,gldzero)
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
         endif
         xmid = sum(x)/size(x)
         ymid = sum(y)/size(y)
         zmid = sum(z)/size(z)
         do j=1,EDGES_PER_ELEMENT
            if (edge(element(elem)%edge(j))%assoc_elem == elem) then
               xx = (x(edge_verts(j,1)) + x(edge_verts(j,2)))/2
               yy = (y(edge_verts(j,1)) + y(edge_verts(j,2)))/2
               zz = (z(edge_verts(j,1)) + z(edge_verts(j,2)))/2
               call gllinewidth(2.0_glfloat)
               call glbegin(gl_lines)
               call glVertex3d(.75_gldouble*xx+.25_gldouble*xmid, &
                               .75_gldouble*yy+.25_gldouble*ymid, &
                               .75_gldouble*zz+.25_gldouble*zmid)
               call glVertex3d(xx,yy,zz)
               call glend
               call gllinewidth(glfone)
            endif
         end do
      endif

! draw a line segment from the faces into the associated element

      if (face_associated_element) then
         if (nolight) then
            call glcolor3d(gldzero,gldzero,gldzero)
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
         endif
         xmid = sum(x)/size(x)
         ymid = sum(y)/size(y)
         zmid = sum(z)/size(z)
         do j=1,FACES_PER_ELEMENT
            if (face(element(elem)%face(j))%assoc_elem == elem) then
               xx = sum(x(face_verts(j,:)))/VERTICES_PER_FACE
               yy = sum(y(face_verts(j,:)))/VERTICES_PER_FACE
               zz = sum(z(face_verts(j,:)))/VERTICES_PER_FACE
               call gllinewidth(2.0_glfloat)
               call glbegin(gl_lines)
               call glVertex3d(.75_gldouble*xx+.25_gldouble*xmid, &
                               .75_gldouble*yy+.25_gldouble*ymid, &
                               .75_gldouble*zz+.25_gldouble*zmid)
               call glVertex3d(xx,yy,zz)
               call glend
               call gllinewidth(glfone)
            endif
         end do
      endif

! edge label

      select case (label_edges)
      case (LABEL_NOLABEL)
      case (LABEL_LID)
         if (nolight) then
            call glcolor3d(edge_lc(1),edge_lc(2),edge_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fedge_lc)
         endif
         do j=1,EDGES_PER_ELEMENT
            xx = (x(edge_verts(j,1)) + x(edge_verts(j,2)))/2
            yy = (y(edge_verts(j,1)) + y(edge_verts(j,2)))/2
            zz = (z(edge_verts(j,1)) + z(edge_verts(j,2)))/2
            h = 2.*sqrt((x(edge_verts(j,1)) - x(edge_verts(j,2)))**2 + &
                        (y(edge_verts(j,1)) - y(edge_verts(j,2)))**2 + &
                        (z(edge_verts(j,1)) - z(edge_verts(j,2)))**2)
            call number(element(elem)%edge(j),xx,yy,zz,h)
         end do
      case (LABEL_GID)
         if (nolight) then
            call glcolor3d(edge_lc(1),edge_lc(2),edge_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fedge_lc)
         endif
         do j=1,EDGES_PER_ELEMENT
            xx = (x(edge_verts(j,1)) + x(edge_verts(j,2)))/2
            yy = (y(edge_verts(j,1)) + y(edge_verts(j,2)))/2
            zz = (z(edge_verts(j,1)) + z(edge_verts(j,2)))/2
            h = 2.*sqrt((x(edge_verts(j,1)) - x(edge_verts(j,2)))**2 + &
                        (y(edge_verts(j,1)) - y(edge_verts(j,2)))**2 + &
                        (z(edge_verts(j,1)) - z(edge_verts(j,2)))**2)
            call hash_pack_key(edge(element(elem)%edge(j))%gid,int_key,1)
            call number(int_key(1),xx,yy,zz,h)
         end do
      end select

! face label

      select case (label_faces)
      case (LABEL_NOLABEL)
      case (LABEL_LID)
         call glcolor3d(face_lc(1),face_lc(2),face_lc(3))
         if (nolight) then
            call glcolor3d(face_lc(1),face_lc(2),face_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fface_lc)
         endif
         do j=1,FACES_PER_ELEMENT
            xx = sum(x(face_verts(j,:)))/VERTICES_PER_FACE
            yy = sum(y(face_verts(j,:)))/VERTICES_PER_FACE
            zz = sum(z(face_verts(j,:)))/VERTICES_PER_FACE
            h = 2.*sqrt((x(face_verts(j,1)) - x(face_verts(j,2)))**2 + &
                        (y(face_verts(j,1)) - y(face_verts(j,2)))**2 + &
                        (z(face_verts(j,1)) - z(face_verts(j,2)))**2)
            call number(element(elem)%face(j),xx,yy,zz,h)
         end do
      case (LABEL_GID)
         if (nolight) then
            call glcolor3d(face_lc(1),face_lc(2),face_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fface_lc)
         endif
         do j=1,FACES_PER_ELEMENT
            xx = sum(x(face_verts(j,:)))/VERTICES_PER_FACE
            yy = sum(y(face_verts(j,:)))/VERTICES_PER_FACE
            zz = sum(z(face_verts(j,:)))/VERTICES_PER_FACE
            h = 2.*sqrt((x(face_verts(j,1)) - x(face_verts(j,2)))**2 + &
                        (y(face_verts(j,1)) - y(face_verts(j,2)))**2 + &
                        (z(face_verts(j,1)) - z(face_verts(j,2)))**2)
            call hash_pack_key(face(element(elem)%face(j))%gid,int_key,1)
            call number(int_key(1),xx,yy,zz,h)
         end do
      end select

   endif

 endif
endif

! draw the elements in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call draw_subgrid_label(children(i),labeled,nolight)
   endif
end do

return
end subroutine draw_subgrid_label

!          ----------
subroutine draw_faces(nolight)
!          ----------

!----------------------------------------------------
! This routine draws the element faces
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: element_list(size(element)), iperm(size(element))
integer :: nelem,ier,i,j,elem
real(my_real) :: element_dist(size(element))
real(glDouble) :: m(0:15)
type(point) :: camera
!----------------------------------------------------
! Begin executable code

! determine the location of the camera

call glGetDoublev(GL_MODELVIEW_MATRIX,m)
camera%x = -(m(0)*m(12) + m(1)*m(13) + m(2)*m(14))
camera%y = -(m(4)*m(12) + m(5)*m(13) + m(6)*m(14))
camera%z = -(m(8)*m(12) + m(9)*m(13) + m(10)*m(14))

! make a list of all the leaf elements along with the distance from the
! center of the element to the camera, and sort them from back to front

nelem = 0
call dist_elements(REFTREE_ROOT,element_list,element_dist,nelem,camera)
call sort(element_dist,nelem,iperm,-1,ier)

! for each element, in the sorted order, draw the faces

do i=1,nelem
   elem = element_list(iperm(i))
   call draw_element_faces(elem,nolight)
end do

end subroutine draw_faces

!          ------------------
subroutine draw_element_faces(elem,nolight)
!          ------------------

!----------------------------------------------------
! This routine draws the element faces for element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
logical, intent(in) :: nolight
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(my_real) :: x(3), y(3), z(3)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! draw this element only if one of the vertices is inside the croping range

if (inside_crop(element(elem)%vertex)) then

! for each face, set the vertices into x,y,z arrays and draw the face

   x = vertex(element(elem)%vertex((/1,2,3/)))%coord%x
   y = vertex(element(elem)%vertex((/1,2,3/)))%coord%y
   z = vertex(element(elem)%vertex((/1,2,3/)))%coord%z

   call draw_triangle_solid(grid,x,y,z,elem,nolight,color_faces,0)

   x = vertex(element(elem)%vertex((/1,2,4/)))%coord%x
   y = vertex(element(elem)%vertex((/1,2,4/)))%coord%y
   z = vertex(element(elem)%vertex((/1,2,4/)))%coord%z

   call draw_triangle_solid(grid,x,y,z,elem,nolight,color_faces,0)

   x = vertex(element(elem)%vertex((/1,3,4/)))%coord%x
   y = vertex(element(elem)%vertex((/1,3,4/)))%coord%y
   z = vertex(element(elem)%vertex((/1,3,4/)))%coord%z

   call draw_triangle_solid(grid,x,y,z,elem,nolight,color_faces,0)

   x = vertex(element(elem)%vertex((/2,3,4/)))%coord%x
   y = vertex(element(elem)%vertex((/2,3,4/)))%coord%y
   z = vertex(element(elem)%vertex((/2,3,4/)))%coord%z

   call draw_triangle_solid(grid,x,y,z,elem,nolight,color_faces,0)

endif

end subroutine draw_element_faces

!                    -------------------
recursive subroutine draw_triangle_solid(grid,xmr,ymr,zmr,elem,nolight, &
                                         use_color,recur_lev)
!                    -------------------

!----------------------------------------------------
! This routine draws the triangle with vertices (xmr,ymr,zmr) solid.
! If recur_lev is smaller than subelement_resolution, it replaces the triangle
! with four triangles by connecting the edge midpoints.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: xmr(:),ymr(:),zmr(:)
integer, intent(in) :: elem, recur_lev, use_color
logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: xmid(3),ymid(3),zmid(3), &
                 val1(1),ppval,ppssval,ppmin,ppssmin,ppmax,ppssmax
real(gldouble) :: x(VERTICES_PER_ELEMENT),y(VERTICES_PER_ELEMENT), &
                  z(VERTICES_PER_ELEMENT), &
                  normal(3),val2,elem_size,max_val
real(gldouble) :: color(4), &
                  grey(4)  = (/gldhalf,gldhalf,gldhalf,gldone/), &
                  white(4) = (/gldone, gldone, gldone, gldone/)
real(glfloat) :: fcolor(4), &
                 fgrey(4)  = (/glfhalf,glfhalf,glfhalf,glfone/), &
                 fwhite(4) = (/glfone, glfone, glfone, glfone/)
!----------------------------------------------------
! Begin executable code

if (recur_lev < subelement_resolution) then

! draw four children of this triangle for better resolution

! define midpoints of the sides

   xmid(1) = (xmr(2)+xmr(3))/2
   ymid(1) = (ymr(2)+ymr(3))/2
   zmid(1) = (zmr(2)+zmr(3))/2
   xmid(2) = (xmr(3)+xmr(1))/2
   ymid(2) = (ymr(3)+ymr(1))/2
   zmid(2) = (zmr(3)+zmr(1))/2
   xmid(3) = (xmr(1)+xmr(2))/2
   ymid(3) = (ymr(1)+ymr(2))/2
   zmid(3) = (zmr(1)+zmr(2))/2

! draw the children

   call draw_triangle_solid(grid, &
                           (/xmr(1),xmid(3),xmid(2)/), &
                           (/ymr(1),ymid(3),ymid(2)/), &
                           (/zmr(1),zmid(3),zmid(2)/), &
                           elem,nolight,use_color,recur_lev+1)
   call draw_triangle_solid(grid, &
                           (/xmr(2),xmid(1),xmid(3)/), &
                           (/ymr(2),ymid(1),ymid(3)/), &
                           (/zmr(2),zmid(1),zmid(3)/), &
                           elem,nolight,use_color,recur_lev+1)
   call draw_triangle_solid(grid, &
                           (/xmr(3),xmid(2),xmid(1)/), &
                           (/ymr(3),ymid(2),ymid(1)/), &
                           (/zmr(3),zmid(2),zmid(1)/), &
                           elem,nolight,use_color,recur_lev+1)
   call draw_triangle_solid(grid, &
                           (/xmid(1),xmid(2),xmid(3)/), &
                           (/ymid(1),ymid(2),ymid(3)/), &
                           (/zmid(1),zmid(2),zmid(3)/), &
                           elem,nolight,use_color,recur_lev+1)

else ! draw the triangle

! (xmr,ymr,zmr) are used for evaluating solution; (x,y,z) are used as vertices

   x = xmr + explode_shift(elem_owner(elem))%x + explode_elem_shift(elem)%x
   y = ymr + explode_shift(elem_owner(elem))%y + explode_elem_shift(elem)%y
   z = zmr + explode_shift(elem_owner(elem))%z + explode_elem_shift(elem)%z

! if coloring elements by size, determine the length of the longest side
! of the parent element

   if (use_color == COLOR_SIZE) then
      elem_size = maxval( &
         (/ sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2 + (z(1)-z(2))**2), &
            sqrt((x(2)-x(3))**2 + (y(2)-y(3))**2 + (z(2)-z(3))**2), &
            sqrt((x(3)-x(1))**2 + (y(3)-y(1))**2 + (z(3)-z(1))**2) /) )
      elem_size = elem_size*2**subelement_resolution
      elem_size = log(elem_size)
   endif

   normal = normcrossprod(xmr,ymr,zmr)
   call glnormal3dv(normal)

! find the color to draw the triangle

   do i=1,3
      select case (use_color)
      case (COLOR_OWNER)
         if (elem_owner(elem) /= NO_OWNER) then
            if(nolight) then
               color = owner_color(:,elem_owner(elem))
               color(4) = face_transparency
               call glcolor4dv(color)
            else
               fcolor = owner_color(:,elem_owner(elem))
               fcolor(4) = face_transparency
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                                 fcolor)
            endif
         else
            if (nolight) then
               color = grey
               color(4) = face_transparency
               call glcolor4dv(color)
            else
               fcolor = fgrey
               fcolor(4) = face_transparency
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                                    fcolor)
            endif
         endif

      case (COLOR_WHITE)
         if (nolight) then
            call glcolor3d(gldone,gldone,gldone)
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fwhite)
         endif

      case (COLOR_SOLUTION)
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),zmr(i:i),elem,val1)
         call preproc_and_scale(use_color,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble),color)
         color(4) = face_transparency
         if (nolight) then
            call glcolor4dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_TRUE)
         val1(1) = graphics_trues(xmr(i),ymr(i),zmr(i))
         call preproc_and_scale(use_color,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble),color)
         color(4) = face_transparency
         if (nolight) then
            call glcolor4dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_SIZE)
         call get_rainbow(elem_size,real(minsize,gldouble), &
                          real(maxsize,gldouble),color)
         color(4) = face_transparency
         if (nolight) then
            call glcolor4dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_DEGREE)
         if (color_scheme == SCHEME_STEP_SEQ) then
            max_val = step_scheme_steps*step_scheme_hues
         else
            max_val = MAX_DEGREE
         endif
         call get_rainbow(real(element(elem)%degree,gldouble),1.0_gldouble, &
                          max_val,color)
         color(4) = face_transparency
         if (nolight) then
            call glcolor4dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_ERROR)
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),zmr(i:i),elem,val1)
         val1(1) = val1(1) - graphics_trues(xmr(i),ymr(i),zmr(i))
         call preproc_and_scale(use_color,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble),color)
         color(4) = face_transparency
         if (nolight) then
            call glcolor4dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      end select

! add this vertex to the drawing of the triangle

      call glvertex3d(x(i),y(i),z(i))

   end do ! next vertex

endif ! check for recursion

end subroutine draw_triangle_solid

!          ----------------
subroutine draw_isosurfaces(nolight)
!          ----------------

!----------------------------------------------------
! This routine draws the isosurfaces
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: element_list(size(element)), iperm(size(element))
integer :: nelem,ier,i,j,elem
real(my_real) :: element_dist(size(element)),x(VERTICES_PER_ELEMENT), &
                 y(VERTICES_PER_ELEMENT),z(VERTICES_PER_ELEMENT)
real(my_real) :: isosurface_value(num_isosurface), &
                 isosurface_color(4,num_isosurface)
real(my_real) :: minz, maxz, pval, psval, psmin, psmax
real(glDouble) :: m(0:15)
type(point) :: camera
!----------------------------------------------------
! Begin executable code

! determine the contour values and set the colors

call preproc_and_scale(draw_isosurface,1.0_my_real,pval,psval,minz,psmin, &
                       maxz,psmax)
do i=1,num_isosurface
  isosurface_value(i) = minz+(maxz-minz)*(i-1)/real(num_isosurface-1,my_real)
  call get_rainbow(isosurface_value(i),minz,maxz,isosurface_color(:,i))
end do
isosurface_color(4,:) = isosurface_transparency

! determine the location of the camera

call glGetDoublev(GL_MODELVIEW_MATRIX,m)
camera%x = -(m(0)*m(12) + m(1)*m(13) + m(2)*m(14))
camera%y = -(m(4)*m(12) + m(5)*m(13) + m(6)*m(14))
camera%z = -(m(8)*m(12) + m(9)*m(13) + m(10)*m(14))

! make a list of all the leaf elements along with the distance from the
! center of the element to the camera, and sort them from back to front

nelem = 0
call dist_elements(REFTREE_ROOT,element_list,element_dist,nelem,camera)
call sort(element_dist,nelem,iperm,-1,ier)

! for each element, in the sorted order, draw the isosurfaces

do i=1,nelem
   elem = element_list(iperm(i))
   do j=1,num_isosurface
      x = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x
      y = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y
      z = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%z
      call draw_element_isosurface(grid,x,y,z,isosurface_value, &
                                   isosurface_color,elem,nolight,0,j)
   end do
end do

end subroutine draw_isosurfaces

!                    -------------
recursive subroutine dist_elements(elem,element_list,element_dist,nelem,camera)
!                    -------------

!----------------------------------------------------
! This routine makes a list of the leaf elements and the distance from the
! center of each element to the camera
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
integer, intent(inout) :: element_list(:)
real(my_real), intent(inout) :: element_dist(:)
integer, intent(inout) :: nelem
type(point), intent(in) :: camera
!----------------------------------------------------
! Local variables:

integer :: i
type(point) :: center
integer :: allc(MAX_CHILD), children(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call dist_elements(i,element_list,element_dist,nelem,camera)
      i = element(i)%next
   end do
   return
endif

! if this is a leaf, add it to the list and compute the distance to the camera

if (element(elem)%isleaf) then
   nelem = nelem + 1
   element_list(nelem) = elem
   center = vertex(element(elem)%vertex(1))%coord + &
            vertex(element(elem)%vertex(2))%coord + &
            vertex(element(elem)%vertex(3))%coord + &
            vertex(element(elem)%vertex(4))%coord
   center = center/4.0_my_real
   element_dist(nelem) = dot_point(center-camera,center-camera)

! if this is not a leaf, do the children below it

else

   allc = ALL_CHILDREN
   children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call dist_elements(children(i),element_list,element_dist,nelem,camera)
      endif
   end do

endif

end subroutine dist_elements

!                    -----------------------
recursive subroutine draw_element_isosurface(grid,xmr,ymr,zmr,isosurface_value,&
                                             isosurface_color,elem,nolight, &
                                             recur_lev,surf)
!                    -----------------------

!----------------------------------------------------
! This routine draws the isosurfaces in an element with vertices (xmr,ymr,zmr).
! TEMP Does recur_lev do anything here? and do I need recursive?
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: xmr(:), ymr(:), zmr(:), isosurface_value(:), &
                             isosurface_color(:,:)
integer, intent(in) :: elem, recur_lev, surf
logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: i, j, i1, i2
real(glfloat) :: fcolor(4)
type(point) :: verts(VERTICES_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

verts = func_intersect_tet(grid,draw_isosurface,isosurface_value(surf),xmr, &
                           ymr,zmr,elem)

! if there are three intersections, draw the triangle formed by them

if (.not. verts(3) == NULL_POINT) then

   verts(1) = verts(1) + explode_shift(elem_owner(elem)) + &
                         explode_elem_shift(elem)
   verts(2) = verts(2) + explode_shift(elem_owner(elem)) + &
                         explode_elem_shift(elem)
   verts(3) = verts(3) + explode_shift(elem_owner(elem)) + &
                         explode_elem_shift(elem)

   if (nolight) then
      call glcolor4dv(isosurface_color(:,surf))
   else
      fcolor = isosurface_color(:,surf)
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                        fcolor)
   endif
   do j=1,3
      call glvertex3d(verts(j)%x,verts(j)%y,verts(j)%z)
   end do

! if there are four intersections, also draw triangle(1,3,4)

   if (.not. verts(4) == NULL_POINT) then

      verts(4) = verts(4) + explode_shift(elem_owner(elem)) + &
                            explode_elem_shift(elem)
      do j=1,4
         if (j==2) cycle
         call glvertex3d(verts(j)%x,verts(j)%y,verts(j)%z)
      end do
   endif

endif

end subroutine draw_element_isosurface

!          -------------------
subroutine draw_cutting_planes(grid,nolight)
!          -------------------

!----------------------------------------------------
! This routine draws the enabled cutting planes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: plane, elem, i
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! do all three planes

do i=1,3

! pick the plane and see if it is active

   select case (i)
   case (1)
      if (yz_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = YZ_PLANE
      value = yz_cutting_plane_value
   case (2)
      if (xz_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = XZ_PLANE
      value = xz_cutting_plane_value
   case (3)
      if (xy_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = XY_PLANE
      value = xy_cutting_plane_value
   end select

! call the worker routine with each element of the initial grid

   elem = grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      call recur_draw_cutting_plane(grid,nolight,elem,plane,value)
      elem = element(elem)%next
   end do

end do

end subroutine draw_cutting_planes

!                    ------------------------
recursive subroutine recur_draw_cutting_plane(grid,nolight,elem,plane,value)
!                    ------------------------

!----------------------------------------------------
! This routine does the work of drawing the enabled cutting planes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical, intent(in) :: nolight
integer, intent(in) :: elem, plane
real(my_real), intent(in) :: value
!----------------------------------------------------
! Local variables:

integer :: allc(MAX_CHILD), children(MAX_CHILD), i
type(point) :: verts(4)
!----------------------------------------------------
! Begin executable code

! only draw the plane through leaf elements

if (element(elem)%isleaf) then

! get the intersection of the cutting plane and element

   verts = func_intersect_tet(grid,plane,value, &
                              vertex(element(elem)%vertex)%coord%x, &
                              vertex(element(elem)%vertex)%coord%y, &
                              vertex(element(elem)%vertex)%coord%z, &
                              elem)

! draw the triangle(s) if an intersection was found

   if (.not. verts(3) == NULL_POINT) then
      call draw_triangle_solid(grid, &
                               (/verts(1)%x,verts(2)%x,verts(3)%x/), &
                               (/verts(1)%y,verts(2)%y,verts(3)%y/), &
                               (/verts(1)%z,verts(2)%z,verts(3)%z/), &
                               elem,nolight,color_cutting_plane,0)
   endif
   if (.not. verts(4) == NULL_POINT) then
      call draw_triangle_solid(grid, &
                               (/verts(3)%x,verts(4)%x,verts(1)%x/), &
                               (/verts(3)%y,verts(4)%y,verts(1)%y/), &
                               (/verts(3)%z,verts(4)%z,verts(1)%z/), &
                               elem,nolight,color_cutting_plane,0)
   endif

else ! not a leaf, do the children

   allc = ALL_CHILDREN
   children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call recur_draw_cutting_plane(grid,nolight,children(i),plane,value)
      endif
   end do

endif

end subroutine recur_draw_cutting_plane

!          ------------------------
subroutine draw_cutting_plane_lines(grid,nolight)
!          ------------------------

!----------------------------------------------------
! This routine draws grid lines on the enabled cutting planes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: plane, elem, i
real(my_real) :: value
!----------------------------------------------------
! Begin executable code

! do all three planes

do i=1,3

! pick the plane and see if it is active

   select case (i)
   case (1)
      if (yz_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = YZ_PLANE
      value = yz_cutting_plane_value
   case (2)
      if (xz_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = XZ_PLANE
      value = xz_cutting_plane_value
   case (3)
      if (xy_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = XY_PLANE
      value = xy_cutting_plane_value
   end select

! call the worker routine with each element of the initial grid

   elem = grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      call recur_draw_cutting_plane_lines(grid,nolight,elem,plane,value)
      elem = element(elem)%next
   end do

end do

end subroutine draw_cutting_plane_lines

!                    ------------------------------
recursive subroutine recur_draw_cutting_plane_lines(grid,nolight,elem,plane, &
                                                    value)
!                    ------------------------------

!----------------------------------------------------
! This routine does the work of drawing the grid lines on the enabled cutting
! planes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical, intent(in) :: nolight
integer, intent(in) :: elem, plane
real(my_real), intent(in) :: value
!----------------------------------------------------
! Local variables:

integer :: allc(MAX_CHILD), children(MAX_CHILD), i
type(point) :: verts(4)
!----------------------------------------------------
! Begin executable code

! only draw the plane through leaf elements

if (element(elem)%isleaf) then

! get the intersection of the cutting plane and element

   verts = func_intersect_tet(grid,plane,value, &
                              vertex(element(elem)%vertex)%coord%x, &
                              vertex(element(elem)%vertex)%coord%y, &
                              vertex(element(elem)%vertex)%coord%z, &
                              elem)

! draw the edges if an intersection was found
! TEMP hard coded black

   if (.not. verts(3) == NULL_POINT) then
      call glcolor4dv((/gldzero,gldzero,gldzero,gldone/))
      verts(1) = verts(1) + explode_shift(elem_owner(elem)) + &
                            explode_elem_shift(elem)
      verts(2) = verts(2) + explode_shift(elem_owner(elem)) + &
                            explode_elem_shift(elem)
      verts(3) = verts(3) + explode_shift(elem_owner(elem)) + &
                            explode_elem_shift(elem)
      call glvertex3d(verts(1)%x,verts(1)%y,verts(1)%z)
      call glvertex3d(verts(2)%x,verts(2)%y,verts(2)%z)
      call glvertex3d(verts(2)%x,verts(2)%y,verts(2)%z)
      call glvertex3d(verts(3)%x,verts(3)%y,verts(3)%z)
      call glvertex3d(verts(3)%x,verts(3)%y,verts(3)%z)
      if (.not. verts(4) == NULL_POINT) then
         verts(4) = verts(4) + explode_shift(elem_owner(elem)) + &
                               explode_elem_shift(elem)
         call glvertex3d(verts(4)%x,verts(4)%y,verts(4)%z)
         call glvertex3d(verts(4)%x,verts(4)%y,verts(4)%z)
      endif
      call glvertex3d(verts(1)%x,verts(1)%y,verts(1)%z)
   endif

else ! not a leaf, do the children

   allc = ALL_CHILDREN
   children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call recur_draw_cutting_plane_lines(grid,nolight,children(i),plane, &
                                             value)
      endif
   end do

endif

end subroutine recur_draw_cutting_plane_lines

!          ---------------------------
subroutine draw_cutting_plane_contours(nolight)
!          ---------------------------

!----------------------------------------------------
! This routine draws contours on the enabled cutting planes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: plane, elem, i
real(my_real) :: contour_value(num_contour), minz, maxz, pval, psval, psmin, &
                 psmax, value
!----------------------------------------------------
! Begin executable code

! draw them black

call glcolor3d(0.0_gldouble,0.0_gldouble,0.0_gldouble)

! set the contour values

if (contour_values_given) then
   contour_value = actual_contours
else
   call preproc_and_scale(color_cutting_plane,1.0_my_real,pval,psval,minz, &
                          psmin,maxz,psmax)
   do i=1,num_contour
      contour_value(i) = minz+(maxz-minz)*(i-1)/real(num_contour-1,my_real)
   end do
endif

! do all three planes

do i=1,3

! pick the plane and see if it is active

   select case (i)
   case (1)
      if (yz_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = YZ_PLANE
      value = yz_cutting_plane_value
   case (2)
      if (xz_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = XZ_PLANE
      value = xz_cutting_plane_value
   case (3)
      if (xy_cutting_plane_value == CUTTING_PLANE_OFF) cycle
      plane = XY_PLANE
      value = xy_cutting_plane_value
   end select

! call the worker routine with each element of the initial grid

   elem = grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      call recur_draw_cut_plane_contours(nolight,elem,plane,contour_value, &
                                         value)
      elem = element(elem)%next
   end do

end do

end subroutine draw_cutting_plane_contours

!                    -----------------------------
recursive subroutine recur_draw_cut_plane_contours(nolight,elem,plane, &
                                                   contour_value,value)
!                    -----------------------------

!----------------------------------------------------
! This routine does the work of drawing the contours on the enabled cutting
! planes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical, intent(in) :: nolight
integer, intent(in) :: elem, plane
real(my_real), intent(in) :: contour_value(:)
real(my_real), intent(in) :: value
!----------------------------------------------------
! Local variables:

integer :: allc(MAX_CHILD), children(MAX_CHILD), i
type(point) :: verts(4)
real(my_real) :: x(3), y(3),z(3)
!----------------------------------------------------
! Begin executable code

! only draw the plane through leaf elements

if (element(elem)%isleaf) then

! get the intersection of the cutting plane and element

   verts = func_intersect_tet(grid,plane,value, &
                              vertex(element(elem)%vertex)%coord%x, &
                              vertex(element(elem)%vertex)%coord%y, &
                              vertex(element(elem)%vertex)%coord%z, &
                              elem)

! draw the edges if an intersection was found

   if (.not. verts(3) == NULL_POINT) then
      x = (/verts(1)%x,verts(2)%x,verts(3)%x/)
      y = (/verts(1)%y,verts(2)%y,verts(3)%y/)
      z = (/verts(1)%z,verts(2)%z,verts(3)%z/)
      call draw_triangle_contour(x,y,z,contour_value,color_cutting_plane,elem,0)
   endif
   if (.not. verts(4) == NULL_POINT) then
      x = (/verts(1)%x,verts(3)%x,verts(4)%x/)
      y = (/verts(1)%y,verts(3)%y,verts(4)%y/)
      z = (/verts(1)%z,verts(3)%z,verts(4)%z/)
      call draw_triangle_contour(x,y,z,contour_value,color_cutting_plane,elem,0)
   endif

else ! not a leaf, do the children

   allc = ALL_CHILDREN
   children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call recur_draw_cut_plane_contours(nolight,children(i),plane, &
                                            contour_value,value)
      endif
   end do

endif

end subroutine recur_draw_cut_plane_contours

!        ------------------
function func_intersect_tet(grid,func,value,x,y,z,elem)
!        ------------------

!----------------------------------------------------
! This routine returns the points where the preprocessed func equals
! value on the edges of the tetrahedon given by (x,y,z).  elem is an
! element that contains the tetrahedron.  There may 0, 3 or 4 intersections;
! excess entries in the result are NULL_POINT.  If there are 4 intersections
! the order is such that triangle(1,2,3) intersect triangle(1,3,4) is line(1,3)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: func
real(my_real), intent(in) :: value
real(my_real), intent(in) :: x(VERTICES_PER_ELEMENT),y(VERTICES_PER_ELEMENT), &
                             z(VERTICES_PER_ELEMENT)
integer, intent(in) :: elem
type(point) :: func_intersect_tet(4)
!----------------------------------------------------
! Local variables:

integer :: nfound, i, j,i1, i2
logical :: found_vert(VERTICES_PER_ELEMENT)
real(my_real) :: f_at_vert(VERTICES_PER_ELEMENT), val, ppval, ppssval, ppmin, &
                 ppssmin, ppmax, ppssmax, frac
!----------------------------------------------------
! Begin executable code

! initialize as none found

func_intersect_tet = NULL_POINT
nfound = 0
found_vert = .false.

! evaluate the function at the vertices

select case (func)

case (YZ_PLANE)
   f_at_vert = x
case (XZ_PLANE)
   f_at_vert = y
case (XY_PLANE)
   f_at_vert = z
case (DRAW_SOLUTION)
   call graphics_evaluate_soln(x,y,z,elem,f_at_vert)
case (DRAW_TRUE)
   do i=1,VERTICES_PER_ELEMENT
      f_at_vert(i) = graphics_trues(x(i),y(i),z(i))
   end do
case (DRAW_ERROR)
   call graphics_evaluate_soln(x,y,z,elem,f_at_vert)
   do i=1,VERTICES_PER_ELEMENT
      f_at_vert(i) = f_at_vert(i) - graphics_trues(x(i),y(i),z(i))
   end do
end select

! preprocess and scale the function values

if (func==DRAW_SOLUTION .or. &
    func==DRAW_TRUE .or. &
    func==DRAW_ERROR) then
   do i=1,VERTICES_PER_ELEMENT
      val = f_at_vert(i)
      call preproc_and_scale(func,val,ppval,ppssval,ppmin,ppssmin, &
                             ppmax,ppssmax)
      f_at_vert(i) = ppval
   end do
endif

! look for places where the function takes on value on an edge

! for each edge ...

do j=1,EDGES_PER_ELEMENT

! i1, i2 are the indices of the endpoints of the edge

   select case(j)
   case (1); i1 = 1; i2 = 2
   case (2); i1 = 1; i2 = 3
   case (3); i1 = 4; i2 = 1
   case (4); i1 = 4; i2 = 3
   case (5); i1 = 4; i2 = 2
   case (6); i1 = 2; i2 = 3
   end select

! see if the value is on this edge

   if (value >= f_at_vert(i1) .and.  value <= f_at_vert(i2)) then

! see if it is at one or both vertices; don't add a vertex more than once

      if (value == f_at_vert(i1) .or. value == f_at_vert(i2)) then
         if (value == f_at_vert(i1)) then
            if (.not. found_vert(i1)) then
               nfound = nfound+1
               func_intersect_tet(nfound)%x = x(i1)
               func_intersect_tet(nfound)%y = y(i1)
               func_intersect_tet(nfound)%z = z(i1)
               found_vert(i1) = .true.
            endif
         endif
         if (value == f_at_vert(i2)) then
            if (.not. found_vert(i2)) then
               nfound = nfound+1
               func_intersect_tet(nfound)%x = x(i2)
               func_intersect_tet(nfound)%y = y(i2)
               func_intersect_tet(nfound)%z = z(i2)
               found_vert(i2) = .true.
            endif
         endif

! not a vertex; find the interpolation point

      else
         frac = (value-f_at_vert(i1))/(f_at_vert(i2)-f_at_vert(i1))
         nfound = nfound+1
         func_intersect_tet(nfound)%x = (1.0_gldouble-frac)*x(i1) + frac*x(i2)
         func_intersect_tet(nfound)%y = (1.0_gldouble-frac)*y(i1) + frac*y(i2)
         func_intersect_tet(nfound)%z = (1.0_gldouble-frac)*z(i1) + frac*z(i2)
      endif
      if (nfound == 4) exit

   elseif (value >= f_at_vert(i2) .and.  value <= f_at_vert(i1)) then

      if (value == f_at_vert(i1) .or. value == f_at_vert(i2)) then
         if (value == f_at_vert(i1)) then
            if (.not. found_vert(i1)) then
               nfound = nfound+1
               func_intersect_tet(nfound)%x = x(i1)
               func_intersect_tet(nfound)%y = y(i1)
               func_intersect_tet(nfound)%z = z(i1)
               found_vert(i1) = .true.
            endif
         endif
         if (value == f_at_vert(i2)) then
            if (.not. found_vert(i2)) then
               nfound = nfound+1
               func_intersect_tet(nfound)%x = x(i2)
               func_intersect_tet(nfound)%y = y(i2)
               func_intersect_tet(nfound)%z = z(i2)
               found_vert(i2) = .true.
            endif
         endif
      else
         frac = (value-f_at_vert(i2))/(f_at_vert(i1)-f_at_vert(i2))
         nfound = nfound+1
         func_intersect_tet(nfound)%x = (1.0_gldouble-frac)*x(i2) + frac*x(i1)
         func_intersect_tet(nfound)%y = (1.0_gldouble-frac)*y(i2) + frac*y(i1)
         func_intersect_tet(nfound)%z = (1.0_gldouble-frac)*z(i2) + frac*z(i1)
      endif
      if (nfound == 4) exit

   endif

end do

end function func_intersect_tet

!                    ---------------------
recursive subroutine draw_triangle_contour(xmr,ymr,zmr,contour_value, &
                                           use_color,elem,recur_lev)
!                    ---------------------

!----------------------------------------------------
! This routine draws the contours in a triangle with vertices (xmr,ymr,zmr)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: xmr(:), ymr(:), zmr(:), contour_value(:)
integer, intent(in) :: use_color, elem, recur_lev
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: xmid(3),ymid(3),zmid(3),fmr(3)
real(gldouble) :: x(3),y(3),z(3),f(3)
real(gldouble) :: xt, yt, zt, ft, frac, xcross1, xcross2, ycross1, ycross2, &
                  zcross1, zcross2, fcross
real(my_real) :: val,ppval,ppssval,ppmin,ppssmin,ppmax,ppssmax
real(my_real) :: delta=.05_my_real
!----------------------------------------------------
! Begin executable code

if (recur_lev < subelement_resolution) then

! draw four children of this triangle for better resolution

! define midpoints of the sides

   xmid(1) = (xmr(2)+xmr(3))/2
   ymid(1) = (ymr(2)+ymr(3))/2
   zmid(1) = (zmr(2)+zmr(3))/2
   xmid(2) = (xmr(3)+xmr(1))/2
   ymid(2) = (ymr(3)+ymr(1))/2
   zmid(2) = (zmr(3)+zmr(1))/2
   xmid(3) = (xmr(1)+xmr(2))/2
   ymid(3) = (ymr(1)+ymr(2))/2
   zmid(3) = (zmr(1)+zmr(2))/2

! draw the children

   call draw_triangle_contour((/xmr(1),xmid(3),xmid(2)/), &
                              (/ymr(1),ymid(3),ymid(2)/), &
                              (/zmr(1),zmid(3),zmid(2)/), &
                              contour_value,use_color,elem,recur_lev+1)
   call draw_triangle_contour((/xmr(2),xmid(1),xmid(3)/), &
                              (/ymr(2),ymid(1),ymid(3)/), &
                              (/zmr(2),zmid(1),zmid(3)/), &
                              contour_value,use_color,elem,recur_lev+1)
   call draw_triangle_contour((/xmr(3),xmid(2),xmid(1)/), &
                              (/ymr(3),ymid(2),ymid(1)/), &
                              (/zmr(3),zmid(2),zmid(1)/), &
                              contour_value,use_color,elem,recur_lev+1)
   call draw_triangle_contour((/xmid(1),xmid(2),xmid(3)/), &
                              (/ymid(1),ymid(2),ymid(3)/), &
                              (/zmid(1),zmid(2),zmid(3)/), &
                              contour_value,use_color,elem,recur_lev+1)

else ! draw the triangle


! (xmr,ymr,zmr) are used for evaluating solution; (x,y,z) are used as vertices

   x = xmr + explode_shift(elem_owner(elem))%x + explode_elem_shift(elem)%x
   y = ymr + explode_shift(elem_owner(elem))%y + explode_elem_shift(elem)%y
   z = zmr + explode_shift(elem_owner(elem))%z + explode_elem_shift(elem)%z

! get the function value at the vertices

   select case (use_color)
   case (COLOR_SOLUTION)
      call graphics_evaluate_soln(xmr,ymr,zmr,elem,fmr)
      f = fmr
   case (COLOR_TRUE)
      do i=1,VERTICES_PER_ELEMENT
         f(i) = graphics_trues(xmr(i),ymr(i),zmr(i))
      end do
   case (COLOR_ERROR)
      call graphics_evaluate_soln(xmr,ymr,zmr,elem,fmr)
      do i=1,VERTICES_PER_ELEMENT
         f(i) = fmr(i) - graphics_trues(xmr(i),ymr(i),zmr(i))
      end do
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("draw_triangle_contour: not set up for contour with function choice",intlist=(/use_color/))
      stop
   end select
   do i=1,VERTICES_PER_ELEMENT
      val = f(i)
      call preproc_and_scale(use_color,val,ppval,ppssval,ppmin,ppssmin,ppmax, &
                             ppssmax)
      f(i) = ppval
   end do

! order the vertices by f value

   xt = x(1); yt = y(1); zt = z(1); ft = f(1)
   if (f(2) < f(1)) then
      xt = x(1); yt = y(1); zt = z(1); ft = f(1)
      if (f(3) < f(1)) then
         if (f(2) < f(3)) then
            x(1) = x(2); y(1) = y(2); z(1) = z(2); f(1) = f(2)
            x(2) = x(3); y(2) = y(3); z(2) = z(3); f(2) = f(3)
         else
            x(1) = x(3); y(1) = y(3); z(1) = z(3); f(1) = f(3)
         endif
         x(3) = xt; y(3) = yt; z(3) = zt; f(3) = ft
      else
         x(1) = x(2); y(1) = y(2); z(1) = z(2); f(1) = f(2)
         x(2) = xt; y(2) = yt; z(2) = zt; f(2) = ft
      endif
   elseif (f(3) < f(1)) then
      x(1) = x(3); y(1) = y(3); z(1) = z(3); f(1) = f(3)
      x(3) = x(2); y(3) = y(2); z(3) = z(2); f(3) = f(2)
      x(2) = xt; y(2) = yt; z(2) = zt; f(2) = ft
   elseif (f(3) < f(2)) then
      xt = x(2); yt = y(2); zt = z(2); ft = f(2)
      x(2) = x(3); y(2) = y(3); z(2) = z(3); f(2) = f(3)
      x(3) = xt; y(3) = yt; z(3) = zt; f(3) = ft
   endif

! if f(1)==f(3), the function is flat in the element and has no contours

   if (f(1) /= f(3)) then

! for each contour value

      do i = 1,num_contour

! see if it passes through this element

         if (contour_value(i) < f(1)) cycle
         if (contour_value(i) > f(3)) exit

! see where it crosses the 1-3 edge

         frac = (contour_value(i)-f(1))/(f(3)-f(1))
         xcross1 = (1.0_gldouble - frac)*x(1) + frac*x(3)
         ycross1 = (1.0_gldouble - frac)*y(1) + frac*y(3)
         zcross1 = (1.0_gldouble - frac)*z(1) + frac*z(3)

! see where it crosses one of the other edges

         if (contour_value(i) == f(2)) then
            xcross2 = x(2)
            ycross2 = y(2)
            zcross2 = z(2)
         elseif (contour_value(i) < f(2)) then
            frac = (contour_value(i)-f(1))/(f(2)-f(1))
            xcross2 = (1.0_gldouble - frac)*x(1) + frac*x(2)
            ycross2 = (1.0_gldouble - frac)*y(1) + frac*y(2)
            zcross2 = (1.0_gldouble - frac)*z(1) + frac*z(2)
         else
            frac = (contour_value(i)-f(2))/(f(3)-f(2))
            xcross2 = (1.0_gldouble - frac)*x(2) + frac*x(3)
            ycross2 = (1.0_gldouble - frac)*y(2) + frac*y(3)
            zcross2 = (1.0_gldouble - frac)*z(2) + frac*z(3)
         endif

! add the line segment to the display list

         call glvertex3d(xcross1,ycross1,zcross1)
         call glvertex3d(xcross2,ycross2,zcross2)

      end do ! next contour
   endif ! not flat element

endif ! check for recursion

end subroutine draw_triangle_contour

!                    --------
recursive subroutine draw_sfc(elem)
!                    --------

!----------------------------------------------------
! This routine draws the space filling curve used for partitioning
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(gldouble) :: x, y, z
integer :: i, j, nelem, allocstat
integer :: allc(MAX_CHILD), children(MAX_CHILD)
integer, allocatable :: initial_elements(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

! the list of elements is in reverse order; make a new list to reverse it

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   nelem = 0
   do while (i /= END_OF_LIST)
      nelem = nelem+1
      i = element(i)%next
   end do
   allocate(initial_elements(nelem),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed in draw_sfc",intlist=(/allocstat/))
      stop
   endif
   i = grid%head_level_elem(1)
   j = 0
   do while (i /= END_OF_LIST)
      j = j+1
      initial_elements(j) = i
      i = element(i)%next
   end do
   do i=nelem,1,-1
      call draw_sfc(initial_elements(i))
   end do
   deallocate(initial_elements)
   return
endif

! add the center of this element to the SFC only if it is a leaf

if (element(elem)%isleaf) then

! set the midpoint

   x = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%x + explode_elem_shift(elem)%x
   y = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%y + explode_elem_shift(elem)%y
   z = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%z)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%z + explode_elem_shift(elem)%z

! add the midpoint to the display list

   call glvertex3d(x,y,z)

endif

! draw the sfc in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(element(elem)%order(i)) /= NO_CHILD) then
      call draw_sfc(children(element(elem)%order(i)))
   endif
end do

return
end subroutine draw_sfc

!                    ----------
recursive subroutine draw_sfcio(elem,h)
!                    ----------

!----------------------------------------------------
! This routine labels the in and out vertices for the space filling curve
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
real(gldouble), intent(in) :: h
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(gldouble) :: x, y, z, hh
real(gldouble), parameter :: shiftx = .1_gldouble, shifty = .01_gldouble
integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      hh = 8*elem_diam(i)
      call draw_sfcio(i,hh)
      i = element(i)%next
   end do
   return
endif

! add the labels to this element only if it is a leaf

if (element(elem)%isleaf) then

! set the midpoint

   x = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%x + explode_elem_shift(elem)%x
   y = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%y + explode_elem_shift(elem)%y
   z = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%z)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%z + explode_elem_shift(elem)%z

! label the in vertex

   call number(1,(x+vertex(element(elem)%in)%coord%x)/2 - h*shiftx, &
                 (y+vertex(element(elem)%in)%coord%y)/2 - h*shifty, &
                 z,h)

! label the out vertex

   call number(0,(x+vertex(element(elem)%out)%coord%x)/2 - h*shiftx, &
                 (y+vertex(element(elem)%out)%coord%y)/2 - h*shifty, &
                 z,h)

endif

! label the elements in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(element(elem)%order(i)) /= NO_CHILD) then
     call draw_sfcio(children(element(elem)%order(i)),h/sqrt(2.0_gldouble))
   endif
end do

return
end subroutine draw_sfcio

!          --------------
subroutine get_edge_verts(grid,elem,edge_verts)
!          --------------

!----------------------------------------------------
! This routine returns the index in grid%element(elem)%vertex of the
! endpoints of the edges of elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
integer, intent(out) :: edge_verts(:,:)
!----------------------------------------------------
! Local variables:

integer :: i,j,k,edge,vert
!----------------------------------------------------
! Begin executable code

do i=1,EDGES_PER_ELEMENT
   edge = grid%element(elem)%edge(i)
   do j=1,2
      vert = grid%edge(edge)%vertex(j)
      edge_verts(i,j) = -1
      do k=1,VERTICES_PER_ELEMENT
         if (grid%element(elem)%vertex(k) == vert) then
            edge_verts(i,j) = k
            exit
         endif
      end do
      if (edge_verts(i,j) == -1) then
         call fatal("get_edge_verts: didn't find edge vertex in the element")
         stop
      endif
   end do
end do

end subroutine get_edge_verts

!          --------------
subroutine get_face_verts(grid,elem,face_verts)
!          --------------

!----------------------------------------------------
! This routine returns the index in grid%element(elem)%vertex of the
! vertices of the faces of elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
integer, intent(out) :: face_verts(:,:)
!----------------------------------------------------
! Local variables:

integer :: i,j,k,face,vert
!----------------------------------------------------
! Begin executable code

do i=1,FACES_PER_ELEMENT
   face = grid%element(elem)%face(i)
   do j=1,VERTICES_PER_FACE
      vert = grid%face(face)%vertex(j)
      face_verts(i,j) = -1
      do k=1,VERTICES_PER_ELEMENT
         if (grid%element(elem)%vertex(k) == vert) then
            face_verts(i,j) = k
            exit
         endif
      end do
      if (face_verts(i,j) == -1) then
         call fatal("get_face_verts: didn't find face vertex in the element")
         stop
      endif
   end do
end do

end subroutine get_face_verts

!          ----------------------
subroutine graphics_evaluate_soln(x,y,z,elem,soln)
!          ----------------------

!----------------------------------------------------
! This routine calls evaluate solution with the current eigen and
! compnt taking care when compnt is the L1 or L2 sum
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:), y(:), z(:)
integer, intent(in) :: elem
real(my_real), intent(out) :: soln(:)
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: solnloc(ncompnt,1,size(soln)), soln1(1,1,size(soln))
!----------------------------------------------------
! Begin executable code

select case (compnt)

case(-1)
   call evaluate_soln_local(grid,x,y,elem,(/(i,i=1,ncompnt)/),(/eigen/), &
                            solnloc,z=z)
   soln = 0.0_my_real
   do i=1,ncompnt
      soln = soln + abs(solnloc(i,1,:))
   end do

case(-2)
   call evaluate_soln_local(grid,x,y,elem,(/(i,i=1,ncompnt)/),(/eigen/), &
                            solnloc,z=z)
   soln = 0.0_my_real
   do i=1,ncompnt
      soln = soln + solnloc(i,1,:)**2
   end do

case default
   call evaluate_soln_local(grid,x,y,elem,(/compnt/),(/eigen/),soln1,z=z)
   soln = soln1(1,1,:)

end select

end subroutine graphics_evaluate_soln

!        --------------
function graphics_trues(x,y,z)
!        --------------

!----------------------------------------------------
! This routine calls trues for the current eigen and compnt taking
! care when compnt is the L1 or L2 sum
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
real(my_real) :: graphics_trues
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: trueloc(ncompnt)
!----------------------------------------------------
! Begin executable code

select case (compnt)

case(-1)
   do i=1,ncompnt
      trueloc(i) = my_trues(x,y,z,i,eigen)
   end do
   graphics_trues = sum(abs(trueloc))

case(-2)
   do i=1,ncompnt
      trueloc(i) = my_trues(x,y,z,i,eigen)
   end do
   graphics_trues = sum(trueloc**2)

case default
   graphics_trues = my_trues(x,y,z,compnt,eigen)

end select

end function graphics_trues

!        ---------
function elem_diam(elem)
!        ---------

!----------------------------------------------------
! This routine computes the diameter of an element, defined to be
! the minimum distance from the center to one of the vertices.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
real(gldouble) :: elem_diam
!----------------------------------------------------
! Local variables:

real(gldouble) :: x, y, mindist
integer :: i
!----------------------------------------------------
! Begin executable code

! set the midpoint

x = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x)/VERTICES_PER_ELEMENT
y = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y)/VERTICES_PER_ELEMENT

! compute the distance to each vertex, keeping the minimum

mindist = huge(0.0_gldouble)

do i=1,VERTICES_PER_ELEMENT
   mindist = min(mindist,sqrt((vertex(element(elem)%vertex(i))%coord%x-x)**2 + &
                              (vertex(element(elem)%vertex(i))%coord%y-y)**2))
end do

elem_diam = mindist

end function elem_diam

!          ---------
subroutine get_rainbow(val,minval,maxval,c)
!          ---------

!----------------------------------------------------
! This routine sets the color for rainbow plots
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments:

real(gldouble), intent(in) :: val,maxval,minval
real(gldouble), intent(out) :: c(4)

!----------------------------------------------------
! Local variables:

real(gldouble) :: f,h,s,v,sbase
integer :: hbase

!----------------------------------------------------
! Begin executable code

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = gldhalf
endif
if (f < 0.0) f = 0
if (f > 1.0) f = 1

if (color_scheme == SCHEME_STEP_SEQ .and. f==1.0_gldouble) f=0.9999_gldouble

select case(color_scheme)

case (SCHEME_RAINBOW)
   if (f < .25) then
      c(1) = 0.0_gldouble
      c(2) = 4.0_gldouble * f
      c(3) = 1.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < .5) then
      c(1) = 0.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 2.0_gldouble - 4.0_gldouble*f
      c(4) = 1.0_gldouble
   elseif (f < .75) then
      c(1) = 4.0_gldouble * f - 2.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   else
      c(1) = 1.0_gldouble
      c(2) = 4.0_gldouble - 4.0_gldouble*f
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   endif

case (SCHEME_DOUBLE_RAINBOW)
   f = 2*f
   if (f < .25) then
      c(1) = 0.0_gldouble
      c(2) = 4.0_gldouble * f
      c(3) = 1.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < .5) then
      c(1) = 0.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 2.0_gldouble - 4.0_gldouble*f
      c(4) = 1.0_gldouble
   elseif (f < .75) then
      c(1) = 4.0_gldouble * f - 2.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < 1.0) then
      c(1) = 1.0_gldouble
      c(2) = 4.0_gldouble - 4.0_gldouble*f
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < 1.25) then
      c(1) = 0.5_gldouble
      c(2) = 0.5_gldouble + 2.0_gldouble * (f-1)
      c(3) = 1.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < 1.5) then
      c(1) = 0.5_gldouble
      c(2) = 1.0_gldouble
      c(3) = 0.5_gldouble + (2.0_gldouble - 4.0_gldouble*(f-1))/2
      c(4) = 1.0_gldouble
   elseif (f < 1.75) then
      c(1) = .5_gldouble+(4.0_gldouble * (f-1) - 2.0_gldouble)/2
      c(2) = 1.0_gldouble
      c(3) = 0.5_gldouble
      c(4) = 1.0_gldouble
   else
      c(1) = 1.0_gldouble
      c(2) = .2_gldouble+(4.0_gldouble - 4.0_gldouble*(f-1))/2
      c(3) = 0.5_gldouble
      c(4) = 1.0_gldouble
   endif

case (SCHEME_GRAY)
   c = f
   c(4) = 1.0_gldouble

case (SCHEME_STEP_SEQ)
   c(4) = 1.0_gldouble
   hbase = floor(f*step_scheme_hues)
   sbase = f*step_scheme_hues - hbase
   h = (360.0_gldouble*hbase)/step_scheme_hues
   s = 0.2_gldouble + 0.8_gldouble*(1-sbase)
   v = 0.8_gldouble + 0.2_gldouble*sbase
   call hsv2rgb(h,s,v,c(1),c(2),c(3))

end select

return
end subroutine get_rainbow

!          -------
subroutine hsv2rgb(h,s,v,r,g,b)
!          -------

!----------------------------------------------------
! This routine converts (h,s,v) color representation to (r,g,b).
! The formulas come from the Wikipedia page on HSL and HSV
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(gldouble), intent(in) :: h,s,v
real(gldouble), intent(out) :: r,g,b
!----------------------------------------------------
! Local variables:

real(gldouble) :: f, p, q, t
!----------------------------------------------------
! Begin executable code

f = h/60 - floor(h/60)
p = v*(1-s)
q = v*(1-f*s)
t = v*(1-(1-f)*s)
select case (mod(floor(h/60),6))
case (0)
   r=v; g=t; b=p
case (1)
   r=q; g=v; b=p
case (2)
   r=p; g=v; b=t
case (3)
   r=p; g=q; b=v
case (4)
   r=t; g=p; b=v
case (5)
   r=v; g=p; b=q
case default
   print *,"bad case in hsv2rgb"
   r=0; g=0; b=0
end select

end subroutine hsv2rgb

!        -------------
function normcrossprod(x,y,z)
!        -------------

!----------------------------------------------------
! The function returns the normalized cross product of the vectors
! (p1,p2) and (p1,p3) with the coordinates of the p's in x,y,z.  If the
! points are counterclockwise, the cross product is negated to point in
! the direction of a clockwise orientation of the points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(gldouble), dimension(:) :: x,y,z
real(gldouble), dimension(3) :: normcrossprod

!----------------------------------------------------
! Local variables

real(gldouble) :: t1(3),t2(3),norm

!----------------------------------------------------
! Begin executable code

t1(1) = x(2) - x(1)
t1(2) = y(2) - y(1)
t1(3) = z(2) - z(1)
t2(1) = x(3) - x(1)
t2(2) = y(3) - y(1)
t2(3) = z(3) - z(1)

normcrossprod(1) = t1(2)*t2(3) - t1(3)*t2(2)
normcrossprod(2) = t1(3)*t2(1) - t1(1)*t2(3)
normcrossprod(3) = t1(1)*t2(2) - t1(2)*t2(1)

if (normcrossprod(3) > 0.0_gldouble) normcrossprod = -normcrossprod

norm = sqrt(dot_product(normcrossprod,normcrossprod))
if (norm /= gldzero) normcrossprod = normcrossprod/norm

return
end function normcrossprod

!                -----------
logical function inside_crop(vert)
!                -----------

!----------------------------------------------------
! This routine determines if at least one of the vertices is inside
! the croping range
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: vert(:)

!----------------------------------------------------
! Local variables:
integer :: i
real(my_real) :: x,y,z

!----------------------------------------------------
! Begin executable code

inside_crop = .false.
do i=1,VERTICES_PER_ELEMENT
   x = vertex(vert(i))%coord%x
   y = vertex(vert(i))%coord%y
   z = vertex(vert(i))%coord%z
   if (x>xcrop1 .and. x<xcrop2 .and. y>ycrop1 .and. y<ycrop2 .and. &
       z>zcrop1 .and. z<zcrop2) inside_crop = .true.
end do

return
end function inside_crop

!          ---------
subroutine draw_axes
!          ---------

!----------------------------------------------------
! This routine draws x, y and z axes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

real(gldouble), parameter :: delta = .05_gldouble, half=.5_gldouble
real(gldouble) :: xstart,xend,ystart,yend,zstart,zend
real(gldouble) :: xmind,xmaxd,ymind,ymaxd,zmind,zmaxd
real(gldouble) :: sdelta
real(my_real) :: hold_xmin, hold_ymin
!----------------------------------------------------
! Begin executable code

if (convert_cylindrical_to_cartesian) then
   hold_xmin = xmin
   hold_ymin = ymin
   xmin = 0
   ymin = 0
endif

xmind = xmin
xmaxd = xmax
ymind = ymin
ymaxd = ymax
zmind = zmin
zmaxd = zmax

sdelta = delta*min(min(xmax-xmin,ymax-ymin),zmax-zmin)

! set the endpoints for the axes

xstart = xmin - delta*(xmax-xmin)
xend   = xmax + delta*(xmax-xmin)
ystart = ymin - delta*(ymax-ymin)
yend   = ymax + delta*(ymax-ymin)
zstart = zmin - delta*(zmax-zmin)
zend   = zmax + delta*(zmax-zmin)

if (convert_cylindrical_to_cartesian) then
   xstart = xmin
   ystart = ymin
   zstart = zmin
endif

! draw the axes

call glbegin(gl_lines)
call glvertex3d(xstart,ystart,zstart) ! x axis
call glvertex3d(xend,ystart,zstart)
call glvertex3d(xstart,ystart,zstart) ! y axis
call glvertex3d(xstart,yend,zstart)
call glvertex3d(xstart,ystart,zstart) ! z axis
call glvertex3d(xstart,ystart,zend)

! add ticks

call glvertex3d(xmind,ystart-half*sdelta,zstart) ! x ticks
call glvertex3d(xmind,ystart+half*sdelta,zstart)
call glvertex3d(half*(xmind+xmaxd),ystart-half*sdelta,zstart)
call glvertex3d(half*(xmind+xmaxd),ystart+half*sdelta,zstart)
call glvertex3d(xmaxd,ystart-half*sdelta,zstart)
call glvertex3d(xmaxd,ystart+half*sdelta,zstart)

call glvertex3d(xstart-half*sdelta,ymind,zstart) ! y ticks
call glvertex3d(xstart+half*sdelta,ymind,zstart)
call glvertex3d(xstart-half*sdelta,half*(ymind+ymaxd),zstart)
call glvertex3d(xstart+half*sdelta,half*(ymind+ymaxd),zstart)
call glvertex3d(xstart-half*sdelta,ymaxd,zstart)
call glvertex3d(xstart+half*sdelta,ymaxd,zstart)

call glvertex3d(xstart,ystart-half*sdelta,zmind) ! z ticks
call glvertex3d(xstart,ystart+half*sdelta,zmind)
call glvertex3d(xstart,ystart-half*sdelta,half*(zmind+zmaxd))
call glvertex3d(xstart,ystart+half*sdelta,half*(zmind+zmaxd))
call glvertex3d(xstart,ystart-half*sdelta,zmaxd)
call glvertex3d(xstart,ystart+half*sdelta,zmaxd)

call glend

! label the axes

call real_number(xmin,xmind,ystart-sdelta,zstart,10.*sdelta)
call real_number((xmin+xmax)/2.,half*(xmind+xmaxd),ystart-sdelta, &
                 zstart,10.*sdelta)
call real_number(xmax,xmaxd,ystart-sdelta,zstart,10.*sdelta)

call real_number(ymin,xstart-sdelta,ymind,zstart,10.*sdelta)
call real_number((ymin+ymax)/2.,xstart-sdelta,half*(ymind+ymaxd), &
                 zstart,10.*sdelta)
call real_number(ymax,xstart-sdelta,ymaxd,zstart,10.*sdelta)

call real_number(zmin,xstart,ystart-sdelta,zmind,10.*sdelta)
call real_number((zmin+zmax)/2.,xstart,ystart-sdelta, &
                 half*(zmind+zmaxd),10.*sdelta)
call real_number(zmax,xstart,ystart-sdelta,zmaxd,10.*sdelta)

if (convert_cylindrical_to_cartesian) then
   xmin = hold_xmin
   ymin = hold_ymin
endif

return
end subroutine draw_axes

!          ------
subroutine number(n,x,y,z,h)
!          ------

!----------------------------------------------------
! This routine draws the number n at the point (x,y,z) using h as a
! scaling factor for the height of the number
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: n
real(gldouble), intent(in) :: x,y,z,h

!----------------------------------------------------
! Local variables

character(len=9) :: chn
integer :: i

!----------------------------------------------------
! Begin executable code

write(chn,'(i9)') n
if (chn /= " ") then
   do while (chn(1:1) == " ")
      chn = chn(2:9)
   end do
endif
call glmatrixmode(gl_modelview)
call glpushmatrix
call gltranslated(x,y,z)
! shift from bottom left corner at center of element to something more central;
! the divisor was determined experimentally and may not always work
call gltranslated(-h/35,0.0_gldouble,-h/35)
! the next line rotates the number 90 degrees around the x axis
!call glrotated(90.0_gldouble,1.0_gldouble, 0.0_gldouble, 0.0_gldouble)
call glscaled(h/2500.0_gldouble,h/2500.0_gldouble,1._gldouble)
do i=1,9
   call glutstrokecharacter(glut_stroke_roman,ichar(chn(i:i)))
end do
call glpopmatrix

return
end subroutine number

!          -----------
subroutine real_number(n,x,y,z,h)
!          -----------

!----------------------------------------------------
! This routine draws the real number n at the point (x,y,z) using h as a
! scaling factor for the height of the number
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: n
real(gldouble), intent(in) :: x,y,z,h

!----------------------------------------------------
! Local variables

character(len=8) chn
integer :: i

!----------------------------------------------------
! Begin executable code

write(chn,'(g8.2)') n
call glmatrixmode(gl_modelview)
call glpushmatrix
call gltranslated(x-h/8.,y,z)
call glscaled(h/2500.0_gldouble,h/2500.0_gldouble,1._gldouble)
do i=1,8
   call glutstrokecharacter(glut_stroke_roman,ichar(chn(i:i)))
end do
call glpopmatrix

return
end subroutine real_number

!          ------------
subroutine grid_display
!          ------------
if (.not. window_initialized) return
call glutsetwindow(grid_win)
call reset_view
if (grid_newview) call recalcview
if (explode_factor /= old_explode_factor) then
   call draw_grid
   old_explode_factor = explode_factor
endif
if (explode_elem_factor /= old_explode_elem_factor) then
   call draw_grid
   old_explode_elem_factor = explode_elem_factor
endif
call glclear(ior(gl_color_buffer_bit,gl_depth_buffer_bit))
call glcalllist(grid_list)
call glutswapbuffers
end subroutine grid_display

!          ----------
subroutine recalcview
!          ----------
grid_newview = .false.
end subroutine recalcview

!          --------------
subroutine set_draw_lines(selection)
!          --------------
integer(glcint), intent(in out) :: selection
color_lines = selection
call draw_grid
return
end subroutine set_draw_lines

!          -------------
subroutine set_draw_face(selection)
!          -------------
integer(glcint), intent(in out) :: selection
color_faces = selection
call draw_grid
return
end subroutine set_draw_face

!          ---------------------
subroutine set_edge_transparency(selection)
!          ---------------------
integer(glcint), intent(in out) :: selection
select case(selection)
case (-1)
   edge_transparency = min(100.0_gldouble,edge_transparency+.01_gldouble)
case (-2)
   edge_transparency = min(100.0_gldouble,edge_transparency+.05_gldouble)
case (-3)
   edge_transparency = max(0.0_gldouble,edge_transparency-.01_gldouble)
case (-4)
   edge_transparency = max(0.0_gldouble,edge_transparency-.05_gldouble)
case default
   edge_transparency = 1.0_gldouble - selection/100.0_gldouble
end select
call draw_grid
return
end subroutine set_edge_transparency

!          ---------------------
subroutine set_face_transparency(selection)
!          ---------------------
integer(glcint), intent(in out) :: selection
select case(selection)
case (-1)
   face_transparency = min(100.0_gldouble,face_transparency+.01_gldouble)
case (-2)
   face_transparency = min(100.0_gldouble,face_transparency+.05_gldouble)
case (-3)
   face_transparency = max(0.0_gldouble,face_transparency-.01_gldouble)
case (-4)
   face_transparency = max(0.0_gldouble,face_transparency-.05_gldouble)
case default
   face_transparency = 1.0_gldouble - selection/100.0_gldouble
end select
call draw_grid
return
end subroutine set_face_transparency

!          ---------------------------
subroutine set_isosurface_transparency(selection)
!          ---------------------------
integer(glcint), intent(in out) :: selection
select case(selection)
case (-1)
   isosurface_transparency = &
      min(100.0_gldouble,isosurface_transparency+.01_gldouble)
case (-2)
   isosurface_transparency = &
      min(100.0_gldouble,isosurface_transparency+.05_gldouble)
case (-3)
   isosurface_transparency = &
      max(0.0_gldouble,isosurface_transparency-.01_gldouble)
case (-4)
   isosurface_transparency = &
      max(0.0_gldouble,isosurface_transparency-.05_gldouble)
case default
   isosurface_transparency = 1.0_gldouble - selection/100.0_gldouble
end select
call draw_grid
return
end subroutine set_isosurface_transparency

!          -----------------------
subroutine set_cutting_plane_color(selection)
!          -----------------------
integer(glcint), intent(in out) :: selection
color_cutting_plane = selection
call draw_grid
return
end subroutine set_cutting_plane_color

!          -----------------------
subroutine set_cutting_plane_edge(selection)
!          -----------------------
integer(glcint), intent(in out) :: selection
color_cutting_plane_lines = selection
call draw_grid
return
end subroutine set_cutting_plane_edge

!          -----------------
subroutine set_cutting_plane(selection)
!          -----------------
integer(glcint), intent(in out) :: selection
select case (selection)
case (YZ_PLANE)
   if (yz_cutting_plane_value == CUTTING_PLANE_OFF) then
      yz_cutting_plane_value = (xmin+xmax)/2
   else
      yz_cutting_plane_value = CUTTING_PLANE_OFF
   endif
case (XZ_PLANE)
   if (xz_cutting_plane_value == CUTTING_PLANE_OFF) then
      xz_cutting_plane_value = (ymin+ymax)/2
   else
      xz_cutting_plane_value = CUTTING_PLANE_OFF
   endif
case (XY_PLANE)
   if (xy_cutting_plane_value == CUTTING_PLANE_OFF) then
      xy_cutting_plane_value = (zmin+zmax)/2
   else
      xy_cutting_plane_value = CUTTING_PLANE_OFF
   endif
case (1_glcint)
   cutting_plane_speed = 2*cutting_plane_speed
case (2_glcint)
   cutting_plane_speed = cutting_plane_speed/2
end select
call draw_grid
return
end subroutine set_cutting_plane

!          ----------------
subroutine set_color_scheme(selection)
!          ----------------
integer(glcint), intent(in out) :: selection
color_scheme = selection
call set_owner_color
call draw_grid
return
end subroutine set_color_scheme

!          ---------------
subroutine set_step_scheme(selection)
!          ---------------
integer(glcint), intent(in out) :: selection
integer(glcint) :: hold
character(len=32) :: str
select case(selection)
case(1)
   step_scheme_steps=5
case(2)
   step_scheme_steps=10
case(3)
   step_scheme_steps=16
case(4)
   step_scheme_steps=max(1,step_scheme_steps-1)
case(5)
   step_scheme_steps=step_scheme_steps+1
case(6)
   step_scheme_hues=5
case(7)
   step_scheme_hues=10
case(8)
   step_scheme_hues=16
case(9)
   step_scheme_hues=max(1,step_scheme_hues-1)
case(10)
   step_scheme_hues=step_scheme_hues+1
end select
hold = glutgetmenu()
call glutsetmenu(scheme_menu)
write(str,"(A18,2I4)") "stepped sequential",step_scheme_steps,step_scheme_hues
call glutchangetomenuentry(5,trim(str),SCHEME_STEP_SEQ)
call glutsetmenu(hold)
call set_owner_color
call draw_grid
end subroutine set_step_scheme

!          -----------------
subroutine set_contour_param(selection)
!          -----------------
integer(glcint), intent(in out) :: selection
integer :: allocstat

select case(selection)
   case(1)
      print *,"enter number of uniform contour lines"
      read *,num_contour
      contour_values_given = .false.
   case(2)
      print *,"enter number of contours:"
      read *, num_contour
      if (allocated(actual_contours)) deallocate(actual_contours,stat=allocstat)
      allocate(actual_contours(num_contour),stat=allocstat)
      if (allocstat /= 0) then
         call fatal("allocation failed for actual_contours", &
                    intlist=(/allocstat/))
         stop
      endif
      print *,"enter ",num_contour," contour values:"
      read *,actual_contours
      contour_values_given = .true.
   case(3)
      num_contour = num_contour + 1
      contour_values_given = .false.
   case(4)
      num_contour = num_contour - 1
      if (num_contour < 2) num_contour = 2
      contour_values_given = .false.
   case(5)
      num_contour = num_contour + 10
      contour_values_given = .false.
   case(6)
      num_contour = num_contour - 10
      if (num_contour < 2) num_contour = 2
      contour_values_given = .false.
   case(7)
      num_contour = num_contour*2 - 1
      contour_values_given = .false.
   case(8)
      num_contour = (num_contour+1)/2
      if (num_contour < 2) num_contour = 2
      contour_values_given = .false.
   case(9)
      cutting_plane_contour = .not. cutting_plane_contour
end select

grid_newview = .true.
call draw_grid
end subroutine set_contour_param

!          --------------
subroutine set_isosurface(selection)
!          --------------
integer(glcint), intent(in out) :: selection
draw_isosurface = selection
grid_newview = .true.
call draw_grid
end subroutine set_isosurface

!          --------------------
subroutine set_isosurface_param(selection)
!          --------------------
integer(glcint), intent(in out) :: selection
integer :: allocstat

select case(selection)
   case(3)
      num_isosurface = num_isosurface + 1
   case(4)
      num_isosurface = num_isosurface - 1
      if (num_isosurface < 2) num_isosurface = 2
   case(5)
      num_isosurface = num_isosurface + 10
   case(6)
      num_isosurface = num_isosurface - 10
      if (num_isosurface < 2) num_isosurface = 2
   case(7)
      num_isosurface = num_isosurface*2 - 1
   case(8)
      num_isosurface = (num_isosurface+1)/2
      if (num_isosurface < 2) num_isosurface = 2
end select

grid_newview = .true.
call draw_grid
end subroutine set_isosurface_param

!          ----------------
subroutine set_preproc_func(selection)
!          ----------------
integer(glcint), intent(in out) :: selection
preproc_func = selection
grid_newview = .true.
call draw_grid
return
end subroutine set_preproc_func

!          -------
subroutine set_sfc(selection)
!          -------
integer(glcint), intent(in out) :: selection

select case(selection)
   case(0)
      draw_sfc_flag = .not. draw_sfc_flag
   case(1)
      draw_sfcio_flag = .not. draw_sfcio_flag
end select

grid_newview = .true.
call draw_grid
end subroutine set_sfc

!          ----------
subroutine set_offset(selection)
!          ----------
integer(glcint), intent(in out) :: selection

offset = offset + selection
if (offset < 0.0_glfloat) offset = 0.0_glfloat
call glpolygonoffset(offset,0.0_glfloat)

grid_newview = .true.
call draw_grid
end subroutine set_offset

!          -------------
subroutine toggle_lights(selection)
!          -------------
integer(glcint), intent(in out) :: selection

select case(selection)
   case(0)
      if (movelighton) then
         call gldisable(gl_light0)
         movelighton = .false.
      else
         call glenable(gl_light0)
         movelighton = .true.
      endif
   case(1)
      if (rightlighton) then
         call gldisable(gl_light1)
         rightlighton = .false.
      else
         call glenable(gl_light1)
         rightlighton = .true.
      endif
   case(2)
      if (leftlighton) then
         call gldisable(gl_light2)
         leftlighton = .false.
      else
         call glenable(gl_light2)
         leftlighton = .true.
      endif
   case(3)
      if (toplighton) then
         call gldisable(gl_light3)
         toplighton = .false.
      else
         call glenable(gl_light3)
         toplighton = .true.
      endif
   case(4)
      if (bottomlighton) then
         call gldisable(gl_light4)
         bottomlighton = .false.
      else
         call glenable(gl_light4)
         bottomlighton = .true.
      endif
end select

return
end subroutine toggle_lights

!          --------------
subroutine set_elem_label(selection)
!          --------------
integer(glcint), intent(in out) :: selection
label_elements = selection
call draw_grid
return
end subroutine set_elem_label

!          -------------
subroutine set_vert_label(selection)
!          -------------
integer(glcint), intent(in out) :: selection
label_verts = selection
call draw_grid
return
end subroutine set_vert_label

!          -------------
subroutine set_edge_label(selection)
!          -------------
integer(glcint), intent(in out) :: selection
label_edges = selection
call draw_grid
return
end subroutine set_edge_label

!          -------------
subroutine set_face_label(selection)
!          -------------
integer(glcint), intent(in out) :: selection
label_faces = selection
call draw_grid
return
end subroutine set_face_label

!          ----------------
subroutine set_compnt_scale(selection)
!          ----------------
integer(glcint), intent(in out) :: selection
select case(selection)
case(0)
   indiv_compnt_scale = .true.
case(1)
   indiv_compnt_scale = .false.
end select
call set_max_min
call draw_grid
return
end subroutine set_compnt_scale

!          -------------
subroutine set_assoc_elem(selection)
!          -------------
integer(glcint), intent(in out) :: selection
select case(selection)
case(1) ! toggle vertex associated element
   vert_associated_element = .not. vert_associated_element
case(2) ! toggle edge associated element
   edge_associated_element = .not. edge_associated_element
case(3) ! toggle face associated element
   face_associated_element = .not. face_associated_element
end select
call draw_grid
return
end subroutine set_assoc_elem

!          ------------
subroutine select_eigen(selection)
!          ------------
integer(glcint), intent(in out) :: selection
eigen = selection
call set_max_min
call draw_grid
return
end subroutine select_eigen

!          -------------
subroutine select_compnt(selection)
!          -------------
integer(glcint), intent(in out) :: selection
compnt = selection
call set_max_min
call draw_grid
return
end subroutine select_compnt

!          ----------------
subroutine write_postscript(selection)
!          ----------------
use rendereps
integer(glcint), intent(in out) :: selection
integer :: bsize
character(len=11) :: filename

! just a colored triangle plus triangle edges seems to require about
! 69*number of leaf elements for the buffer size.  100*total elements
! should be plenty.  If it is not, then the output postscript file
! contains just a white square.

filename = "render .eps"
write(filename(7:7),"(i1)") my_processor
bsize = 100*size(element)
call glutsetcursor(glut_cursor_wait)
select case(selection)
case(0)
   call outputeps(bsize,.true.,filename)
case(1)
   call outputeps(bsize,.false.,filename)
end select
call glutsetcursor(glut_cursor_inherit)

return
end subroutine write_postscript

!          ------------
subroutine menu_handler(selection)
!          ------------
integer(glcint), intent(in out) :: selection

select case(selection)
case(12)
   draw_axes_flag = .not. draw_axes_flag
   grid_newview = .true.
   call draw_grid
case(13)
   print *,"enter crop region as xmin, xmax, ymin, ymax, zmin, zmax"
   read *,xcrop1, xcrop2, ycrop1, ycrop2, zcrop1, zcrop2
   grid_newview = .true.
   call draw_grid
case(14)
   call flip_cart_to_cyl
   grid_newview = .true.
   call draw_grid

end select

return
end subroutine menu_handler

!          ----------------
subroutine flip_cart_to_cyl()
!          ----------------

logical :: visited_vert(grid%biggest_vert)
real(gldouble) :: lookfrom_x, lookfrom_y, lookfrom_z

convert_cylindrical_to_cartesian = .not. convert_cylindrical_to_cartesian
visited_vert = .false.
if (convert_cylindrical_to_cartesian) then
   call cylindrical_vertices_to_cartesian(REFTREE_ROOT,visited_vert)
   lookat_x = 0
   lookat_y = 0
   lookat_z = 0
   lookfrom_x = lookat_x + 3.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   lookfrom_y = lookat_y - 9.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   lookfrom_z = lookat_z + 4.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   call reset_look_from_at(lookfrom_x,lookfrom_y,lookfrom_z, &
                           lookat_x,lookat_y,lookat_z)
else
   call cartesian_vertices_to_cylindrical(REFTREE_ROOT,visited_vert)
   lookat_x = (xmax+xmin)/2
   lookat_y = (ymax+ymin)/2
   lookat_z = (zmax+zmin)/2
   lookfrom_x = lookat_x + 3.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   lookfrom_y = lookat_y - 9.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   lookfrom_z = lookat_z + 4.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   call reset_look_from_at(lookfrom_x,lookfrom_y,lookfrom_z, &
                           lookat_x,lookat_y,lookat_z)
endif
end subroutine flip_cart_to_cyl

!                    ---------------------------------
recursive subroutine cartesian_vertices_to_cylindrical(elem,visited_vert)
!                    ---------------------------------
integer, intent(in) :: elem
logical, intent(inout) :: visited_vert(:)
integer :: i, vert
real(my_real) :: x, y, rho, phi
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(my_real), parameter :: small=.001_my_real

if (elem == REFTREE_ROOT) then
   xmin = huge(0.0_my_real)
   xmax = -huge(0.0_my_real)
   ymin = huge(0.0_my_real)
   ymax = -huge(0.0_my_real)
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call cartesian_vertices_to_cylindrical(i,visited_vert)
      i = element(i)%next
   end do
   return
endif

do i=1,VERTICES_PER_ELEMENT
   vert = element(elem)%vertex(i)
   if (.not. visited_vert(vert)) then
      visited_vert(vert) = .true.
      x = vertex(vert)%coord%x
      y = vertex(vert)%coord%y
      rho = max(0.0_my_real,sqrt(x**2 + y**2)-small)
      if (x == 0.0_my_real .and. y == 0.0_my_real) then
         phi = 0.0_my_real
      else
         phi = atan2(y,x)
         if (phi < 0.0_my_real) phi = phi + 8.0_my_real*atan(1.0_my_real)
      endif
      vertex(vert)%coord%x = rho
      vertex(vert)%coord%y = phi
      xmin = min(xmin,rho)
      xmax = max(xmax,rho)
      ymin = min(ymin,phi)
      ymax = max(ymax,phi)
   endif
end do

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call cartesian_vertices_to_cylindrical(children(i),visited_vert)
   endif
end do

end subroutine cartesian_vertices_to_cylindrical

!                    ---------------------------------
recursive subroutine cylindrical_vertices_to_cartesian(elem,visited_vert)
!                    ---------------------------------
integer, intent(in) :: elem
logical, intent(inout) :: visited_vert(:)
integer :: i, vert
real(my_real) :: x, y, rho, phi
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(my_real), parameter :: small=.001_my_real

if (elem == REFTREE_ROOT) then
   xmin = huge(0.0_my_real)
   xmax = -huge(0.0_my_real)
   ymin = huge(0.0_my_real)
   ymax = -huge(0.0_my_real)
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call cylindrical_vertices_to_cartesian(i,visited_vert)
      i = element(i)%next
   end do
   return
endif

! Increment rho by a small number so we don't have the whole rho=0 plane
! mapped to the z axis, resulting in degenerate tetrahedra.
! Bound phi away from 2*pi by a small number so that the phi=2*pi plane is
! not mapped to the phi=0 plane in the inverse mapping.

do i=1,VERTICES_PER_ELEMENT
   vert = element(elem)%vertex(i)
   if (.not. visited_vert(vert)) then
      visited_vert(vert) = .true.
      rho = vertex(vert)%coord%x + small
      phi = min(vertex(vert)%coord%y,8*atan(1.0_my_real)-small)
      x = rho*cos(phi)
      y = rho*sin(phi)
      vertex(vert)%coord%x = x
      vertex(vert)%coord%y = y
      xmin = min(xmin,x)
      xmax = max(xmax,x)
      ymin = min(ymin,y)
      ymax = max(ymax,y)
   endif
end do

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call cylindrical_vertices_to_cartesian(children(i),visited_vert)
   endif
end do

end subroutine cylindrical_vertices_to_cartesian

!          ---------
subroutine menu_nada(selection)
!          ---------
integer(glcint), intent(in out) :: selection

return
end subroutine menu_nada

!          -------------
subroutine set_color_key(selection)
!          -------------
integer(glcint), intent(in out) :: selection
integer(glcint) :: xsize, ysize, scr_width
logical, save :: full_screen = .false.

xsize = glutGet(GLUT_WINDOW_WIDTH)
ysize = glutGet(GLUT_WINDOW_HEIGHT)
scr_width = glutGet(GLUT_SCREEN_WIDTH)

if (selection == NO_KEY) then
   if (show_key_flag) then
      if (.not. full_screen) then
         xsize = 2.0d0*xsize/3.0d0
      endif
      show_key_flag = .false.
      call glutreshapewindow(xsize,ysize)
      call gluttimerfunc(20_glcuint,finish_toggle_key,1_glcint)
   else
      call gluttimerfunc(20_glcuint,finish_toggle_key,0_glcint)
   endif
else
   select case (selection)
   case (EDGE_KEY)
      color_key_color = color_lines
   case (FACE_KEY)
      color_key_color = color_faces
   case (CUTTING_KEY)
      color_key_color = color_cutting_plane
   end select
   if (.not. show_key_flag) then
      xsize = 3.0d0*xsize/2.0d0
      if (xsize >= scr_width .and. scr_width /= 0) then
         xsize = scr_width
         full_screen = .true.
      else
         full_screen = .false.
      endif
      show_key_flag = .true.
      call glutreshapewindow(xsize,ysize)
      call gluttimerfunc(20_glcuint,finish_toggle_key,1_glcint)
   else
      call gluttimerfunc(20_glcuint,finish_toggle_key,0_glcint)
   endif
endif

end subroutine set_color_key

!          -----------------
subroutine finish_toggle_key(selection)
!          -----------------
integer(glcint), intent(in out) :: selection

! This routine finishes the process of toggling the color key.  It seems that
! the window is not really resized as far as glutGet is concerned until we
! return to the main loop.  So we return to the main loop with a timer set
! to call this routine after a very short amount of time.

if (selection == 1_glcint) then
   call reset_view(.true.)
else
   call reset_view
endif
grid_newview = .true.
call draw_grid

end subroutine finish_toggle_key

!          --------
subroutine draw_key
!          --------

! draw the color key

real(gldouble) :: color(4), xmin,xmax,ymin,ymax, dval, dminval, dmaxval, dnum
real(my_real) :: rminval, rmaxval, ppval, ppssval, ppssmin, ppssmax
integer :: i, iminval, imaxval, inum, nbox, nlabel, ival
integer(glcint) :: xsize, ysize, xstart, scr_width

! set the view

call gldisable(gl_lighting)
call glMatrixMode(GL_MODELVIEW)
call glpushmatrix()
call glloadidentity()
call glMatrixMode(GL_PROJECTION)
call glpushmatrix
call glloadidentity()
call glortho(0._gldouble,1.0_gldouble,0._gldouble,1._gldouble,-1._gldouble,1._gldouble)
scr_width = glutGet(GLUT_SCREEN_WIDTH)
ysize = glutGet(GLUT_WINDOW_HEIGHT)
xsize = glutGet(GLUT_WINDOW_WIDTH)
xstart = 2.0d0*xsize/3.0d0
call glViewport(xstart,0,xsize,ysize)

! determine how many color boxes to draw

select case(color_key_color)
case (COLOR_SOLUTION, COLOR_TRUE, COLOR_ERROR, COLOR_SIZE)
   nbox = 256
case (COLOR_OWNER, COLOR_VERT_OWNER, COLOR_DEGREE)
   select case(color_scheme)
   case (SCHEME_STEP_SEQ)
      if (color_key_color == COLOR_OWNER .or. &
          color_key_color == COLOR_VERT_OWNER) then
         nbox = max(2,nproc)
      else
         nbox = step_scheme_steps*step_scheme_hues
      endif
   case (SCHEME_RAINBOW, SCHEME_GRAY, SCHEME_DOUBLE_RAINBOW)
      select case(color_key_color)
      case (COLOR_OWNER, COLOR_VERT_OWNER)
         nbox = max(2,nproc)
      case (COLOR_DEGREE)
         nbox = MAX_DEGREE
      end select
   end select
case (COLOR_TRANSPARENT, COLOR_WHITE, COLOR_BLACK, COLOR_PART_BOUND)
   nbox = 0
end select

! determine the minumum and maximum values for the key

select case(color_key_color)
case(COLOR_SOLUTION, COLOR_TRUE, COLOR_ERROR)
   call preproc_and_scale(color_key_color,1.0_my_real,ppval,ppssval,rminval, &
                          ppssmin,rmaxval,ppssmax)
   dminval = rminval
   dmaxval = rmaxval
case(COLOR_SIZE)
   dminval = minsize
   dmaxval = maxsize
case(COLOR_OWNER, COLOR_VERT_OWNER)
   iminval = 1
   imaxval = max(2,nproc)
   dminval = iminval
   dmaxval = imaxval
case(COLOR_DEGREE)
   if (color_scheme == SCHEME_STEP_SEQ) then
      iminval = 1
      imaxval = step_scheme_steps*step_scheme_hues
   else
      iminval = 1
      imaxval = MAX_DEGREE
   endif
   dminval = iminval
   dmaxval = imaxval
end select

! draw the color boxes

xmin = 0.05_gldouble
xmax = 0.1_gldouble
do i=1,nbox
   if (color_key_color == COLOR_OWNER .or. &
       color_key_color == COLOR_VERT_OWNER) then
      ival = iminval + (i-1)*(imaxval-iminval)/(nbox-1)
      color = owner_color(:,ival)
   else
      dval = dminval + (i-1)*(dmaxval-dminval)/(nbox-1)
      call get_rainbow(dval,dminval,dmaxval,color)
   endif
   call glpolygonmode(gl_front_and_back, gl_fill)
   call glcolor4dv(color)
   ymin = .1_gldouble + (i-1)*.8_gldouble/nbox
   ymax = .1_gldouble +     i*.8_gldouble/nbox
   call glrectd(xmin,ymin,xmax,ymax)
   if (nbox /= 256) then
     call glpolygonmode(gl_front_and_back, gl_line)
     call glcolor4d(0.0_gldouble,0.0_gldouble,0.0_gldouble,0.0_gldouble)
     call glrectd(xmin,ymin,xmax,ymax)
     call glpolygonmode(gl_front_and_back, gl_fill)
   endif
end do

! label some of the color boxes

call glcolor4d(0.0_gldouble,0.0_gldouble,0.0_gldouble,0.0_gldouble)
select case(color_key_color)
case (COLOR_SOLUTION, COLOR_TRUE, COLOR_ERROR, COLOR_SIZE)
   xmin = 0.1875_gldouble
   nlabel = 5
case (COLOR_OWNER, COLOR_VERT_OWNER, COLOR_DEGREE)
   xmin = 0.125_gldouble
   nlabel = min(5,nbox)
case (COLOR_TRANSPARENT, COLOR_WHITE, COLOR_BLACK, COLOR_PART_BOUND)
   nlabel = 0
end select
do i=1,nlabel
   ymin = .1_gldouble + (i-1)*.8_gldouble*(nbox-1)/((nlabel-1.0d0)*nbox)
   select case(color_key_color)
   case (COLOR_SOLUTION, COLOR_TRUE, COLOR_ERROR, COLOR_SIZE)
      dnum = dminval + ((i-1)*(dmaxval-dminval))/(nlabel-1)
      if (color_key_color == COLOR_SIZE) then
         dnum = exp(dnum)
      endif
      call real_number(dnum,xmin,ymin,0.0_gldouble,0.5_gldouble)
   case (COLOR_OWNER, COLOR_VERT_OWNER, COLOR_DEGREE)
      inum = iminval + ((i-1)*(imaxval-1))/(nlabel-1)
      call number(inum,xmin,ymin,0.0_gldouble,0.5_gldouble)
   end select
end do

! restore the grid view

call glMatrixMode(GL_PROJECTION)
call glpopmatrix()
call glMatrixMode(GL_MODELVIEW)
call glpopmatrix()

end subroutine draw_key

!          ------------
subroutine init_windows
!          ------------
real(gldouble) :: lookfrom_x, lookfrom_y, lookfrom_z
integer(glcint) :: lines, face, lights, elem_label, vert_label, &
                   elem_edge, elem_face, &
                   edge_label, face_label, mark_av, color_key_menu, &
                   isosurface, isosurface_niso, &
                   postscript, sfc, off, preproc, &
                   step_scheme_params, component_scale, &
                   cutplane, cutplanecolor, edge_trans, face_trans, &
                   isosurface_trans, cutplaneedge, cutplanecontour
integer(glcint) :: menuid
character(len=32) :: str

! set the grid window view and callback functions

   call glutsetwindow(grid_win)
   call glpolygonoffset(offset,0.0_glfloat)
   lookat_x = (xmax+xmin)/2
   lookat_y = (ymax+ymin)/2
   lookat_z = (zmax+zmin)/2
   if (convert_cylindrical_to_cartesian) then
      lookat_x = 0
      lookat_y = 0
   endif
   lookfrom_x = lookat_x + 3.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   lookfrom_y = lookat_y - 9.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   lookfrom_z = lookat_z + 4.0_gldouble*max(max(xmax-xmin,ymax-ymin), &
                                            zmax-zmin)
   menuid = view_modifier_init(lookfrom_x, lookfrom_y, lookfrom_z, &
                               lookat_x,   lookat_y,   lookat_z, &
                               real(maxdomain,gldouble))

! create menu for grid window

   lines = glutcreatemenu(set_draw_lines)
   call glutaddmenuentry("no lines",COLOR_TRANSPARENT)
   call glutaddmenuentry("black",COLOR_BLACK)
   call glutaddmenuentry("edge owner",COLOR_OWNER)
   call glutaddmenuentry("vertex owner",COLOR_VERT_OWNER)
   call glutaddmenuentry("computed solution",COLOR_SOLUTION)
   call glutaddmenuentry("true solution",COLOR_TRUE)
   call glutaddmenuentry("error",COLOR_ERROR)
   call glutaddmenuentry("size",COLOR_SIZE)
   call glutaddmenuentry("degree",COLOR_DEGREE)
   edge_trans = glutcreatemenu(set_edge_transparency)
   call glutaddmenuentry("0 (opaque)",0_glcint)
   call glutaddmenuentry("25",25_glcint)
   call glutaddmenuentry("50",50_glcint)
   call glutaddmenuentry("75",75_glcint)
   call glutaddmenuentry("90",90_glcint)
   call glutaddmenuentry("95",95_glcint)
   call glutaddmenuentry("99",99_glcint)
   call glutaddmenuentry("decrease by 1",-1_glcint)
   call glutaddmenuentry("decrease by 5",-2_glcint)
   call glutaddmenuentry("increase by 1",-3_glcint)
   call glutaddmenuentry("increase by 5",-4_glcint)
   elem_edge = glutcreatemenu(menu_nada)
   call glutaddsubmenu("color",lines)
   call glutaddsubmenu("transparency",edge_trans)
   face = glutcreatemenu(set_draw_face)
   call glutaddmenuentry("transparent",COLOR_TRANSPARENT)
   call glutaddmenuentry("owner",COLOR_OWNER)
   call glutaddmenuentry("computed solution",COLOR_SOLUTION)
   call glutaddmenuentry("true solution",COLOR_TRUE)
   call glutaddmenuentry("error",COLOR_ERROR)
   call glutaddmenuentry("size",COLOR_SIZE)
   call glutaddmenuentry("degree",COLOR_DEGREE)
   face_trans = glutcreatemenu(set_face_transparency)
   call glutaddmenuentry("0 (opaque)",0_glcint)
   call glutaddmenuentry("25",25_glcint)
   call glutaddmenuentry("50",50_glcint)
   call glutaddmenuentry("75",75_glcint)
   call glutaddmenuentry("90",90_glcint)
   call glutaddmenuentry("95",95_glcint)
   call glutaddmenuentry("99",99_glcint)
   call glutaddmenuentry("decrease by 1",-1_glcint)
   call glutaddmenuentry("decrease by 5",-2_glcint)
   call glutaddmenuentry("increase by 1",-3_glcint)
   call glutaddmenuentry("increase by 5",-4_glcint)
   elem_face = glutcreatemenu(menu_nada)
   call glutaddsubmenu("color",face)
   call glutaddsubmenu("transparency",face_trans)
   cutplanecolor = glutcreatemenu(set_cutting_plane_color)
   call glutaddmenuentry("computed solution",COLOR_SOLUTION)
   call glutaddmenuentry("true solution",COLOR_TRUE)
   call glutaddmenuentry("error",COLOR_ERROR)
   cutplaneedge = glutcreatemenu(set_cutting_plane_edge)
   call glutaddmenuentry("no lines",COLOR_TRANSPARENT)
   call glutaddmenuentry("black",COLOR_BLACK)
   cutplanecontour = glutcreatemenu(set_contour_param)
   call glutaddmenuentry("toggle contour plot",9_glcint)
   call glutaddmenuentry("increment number by 1",3_glcint)
   call glutaddmenuentry("decrement number by 1",4_glcint)
   call glutaddmenuentry("increment number by 10",5_glcint)
   call glutaddmenuentry("decrement number by 10",6_glcint)
   call glutaddmenuentry("double number",7_glcint)
   call glutaddmenuentry("cut number in half",8_glcint)
   call glutaddmenuentry("enter number in debug window",1_glcint)
   call glutaddmenuentry("enter values in debug window",2_glcint)
   cutplane = glutcreatemenu(set_cutting_plane)
   call glutaddmenuentry("toggle xz cutting plane",XZ_PLANE)
   call glutaddmenuentry("toggle yz cutting plane",YZ_PLANE)
   call glutaddmenuentry("toggle xy cutting plane",XY_PLANE)
   call glutaddsubmenu("cutting plane color",cutplanecolor)
   call glutaddsubmenu("edge color",cutplaneedge)
   call glutaddsubmenu("contour plot",cutplanecontour)
   call glutaddmenuentry("speed up movement",1_glcint)
   call glutaddmenuentry("slow down movement",2_glcint)
   isosurface_trans = glutcreatemenu(set_isosurface_transparency)
   call glutaddmenuentry("0 (opaque)",0_glcint)
   call glutaddmenuentry("25",25_glcint)
   call glutaddmenuentry("50",50_glcint)
   call glutaddmenuentry("75",75_glcint)
   call glutaddmenuentry("90",90_glcint)
   call glutaddmenuentry("95",95_glcint)
   call glutaddmenuentry("99",99_glcint)
   call glutaddmenuentry("decrease by 1",-1_glcint)
   call glutaddmenuentry("decrease by 5",-2_glcint)
   call glutaddmenuentry("increase by 1",-3_glcint)
   call glutaddmenuentry("increase by 5",-4_glcint)
   isosurface_niso = glutcreatemenu(set_isosurface_param)
   call glutaddmenuentry("decrement by 1",4_glcint)
   call glutaddmenuentry("increment by 1",3_glcint)
   call glutaddmenuentry("decrement by 10",6_glcint)
   call glutaddmenuentry("increment by 10",5_glcint)
   call glutaddmenuentry("cut in half",8_glcint)
   call glutaddmenuentry("double",7_glcint)
   isosurface = glutcreatemenu(set_isosurface)
   call glutaddmenuentry("no isosurface plot",DRAW_NO_FUNCTION)
   call glutaddmenuentry("computed solution",DRAW_SOLUTION)
   call glutaddmenuentry("true solution",DRAW_TRUE)
   call glutaddmenuentry("error",DRAW_ERROR)
   call glutaddsubmenu("isosurface transparency",isosurface_trans)
   call glutaddsubmenu("number of isosurfaces",isosurface_niso)
   preproc = glutcreatemenu(set_preproc_func)
   call glutaddmenuentry("none",PREPROC_NONE)
   call glutaddmenuentry("-f",PREPROC_NEG)
   call glutaddmenuentry("abs(f)",PREPROC_ABS)
   call glutaddmenuentry("f**2",PREPROC_SQ)
   call glutaddmenuentry("log(abs(f))",PREPROC_LOG)
!   lights = glutcreatemenu(toggle_lights)
!   call glutaddmenuentry("movable light",0_glcint)
!   call glutaddmenuentry("right light",1_glcint)
!   call glutaddmenuentry("left light",2_glcint)
!   call glutaddmenuentry("top light",3_glcint)
!   call glutaddmenuentry("bottom light",4_glcint)
   elem_label = glutcreatemenu(set_elem_label)
   call glutaddmenuentry("none",LABEL_NOLABEL)
   call glutaddmenuentry("local id",LABEL_LID)
   call glutaddmenuentry("global id",LABEL_GID)
   vert_label = glutcreatemenu(set_vert_label)
   call glutaddmenuentry("none",LABEL_NOLABEL)
   call glutaddmenuentry("local id",LABEL_LID)
   call glutaddmenuentry("global id",LABEL_GID)
   edge_label = glutcreatemenu(set_edge_label)
   call glutaddmenuentry("none",LABEL_NOLABEL)
   call glutaddmenuentry("local id",LABEL_LID)
   call glutaddmenuentry("global id",LABEL_GID)
   face_label = glutcreatemenu(set_face_label)
   call glutaddmenuentry("none",LABEL_NOLABEL)
   call glutaddmenuentry("local id",LABEL_LID)
   call glutaddmenuentry("global id",LABEL_GID)
   step_scheme_params = glutcreatemenu(set_step_scheme)
   call glutaddmenuentry("steps=5",1_glcint)
   call glutaddmenuentry("steps=10",2_glcint)
   call glutaddmenuentry("steps=16",3_glcint)
   call glutaddmenuentry("decrement steps",4_glcint)
   call glutaddmenuentry("increment steps",5_glcint)
   call glutaddmenuentry("hues=5",6_glcint)
   call glutaddmenuentry("hues=10",7_glcint)
   call glutaddmenuentry("hues=16",8_glcint)
   call glutaddmenuentry("decrement hues",9_glcint)
   call glutaddmenuentry("increment hues",10_glcint)
   scheme_menu = glutcreatemenu(set_color_scheme)
   call glutaddmenuentry("rainbow",SCHEME_RAINBOW)
   call glutaddmenuentry("double rainbow",SCHEME_DOUBLE_RAINBOW)
   call glutaddmenuentry("gray scale",SCHEME_GRAY)
! if this is changed so that stepped sequential is not the 5th entry, need to also
! make a change in set_step_scheme
   write(str,"(A18,2I4)") "stepped sequential",step_scheme_steps,step_scheme_hues
   call glutaddmenuentry(trim(str),SCHEME_STEP_SEQ)
   call glutaddsubmenu("change stepped sequential",step_scheme_params)
   mark_av = glutcreatemenu(set_assoc_elem)
   call glutaddmenuentry("toggle vertex associated element",1_glcint)
   call glutaddmenuentry("toggle edge associated element",2_glcint)
   call glutaddmenuentry("toggle face associated element",3_glcint)
   color_key_menu = glutcreatemenu(set_color_key)
   call glutaddmenuentry("no color key",NO_KEY)
   call glutaddmenuentry("key for edge color",EDGE_KEY)
   call glutaddmenuentry("key for face color",FACE_KEY)
   call glutaddmenuentry("key for cutting plane",CUTTING_KEY)
   postscript = glutcreatemenu(write_postscript)
   call glutaddmenuentry("Write out Encapsulated PS (sorted)",0_glcint)
   call glutaddmenuentry("Write out Encapsulated PS (unsorted)",1_glcint)
!   sfc = glutcreatemenu(set_sfc)
!   call glutaddmenuentry("toggle space filling curve",0_glcint)
!   call glutaddmenuentry("toggle in/out vertex label",1_glcint)
   off = glutcreatemenu(set_offset)
   call glutaddmenuentry("decrease by 10",-10_glcint)
   call glutaddmenuentry("decrease by 1",-1_glcint)
   call glutaddmenuentry("increase by 1",1_glcint)
   call glutaddmenuentry("increase by 10",10_glcint)
   eigen_menu = glutcreatemenu(select_eigen)
   menu_neigen = 1
   write(str,"(a14,i3)") "eigenfunction ",1
   call glutaddmenuentry(trim(str),1_glcint)
   component_menu = glutcreatemenu(select_compnt)
   menu_ncompnt = 3
   write(str,"(a6)") "L1 sum"
   call glutaddmenuentry(trim(str),-1_glcint)
   write(str,"(a6)") "L2 sum"
   call glutaddmenuentry(trim(str),-2_glcint)
   write(str,"(a10,i3)") "component ",1
   call glutaddmenuentry(trim(str),1_glcint)
   component_scale = glutcreatemenu(set_compnt_scale)
   call glutaddmenuentry("individual",0_glcint)
   call glutaddmenuentry("all the same",1_glcint)
   grid_menu = glutcreatemenu(menu_handler)
   call glutaddsubmenu("view modifier",menuid)
   call glutaddsubmenu("element edge",elem_edge)
   call glutaddsubmenu("element face",elem_face)
   call glutaddsubmenu("cutting planes",cutplane)
   call glutaddsubmenu("isosurface",isosurface)
   call glutaddsubmenu("preprocess function",preproc)
   call glutaddsubmenu("color scheme",scheme_menu)
!   call glutaddsubmenu("toggle lights",lights)
   call glutaddsubmenu("element label",elem_label)
   call glutaddsubmenu("face label", face_label)
   call glutaddsubmenu("edge label", edge_label)
   call glutaddsubmenu("vertex label", vert_label)
   call glutaddsubmenu("associated element", mark_av)
   call glutaddsubmenu("eigenfunction to use",eigen_menu)
   call glutaddsubmenu("component to use",component_menu)
   call glutaddsubmenu("component scale",component_scale)
!   call glutaddsubmenu("space filling curve",sfc)
   call glutaddsubmenu("grid offset",off)
!   call glutaddmenuentry("crop (debug window)",13_glcint)
   call glutaddmenuentry("toggle axes",12_glcint)
   call glutaddsubmenu("color key",color_key_menu)
   call glutaddsubmenu("write postscript",postscript)
   call glutattachmenu(glut_right_button)
! secret conversion of cylindrical coordinates to cartesian
   call glutaddmenuentry("toggle cylindrical to cartesian",14_glcint)

! set up lighting conditions for grid

   call glclearcolor(1.0_glclampf, 1.0_glclampf, 1.0_glclampf, 1.0_glclampf)
   call gllightfv(gl_light0, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light1, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light2, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light3, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light4, gl_diffuse, (/glfone,glfone,glfone,glfone/))
!   call gllightfv(gl_light1, gl_position, (/-2.5, -.5, -2.0, 0.0/))
   call gllightfv(gl_light1, gl_position, (/-7.1268,1.68,-6.81, 0.0/))
   call gllightfv(gl_light2, gl_position, (/1.5,-.5,-2.,0./))
   call gllightfv(gl_light3, gl_position, (/-.5,-.5,-2.,0./))
   call gllightfv(gl_light4, gl_position, (/-.5,-.5, 2.,0./))
   call glenable(gl_lighting)
   if (movelighton) call glenable(gl_light0)
   if (rightlighton) call glenable(gl_light1)
   if (leftlighton) call glenable(gl_light2)
   if (toplighton) call glenable(gl_light3)
   if (bottomlighton) call glenable(gl_light4)
   call gllightmodelfv(gl_light_model_ambient, (/glfhalf,glfhalf,glfhalf,glfone/))
   call gldepthfunc(gl_lequal)
   call glenable(gl_depth_test)

   window_initialized = .true.

end subroutine init_windows

! some compiler doesn't let timer pass itself, so use 2 timers
! that pass each other

!                    -----
recursive subroutine timer(selection)
!                    -----
integer(glcint), intent(in out) :: selection
integer :: i

i = selection ! just to shut up compilers that warn about unused arguments

call process_message
call gluttimerfunc(10_glcuint,timer2,0_glcint)

return
end subroutine timer

!                    -----
recursive subroutine timer2(selection)
!                    -----
integer(glcint), intent(in out) :: selection
integer :: i

i = selection ! just to shut up compilers that warn about unused arguments

call process_message
call gluttimerfunc(10_glcuint,timer,0_glcint)

return
end subroutine timer2

end module graphics_mod

!-------------------------------------------
! The main program
!-------------------------------------------

!          --------------
subroutine phaml_graphics
!          --------------
use graphics_mod
use opengl_gl
use opengl_glut
implicit none

integer :: spawn_form=0, allocstat, system_size, eq_type, ijunk
character(len=32) :: dummy_char
logical :: junk1,junk2,update_umod
type(phaml_solution_type) :: phaml_solution
integer :: proc,ni,nr
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
!----------------------------------------------------
! Begin executable code

! nullify grid data structures pointers so they can be tested for associated

nullify (grid%element,grid%edge,grid%vertex)

! initialize communication

junk1 = .false.
junk2 = .false.
allocate(procs,stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for procs",intlist=(/allocstat/))
   stop
endif
call init_comm(procs,spawn_form,junk1,junk2,dummy_char,outunit, &
               errunit,system_size,eq_type,my_pde_id,nproc,ijunk,grid%max_blen,&
               grid%triangle_files,update_umod)
if (update_umod) then
   if (PARALLEL /= SEQUENTIAL) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,101)
      if (recv_int(1) /= 9) then
         call fatal("In graphics start up, received code that is not update_usermod.")
         stop
      endif
      phaml_solution%procs = procs
      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)
      call update_usermod(phaml_solution)
   endif
endif

! use unit 6 for all output from the graphics process

outunit = 6
errunit = 6

! initialize GLUT

call glutinit
call glutinitdisplaymode(ior(ior(glut_double,glut_rgb),glut_depth))
call glutinitwindowsize(500_glcint,500_glcint)
call gluttimerfunc(10_glcuint,timer,0_glcint)

! create the window

grid_win = glutcreatewindow("PHAML grid")
call glutdisplayfunc(grid_display)

! go into the GLUT infinite loop

call glutmainloop

end subroutine phaml_graphics
