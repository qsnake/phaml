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

module grid_io

!----------------------------------------------------
! This module contains printed and graphical i/o for grids
!
! communication tags in this module are of the form 7xx
!----------------------------------------------------

use global
use message_passing
use hash_mod
use gridtype_mod
use grid_util
use stopwatch
use zoltan_interf
use error_estimators
use sort_mod
!----------------------------------------------------

implicit none
private
public draw_grid, print_element, print_vertex, print_edge, print_face, &
       print_all_elements, print_all_vertices, print_all_edges, &
       print_all_faces, print_entire_grid, print_grid_info, store_grid

contains

!          ------------------
subroutine print_all_elements(grid)
!          ------------------

!----------------------------------------------------
! This subroutine prints the data structure for all the elements
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables

integer :: lev, elem

!----------------------------------------------------
! Begin executable code

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      call print_element(grid,elem)
      elem = grid%element(elem)%next
   end do
end do

end subroutine print_all_elements

!          ---------------
subroutine print_all_faces(grid)
!          ---------------

!----------------------------------------------------
! This subroutine prints the data structure for all the faces
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables

integer :: lev, elem, face, i
logical :: printed(grid%biggest_face)

!----------------------------------------------------
! Begin executable code

! go through elements and print any faces that haven't already been printed

printed = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,FACES_PER_ELEMENT
         face = grid%element(elem)%face(i)
         if (.not. printed(face)) then
            call print_face(grid,face)
            printed(face) = .true.
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

end subroutine print_all_faces

!          ---------------
subroutine print_all_edges(grid)
!          ---------------

!----------------------------------------------------
! This subroutine prints the data structure for all the edges
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables

integer :: lev, elem, edge, i
logical :: printed(grid%biggest_edge)

!----------------------------------------------------
! Begin executable code

! go through elements and print any edges that haven't already been printed

printed = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)
         if (.not. printed(edge)) then
            call print_edge(grid,edge)
            printed(edge) = .true.
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

end subroutine print_all_edges

!          ------------------
subroutine print_all_vertices(grid)
!          ------------------

!----------------------------------------------------
! This subroutine prints the data structure for all the vertices
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables

integer lev, vert, elem, ivert
logical :: visited_vert(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

visited_vert = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do ivert=1,VERTICES_PER_ELEMENT
         vert = grid%element(elem)%vertex(ivert)
         if (visited_vert(vert)) cycle
         visited_vert(vert) = .true.
         call print_vertex(grid,vert)
      end do
      elem = grid%element(elem)%next
   end do
end do

end subroutine print_all_vertices

!          -----------------
subroutine print_entire_grid(grid)
!          -----------------

!----------------------------------------------------
! This subroutine prints the data structures for the whole grid
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Begin executable code

call print_all_vertices(grid)
call print_all_edges(grid)
if (global_element_kind == TETRAHEDRAL_ELEMENT) call print_all_faces(grid)
call print_all_elements(grid)

end subroutine print_entire_grid

!          -------------
subroutine print_element(grid,elem)
!          -------------

!----------------------------------------------------
! This subroutine prints the data structure for element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: elem

!----------------------------------------------------
! Begin executable code

write(outunit,"(A)") '--------------------'
write(outunit,"(A,I11)") 'Element number ',elem
write(outunit,"(A)") '--------------------'
call hash_print_key(grid%element(elem)%gid,outunit,"Global ID: ")
write(outunit,"(A,4I11)") 'Vertices:  ',grid%element(elem)%vertex
write(outunit,"(A,6I11)") 'Edges:     ',grid%element(elem)%edge
if (global_element_kind == TETRAHEDRAL_ELEMENT) then
   write(outunit,"(A,4I11)") 'Faces:     ',grid%element(elem)%face
   call hash_print_key(grid%element(elem)%neighbor_hint,outunit,"Neighbor hint: ")
   write(outunit,"(A,I11)") 'Refinement edge: ',grid%element(elem)%refinement_edge
   write(outunit,"(2A)") 'Type: ',grid%element(elem)%type
end if
if (global_element_kind == TRIANGULAR_ELEMENT) then
   call hash_print_key(grid%element(elem)%mate,outunit,"Mate GID: ")
endif
write(outunit,"(A,I11)") 'Level:     ',grid%element(elem)%level
write(outunit,"(A,I11)") 'Degree:    ',grid%element(elem)%degree
if (associated(grid%element(elem)%solution)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Solution:  ',grid%element(elem)%solution
else
   write(outunit,"(A)") 'Solution:  disassociated'
endif
if (associated(grid%element(elem)%exact)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Exact:     ',grid%element(elem)%exact
else
   write(outunit,"(A)") 'Exact:     disassociated'
endif
if (associated(grid%element(elem)%oldsoln)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Oldsoln:   ',grid%element(elem)%oldsoln
else
   write(outunit,"(A)") 'Oldsoln:   disassociated'
endif
write(outunit,"(A,L1)") 'Owner:     ',grid%element(elem)%iown
write(outunit,"(A,2I11)") 'Next/Prev: ',grid%element(elem)%next,grid%element(elem)%previous
write(outunit,"(A,2I11)") 'In/Out:    ',grid%element(elem)%in,grid%element(elem)%out
if (grid%element(elem)%level == 1) then
   write(outunit,"(A,4I11)") 'Neighbors: ',grid%initial_neighbor(:,elem)
endif
write(outunit,"(A,100I11)") 'Tags:      ',grid%element(elem)%tags

end subroutine print_element

!          ----------
subroutine print_face(grid,face)
!          ----------

!----------------------------------------------------
! This subroutine prints the data structure for face face
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: face

!----------------------------------------------------
! Begin executable code

write(outunit,"(A)") '------------------'
write(outunit,"(A,I11)") 'Face number ',face
write(outunit,"(A)") '------------------'
call hash_print_key(grid%face(face)%gid,outunit,"Global ID: ")
write(outunit,"(A,3I11)") 'Vertices:           ',grid%face(face)%vertex
write(outunit,"(A,3I11)") 'Edges:              ',grid%face(face)%edge
if (associated(grid%face(face)%solution)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Solution:           ',grid%face(face)%solution
else
   write(outunit,"(A)") 'Solution:           disassociated'
endif
if (associated(grid%face(face)%exact)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Exact:              ',grid%face(face)%exact
else
   write(outunit,"(A)") 'Exact:              disassociated'
endif
if (associated(grid%face(face)%oldsoln)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Oldsoln:            ',grid%face(face)%oldsoln
else
   write(outunit,"(A)") 'Oldsoln:            disassociated'
endif
write(outunit,"(A,I11)") 'Bmark:              ',grid%face(face)%bmark
write(outunit,"(A,100I11)") 'Type:               ',grid%face_type(face,:)
write(outunit,"(A,I11)") 'Degree:             ',grid%face(face)%degree
write(outunit,"(A,I11)") 'Associated element: ',grid%face(face)%assoc_elem
write(outunit,"(A,I11)") 'Marked edge:        ',grid%face(face)%marked_edge
write(outunit,"(A,I11)") 'Next:               ',grid%face(face)%next
write(outunit,"(A,100I11)") 'Tags:               ',grid%face(face)%tags

end subroutine print_face

!          ----------
subroutine print_edge(grid,edge)
!          ----------

!----------------------------------------------------
! This subroutine prints the data structure for edge edge
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: edge

!----------------------------------------------------
! Begin executable code

write(outunit,"(A)") '------------------'
write(outunit,"(A,I11)") 'Edge number ',edge
write(outunit,"(A)") '------------------'
call hash_print_key(grid%edge(edge)%gid,outunit,"Global ID: ")
write(outunit,"(A,2I11)") 'Vertices:           ',grid%edge(edge)%vertex
if (associated(grid%edge(edge)%solution)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Solution:           ',grid%edge(edge)%solution
else
   write(outunit,"(A)") 'Solution:           disassociated'
endif
if (associated(grid%edge(edge)%exact)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Exact:              ',grid%edge(edge)%exact
else
   write(outunit,"(A)") 'Exact:              disassociated'
endif
if (associated(grid%edge(edge)%oldsoln)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Oldsoln:            ',grid%edge(edge)%oldsoln
else
   write(outunit,"(A)") 'Oldsoln:            disassociated'
endif
write(outunit,"(A,I11)") 'Bmark:              ',grid%edge(edge)%bmark
write(outunit,"(A,100I11)") 'Type:               ',grid%edge_type(edge,:)
write(outunit,"(A,I11)") 'Degree:             ',grid%edge(edge)%degree
write(outunit,"(A,I11)") 'Associated element: ',grid%edge(edge)%assoc_elem
write(outunit,"(A,I11)") 'Next:               ',grid%edge(edge)%next
write(outunit,"(A,100I11)") 'Tags:               ',grid%edge(edge)%tags

end subroutine print_edge

!          ------------
subroutine print_vertex(grid,vert)
!          ------------

!----------------------------------------------------
! This subroutine prints the data structure for vertex vertex
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: vert

!----------------------------------------------------
! Begin executable code

write(outunit,"(A)") '--------------------'
write(outunit,"(A,I11)") 'Vertex number ',vert
write(outunit,"(A)") '--------------------'
call hash_print_key(grid%vertex(vert)%gid,outunit,"Global ID: ")
write(outunit,"(SS,1P,A,3E18.10E2)") 'Coordinates:        ',grid%vertex(vert)%coord
write(outunit,"(SS,1P,A,100E18.10E2)") 'Solution:           ',grid%vertex_solution(vert,:,:)
if (associated(grid%vertex_exact)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Exact:              ',grid%vertex_exact(vert,:,:)
else
   write(outunit,"(A)") 'Exact:              disassociated'
endif
if (associated(grid%vertex_oldsoln)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Oldsoln:            ',grid%vertex_oldsoln(vert,:,:)
else
   write(outunit,"(A)") 'Oldsoln:            disassociated'
endif
write(outunit,"(A,I11)") 'Bmark:              ',grid%vertex(vert)%bmark
write(outunit,"(A,100I11)") 'Type:               ',grid%vertex_type(vert,:)
write(outunit,"(A,I11)") 'Associated element: ',grid%vertex(vert)%assoc_elem
write(outunit,"(A,I11)") 'Next:               ',grid%vertex(vert)%next
write(outunit,"(A,100I11)") 'Tags:               ',grid%vertex(vert)%tags

return
end subroutine print_vertex

!          ---------------
subroutine print_grid_info(grid,procs,io_control,still_sequential,this_time, &
                           tag)
!          ---------------

!----------------------------------------------------
! This routine prints information about the grid
!----------------------------------------------------
 
implicit none
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
type (io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
integer, intent(in) :: this_time(:),tag
!----------------------------------------------------
 
!----------------------------------------------------
! Local variables:
 
integer :: i,proc,np,who,when,astat
integer :: nelem,nlev,nvert,nvert_own,nelem_own,nelem_leaf,nelem_leaf_own, &
           mindeg,maxdeg,dof,dof_own
integer, allocatable :: ivert(:),ielem(:),ilev(:),ielem_leaf(:), &
                        ivert_own(:),ielem_own(:),ielem_leaf_own(:), &
                        imindeg(:),imaxdeg(:),idof(:),idof_own(:)
integer :: send_int(11),ni,nr
real (my_real) :: no_reals(1)
integer, pointer :: recv_int(:)
real (my_real), pointer :: recv_real(:)
 
!----------------------------------------------------
 
!----------------------------------------------------
! Begin executable code
 
! If this is not the right time to print, return

who = io_control%print_grid_who
when = io_control%print_grid_when

if (.not. any(this_time == when) .or. who == NO_ONE) return

! If I'm the master and only slaves print the grid, return

if (my_proc(procs) == MASTER .and. who == SLAVES) return

! stop the clocks

call pause_watch(all_watches)

! If the master will print, wait for the master to get here so the slaves
! don't run away from it

if (who == MASTER .or. who == MASTER_ALL .or. who == EVERYONE) then
   if (my_proc(procs) == MASTER) then
      do proc=1,num_proc(procs)
         call phaml_send(procs,proc,(/1/),1,no_reals,0,711)
      end do
   else
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,711)
      if (associated(recv_int)) deallocate(recv_int)
   endif
endif

! If I'm not the master, collect my grid info

if (my_proc(procs) /= MASTER) then
   call get_grid_info(grid,procs,still_sequential,710,nelem,nlev,nvert, &
                      nvert_own,nelem_own,nelem_leaf,nelem_leaf_own, &
                      dof,dof_own,mindeg=mindeg,maxdeg=maxdeg,no_master=.true.)
endif

! If the master will print, get it the info

if (who == MASTER .or. who == MASTER_ALL .or. who == EVERYONE) then
   if (my_proc(procs) /= MASTER) then
      send_int(1) = nvert
      send_int(2) = nelem
      send_int(3) = nlev
      send_int(4) = nelem_leaf
      send_int(5) = nvert_own
      send_int(6) = nelem_own
      send_int(7) = nelem_leaf_own
      send_int(8) = mindeg
      send_int(9) = maxdeg
      send_int(10) = dof
      send_int(11) = dof_own
      ni = 11
      nr = 0
      call phaml_send(procs,MASTER,send_int,ni,no_reals,nr,tag)
   else
      np = num_proc(procs)
      allocate(ivert(np),ielem(np),ilev(np),ielem_leaf(np),ivert_own(np), &
               ielem_own(np),ielem_leaf_own(np),imindeg(np),imaxdeg(np), &
               idof(np),idof_own(np),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in print_grid_info",procs=procs)
         return
      endif
      do i=1,np
         call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag)
         ivert(proc) = recv_int(1)
         ielem(proc) = recv_int(2)
         ilev(proc) = recv_int(3)
         ielem_leaf(proc) = recv_int(4)
         ivert_own(proc) = recv_int(5)
         ielem_own(proc) = recv_int(6)
         ielem_leaf_own(proc) = recv_int(7)
         imindeg(proc) = recv_int(8)
         imaxdeg(proc) = recv_int(9)
         idof(proc) = recv_int(10)
         idof_own(proc) = recv_int(11)
         deallocate(recv_int,stat=astat)
      end do
      if (still_sequential) then
         nvert = ivert_own(1)
         nelem = ielem_own(1)
         nlev = ilev(1)
         nelem_leaf = ielem_leaf_own(1)
         mindeg = imindeg(1)
         maxdeg = imaxdeg(1)
         dof = idof(1)
      else
         nvert = sum(ivert_own)
         nelem = sum(ielem_own)
         nlev = maxval(ilev)
         nelem_leaf = sum(ielem_leaf_own)
         mindeg = minval(imindeg)
         maxdeg = maxval(imaxdeg)
         dof = sum(idof_own)
      endif
   endif
endif

! print the info, if requested

if (my_proc(procs) /= MASTER) then
 if (who == SLAVES .or. who == EVERYONE) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'My Grid:'
   write(outunit,"(A,I11)") '   number of vertices       = ',nvert
   write(outunit,"(A,I11)") '   number of vertices I own = ',nvert_own
   write(outunit,"(A,I11)") '   number of elements       = ',nelem
   write(outunit,"(A,I11)") '   number of leaf elements  = ',nelem_leaf
   write(outunit,"(A,I11)") '   leaf elements I own      = ',nelem_leaf_own
   write(outunit,"(A,I11)") '   number of levels         = ',nlev
   write(outunit,"(A,I11)") '   degrees of freedom       = ',dof
   write(outunit,"(A,I11)") '   degrees of freedom I own = ',dof_own
   write(outunit,"(A,2I11)") '   min and max degree       = ',mindeg,maxdeg
 endif
else
 if (who == MASTER_ALL) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'Individual Grids:'
   write(outunit,"(A,128I11)") '   number of vertices       = ',ivert
   write(outunit,"(A,128I11)") '   number of vertices I own = ',ivert_own
   write(outunit,"(A,128I11)") '   number of elements       = ',ielem
   write(outunit,"(A,128I11)") '   number of leaf elements  = ',ielem_leaf
   write(outunit,"(A,128I11)") '   leaf elements I own      = ',ielem_leaf_own
   write(outunit,"(A,128I11)") '   number of levels         = ',ilev
   write(outunit,"(A,128I11)") '   degrees of freedom       = ',idof
   write(outunit,"(A,128I11)") '   degrees of freedom I own = ',idof_own
   write(outunit,"(A,128I11)") '   min degree               = ',imindeg
   write(outunit,"(A,128I11)") '   max degree               = ',imaxdeg
 endif
 if (who == MASTER .or. who == EVERYONE .or. who == MASTER_ALL) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'Total Grid:'
   write(outunit,"(A,I11)") '   number of vertices      = ',nvert
   write(outunit,"(A,I11)") '   number of leaf elements = ',nelem_leaf
   write(outunit,"(A,I11)") '   number of levels        = ',nlev
   write(outunit,"(A,I11)") '   degrees of freedom      = ',dof
   write(outunit,"(A,2I11)") '   min and max degree      = ',mindeg,maxdeg
 endif
endif

if (allocated(ivert)) then
   deallocate(ivert,ielem,ilev,ielem_leaf,ivert_own,ielem_own,ielem_leaf_own, &
              imindeg,imaxdeg,idof,idof_own,stat=astat)
endif

call end_pause_watch(all_watches)

return
end subroutine print_grid_info

!          ---------
subroutine draw_grid(grid,procs,io_control,ref_control,i_draw,master_draws, &
                     still_sequential,this_time,partition_method,lb)
!          ---------

!----------------------------------------------------
! This subroutine sends a message to the graphics process to redraw the grid.
! If the grid data has changed, it also sends the new data.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type (io_options), intent(in) :: io_control
type (refine_options), intent(in) :: ref_control
logical, intent(in) :: i_draw, master_draws, still_sequential
integer, intent(in) :: this_time(:)
integer, intent(in) :: partition_method
type(Zoltan_Struct), pointer :: lb
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: ni, nr, astat, markers
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
integer :: imess1(1)
real (my_real) :: noreals(1)
integer, allocatable :: hold_next(:), hold_prev(:)
integer :: hold_head
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! see if I should do anything

if (.not. (i_draw .or. master_draws)) return
if (.not. any(this_time == io_control%draw_grid_when)) return

call pause_watch(all_watches)

! slave

if (my_proc(procs) /= MASTER) then

! if the grid has changed since last sent to the graphics process ...

   if (grid_changed) then

! if using ZOLTAN_REFTREE,  preserve the level-by-level element linked list
! and get the child order from Zoltan

      if (partition_method == ZOLTAN_REFTREE) then

         allocate(hold_next(grid%biggest_elem), &
                  hold_prev(grid%biggest_elem),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in draw_grid",procs=procs)
            return
         endif
         hold_next = grid%element(1:grid%biggest_elem)%next
         hold_prev = grid%element(1:grid%biggest_elem)%previous
         hold_head = grid%head_level_elem(1)

         call child_order_from_zoltan(grid,procs,lb,partition_method, &
                                      still_sequential)
      endif

! if slaves are plotting, send the grid to the graphics process

      select case (global_element_kind)
      case (TRIANGULAR_ELEMENT)
         markers = 3
      case (TETRAHEDRAL_ELEMENT)
         markers = 4
      end select
      if (i_draw) then
         allocate(imess(markers+scalar_imess_size(grid) + &
                        elem_imess_size(grid)*grid%nelem + &
                        face_imess_size(grid)*grid%nface + &
                        edge_imess_size(grid)*grid%nedge + &
                        vert_imess_size(grid)*grid%nvert), &
                  rmess(scalar_rmess_size(grid) + &
                        elem_rmess_size(grid)*grid%nelem + &
                        face_rmess_size(grid)*grid%nface + &
                        edge_rmess_size(grid)*grid%nedge + &
                        vert_rmess_size(grid)*grid%nvert), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in draw_grid",procs=procs)
            return
         endif
         call pack_grid(grid,procs,ref_control,.false.,still_sequential, &
                        imess,rmess,ni,nr)
         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess(1:ni),ni,rmess(1:nr),nr)
         else
            call phaml_send(procs,graphics_proc(procs),imess,ni,rmess,nr,101)
         endif
         deallocate(imess,rmess,stat=astat)
      endif

! if master is plotting, send the grid to the master

      if (master_draws) then
         allocate(imess(markers+scalar_imess_size(grid) + &
                        elem_imess_size(grid)*grid%nelem + &
                        face_imess_size(grid)*grid%nface + &
                        edge_imess_size(grid)*grid%nedge + &
                        vert_imess_size(grid)*grid%nvert), &
                  rmess(scalar_rmess_size(grid) + &
                        elem_rmess_size(grid)*grid%nelem + &
                        face_rmess_size(grid)*grid%nface + &
                        edge_rmess_size(grid)*grid%nedge + &
                        vert_rmess_size(grid)*grid%nvert), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in draw_grid",procs=procs)
            return
         endif
         call pack_grid(grid,procs,ref_control,.true.,still_sequential, &
                        imess,rmess,ni,nr)
         if (.not.still_sequential .or. my_proc(procs) == 1) then
            call phaml_send(procs,MASTER,imess,ni,rmess,nr,720)
         endif
         deallocate(imess,rmess,stat=astat)
      endif

! if using ZOLTAN_REFTREE, restore the level-by-level element linked list

      if (partition_method == ZOLTAN_REFTREE) then
         grid%element(1:grid%biggest_elem)%next = hold_next
         grid%element(1:grid%biggest_elem)%previous = hold_prev
         grid%head_level_elem(1) = hold_head
         deallocate(hold_next, hold_prev, stat=astat)
      endif

! if the grid has not changed, just send a message to redraw

   else

      if (i_draw) then
         imess1(1) = GRAPHICS_GRID
         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess1,1,noreals,0)
         else
            call phaml_send(procs,graphics_proc(procs),imess1,1,noreals,0,101)
         endif
      endif

   endif

! master

else

   if (master_draws) then

! if the grid has changed since last sent to the graphics process ...

      if (grid_changed) then

! receive the individual grids from all the slaves and build the composite grid

         call merge_grids(grid,procs,imess,ni,rmess,nr,still_sequential)

! send the composite grid to the master's graphics process

         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess(1:ni),ni,rmess(1:nr),nr)
         else
            call phaml_send(procs,graphics_proc(procs),imess,ni,rmess,nr,101)
         endif
         deallocate(imess,rmess,stat=astat)

! if the grid has not changed, just send a message to redraw

      else

         imess1(1) = GRAPHICS_GRID
         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess1,1,noreals,0)
         else
            call phaml_send(procs,graphics_proc(procs),imess1,1,noreals,0,101)
         endif

      endif
   endif

endif ! slave or master

grid_changed = .false. ! says the graphics data has been sent since last change

! pause, if requested

call pause_until_enter(procs,io_control%pause_after_draw,dont_pausewatch=.true.)

call end_pause_watch(all_watches)

end subroutine draw_grid

!          ---------
subroutine pack_grid(grid,procs,ref_control,for_master,still_sequential, &
                     imess,rmess,ni,nr)
!          ---------

!----------------------------------------------------
! This subroutine packs the message containing the grid and flags for
! the graphics processor
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type (refine_options), intent(in) :: ref_control
logical, intent(in) :: for_master, still_sequential
integer, intent(out) :: imess(:)
real (my_real), intent(out) :: rmess(:)
integer, intent(out) :: ni,nr

!----------------------------------------------------
! Local variables:

logical(small_logical), allocatable :: sendit(:), sendite(:), senditf(:), &
                                       seen_face(:), seen_edge(:)
character(len=HOSTLEN) :: host
integer :: i,k,ind,lev,elem,edge,rind,vert,elem_per_proc,edge_per_proc, &
           vert_per_proc,nproc,my_processor,melem,medge,mvert,astat,ssize, &
           neigh_lid(NEIGHBORS_PER_ELEMENT),mface,face_per_proc,face,ivert
type(hash_key) :: boundary_gid
logical :: visited_vert(grid%biggest_vert)

!----------------------------------------------------
! Begin executable code

if (.not. grid%errind_up2date) then
   call all_error_indicators(grid,ref_control%error_estimator)
endif

boundary_gid = BOUNDARY
my_processor = my_proc(procs)
nproc = num_proc(procs)

if (still_sequential .and. for_master .and. my_processor/=1) return

! if packing for the master, determine which elements and edges to send to the
! master as all that I own or have a descendent that I own

if (for_master) then
   allocate(sendit(grid%biggest_elem),sendite(grid%biggest_edge),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in pack grid",procs=procs)
      return
   endif
   if (global_element_kind == TETRAHEDRAL_ELEMENT) then
      allocate(senditf(grid%biggest_face),stat=astat)
   else
      allocate(senditf(0),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in pack grid",procs=procs)
         return
      endif
   endif
   call set_sendit(grid,sendit,sendite,senditf,my_processor)
endif

! pack the command number

imess(1) = GRAPHICS_GRID

! determine the dimension for elements, faces, edges and vertices needed on
! the graphics server

melem = 0
mface = 0
medge = 0
mvert = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      melem = max(melem,elem)
      do i=1,VERTICES_PER_ELEMENT
         mvert = max(mvert,grid%element(elem)%vertex(i))
      end do
      do i=1,EDGES_PER_ELEMENT
         medge = max(medge,grid%element(elem)%edge(i))
      end do
      do i=1,FACES_PER_ELEMENT
         mface = max(mface,grid%element(elem)%face(i))
      end do
      elem = grid%element(elem)%next
   end do
end do

if (for_master .and. nproc > 1) then
   if (still_sequential) then
      elem_per_proc=10**int(log10(real(melem,my_real))+1)
      edge_per_proc=10**int(log10(real(medge,my_real))+1)
      vert_per_proc=10**int(log10(real(mvert,my_real))+1)
      select case (global_element_kind)
      case (TRIANGULAR_ELEMENT)
         face_per_proc = 0
      case (TETRAHEDRAL_ELEMENT)
         face_per_proc=10**int(log10(real(mface,my_real))+1)
      end select
   else
      elem_per_proc = &
              10**int(log10(phaml_global_max(procs,real(melem,my_real),730))+1)
      edge_per_proc = &
              10**int(log10(phaml_global_max(procs,real(medge,my_real),735))+1)
      vert_per_proc = &
              10**int(log10(phaml_global_max(procs,real(mvert,my_real),740))+1)
      select case (global_element_kind)
      case (TRIANGULAR_ELEMENT)
         face_per_proc = 0
      case (TETRAHEDRAL_ELEMENT)
         face_per_proc = &
              10**int(log10(phaml_global_max(procs,real(mface,my_real),736))+1)
      end select
   endif
   imess(2) = (nproc+1)*elem_per_proc
   imess(3) = (nproc+1)*edge_per_proc
   imess(4) = (nproc+1)*vert_per_proc
else
   imess(2) = melem
   imess(3) = medge
   imess(4) = mvert
endif
imess(5) = size(grid%vertex_solution,2)
imess(6) = size(grid%vertex_solution,3)

! pack the processor info

imess(7) = nproc
imess(8) = my_processor

select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   ind = 8
case (TETRAHEDRAL_ELEMENT)
   if (for_master .and. nproc > 1) then
      imess(9) = (nproc+1)*face_per_proc
   else
      imess(9) = mface
   endif
   ind = 9
end select

host = hostname(procs)
do i=1,HOSTLEN
   imess(ind+ i) = ichar(host(i:i))
end do
ind = ind+HOSTLEN

! pack grid_type variables

imess(ind+1) = grid%nsoln/grid%system_size ! number of solutions or eigenvectors
imess(ind+2) = grid%system_size
ind = ind+2

! pack number of children of the root

ind = scalar_imess_size(grid)
rind = 0

! pack the elements

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (for_master) then
         if (.not. sendit(elem)) then
            elem = grid%element(elem)%next
            cycle
         endif
         if (nproc > 1) then
            imess(ind+1) = elem + my_processor*elem_per_proc
         else
            imess(ind+1) = elem
         endif
      else
         imess(ind+1) = elem
      endif
      ind = ind + 1
      if (associated(grid%element(elem)%solution)) then
         ssize = size(grid%element(elem)%solution)
      else
         ssize = 0
      endif
      imess(ind+1) = ssize
      ind = ind + 1
      call hash_pack_key(grid%element(elem)%gid,imess,ind+1)
      ind = ind + KEY_SIZE
      if (global_element_kind == TRIANGULAR_ELEMENT) then
         rmess(rind+1) = maxval(grid%element_errind(elem,:))
         rmess(rind+2) = grid%element(elem)%work
         rind = rind + 2
      endif
      if (ssize /= 0) then
         rmess(rind+1:rind+ssize) = reshape(grid%element(elem)%solution,(/ssize/))
         rind = rind + ssize
         if (grid%have_true) then
            rmess(rind+1:rind+ssize) = reshape(grid%element(elem)%exact,(/ssize/))
         else
            rmess(rind+1:rind+ssize) = 0.0_my_real
         endif
         where (rmess(rind+1:rind+ssize) == huge(0.0_my_real))
            rmess(rind+1:rind+ssize) = 0.0_my_real
         endwhere
         rind = rind + ssize
      endif
      do i=1,VERTICES_PER_ELEMENT
         call hash_pack_key(grid%vertex(grid%element(elem)%vertex(i))%gid, &
                            imess,ind+1+(i-1)*KEY_SIZE)
      end do
      ind = ind + VERTICES_PER_ELEMENT*KEY_SIZE
      do i=1,EDGES_PER_ELEMENT
         call hash_pack_key(grid%edge(grid%element(elem)%edge(i))%gid,imess, &
                           ind+1+(i-1)*KEY_SIZE)
      end do
      ind = ind + EDGES_PER_ELEMENT*KEY_SIZE
      do i=1,FACES_PER_ELEMENT
         call hash_pack_key(grid%face(grid%element(elem)%face(i))%gid,imess, &
                           ind+1+(i-1)*KEY_SIZE)
      end do
      ind = ind + FACES_PER_ELEMENT*KEY_SIZE
      imess(ind+1) = grid%element(elem)%degree
      ind = ind + 1
      imess(ind+1) = grid%element(elem)%level
      ind = ind + 1
      if (global_element_kind == TRIANGULAR_ELEMENT) then
         call hash_pack_key(grid%vertex(grid%element(elem)%in)%gid,imess, &
                            ind+1)
         ind = ind + KEY_SIZE
         call hash_pack_key(grid%vertex(grid%element(elem)%out)%gid,imess, &
                            ind+1)
         ind = ind + KEY_SIZE
         imess(ind+1:ind+MAX_CHILD) = grid%element(elem)%order
         ind = ind + MAX_CHILD
      endif
      if (grid%element(elem)%isleaf) then
         imess(ind+1) = 1
      else
         imess(ind+1) = 0
      endif
      ind = ind + 1
      if (grid%element(elem)%iown) then
         imess(ind+1) = my_proc(procs)
      else
         imess(ind+1) = 0
      endif
      ind = ind + 1
      if (global_element_kind == TRIANGULAR_ELEMENT) then
         neigh_lid = get_neighbors(elem,grid)
         do i=1,NEIGHBORS_PER_ELEMENT
            if (neigh_lid(i) == BOUNDARY) then
               call hash_pack_key(boundary_gid,imess,ind+1+(i-1)*KEY_SIZE)
            else
               call hash_pack_key(grid%element(neigh_lid(i))%gid,imess, &
                                  ind+1+(i-1)*KEY_SIZE)
            endif
         end do
         ind = ind + NEIGHBORS_PER_ELEMENT*KEY_SIZE
      endif
      elem = grid%element(elem)%next
   end do
end do

imess(ind+1) = END_OF_ELEMENTS
ind = ind+1

! pack the vertices

visited_vert = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do ivert=1,VERTICES_PER_ELEMENT
         vert = grid%element(elem)%vertex(ivert)
         if (visited_vert(vert)) cycle
         visited_vert(vert) = .true.
         if (for_master) then
            if (.not. grid%element(grid%vertex(vert)%assoc_elem)%iown) cycle
            if (nproc > 1) then
               imess(ind+1) = vert + vert_per_proc*my_processor
               imess(ind+2) = grid%vertex(vert)%assoc_elem + &
                              my_processor*elem_per_proc
            else
               imess(ind+1) = vert
               imess(ind+2) = grid%vertex(vert)%assoc_elem
            endif
         else
            imess(ind+1) = vert
            imess(ind+2) = grid%vertex(vert)%assoc_elem
         endif
         ind = ind + 2
         call hash_pack_key(grid%vertex(vert)%gid,imess,ind+1)
         ind = ind + KEY_SIZE
         rmess(rind+1) = grid%vertex(vert)%coord%x
         rmess(rind+2) = grid%vertex(vert)%coord%y
         rind = rind + 2
         if (global_element_kind == TETRAHEDRAL_ELEMENT) then
            rmess(rind+1) = zcoord(grid%vertex(vert)%coord)
            rind = rind + 1
         endif
         rmess(rind+1:rind+grid%nsoln) = reshape(grid%vertex_solution(vert,:,:), &
                                                 (/grid%nsoln/))
         rind = rind + grid%nsoln
         if (associated(grid%vertex_exact)) then
            where (reshape(grid%vertex_exact(vert,:,:),(/grid%nsoln/)) == huge(0.0_my_real))
               rmess(rind+1:rind+grid%nsoln) = 0.0_my_real
            elsewhere
               rmess(rind+1:rind+grid%nsoln) = reshape(grid%vertex_exact(vert,:,:), &
                                                       (/grid%nsoln/))
            endwhere
         else
            rmess(rind+1:rind+grid%nsoln) = 0.0_my_real
         endif
         rind = rind + grid%nsoln
      end do
      elem = grid%element(elem)%next
   end do
end do
imess(ind+1) = END_OF_VERTICES
ind = ind+1

! pack the faces

if (global_element_kind == TETRAHEDRAL_ELEMENT) then
   allocate(seen_face(mface),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in pack grid",procs=procs)
      return
   endif
   seen_face = .false.

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         do k=1,FACES_PER_ELEMENT
            face = grid%element(elem)%face(k)
            if (seen_face(face)) cycle
            seen_face(face) = .true.
            if (associated(grid%face(face)%solution)) then
               ssize = size(grid%face(face)%solution)
            else
               ssize = 0
            endif
            if (for_master) then
               if (.not. senditf(face)) cycle
               if (nproc > 1) then
                  imess(ind+1) = face + face_per_proc*my_processor
                  imess(ind+2) = ssize
                  imess(ind+3) = grid%face(face)%assoc_elem + &
                                 my_processor*elem_per_proc
               else
                  imess(ind+1) = face
                  imess(ind+2) = ssize
                  imess(ind+3) = grid%face(face)%assoc_elem
               endif
            else
               imess(ind+1) = face
               imess(ind+2) = ssize
               imess(ind+3) = grid%face(face)%assoc_elem
            endif
            ind = ind + 3
            call hash_pack_key(grid%face(face)%gid,imess,ind+1)
            ind = ind + KEY_SIZE
            call hash_pack_key(grid%vertex(grid%face(face)%vertex(1))%gid, &
                               imess,ind+1)
            ind = ind + KEY_SIZE
            call hash_pack_key(grid%vertex(grid%face(face)%vertex(2))%gid, &
                               imess,ind+1)
            ind = ind + KEY_SIZE
            call hash_pack_key(grid%vertex(grid%face(face)%vertex(3))%gid, &
                               imess,ind+1)
            ind = ind + KEY_SIZE
            imess(ind+1) = grid%face(face)%degree
            ind = ind + 1
            if (ssize /= 0) then
               rmess(rind+1:rind+ssize) = reshape(grid%face(face)%solution,(/ssize/))
               rind = rind + ssize
               if (grid%have_true) then
                  rmess(rind+1:rind+ssize) = reshape(grid%face(face)%exact,(/ssize/))
               else
                  rmess(rind+1:rind+ssize) = 0.0_my_real
               endif
               where (rmess(rind+1:rind+ssize) == huge(0.0_my_real))
                  rmess(rind+1:rind+ssize) = 0.0_my_real
               endwhere
               rind = rind + ssize
            endif
         end do
         elem = grid%element(elem)%next
      end do
   end do

   deallocate(seen_face,stat=astat)
   imess(ind+1) = END_OF_FACES
   ind = ind+1

endif

! pack the edges

allocate(seen_edge(medge),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in pack grid",procs=procs)
   return
endif
seen_edge = .false.

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do k=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(k)
         if (seen_edge(edge)) cycle
         seen_edge(edge) = .true.
         if (associated(grid%edge(edge)%solution)) then
            ssize = size(grid%edge(edge)%solution)
         else
            ssize = 0
         endif
         if (for_master) then
            if (.not. sendite(edge)) cycle
            if (nproc > 1) then
               imess(ind+1) = edge + edge_per_proc*my_processor
               imess(ind+2) = ssize
               imess(ind+3) = grid%edge(edge)%assoc_elem + &
                              my_processor*elem_per_proc
            else
               imess(ind+1) = edge
               imess(ind+2) = ssize
               imess(ind+3) = grid%edge(edge)%assoc_elem
            endif
         else
            imess(ind+1) = edge
            imess(ind+2) = ssize
            imess(ind+3) = grid%edge(edge)%assoc_elem
         endif
         ind = ind + 3
         call hash_pack_key(grid%edge(edge)%gid,imess,ind+1)
         ind = ind + KEY_SIZE
         call hash_pack_key(grid%vertex(grid%edge(edge)%vertex(1))%gid, &
                            imess,ind+1)
         ind = ind + KEY_SIZE
         call hash_pack_key(grid%vertex(grid%edge(edge)%vertex(2))%gid, &
                            imess,ind+1)
         ind = ind + KEY_SIZE
         imess(ind+1) = grid%edge(edge)%degree
         ind = ind + 1
         if (ssize /= 0) then
            rmess(rind+1:rind+ssize) = reshape(grid%edge(edge)%solution,(/ssize/))
            rind = rind + ssize
            if (grid%have_true) then
               rmess(rind+1:rind+ssize) = reshape(grid%edge(edge)%exact,(/ssize/))
            else
               rmess(rind+1:rind+ssize) = 0.0_my_real
            endif
            where (rmess(rind+1:rind+ssize) == huge(0.0_my_real))
               rmess(rind+1:rind+ssize) = 0.0_my_real
            endwhere
            rind = rind + ssize
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

deallocate(seen_edge,stat=astat)
imess(ind+1) = END_OF_EDGES

ni = ind+1
nr = rind

if (ni > size(imess) .or. nr > size(rmess)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Incorrect allocation for message to graphics process.", &
              intlist=(/ni,size(imess),nr,size(rmess)/))
   stop
endif

if (for_master) deallocate(sendit,sendite,senditf,stat=astat)

end subroutine pack_grid

!          -----------------
! function scalar_imess_size and friends
!          -----------------
! These functions give the size of the message containing the grid.  They
! give the amount of integer and real space used for scalar information and
! for each element, edge and vertex.

integer function scalar_imess_size(grid)
type(grid_type), intent(in) :: grid
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   scalar_imess_size = 10 + HOSTLEN
case (TETRAHEDRAL_ELEMENT)
   scalar_imess_size = 11 + HOSTLEN
end select
end function scalar_imess_size

integer function elem_imess_size(grid)
type(grid_type), intent(in) :: grid
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   elem_imess_size = 6 + MAX_CHILD + &
                    (3 + EDGES_PER_ELEMENT + VERTICES_PER_ELEMENT + &
                     NEIGHBORS_PER_ELEMENT)*KEY_SIZE
case (TETRAHEDRAL_ELEMENT)
   elem_imess_size = 6 + &
                    (1 + EDGES_PER_ELEMENT + VERTICES_PER_ELEMENT + &
                     FACES_PER_ELEMENT)*KEY_SIZE
end select
end function elem_imess_size

integer function face_imess_size(grid)
type(grid_type), intent(in) :: grid
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   face_imess_size = 0
case (TETRAHEDRAL_ELEMENT)
   face_imess_size = 4 + 4*KEY_SIZE
end select
end function face_imess_size

integer function edge_imess_size(grid)
type(grid_type), intent(in) :: grid
edge_imess_size = 4 + 3*KEY_SIZE
end function edge_imess_size

integer function vert_imess_size(grid)
type(grid_type), intent(in) :: grid
vert_imess_size = 2 + KEY_SIZE
end function vert_imess_size

integer function scalar_rmess_size(grid)
type(grid_type), intent(in) :: grid
logical(small_logical) :: seen_face(grid%biggest_face), &
                          seen_edge(grid%biggest_edge)
integer :: lev, elem, face, edge, i
seen_face = .false.
seen_edge = .false.
scalar_rmess_size = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (associated(grid%element(elem)%solution)) then
         scalar_rmess_size = scalar_rmess_size + &
                             2*size(grid%element(elem)%solution)
      endif
      do i=1,FACES_PER_ELEMENT
         face = grid%element(elem)%face(i)
         if (seen_face(face)) cycle
         seen_face(face) = .true.
         if (associated(grid%face(face)%solution)) then
            scalar_rmess_size = scalar_rmess_size + &
                                2*size(grid%face(face)%solution)
         endif
      end do
      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)
         if (seen_edge(edge)) cycle
         seen_edge(edge) = .true.
         if (associated(grid%edge(edge)%solution)) then
            scalar_rmess_size = scalar_rmess_size + &
                                2*size(grid%edge(edge)%solution)
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do
end function scalar_rmess_size

integer function elem_rmess_size(grid)
type(grid_type), intent(in) :: grid
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   elem_rmess_size = 2
case (TETRAHEDRAL_ELEMENT)
   elem_rmess_size = 0
end select
end function elem_rmess_size

integer function face_rmess_size(grid)
type(grid_type), intent(in) :: grid
face_rmess_size = 0
end function face_rmess_size

integer function edge_rmess_size(grid)
type(grid_type), intent(in) :: grid
edge_rmess_size = 0
end function edge_rmess_size

integer function vert_rmess_size(grid)
type(grid_type), intent(in) :: grid
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   vert_rmess_size = 2 + 2*grid%nsoln
case (TETRAHEDRAL_ELEMENT)
   vert_rmess_size = 3 + 2*grid%nsoln
end select
end function vert_rmess_size

!          -----------
subroutine merge_grids(grid,procs,imess,ni,rmess,nr,still_sequential)
!          -----------

!----------------------------------------------------
! This routine is called by the master to merge the grids from each slave
! into a single grid for the master's graphics processor
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
integer, intent(out) :: ni,nr
logical, intent(in) :: still_sequential
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

character(len=HOSTLEN) :: host
integer :: i, j, iind, rind, piind, prind, np, isize, rsize, astat, ssize
integer, allocatable :: ord(:)
type recv_type
   integer, pointer :: imess(:)
   real(my_real), pointer :: rmess(:)
   integer :: ni,nr,proc,iind,rind
end type recv_type
type (recv_type), allocatable :: recv(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! receive all the slave messages

if (still_sequential) then
   np = 1
else
   np = num_proc(procs)
endif

allocate(recv(np),ord(np),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in merge grids",procs=procs)
   return
endif
do i=1,np
   call phaml_recv(procs,recv(i)%proc,recv(i)%imess,recv(i)%ni,recv(i)%rmess, &
                   recv(i)%nr,720)
   ord(recv(i)%proc) = i
end do

allocate(imess(sum(recv%ni)),rmess(sum(recv%nr)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in merge grids",procs=procs)
   return
endif

! pack the scalar information from processor 1; the others should be the same

imess(1:scalar_imess_size(grid)) = recv(1)%imess(1:scalar_imess_size(grid))

! change the processor to this processor

imess(8) = my_proc(procs)
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   iind = 8
case (TETRAHEDRAL_ELEMENT)
   iind = 9
end select
host = hostname(procs)
do i=1,HOSTLEN
   imess(iind+ i) = ichar(host(i:i))
end do
iind = scalar_imess_size(grid)
rind = 0

! pack elements from each processor

isize = elem_imess_size(grid)
rsize = elem_rmess_size(grid)

do j=1,np
   i = ord(j)
   piind = scalar_imess_size(grid)
   prind = 0
   do while (recv(i)%imess(piind+1) /= END_OF_ELEMENTS)
      imess(iind+1:iind+isize) = recv(i)%imess(piind+1:piind+isize)
      ssize = imess(iind+2)
      iind = iind   + isize
      piind = piind + isize
      rmess(rind+1:rind+rsize+2*ssize) = recv(i)%rmess(prind+1:prind+rsize+2*ssize)
      rind = rind   + rsize + 2*ssize
      prind = prind + rsize + 2*ssize
   end do
   recv(i)%iind = piind + 1 ! bookmark for when I copy the vertices
   recv(i)%rind = prind
end do
imess(iind+1) = END_OF_ELEMENTS
iind = iind + 1

! pack vertices from each processor

isize = vert_imess_size(grid)
rsize = vert_rmess_size(grid)

do j=1,np
   i = ord(j)
   piind = recv(i)%iind
   prind = recv(i)%rind
   do while (recv(i)%imess(piind+1) /= END_OF_VERTICES)
      imess(iind+1:iind+isize) = recv(i)%imess(piind+1:piind+isize)
      iind = iind   + isize
      piind = piind + isize
      rmess(rind+1:rind+rsize) = recv(i)%rmess(prind+1:prind+rsize)
      rind = rind   + rsize
      prind = prind + rsize
   end do
   recv(i)%iind = piind + 1 ! bookmark for when I copy the faces or edges
   recv(i)%rind = prind
end do
imess(iind+1) = END_OF_VERTICES
iind = iind + 1

! pack faces from each processor

if (global_element_kind == TETRAHEDRAL_ELEMENT) then

   isize = face_imess_size(grid)
   rsize = face_rmess_size(grid)

   do j=1,np
      i = ord(j)
      piind = recv(i)%iind
      prind = recv(i)%rind
      do while (recv(i)%imess(piind+1) /= END_OF_FACES)
         imess(iind+1:iind+isize) = recv(i)%imess(piind+1:piind+isize)
         ssize = imess(iind+2)
         iind = iind   + isize
         piind = piind + isize
         rmess(rind+1:rind+rsize+2*ssize) = recv(i)%rmess(prind+1:prind+rsize+2*ssize)
         rind = rind   + rsize + 2*ssize
         prind = prind + rsize + 2*ssize
      end do
      recv(i)%iind = piind + 1 ! bookmark for when I copy the edges
      recv(i)%rind = prind
   end do
   imess(iind+1) = END_OF_FACES
   iind = iind + 1

endif

! pack edges from each processor

isize = edge_imess_size(grid)
rsize = edge_rmess_size(grid)

do j=1,np
   i = ord(j)
   piind = recv(i)%iind
   prind = recv(i)%rind
   do while (recv(i)%imess(piind+1) /= END_OF_EDGES)
      imess(iind+1:iind+isize) = recv(i)%imess(piind+1:piind+isize)
      ssize = imess(iind+2)
      iind = iind   + isize
      piind = piind + isize
      rmess(rind+1:rind+rsize+2*ssize) = recv(i)%rmess(prind+1:prind+rsize+2*ssize)
      rind = rind   + rsize + 2*ssize
      prind = prind + rsize + 2*ssize
   end do
end do
imess(iind+1) = END_OF_EDGES

ni = iind + 1
nr = rind

do i=1,np
   if (recv(i)%ni /= 0) deallocate(recv(i)%imess,stat=astat)
   if (recv(i)%nr /= 0) deallocate(recv(i)%rmess,stat=astat)
end do
deallocate(recv,ord,stat=astat)

return
end subroutine merge_grids

!          ----------
subroutine set_sendit(grid,sendit,sendite,senditf,my_processor)
!          ----------

!----------------------------------------------------
! This routine sets flags to indicate whether or not to send an element
! to the master for the composite grid graphics.  Send it if it or any
! of its descendents are owned by my_processor
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical(small_logical), intent(inout) :: sendit(:), sendite(:), senditf(:)
integer, intent(in) :: my_processor
!----------------------------------------------------
! Local variables:

integer :: elem
!----------------------------------------------------
! Begin executable code

! TEMP3D as long as it is not parallel we send everything

if (global_element_kind == TETRAHEDRAL_ELEMENT) then
   if (my_processor == 0 .or. my_processor == 1) then
      sendit = .true.; sendite = .true.; senditf = .true.
      return
   else
      call fatal("need to finish writing set_sendit_recur for 3D")
      stop
   endif
endif

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call set_sendit_recur(grid,sendit,sendite,senditf,my_processor,elem)
   elem = grid%element(elem)%next
end do
end subroutine set_sendit

!                    ----------------
recursive subroutine set_sendit_recur(grid,sendit,sendite,senditf, &
                                      my_processor,subroot)
!                    ----------------

!----------------------------------------------------
! This routine recursively traverses the tree to do the work of set_sendit
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical(small_logical), intent(inout) :: sendit(:), sendite(:), senditf(:)
integer, intent(in) :: my_processor, subroot
!----------------------------------------------------
! Local variables:

integer :: child(MAX_CHILD), i, allc(MAX_CHILD), face, edge, cedge
!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
child = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)

if (child(1) == NO_CHILD) then

! for a leaf, return whether or not I own it

   sendit(subroot) = grid%element(subroot)%iown
   do i=1,FACES_PER_ELEMENT
      face = grid%element(subroot)%face(i)
      senditf(face) = grid%element(grid%face(face)%assoc_elem)%iown
   end do
   do i=1,EDGES_PER_ELEMENT
      edge = grid%element(subroot)%edge(i)
      sendite(edge) = grid%element(grid%edge(edge)%assoc_elem)%iown
   end do

else

! otherwise, see if any of the children are owned
! RESTRICTION bisected triangles; only the third edge needs children checked
! TEMP 3D for tetrahedra, I need to know which edge is the bisection edge and
!         which two faces are bisected and which edges and faces are the
!         children of them
   if (global_element_kind == TETRAHEDRAL_ELEMENT) then
      call fatal("need to finish writing set_sendit_recur for 3D")
      stop
   endif

   sendit(subroot) = .false.
   edge = grid%element(subroot)%edge(3)
   sendite(edge) = .false.
   do i=1,MAX_CHILD
      if (child(i) /= NO_CHILD) then
         call set_sendit_recur(grid,sendit,sendite,senditf,my_processor, &
                               child(i))
         sendit(subroot) = sendit(subroot) .or. sendit(child(i))
         cedge = grid%element(child(i))%edge(2)
         if (sendite(cedge)) sendite(edge) = .true.
      endif
   end do

endif

end subroutine set_sendit_recur

!          -----------------------
subroutine child_order_from_zoltan(grid,procs,lb,partition_method, &
                                   still_sequential)
!          -----------------------

!----------------------------------------------------
! This routine gets the child order (for the reftree sfc) from Zoltan.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: partition_method
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: i, nelem_init, maxelem_init, child1, child2, parent, &
           astat, jerr, allc(MAX_CHILD), children(MAX_CHILD)
integer, allocatable :: order(:)
type(hash_key) :: tempkey
!----------------------------------------------------
! Begin executable code

! must be using Zoltan REFTREE as the partition method

if (partition_method /= ZOLTAN_REFTREE) return

! count number of elements in the initial grid

nelem_init = 0
maxelem_init = 0
i = grid%head_level_elem(1)
do while (i /= END_OF_LIST)
   nelem_init = nelem_init + 1
   maxelem_init = max(maxelem_init,i)
   i = grid%element(i)%next
end do

! space for order returned by zoltan and retaining element linked list

allocate(order(KEY_SIZE*(1 + 3*nelem_init + 7*(grid%nelem-grid%nelem_leaf))), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in child_order_from_zoltan",procs=procs)
   return
endif

! get the order from zoltan

order = -1
call zoltan_child_order(order,jerr,lb,grid,procs,still_sequential)
if (jerr /= 0) then
   call warning("zoltan_child_order returned error.", &
                "Retaining PHAML's child order.",(/jerr/))
   deallocate(order,stat=astat)
   return
endif

! the first 1+nelem_init are the root and it's children (initial grid). Skip
! the root and reorder the initial grid by resetting previous/next

i = KEY_SIZE+1
tempkey = hash_unpack_key(order,i)
child1 = hash_decode_key(tempkey, grid%elem_hash)
if (child1 == HASH_NOT_FOUND) then
   call fatal("zoltan child order returned unknown GID",procs=procs)
   stop
endif
grid%head_level_elem(1) = child1
grid%element(child1)%previous = END_OF_LIST
tempkey = hash_unpack_key(order,i+KEY_SIZE)
grid%element(child1)%in = hash_decode_key(tempkey, grid%vert_hash)
if (grid%element(child1)%in == HASH_NOT_FOUND) then
   call fatal("zoltan child order returned unknown GID for in",procs=procs)
   stop
endif
tempkey = hash_unpack_key(order,i+2*KEY_SIZE)
grid%element(child1)%out = hash_decode_key(tempkey, grid%vert_hash)
if (grid%element(child1)%out == HASH_NOT_FOUND) then
   call fatal("zoltan child order returned unknown GID for out",procs=procs)
   stop
endif
child2 = child1
do i=4*KEY_SIZE+1,(3*nelem_init+1)*KEY_SIZE,3*KEY_SIZE
   tempkey = hash_unpack_key(order,i)
   child1 = hash_decode_key(tempkey, grid%elem_hash)
   if (child1 == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID",procs=procs)
      stop
   endif
   grid%element(child2)%next = child1
   grid%element(child1)%previous = child2
   tempkey = hash_unpack_key(order,i+KEY_SIZE)
   grid%element(child1)%in = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%in == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for in",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+2*KEY_SIZE)
   grid%element(child1)%out = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%out == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for out",procs=procs)
      stop
   endif
   child2 = child1
end do
grid%element(child2)%next = END_OF_LIST
i = (3*nelem_init+1)*KEY_SIZE + 1

! don't get refinements if still sequential and not processor 1
if (still_sequential .and. my_proc(procs) /= 1) i = size(order)+1

! go through the rest of the list returned by Zoltan and reset order
! to correspond to the order of the children

do while (i < size(order))
! RESTRICTION bisection
   tempkey = hash_unpack_key(order,i)
   if (tempkey == -1) then
      i = i+7*KEY_SIZE
      cycle
   endif
   parent = hash_decode_key(tempkey, grid%elem_hash)
   if (parent == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for parent",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+KEY_SIZE)
   child1 = hash_decode_key(tempkey, grid%elem_hash)
   if (child1 == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for child1",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+2*KEY_SIZE)
   grid%element(child1)%in = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%in == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for in",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+3*KEY_SIZE)
   grid%element(child1)%out = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%out == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for out",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+4*KEY_SIZE)
   child2 = hash_decode_key(tempkey, grid%elem_hash)
   if (child2 == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for child2",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+5*KEY_SIZE)
   grid%element(child2)%in = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child2)%in == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for in",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+6*KEY_SIZE)
   grid%element(child2)%out = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child2)%out == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for out",procs=procs)
      stop
   endif
   allc = ALL_CHILDREN
   children = get_child_lid(grid%element(parent)%gid,allc,grid%elem_hash)
   if (child1 == children(1) .and. child2 == children(2)) then
      grid%element(parent)%order = (/1,2/)
   elseif (child2 == children(1) .and. child1 == children(2)) then
      grid%element(parent)%order = (/2,1/)
   else
      call warning("in draw_grid zoltan did not return the right children", &
                   "parent, children, zoltan", &
                   intlist=(/parent,children,child1,child2/))
   endif
   i = i+7*KEY_SIZE
end do

deallocate(order,stat=astat)

end subroutine child_order_from_zoltan

!          ----------
subroutine store_grid(grid,procs,still_sequential,unit,fmt,comp,eigen)
!          ----------

!----------------------------------------------------
! This routine calls the appropriate routine to save the grid to a file
! in the requested format
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: unit, fmt, comp, eigen
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

select case(fmt)
case (GRIDFILE_POLY, GRIDFILE_POLY_SOLN)
   if (global_element_kind == TRIANGULAR_ELEMENT) then
      call store_grid_poly(grid,procs,still_sequential,unit,fmt,comp,eigen)
   else
      call warning("can only save a .poly file with triangular elements; grid not saved")
      return
   endif
case (GRIDFILE_MSH, GRIDFILE_MSH_SOLN)
   call store_grid_msh(grid,procs,still_sequential,unit,fmt,comp,eigen)
case default
   ierr = USER_INPUT_ERROR
   call fatal("unknown format requested in store_grid")
   stop
end select

end subroutine store_grid

!          ---------------
subroutine store_grid_poly(grid,procs,still_sequential,unit,fmt,comp,eigen)
!          ---------------

!----------------------------------------------------
! This routine stores a triangle grid in triangle's .poly format.
! If fmt is GRIDFILE_POLY_SOLN then one or more solutions are saved as
! attributes of the vertices.  If comp is -1, all components are saved,
! otherwise comp tells which component to save.  Same with eigen.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: unit, fmt, comp, eigen
!----------------------------------------------------
! Local variables:

integer :: myproc, ni, nr, astat, i, lev, vert, errcode, proc, nsolut, &
           ipoint, rpoint, elem, edge, nhole, hole, k, piece, ivert
integer, allocatable :: vertlist(:), iperm(:), end_hole(:)
logical :: iprint
logical, allocatable :: seen_edge(:)
character(len=28) :: fmat
integer, allocatable, target :: isend(:)
integer, pointer :: irecv(:)
real(my_real) :: x1, y1, x2, y2
real(my_real), allocatable, target :: rsend(:)
real(my_real), pointer :: rrecv(:)
logical :: visited_vert(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

myproc = my_proc(procs)

! TEMP requires one processor, but it is OK if it is still sequential

if (num_proc(procs) > 1 .and. .not. still_sequential) then
   print *,"only one processor if saving the grid as .poly"
   stop
endif
if (myproc > 1) return

! Slave 1 has the grid.  For each section, pack it into a message.  If the
! compilation is sequential, then slave 1 sets recv to point to send and
! prints.  If not, then slave 1 sends send to recv on the master and the
! master prints.

iprint = (parallel == SEQUENTIAL .and. myproc == 1) .or. &
         (parallel /= SEQUENTIAL .and. myproc == MASTER)
nullify(irecv,rrecv)

! The description of triangle .poly and .node files says the vertices must
! be in increasing order with no gaps.  Experimenting confirmed this must
! be the case for showme and triangle to work correctly.  Make sure this
! is met by 1) make a list of the vertices, 2) sort the vertices, and
! 3) make sure the last vertex in the list is the length of the list

if (myproc == 1) then

   allocate(vertlist(grid%nvert),iperm(grid%nvert),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_grid_poly",procs=procs)
      stop
   endif

   i = 0
   visited_vert = .false.
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         do ivert=1,VERTICES_PER_ELEMENT
            vert = grid%element(elem)%vertex(ivert)
            if (visited_vert(vert)) cycle
            visited_vert(vert) = .true.
            i = i + 1
            vertlist(i) = vert
         end do
         elem = grid%element(elem)%next
      end do
   end do

   call sort(vertlist,grid%nvert,iperm,1,errcode)

   if (vertlist(iperm(grid%nvert)) == grid%nvert) then
      call phaml_send(procs,MASTER,(/0/),1,(/0.0_my_real/),0,760)
   else
      call phaml_send(procs,MASTER,(/1/),1,(/0.0_my_real/),0,760)
      deallocate(vertlist,iperm)
      return
   endif

   deallocate(vertlist,iperm)

else ! master

   call phaml_recv(procs,proc,irecv,ni,rrecv,nr,760)
   if (irecv(1) == 1) then
      call warning("The list of vertex indices has gaps; cannot create a .poly file")
      deallocate(irecv)
      if (associated(rrecv)) deallocate(rrecv)
      return
   else
      deallocate(irecv)
      if (associated(rrecv)) deallocate(rrecv)
   endif

endif

! number of solutions to include as attributes

if (fmt == GRIDFILE_POLY_SOLN) then
   nsolut = 1
   if (comp == -1)  nsolut = nsolut*size(grid%vertex_solution,dim=2)
   if (eigen == -1) nsolut = nsolut*size(grid%vertex_solution,dim=3)
   if (nsolut > 997) then
      call warning("Too many solutions per vertex in store_grid_msh", &
                   "Solutions not included in grid file")
      nsolut = 0
   endif
else
   nsolut = 0
endif

! write the vertices

if (myproc == 1) then
   allocate(isend(2*grid%nvert+1),rsend((2+nsolut)*grid%nvert),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_grid_poly",procs=procs)
      stop
   endif
   isend(1) = grid%nvert
   ipoint = 1
   rpoint = 0
   do vert=1,grid%nvert
      isend(ipoint+1) = vert
      isend(ipoint+2) = grid%vertex(vert)%bmark
      ipoint = ipoint + 2
      rsend(rpoint+1) = grid%vertex(vert)%coord%x
      rsend(rpoint+2) = grid%vertex(vert)%coord%y
      rpoint = rpoint + 2
      if (nsolut /= 0) then
         if (comp == -1) then
            if (eigen == -1) then
               rsend(rpoint+1:rpoint+nsolut) = &
                  reshape(grid%vertex_solution(vert,:,:),(/nsolut/))
            else
               rsend(rpoint+1:rpoint+nsolut) = &
                  grid%vertex_solution(vert,:,eigen)
            endif
         else
            if (eigen == -1) then
               rsend(rpoint+1:rpoint+nsolut) = &
                  grid%vertex_solution(vert,comp,:)
            else
               rsend(rpoint+1) = grid%vertex_solution(vert,comp,eigen)
            endif
         endif
         rpoint = rpoint + nsolut
      endif
   end do
   ni = ipoint
   nr = rpoint
   if (parallel == SEQUENTIAL) then
      irecv => isend
      rrecv => rsend
   else
      call phaml_send(procs,MASTER,isend,ni,rsend,nr,761)
   endif
else ! MASTER
   call phaml_recv(procs,proc,irecv,ni,rrecv,nr,761)
endif

if (iprint) then
   fmat = "(SS,1P,I11,    E21.13E2,I11)"
   write(fmat(12:15),"(I3)") 2+nsolut
   write(unit,"(I11,I2,I4,I2)") irecv(1),2,nsolut,1
   ipoint = 1
   rpoint = 0
   do i=1,irecv(1)
      write(unit,fmat) irecv(ipoint+1),rrecv(rpoint+1:rpoint+2+nsolut), &
                       irecv(ipoint+2)
      ipoint = ipoint + 2
      rpoint = rpoint + 2+nsolut
   end do
endif

if (parallel == SEQUENTIAL) then
   nullify(irecv,rrecv)
   deallocate(isend,rsend)
else
   if (allocated(isend)) deallocate(isend)
   if (allocated(rsend)) deallocate(rsend)
   if (associated(irecv)) deallocate(irecv)
   if (associated(rrecv)) deallocate(rrecv)
endif

! The description of triangle .poly files does not say the edges must be
! in increasing order with no gaps.  Experimenting with out-of-order edges
! and gaps in the edge numbers showed that showme gives a warning but
! displays correctly, and triangle works correctly but renumbers the edges
! to be increasing and with no gaps.  So we might just as well renumber the
! edges here by IDing the edge with a counter instead of the actual edge ID.

! write the edges by going through the elements and writing any edge not
! already written

if (myproc == 1) then

   allocate(isend(4*grid%nedge+1),seen_edge(grid%biggest_edge),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_grid_poly",procs=procs)
      stop
   endif
   seen_edge = .false.

   isend(1) = 0
   ipoint = 1
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf) then
            do i=1,EDGES_PER_ELEMENT
               edge = grid%element(elem)%edge(i)
               if (seen_edge(edge)) cycle
               seen_edge(edge) = .true.
               isend(1) = isend(1) + 1
               isend(ipoint+1) = isend(1) ! new ID for the edge
               isend(ipoint+2) = grid%edge(edge)%vertex(1)
               isend(ipoint+3) = grid%edge(edge)%vertex(2)
               isend(ipoint+4) = grid%edge(edge)%bmark
               ipoint = ipoint + 4
            end do
         endif
         elem = grid%element(elem)%next
      end do
   end do
   deallocate(seen_edge)

   ni = ipoint
   if (parallel == SEQUENTIAL) then
      irecv => isend
   else
      call phaml_send(procs,MASTER,isend,ni,(/0.0_my_real/),0,762)
   endif
else ! MASTER
   call phaml_recv(procs,proc,irecv,ni,rrecv,nr,762)
endif

if (iprint) then
   write(unit,"(I11,I2)") irecv(1),1
   ipoint = 1
   do i=1,irecv(1)
      write(unit,"(4I11)") irecv(ipoint+1:ipoint+4)
      ipoint = ipoint + 4
   end do
endif

if (parallel == SEQUENTIAL) then
   nullify(irecv)
   deallocate(isend)
else
   if (allocated(isend)) deallocate(isend)
   if (allocated(rsend)) deallocate(rsend)
   if (associated(irecv)) deallocate(irecv)
   if (associated(rrecv)) deallocate(rrecv)
endif

! write the holes

! This requires that boundary_point and boundary_npiece are available.
! For each hole compute a point inside the hole as the average of the
! piece endpoints and midpoints.  This can fail if the hole is too concave.

if (iprint) then

! count the number of holes

   nhole = 1
   do while (my_boundary_npiece(nhole) > 0)
      nhole = nhole + 1
   end do
   nhole = nhole - 1
   allocate(end_hole(-1:nhole),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in store_grid_poly")
      stop
   endif
   write(unit,"(I11)") nhole

! get the piece that ends each hole, and
! parameter limits

   end_hole(-1) = 0
   do hole=0,nhole
      end_hole(hole) = end_hole(hole-1) + my_boundary_npiece(hole)
   end do

! compute the point for each hole

   do hole=1,nhole
      x2 = 0.0_my_real
      y2 = 0.0_my_real
      k = 0
      do piece=end_hole(hole-1)+1,end_hole(hole)
         call my_boundary_point(piece,grid%bp_start(piece),x1,y1)
         x2 = x2 + x1
         y2 = y2 + y1
         call my_boundary_point(piece, &
                                (grid%bp_start(piece)+grid%bp_finish(piece))/2,&
                                x1,y1)
         x2 = x2 + x1
         y2 = y2 + y1
         k = k+2
      end do
      x2 = x2/k
      y2 = y2/k
      write(unit,"(I11,SS,1P,2E21.13E2)") hole,x2,y2
   end do

   deallocate(end_hole)

endif

end subroutine store_grid_poly

!          --------------
subroutine store_grid_msh(grid,procs,still_sequential,unit,fmt,comp,eigen)
!          --------------

!----------------------------------------------------
! This routine saves the grid in Gmsh's .msh format
! This routine stores a tetrahedral grid in Gmsh's .msh format.
! If fmt is GRIDFILE_MSH_SOLN then one or more solutions are saved as
! a NodeData section.  If comp is -1, all components are saved,
! otherwise comp tells which component to save.  Same with eigen.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: unit, fmt, comp, eigen
!----------------------------------------------------
! Local variables:

integer :: i, myproc, proc, lev, vert, edge, face, elem, ntag, astat, ni, nr, &
           ipoint, rpoint, nelement, nval, nsolut, itemp, ivert
integer, allocatable, target :: isend(:)
integer, pointer :: irecv(:)
real(my_real), allocatable, target :: rsend(:)
real(my_real), pointer :: rrecv(:)
character(len=24) :: fmat
logical :: iprint
logical, allocatable :: seen_edge(:), seen_face(:)
logical :: visited_vert(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

myproc = my_proc(procs)

! TEMP requires one processor, but it is OK if it is still sequential

if (num_proc(procs) > 1 .and. .not. still_sequential) then
   print *,"only one processor if saving the grid as .msh"
   stop
endif
if (myproc > 1) return

! Slave 1 has the grid.  For each section, pack it into a message.  If the
! compilation is sequential, then slave 1 sets recv to point to send and
! prints.  If not, then slave 1 sends send to recv on the master and the
! master prints.

iprint = (parallel == SEQUENTIAL .and. myproc == 1) .or. &
         (parallel /= SEQUENTIAL .and. myproc == MASTER)
nullify(irecv,rrecv)

! mesh format

if (iprint) then
   write(unit,"(A)") "$MeshFormat"
   write(unit,"(A)") "2.2 0 8"
   write(unit,"(A)") "$EndMeshFormat"
endif

! nodes

if (myproc == 1) then
   allocate(isend(grid%nvert+1),rsend(3*grid%nvert),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_grid_msh",procs=procs)
      stop
   endif
   isend(1) = grid%nvert
   ipoint = 1
   rpoint = 0
   visited_vert = .false.
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         do ivert=1,VERTICES_PER_ELEMENT
            vert = grid%element(elem)%vertex(ivert)
            if (visited_vert(vert)) cycle
            visited_vert(vert) = .true.
            isend(ipoint+1) = vert
            ipoint = ipoint + 1
            rsend(rpoint+1) = grid%vertex(vert)%coord%x
            rsend(rpoint+2) = grid%vertex(vert)%coord%y
            rsend(rpoint+3) = zcoord(grid%vertex(vert)%coord)
            rpoint = rpoint + 3
         end do
         elem = grid%element(elem)%next
      end do
   end do
   ni = ipoint
   nr = rpoint
   if (parallel == SEQUENTIAL) then
      irecv => isend
      rrecv => rsend
   else
      call phaml_send(procs,MASTER,isend,ni,rsend,nr,750)
   endif
else ! MASTER
   call phaml_recv(procs,proc,irecv,ni,rrecv,nr,750)
endif

if (iprint) then
   write(unit,"(A)") "$Nodes"
   write(unit,"(I11)") irecv(1) ! grid%nvert
   ipoint = 1
   rpoint = 0
   do i=2,irecv(1)+1
      write(unit,"(SS,1P,I11,3E21.13E2)") irecv(ipoint+1),rrecv(rpoint+1:rpoint+3)
      ipoint = ipoint + 1
      rpoint = rpoint + 3
   end do
   write(unit,"(A)") "$EndNodes"
endif

if (parallel == SEQUENTIAL) then
   nullify(irecv,rrecv)
   deallocate(isend,rsend)
else
   if (allocated(isend)) deallocate(isend)
   if (allocated(rsend)) deallocate(rsend)
   if (associated(irecv)) deallocate(irecv)
   if (associated(rrecv)) deallocate(rrecv)
endif

! elements

if (MAX_TAG > 99) then
   call warning("Only saving the first 99 tags in .msh file.")
   ntag = 99
elseif (MAX_TAG == 0) then
   ntag = 1 ! if no tags, use 1 tag for bmark
else
   ntag = MAX_TAG
endif

if (myproc == 1) then
   select case (global_element_kind)
   case (TRIANGULAR_ELEMENT)
      allocate(isend(grid%nvert*(4+ntag) + grid%nedge*(5+ntag) + &
                     grid%nelem_leaf*(6+ntag)), stat=astat)
   case (TETRAHEDRAL_ELEMENT)
      allocate(isend(grid%nvert*(4+ntag) + grid%nedge*(5+ntag) + &
                     grid%nface*(6+ntag) + grid%nelem_leaf*(7+ntag)), &
               stat=astat)
   end select
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_grid_msh",procs=procs)
      stop
   endif

   nelement = 0
   ipoint = 1

! first the boundary points

   visited_vert = .false.
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         do ivert=1,VERTICES_PER_ELEMENT
            vert = grid%element(elem)%vertex(ivert)
            if (visited_vert(vert)) cycle
            visited_vert(vert) = .true.
            if (grid%vertex(vert)%bmark /= 0) then
               nelement = nelement + 1
               isend(ipoint+1) = vert
               isend(ipoint+2) = 15   ! element type is 1-node point
               isend(ipoint+3) = ntag
               if (ntag /= 1) then
                  isend(ipoint+4:ipoint+3+ntag) = grid%vertex(vert)%tags(1:ntag)
               endif
               isend(ipoint+4) = grid%vertex(vert)%bmark
               ipoint = ipoint + 3+ntag
               isend(ipoint+1) = vert
               ipoint = ipoint + 1
            endif
         end do
         elem = grid%element(elem)%next
      end do
   end do

! then go through the elements adding the leaf elements and any boundary edges
! and, in 3D, boundary faces not already added

   allocate(seen_edge(grid%biggest_edge),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_grid_msh",procs=procs)
      stop
   endif
   seen_edge = .false.
   if (global_element_kind == TETRAHEDRAL_ELEMENT) then
      allocate(seen_face(grid%biggest_face),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in store_grid_msh",procs=procs)
         stop
      endif
      seen_face = .false.
   endif

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf) then
! element
            nelement = nelement + 1
            isend(ipoint+1) = elem
            select case (global_element_kind)
            case (TRIANGULAR_ELEMENT)
               isend(ipoint+2) = 2    ! 3-node triangle
            case (TETRAHEDRAL_ELEMENT)
               isend(ipoint+2) = 4    ! 4-node tetrahedron
            end select
            if (MAX_TAG == 0) then
               isend(ipoint+3) = 0
               ipoint = ipoint + 3
            else
               isend(ipoint+3) = ntag
               isend(ipoint+4:ipoint+3+ntag) = grid%element(elem)%tags(1:ntag)
               ipoint = ipoint + 3+ntag
            endif
            select case (global_element_kind)
            case (TRIANGULAR_ELEMENT)
               isend(ipoint+1:ipoint+3) = grid%element(elem)%vertex
! orient triangles clockwise
               if ((grid%vertex(isend(ipoint+2))%coord%x - &
                    grid%vertex(isend(ipoint+1))%coord%x) * &
                   (grid%vertex(isend(ipoint+3))%coord%y - &
                    grid%vertex(isend(ipoint+1))%coord%y) - &
                   (grid%vertex(isend(ipoint+3))%coord%x - &
                    grid%vertex(isend(ipoint+1))%coord%x) * &
                   (grid%vertex(isend(ipoint+2))%coord%y - &
                    grid%vertex(isend(ipoint+1))%coord%y) > 0.0_my_real) then
                  itemp = isend(ipoint+1)
                  isend(ipoint+1) = isend(ipoint+2)
                  isend(ipoint+2) = itemp
               endif
               ipoint = ipoint + 3
            case (TETRAHEDRAL_ELEMENT)
               isend(ipoint+1:ipoint+4) = grid%element(elem)%vertex
               ipoint = ipoint + 4
            end select
! edges
            do i=1,EDGES_PER_ELEMENT
               edge = grid%element(elem)%edge(i)
               if (seen_edge(edge)) cycle
               seen_edge(edge) = .true.
               if (grid%edge(edge)%bmark == 0) cycle
               nelement = nelement + 1
               isend(ipoint+1) = edge
               isend(ipoint+2) = 1   ! 2-node line
               isend(ipoint+3) = ntag
               if (ntag /= 1) then
                  isend(ipoint+4:ipoint+3+ntag) = grid%edge(edge)%tags(1:ntag)
               endif
               isend(ipoint+4) = grid%edge(edge)%bmark
               ipoint = ipoint + 3+ntag
               isend(ipoint+1:ipoint+2) = grid%edge(edge)%vertex
               ipoint = ipoint + 2
            end do
! faces
            do i=1,FACES_PER_ELEMENT
               face = grid%element(elem)%face(i)
               if (seen_face(face)) cycle
               seen_face(face) = .true.
               if (grid%face(face)%bmark == 0) cycle
               nelement = nelement + 1
               isend(ipoint+1) = face
               isend(ipoint+2) = 2    ! 3-node triangle
               isend(ipoint+3) = ntag
               if (ntag /= 1) then
                  isend(ipoint+4:ipoint+3+ntag) = grid%face(face)%tags(1:ntag)
               endif
               isend(ipoint+4) = grid%face(face)%bmark
               ipoint = ipoint + 3+ntag
               isend(ipoint+1:ipoint+3) = grid%face(face)%vertex
               ipoint = ipoint + 3
            end do
         endif
         elem = grid%element(elem)%next
      end do
   end do
   deallocate(seen_edge)
   if (global_element_kind == TETRAHEDRAL_ELEMENT) deallocate(seen_face)

   isend(1) = nelement
   ni = ipoint
   if (parallel == SEQUENTIAL) then
      irecv => isend
   else
      call phaml_send(procs,MASTER,isend,ni,(/0.0_my_real/),0,751)
   endif
else ! MASTER
   call phaml_recv(procs,proc,irecv,ni,rrecv,nr,751)
endif

if (iprint) then
   write(unit,"(A)") "$Elements"
   write(unit,"(I11)") irecv(1)
   fmat = "(    I11)"

   ipoint = 1
   do i=1,irecv(1)
      select case (irecv(ipoint+2)) ! element type
      case (15) ! point
         nval = 4+irecv(ipoint+3)
      case (1) ! line
         nval = 5+irecv(ipoint+3)
      case (2) ! triangle
         nval = 6+irecv(ipoint+3)
      case (4) ! tetrahedron
         nval = 7+irecv(ipoint+3)
      end select
      write(fmat(2:5),"(I3)") nval
      write(unit,fmat) irecv(ipoint+1:ipoint+nval)
      ipoint = ipoint + nval
   end do

   write(unit,"(A)") "$EndElements"
endif

if (parallel == SEQUENTIAL) then
   nullify(irecv)
   deallocate(isend)
else
   if (allocated(isend)) deallocate(isend)
   if (associated(irecv)) deallocate(irecv)
   if (associated(rrecv)) deallocate(rrecv)
endif

! output solution at vertices

if (fmt == GRIDFILE_MSH_SOLN) then

   nsolut = 1
   if (comp == -1)  nsolut = nsolut*size(grid%vertex_solution,dim=2)
   if (eigen == -1) nsolut = nsolut*size(grid%vertex_solution,dim=3)
   if (nsolut > 999) then
      call warning("Too many solutions per vertex in store_grid_msh", &
                   "Solutions not included in grid file")
      return
   endif

   if (myproc == 1) then
      allocate(isend(grid%nvert+1),rsend(nsolut*grid%nvert),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in store_grid_msh",procs=procs)
         stop
      endif
      isend(1) = grid%nvert
      ipoint = 1
      rpoint = 0
      visited_vert = .false.
      do lev=1,grid%nlev
         elem = grid%head_level_elem(lev)
         do while (elem /= END_OF_LIST)
            do ivert=1,VERTICES_PER_ELEMENT
               vert = grid%element(elem)%vertex(ivert)
               if (visited_vert(vert)) cycle
               visited_vert(vert) = .true.
               isend(ipoint+1) = vert
               ipoint = ipoint + 1
               if (comp == -1) then
                  if (eigen == -1) then
                     rsend(rpoint+1:rpoint+nsolut) = &
                        reshape(grid%vertex_solution(vert,:,:),(/nsolut/))
                  else
                     rsend(rpoint+1:rpoint+nsolut) = &
                        grid%vertex_solution(vert,:,eigen)
                  endif
               else
                  if (eigen == -1) then
                     rsend(rpoint+1:rpoint+nsolut) = &
                        grid%vertex_solution(vert,comp,:)
                  else
                     rsend(rpoint+1) = grid%vertex_solution(vert,comp,eigen)
                  endif
               endif
               rpoint = rpoint + nsolut
            end do
            elem = grid%element(elem)%next
         end do
      end do
      ni = ipoint
      nr = rpoint
      if (parallel == SEQUENTIAL) then
         irecv => isend
         rrecv => rsend
      else
         call phaml_send(procs,MASTER,isend,ni,rsend,nr,750)
      endif
   else ! MASTER
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,750)
   endif

   if (iprint) then
      fmat = "(SS,1P,I11,    E21.13E2)"
      write(fmat(12:15),"(I3)") nsolut

      write(unit,"(A)") "$NodeData"
      write(unit,"(I2)") 1
      write(unit,"(A)") "Solution"
      write(unit,"(I2)") 1
      write(unit,"(F4.1)") 0.0
      write(unit,"(I2)") 3
      write(unit,"(I2)") 0
      write(unit,"(I4)") nsolut
      write(unit,"(I11)") irecv(1)

      ipoint = 1
      rpoint = 0
      do i=2,irecv(1)+1
         write(unit,fmat) irecv(ipoint+1), rrecv(rpoint+1:rpoint+nsolut)
         ipoint = ipoint + 1
         rpoint = rpoint + nsolut
      end do
      write(unit,"(A)") "$EndNodeData"
   endif

   if (parallel == SEQUENTIAL) then
      nullify(irecv,rrecv)
      deallocate(isend,rsend)
   else
      if (allocated(isend)) deallocate(isend)
      if (allocated(rsend)) deallocate(rsend)
      if (associated(irecv)) deallocate(irecv)
      if (associated(rrecv)) deallocate(rrecv)
   endif

endif

end subroutine store_grid_msh

end module grid_io
