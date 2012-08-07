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

module evaluate

!----------------------------------------------------
! This module contains routines for evaluation of the solution.
!
! communication tags in this module are of the form 9xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use gridtype_mod
use grid_util
use basis_functions
!----------------------------------------------------

implicit none
private
public evaluate_soln, evaluate_soln_slave, evaluate_soln_local, &
       evaluate_oldsoln_local, set_grid_for_old_soln

!----------------------------------------------------
! The following variables are defined:

type(grid_type), pointer, save :: grid_for_old_soln
integer, parameter :: roundoff_fudge = 1000

!----------------------------------------------------

contains

!          -------------
subroutine evaluate_soln(procs,x,y,u,ux,uy,uxx,uyy,comp,eigen,z,uz,uzz)
!          -------------

!----------------------------------------------------
! This routine evaluates the solution and/or derivatives at the points in the
! arrays (x,y,z) and returns them in the u arrays.  The return value is 0.0 for
! points that are outside the domain.  This would only be called by
! the master or a slave requesting an evaluation by another slave.
! If comp and/or eigen is present, it tells which component and/or which
! eigenfunction to return.  The default is 1 for each.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(in) :: procs
real(my_real), intent(in) :: x(:),y(:)
real(my_real), optional, intent(out) :: u(:),ux(:),uy(:),uxx(:),uyy(:)
integer, optional, intent(in) :: comp,eigen
real(my_real), optional, intent(in) :: z(:)
real(my_real), optional, intent(out) :: uz(:),uzz(:)

!----------------------------------------------------
! Local variables:

integer, allocatable :: send_int(:)
real(my_real), allocatable :: send_real(:)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
integer :: i, j, ni, nr, proc, first_int, first_real, allocstat

!----------------------------------------------------
! Begin executable code

if (present(z) .neqv. global_element_kind == TETRAHEDRAL_ELEMENT) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("evaluate_soln: present(z) and TETRAHEDRAL_ELEMENT don't agree")
   stop
endif

! Invoke the slaves to evaluate the solution at the points in
! elements they own, and merge the results into a single array

! send the message to evaluate the solution and the points at which to evaluate

call pack_procs_size(this_processors_procs,ni,nr)
allocate(send_int(10+ni),send_real(nr+3*size(x)),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in evaluate_soln",procs=procs)
   return
endif
send_int(1) = 3 ! code for "evaluate"
if (present(comp)) then
   send_int(2) = comp
else
   send_int(2) = 1
endif
if (present(eigen)) then
   send_int(3) = eigen
else
   send_int(3) = 1
endif
send_int(4:10) = 0
if (present(u  )) send_int(4) = 1
if (present(ux )) send_int(5) = 1
if (present(uy )) send_int(6) = 1
if (present(uxx)) send_int(7) = 1
if (present(uyy)) send_int(8) = 1
if (present(uz )) send_int(9) = 1
if (present(uzz)) send_int(10)= 1
ni = ni+10
! a change in first_int or first_real requires same change in phaml_slave
first_int = 11
first_real = 1
call pack_procs(this_processors_procs,send_int,first_int,send_real,first_real)
send_real(nr+1:nr+size(x)) = x
send_real(nr+size(x)+1:nr+2*size(x)) = y
nr = nr+2*size(x)
if (present(z)) then
   send_real(nr+1:nr+size(x)) = z
   nr = nr + size(x)
endif
do proc=1,num_proc(procs)
   call phaml_send(procs,proc,send_int,ni,send_real,nr,101)
end do
deallocate(send_real,stat=allocstat)

! receive the solution and flags indicating which points were evaluated
! from each slave

if (present(u  )) u   = 0.0_my_real
if (present(ux )) ux  = 0.0_my_real
if (present(uy )) uy  = 0.0_my_real
if (present(uxx)) uxx = 0.0_my_real
if (present(uyy)) uyy = 0.0_my_real
if (present(uz )) uz  = 0.0_my_real
if (present(uzz)) uzz = 0.0_my_real

do i=1,num_proc(procs)
   call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,910)
   j = 0
   if (present(u  )) then
      where (recv_int == 1) u   = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(ux )) then
      where (recv_int == 1) ux  = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uy )) then
      where (recv_int == 1) uy  = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uxx)) then
      where (recv_int == 1) uxx = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uyy)) then
      where (recv_int == 1) uyy = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uz )) then
      where (recv_int == 1) uz  = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uzz)) then
      where (recv_int == 1) uzz = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   deallocate(recv_int,recv_real,stat=allocstat)
end do

end subroutine evaluate_soln

!          -------------------
subroutine evaluate_soln_slave(grid,invoker_procs,x,y,comp,eigen, &
                               u,ux,uy,uxx,uyy,which,z,uz,uzz)
!          -------------------

!----------------------------------------------------
! This is the slaves version of evaluate_soln.  If the u's are present, the
! solution is returned in them.  Otherwise, it sends the solutions indicated
! by which to the invoker.  Use of u's or which should be exclusive.
! comp and eigen tell which component and eigenfunction to return.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: invoker_procs
real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: comp, eigen
real(my_real), intent(out), optional :: u(:),ux(:),uy(:),uxx(:),uyy(:)
integer, optional, intent(in) :: which(:)
real(my_real), optional, intent(in) :: z(:)
real(my_real), optional, intent(out) :: uz(:), uzz(:)

!----------------------------------------------------
! Local variables:

real(my_real) :: loc_u(size(x)), loc_ux(size(x)), loc_uy(size(x)), &
                 loc_uz(size(x)), loc_uxx(size(x)),loc_uyy(size(x)), &
                 loc_uzz(size(x)), xc(VERTICES_PER_ELEMENT), &
                 yc(VERTICES_PER_ELEMENT), zc(VERTICES_PER_ELEMENT)
integer :: have_it(size(x))
integer :: i,j,k,elem,nbasis,astat,isub,nr
real(my_real), allocatable :: basis(:),basisx(:),basisy(:),basisz(:), &
                              basisxx(:),basisyy(:),basiszz(:)
real(my_real), allocatable :: send_real(:)
integer :: loc_which(7)

!----------------------------------------------------
! Begin executable code

if (present(z) .neqv. global_element_kind == TETRAHEDRAL_ELEMENT) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("evaluate_soln_slave: present(z) and TETRAHEDRAL_ELEMENT don't agree")
   stop
endif

! determine which solutions to evaluate

if (present(which)) then
   loc_which = which
else
   loc_which = 0
   if (present(u  )) loc_which(1) = 1
   if (present(ux )) loc_which(2) = 1
   if (present(uy )) loc_which(3) = 1
   if (present(uxx)) loc_which(4) = 1
   if (present(uyy)) loc_which(5) = 1
   if (present(uz )) loc_which(6) = 1
   if (present(uzz)) loc_which(7) = 1
endif

have_it = 0

! look for points in my partition

! For each point ...

do i=1,size(x)

! find the element it is in.  have_it==1 if I own the element.
! elem is 0 if I do not own it.

   select case (global_element_kind)
   case (TRIANGULAR_ELEMENT)
      call find_containing_element(x(i),y(i),have_it(i),elem,grid)
   case (TETRAHEDRAL_ELEMENT)
      call find_containing_element(x(i),y(i),have_it(i),elem,grid,z=z(i))
   end select

   if (have_it(i)==1) then

! evaluate the basis functions at this point

      nbasis = VERTICES_PER_ELEMENT
      do j=1,EDGES_PER_ELEMENT
         if (grid%edge(grid%element(elem)%edge(j))%degree > 1) then
            nbasis = nbasis + grid%edge(grid%element(elem)%edge(j))%degree - 1
         endif
      end do
      do j=1,FACES_PER_ELEMENT
         if (grid%face(grid%element(elem)%face(j))%degree > 2) then
            nbasis = nbasis + &
               ((grid%face(grid%element(elem)%face(j))%degree-2)* &
                (grid%face(grid%element(elem)%face(j))%degree-1))/2
         endif
      end do
      nbasis = nbasis + element_dof(grid%element(elem)%degree)
      allocate(basis(nbasis),basisx(nbasis),basisy(nbasis),basisz(nbasis), &
               basisxx(nbasis),basisyy(nbasis),basiszz(nbasis),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in evaluate_soln_slave")
         stop
      endif
      xc = grid%vertex(grid%element(elem)%vertex)%coord%x
      yc = grid%vertex(grid%element(elem)%vertex)%coord%y
      zc = zcoord(grid%vertex(grid%element(elem)%vertex)%coord)
      select case (global_element_kind)
      case (TRIANGULAR_ELEMENT)
         if (loc_which(4)==1 .or. loc_which(5)==1) then
            call p_hier_basis_func(x(i),y(i),xc,yc, &
                                   (/grid%edge(grid%element(elem)%edge)%degree,&
                                     grid%element(elem)%degree /), "a", &
                                     basis,basisx,basisy,basisxx,basisyy)
         elseif (loc_which(2)==1 .or. loc_which(3)==1) then
            call p_hier_basis_func(x(i),y(i),xc,yc, &
                                   (/grid%edge(grid%element(elem)%edge)%degree,&
                                     grid%element(elem)%degree /), "a", &
                                     basis,basisx,basisy)
         else
            call p_hier_basis_func(x(i),y(i),xc,yc, &
                                   (/grid%edge(grid%element(elem)%edge)%degree,&
                                     grid%element(elem)%degree /), "a", &
                                     basis)
         endif
      case (TETRAHEDRAL_ELEMENT)
         if (loc_which(4)==1 .or. loc_which(5)==1) then
            call p_hier_basis_func(x(i),y(i),z(i),xc,yc,zc, &
                                   (/grid%edge(grid%element(elem)%edge)%degree,&
                                     grid%face(grid%element(elem)%face)%degree,&
                                     grid%element(elem)%degree /), "a", &
                                     basis,basisx,basisy,basisz,basisxx, &
                                     basisyy,basiszz)
         elseif (loc_which(2)==1 .or. loc_which(3)==1) then
            call p_hier_basis_func(x(i),y(i),z(i),xc,yc,zc, &
                                   (/grid%edge(grid%element(elem)%edge)%degree,&
                                     grid%face(grid%element(elem)%face)%degree,&
                                     grid%element(elem)%degree /), "a", &
                                     basis,basisx,basisy,basisz)
         else
            call p_hier_basis_func(x(i),y(i),z(i),xc,yc,zc, &
                                   (/grid%edge(grid%element(elem)%edge)%degree,&
                                     grid%face(grid%element(elem)%face)%degree,&
                                     grid%element(elem)%degree /), "a", &
                                     basis)
         endif
      end select

! compute the solution at this point

      loc_u  (i) = 0.0_my_real
      loc_ux (i) = 0.0_my_real
      loc_uy (i) = 0.0_my_real
      loc_uz (i) = 0.0_my_real
      loc_uxx(i) = 0.0_my_real
      loc_uyy(i) = 0.0_my_real
      loc_uzz(i) = 0.0_my_real

      do j=1,VERTICES_PER_ELEMENT
         if (loc_which(1)==1) loc_u  (i) = loc_u  (i) + &
              basis  (j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(2)==1) loc_ux (i) = loc_ux (i) + &
              basisx (j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(3)==1) loc_uy (i) = loc_uy (i) + &
              basisy (j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(4)==1) loc_uxx(i) = loc_uxx(i) + &
              basisxx(j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(5)==1) loc_uyy(i) = loc_uyy(i) + &
              basisyy(j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(6)==1) loc_uz(i) = loc_uz(i) + &
              basisz (j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(7)==1) loc_uzz(i) = loc_uzz(i) + &
              basiszz(j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
      end do
      isub = VERTICES_PER_ELEMENT
      do j=1,EDGES_PER_ELEMENT
         if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
         do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
            isub = isub + 1
            if (loc_which(1)==1) loc_u  (i) = loc_u  (i) + &
               basis  (isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(2)==1) loc_ux (i) = loc_ux (i) + &
               basisx (isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(3)==1) loc_uy (i) = loc_uy (i) + &
               basisy (isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(4)==1) loc_uxx(i) = loc_uxx(i) + &
               basisxx(isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(5)==1) loc_uyy(i) = loc_uyy(i) + &
               basisyy(isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(6)==1) loc_uz (i) = loc_uz (i) + &
               basisz (isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(7)==1) loc_uzz(i) = loc_uzz(i) + &
               basiszz(isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
         end do
      end do
      do j=1,FACES_PER_ELEMENT
         if (grid%face(grid%element(elem)%face(j))%degree < 3) cycle
         do k=1,((grid%face(grid%element(elem)%face(j))%degree-2) * &
                 (grid%face(grid%element(elem)%face(j))%degree-1))/2
            isub = isub + 1
            if (loc_which(1)==1) loc_u  (i) = loc_u  (i) + &
               basis  (isub)*grid%face(grid%element(elem)%face(j))%solution(k, &
               comp,eigen)
            if (loc_which(2)==1) loc_ux (i) = loc_ux (i) + &
               basisx (isub)*grid%face(grid%element(elem)%face(j))%solution(k, &
               comp,eigen)
            if (loc_which(3)==1) loc_uy (i) = loc_uy (i) + &
               basisy (isub)*grid%face(grid%element(elem)%face(j))%solution(k, &
               comp,eigen)
            if (loc_which(4)==1) loc_uxx(i) = loc_uxx(i) + &
               basisxx(isub)*grid%face(grid%element(elem)%face(j))%solution(k, &
               comp,eigen)
            if (loc_which(5)==1) loc_uyy(i) = loc_uyy(i) + &
               basisyy(isub)*grid%face(grid%element(elem)%face(j))%solution(k, &
               comp,eigen)
            if (loc_which(6)==1) loc_uz (i) = loc_uz (i) + &
               basisz (isub)*grid%face(grid%element(elem)%face(j))%solution(k, &
               comp,eigen)
            if (loc_which(7)==1) loc_uzz(i) = loc_uzz(i) + &
               basiszz(isub)*grid%face(grid%element(elem)%face(j))%solution(k, &
               comp,eigen)
         end do
      end do
      if (grid%element(elem)%degree > 2) then
         do k=1,element_dof(grid%element(elem)%degree)
            isub = isub + 1
            if (loc_which(1)==1) loc_u  (i) = loc_u  (i) + &
               basis  (isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(2)==1) loc_ux (i) = loc_ux (i) + &
               basisx (isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(3)==1) loc_uy (i) = loc_uy (i) + &
               basisy (isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(4)==1) loc_uxx(i) = loc_uxx(i) + &
               basisxx(isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(5)==1) loc_uyy(i) = loc_uyy(i) + &
               basisyy(isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(6)==1) loc_uz (i) = loc_uz (i) + &
               basisz (isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(7)==1) loc_uzz(i) = loc_uzz(i) + &
               basiszz(isub)*grid%element(elem)%solution(k,comp,eigen)
         end do
      endif

   else

      loc_u  (i) = 0.0_my_real
      loc_ux (i) = 0.0_my_real
      loc_uy (i) = 0.0_my_real
      loc_uxx(i) = 0.0_my_real
      loc_uyy(i) = 0.0_my_real
      loc_uz (i) = 0.0_my_real
      loc_uzz(i) = 0.0_my_real

   endif
   deallocate(basis,basisx,basisy,basisz,basisxx,basisyy,basiszz,stat=astat)
end do

! Copy the solution to the u's or send the solution to the invoker.
! Note that my_proc(invoker_procs) is not me, it is the invoker.

if (present(which)) then
   nr = size(x)*count(loc_which==1)
   allocate(send_real(nr),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in evaluate_soln_slave")
      stop
   endif
   j = 0
   if (loc_which(1) == 1) then
      send_real(j+1:j+size(x)) = loc_u
      j = j + size(x)
   endif
   if (loc_which(2) == 1) then
      send_real(j+1:j+size(x)) = loc_ux
      j = j + size(x)
   endif
   if (loc_which(3) == 1) then
      send_real(j+1:j+size(x)) = loc_uy
      j = j + size(x)
   endif
   if (loc_which(4) == 1) then
      send_real(j+1:j+size(x)) = loc_uxx
      j = j + size(x)
   endif
   if (loc_which(5) == 1) then
      send_real(j+1:j+size(x)) = loc_uyy
      j = j + size(x)
   endif
   if (loc_which(6) == 1) then
      send_real(j+1:j+size(x)) = loc_uz
      j = j + size(x)
   endif
   if (loc_which(7) == 1) then
      send_real(j+1:j+size(x)) = loc_uzz
      j = j + size(x)
   endif
   call phaml_send(invoker_procs,my_proc(invoker_procs),have_it,size(x), &
                   send_real,nr,910)
   deallocate(send_real)
else
   if (present(u  )) u   = loc_u
   if (present(ux )) ux  = loc_ux
   if (present(uy )) uy  = loc_uy
   if (present(uxx)) uxx = loc_uxx
   if (present(uyy)) uyy = loc_uyy
   if (present(uz )) uz  = loc_uz
   if (present(uzz)) uzz = loc_uzz
endif

end subroutine evaluate_soln_slave

!          -------------------
subroutine evaluate_soln_local(grid,x,y,elem,comp,eigen,u,ux,uy,uxx,uyy,uxy, &
                               z,uz,uzz,uxz,uyz)
!          -------------------

!----------------------------------------------------
! This routine evaluates the solution and/or derivatives of the solution
! on this processor.  The points (x,y) must all be in element elem.
! comp and eigen tell which components and eigenfunctions to return.
! The three dimensions of u correspond to comp, eigen and point.
! If ux or uy is present, then both must be present.  If uxx or uyy is
! present, then all 5 must be present.  This is only to reduce
! the number of forms of the call to basis_function, so this restriction can
! easily be removed.  Likewise, if uxy is present then u, ux and uy must be
! present and uxx and uyy either both present or not present.  In 3D, these
! rules are extended to include the z derivatives.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: elem
integer, intent(in) :: comp(:),eigen(:)
real(my_real), optional, intent(out) :: u(:,:,:),ux(:,:,:),uy(:,:,:), &
                                        uxx(:,:,:),uyy(:,:,:),uxy(:,:,:)
! dimensions are (comp,eigen,point)
real(my_real), optional, intent(in) :: z(:)
real(my_real), optional, intent(out) :: uz(:,:,:),uzz(:,:,:),uxz(:,:,:), &
                                        uyz(:,:,:)

!----------------------------------------------------
! Local variables:

real(my_real) :: xc(VERTICES_PER_ELEMENT), yc(VERTICES_PER_ELEMENT), &
                 zc(VERTICES_PER_ELEMENT)
real(my_real), allocatable :: basis(:,:),basisx(:,:),basisy(:,:),basisxx(:,:), &
                              basisyy(:,:),basisxy(:,:),solnvect(:,:), &
                              basisz(:,:),basiszz(:,:),basisxz(:,:), &
                              basisyz(:,:)
integer :: i,j,k,l,p,nbasis,astat,isub
logical :: useblas

!----------------------------------------------------
! Begin executable code

if (present(z) .neqv. global_element_kind == TETRAHEDRAL_ELEMENT) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("evaluate_soln_local: present(z) and TETRAHEDRAL_ELEMENT don't agree")
   stop
endif

! evaluate the basis functions at these points

nbasis = VERTICES_PER_ELEMENT
do j=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(j))%degree > 1) then
      nbasis = nbasis + grid%edge(grid%element(elem)%edge(j))%degree - 1
   endif
end do
do j=1,FACES_PER_ELEMENT
   if (grid%face(grid%element(elem)%face(j))%degree > 2) then
      nbasis = nbasis + ((grid%edge(grid%element(elem)%edge(j))%degree-2) * &
                         (grid%edge(grid%element(elem)%edge(j))%degree-1))/2
   endif
end do
if (grid%element(elem)%degree > 2) then
   nbasis = nbasis + element_dof(grid%element(elem)%degree)
endif
allocate(basis(nbasis,size(x)),basisx(nbasis,size(x)),basisy(nbasis,size(x)), &
         basisxx(nbasis,size(x)),basisyy(nbasis,size(x)), &
         basisxy(nbasis,size(x)),basisz(nbasis,size(x)), &
         basiszz(nbasis,size(x)),basisxz(nbasis,size(x)), &
         basisyz(nbasis,size(x)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in evaluate_soln_local")
   stop
endif

select case (global_element_kind)

case (TRIANGULAR_ELEMENT)

if (present(ux) .neqv. present(uy)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("need both or neither of ux and uy in evaluate_soln_local")
   stop
endif
if (present(uxx) .or. present(uyy)) then
   if (.not. present(u) .or. .not. present(ux) .or. .not. present(uy) .or. &
       .not. present(uxx) .or. .not. present(uyy)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("if uxx or uyy is present in evaluate_soln_local then all args except uxy must be present")
      stop
   endif
endif
if (present(uxy)) then
   if (.not. present(u) .or. .not. present(ux) .or. .not. present(uy) .or. &
       (present(uxx) .neqv. present(uyy))) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("evaluate_soln_local: uxy present requires u, ux and uy present and uxx and uyy both or neither present")
      stop
   endif
endif

case (TETRAHEDRAL_ELEMENT)

if ((present(ux) .neqv. present(uy)) .or. &
    (present(ux) .neqv. present(uz))) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("need all or none of ux, uy and uz in evaluate_soln_local")
   stop
endif
if (present(uxx) .or. present(uyy) .or. present(uxx)) then
   if (.not. present(u) .or. .not. present(ux) .or. .not. present(uy) .or. &
       .not. present(uxx) .or. .not. present(uyy) .or. .not. present(uz) .or. &
       .not. present(uzz)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("if uxx, uyy or uzz is present in evaluate_soln_local then all args except mixed derivatives must be present")
      stop
   endif
endif
if (present(uxy) .or. present(uxz) .or. present(uyz)) then
   if (.not. present(u) .or. .not. present(ux) .or. .not. present(uy) .or. &
       (present(uxx) .neqv. present(uyy)) .or. &
       (present(uxx) .neqv. present(uzz)) .or. .not. present(uxy) .or. &
       .not. present(uxz) .or. .not. present(uyz)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("evaluate_soln_local: uxy, uxz or uyz present requires u,"// &
        " ux, uy and uz present and uxx, uyy and uzz all present or all absent")
      stop
   endif
endif

end select

xc = grid%vertex(grid%element(elem)%vertex)%coord%x
yc = grid%vertex(grid%element(elem)%vertex)%coord%y
zc = zcoord(grid%vertex(grid%element(elem)%vertex)%coord)
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   if (present(u) .and. .not. present(ux)) then
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis)
   elseif (.not. present(u) .and. present(ux)) then
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basisx=basisx,basisy=basisy)
   elseif (present(uxy) .and. present(uxx)) then
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisxx,basisyy,basisxy)
   elseif (present(uxy) .and. .not. present(uxx)) then
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisxy=basisxy)
   elseif (present(uxx)) then
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisxx,basisyy)
   else
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                              grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy)
   endif
case (TETRAHEDRAL_ELEMENT)
   if (present(u) .and. .not. present(ux)) then
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis)
   elseif (.not. present(u) .and. present(ux)) then
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basisx=basisx,basisy=basisy,basisz=basisz)
   elseif (present(uxy) .and. present(uxx)) then
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisz,basisxx,basisyy, &
                             basiszz,basisxy,basisxz,basisyz)
   elseif (present(uxy) .and. .not. present(uxx)) then
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisz,basisxy=basisxy, &
                             basisxz=basisxz,basisyz=basisyz)
   elseif (present(uxx)) then
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisz,basisxx,basisyy, &
                             basiszz)
   else
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                              grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisz)
   endif
end select

! compute the solution at this point

if (present(u)) u = 0.0_my_real
if (present(ux)) ux = 0.0_my_real
if (present(uy)) uy = 0.0_my_real
if (present(uxx)) uxx = 0.0_my_real
if (present(uyy)) uyy = 0.0_my_real
if (present(uxy)) uxy = 0.0_my_real
if (present(uz)) uz = 0.0_my_real
if (present(uzz)) uzz = 0.0_my_real
if (present(uxz)) uxz = 0.0_my_real
if (present(uyz)) uyz = 0.0_my_real

! special BLAS code for size(comp)==size(eigen)==1.  Copy solution values to
! an array and use GEMM.  TEMP Probably can do something similar with /=1.
! Also need the first 2 dimensions of the u's to be 1 since we're
! (illegally!) passing a rank 3 array to a rank 2 dummy argument.

useblas = size(comp)==1 .and. size(eigen)==1
if (present(u)) useblas = useblas .and. size(u,dim=1) == 1 .and. &
                          size(u,dim=2) == 1
if (present(ux)) useblas = useblas .and. size(ux,dim=1) == 1 .and. &
                          size(ux,dim=2) == 1
if (present(uy)) useblas = useblas .and. size(uy,dim=1) == 1 .and. &
                          size(uy,dim=2) == 1
if (present(uxx)) useblas = useblas .and. size(uxx,dim=1) == 1 .and. &
                          size(uxx,dim=2) == 1
if (present(uyy)) useblas = useblas .and. size(uyy,dim=1) == 1 .and. &
                          size(uyy,dim=2) == 1
if (present(uxy)) useblas = useblas .and. size(uxy,dim=1) == 1 .and. &
                          size(uxy,dim=2) == 1
if (present(uz)) useblas = useblas .and. size(uz,dim=1) == 1 .and. &
                          size(uz,dim=2) == 1
if (present(uzz)) useblas = useblas .and. size(uzz,dim=1) == 1 .and. &
                          size(uzz,dim=2) == 1
if (present(uxz)) useblas = useblas .and. size(uxz,dim=1) == 1 .and. &
                          size(uxz,dim=2) == 1
if (present(uyz)) useblas = useblas .and. size(uyz,dim=1) == 1 .and. &
                          size(uyz,dim=2) == 1

if (useblas) then
   allocate(solnvect(1,nbasis),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in evaluate_soln_local")
      stop
   endif

   do j=1,VERTICES_PER_ELEMENT
      solnvect(1,j) = grid%vertex_solution(grid%element(elem)%vertex(j), &
                                           comp(1),eigen(1))
   end do
   isub = VERTICES_PER_ELEMENT
   do j=1,EDGES_PER_ELEMENT
      if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
         isub = isub + 1
         solnvect(1,isub) = grid%edge(grid%element(elem)%edge(j))%solution(k, &
                            comp(1),eigen(1))
      end do
   end do
   do j=1,FACES_PER_ELEMENT
      if (grid%face(grid%element(elem)%face(j))%degree < 3) cycle
      do k=1,((grid%face(grid%element(elem)%face(j))%degree-2) * &
              (grid%face(grid%element(elem)%face(j))%degree-1))/2
         isub = isub + 1
         solnvect(1,isub) = grid%face(grid%element(elem)%face(j))%solution(k, &
                            comp(1),eigen(1))
      end do
   end do
   if (grid%element(elem)%degree > 2) then
      do k=1,element_dof(grid%element(elem)%degree)
         isub = isub + 1
         solnvect(1,isub) = grid%element(elem)%solution(k,comp(1),eigen(1))
      end do
   endif

   if (my_real == kind(1.0)) then
      if (present(u)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                 basis,nbasis,1.0,u,1)
      if (present(ux)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisx,nbasis,1.0,ux,1)
      if (present(uy)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisy,nbasis,1.0,uy,1)
      if (present(uxx)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisxx,nbasis,1.0,uxx,1)
      if (present(uyy)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisyy,nbasis,1.0,uyy,1)
      if (present(uxy)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisxy,nbasis,1.0,uxy,1)
      if (present(uz)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisz,nbasis,1.0,uz,1)
      if (present(uzz)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basiszz,nbasis,1.0,uzz,1)
      if (present(uxz)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisxz,nbasis,1.0,uxz,1)
      if (present(uyz)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisyz,nbasis,1.0,uyz,1)
   elseif (my_real == kind(1.0d0)) then
      if (present(u)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                 basis,nbasis,1.0d0,u,1)
      if (present(ux)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisx,nbasis,1.0d0,ux,1)
      if (present(uy)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisy,nbasis,1.0d0,uy,1)
      if (present(uxx)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisxx,nbasis,1.0d0,uxx,1)
      if (present(uyy)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisyy,nbasis,1.0d0,uyy,1)
      if (present(uxy)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisxy,nbasis,1.0d0,uxy,1)
      if (present(uz)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisz,nbasis,1.0d0,uz,1)
      if (present(uzz)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basiszz,nbasis,1.0d0,uzz,1)
      if (present(uxz)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisxz,nbasis,1.0d0,uxz,1)
      if (present(uyz)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisyz,nbasis,1.0d0,uyz,1)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call GEMM")
      stop
   endif

   deallocate(solnvect,stat=astat)

else ! general case with size(comp)/=1 .or. size(eigen)/=1

   do j=1,VERTICES_PER_ELEMENT
      do p=1,size(x)
       do i=1,size(comp)
        do l=1,size(eigen)
         if (present(u)) u(i,l,p) = u(i,l,p) + &
            basis(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                       comp(i),eigen(l))
         if (present(ux)) ux(i,l,p) = ux(i,l,p) + &
            basisx(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uy)) uy(i,l,p) = uy(i,l,p) + &
            basisy(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uxx)) uxx(i,l,p) = uxx(i,l,p) + &
            basisxx(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uyy)) uyy(i,l,p) = uyy(i,l,p) + &
            basisyy(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uxy)) uxy(i,l,p) = uxy(i,l,p) + &
            basisxy(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uz)) uz(i,l,p) = uz(i,l,p) + &
            basisz(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uzz)) uzz(i,l,p) = uzz(i,l,p) + &
            basiszz(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uxz)) uxz(i,l,p) = uxz(i,l,p) + &
            basisxz(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uyz)) uyz(i,l,p) = uyz(i,l,p) + &
            basisyz(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
        end do
       end do
      end do
   end do
   isub = VERTICES_PER_ELEMENT
   do j=1,EDGES_PER_ELEMENT
      if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
         isub = isub + 1
         do p=1,size(x)
          do i=1,size(comp)
           do l=1,size(eigen)
            if (present(u)) u(i,l,p) = u(i,l,p) + &
               basis(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
                             comp(i),eigen(l))
            if (present(ux)) ux(i,l,p) = ux(i,l,p) + &
               basisx(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uy)) uy(i,l,p) = uy(i,l,p) + &
               basisy(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxx)) uxx(i,l,p) = uxx(i,l,p) + &
              basisxx(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uyy)) uyy(i,l,p) = uyy(i,l,p) + &
              basisyy(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxy)) uxy(i,l,p) = uxy(i,l,p) + &
              basisxy(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uz)) uz(i,l,p) = uz(i,l,p) + &
               basisz(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uzz)) uzz(i,l,p) = uzz(i,l,p) + &
              basiszz(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxz)) uxz(i,l,p) = uxz(i,l,p) + &
              basisxz(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uyz)) uyz(i,l,p) = uyz(i,l,p) + &
              basisyz(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
           end do
          end do
         end do
      end do
   end do
   do j=1,FACES_PER_ELEMENT
      if (grid%face(grid%element(elem)%face(j))%degree < 3) cycle
      do k=1,((grid%face(grid%element(elem)%face(j))%degree-2) * &
              (grid%face(grid%element(elem)%face(j))%degree-1))/2
         isub = isub + 1
         do p=1,size(x)
          do i=1,size(comp)
           do l=1,size(eigen)
            if (present(u)) u(i,l,p) = u(i,l,p) + &
               basis(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k, &
                             comp(i),eigen(l))
            if (present(ux)) ux(i,l,p) = ux(i,l,p) + &
               basisx(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uy)) uy(i,l,p) = uy(i,l,p) + &
               basisy(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxx)) uxx(i,l,p) = uxx(i,l,p) + &
              basisxx(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uyy)) uyy(i,l,p) = uyy(i,l,p) + &
              basisyy(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxy)) uxy(i,l,p) = uxy(i,l,p) + &
              basisxy(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uz)) uz(i,l,p) = uz(i,l,p) + &
               basisz(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uzz)) uzz(i,l,p) = uzz(i,l,p) + &
              basiszz(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxz)) uxz(i,l,p) = uxz(i,l,p) + &
              basisxz(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uyz)) uyz(i,l,p) = uyz(i,l,p) + &
              basisyz(isub,p)*grid%face(grid%element(elem)%face(j))%solution(k,&
                              comp(i),eigen(l))
           end do
          end do
         end do
      end do
   end do
   if (grid%element(elem)%degree > 2) then
      do k=1,element_dof(grid%element(elem)%degree)
         isub = isub + 1
         do p=1,size(x)
          do i=1,size(comp)
           do l=1,size(eigen)
            if (present(u)) u(i,l,p) = u(i,l,p) + &
               basis(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(ux)) ux(i,l,p) = ux(i,l,p) + &
               basisx(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uy)) uy(i,l,p) = uy(i,l,p) + &
               basisy(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uxx)) uxx(i,l,p) = uxx(i,l,p) + &
               basisxx(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uyy)) uyy(i,l,p) = uyy(i,l,p) + &
               basisyy(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uxy)) uxy(i,l,p) = uxy(i,l,p) + &
               basisxy(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uz)) uz(i,l,p) = uz(i,l,p) + &
               basisz(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uzz)) uzz(i,l,p) = uzz(i,l,p) + &
               basiszz(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uxz)) uxz(i,l,p) = uxz(i,l,p) + &
               basisxz(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uyz)) uyz(i,l,p) = uyz(i,l,p) + &
               basisyz(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
           end do
          end do
         end do
      end do
   endif

endif ! size(comp)==size(eigen)==1

deallocate(basis,basisx,basisy,basisxx,basisyy,basisxy,basisz,basiszz, &
           basisxz,basisyz,stat=astat)

end subroutine evaluate_soln_local

!          ----------------------
subroutine evaluate_oldsoln_local(x,y,u,ux,uy,uxx,uyy,comp,eigen,z,uz,uzz)
!          ----------------------

!----------------------------------------------------
! This routine evaluates the old solution and/or derivatives of the old solution
! on this processor.
! If present, comp and/or eigen tell which component and eigenfunction to return
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
real(my_real), optional, intent(out) :: u,ux,uy,uxx,uyy
integer, optional, intent(in) :: comp, eigen
real(my_real), optional, intent(in) :: z
real(my_real), optional, intent(out) :: uz,uzz

!----------------------------------------------------
! Local variables:

integer :: j,k,nbasis,astat,isub,elem,loc_comp,loc_eigen
real(my_real), allocatable :: basis(:),basisx(:),basisy(:),basisxx(:), &
                              basisyy(:),basisz(:),basiszz(:)
real(my_real) :: xc(VERTICES_PER_ELEMENT), yc(VERTICES_PER_ELEMENT), &
                 zc(VERTICES_PER_ELEMENT), solnval
type(grid_type), pointer :: grid

!----------------------------------------------------
! Begin executable code

if (present(z) .neqv. global_element_kind == TETRAHEDRAL_ELEMENT) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("evaluate_oldsoln_local: present(z) and TETRAHEDRAL_ELEMENT don't agree")
   stop
endif

! make sure the old solution exists

grid => grid_for_old_soln
if (.not. grid%oldsoln_exists) then
   call warning("Old solution has not been saved, and hence cannot be evaluated", &
                "Returning 0.0 as old solution.")
   if (present(u)) u=0.0_my_real
   if (present(ux)) ux=0.0_my_real
   if (present(uy)) uy=0.0_my_real
   if (present(uxx)) uxx=0.0_my_real
   if (present(uyy)) uyy=0.0_my_real
   if (present(uz)) uz=0.0_my_real
   if (present(uzz)) uzz=0.0_my_real
   return
endif

! determine the component and eigenfunction

if (present(comp)) then
   if (comp <= 0) then
      call fatal("phaml_evaluate_oldsoln: comp must be a positive integer.  It is ", &
                  intlist=(/comp/))
      stop
   endif
   if (comp > grid%system_size) then
      call fatal("phaml_evaluate_oldsoln: comp must be no larger than system size.  They are ", &
                 intlist=(/comp,grid%system_size/))
      stop
   endif
   loc_comp = comp
else
   loc_comp = 1
endif
if (present(eigen)) then
   if (eigen <= 0) then
      call fatal("phaml_evaluate_oldsoln: eigen must be a positive integer.  It is ", &
                 intlist=(/eigen/))
      stop
   endif
   if (eigen > max(1,grid%num_eval)) then
      call fatal("phaml_evaluate_oldsoln: eigen must be no larger than the number of eigenvalues computed.  They are ", &
                 intlist=(/eigen,grid%num_eval/))
      stop
   endif
   loc_eigen = eigen
else
   loc_eigen = 1 
endif

! find the element containing (x,y)

call find_old_containing_element(x,y,elem,grid,z)
if (elem == 0) then
   if (present(u)) u=0.0_my_real
   if (present(ux)) ux=0.0_my_real
   if (present(uy)) uy=0.0_my_real
   if (present(uxx)) uxx=0.0_my_real
   if (present(uyy)) uyy=0.0_my_real
   if (present(uz)) uz=0.0_my_real
   if (present(uzz)) uzz=0.0_my_real
   return
endif

! evaluate the basis functions at the point

nbasis = VERTICES_PER_ELEMENT
do j=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(j))%degree > 1) then
      nbasis = nbasis + grid%edge(grid%element(elem)%edge(j))%degree - 1
   endif
end do
do j=1,FACES_PER_ELEMENT
   if (grid%face(grid%element(elem)%face(j))%degree > 2) then
      nbasis = nbasis + ((grid%face(grid%element(elem)%face(j))%degree-2) * &
                         (grid%face(grid%element(elem)%face(j))%degree-1))/2
   endif
end do
if (grid%element(elem)%degree > 2) then
   nbasis = nbasis + element_dof(grid%element(elem)%degree)
endif
allocate(basis(nbasis),basisx(nbasis),basisy(nbasis),basisxx(nbasis), &
         basisyy(nbasis),basisz(nbasis),basiszz(nbasis),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in evaluate_oldsoln_local")
   stop
endif
xc = grid%vertex(grid%element(elem)%vertex)%coord%x
yc = grid%vertex(grid%element(elem)%vertex)%coord%y
zc = zcoord(grid%vertex(grid%element(elem)%vertex)%coord)
select case (global_element_kind)
case (TRIANGULAR_ELEMENT)
   if (present(uxx) .or. present(uyy) .or. present(uzz)) then
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                              grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisxx,basisyy)
   elseif (present(ux) .or. present(uy) .or. present(uz)) then
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy)
   else
      call p_hier_basis_func(x,y,xc,yc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                              grid%element(elem)%degree /),"a", &
                             basis)
   endif
case (TETRAHEDRAL_ELEMENT)
   if (present(uxx) .or. present(uyy) .or. present(uzz)) then
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                              grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisz,basisxx,basisyy, &
                             basiszz)
   elseif (present(ux) .or. present(uy) .or. present(uz)) then
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                               grid%element(elem)%degree /),"a", &
                             basis,basisx,basisy,basisz)
   else
      call p_hier_basis_func(x,y,z,xc,yc,zc, &
                             (/grid%edge(grid%element(elem)%edge)%degree, &
                               grid%face(grid%element(elem)%face)%degree, &
                              grid%element(elem)%degree /),"a", &
                             basis)
   endif
end select

! compute the solution at this point

if (present(u)) u = 0.0_my_real
if (present(ux)) ux = 0.0_my_real
if (present(uy)) uy = 0.0_my_real
if (present(uxx)) uxx = 0.0_my_real
if (present(uyy)) uyy = 0.0_my_real
if (present(uz)) uz = 0.0_my_real
if (present(uzz)) uzz = 0.0_my_real

do j=1,VERTICES_PER_ELEMENT
   solnval = grid%vertex_oldsoln(grid%element(elem)%vertex(j),loc_comp,loc_eigen)
   if (present(u)) u = u + solnval*basis(j)
   if (present(ux)) ux = ux + solnval*basisx(j)
   if (present(uy)) uy = uy + solnval*basisy(j)
   if (present(uxx)) uxx = uxx + solnval*basisxx(j)
   if (present(uyy)) uyy = uyy + solnval*basisyy(j)
   if (present(uz)) uz = uz + solnval*basisz(j)
   if (present(uzz)) uzz = uzz + solnval*basiszz(j)
end do
isub = VERTICES_PER_ELEMENT
do j=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
   do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
      isub = isub + 1
      if (.not. associated(grid%edge(grid%element(elem)%edge(j))%oldsoln)) cycle
      if (k > size(grid%edge(grid%element(elem)%edge(j))%oldsoln,dim=1)) cycle
      solnval = grid%edge(grid%element(elem)%edge(j))%oldsoln(k,loc_comp,loc_eigen)
      if (present(u)) u = u + solnval*basis(isub)
      if (present(ux)) ux = ux + solnval*basisx(isub)
      if (present(uy)) uy = uy + solnval*basisy(isub)
      if (present(uxx)) uxx = uxx + solnval*basisxx(isub)
      if (present(uyy)) uyy = uyy + solnval*basisyy(isub)
      if (present(uz)) uz = uz + solnval*basisz(isub)
      if (present(uzz)) uzz = uzz + solnval*basiszz(isub)
   end do
end do
do j=1,FACES_PER_ELEMENT
   if (grid%face(grid%element(elem)%face(j))%degree < 3) cycle
   do k=1,((grid%face(grid%element(elem)%face(j))%degree-2) * &
           (grid%face(grid%element(elem)%face(j))%degree-1))/2
      isub = isub + 1
      if (.not. associated(grid%face(grid%element(elem)%face(j))%oldsoln)) cycle
      if (k > size(grid%face(grid%element(elem)%face(j))%oldsoln,dim=1)) cycle
      solnval = grid%face(grid%element(elem)%face(j))%oldsoln(k,loc_comp,loc_eigen)
      if (present(u)) u = u + solnval*basis(isub)
      if (present(ux)) ux = ux + solnval*basisx(isub)
      if (present(uy)) uy = uy + solnval*basisy(isub)
      if (present(uxx)) uxx = uxx + solnval*basisxx(isub)
      if (present(uyy)) uyy = uyy + solnval*basisyy(isub)
      if (present(uz)) uz = uz + solnval*basisz(isub)
      if (present(uzz)) uzz = uzz + solnval*basiszz(isub)
   end do
end do
if (grid%element(elem)%degree > 2) then
   do k=1,element_dof(grid%element(elem)%degree)
      isub = isub + 1
      if (.not. associated(grid%element(elem)%oldsoln)) cycle
      solnval = grid%element(elem)%oldsoln(k,loc_comp,loc_eigen)
      if (present(u)) u = u + solnval*basis(isub)
      if (present(ux)) ux = ux + solnval*basisx(isub)
      if (present(uy)) uy = uy + solnval*basisy(isub)
      if (present(uxx)) uxx = uxx + solnval*basisxx(isub)
      if (present(uyy)) uyy = uyy + solnval*basisyy(isub)
      if (present(uz)) uz = uz + solnval*basisz(isub)
      if (present(uzz)) uzz = uzz + solnval*basiszz(isub)
   end do
endif

deallocate(basis,basisx,basisy,basisxx,basisyy,basisz,basiszz,stat=astat)

end subroutine evaluate_oldsoln_local

!          ---------------------------
subroutine find_old_containing_element(x,y,elem,grid,z)
!          ---------------------------

!----------------------------------------------------
! This routine looks for an "effective old leaf" element containing the point
! (x,y).  An element is an effective old leaf if oldsoln is allocated.
! If the point falls on an element boundary, it is indeterminant as to which
! of the containing elements it returns.
! If the point does not fall in any element, 0 is returned.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(out) :: elem
type(grid_type), intent(in) :: grid
real(my_real), optional, intent(in) :: z

!----------------------------------------------------
! Local variables:

integer :: e,root,i,neigh(NEIGHBORS_PER_ELEMENT)
real(my_real) :: bc(VERTICES_PER_ELEMENT,1)
!----------------------------------------------------
! Begin executable code

! find an element in the initial grid that contains (x,y)

! first look for it by moving from a triangle to a neighbor in the direction
! of the point, which is characterized by a negative barycentric coordinate

root = -10
e = grid%head_level_elem(1)
do
! if all barycentric coordinates are positive, it's in there
   select case (global_element_kind)
   case (TRIANGULAR_ELEMENT)
      call barycentric((/x/),(/y/), &
                       grid%vertex(grid%element(e)%vertex)%coord%x, &
                       grid%vertex(grid%element(e)%vertex)%coord%y, &
                       bc,no_det=.true.)
   case (TETRAHEDRAL_ELEMENT)
      call barycentric_tet((/x/),(/y/),(/z/), &
                           grid%vertex(grid%element(e)%vertex)%coord%x, &
                           grid%vertex(grid%element(e)%vertex)%coord%y, &
                           zcoord(grid%vertex(grid%element(e)%vertex)%coord), &
                           bc)
   end select
   if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
      root = e
      exit
   endif
   neigh = grid%initial_neighbor(:,e)
   do i=1,VERTICES_PER_ELEMENT
      if (bc(i,1) < 0) then
         e = neigh(i)
         if (e /= BOUNDARY) exit
      endif
   end do
   if (e == BOUNDARY) exit
end do

! if that failed, go through all the initial elements

if (root == -10) then
   e = grid%head_level_elem(1)
   do while (e /= END_OF_LIST)
      select case (global_element_kind)
      case (TRIANGULAR_ELEMENT)
         call barycentric((/x/),(/y/), &
                          grid%vertex(grid%element(e)%vertex)%coord%x, &
                          grid%vertex(grid%element(e)%vertex)%coord%y, &
                          bc,no_det=.true.)
      case (TETRAHEDRAL_ELEMENT)
         call barycentric_tet((/x/),(/y/),(/z/), &
                              grid%vertex(grid%element(e)%vertex)%coord%x, &
                              grid%vertex(grid%element(e)%vertex)%coord%y, &
                            zcoord(grid%vertex(grid%element(e)%vertex)%coord), &
                              bc)
      end select
      if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
         root = e
         exit
      endif
      e = grid%element(e)%next
   end do
endif

! go down the refinement tree to find the effective leaf element that
! contains (x,y)

if (root /= -10) then
   call find_old_containing_leaf(x,y,root,elem,grid,z)
else
   elem = 0
endif

end subroutine find_old_containing_element

!                    ------------------------
recursive subroutine find_old_containing_leaf(x,y,root,elem,grid,z)
!                    ------------------------

!----------------------------------------------------
! This routine recursively goes down the refinement tree, starting at root,
! to find a element containing (x,y) that is marked as oldleaf
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: root
integer, intent(out) :: elem
type(grid_type), intent(in) :: grid
real(my_real), optional, intent(in) :: z

!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
real(my_real) :: bc(VERTICES_PER_ELEMENT,1)
!----------------------------------------------------
! Begin executable code

! if root is oldleaf, we've found it

if (grid%element(root)%oldleaf) then
   elem = root
   return
endif

! otherwise, check the children

allc = ALL_CHILDREN
children = get_child_lid(grid%element(root)%gid,allc,grid%elem_hash)
if (children(1) == NO_CHILD) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("evaluate_oldsoln could not find a containing element marked as oldleaf.")
   stop
endif

! look at the barycentric coordinates of (x,y) in each child,
! until one is found where they are all positive

elem = 0
do i=1,MAX_CHILD
      select case(global_element_kind)
      case (TRIANGULAR_ELEMENT)
         call barycentric((/x/),(/y/), &
                       grid%vertex(grid%element(children(i))%vertex)%coord%x, &
                       grid%vertex(grid%element(children(i))%vertex)%coord%y, &
                       bc,no_det=.true.)
      case (TETRAHEDRAL_ELEMENT)
         call barycentric_tet((/x/),(/y/),(/z/), &
                       grid%vertex(grid%element(children(i))%vertex)%coord%x, &
                       grid%vertex(grid%element(children(i))%vertex)%coord%y, &
                  zcoord(grid%vertex(grid%element(children(i))%vertex)%coord), &
                       bc)
      end select
   if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
      call find_old_containing_leaf(x,y,children(i),elem,grid,z)
      exit
   endif
end do

end subroutine find_old_containing_leaf

!          ---------------------
subroutine set_grid_for_old_soln(grid)
!          ---------------------

!----------------------------------------------------
! This routine sets grid_for_old_soln so that the grid is accessible even
! though it cannot be passed through the user routines pdecoef, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), target :: grid

!----------------------------------------------------
! Begin executable code

grid_for_old_soln => grid

end subroutine set_grid_for_old_soln

end module evaluate
