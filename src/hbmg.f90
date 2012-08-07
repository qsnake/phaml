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

module hbmg

!----------------------------------------------------
! This module contains routines for the hierarchical basis multigrid method
!
! communication tags in this module are of the form 14xx and 16xx
! TEMP090127 new communication approach
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use gridtype_mod
use linsystype_mod
use linsys_util
use lapack_solve
use linsys_io
use error_estimators
use hash_mod
use omp_lib, only: omp_get_num_threads

!----------------------------------------------------

implicit none
private
public multigrid, mg_precon

!----------------------------------------------------
! The following variables are defined:

! When running the phaml tests, you want IDENTICAL results to a sequential run.
logical, parameter :: global_need_identical = .false.
!----------------------------------------------------

contains

!          ---------
subroutine multigrid(grid,procs,linear_system,io_cntl,solver_cntl, &
                     still_sequential,maxlev,no_master)
!          ---------

!----------------------------------------------------
! This routine solves the linear system via multigrid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(linsys_type), intent(inout) :: linear_system
type(io_options), intent(in) :: io_cntl
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
integer, intent(in), optional :: maxlev
logical, intent(in), optional :: no_master

!----------------------------------------------------
! Local variables:

integer :: ncyc, cycle, lev, nu1, nu2, nu1ho, nu2ho, nlev, proc, ni, nr
real(my_real) :: mg_tol, resid
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
logical :: fudop_comm, conventional_comm, yes_master, use_simple_relax_ho
integer, allocatable :: neqlist(:,:,:), eqlist(:,:,:,:,:)

!----------------------------------------------------
! Begin executable code

if (global_element_kind == TETRAHEDRAL_ELEMENT) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("3D version of PHAML multigrid has not been written")
   stop
endif

if (present(no_master)) then
   yes_master = .not. no_master
else
   yes_master = .true.
endif

fudop_comm = solver_cntl%mg_comm == MGCOMM_FUDOP .and. &
             num_proc(procs) > 1 .and. &
             .not. still_sequential
conventional_comm = solver_cntl%mg_comm == MGCOMM_CONVENTIONAL .and. &
             num_proc(procs) > 1 .and. &
             .not. still_sequential
ncyc = solver_cntl%ncycle
nu1 = solver_cntl%prerelax
nu2 = solver_cntl%postrelax
nu1ho = solver_cntl%prerelax_ho
nu2ho = solver_cntl%postrelax_ho
if (present(maxlev)) then
   nlev = maxlev
else
   if (linear_system%maxdeg == 1) then
      nlev = linear_system%nlev_vert
   else
      nlev = linear_system%edge_level
   endif
endif

! determine tolerance of residual for termination

mg_tol = solver_cntl%mg_tol
if (mg_tol == MG_ERREST_TOL) then
   if (my_proc(procs) == MASTER) then
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1400)
      mg_tol = rrecv(1)
      deallocate(rrecv)
   else
! TEMP for systems of equations, should be a combination not the max
      call error_estimate(grid,procs,HIERARCHICAL_COEFFICIENT,errest_L2=mg_tol)
      mg_tol = sqrt(phaml_global_sum(procs,mg_tol**2,1401))
      mg_tol = max(100*epsilon(1.0_my_real),mg_tol/1000)
      if (my_proc(procs) == 1 .and. yes_master) then
         call phaml_send(procs,MASTER,(/0/),0,(/mg_tol/),1,1400)
      endif
   endif
endif

resid = 0.0_my_real

! master process just prints the error (if requested) and returns

if (my_proc(procs) == MASTER) then
   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      call linsys_residual(linear_system,procs,still_sequential,0,1410,.true., &
                           .true.)
   endif
   do cycle=1,ncyc
      if (io_cntl%print_error_when == FREQUENTLY .or. &
          io_cntl%print_error_when == TOO_MUCH) then
         call linsys_residual(linear_system,procs,still_sequential,cycle, &
                              1410+cycle,.true.,.false.,relresid=resid)
      else
         if (mg_tol /= MG_NO_TOL) then
            call linsys_residual(linear_system,procs,still_sequential,cycle, &
                                 1410+cycle,.false.,.false.,relresid=resid)
         endif
      endif
      if (resid < mg_tol) exit
   end do
   if (solver_cntl%solver == MG_SOLVER .and. mg_tol /= MG_NO_TOL .and. &
       resid > mg_tol) then
      call warning("multigrid did not reach tolerance; ncyc, tol, resid", &
                   intlist=(/ncyc/),reallist=(/mg_tol,resid/))
   endif
   linear_system%solver_niter = min(cycle,ncyc)
   if (io_cntl%print_error_when == FREQUENTLY .or. &
          io_cntl%print_error_when == TOO_MUCH .or. mg_tol /= MG_NO_TOL) then
      linear_system%relresid = resid
   else
      linear_system%relresid = -huge(0.0_my_real)
   endif
   return
endif

! For p multigrid, create independent sets of edge equations so that the
! relaxations can be performed in OpenMP parallel without conflicts and
! always in the same order.

use_simple_relax_ho = (num_proc(procs) == 1 .and. omp_get_num_threads() == 1)

if (linear_system%maxdeg > 1 .and. .not. use_simple_relax_ho) then
   call make_pmg_sets
endif

! make sure the solution at unowned points is current

if (fudop_comm .or. conventional_comm) then
   call exchange_fudop_vect(linear_system%solution(1:),procs, &
                            linear_system,1402,1403,1404)
endif

! print the error before solution (if requested)

if (io_cntl%print_error_when == FREQUENTLY .or. &
    io_cntl%print_error_when == TOO_MUCH) then
   call linsys_residual(linear_system,procs,still_sequential,0,1410,.true., &
                        .true.,no_master=no_master)
endif

!$omp parallel &
!$omp default(shared) &
!$omp private(cycle,lev)

! initialize r_other

if (fudop_comm .or. conventional_comm) then
   call init_r_other(linear_system%solution(1:),procs,still_sequential)
else
!$omp master
   linear_system%r_mine = 0.0_my_real
   linear_system%r_others = 0.0_my_real
!$omp end master
!$omp barrier
endif

! repeat ncycle times or until L2 norm of residual < mg_tol

do cycle=1,ncyc

! for each level from finest to coarsest, relaxation and restriction

   do lev=nlev,2,-1
      if (lev > linear_system%nlev_vert) then
         if (use_simple_relax_ho) then
            call singleproc_relax_ho(lev,nu1ho,linear_system,-1)
         else
            call relax_ho(lev,nu1ho,linear_system,-1,procs, &
                          fudop_comm.or.conventional_comm,still_sequential, &
                          neqlist,eqlist)
         endif
      else
         call relax(lev,nu1,linear_system,conventional_comm)
      endif
!$omp master
      call basis_change(lev,TO_HIER,linear_system)
!$omp end master
!$omp barrier
      call rhs_minus_Axh(linear_system,lev,still_sequential,procs)
      if (conventional_comm) then
!$omp master
         call exchange_fudop_soln_residual(procs,linear_system, &
                                           1420+cycle+lev,1430+cycle+lev, &
                                           1440+cycle+lev,max_lev=lev)
!$omp end master
!$omp barrier
      endif
   end do

! fix solution values and residuals via information from other processors

   if (fudop_comm) then
!$omp master
      call exchange_fudop_soln_residual(procs,linear_system, &
                                        1420+cycle,1430+cycle,1440+cycle)
!$omp end master
!$omp barrier
   endif

! solve on coarsest grid

!$omp master
   call coarse_solve(linear_system)
!$omp end master
!$omp barrier

! for each level from coarsest to finest, prolongation and relaxation

   do lev=2,nlev
      call rhs_plus_Axh(linear_system,lev,still_sequential,procs)
!$omp master
      if (conventional_comm) then
         call exchange_fudop_soln_residual(procs,linear_system, &
                                           1450+cycle+lev,1460+cycle+lev, &
                                           1470+cycle+lev,max_lev=lev)
      endif
      call basis_change(lev,TO_NODAL,linear_system)
!$omp end master
!$omp barrier
      if (lev > linear_system%nlev_vert) then
         if (use_simple_relax_ho) then
            call singleproc_relax_ho(lev,nu2ho,linear_system,1)
         else
            call relax_ho(lev,nu2ho,linear_system,1,procs, &
                          fudop_comm.or.conventional_comm,still_sequential, &
                          neqlist,eqlist)
         endif
      else
         call relax(lev,nu2,linear_system,conventional_comm)
      endif
   end do

! fix solution values via information from other processors

   if (fudop_comm .or. conventional_comm) then
!$omp master
      call exchange_fudop_vect(linear_system%solution(1:),procs, &
                               linear_system,1450+cycle,1460+cycle, &
                               1470+cycle)
!$omp end master
!$omp barrier
   endif

! print error after this cycle (if requested) and test for small residual

!$omp master
   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      call linsys_residual(linear_system,procs,still_sequential,cycle, &
                           1410+cycle,.true.,.false.,relresid=resid,no_master=no_master)
   else
      if (mg_tol /= MG_NO_TOL) then
         call linsys_residual(linear_system,procs,still_sequential,cycle, &
                              1410+cycle,.false.,.false.,relresid=resid,no_master=no_master)
      endif
   endif
!$omp end master
!$omp barrier

   if (resid < mg_tol) exit

end do ! next cycle

!$omp master
linear_system%solver_niter = min(cycle,ncyc)
!$omp end master
!$omp barrier

!$omp end parallel

if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH .or. mg_tol /= MG_NO_TOL) then
   linear_system%relresid = resid
else
   linear_system%relresid = -huge(0.0_my_real)
endif

if (allocated(eqlist)) deallocate(neqlist, eqlist)

return

contains

!          ------------
subroutine init_r_other(solution,procs,still_sequential)
!          ------------

!----------------------------------------------------
! This routine initializes r_other before the first V-cycle
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

real(my_real) :: hold_soln(size(solution))
integer :: lev
!----------------------------------------------------
! Begin executable code

! clear out residuals

!$omp single
linear_system%r_mine = 0.0_my_real
linear_system%r_others = 0.0_my_real
!$omp end single

! hold the solution so that this doesn't change it

!$omp master
hold_soln = solution
!$omp end master
!$omp barrier

! do the first half of a V cycle to generate the residual

do lev=nlev,2,-1
   if (lev > linear_system%nlev_vert) then
      if (use_simple_relax_ho) then
         call singleproc_relax_ho(lev,nu1ho,linear_system,-1)
      else
         call relax_ho(lev,nu1ho,linear_system,-1,procs, &
                       fudop_comm.or.conventional_comm,still_sequential, &
                       neqlist,eqlist)
      endif
   else
      call relax(lev,nu1,linear_system,conventional_comm)
   endif
!$omp single
   call basis_change(lev,TO_HIER,linear_system)
!$omp end single
   call rhs_minus_Axh(linear_system,lev,still_sequential,procs)
!$omp single
   if (conventional_comm) then
      call exchange_fudop_soln_residual(procs,linear_system,1600+lev,1610+lev, &
                                        1620+lev,max_lev=lev)
   endif
!$omp end single
end do

! exchange the residual

!$omp single
if (fudop_comm) then
   call exchange_fudop_soln_residual(procs,linear_system,1600,1610,1620)
endif
!$omp end single

! second half of V cycle to restore to nodal form, but don't need relaxation

do lev=2,nlev
   call rhs_plus_Axh(linear_system,lev,still_sequential,procs)
!$omp single
   if (conventional_comm) then
      call exchange_fudop_soln_residual(procs,linear_system, &
                                        1630+lev,1640+lev,1650+lev, &
                                        no_soln=.true.,max_lev=lev)
   endif
   call basis_change(lev,TO_NODAL,linear_system)
!$omp end single
end do

! restore the solution

!$omp master
solution = hold_soln
!$omp end master
!$omp barrier

end subroutine init_r_other

!          -------------
subroutine make_pmg_sets()
!          -------------

!----------------------------------------------------
! This routine creates independent sets of edge equations so that the
! p-multigrid relaxations can be performed in OpenMP parallel without
! conflicts and always in the same order.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

! Uses grid, linear_system, neqlist, and eqlist from the host routine.

!----------------------------------------------------
! Local variables:

integer :: ss, astat, lev, elem, side, edge, mate, degm1, sys, eq
logical :: compat_div
logical(small_logical) :: on_list(grid%biggest_edge)
!----------------------------------------------------
! Begin executable code

! For each level, pass through the triangles creating one set with first edges,
! one set with second edges, and one set with third edges iff the triangle is
! compatibly divisible.  Don't add an edge more than once.  The sets are
! further separated by degree.

ss = linear_system%system_size

allocate(neqlist(grid%nlev,linear_system%maxdeg-1,3),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_pmg_sets")
   stop
endif
neqlist = 0

! First pass, count how many equations will be in each set.

on_list = .false.

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         compat_div = .false.
         do side = 1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(side)
            if (.not. on_list(edge)) then
               if (side == 3) then
                  if (grid%element(elem)%mate == BOUNDARY) then
                     compat_div = .true.
                  else
                     mate = hash_decode_key(grid%element(elem)%mate, &
                                            grid%elem_hash)
                     compat_div = (mate /= HASH_NOT_FOUND)
                  endif
               endif
               if (side < 3 .or. compat_div) then
                  neqlist(lev,1:grid%edge(edge)%degree-1,side) = &
                            neqlist(lev,1:grid%edge(edge)%degree-1,side) + 1
                  on_list(edge) = .true.
               endif
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

! Allocate results

allocate(eqlist(maxval(neqlist),grid%nlev,linear_system%maxdeg-1,3,ss), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_pmg_sets")
   stop
endif

! Second pass, fill in the equations

neqlist = 0
on_list = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         compat_div = .false.
         do side = 1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(side)
            if (.not. on_list(edge)) then
               if (side == 3) then
                  if (grid%element(elem)%mate == BOUNDARY) then
                     compat_div = .true.
                  else
                     mate = hash_decode_key(grid%element(elem)%mate, &
                                            grid%elem_hash)
                     compat_div = (mate /= HASH_NOT_FOUND)
                  endif
               endif
               if (side < 3 .or. compat_div) then
                  do degm1=1,grid%edge(edge)%degree-1
                     neqlist(lev,degm1,side) = neqlist(lev,degm1,side) + 1
                     do sys=1,ss
                        call grid_to_eq(grid,linear_system,EDGE_ID,degm1,sys, &
                                        grid%edge(edge)%gid,eq)
                        eqlist(neqlist(lev,degm1,side),lev,degm1,side,sys) = eq
                     end do
                  end do
                  on_list(edge) = .true.
               endif
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

end subroutine make_pmg_sets

end subroutine multigrid

!          -----
subroutine relax(lev,nu,matrix,conventional_comm)
!          -----

!----------------------------------------------------
! This routine perform nu/2 iterations of red-(local)black relaxation, where
! the red points are the level lev points.
! Half an iteration means red relaxation.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: lev, nu
type(linsys_type), intent(inout) :: matrix
logical, intent(in) :: conventional_comm

!----------------------------------------------------
! Local variables:

integer :: iter, eq, i, j, neigh, allocstat, ss, start_block, in_block
integer, save :: nblack(4)
logical :: red_only
integer, allocatable, save :: black_list(:,:)
logical(small_logical), allocatable, save :: on_black_list(:)
!----------------------------------------------------
! Begin executable code

ss = matrix%system_size

! Create three independent sets of black neighbors of red equations.  The
! first is those with level lev-1, the second lev-2, and third has levels
! below lev-2.  Each has no neighbors within the set, allowing them to be
! processed in parallel with no concerns.

!$omp single
if (nu > 1) then
   allocate(on_black_list(matrix%neq),black_list(matrix%neq,4),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed for black_list in relax")
      stop
   endif
   nblack = 0
   on_black_list = .false.
   do eq = matrix%begin_level(lev), matrix%begin_level(lev+1)-1
      if (matrix%equation_type(eq) == DIRICHLET) cycle
      do i=matrix%begin_row(eq)+1,matrix%end_row(eq)
         neigh = matrix%column_index(i)
         if (neigh == NO_ENTRY) cycle
         if (on_black_list(neigh)) cycle
         if (neigh >= matrix%begin_level(lev)) cycle ! must be a lower level
         if (neigh < matrix%begin_level(2)) then
            on_black_list(neigh) = .true.
            nblack(4) = nblack(4) + 1
            black_list(nblack(4),4) = neigh
         elseif (neigh >= matrix%begin_level(lev-1)) then
            on_black_list(neigh) = .true.
            nblack(1) = nblack(1) + 1
            black_list(nblack(1),1) = neigh
         elseif (neigh >= matrix%begin_level(lev-2)) then
            on_black_list(neigh) = .true.
            nblack(2) = nblack(2) + 1
            black_list(nblack(2),2) = neigh
         else
            on_black_list(neigh) = .true.
            nblack(3) = nblack(3) + 1
            black_list(nblack(3),3) = neigh
         endif
      end do
   end do
endif
!$omp end single

! repeat nu/2 times, but make sure the red get relaxed when nu is odd

do iter = 1, (nu+1)/2

! will this be red relaxation only?

!                nu is odd           last iteration
   red_only = (2*(nu/2) /= nu .and. iter == (nu+1)/2)

! for each non-Dirichlet equation of level lev

!$omp do
   do start_block = matrix%begin_level(lev), matrix%begin_level(lev+1)-1, ss
      do in_block = 0, ss-1
         eq = start_block + in_block

         if (matrix%equation_type(eq) == DIRICHLET) cycle

! red relaxation

         if (.not. conventional_comm .or. matrix%iown(eq)) then
            matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                                  matrix%r_others(eq)
            do i=matrix%begin_row(eq)+1,matrix%end_row(eq)
               matrix%solution(eq) = matrix%solution(eq) - &
                   matrix%matrix_val(i)*matrix%solution(matrix%column_index(i))
            end do
            matrix%solution(eq) = matrix%solution(eq) / &
                                  matrix%matrix_val(matrix%begin_row(eq))
         endif

      end do
   end do
!$omp end do

! black relaxation, if desired

   if (red_only) exit

   do j=1,3
!$omp do
      do start_block = 1, nblack(j), ss
         do in_block = 0, ss-1
            eq = black_list(start_block+in_block,j)
            if (matrix%equation_type(eq) == DIRICHLET) cycle
            if (.not. conventional_comm .or. matrix%iown(eq)) then
               matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                                     matrix%r_others(eq)
               do i=matrix%begin_row(eq)+1,matrix%end_row(eq)
                  matrix%solution(eq) = matrix%solution(eq) - &
                   matrix%matrix_val(i)*matrix%solution(matrix%column_index(i))
               end do
               matrix%solution(eq) = &
                  matrix%solution(eq) / matrix%matrix_val(matrix%begin_row(eq))
            endif
         end do
      end do
!$omp end do

   end do
! do level 1 equations sequentially because they can be neighbors of each other
!$omp single
   do start_block = 1, nblack(4), ss
      do in_block = 0, ss-1
         eq = black_list(start_block+in_block,4)
         if (matrix%equation_type(eq) == DIRICHLET) cycle
         matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                               matrix%r_others(eq)
         do i=matrix%begin_row(eq)+1,matrix%end_row(eq)
            matrix%solution(eq) = matrix%solution(eq) - &
                matrix%matrix_val(i)*matrix%solution(matrix%column_index(i))
         end do
         matrix%solution(eq) = &
               matrix%solution(eq) / matrix%matrix_val(matrix%begin_row(eq))
      end do
   end do
!$omp end single

end do ! next iteration

! free memory

!$omp barrier
!$omp single
if (nu > 1) then
   if (allocated(black_list)) then
      deallocate(on_black_list,black_list,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed in relax")
      endif
   endif
endif
!$omp end single

end subroutine relax

!          --------
subroutine relax_ho(edgelev,nu,matrix,dir,procs,communicate,still_sequential, &
                    neqlist,eqlist)
!          --------

!----------------------------------------------------
! This routine performs one direction of the p-multigrid cycle.  If dir==-1
! it performs nu Gauss-Seidel iterations with the entire (condensed) matrix,
! then with the matrix up to maxdeg-1, then maxdeg-2, ... 2.  If dir==1 it
! starts at 2 and works up to maxdeg.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: edgelev, nu
type(linsys_type), intent(inout) :: matrix
integer, intent(in) :: dir
type(proc_info), intent(in) :: procs
logical, intent(in) :: communicate
logical, intent(in) :: still_sequential
integer, intent(in) :: neqlist(:,:,:), eqlist(:,:,:,:,:)

!----------------------------------------------------
! Local variables:

integer :: eq, i, j, gsiter, loop, astat, plev, degm1, lev, set, ss, sys
logical, save, allocatable :: isneigh(:)
real(my_real) :: temp
!----------------------------------------------------
! Begin executable code
 
ss = matrix%system_size

!$omp single
allocate(isneigh(matrix%begin_level(edgelev+1)-1),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in relax_ho")
   stop
endif
!$omp end single

! Downward half of V-cycle

if (dir == -1) then

! For each p level

   do plev = matrix%maxdeg,2,-1

! isneigh is used to limit the equations with level less than plev to those
! that are neighbors of equations of level plev

!$omp single
      isneigh = .false.
!$omp end single

! nu Gauss-Siedel iterations

      do gsiter=1,nu

! With new communication, first loop does interior equations and second
! loop does equations by the partition boundary

         do loop=1,2
            if (.not. new_comm .and. loop==2) cycle ! TEMP090127

! Communicate solution values near the partition boundary

            if (communicate) then
!$omp single
               if (new_comm) then ! TEMP090127
                  if (loop==1) then
                     call send_neigh_vect(matrix%solution(1:),procs,matrix, &
                                          1660+plev,min(plev+1,matrix%maxdeg))
                  else
                     call recv_neigh_vect(matrix%solution(1:),procs,matrix, &
                                          1660+plev)
                  endif
               else ! TEMP090127
                  call exchange_neigh_vect(matrix%solution(1:),procs, & ! TEMP090127
                                           matrix,1660+plev,1670+plev,1680+plev) ! TEMP090127
               endif ! TEMP090127
!$omp end single
            endif

! Perform relaxation on all equations (with degree plev or neighbors of degree
! plev) up to degree plev, in independent sets.  On the outer most loop, do
! degree plev down to degree 2, followed by a separate loop for degree 1 because
! it has to be done sequentially.  Inside of that do each h-level from finest
! to coarsest, and inside that are the sets associated with first, second
! and third sides.

            do degm1=plev-1,1,-1
               do lev=size(neqlist,dim=1),1,-1
                  do set=1,3
!$omp do
                     do i = 1,neqlist(lev,degm1,set)
                        do sys=1,ss
                           eq = eqlist(i,lev,degm1,set,sys)
                           if (new_comm) then ! TEMP090127
                              if ((loop==1 .and. matrix%nn_comm_remote_neigh(eq)) .or. &
                                  (loop==2 .and. .not. matrix%nn_comm_remote_neigh(eq))) then
                                 cycle
                              endif
                           endif ! TEMP090127
                           if (matrix%equation_type(eq) == DIRICHLET) cycle
                           if (degm1 < plev-1 .and. .not. isneigh(eq)) cycle
                           if (.not. communicate .or. matrix%iown(eq)) then
                              temp = matrix%rhs(eq) + matrix%r_mine(eq) + matrix%r_others(eq)
                              do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
                                 if (matrix%column_index(j) == NO_ENTRY) cycle
                                 if (degm1 == plev-1) then
                                    isneigh(matrix%column_index(j)) = .true.
                                 endif
                                 temp = temp - &
                                        matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
                              end do
                              matrix%solution(eq) = temp / &
                                         matrix%matrix_val(matrix%begin_row(eq))
                           endif
                        end do
                     end do
!$omp end do
                  end do ! set
               end do ! h-level
            end do ! degree
!$omp single
            do eq=1,matrix%begin_level(matrix%nlev_vert+1)-1
               if (new_comm) then ! TEMP090127
                  if ((loop==1 .and. matrix%nn_comm_remote_neigh(eq)) .or. &
                      (loop==2 .and. .not. matrix%nn_comm_remote_neigh(eq))) then
                     cycle
                  endif
               endif ! TEMP090127
               if (matrix%equation_type(eq) == DIRICHLET) cycle
               if (.not. isneigh(eq)) cycle
               if (.not. communicate .or. matrix%iown(eq)) then
                  matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                                        matrix%r_others(eq)
                  do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
                     if (matrix%column_index(j) == NO_ENTRY) cycle
                     matrix%solution(eq) = matrix%solution(eq) - &
                            matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
                  end do
                  matrix%solution(eq) = matrix%solution(eq) / &
                                        matrix%matrix_val(matrix%begin_row(eq))
               endif
            end do
!$omp end single
         end do ! loop for communication
      end do ! next Gauss-Seidel iteration
   end do ! p-level

   if (communicate) then
!$omp single
      call exchange_fudop_vect(matrix%solution(1:),procs,matrix,1661,1671,1681)
!$omp end single
   endif

! Upward half of V-cycle
! See the downward half for comments.  This is the same except reversing plev.

elseif (dir == 1) then

   do plev = 2,matrix%maxdeg
!$omp single
      isneigh = .false.
!$omp end single
      do gsiter=1,nu
         do loop=1,2
            if (.not. new_comm .and. loop==2) cycle ! TEMP090127
            if (communicate) then
!$omp single
               if (new_comm) then ! TEMP090127
                  if (loop==1) then
                     call send_neigh_vect(matrix%solution(1:),procs,matrix, &
                                          1660+plev,min(plev+1,matrix%maxdeg))
                  else
                     call recv_neigh_vect(matrix%solution(1:),procs,matrix, &
                                          1660+plev)
                  endif
               else ! TEMP090127
                  call exchange_neigh_vect(matrix%solution(1:),procs, & ! TEMP090127
                                           matrix,1480+plev,1490+plev,1690+plev) ! TEMP090127
               endif ! TEMP090127
!$omp end single
            endif
            do degm1=plev-1,1,-1
               do lev=size(neqlist,dim=1),1,-1
                  do set=1,3
!$omp do
                     do i = 1,neqlist(lev,degm1,set)
                        do sys=1,ss
                           eq = eqlist(i,lev,degm1,set,sys)
                           if (new_comm) then ! TEMP090127
                              if ((loop==1 .and. matrix%nn_comm_remote_neigh(eq)) .or. &
                                  (loop==2 .and. .not. matrix%nn_comm_remote_neigh(eq))) then
                                 cycle
                              endif
                           endif ! TEMP090127
                           if (matrix%equation_type(eq) == DIRICHLET) cycle
                           if (degm1 < plev-1 .and. .not. isneigh(eq)) cycle
                           if (.not. communicate .or. matrix%iown(eq)) then
                              temp = matrix%rhs(eq) + matrix%r_mine(eq) + matrix%r_others(eq)
                              do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
                                 if (matrix%column_index(j) == NO_ENTRY) cycle
                                 if (degm1 == plev-1) then
                                    isneigh(matrix%column_index(j)) = .true.
                                 endif
                                 temp = temp - &
                                        matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
                              end do
                              matrix%solution(eq) = temp / &
                                         matrix%matrix_val(matrix%begin_row(eq))
                           endif
                        end do
                     end do
!$omp end do
                  end do ! set
               end do ! h-level
            end do ! degree
!$omp single
            do eq=1,matrix%begin_level(matrix%nlev_vert+1)-1
               if (new_comm) then ! TEMP090127
                  if ((loop==1 .and. matrix%nn_comm_remote_neigh(eq)) .or. &
                      (loop==2 .and. .not. matrix%nn_comm_remote_neigh(eq))) then
                     cycle
                  endif
               endif ! TEMP090127
               if (matrix%equation_type(eq) == DIRICHLET) cycle
               if (.not. isneigh(eq)) cycle
               if (.not. communicate .or. matrix%iown(eq)) then
                  matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                                        matrix%r_others(eq)
                  do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
                     if (matrix%column_index(j) == NO_ENTRY) cycle
                     matrix%solution(eq) = matrix%solution(eq) - &
                            matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
                  end do
                  matrix%solution(eq) = matrix%solution(eq) / &
                                        matrix%matrix_val(matrix%begin_row(eq))
               endif
            end do
!$omp end single
         end do ! loop for communication
      end do ! next Gauss-Seidel iteration
   end do ! p-level

else

   call fatal("bad value for dir in relax_ho",intlist=(/dir/))
   stop

endif

!$omp barrier
!$omp single
deallocate(isneigh)
!$omp end single

end subroutine relax_ho

!          -------------------
subroutine singleproc_relax_ho(edgelev,nu,matrix,dir)
!          -------------------

!----------------------------------------------------
! This routine performs one direction of the p-multigrid cycle when we are
! running one thread on one processor.    If dir==-1
! it performs nu Gauss-Seidel iterations with the entire (condensed) matrix,
! then with the matrix up to maxdeg-1, then maxdeg-2, ... 2.  If dir==1 it
! starts at 2 and works up to maxdeg.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: edgelev, nu
type(linsys_type), intent(inout) :: matrix
integer, intent(in) :: dir
!----------------------------------------------------
! Local variables:

logical :: isneigh(matrix%begin_level(edgelev+1)-1)
integer :: plev, gsiter, eq, j
real(my_real) :: temp
!----------------------------------------------------
! Begin executable code

! Downward half of V-cycle

if (dir == -1) then

! For each p level

   do plev = matrix%maxdeg,2,-1

! isneigh is used to limit the equations with level less than plev to those
! that are neighbors of equations of level plev

      isneigh = .false.

! nu Gauss-Siedel iterations

      do gsiter=1,nu

! Perform relaxation on equations with degree plev or neighbors of degree plev

! Do level plev first to identify neighbors

         do eq=matrix%edge_begin_degree(plev),matrix%edge_begin_degree(plev+1)-1
            if (matrix%equation_type(eq) == DIRICHLET) cycle
            temp = matrix%rhs(eq) + matrix%r_mine(eq) + matrix%r_others(eq)
            do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
               if (matrix%column_index(j) == NO_ENTRY) cycle
               isneigh(matrix%column_index(j)) = .true.
               temp = temp - &
                    matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
            end do
            matrix%solution(eq) = temp / matrix%matrix_val(matrix%begin_row(eq))
         end do

! Then do level 1 up to level plev-1

         do eq=1,matrix%edge_begin_degree(plev)-1
            if (matrix%equation_type(eq) == DIRICHLET) cycle
            if (.not. isneigh(eq)) cycle
            temp = matrix%rhs(eq) + matrix%r_mine(eq) + matrix%r_others(eq)
            do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
               if (matrix%column_index(j) == NO_ENTRY) cycle
               temp = temp - &
                    matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
            end do
            matrix%solution(eq) = temp / matrix%matrix_val(matrix%begin_row(eq))
         end do

      end do
   end do
 
! Upward half of V-cycle
! See the downward half for comments.  This is the same except reversing plev.

elseif (dir == 1) then
 
   do plev = 2,matrix%maxdeg
      isneigh = .false.
      do gsiter=1,nu
         do eq=matrix%edge_begin_degree(plev),matrix%edge_begin_degree(plev+1)-1
            if (matrix%equation_type(eq) == DIRICHLET) cycle
            temp = matrix%rhs(eq) + matrix%r_mine(eq) + matrix%r_others(eq)
            do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
               if (matrix%column_index(j) == NO_ENTRY) cycle
               isneigh(matrix%column_index(j)) = .true.
               temp = temp - &
                    matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
            end do
            matrix%solution(eq) = temp / matrix%matrix_val(matrix%begin_row(eq))
         end do
         do eq=1,matrix%edge_begin_degree(plev)-1
            if (matrix%equation_type(eq) == DIRICHLET) cycle
            if (.not. isneigh(eq)) cycle
            temp = matrix%rhs(eq) + matrix%r_mine(eq) + matrix%r_others(eq)
            do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
               if (matrix%column_index(j) == NO_ENTRY) cycle
               temp = temp - &
                    matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
            end do
            matrix%solution(eq) = temp / matrix%matrix_val(matrix%begin_row(eq))
         end do
      end do
   end do

endif

end subroutine singleproc_relax_ho

!          ------------
subroutine rhs_plus_Axh(matrix,lev,still_sequential,procs)
!          ------------

!----------------------------------------------------
! This routine computes the addition of the coarse-fine block of the
! hierarchical matrix times the fine part of the hierarchical solution
! vector and removes it from r_mine and r_others
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: matrix
integer, intent(in) :: lev
logical, intent(in) :: still_sequential
type(proc_info), intent(in) :: procs

!----------------------------------------------------
! Local variables:

integer :: eq, col
logical :: need_identical

real(my_real) :: temp
!----------------------------------------------------
! Begin executable code

! The results must be identical to that of the other processors if there
! is more than one processor and it is still sequential

if (still_sequential .and. num_proc(procs) > 1) then
   need_identical = .true.
else
   need_identical = global_need_identical
endif

! To get identical results, the floating point additions that create r_mine
! and r_others need to be in the same order, which sequentiallizes the loop.

if (need_identical) then

!$omp single
   do eq = matrix%begin_level(lev),matrix%begin_level(lev+1)-1
      do col = matrix%begin_row(eq)+1,matrix%end_row(eq)
         if (matrix%column_index(col) == NO_ENTRY) cycle
         if (matrix%column_index(col) >= matrix%begin_level(lev)) cycle
         if (matrix%iown(eq)) then
            matrix%r_mine(matrix%column_index(col)) = &
                                     matrix%r_mine(matrix%column_index(col)) + &
                                     matrix%matrix_val(col)*matrix%solution(eq)
         else
            matrix%r_others(matrix%column_index(col)) = &
                                   matrix%r_others(matrix%column_index(col)) + &
                                   matrix%matrix_val(col)*matrix%solution(eq)
         endif
      end do
   end do
!$omp end single

else

!$omp do
   do eq = matrix%begin_level(lev),matrix%begin_level(lev+1)-1
      do col = matrix%begin_row(eq)+1,matrix%end_row(eq)
         if (matrix%column_index(col) == NO_ENTRY) cycle
         if (matrix%column_index(col) >= matrix%begin_level(lev)) cycle
         temp = matrix%matrix_val(col)*matrix%solution(eq)
         if (matrix%iown(eq)) then
!$omp atomic
            matrix%r_mine(matrix%column_index(col)) = &
                                  matrix%r_mine(matrix%column_index(col)) + temp
         else
!$omp atomic
            matrix%r_others(matrix%column_index(col)) = &
                                matrix%r_others(matrix%column_index(col)) + temp
         endif
      end do
   end do
!$omp end do

endif

end subroutine rhs_plus_Axh

!          -------------
subroutine rhs_minus_Axh(matrix,lev,still_sequential,procs)
!          -------------

!----------------------------------------------------
! This routine computes the subtraction of the coarse-fine block of the
! hierarchical matrix times the fine part of the hierarchical solution
! vector and stores it in r_mine and r_others
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: matrix
integer, intent(in) :: lev
logical, intent(in) :: still_sequential
type(proc_info), intent(in) :: procs

!----------------------------------------------------
! Local variables:

integer :: eq, col
logical :: need_identical

real(my_real) :: temp
!----------------------------------------------------
! Begin executable code

! The results must be identical to that of the other processors if there
! is more than one processor and it is still sequential

if (still_sequential .and. num_proc(procs) > 1) then
   need_identical = .true.
else
   need_identical = global_need_identical
endif

! To get identical results, the floating point additions that create r_mine
! and r_others need to be in the same order, which sequentiallizes the loop.

if (need_identical) then

!$omp single
   do eq = matrix%begin_level(lev),matrix%begin_level(lev+1)-1
      do col = matrix%begin_row(eq)+1,matrix%end_row(eq)
         if (matrix%column_index(col) == NO_ENTRY) cycle
         if (matrix%column_index(col) >= matrix%begin_level(lev)) cycle
         if (matrix%iown(eq)) then
            matrix%r_mine(matrix%column_index(col)) = &
                                     matrix%r_mine(matrix%column_index(col)) - &
                                     matrix%matrix_val(col)*matrix%solution(eq)
         else
            matrix%r_others(matrix%column_index(col)) = &
                                   matrix%r_others(matrix%column_index(col)) - &
                                   matrix%matrix_val(col)*matrix%solution(eq)
         endif
      end do
   end do
!$omp end single

else

!$omp do
   do eq = matrix%begin_level(lev),matrix%begin_level(lev+1)-1
      do col = matrix%begin_row(eq)+1,matrix%end_row(eq)
         if (matrix%column_index(col) == NO_ENTRY) cycle
         if (matrix%column_index(col) >= matrix%begin_level(lev)) cycle
         temp = matrix%matrix_val(col)*matrix%solution(eq)
         if (matrix%iown(eq)) then
!$omp atomic
            matrix%r_mine(matrix%column_index(col)) = &
                                  matrix%r_mine(matrix%column_index(col)) - temp
         else
!$omp atomic
            matrix%r_others(matrix%column_index(col)) = &
                                matrix%r_others(matrix%column_index(col)) - temp
         endif
      end do
   end do
!$omp end do

endif

end subroutine rhs_minus_Axh

!          ------------
subroutine coarse_solve(linear_system)
!          ------------

!----------------------------------------------------
! This routine solves the coarse grid problem using the LAPACK routines for
! symmetric positive definite band matricies.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! All partitions have the same coarse grid.  Let each processor do the
! coarse grid for itself.

! If the coarse grid matrix has not been put into band form and
! factorized, do so.

if (.not. linear_system%coarse_band_exists) then
   call make_lapack_symm_band(1,linear_system,linear_system%coarse_matrix)
   if (ierr /= NO_ERROR) return
   linear_system%coarse_band_exists = .true.
endif

! Solve the coarse grid system

call lapack_spd(1,linear_system,linear_system%coarse_matrix)

end subroutine coarse_solve

!          ---------
subroutine mg_precon(invec,outvec,matrix,grid,procs,io_cntl, &
                     solver_cntl,still_sequential)
!          ---------

!----------------------------------------------------
! This routine applies a multigrid V cycle as a preconditioner.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: invec(:)
real(my_real), intent(out) :: outvec(:)
type(linsys_type), intent(inout) :: matrix
type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(io_options), intent(in) :: io_cntl
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

real(my_real) :: holdrhs(size(matrix%rhs)), holdsoln(0:size(matrix%rhs))
!----------------------------------------------------
! Begin executable code

! Keep rhs and solution

holdrhs = matrix%rhs
holdsoln = matrix%solution

! Copy the invec to rhs; the size of invec should be the same as rhs

matrix%rhs = invec

! Set the initial guess to 0.0

matrix%solution(1:) = 0.0_my_real

! Set Dirichlet points

where (matrix%equation_type == DIRICHLET) matrix%solution(1:) = matrix%rhs

! Apply a multigrid cycle

call multigrid(grid,procs,matrix,io_cntl,solver_cntl,still_sequential)

! Copy solution (which now contains the preconditioner times invec) to outvec

outvec = matrix%solution(1:)

! Restore rhs and solution

matrix%rhs = holdrhs
matrix%solution = holdsoln

end subroutine mg_precon

end module hbmg
