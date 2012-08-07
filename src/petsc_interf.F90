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

module petsc_interf

!----------------------------------------------------
! This module contains routines that interface to PETSc.  These are
! in a separate file, rather than being part of module linear_system,
! because PETSc requires the use of a C-type preprocessor (by using
! F90 instead of f90 as the suffix, with most compilers).
!
! communication tags in this module are of the form 16xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use hash_mod
use hash_eq_mod
use message_passing
use petsc_type_mod
use gridtype_mod
use linsystype_mod
use hbmg
use lapack_solve
use linsys_util

!----------------------------------------------------

implicit none
private
public create_petsc_linear_system, create_petsc_linear_system_mf, &
       change_petsc_rhs, petsc_solve, destroy_petsc_linear_system, &
       phamlvec_to_petsc, petscvec_to_phaml, set_petsc_hold, &
       set_petsc_solver, set_petsc_preconditioner, set_petsc_options

!----------------------------------------------------
! The PETSc include files.  Note the use of preprocessor #include instead of
! the Fortran include statement, because the include files contain
! preprocessor directives.

! At PETSc 3.1 petsc.h was changed to just include the other .h's, so we
! no longer need the others (in fact, they causes duplicate declarations)

#include "include/petscversion.h"
#include "include/finclude/petsc.h"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0))
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#endif

!----------------------------------------------------
! The following parameters are defined:

!----------------------------------------------------

!----------------------------------------------------
! The following types are defined:

type petsc_hold_data
   type(petsc_matrix_type), pointer :: petsc_matrix
   type(linsys_type), pointer :: phaml_matrix
   type(linsys_type), pointer :: phaml_full_matrix
   type(grid_type), pointer :: grid
   type(proc_info), pointer :: procs
   type(io_options), pointer :: io_cntl
   type(solver_options), pointer :: solver_cntl
   logical :: still_sequential
end type petsc_hold_data

!----------------------------------------------------

!----------------------------------------------------
! The following variables are defined:

! Pointers to maintain access to variables from callbacks

type(petsc_hold_data) :: petsc_hold

integer :: null_data ! no data to pass

!----------------------------------------------------

contains

!          --------------
subroutine set_petsc_hold(petsc_matrix,phaml_matrix,grid,procs,io_cntl, &
                          solver_cntl,still_sequential)
!          --------------

!----------------------------------------------------
! This routine sets selected components of petsc_hold
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(petsc_matrix_type), target, optional :: petsc_matrix
type(linsys_type), target, optional :: phaml_matrix
type(grid_type), target, optional :: grid
type(proc_info), target, optional :: procs
type(io_options), target, optional :: io_cntl
type(solver_options), target, optional :: solver_cntl
logical, optional :: still_sequential

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (present(petsc_matrix)) petsc_hold%petsc_matrix => petsc_matrix
if (present(phaml_matrix)) petsc_hold%phaml_matrix => phaml_matrix
if (present(grid)) petsc_hold%grid => grid
if (present(procs)) petsc_hold%procs => procs
if (present(io_cntl)) petsc_hold%io_cntl => io_cntl
if (present(solver_cntl)) petsc_hold%solver_cntl => solver_cntl
if (present(still_sequential)) petsc_hold%still_sequential = still_sequential

end subroutine set_petsc_hold

!          --------------------------
subroutine create_petsc_linear_system(phaml_matrix,petsc_matrix, &
                                      solver_cntl,still_sequential,procs)
!          --------------------------

!----------------------------------------------------
! This routine creates a PETSc linear system (matrix A and vector b) from
! a PHAML linear system.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), target :: phaml_matrix
type(petsc_matrix_type), intent(out) :: petsc_matrix
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: astat, i, j, k, pi
PetscErrorCode jerr
PetscInt, allocatable :: d_nnz(:), o_nnz(:), cvalues(:)
PetscScalar, allocatable :: avalues(:), bvalues(:)
!----------------------------------------------------
! Begin executable code

petsc_matrix%my_total_eq = phaml_matrix%neq

! Set the equation hash table

petsc_matrix%eq_hash => phaml_matrix%eq_hash

! Determine which equations I own and which are Dirichlet rows


allocate(petsc_matrix%iown(petsc_matrix%my_total_eq), &
         petsc_matrix%dirich(petsc_matrix%my_total_eq), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

petsc_matrix%iown = phaml_matrix%iown(1:petsc_matrix%my_total_eq)
! TEMP some preconditioners don't have the Dirichlet rows eliminated
if (solver_cntl%preconditioner == MG_PRECONDITION .or. &
    solver_cntl%preconditioner == FUDOP_DD_PRECONDITION .or. &
    solver_cntl%preconditioner == COARSE_GRID_PRECONDITION) then
   petsc_matrix%dirich = .false.
else
   petsc_matrix%dirich = &
      phaml_matrix%equation_type(1:petsc_matrix%my_total_eq) == DIRICHLET
endif
petsc_matrix%my_own_eq = count(petsc_matrix%iown .and. .not.petsc_matrix%dirich)

! Determine how many nonzeroes in each row.

allocate(d_nnz(petsc_matrix%my_own_eq), o_nnz(petsc_matrix%my_own_eq), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

j = 0
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i) .and. .not.petsc_matrix%dirich(i)) then
      j = j+1
      d_nnz(j) = 0
      o_nnz(j) = 0
      do k=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (petsc_matrix%dirich(phaml_matrix%column_index(k))) cycle
         if (phaml_matrix%matrix_val(k) == 0.0_my_real) cycle
         if (petsc_matrix%iown(phaml_matrix%column_index(k))) then
            d_nnz(j) = d_nnz(j) + 1
         else
            o_nnz(j) = o_nnz(j) + 1
         endif
      end do
   endif
end do

! Create the matrix and rhs data structures

if (still_sequential) then
   call VecCreateSeq(PETSC_COMM_SELF,petsc_matrix%my_own_eq,petsc_matrix%b, &
                     jerr);CHKERRQ(jerr)
else
   call VecCreateMPI(PETSC_COMM_WORLD,petsc_matrix%my_own_eq,PETSC_DETERMINE, &
                     petsc_matrix%b,jerr);CHKERRQ(jerr)
endif

! Get the owned range for this processor and total number of equations.

call VecGetOwnershipRange(petsc_matrix%b,petsc_matrix%my_global_low, &
                          petsc_matrix%my_global_hi,jerr);CHKERRQ(jerr)
call VecGetSize(petsc_matrix%b,petsc_matrix%global_eq,jerr);CHKERRQ(jerr)

! Create the matrix data structure

if (still_sequential) then
   call MatCreateSeqAIJ(PETSC_COMM_SELF,petsc_matrix%my_own_eq, &
                        petsc_matrix%my_own_eq,0,d_nnz,petsc_matrix%A, &
                        jerr);CHKERRQ(jerr)
else
   call MatCreateMPIAIJ(PETSC_COMM_WORLD,petsc_matrix%my_own_eq, &
                        petsc_matrix%my_own_eq,PETSC_DETERMINE,PETSC_DETERMINE,&
                        0,d_nnz,0,o_nnz,petsc_matrix%A,jerr);CHKERRQ(jerr)
endif

! Determine the PETSc index for each equation on this processor by
! setting the ones that this processor owns and requesting the
! others from other processors.

allocate(petsc_matrix%petsc_index(petsc_matrix%my_total_eq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

! Set the ones this processor owns

petsc_matrix%petsc_index = 0
pi = petsc_matrix%my_global_low
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i) .and. .not.petsc_matrix%dirich(i)) then
      petsc_matrix%petsc_index(i) = pi
      pi = pi + 1
   endif
end do

if (pi /= petsc_matrix%my_global_hi) then
   call warning("count of owned equations not equal to size of petsc range", &
                intlist=(/pi,petsc_matrix%my_global_hi/))
endif

! Get the unowned ones from other processors

call update_shadows(petsc_matrix,phaml_matrix%gid,procs,1601,1602,1603, &
                    idata=petsc_matrix%petsc_index)

! copy the values from rhs to bvalues and remove Dirichlet boundary conditions

allocate(bvalues(petsc_matrix%my_own_eq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

j = 0
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i) .and. .not.petsc_matrix%dirich(i)) then
      j = j+1
      bvalues(j) = phaml_matrix%rhs(i)
      do k=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (petsc_matrix%dirich(phaml_matrix%column_index(k))) then
            bvalues(j) = bvalues(j) - phaml_matrix%matrix_val(k) * &
                                  phaml_matrix%rhs(phaml_matrix%column_index(k))
         endif
      end do
   endif
end do

! put bvalues in b

call VecSetValues(petsc_matrix%b,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  bvalues,INSERT_VALUES,jerr);CHKERRQ(jerr)
call VecAssemblyBegin(petsc_matrix%b,jerr);CHKERRQ(jerr)

! Copy the values in matrix_val into A

allocate(avalues(maxval(d_nnz)+maxval(o_nnz)), &
         cvalues(maxval(d_nnz)+maxval(o_nnz)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

do i=1,petsc_matrix%my_total_eq
   if (.not. petsc_matrix%iown(i) .or. petsc_matrix%dirich(i)) cycle
   j = 0
! TEMP for preconditioners that don't have Dirichlet rows eliminated
   if (phaml_matrix%equation_type(i) == DIRICHLET) then
      j = 1
      avalues(j) = 1
      cvalues(j) = petsc_matrix%petsc_index(i)
   else
      do k=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (petsc_matrix%dirich(phaml_matrix%column_index(k))) cycle
         if (phaml_matrix%matrix_val(k) == 0.0_my_real) cycle
         j = j+1
         avalues(j) = phaml_matrix%matrix_val(k)
         cvalues(j) = petsc_matrix%petsc_index(phaml_matrix%column_index(k))
      end do
   endif
   call MatSetValues(petsc_matrix%A,1,(/petsc_matrix%petsc_index(i)/), &
                     j,cvalues,avalues,INSERT_VALUES,jerr);CHKERRQ(jerr)
end do
call MatAssemblyBegin(petsc_matrix%A,MAT_FINAL_ASSEMBLY,jerr);CHKERRQ(jerr)

! finish messages associated with assembling b and A

call VecAssemblyEnd(petsc_matrix%b,jerr);CHKERRQ(jerr)
call MatAssemblyEnd(petsc_matrix%A,MAT_FINAL_ASSEMBLY,jerr);CHKERRQ(jerr)

deallocate(d_nnz,o_nnz,avalues,bvalues,cvalues)

end subroutine create_petsc_linear_system

!          -----------------------------
subroutine create_petsc_linear_system_mf(phaml_matrix,petsc_matrix, &
                                         still_sequential)
!          -----------------------------

!----------------------------------------------------
! This routine creates a PETSc linear system shell for matrix free solvers.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), target :: phaml_matrix
type(petsc_matrix_type), intent(out) :: petsc_matrix
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: astat,i,jerr
!----------------------------------------------------
! Begin executable code

petsc_matrix%my_total_eq = phaml_matrix%neq

! Set the equation hash table

petsc_matrix%eq_hash => phaml_matrix%eq_hash

! Determine which equations I own

allocate(petsc_matrix%iown(petsc_matrix%my_total_eq), &
         petsc_matrix%dirich(petsc_matrix%my_total_eq), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

petsc_matrix%iown = phaml_matrix%iown(1:petsc_matrix%my_total_eq)
petsc_matrix%dirich = .false.

! Create the matrix and rhs data structures

petsc_matrix%my_own_eq = count(petsc_matrix%iown)

if (still_sequential) then
   call MatCreateShell(PETSC_COMM_SELF,petsc_matrix%my_own_eq, &
                       petsc_matrix%my_own_eq,PETSC_DECIDE, &
                       PETSC_DECIDE,null_data,petsc_matrix%A,jerr);CHKERRQ(jerr)
   call VecCreateSeq(PETSC_COMM_SELF,petsc_matrix%my_own_eq,petsc_matrix%b, &
                     jerr);CHKERRQ(jerr)
else
   call MatCreateShell(PETSC_COMM_WORLD,petsc_matrix%my_own_eq, &
                       petsc_matrix%my_own_eq,PETSC_DECIDE, &
                       PETSC_DECIDE,null_data,petsc_matrix%A,jerr);CHKERRQ(jerr)
   call VecCreateMPI(PETSC_COMM_WORLD,petsc_matrix%my_own_eq,PETSC_DECIDE, &
                     petsc_matrix%b,jerr);CHKERRQ(jerr)
endif

! Get the owned range for this processor and total number of equations.

call VecGetOwnershipRange(petsc_matrix%b,petsc_matrix%my_global_low, &
                          petsc_matrix%my_global_hi,jerr);CHKERRQ(jerr)
call VecGetSize(petsc_matrix%b,petsc_matrix%global_eq,jerr);CHKERRQ(jerr)

! copy the values from rhs to b

call VecSetValues(petsc_matrix%b,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  pack(phaml_matrix%rhs(1:petsc_matrix%my_total_eq),mask=petsc_matrix%iown), &
                  INSERT_VALUES,jerr);CHKERRQ(jerr)
call VecAssemblyBegin(petsc_matrix%b,jerr);CHKERRQ(jerr)
call VecAssemblyEnd(petsc_matrix%b,jerr);CHKERRQ(jerr)

! set the matrix multiply callback routine

call MatShellSetOperation(petsc_matrix%A,MATOP_MULT,matmult_mf, &
                          jerr);CHKERRQ(jerr)

nullify(petsc_matrix%petsc_index)

end subroutine create_petsc_linear_system_mf

!          ----------------
subroutine change_petsc_rhs(rhs,petsc_matrix)
!          ----------------

!----------------------------------------------------
! This routine copies a PHAML right hand side to a PETSc right hand side.
! This requires that the Dirichlet boundary conditions be homogeneous because
! it does not know the values in the Dirichlet columns of the matrix.  At the
! time of the writing (11/18/11) this routine is only called by the
! eigensolvers, which require homogeneous boundary conditions.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: rhs(:)
type(petsc_matrix_type), intent(inout) :: petsc_matrix
!----------------------------------------------------
! Local variables:

integer :: i
PetscErrorCode jerr
!----------------------------------------------------
! Begin executable code

! copy the values from rhs to b

call VecSetValues(petsc_matrix%b,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  pack(rhs,mask=(petsc_matrix%iown(:size(rhs)) .and. &
                                 .not.petsc_matrix%dirich(:size(rhs)))), &
                  INSERT_VALUES,jerr);CHKERRQ(jerr)
call VecAssemblyBegin(petsc_matrix%b,jerr);CHKERRQ(jerr)
call VecAssemblyEnd(petsc_matrix%b,jerr);CHKERRQ(jerr)

end subroutine change_petsc_rhs

!          -----------
subroutine petsc_solve(phaml_matrix,petsc_matrix,solver_cntl,io_cntl, &
                       still_sequential,grid,procs)
!          -----------

!----------------------------------------------------
! This routine solves the linear system by the specified method in PETSc
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout), target :: phaml_matrix
type(petsc_matrix_type), intent(inout), target :: petsc_matrix
type(solver_options), intent(in), target :: solver_cntl
type(io_options), intent(in), target :: io_cntl
logical, intent(in) :: still_sequential
type(grid_type), intent(in), target :: grid
type(proc_info), intent(in), target :: procs
!----------------------------------------------------
!----------------------------------------------------
! Local variables:

integer :: eq, i
PC :: pc
KSP :: ksp
PetscErrorCode jerr
Vec :: x
PetscInt :: niter
!----------------------------------------------------
! Begin executable code

! If there is only one equation, just set the solution

if (petsc_matrix%global_eq == 1) then
   eq = 1
   do while (phaml_matrix%equation_type(eq) == DIRICHLET)
      eq = eq + 1
   end do
   phaml_matrix%solution(eq) = phaml_matrix%rhs(eq)
   do i=phaml_matrix%begin_row(eq)+1,phaml_matrix%end_row(eq)
      phaml_matrix%solution(eq) = phaml_matrix%solution(eq) - &
                             phaml_matrix%matrix_val(i) * &
                             phaml_matrix%solution(phaml_matrix%column_index(i))
   end do
   phaml_matrix%solution(eq) = phaml_matrix%solution(eq) / &
                             phaml_matrix%matrix_val(phaml_matrix%begin_row(eq))
   return
endif

! For later access

call set_petsc_hold(petsc_matrix,phaml_matrix,grid,procs,io_cntl,solver_cntl, &
                    still_sequential)

! set up the solver context

if (still_sequential) then
   call KSPCreate(PETSC_COMM_SELF,ksp,jerr);CHKERRQ(jerr)
else
   call KSPCreate(PETSC_COMM_WORLD,ksp,jerr);CHKERRQ(jerr)
endif
call KSPSetOperators(ksp,petsc_matrix%A,petsc_matrix%A, &
                      DIFFERENT_NONZERO_PATTERN,jerr);CHKERRQ(jerr)
call VecDuplicate(petsc_matrix%b,x,jerr);CHKERRQ(jerr)

! select the solver method

call set_petsc_solver(ksp,solver_cntl)

! select the preconditioning method

call set_petsc_preconditioner(ksp,solver_cntl)

! set PETSc options

call set_petsc_options(ksp,solver_cntl,io_cntl)

! solve the system

call KSPSolve(ksp,petsc_matrix%b,x,jerr);CHKERRQ(jerr)
if (jerr /= 0) call warning("PETSc solver returned error code", &
                            intlist=(/jerr/))

! extract the solution

call petscvec_to_phaml(x,phaml_matrix%solution(1:),petsc_matrix)

! get the number of iterations

call KSPGetIterationNumber(ksp,niter,jerr);CHKERRQ(jerr)
phaml_matrix%solver_niter = niter

! free memory

call KSPDestroy(ksp,jerr);CHKERRQ(jerr)
call VecDestroy(x,jerr);CHKERRQ(jerr)

end subroutine petsc_solve

!          ----------------
subroutine set_petsc_solver(ksp,solver_cntl)
!          ----------------

!----------------------------------------------------
! This routine sets the PETSc solver in ksp
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

KSP ksp
type(solver_options), intent(in) :: solver_cntl
!----------------------------------------------------
! Local variables:

PetscErrorCode jerr
!----------------------------------------------------
! Begin executable code

select case (solver_cntl%solver)
case (PETSC_RICHARDSON_SOLVER)
   call KSPSetType(ksp,KSPRICHARDSON,jerr);CHKERRQ(jerr)
case (PETSC_CHEBYCHEV_SOLVER)
   call KSPSetType(ksp,KSPCHEBYCHEV,jerr);CHKERRQ(jerr)
case (PETSC_CG_SOLVER)
   call KSPSetType(ksp,KSPCG,jerr);CHKERRQ(jerr)
case (PETSC_GMRES_SOLVER)
   call KSPSetType(ksp,KSPGMRES,jerr);CHKERRQ(jerr)
case (PETSC_TCQMR_SOLVER)
   call KSPSetType(ksp,KSPTCQMR,jerr);CHKERRQ(jerr)
case (PETSC_BCGS_SOLVER)
   call KSPSetType(ksp,KSPBCGS,jerr);CHKERRQ(jerr)
case (PETSC_CGS_SOLVER)
   call KSPSetType(ksp,KSPCGS,jerr);CHKERRQ(jerr)
case (PETSC_TFQMR_SOLVER)
   call KSPSetType(ksp,KSPTFQMR,jerr);CHKERRQ(jerr)
case (PETSC_CR_SOLVER)
   call KSPSetType(ksp,KSPCR,jerr);CHKERRQ(jerr)
case (PETSC_LSQR_SOLVER)
   call KSPSetType(ksp,KSPLSQR,jerr);CHKERRQ(jerr)
case (PETSC_BICG_SOLVER)
   call KSPSetType(ksp,KSPBICG,jerr);CHKERRQ(jerr)
case (PETSC_MUMPS_GEN_SOLVER, PETSC_MUMPS_SPD_SOLVER, PETSC_SUPERLU_SOLVER)
   call KSPSetType(ksp,KSPPREONLY,jerr);CHKERRQ(jerr)
case default
   ierr = USER_INPUT_ERROR
   call fatal("illegal solver choice in set_petsc_solver", &
              intlist=(/solver_cntl%solver/))
   return
end select

end subroutine set_petsc_solver

!          ------------------------
subroutine set_petsc_preconditioner(ksp,solver_cntl)
!          ------------------------

!----------------------------------------------------
! This routine set the PETSc preconditioner in ksp
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

KSP ksp
type(solver_options), intent(in) :: solver_cntl
!----------------------------------------------------
! Local variables:

PC pc
PetscErrorCode jerr
Mat A, F
PetscInt ival, icntl
!----------------------------------------------------
! Begin executable code

call KSPGetPC(ksp,pc,jerr);CHKERRQ(jerr)

if (solver_cntl%solver == PETSC_MUMPS_GEN_SOLVER .or. &
    solver_cntl%solver == PETSC_MUMPS_SPD_SOLVER) then
#ifdef PETSC_HAVE_MUMPS
   if (solver_cntl%solver == PETSC_MUMPS_GEN_SOLVER) then
      call PCSetType(pc,PCLU,jerr);CHKERRQ(jerr)
   else
      call KSPGetOperators(ksp,A,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                           jerr);CHKERRQ(jerr)
      call MatSetOption(A,MAT_SPD,PETSC_TRUE,jerr);CHKERRQ(jerr)
      call PCSetType(pc,PCCHOLESKY,jerr);CHKERRQ(jerr)
   endif
   call PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS,jerr);CHKERRQ(jerr)
   call PCFactorSetUpMatSolverPackage(pc,jerr);CHKERRQ(jerr)
   call PCFactorGetMatrix(pc,F,jerr);CHKERRQ(jerr)
   icntl=7
   ival=2
   call MatMumpsSetIcntl(F,icntl,ival,jerr);CHKERRQ(ierr)
#else
   ierr = USER_INPUT_ERROR
   call fatal("Cannot use PETSC_MUMPS solver because PETSc was not built with MUMPS support")
   stop
#endif

elseif (solver_cntl%solver == PETSC_SUPERLU_SOLVER) then
   if (PARALLEL == SEQUENTIAL) then
#ifdef PETSC_HAVE_SUPERLU
      call PCSetType(pc,PCLU,jerr);CHKERRQ(jerr)
! SuperLU cannot handle a system with only one equation
      call KSPGetOperators(ksp,A,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                           jerr);CHKERRQ(jerr)
      call MatGetSize(A,ival,PETSC_NULL_INTEGER,jerr);CHKERRQ(jerr)
      if (ival <= 1) then
         call PCFactorSetMatSolverPackage(pc,MATSOLVERPETSC,jerr);CHKERRQ(jerr)
      else
         call PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU, &
                                          jerr);CHKERRQ(jerr)
      endif
      call PCFactorSetUpMatSolverPackage(pc,jerr);CHKERRQ(jerr)
      call PCFactorGetMatrix(pc,F,jerr);CHKERRQ(jerr)
#else
      ierr = USER_INPUT_ERROR
      call fatal("Cannot use PETSC_SUPERLU solver because PETSc was not built with SUPERLU support")
      stop
#endif
   else
#ifdef PETSC_HAVE_SUPERLU_DIST
      call PCSetType(pc,PCLU,jerr);CHKERRQ(jerr)
! SuperLU cannot handle a system with only one equation
      call KSPGetOperators(ksp,A,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                           jerr);CHKERRQ(jerr)
      call MatGetSize(A,ival,PETSC_NULL_INTEGER,jerr);CHKERRQ(jerr)
      if (ival <= 1) then
         call PCFactorSetMatSolverPackage(pc,MATSOLVERPETSC,jerr);CHKERRQ(jerr)
      else
         call PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST, &
                                          jerr);CHKERRQ(jerr)
      endif
      call PCFactorSetUpMatSolverPackage(pc,jerr);CHKERRQ(jerr)
      call PCFactorGetMatrix(pc,F,jerr);CHKERRQ(jerr)
#else
      ierr = USER_INPUT_ERROR
      call fatal("Cannot use PETSC_SUPERLU solver because PETSc was not built with SUPERLU support")
      stop
#endif
   endif

else

   select case (solver_cntl%preconditioner)
   case (NO_PRECONDITION)
      call PCSetType(pc,PCSHELL,jerr);CHKERRQ(jerr)
      call PCShellSetApply(pc,precon_no,null_data,jerr);CHKERRQ(jerr)
   case (MG_PRECONDITION)
      call PCSetType(pc,PCSHELL,jerr);CHKERRQ(jerr)
      call PCShellSetApply(pc,precon_mg,null_data,jerr);CHKERRQ(jerr)
   case (FUDOP_DD_PRECONDITION)
      call PCSetType(pc,PCSHELL,jerr);CHKERRQ(jerr)
      call PCShellSetApply(pc,precon_fudop_dd,null_data,jerr);CHKERRQ(jerr)
   case (COARSE_GRID_PRECONDITION)
      call PCSetType(pc,PCSHELL,jerr);CHKERRQ(jerr)
      call PCShellSetApply(pc,precon_coarse,null_data,jerr);CHKERRQ(jerr)
   case (TEST_PRECONDITION)
      call PCSetType(pc,PCSHELL,jerr);CHKERRQ(jerr)
      call PCShellSetApply(pc,precon_test,null_data,jerr);CHKERRQ(jerr)
   case (PETSC_JACOBI_PRECONDITION)
      call PCSetType(pc,PCJACOBI,jerr);CHKERRQ(jerr)
   case (PETSC_BJACOBI_PRECONDITION)
      call PCSetType(pc,PCBJACOBI,jerr);CHKERRQ(jerr)
   case (PETSC_SOR_PRECONDITION)
      call PCSetType(pc,PCSOR,jerr);CHKERRQ(jerr)
   case (PETSC_EISENSTAT_PRECONDITION)
      call PCSetType(pc,PCEISENSTAT,jerr);CHKERRQ(jerr)
   case (PETSC_ICC_PRECONDITION)
      call PCSetType(pc,PCICC,jerr);CHKERRQ(jerr)
   case (PETSC_ILU_PRECONDITION)
      call PCSetType(pc,PCILU,jerr);CHKERRQ(jerr)
   case (PETSC_ASM_PRECONDITION)
      call PCSetType(pc,PCASM,jerr);CHKERRQ(jerr)
   case default
      ierr = USER_INPUT_ERROR
      call fatal("illegal preconditioner choice in set_petsc_preconditioner", &
                 intlist=(/solver_cntl%preconditioner/))
      return
   end select

endif

end subroutine set_petsc_preconditioner

!          -----------------
subroutine set_petsc_options(ksp,solver_cntl,io_cntl)
!          -----------------

!----------------------------------------------------
! This routine set the PETSc options in ksp
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

KSP ksp
type(solver_options), intent(in) :: solver_cntl
type(io_options), intent(in) :: io_cntl
!----------------------------------------------------
! Local variables:

real(kind(0.0d0)) :: temp1, temp2, temp3
integer :: temp4, temp5
PetscErrorCode jerr
PC pc

! TEMP101110 gfortran bug?
! Rediscovered 7/2/2012; GNU Fortran (GCC) 4.1.2 20080704 (Red Hat 4.1.2-52)
! It thinks these two should be functions in this module, i.e.
! undefined reference to `__petsc_interf__kspmonitortrueresidualnorm'
! These external statements are in petscksp.h and petscsys.h respectively,
! but the external statements need to be in this subroutine, not at the
! module level.

!external KSPMonitorTrueResidualNorm
!external PETSC_NULL_FUNCTION
!----------------------------------------------------
! Begin executable code

call KSPGetPC(ksp,pc,jerr);CHKERRQ(jerr)

if (solver_cntl%petsc_cntl%petsc_richardson_damping_factor /= huge(0.0d0)) then
   call KSPRichardsonSetScale(ksp, &
                      solver_cntl%petsc_cntl%petsc_richardson_damping_factor, &
                      jerr);CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_chebychev_emin /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_chebychev_emin /= huge(0.0d0)) then
   temp1 = solver_cntl%petsc_cntl%petsc_chebychev_emin
   if (temp1 == huge(0.0d0)) temp1 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp2 = solver_cntl%petsc_cntl%petsc_chebychev_emax
   if (temp2 == huge(0.0d0)) temp2 = PETSC_DEFAULT_DOUBLE_PRECISION
   call KSPChebychevSetEigenvalues(ksp,temp1,temp2,jerr);CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_gmres_max_steps /= huge(0)) then
   call KSPGMRESSetRestart(ksp,solver_cntl%petsc_cntl%petsc_gmres_max_steps, &
                           jerr);CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_rtol /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_atol /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_dtol /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_maxits /= huge(0)) then
   temp1 = solver_cntl%petsc_cntl%petsc_rtol
   if (temp1 == huge(0.0d0)) temp1 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp2 = solver_cntl%petsc_cntl%petsc_atol
   if (temp2 == huge(0.0d0)) temp2 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp3 = solver_cntl%petsc_cntl%petsc_dtol
   if (temp3 == huge(0.0d0)) temp3 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp4 = solver_cntl%petsc_cntl%petsc_maxits
   if (temp4 == huge(0)) temp4 = PETSC_DEFAULT_INTEGER
   call KSPSetTolerances(ksp,temp1,temp2,temp3,temp4,jerr);CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_ilu_levels /= huge(0)) then

! For PETSc versions before 2.3.1
!   call PCILUSetLevels(pc,solver_cntl%petsc_cntl%petsc_ilu_levels,jerr);CHKERRQ(jerr)

! For PETSc version 2.3.1 and later
   call PCFactorSetLevels(pc,solver_cntl%petsc_cntl%petsc_ilu_levels,jerr)
   CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_icc_levels /= huge(0)) then

! For PETSc versions before 2.3.1
!   call PCICCSetLevels(pc,solver_cntl%petsc_cntl%petsc_icc_levels,jerr);CHKERRQ(jerr)

! For PETSc version 2.3.1 and later
   call PCFactorSetLevels(pc,solver_cntl%petsc_cntl%petsc_icc_levels,jerr)
   CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_sor_omega /= huge(0.0d0)) then
   call PCSORSetOmega(pc,solver_cntl%petsc_cntl%petsc_sor_omega,jerr);CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_sor_its /= huge(0) .or. &
    solver_cntl%petsc_cntl%petsc_sor_lits /= huge(0)) then
   temp4 = solver_cntl%petsc_cntl%petsc_sor_its
   if (temp4 == huge(0)) temp4 = PETSC_DEFAULT_INTEGER
   temp5 = solver_cntl%petsc_cntl%petsc_sor_lits
   if (temp5 == huge(0)) temp5 = PETSC_DEFAULT_INTEGER
   call PCSORSetIterations(pc,temp4,temp5,jerr);CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_eisenstat_nodiagscaling) then
   call PCEisenstatNoDiagonalScaling()
endif

if (solver_cntl%petsc_cntl%petsc_eisenstat_omega /= huge(0.0d0)) then
   call PCEisenstatSetOmega(pc,solver_cntl%petsc_cntl%petsc_eisenstat_omega, &
                            jerr);CHKERRQ(jerr)
endif

if (solver_cntl%petsc_cntl%petsc_asm_overlap /= huge(0)) then
   call PCASMSetOverlap(pc,solver_cntl%petsc_cntl%petsc_asm_overlap,jerr)
   CHKERRQ(jerr)
endif

! Print L2 norm of residual after each iteration

if (io_cntl%print_error_when == FREQUENTLY .or. &
    io_cntl%print_error_when == TOO_MUCH) then

! For PETSc versions before 2.3.3
!   call KSPSetMonitor(ksp,KSPTrueMonitor,PETSC_NULL_OBJECT, &
!                       PETSC_NULL_FUNCTION,jerr);CHKERRQ(jerr)

! For PETSc version 2.3.3 and later
  call KSPMonitorSet(ksp,KSPMonitorTrueResidualNorm,PETSC_NULL_OBJECT, &
                     PETSC_NULL_FUNCTION,jerr);CHKERRQ(jerr)

endif
end subroutine set_petsc_options

!          ---------------------------
subroutine destroy_petsc_linear_system(petsc_matrix)
!          ---------------------------

!----------------------------------------------------
! This routine frees the space for a PETSc matrix and right hand side
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(petsc_matrix_type), intent(inout) :: petsc_matrix
!----------------------------------------------------
! Local variables:

integer :: astat
PetscErrorCode jerr
!----------------------------------------------------
! Begin executable code

call MatDestroy(petsc_matrix%A,jerr);CHKERRQ(jerr)
call VecDestroy(petsc_matrix%b,jerr);CHKERRQ(jerr)
if (associated(petsc_matrix%iown)) deallocate(petsc_matrix%iown,stat=astat)
if (associated(petsc_matrix%dirich)) deallocate(petsc_matrix%dirich,stat=astat)
if (associated(petsc_matrix%petsc_index)) deallocate(petsc_matrix%petsc_index,stat=astat)

end subroutine destroy_petsc_linear_system

!          -----------------
subroutine petscvec_to_phaml(petscvec,phamlvec,petsc_matrix)
!          -----------------

!----------------------------------------------------
! This routine copies a PETSc Vec to a PHAML distributed vector
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real) :: phamlvec(:)
Vec :: petscvec
type(petsc_matrix_type), intent(in) :: petsc_matrix
!----------------------------------------------------
! Local variables:

PetscScalar :: x_array(1)
PetscOffset :: i_x
PetscErrorCode jerr
integer :: pi, i
!----------------------------------------------------
! Begin executable code

! Start with the phaml solution vector to get Dirichlet boundary values

phamlvec = petsc_hold%phaml_matrix%solution(1:size(phamlvec))

! Extract the local values from petscvec

call VecGetArray(petscvec,x_array,i_x,jerr);CHKERRQ(jerr)

! Copy the ones I own to the PHAML vector

pi = 0
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i) .and. .not.petsc_matrix%dirich(i)) then
      pi = pi + 1
      phamlvec(i) = x_array(i_x + pi)
   endif
end do

! Request the unowned ones from other processors

if (.not. petsc_hold%still_sequential) then
   call update_shadows(petsc_matrix,petsc_hold%phaml_matrix%gid, &
                       petsc_hold%procs, 1604, 1605, 1606, rdata=phamlvec)
endif

! free memory

call VecRestoreArray(petscvec,x_array,i_x,jerr);CHKERRQ(jerr)

end subroutine petscvec_to_phaml

!          -----------------
subroutine phamlvec_to_petsc(phamlvec,petscvec,petsc_matrix)
!          -----------------

!----------------------------------------------------
! This routine copies a PHAML distributed vector to a PETSC vec
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real) :: phamlvec(:)
Vec :: petscvec
type(petsc_matrix_type), intent(in) :: petsc_matrix
!----------------------------------------------------
! Local variables:

integer :: i
PetscErrorCode jerr
!----------------------------------------------------
! Begin executable code

call VecSetValues(petscvec,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  pack(phamlvec,mask=(petsc_matrix%iown .and. &
                                      .not. petsc_matrix%dirich)), &
                  INSERT_VALUES,jerr);CHKERRQ(jerr)
call VecAssemblyBegin(petscvec,jerr);CHKERRQ(jerr)
call VecAssemblyEnd(petscvec,jerr);CHKERRQ(jerr)

end subroutine phamlvec_to_petsc

!          --------------
subroutine update_shadows(petsc_matrix,gid,procs,tag1,tag2,tag3,rdata,idata)
!          --------------

!----------------------------------------------------
! This routine updates the values of rdata and/or idata at shadow points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(petsc_matrix_type), intent(in) :: petsc_matrix
type(hash_key_eq), intent(in) :: gid(:)
type(proc_info), intent(in) :: procs
integer, intent(in) :: tag1, tag2, tag3
real(my_real), intent(inout), optional :: rdata(:)
integer, intent(inout), optional :: idata(:)
!----------------------------------------------------
! Local variables:

integer :: counter, i, p, lid, astat, nproc, my_processor, isub, rsub, &
           oisub, orsub, limit, KEY_SIZE_EQ
! newcomm
integer :: nisend
integer, allocatable :: isend(:), nisendv(:), nirecv(:), nrsend(:), nrrecv(:)
integer, pointer :: irecv(:)
real(my_real), allocatable :: rsend(:)
real(my_real), pointer :: rrecv(:)
!----------------------------------------------------
! Begin executable code

KEY_SIZE_EQ = KEY_SIZE+1

! allocate space for received messages

nproc = num_proc(procs)
my_processor = my_proc(procs)
allocate(nisendv(nproc),nirecv(nproc),nrsend(nproc),nrrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in update_shadows")
   return
endif

! this is overallocated if Dirichlet rows are not included
nisend = (petsc_matrix%my_total_eq-petsc_matrix%my_own_eq)*KEY_SIZE_EQ
allocate(isend(nisend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in update_shadows")
   return
endif

! Make a list of the ones this processor doesn't own and are not Dirichlet

counter=1
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i)) cycle
   if (petsc_matrix%dirich(i)) cycle
   call hash_pack_key(gid(i),isend,counter)
   counter = counter + KEY_SIZE_EQ
end do
nisend = counter-1

! Send the request

call phaml_alltoall(procs,isend,nisend,irecv,nirecv,tag1)

! Reply with ones I own

! Count the number of responses

counter = 0
isub = 1
do p=1,nproc
   if (p == my_processor) then
      isub = isub + nirecv(p)
      cycle
   endif
   do i=1,nirecv(p)/KEY_SIZE_EQ
      lid = hash_decode_key(hash_unpack_key(irecv,isub,.true.), &
                            petsc_matrix%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (petsc_matrix%iown(lid)) then
            counter = counter + 1
         endif
      endif
      isub = isub + KEY_SIZE_EQ
   end do
end do

! allocate memory

deallocate(isend)
if (present(idata)) then
   allocate(isend(counter*(KEY_SIZE_EQ+1)), stat=astat)
else
   allocate(isend(counter*KEY_SIZE_EQ), stat=astat)
endif
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in update_shadows")
   return
endif
if (present(rdata)) then
   allocate(rsend(counter), stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in update_shadows")
      return
   endif
endif

! Make the arrays with the responses

counter = 1
isub = 1
oisub = isub
rsub = 1
orsub = rsub
do p=1,nproc
   if (p == my_processor) then
      nisendv(p) = 0
      nrsend(p) = 0
      counter = counter + nirecv(p)
      cycle
   endif
   do i=1,nirecv(p)/KEY_SIZE_EQ
      lid = hash_decode_key(hash_unpack_key(irecv,counter,.true.), &
                            petsc_matrix%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (petsc_matrix%iown(lid)) then
            isend(isub:isub+KEY_SIZE_EQ-1) = irecv(counter:counter+KEY_SIZE_EQ-1)
            isub = isub + KEY_SIZE_EQ
            if (present(rdata)) then
               rsend(rsub) = rdata(lid)
               rsub = rsub + 1
            endif
            if (present(idata)) then
               isend(isub) = idata(lid)
               isub = isub + 1
            endif
         endif
      endif
      counter = counter + KEY_SIZE_EQ
   end do
   nisendv(p) = isub - oisub
   oisub = isub
   nrsend(p) = rsub - orsub
   orsub = rsub
end do

if (associated(irecv)) deallocate(irecv,stat=astat)

! Send the replies

call phaml_alltoall(procs,isend,nisendv,irecv,nirecv,tag2)
if (present(rdata)) then
   call phaml_alltoall(procs,rsend,nrsend,rrecv,nrrecv,tag3)
else
   nullify(rrecv)
endif

deallocate(isend,rsend,stat=astat)

! Set the shadow values from the replies

isub = 1
rsub = 1
do p=1,nproc
   if (p == my_processor) cycle
   if (present(idata)) then
      limit = nirecv(p)/(KEY_SIZE_EQ+1)
   else
      limit = nirecv(p)/KEY_SIZE_EQ
   endif
   do i=1,limit
      lid = hash_decode_key(hash_unpack_key(irecv,isub,.true.), &
                            petsc_matrix%eq_hash)
      isub = isub+KEY_SIZE_EQ
      if (lid == HASH_NOT_FOUND) then
         call warning("received reply for an equation I don't have in update_shadows")
      else
         if (present(rdata)) then
            rdata(lid) = rrecv(rsub)
            rsub = rsub + 1
         endif
         if (present(idata)) then
            idata(lid) = irecv(isub)
            isub = isub + 1
         endif
      endif
   end do
end do

if (associated(irecv)) deallocate(irecv,stat=astat)
if (associated(rrecv)) deallocate(rrecv,stat=astat)
deallocate(nisendv,nirecv,nrsend,nrrecv,stat=astat)

end subroutine update_shadows

!          ---------
subroutine precon_no(nodata,avec,bvec,jerr)
!          ---------

!----------------------------------------------------
! This routine applies no preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
PetscErrorCode jerr
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call VecCopy(avec,bvec,jerr);CHKERRQ(jerr)

end subroutine precon_no

!          ---------
subroutine precon_mg(nodata,avec,bvec,jerr)
!          ---------

!----------------------------------------------------
! This routine provides the callback for V-cycle multigrid preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
PetscErrorCode jerr
!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: invec(:),outvec(:) 
integer :: isize

!----------------------------------------------------
! Begin executable code

! TEMP 04/04/08 the condensed number of equations should be sufficient (ifort
!      with array bounds checking ran fine), but gfortran gives an error during
!      free unless I use the full number of equations

isize = petsc_hold%phaml_matrix%neq_full

allocate(invec(isize),outvec(isize))

! convert PETSc avec to a normal vector, apply preconditioner, and convert back

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call mg_precon(invec,outvec,petsc_hold%phaml_matrix, &
               petsc_hold%grid,petsc_hold%procs,petsc_hold%io_cntl, &
               petsc_hold%solver_cntl,petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

deallocate(invec,outvec)

end subroutine precon_mg

!          ---------------
subroutine precon_fudop_dd(nodata,avec,bvec,jerr)
!          ---------------

!----------------------------------------------------
! This routine provides the callback for FuDoP domain decomposition
! preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
PetscErrorCode jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call lapack_precon(invec,outvec,FUDOP_DD_PRECONDITION,petsc_hold%phaml_matrix, &
                   petsc_hold%procs,petsc_hold%solver_cntl, &
                   petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine precon_fudop_dd

!          -------------
subroutine precon_coarse(nodata,avec,bvec,jerr)
!          -------------

!----------------------------------------------------
! This routine provides the callback for coarse grid preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
PetscErrorCode jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call lapack_precon(invec,outvec,COARSE_GRID_PRECONDITION, &
                   petsc_hold%phaml_matrix, &
                   petsc_hold%procs,petsc_hold%solver_cntl, &
                   petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine precon_coarse

!          -----------
subroutine precon_test(nodata,avec,bvec,jerr)
!          -----------

!----------------------------------------------------
! This routine provides the callback for testing a new preconditioner.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
PetscErrorCode jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call test_precon(invec,outvec,petsc_hold%phaml_matrix,petsc_hold%grid, &
                 petsc_hold%procs,petsc_hold%solver_cntl, &
                 petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine precon_test

!          -------------
subroutine matmult_mf(matrix,avec,bvec,jerr)
!          -------------

!----------------------------------------------------
! This routine provides the callback for matrix free matmult
! as an interface between PETSc and an external routine in linsys.f90.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Mat :: matrix
Vec :: avec, bvec
PetscErrorCode jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call matrix_times_vector(invec,outvec,petsc_hold%phaml_matrix, &
                         petsc_hold%procs,petsc_hold%still_sequential, &
                         1621,1622,1623,1624,1625,1626)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine matmult_mf

end module petsc_interf
