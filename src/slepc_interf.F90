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
!     Mathematical and Computational Sciences Division                !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!---------------------------------------------------------------------!

module slepc_interf

!----------------------------------------------------
! This module contains routines that interface to SLEPc.
!
! communication tags in this module are of the form ??xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use petsc_interf
use petsc_type_mod
use linsystype_mod
use gridtype_mod
use sort_mod
use message_passing

!----------------------------------------------------

implicit none
private
public eigen_slepc

!----------------------------------------------------
! The PETSc and SLEPc include files.  Note the use of preprocessor #include
! instead of the Fortran include statement, because the include files contain
! preprocessor directives.

#include "petscversion.h"
#include "finclude/petsc.h"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0))
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#endif
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"

contains

!          -----------
subroutine eigen_slepc(linear_system,procs,solver_control,io_control, &
                       still_sequential,eigenvalues,monitor)
!          -----------

!----------------------------------------------------
! This routine solves the eigenvalue problem using SLEPc
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
real(my_real), intent(out) :: eigenvalues(:)
type(eigen_monitor), intent(inout) :: monitor

!----------------------------------------------------
! Local variables:

EPS eps
ST st
KSP ksp
PC pc
type(petsc_matrix_type) :: petsc_stiffness, petsc_mass
Vec xr, xi, init_space(solver_control%num_eval)
PetscScalar kr, ki
PetscInt j, nconv, nev, ncv, maxit, numit, neq, init_space_dim
PetscReal error, tol, norm
PetscScalar target, shift
PetscErrorCode jerr

integer :: i, iperm(solver_control%num_eval), eigensolver
real(my_real), parameter :: near_roundoff = 1.0e-13_my_real
logical :: reverse_M_A

! TEMP120702 apparent gfortran bug, see comments in petsc_interf.F90
!external EPSMONITORCONVERGED
!external EPSMONITORFIRST
!external PETSC_NULL_FUNCTION
!----------------------------------------------------
! Begin executable code

! Master does not participate

if (my_proc(procs) == MASTER) return

! Local copy of eigensolver, so it can be changed if necessary.

eigensolver = solver_control%eigensolver

! The following "held" variables are used in petscvec_to_phaml

call set_petsc_hold(phaml_matrix=linear_system,procs=procs, &
                    still_sequential=still_sequential)

! Create the PETSc version of the stiffness and mass matrices

linear_system%matrix_val => linear_system%mass
call create_petsc_linear_system(linear_system,petsc_mass, &
                                solver_control,still_sequential,procs)
linear_system%matrix_val => linear_system%stiffness
call create_petsc_linear_system(linear_system,petsc_stiffness, &
                                solver_control,still_sequential,procs)

! PETSc version of neq

neq = petsc_stiffness%global_eq

! Create vectors that are compatible with the matrix

call MatGetVecs(petsc_stiffness%A,PETSC_NULL_OBJECT,xr,jerr);CHKERRQ(jerr)
call MatGetVecs(petsc_stiffness%A,PETSC_NULL_OBJECT,xi,jerr);CHKERRQ(jerr)
do i=1,solver_control%num_eval
   call MatGetVecs(petsc_stiffness%A,PETSC_NULL_OBJECT,init_space(i), &
                   jerr);CHKERRQ(jerr)
end do

! Create eigensolver context

if (still_sequential) then
   call EPSCreate(PETSC_COMM_SELF,eps,jerr);CHKERRQ(jerr)
else
   call EPSCreate(PETSC_COMM_WORLD,eps,jerr);CHKERRQ(jerr)
endif

! Set operators

! For methods that can only find the largest eigenvalues,
! solve Mx = (1/lambda) Ax to get smallest magnitude eigenvalues.  Note this
! is not the same as the smallest eigenvalues if there are negative eigenvalues.

if (eigensolver == SLEPC_POWER .or. &
    eigensolver == SLEPC_SUBSPACE) then
   call EPSSetOperators(eps,petsc_mass%A,petsc_stiffness%A,jerr);CHKERRQ(jerr)
   reverse_M_A = .true.
else
   call EPSSetOperators(eps,petsc_stiffness%A,petsc_mass%A,jerr);CHKERRQ(jerr)
   reverse_M_A = .false.
endif

! Set eigenproblem type

! nonsymmetric

if (solver_control%pde_has_first_order_terms .or. &
    solver_control%pde_has_cross_derivative) then
   select case (eigensolver)
   case (SLEPC_LANCZOS, SLEPC_BLOPEX)
      ierr = USER_INPUT_ERROR
      call fatal("the selected eigensolver requires a Hermitian matrix, i.e. no first order or mixed derivative terms in the PDE")
      stop
   case default
      call EPSSetProblemType(eps,EPS_PGNHEP,jerr);CHKERRQ(jerr)
   end select

! symmetric

else
   call EPSSetProblemType(eps,EPS_GHEP,jerr);CHKERRQ(jerr)
endif

! Set the number of eigenvalues to compute and the number of column vectors.
! The fourth argument in SetDimensions is for extremely large number of
! eigenvalues; we'll assume it doesn't matter.

nev = solver_control%num_eval

! For ncv, use the default given in the SLEPc manual, if not user specified

if (solver_control%eigen_cntl%ncv <= 0) then
   select case (eigensolver)
   case (SLEPC_POWER, SLEPC_BLOPEX)
      ncv = nev
   case (SLEPC_SUBSPACE, SLEPC_ARNOLDI, SLEPC_LANCZOS, SLEPC_KRYLOV_SCHUR)
      ncv = max(2*nev,nev+15)
   case (SLEPC_GEN_DAVIDSON, SLEPC_JACOBI_DAVIDSON)
      ncv = max(2*nev,nev+15) + 1
   case (SLEPC_LAPACK)
      ncv = neq
   case (SLEPC_ARPACK)
      ncv = max(20,2*nev+1)
   case default
      ierr = USER_INPUT_ERROR
      call fatal("illegal selection for eigensolver", &
                 intlist=(/eigensolver/))
      stop
   end select
else
   ncv = solver_control%eigen_cntl%ncv
endif

call EPSSetDimensions(eps,nev,ncv,PETSC_DECIDE,jerr);CHKERRQ(jerr)

! Some cases where an eigensolver cannot handle extremely small matrices.
! Use LAPACK instead.

! The Davidson methods require the dimension of the matrix to be at least 2.
   if ((eigensolver == SLEPC_GEN_DAVIDSON .or. &
        eigensolver == SLEPC_JACOBI_DAVIDSON) .and. &
        neq < 2) eigensolver = SLEPC_LAPACK

! ARPACK needs the dimension of the matrix to be at least nev+2.
   if (eigensolver == SLEPC_ARPACK .and. &
       neq < nev+2) eigensolver = SLEPC_LAPACK

! Set the target eigenvalues

! smallest eigenvalues

if (solver_control%lambda0 == -huge(0.0_my_real)) then

   if (reverse_M_A) then
      call EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE,jerr);CHKERRQ(jerr)
   else
      call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL,jerr);CHKERRQ(jerr)
   endif

! closest to lambda0

else

   select case (eigensolver)
   case (SLEPC_POWER, SLEPC_SUBSPACE, SLEPC_BLOPEX)
      ierr = USER_INPUT_ERROR
      call fatal("cannot ask for interior eigenvalues with the selected eigensolver")
      stop
   case default
      target = solver_control%lambda0
      call EPSSetTarget(eps,target,jerr);CHKERRQ(jerr)
      call EPSSetWhichEigenpairs(eps,EPS_TARGET_REAL,jerr);CHKERRQ(jerr)
   end select

endif

! Set the maximum number of iterations and tolerance.  For maxit, use the
! default in the SLEPc users guide, if not user specified.  Note the tolerance
! is on ||r||/|lambda|, and ||r|| is a bound on the error in the eigenvalue, so
! tolerance is on the relative error of the eigenvalue.  This is the default
! in SLEPc and can be changed with EPSSetConvergenceTest.

tol = solver_control%eigen_cntl%tol

if (solver_control%eigen_cntl%maxit <= 0) then
   select case (eigensolver)
   case (SLEPC_POWER)
      maxit = max(2000,100*neq)
   case (SLEPC_SUBSPACE, SLEPC_ARNOLDI, SLEPC_LANCZOS, SLEPC_KRYLOV_SCHUR, &
         SLEPC_GEN_DAVIDSON, SLEPC_JACOBI_DAVIDSON, SLEPC_BLOPEX)
      maxit = max(100,ceiling((2.0_my_real*neq)/ncv))
   case (SLEPC_LAPACK)
      maxit = 0
   case (SLEPC_ARPACK)
      maxit = max(300,ceiling((2.0_my_real*neq)/ncv))
   case default
      ierr = USER_INPUT_ERROR
      call fatal("illegal selection for eigensolver", &
                 intlist=(/eigensolver/))
      stop
   end select
else
   maxit = solver_control%eigen_cntl%maxit
endif

call EPSSetTolerances(eps,tol,maxit,jerr);CHKERRQ(jerr)

! Monitor convergence, if requested

if (io_control%print_error_when == FREQUENTLY) then
   call EPSMonitorSet(eps,EPSMONITORCONVERGED,PETSC_NULL_OBJECT, &
                      PETSC_NULL_FUNCTION,jerr);CHKERRQ(jerr)
elseif (io_control%print_error_when == TOO_MUCH) then
   call EPSMonitorSet(eps,EPSMONITORFIRST,PETSC_NULL_OBJECT, &
                      PETSC_NULL_FUNCTION,jerr);CHKERRQ(jerr)
endif

! Some solvers (e.g., Krylov solvers) can avoid computing the residual by
! using a cheap estimate of the residual norm, but this may sometimes give
! inaccurate results (especially if a spectral transform is being used).
! This option forces them to compute the true residual.

if (solver_control%eigen_cntl%true_residual) then
   call EPSSetTrueResidual(eps,PETSC_TRUE,jerr);CHKERRQ(jerr)
endif

! Set the solver

select case (eigensolver)

case (SLEPC_POWER)
   call EPSSetType(eps,EPSPOWER,jerr);CHKERRQ(jerr)
case (SLEPC_SUBSPACE)
   call EPSSetType(eps,EPSSUBSPACE,jerr);CHKERRQ(jerr)
case (SLEPC_ARNOLDI)
   call EPSSetType(eps,EPSARNOLDI,jerr);CHKERRQ(jerr)
case (SLEPC_LANCZOS)
   call EPSSetType(eps,EPSLANCZOS,jerr);CHKERRQ(jerr)
case (SLEPC_KRYLOV_SCHUR)
   call EPSSetType(eps,EPSKRYLOVSCHUR,jerr);CHKERRQ(jerr)
case (SLEPC_GEN_DAVIDSON)
   call EPSSetType(eps,EPSGD,jerr);CHKERRQ(jerr)
case (SLEPC_JACOBI_DAVIDSON)
   call EPSSetType(eps,EPSJD,jerr);CHKERRQ(jerr)
case (SLEPC_LAPACK)
   call EPSSetType(eps,EPSLAPACK,jerr);CHKERRQ(jerr)
case (SLEPC_ARPACK)
   call EPSSetType(eps,EPSARPACK,jerr);CHKERRQ(jerr)
case (SLEPC_BLOPEX)
   call EPSSetType(eps,EPSBLOPEX,jerr);CHKERRQ(jerr)
case default
   ierr = USER_INPUT_ERROR
   call fatal("illegal selection for eigensolver", &
              intlist=(/eigensolver/))
   stop
end select

! Set harmonic extraction

if (solver_control%eigen_cntl%harmonic_extraction) then
   if (solver_control%lambda0 == -huge(0.0_my_real)) then
      ierr = USER_INPUT_ERROR
      call fatal("harmonic extraction requires lambda0 not be minus infinity")
      stop
   elseif (eigensolver == SLEPC_LANCZOS) then
      ierr = USER_INPUT_ERROR
      call fatal("harmonic extraction is not available for the selected eigensolver")
      stop
   else
      call EPSSetExtraction(eps,EPS_HARMONIC,jerr);CHKERRQ(jerr)
   endif
endif

! Extract the spectral transformation context

   call EPSGetST(eps,st,jerr);CHKERRQ(jerr)

! Set up the spectral transformation, if one is requested.
! The Davidsons and BLOPEX set ST to PRECOND by default.

if (solver_control%eigen_cntl%transformation /= ST_NONE      .and. &
    eigensolver /= SLEPC_GEN_DAVIDSON         .and. &
    eigensolver /= SLEPC_JACOBI_DAVIDSON      .and. &
    eigensolver /= SLEPC_BLOPEX) then

! Set the spectral transform shift

   if (solver_control%eigen_cntl%st_shift /= -huge(0.0_my_real)) then
      shift = solver_control%eigen_cntl%st_shift
      call STSetShift(st,shift,jerr);CHKERRQ(jerr)
   endif

! Set the type of spectral transform

   select case (solver_control%eigen_cntl%transformation)
   case (ST_SHIFT_ORIGIN)
      call STSetType(st,STSHIFT,jerr);CHKERRQ(jerr)
   case (ST_FOLD)
      call STSetType(st,STFOLD,jerr);CHKERRQ(jerr)
   case (ST_SHIFT_INVERT)
      call STSetType(st,STSINVERT,jerr);CHKERRQ(jerr)
   case (ST_CAYLEY)
      call STSetType(st,STCAYLEY,jerr);CHKERRQ(jerr)
      if (solver_control%eigen_cntl%st_antishift /= -huge(0.0_my_real)) then
         call STCayleySetAntishift(st,solver_control%eigen_cntl%st_antishift, &
                                   jerr);CHKERRQ(jerr)
      endif
   case default
      ierr = USER_INPUT_ERROR
      call fatal("illegal value for transformation", &
                 intlist=(/solver_control%eigen_cntl%transformation/))
      stop
   end select

endif ! spectral transformation requested

if (eigensolver == SLEPC_JACOBI_DAVIDSON) then
   call STSetType(st,STPRECOND,jerr);CHKERRQ(jerr)
endif

! Set up the KSP solver and preconditioner, if not using SLEPc's default

if (solver_control%solver /= DEFAULT_SOLVER .or. &
    (eigensolver == SLEPC_GEN_DAVIDSON .and. &
     solver_control%preconditioner /= DEFAULT_PRECONDITION)) then

! Get the KSP context

   call STGetKSP(st,ksp,jerr);CHKERRQ(jerr)

! Set the solver, preconditioner, and options

   if (eigensolver /= SLEPC_GEN_DAVIDSON) then
      call set_petsc_solver(ksp,solver_control)
   endif
   if (solver_control%solver == PETSC_MUMPS_GEN_SOLVER .or. &
       solver_control%solver == PETSC_MUMPS_SPD_SOLVER .or. &
       solver_control%solver == PETSC_SUPERLU_SOLVER) then
      call KSPGetPC(ksp,pc,jerr);CHKERRQ(jerr)
      call PCSetType(pc,PCNONE,jerr);CHKERRQ(jerr)
      call STSetUp(st,jerr);CHKERRQ(jerr)
   endif
   call set_petsc_preconditioner(ksp,solver_control)
   call set_petsc_options(ksp,solver_control,io_control)

endif ! set KSP solver and preconditioner

! Set the initial space to be the previous solution

! copy each nontrivial solution to the initial space

init_space_dim = 0
do i=1,nev
   if (i==1) then
      norm = maxval(abs(linear_system%solution(1:)))
   else
      norm = maxval(abs(linear_system%evecs(:,i-1)))
   endif
   if (norm > near_roundoff) then
      init_space_dim = init_space_dim + 1
      if (i==1) then
         call phamlvec_to_petsc(linear_system%solution(1:), &
                                init_space(init_space_dim), &
                                petsc_stiffness)
      else
         call phamlvec_to_petsc(linear_system%evecs(:,i-1), &
                                init_space(init_space_dim), &
                                petsc_stiffness)
      endif
   endif
end do

! If any solution is nontrivial, set the initial space.  Otherwise SLEPc will
! set a random initial guess.

if (init_space_dim /= 0) then
   call EPSSetInitialSpace(eps,init_space_dim,init_space,jerr);CHKERRQ(jerr)
endif

! I don't expect anyone to set SLEPc options from command line arguments,
! but the SLEPc folks seem to expect this will be called, and I keep running
! into problems if I don't call it.

call EPSSetFromOptions(eps,jerr)

! Solve the eigensystem

call EPSSolve(eps,jerr);CHKERRQ(jerr)

! Check how many eigenpairs converged

call EPSGetConverged(eps,nconv,jerr);CHKERRQ(jerr)
if (nconv < nev) then
   call warning("Number of converged eigenvalues is less than the number requested.", &
                "Remaining eigenpairs are set to zero.", &
                intlist=(/nconv,nev/))
endif
if (nconv > nev) then
   call warning("Number of converged eigenvalues is greater than the number requested.", &
                "Extraneous eigenpairs are dropped, but may have been the desired ones.", &
                intlist=(/nconv,nev/))
   nconv = nev
endif

! Get the eigenvalues

do j=1,nconv
   call EPSGetEigenvalue(eps,j-1,kr,ki,jerr);CHKERRQ(jerr)
   if (abs(ki) > near_roundoff) then
      call warning("Eigenvalue is not real.  Only keeping real part.", &
                   reallist=(/real(ki,my_real)/))
   endif
   if (reverse_M_A) then
      eigenvalues(j) = 1/kr
   else
      eigenvalues(j) = kr
   endif
end do

! Sort so eigenvalues go from smallest to largest

if (nconv > 1) then
   call sort(eigenvalues(1:nconv),nconv,iperm,2,jerr)
else
   iperm(1) = 1
endif

! Space to save the residual norms

if (associated(monitor%eigensolver_l2_resid)) then
   if (size(monitor%eigensolver_l2_resid) /= nev) then
      deallocate(monitor%eigensolver_l2_resid)
      allocate(monitor%eigensolver_l2_resid(nev))
   endif
else
   allocate(monitor%eigensolver_l2_resid(nev))
endif

if (associated(monitor%eigensolver_errbound)) then
   if (size(monitor%eigensolver_errbound) /= nev) then
      deallocate(monitor%eigensolver_errbound)
      allocate(monitor%eigensolver_errbound(nev))
   endif
else
   allocate(monitor%eigensolver_errbound(nev))
endif

! For each converged eigenvalue ...

do i=1,nconv
   j = iperm(i)

! Get the eigenvector, in sorted order

   call EPSGetEigenvector(eps,j-1,xr,xi,jerr);CHKERRQ(jerr)

! Normalize the eigenvector, set orientation so the maximum magnitude is along
! positive real axis, and copy to linear_system solution

   if (i==1) then
      call normalize_and_orient(xr,linear_system%solution(1:),petsc_mass, &
                                solver_control)
   else
      call normalize_and_orient(xr,linear_system%evecs(:,i-1),petsc_mass, &
                                solver_control)
   endif

! Some values for monitoring

   call EPSComputeResidualNorm(eps,j-1,norm,jerr);CHKERRQ(jerr)
   monitor%eigensolver_l2_resid(i) = norm
   call EPSComputeRelativeError(eps,j-1,error,jerr);CHKERRQ(jerr)
   monitor%eigensolver_errbound(i) = error
end do

! Wipe out unconverged eigenpairs

eigenvalues(nconv+1:nev) = 0.0_my_real
if (nconv == 0) then
   linear_system%solution(1:) = 0.0_my_real
endif
if (nev > 1) then
   linear_system%evecs(:,max(1,nconv):nev-1) = 0.0_my_real
endif
monitor%eigensolver_l2_resid(nconv+1:nev) = 0.0_my_real
monitor%eigensolver_errbound(nconv+1:nev) = 0.0_my_real

! Other monitor information

call EPSGetIterationNumber(eps,numit,jerr);CHKERRQ(jerr)
monitor%niter = numit
monitor%nconv = nconv
monitor%ncv = ncv
monitor%maxit = maxit

! Destroy the PETSc form of the matrices and PETSc vectors

call destroy_petsc_linear_system(petsc_stiffness)
call destroy_petsc_linear_system(petsc_mass)
call VecDestroy(xr,jerr);CHKERRQ(jerr)
call VecDestroy(xi,jerr);CHKERRQ(jerr)
do i=1,solver_control%num_eval
   call VecDestroy(init_space(i),jerr);CHKERRQ(jerr)
end do

end subroutine eigen_slepc

!          --------------------
subroutine normalize_and_orient(x,y,petsc_mass,solver_control)
!          --------------------

!----------------------------------------------------
! This routine normalizes the PETSc vector x and orients it such that the
! maximum magnitude is along the positive real axis, and returns it in the
! PHAML vector y
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Vec x
real(my_real), intent(out) :: y(:)
type(petsc_matrix_type), intent(in) :: petsc_mass
type(solver_options), intent(in) :: solver_control

!----------------------------------------------------
! Local variables:

real(my_real) :: norm, xmax, xmin
PetscReal val
PetscScalar mult
integer :: jerr
!----------------------------------------------------
! Begin executable code

! compute the norm of x in which we want to normalize

select case (solver_control%scale_evec)
case (SCALE_LINF)
   norm = Linf_norm(x)
case (SCALE_L2)
   norm = L2_norm(x)
case (SCALE_M)
   norm = M_norm(x,petsc_mass%A)
case default
   ierr = USER_INPUT_ERROR
   call fatal("invalid value for scale_evec")
   stop
end select

! normalize

if (norm /= 1.0_my_real .and. norm /= 0.0_my_real) then
   mult = 1.0_my_real/norm
   call VecScale(x,mult,jerr);CHKERRQ(jerr)
endif

! determine if the maximum magnitude is positive or negative,
! and negate if it is negative

call VecMax(x,PETSC_NULL_INTEGER,val,jerr);CHKERRQ(jerr)
xmax = val
call VecMin(x,PETSC_NULL_INTEGER,val,jerr);CHKERRQ(jerr)
xmin = val

if (abs(xmin) > abs(xmax)) then
   mult = -1.0_my_real
   call VecScale(x,mult,jerr);CHKERRQ(jerr)
endif

! copy the result from PETSc Vec to a PHAML solution vector

call petscvec_to_phaml(x,y,petsc_mass)

end subroutine normalize_and_orient

!        ---------
function Linf_norm(x)
!        ---------

!----------------------------------------------------
! This routine computes the L infinity norm of PETSc vector x
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Vec x
real(my_real) :: Linf_norm
!----------------------------------------------------
! Local variables:

PetscReal norm
integer :: jerr
!----------------------------------------------------
! Begin executable code

call VecNorm(x,NORM_INFINITY,norm,jerr);CHKERRQ(jerr)
Linf_norm = norm

end function Linf_norm

!        -------
function L2_norm(x)
!        -------

!----------------------------------------------------
! This routine computes the L^2 norm of PETSc vector x
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Vec x
real(my_real) :: L2_norm
!----------------------------------------------------
! Local variables:

PetscReal norm
integer :: jerr
!----------------------------------------------------
! Begin executable code

call VecNorm(x,NORM_2,norm,jerr);CHKERRQ(jerr)
L2_norm = norm

end function L2_norm

!        ------
function M_norm(x,M)
!        ------

!----------------------------------------------------
! This routine computes the M norm of PETSc vector x
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Vec x
Mat M
real(my_real) :: M_norm
!----------------------------------------------------
! Local variables:

Vec y
PetscScalar norm
integer :: jerr
!----------------------------------------------------
! Begin executable code

call VecDuplicate(x,y,jerr);CHKERRQ(jerr)
call VecCopy(x,y,jerr);CHKERRQ(jerr)
call MatMult(M,x,y,jerr);CHKERRQ(jerr)
call VecDot(x,y,norm,jerr);CHKERRQ(jerr)
M_norm = sqrt(norm)
call VecDestroy(y,jerr);CHKERRQ(jerr)

end function M_norm

end module slepc_interf
