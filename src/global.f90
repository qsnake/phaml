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

module global

!----------------------------------------------------
! This module contains some global variables
!----------------------------------------------------

use stopwatch
implicit none

!----------------------------------------------------
! Parameters defined are:

logical :: new_comm = .false.  ! TEMP090128 testing new forms of communication

! kind numbers

! an easy way to switch between single and double precision, or any
! other precision supported

! Note: only single and double precision are supported (not quad or other
! precisions) because BLAS and LAPACK only have single and double precision.

integer, parameter :: my_real = kind(1.0d0)

real(my_real), save, allocatable :: charge_density(:) ! TEMP110525

! save on memory by using a compiler specific small logical kind.
! use default kind when in doubt

integer, parameter :: small_logical = kind(.false.)
!integer, parameter :: small_logical = 1   ! Works on many compilers

! quadruple precision real, if available, double precision otherwise.  This
! is used when computing barycentric coordinates while trying to find a
! point in a very small triangle

!integer, parameter :: quad_real = kind(1.0q0)
integer, parameter :: quad_real = kind(1.0d0)

! maximum length of file names

integer, parameter :: FN_LEN = 256

! directory to use as scratch space to write messages to the graphics
! process when running a sequential program.  Be sure the length is
! long enough for the string, and that you end it with a / (or for
! non-unix, whatever fully specifies a path)

character(len=5), parameter :: scratch = "/tmp/"

! version number

character(len=6), parameter :: version_number = "1.11.0"

! termination codes

integer, parameter :: NO_ERROR               =  0, &
                      MAX_VERTICES_ACHIEVED  =  1, &
                      MAX_ELEMENTS_ACHIEVED  =  3, &
                      MAX_EQUATIONS_ACHIEVED =  4, &
                      MAX_LOOP_ACHIEVED      =  5, &
                      TOLERANCE_ACHIEVED     =  6, &
                      DONE_BALANCE_ONLY      =  7, &
                      DONE_SOLVE_ONLY        =  8, &
                      DONE_REFINE_ONLY       =  9, &
                      ENERGY_ERREST_ACHIEVED = 10, &
                      LINF_ERREST_ACHIEVED   = 12, &
                      L2_ERREST_ACHIEVED     = 14, &
                      MAX_LEV_ACHIEVED       = 15, &
                      MAX_DEG_ACHIEVED       = 16, &
                      ALLOC_FAILED           = -1, &
                      USER_INPUT_ERROR       = -2, &
                      PHAML_INTERNAL_ERROR   = -3, &
                      REFINEMENT_STALLED     = -4, &
                      UNCLASSIFIED_ERROR     = -5, &
                      MULTIPLE_ERRORS        = -6

! values for [draw|print]_[grid|error|time]_who

! MASTER is also used in tests for my_proc==MASTER so it must agree
! with the master's process number, usually 0

integer, parameter :: MASTER     = 0, &
                      SLAVES     = 2, &
                      EVERYONE   = 3, &
                      MASTER_ALL = 4, &
                      NO_ONE     = 5

! in addition, MASTER, SLAVES and GRAPHICS identify the proc type

integer, parameter :: GRAPHICS = 6

! values for [draw|print]_[grid|error|time]_when

integer, parameter :: NEVER          =  5, &
                      FINAL          =  6, &
                      PHASES         =  7, &
                      FREQUENTLY     =  8, &
                      TOO_MUCH       =  9, &
                      LAST           = 10, &
                      LAST_AND_FINAL = 11

! values for print_[error|errest]_what

integer, parameter :: ENERGY_ERR              = 21, &
                      LINF_ERR                = 22, &
                      L2_ERR                  = 23, &
                      ENERGY_LINF_ERR         = 24, &
                      ENERGY_L2_ERR           = 25, &
                      LINF_L2_ERR             = 26, &
                      ENERGY_LINF_L2_ERR      = 27, &
                      ENERGY_ERREST           = 31, &
                      LINF_ERREST             = 32, &
                      L2_ERREST               = 33, &
                      ENERGY_LINF_ERREST      = 34, &
                      ENERGY_L2_ERREST        = 35, &
                      LINF_L2_ERREST          = 36, &
                      ENERGY_LINF_L2_ERREST   = 37

! values for errtype

integer, parameter :: ABSOLUTE_ERROR = 1, &
                      RELATIVE_ERROR = 2

! control over printing warning messages

logical :: warn = .true.

! data separators for messages to graphics processes

integer, parameter :: END_OF_ELEMENTS  = -101, &
                      END_OF_FACES     = -102, &
                      END_OF_EDGES     = -103, &
                      END_OF_VERTICES  = -104

! values for selecting which clocks

integer, parameter :: CLOCK_W  = 1, &
                      CLOCK_C  = 2, &
                      CLOCK_CW = 3, &
                      CLOCK_WC = CLOCK_CW

! methods for error indicator

integer, parameter :: HIERARCHICAL_COEFFICIENT = 1, &
                      TRUE_DIFF                = 2, &
                      LOCAL_PROBLEM_H          = 4, &
                      LOCAL_PROBLEM_P          = 5, &
                      INITIAL_CONDITION        = 6, &
                      EXPLICIT_ERRIND          = 7, &
                      EQUILIBRATED_RESIDUAL    = 8

! selection of refinement method

integer, parameter :: H_UNIFORM   = 1, &
                      H_ADAPTIVE  = 2, &
                      P_UNIFORM   = 3, &
                      P_ADAPTIVE  = 4, &
                      HP_ADAPTIVE = 5

! selection of hp adaptive strategy

integer, parameter :: HP_BIGGER_ERRIND  =  1, &
                      HP_APRIORI        =  2, &
                      HP_PRIOR2P_E      =  3, &
                      HP_PRIOR2P_H1     =  4, &
                      HP_T3S            =  5, &
                      HP_ALTERNATE      =  6, &
                      HP_TYPEPARAM      =  7, &
                      HP_COEF_DECAY     =  8, &
                      HP_COEF_ROOT      =  9, &
                      HP_SMOOTH_PRED    = 10, &
                      HP_NEXT3P         = 11, &
                      HP_REFSOLN_EDGE   = 12, &
                      HP_REFSOLN_ELEM   = 13, &
                      HP_NLP            = 14, &
                      HP_STEEPEST_SLOPE = 15

! selection of rule for edge/face degree

integer, parameter :: MINIMUM_RULE = 1, &
                      MAXIMUM_RULE = 2

! choices for type of PDE

integer, parameter :: ELLIPTIC = 1, &
                      EIGENVALUE = 2

! choices for solver

integer, parameter :: MG_SOLVER               = 1, &
                      CG_SOLVER               = 2, &
                      GMRES_SOLVER            = 3, &
                      LAPACK_INDEFINITE_SOLVER= 4, &
                      LAPACK_SPD_SOLVER       = 5, &
                      HYPRE_BOOMERAMG_SOLVER  = 10, &
                      HYPRE_PCG_SOLVER        = 11, &
                      HYPRE_GMRES_SOLVER      = 12, &
                      PETSC_RICHARDSON_SOLVER = 13, &
                      PETSC_CHEBYCHEV_SOLVER  = 14, &
                      PETSC_CG_SOLVER         = 15, &
                      PETSC_GMRES_SOLVER      = 16, &
                      PETSC_TCQMR_SOLVER      = 17, &
                      PETSC_BCGS_SOLVER       = 18, &
                      PETSC_CGS_SOLVER        = 19, &
                      PETSC_TFQMR_SOLVER      = 20, &
                      PETSC_CR_SOLVER         = 21, &
                      PETSC_LSQR_SOLVER       = 22, &
                      PETSC_BICG_SOLVER       = 23, &
                      PETSC_MUMPS_GEN_SOLVER  = 24, &
                      PETSC_MUMPS_SPD_SOLVER  = 25, &
                      PETSC_SUPERLU_SOLVER    = 26, &
                      DEFAULT_SOLVER          = 27

! choices for preconditioner

integer, parameter :: NO_PRECONDITION              = 1, &
                      MG_PRECONDITION              = 2, &
                      COARSE_GRID_PRECONDITION     = 4, &
                      TEST_PRECONDITION            = 5, &
                      FUDOP_DD_PRECONDITION        = 6, &
                      PETSC_JACOBI_PRECONDITION    = 11, &
                      PETSC_BJACOBI_PRECONDITION   = 12, &
                      PETSC_SOR_PRECONDITION       = 13, &
                      PETSC_EISENSTAT_PRECONDITION = 14, &
                      PETSC_ICC_PRECONDITION       = 15, &
                      PETSC_ILU_PRECONDITION       = 16, &
                      PETSC_ASM_PRECONDITION       = 17, &
                      HYPRE_BOOMERAMG_PRECONDITION = 21, &
                      HYPRE_PARASAILS_PRECONDITION = 22, &
                      HYPRE_DS_PRECONDITION        = 23, &
                      DEFAULT_PRECONDITION         = 24

! choices for eigensolver

integer, parameter :: SLEPC_POWER             = 1, &
                      SLEPC_SUBSPACE          = 2, &
                      SLEPC_ARNOLDI           = 3, &
                      SLEPC_LANCZOS           = 4, &
                      SLEPC_KRYLOV_SCHUR      = 5, &
                      SLEPC_GEN_DAVIDSON      = 6, &
                      SLEPC_JACOBI_DAVIDSON   = 7, &
                      SLEPC_LAPACK            = 8, &
                      SLEPC_ARPACK            = 9, &
                      SLEPC_BLOPEX            = 10

! multigrid choices

real(my_real), parameter :: MG_NO_TOL = -1.0_my_real, &
                            MG_ERREST_TOL = -2.0_my_real, &
                            KRYLOV_ERREST_TOL = MG_ERREST_TOL

integer, parameter :: MGCOMM_NONE         = 1, &
                      MGCOMM_FUDOP        = 2, &
                      MGCOMM_CONVENTIONAL = 3

! what to balance during load balancing

integer, parameter :: BALANCE_NONE     = 1, &
                      BALANCE_ELEMENTS = 2, &
                      BALANCE_VERTICES = 3, &
                      BALANCE_EQUATIONS= 4

! refinement termination choices

integer, parameter :: DOUBLE_NVERT        = 1, &
                      DOUBLE_NELEM        = 3, &
                      DOUBLE_NEQ          = 5, &
                      HALVE_ERREST        = 7, &
                      KEEP_NVERT          = 8, &
                      KEEP_NELEM          = 10, &
                      KEEP_NEQ            = 12, &
                      KEEP_ERREST         = 14, &
                      ONE_REF             = 15, &
                      ONE_REF_HALF_ERRIND = 16

! equation/vertex/edge/face types

integer, parameter :: INTERIOR  = 1, &
                      DIRICHLET = -2, &  ! must be negative
                      NATURAL   = 3, &
                      MIXED     = 4, &
                      PERIODIC  = 5, &
                      PERIODIC_MASTER = 6, &
                      PERIODIC_SLAVE = 7, &
                      PERIODIC_MASTER_DIR = 8, &
                      PERIODIC_SLAVE_DIR = 9, &
                      PERIODIC_MASTER_NAT = 10, &
                      PERIODIC_SLAVE_NAT = 11, &
                      PERIODIC_MASTER_MIX = 12, &
                      PERIODIC_SLAVE_MIX = 13

! methods for partitioning algorithm

integer, parameter :: RTK            = 1, &
                      ZOLTAN_RCB     = 2, &
                      ZOLTAN_OCT     = 3, &
                      ZOLTAN_METIS   = 4, &
                      ZOLTAN_REFTREE = 5, &
                      ZOLTAN_RIB     = 6, &
                      ZOLTAN_HSFC    = 7, &
                      ZOLTAN_FILE    = 8

! values for task

integer, parameter :: BALANCE_REFINE_SOLVE = 1, &
                      SET_INITIAL          = 2, &
                      BALANCE_ONLY         = 3, &
                      REFINE_ONLY          = 4, &
                      SOLVE_ONLY           = 5

! forms of parallelism

integer, parameter :: SEQUENTIAL = 1, &
                      MPI1       = 2, &
                      MPI2       = 3

! choices for scaling eigenvectors

integer, parameter :: SCALE_LINF = 1, &
                      SCALE_L2   = 2, &
                      SCALE_M    = 3

! choices for the spectral transformation

integer, parameter :: ST_NONE         = 1, &
                      ST_SHIFT_ORIGIN = 2, &
                      ST_FOLD         = 3, &
                      ST_SHIFT_INVERT = 4, &
                      ST_CAYLEY       = 5

! message tag for error handling

integer, parameter :: ERR_TAG = 16000

! types of elements

integer, parameter :: TRIANGULAR_ELEMENT  = 1, &
                      TETRAHEDRAL_ELEMENT = 2

! format for saving a grid

integer, parameter :: GRIDFILE_POLY      = 1, &
                      GRIDFILE_POLY_SOLN = 2, &
                      GRIDFILE_MSH       = 3, &
                      GRIDFILE_MSH_SOLN  = 4

!----------------------------------------------------
! Types defined are :

type io_options
   integer :: print_grid_when,   print_grid_who,  &
              print_linsys_when, print_linsys_who,&
              print_solver_when, print_solver_who, &
              print_error_when,  print_error_who, &
              print_error_what,  print_errest_what, &
              print_time_when,   print_time_who,  &
              print_eval_when,   print_eval_who, &
              draw_grid_when
   logical :: pause_after_draw
end type io_options

!----------------------------------------------------
! Variables defined are:

integer, save :: ierr = 0               ! an error code flag
logical, save :: grid_changed = .false. ! for determining when to send graphics
integer, save :: my_pde_id              ! id for multiple pdes
integer, save :: outunit, errunit       ! I/O unit numbers
integer, save :: global_element_kind    ! kind of elements, set during create
logical, save :: crank_nicholson=.false.! undocumented feature for using Crank-Nicholson
                                        ! for time-dependent Schroedinger equation
logical, save :: charged_particles=.false. ! undocumented feature for equations
                                           ! with charged particles TEMP080114
logical, save :: save_convergence = .false. ! undocumented feature; see
                                            ! linsys_io.f90
integer, save :: convfileunit = 21      ! unit to attach to convergence file
integer, save :: tds_scale = 0          ! undocumented feature for special
                                        ! scaling of graphics for time
                                        ! dependent Schroedinger
integer, save, allocatable :: tags(:)   ! undocumented feature to make the tags
                                        ! of the current element available to
                                        ! pde_coefs

! Undocumented feature to make the grid and current element lid available to
! pde_coefs to, for example, access tags.  But don't tell anyone!
! pde_coefs need to use both global and grid_util.

integer, save :: secret_elem

! TEMP071217 to keep desired initial grid with battery example.  The battery
!            main program changes it to .true.
logical :: dont_smooth = .false.

! TEMP081103 also for the battery problem, say to move the points slightly
!            inside and outside the element in local_problem_p
logical, save :: battery_cludge = .false.

! watches
type (watchtype), save :: &
! total times
   ttotal, trefine, trecon, tpartition, tdistribute, tassemble, tsolve, &
! current phase times
   ptotal, prefine, precon, ppartition, pdistribute, passemble, psolve, &
! time spent in communication, total
   ctrecon, ctpartition, ctdistribute, ctassemble, ctsolve, &
! time spent in communication, phase
   cprecon, cppartition, cpdistribute, cpassemble, cpsolve

type (watchgroup), save :: all_watches

contains

!          --------
subroutine my_pause(sec)
!          --------

!----------------------------------------------------
! This routine creates a busy wait for sec seconds before returning
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real, intent(in) :: sec
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

type (watchtype) :: w
real :: time_read
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

call create_watch(w)
call start_watch(w)
time_read = 0.
do while (time_read < sec)
   call read_watch(time_read,w,'wall')
end do
call destroy_watch(w)

end subroutine my_pause

!        -------
function get_lun()
!        -------

!----------------------------------------------------
! This routine returns an available logical unit number for i/o
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: get_lun
!----------------------------------------------------
! Local variables:

logical :: exists, opened
!----------------------------------------------------
! Begin executable code

get_lun = 11
do
   inquire(unit=get_lun,exist=exists,opened=opened)
   if (exists .and. .not. opened) exit
   get_lun = get_lun + 1
end do

end function get_lun

end module global
