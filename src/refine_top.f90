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

module refine_top

!----------------------------------------------------
! This module contains the top level refinement routines
!
! communication tags in this module are of the form 34xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use grid_util
use linsystype_mod
use message_passing
use refine_uniform_mod
use refine_adapt_mod
use refine_elements
use hp_strategies
use hash_mod
use zoltan_interf
!----------------------------------------------------

implicit none
private
public refine

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ------
subroutine refine(grid,procs,refine_control,solver_control,io_control, &
                  lb,still_sequential,init_nvert,init_nelem,init_dof,loop, &
                  partition_method,balance_what,predictive,stalled,no_time)
!          ------

!----------------------------------------------------
! Top level refinement routine.  It mainly determines which refinement
! routine to call.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: init_nvert, init_nelem, init_dof, loop, &
                       partition_method,balance_what
logical, intent(in) :: still_sequential, predictive
logical, intent(out) :: stalled
logical, intent(in), optional :: no_time

!----------------------------------------------------
! Local variables:

logical :: timeit
!----------------------------------------------------
! Begin executable code

! start timing the refinement process

if (present(no_time)) then
   timeit = .not. no_time
else
   timeit = .true.
endif
if (timeit) then
   call reset_watch(prefine)
   call start_watch((/prefine,trefine/))
endif

! The REFSOLN_EDGE and REFSOLN_ELEM hp-adaptive strategies have their own
! refine routine

if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_REFSOLN_EDGE .or. &
     refine_control%hp_strategy == HP_REFSOLN_ELEM)) then
   call refine_refsoln(grid,procs,refine_control,solver_control, &
                       io_control,still_sequential)
   stalled = .false.

! NLP also has it's own routine

elseif (refine_control%reftype == HP_ADAPTIVE .and. &
        refine_control%hp_strategy == HP_NLP) then
   call refine_nlp(grid,procs,refine_control,solver_control,still_sequential)
   stalled = .false.

! uniform refinement

elseif (refine_control%reftype == P_UNIFORM) then
   call refine_uniform_p(grid,refine_control)
   stalled = .false.
elseif (refine_control%reftype == H_UNIFORM) then
   call refine_uniform_h(grid,refine_control,solver_control)
   stalled = .false.

! all other adaptive refinements

else
   call refine_adaptive(grid,procs,refine_control,solver_control,io_control, &
                        lb,still_sequential,init_nvert,init_nelem,init_dof, &
                        loop,partition_method,balance_what,predictive, &
                        stalled,no_time)

endif

! error indicators should be recomputed and say we need to resend graphics data

grid%errind_up2date = .false.
grid_changed = .true.

! stop timing the refinement process

if (timeit) call stop_watch((/prefine,trefine/))

end subroutine refine

end module refine_top
