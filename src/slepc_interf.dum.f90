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
! This module contains dummy routines corresponding to those in
! slepc_interf.F90 to satisfy references when SLEPc is not used.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use linsystype_mod
use message_passing
use gridtype_mod
!----------------------------------------------------

implicit none
private
public eigen_slepc

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          -----------
subroutine eigen_slepc(linear_system,procs,solver_control,io_control, &
                       still_sequential,eigenvalues,monitor)
!          -----------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
real(my_real), intent(out) :: eigenvalues(:)
type(eigen_monitor), intent(inout) :: monitor

end subroutine eigen_slepc

end module slepc_interf
