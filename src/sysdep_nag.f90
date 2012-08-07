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

module sysdep

!----------------------------------------------------
! This module contains routines that are not standard Fortran 90 and may
! need to be changed for different compilers.  This version is for NAG Fortran.
!----------------------------------------------------

!----------------------------------------------------

implicit none
private
public my_system, my_flush

!----------------------------------------------------

contains

!          ---------
subroutine my_system(cmd)
!          ---------

!----------------------------------------------------
! This routine executes a system command as if from the system command line.
!----------------------------------------------------

use f90_unix_proc, only: system
!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: cmd
!----------------------------------------------------
! Begin executable code

call system(cmd)

end subroutine my_system

!          --------
subroutine my_flush(unit)
!          --------

!----------------------------------------------------
! This routine flushes the io buffer of the given unit.
!
!----------------------------------------------------

use f90_unix_io, only:flush
!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: unit
!----------------------------------------------------
! Begin executable code

call flush(unit)

end subroutine my_flush

end module sysdep
