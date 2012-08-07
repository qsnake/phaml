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

module slepc_init_mod

!----------------------------------------------------
! This module contains routines that initialize and finalize SLEPc.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use petsc_init_mod
!----------------------------------------------------

implicit none
private
public slepc_init, slepc_finalize

!----------------------------------------------------
! The PETSc include files.  Note the use of preprocessor #include instead of
! the Fortran include statement, because the include files contain
! preprocessor directives.

#include "include/finclude/petsc.h"

!----------------------------------------------------
! The following parameters are defined:

!----------------------------------------------------

!----------------------------------------------------
! The following types are defined:

!----------------------------------------------------

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ----------
subroutine slepc_init(slaves_communicator)
!          ----------

!----------------------------------------------------
! This routine initializes SLEPc
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: slaves_communicator
!----------------------------------------------------
! Local variables:

integer :: jerr
!----------------------------------------------------
! Begin executable code

PETSC_COMM_WORLD = slaves_communicator
call SlepcInitialize(PETSC_NULL_CHARACTER,jerr)
petsc_init_called = .true.
petsc_final_called = .false.

end subroutine slepc_init

!          --------------
subroutine slepc_finalize
!          --------------

!----------------------------------------------------
! This routine finalizes SLEPc
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

integer :: jerr
!----------------------------------------------------
! Begin executable code

call SlepcFinalize(jerr)
petsc_init_called = .false.
petsc_final_called = .true.

end subroutine slepc_finalize

end module slepc_init_mod
