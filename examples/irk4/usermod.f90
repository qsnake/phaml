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

!----------------------------------------------------
! This file contains module phaml_user_mod and subroutine update_usermod
!----------------------------------------------------

module phaml_user_mod

!----------------------------------------------------
! This module contains user global data.
!
! This version is for the Implicite Runge Kutta 4 solution of a parabolic
! equation.  It contains the variables that hold the current time and time step,
! and the Butcher coefficients for the IRK4 formula.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
!----------------------------------------------------

implicit none

!----------------------------------------------------

! The time step

real(my_real) :: t, deltat

! The Butcher variables, given as rational numbers in the Wang & Mavriplis paper

real(my_real), parameter :: butcher(6,6) = reshape( (/ &
0.0_my_real, 0.25_my_real, 0.137776_my_real, & ! a_11, a_21, a_31
0.144636866026982_my_real, .09825878328356477_my_real, 0.157916295161671_my_real, & ! a_41, a_51, a_61
0.0_my_real, 0.25_my_real, -.055776_my_real, & ! a_12, a_22, a_32
-0.223931907613345_my_real, -0.59154424281967_my_real, 0.0_my_real, & ! a_42, a_52, a_62
0.0_my_real, 0.0_my_real, 0.25_my_real, & ! a_13, a_23, a_33
0.449295041586363_my_real, 0.8101210538283_my_real, 0.186758940524001_my_real, & ! a_43, a_53, a_63
0.0_my_real, 0.0_my_real, 0.0_my_real, & ! a_14, a_24, a_34
0.25_my_real, 0.283164405707806_my_real, 0.680565295309335_my_real, & ! a_44, a_54, a_64
0.0_my_real, 0.0_my_real, 0.0_my_real, & ! a_15, a_25, a_35
0.0_my_real, 0.25_my_real, -0.275240530995007_my_real, & ! a_45, a_55, a_65
0.0_my_real, 0.0_my_real, 0.0_my_real, & ! a_16, a_26, a_36
0.0_my_real, 0.0_my_real, 0.25_my_real /), & ! a_46, a_56, a_66
(/6,6/))

end module phaml_user_mod

!          --------------
subroutine update_usermod(phaml_solution) 
!          --------------

!----------------------------------------------------
! This routine updates the module variables on the slave processes by
! sending them from the master process
!----------------------------------------------------

use phaml
use phaml_user_mod
   
!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(in) :: phaml_solution

!----------------------------------------------------
! Local variables:

! Declare these arrays big enough to hold the variables to be sent

integer :: iparam(1)
real(my_real) :: rparam(2)

!----------------------------------------------------
! Begin executable code

! Copy the module variables into the arrays, putting integer variables
! into iparam and real variables into rparam.

   rparam(1) = deltat
   rparam(2) = t

! Call the routine that performs the actual exchange.  Don't change this line.

   call master_to_slaves(phaml_solution,iparam,rparam)

! Copy the arrays into the module variables, using the same correspondence
! between module variable and array index as was used above.
   
   deltat = rparam(1)
   t = rparam(2)
   
end subroutine update_usermod
