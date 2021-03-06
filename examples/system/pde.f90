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
! This file contains the user supplied external subroutines that define
! the PDE(s) to be solved, and other required external subroutines.
!   pdecoefs bconds boundary_point boundary_npiece boundary_param iconds trues
!   truexs trueys update_usermod phaml_integral_kernel
!
! This version illustrates a system of elliptic PDEs solved simultaneously.
!----------------------------------------------------

!          --------
subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
!          --------

!----------------------------------------------------
! This subroutine returns the coefficient and right hand side of the PDEs
! at the points (x,y)
!
! PDEs are
!
!    -( cxx(x,y) * u  )  -( cyy(x,y) * u  ) + c(x,y) * u = rs(x,y)
!                   x  x                y  y
!
! where cxx, cyy and c are matrices in which the (i,j)th entry is the
! coefficient of the jth solution in the ith equation, u is a vector of
! the solutions, and rs is a vector of the right hand sides.
!
! For eigenvalue problems, the right hand side is lambda * u * rs(x,y)
!
! cxy, cx and cy are not currently used and should not be set.  They are
! for future expansion.
!
! NOTE: BE CAREFUL TO GET THE SIGNS RIGHT
! e.g. cxx=cyy=1 means rs=-(uxx+uyy)
!
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)
!----------------------------------------------------

!----------------------------------------------------
! Local variables

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! The particular equations for this example are
!
! -uxx -uyy + v = f1
! -vxx -vyy + u = f2

cxx(1,1) = 1.0_my_real;  cxx(1,2) = 0.0_my_real
cxx(2,1) = 0.0_my_real;  cxx(2,2) = 1.0_my_real

cyy(1,1) = 1.0_my_real;  cyy(1,2) = 0.0_my_real
cyy(2,1) = 0.0_my_real;  cyy(2,2) = 1.0_my_real

  c(1,1) = 0.0_my_real;    c(1,2) = 1.0_my_real
  c(2,1) = 1.0_my_real;    c(2,2) = 0.0_my_real


rs(1) = -(2.0_my_real*exp(x-y) - (x+y)**4/8.0_my_real)
rs(2) = -(3.0_my_real*(x+y)**2 - exp(x-y))

cxy=0; cx=0; cy=0

end subroutine pdecoefs

!          -----
subroutine bconds(x,y,bmark,itype,c,rs)
!          -----

!----------------------------------------------------
! This subroutine returns the boundary conditions at the points (x,y).
!
! Each boundary condition is either
!
!    u = rs(x,y) or  u  + c(x,y)*u = rs(x,y)
!                     n
!
! c is an array in which the (i,j)th component is the coefficient of the
! jth solution in the ith boundary condition.  rs is a vector in which the
! ith component is the right side of the ith boundary condition.
!
! itype indicates which type of condition applies to each point, using
! symbolic constants from module phaml.  It is DIRICHLET for the first
! condition, NATURAL for the second condition with c==0, and MIXED for
! the second condition with c/=0.  The ith component is the type of the
! ith boundary condition.
!
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
!----------------------------------------------------

!----------------------------------------------------
! Non-module procedures used are:

interface

   function trues(x,y,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trues
   end function trues

end interface

!----------------------------------------------------
! Local variables

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! Dirichlet boundary conditions

itype = DIRICHLET
c = 0.0_my_real
rs(1) = trues(x,y,1,1)
rs(2) = trues(x,y,2,1)

end subroutine bconds

!        ------
function iconds(x,y,comp,eigen)
!        ------

!----------------------------------------------------
! This routine returns the initial condition of the comp'th solution
! for a time dependent problem.
! It can also be used for the initial guess for a nonlinear problem, or
! systems of equations solved by sucessive substitution.
! For problems where there are no initial conditions, it is a dummy.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real(my_real) :: iconds

!----------------------------------------------------
! Begin executable code

iconds = 0.0_my_real

end function iconds

!        -----
function trues(x,y,comp,eigen) ! real (my_real)
!        -----

!----------------------------------------------------
! This is the comp'th true solution of the differential equation, if known.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trues
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

select case (comp)

   case (1)
      trues = exp(x-y)

   case (2)
      trues = (x+y)**4/8.0_my_real

   case default
      print *,"ERROR -- bad value ",comp," for comp in trues"
      stop

end select

end function trues

!        ------
function truexs(x,y,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the x derivative of the comp'th true solution of the differential
! equation, if known.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: truexs
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

select case (comp)

   case (1)
      truexs = exp(x-y)

   case (2)
      truexs = (x+y)**3/2.0_my_real

   case default
      print *,"ERROR -- bad value ",comp," for comp in truexs"
      stop

end select

end function truexs

!        ------
function trueys(x,y,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the y derivative of the comp'th true solution of the differential
! equation, if known.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trueys
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

select case (comp)

   case (1)
      trueys = -exp(x-y)

   case (2)
      trueys = (x+y)**3/2.0_my_real

   case default
      print *,"ERROR -- bad value ",comp," for comp in trueys"
      stop

end select

end function trueys

!          --------------
subroutine boundary_point(ipiece,s,x,y)
!          --------------

!----------------------------------------------------
! This routine defines the boundary of the domain.  It returns the point
! (x,y) with parameter s on piece ipiece.
! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
! data files.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: ipiece
real(my_real), intent(in) :: s
real(my_real), intent(out) :: x,y

!----------------------------------------------------
! Begin executable code

! Dummy version

x = 0.0_my_real
y = 0.0_my_real

end subroutine boundary_point

!        ---------------
function boundary_npiece(hole)
!        ---------------

!----------------------------------------------------
! This routine gives the number of pieces in the boundary definition.
! If boundary_npiece <= 0 the domain is given by triangle data files.
!----------------------------------------------------

!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: hole
integer :: boundary_npiece
!----------------------------------------------------
! Begin executable code

boundary_npiece = 0

end function boundary_npiece

!          --------------
subroutine boundary_param(start,finish)
!          --------------

!----------------------------------------------------
! This routine gives the range of parameter values for each piece of the
! boundary.
! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
! data files.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(out) :: start(:), finish(:)
!----------------------------------------------------
! Begin executable code

! Dummy version

start = 0.0_my_real
finish = 0.0_my_real

end subroutine boundary_param

!          --------------
subroutine update_usermod(phaml_solution)
!          --------------

use phaml
type(phaml_solution_type), intent(in) :: phaml_solution

! Dummy version

end subroutine update_usermod

!        ---------------------
function phaml_integral_kernel(kernel,x,y)
!        ---------------------

use phaml
integer, intent(in) :: kernel
real(my_real), intent(in) :: x,y
real(my_real) :: phaml_integral_kernel

! Identity function

phaml_integral_kernel = 1.0

end function phaml_integral_kernel

!        ----------
function regularity(x,y)
!        ----------

use phaml
real(my_real), intent(in) :: x(3),y(3)
real(my_real) :: regularity

! Dummy version, assume infinitely differentiable everywhere.

regularity = huge(0.0_my_real)

end function regularity
