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
! This version illustrates the subroutines for a single PDE, thus all the
! array arguments have dimension 1.  For a coupled system of PDEs, see
! examples/system.
!----------------------------------------------------

module phaml_user_mod
use global
real(my_real) :: prob_param = 6.0_my_real
end module phaml_user_mod

!          --------
subroutine pdecoefs(x,y,z,cxx,cyy,czz,cxy,cxz,cyz,cx,cy,cz,c,rs)
!          --------

!----------------------------------------------------
! This subroutine returns the coefficient and right hand side of the PDE
! at the point (x,y)
!
! The PDE is
!
!    -( cxx(x,y) * u  )  -( cyy(x,y) * u  ) + c(x,y) * u = rs(x,y)
!                   x  x                y  y
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
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
real(my_real), intent(out) :: cxx(:,:),cyy(:,:),czz(:,:),cxy(:,:),cxz(:,:), &
                              cyz(:,:),cx(:,:),cy(:,:),cz(:,:),c(:,:),rs(:)
!----------------------------------------------------

!----------------------------------------------------
! Local variables

real(my_real) :: r
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

cxy=0; cxz=0; cyz=0; cx=0; cy=0; cz=0; c=0

cxx = 1.0_my_real
cyy = 1.0_my_real
czz = 1.0_my_real
r = prob_param*sqrt(x**2+y**2)
rs = z**4*prob_param**2*(sin(r) -cos(r)/r)

end subroutine pdecoefs

!          ------
subroutine bconds(x,y,z,bmark,itype,c,rs)
!          ------

!----------------------------------------------------
! This subroutine returns the boundary conditions at the point (x,y).
!
! Each boundary condition is either
!
!    u = rs(x,y) or  u  + c(x,y)*u = rs(x,y)
!                     n
!
! itype indicates which type of condition applies to each point, using
! symbolic constants from module phaml.  It is DIRICHLET for the first
! condition, NATURAL for the second condition with c==0, and MIXED for
! the second condition with c/=0.
!
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
!----------------------------------------------------

!----------------------------------------------------
! Non-module procedures used are:

interface

   function trues(x,y,z,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y,z
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
rs(1) = trues(x,y,z,1,1)

end subroutine bconds

!        ------
function iconds(x,y,z,comp,eigen)
!        ------

!----------------------------------------------------
! This routine returns the initial condition for a time dependent problem.
! It can also be used for the initial guess for a nonlinear problem, or
! systems of equations solved by sucessive substitution.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
! For problems where there are no initial conditions, it is a dummy.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
integer, intent(in) :: comp,eigen
real(my_real) :: iconds

!----------------------------------------------------
! Begin executable code

iconds = 0.0_my_real

end function iconds

!        -----
function trues(x,y,z,comp,eigen) ! real (my_real)
!        -----

!----------------------------------------------------
! This is the true solution of the differential equation, if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
integer, intent(in) :: comp,eigen
real (my_real) :: trues
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

trues = z**4*sin(prob_param*sqrt(x**2+y**2))

end function trues

!        ------
function truexs(x,y,z,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the x derivative of the true solution of the differential
! equation, if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
integer, intent(in) :: comp,eigen
real (my_real) :: truexs
!----------------------------------------------------

real(my_real) :: r
!----------------------------------------------------
! Begin executable code

r = sqrt(x**2+y**2)
truexs = z**4*cos(prob_param*r)*prob_param*x/r

end function truexs

!        ------
function trueys(x,y,z,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the y derivative of the true solution of the differential
! equation, if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
integer, intent(in) :: comp,eigen
real (my_real) :: trueys
!----------------------------------------------------

real(my_real) :: r
!----------------------------------------------------
! Begin executable code

r = sqrt(x**2+y**2)
trueys = z**4*cos(prob_param*r)*prob_param*y/r

end function trueys

!        ------
function truezs(x,y,z,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the z derivative of the true solution of the differential
! equation, if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y,z
integer, intent(in) :: comp,eigen
real (my_real) :: truezs
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

truezs = 4*z**3*sin(prob_param*sqrt(x**2+y**2))

end function truezs

!          --------------
subroutine update_usermod(phaml_solution)
!          --------------

use phaml
type(phaml_solution_type), intent(in) :: phaml_solution

! Dummy version.  See examples/elliptic/usermod.f90 for a functional example.

end subroutine update_usermod

!        ---------------------
function phaml_integral_kernel(kernel,x,y,z)
!        ---------------------

use phaml
integer, intent(in) :: kernel
real(my_real), intent(in) :: x,y,z
real(my_real) :: phaml_integral_kernel

! Identity function

phaml_integral_kernel = 1.0

end function phaml_integral_kernel

!        ----------
function regularity(x,y,z)
!        ----------

use phaml
real(my_real), intent(in) :: x(*),y(*),z(*)
real(my_real) :: regularity

! Dummy version, assume infinitely differentiable everywhere.

regularity = huge(0.0_my_real)

end function regularity
