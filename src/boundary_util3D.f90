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

module boundary_util

!----------------------------------------------------
! This module contains utilities for dealing with the 3D surface boundary.
!
! communication tags in this module are of the form xxx (currently there is
! no communication in this module)
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use boundtype_mod
use gridtype_mod
use message_passing, only : fatal
!----------------------------------------------------

implicit none
private
public define_boundary_vertex, deallocate_boundary, coord_is_on_surface

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

! these are needed to evaluate the function whose zero gives the
! point at which to bisect the refinement edge
type(point) :: coord1, coord2
type(bvert_type), pointer :: bv1, bv2
type(grid_type), pointer :: hold_grid
integer :: hold_line

! these are needed to evaluate the function for a point on an ellipse
type(point) :: e_C, e_u, e_v, e_P1, e_P2
real(my_real) :: e_A, e_alpha, e_beta

! TEMP120406 until I decide if this is an option
logical :: gmsh_sphere_hack = .true.
!----------------------------------------------------

contains

!          ----------------------
subroutine define_boundary_vertex(grid,refinement_edge,vert_mid)
!          ----------------------

!----------------------------------------------------
! This routine defines a new boundary vertex at the "midpoint" of
! refinement_edge, meaning a point on refinement edge that is equidistant
! from the endpoints of refinement_edge.  vert_mid is the local ID of the
! vertex at the midpoint.  The results are defined in the grid data structure.
! vertex(vert_mid)%coord gets the coordinates of the midpoint.
! vertex(vert_mid)%boundary_vertex gets a linked list of boundary vertices
! associated with vert_mid.  The boundary vertices contain additional
! information about the boundary to be added to the usual vertex information.
! A vertex has a boundary vertex associated to it for each of the boundary
! surface pieces that contain the vertex.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout), target :: grid
integer, intent(in) :: refinement_edge, vert_mid
!----------------------------------------------------
! Local variables:

type(bvert_type), pointer :: bvert1, bvert2, this_bvert, new_bvert
integer :: vert1, vert2
!----------------------------------------------------
! Begin executable code

! convenience variables for the vertices

vert1 = grid%edge(refinement_edge)%vertex(1)
vert2 = grid%edge(refinement_edge)%vertex(2)
nullify(this_bvert,new_bvert)

! for each associated boundary vertex of vert1, see if there is an associated
! boundary vertex of vert2 with the same surface piece, and process if there is

bvert1 => grid%vertex(vert1)%boundary_vertex
do while (associated(bvert1))
   bvert2 => grid%vertex(vert2)%boundary_vertex
   do while (associated(bvert2))
      if (bvert1%surface == bvert2%surface) then
         call new_boundary_vertex(grid,bvert1,bvert2,new_bvert,vert1, &
                                  vert2,vert_mid)
         if (associated(this_bvert)) then
            this_bvert%next => new_bvert
            this_bvert => this_bvert%next
         else
            grid%vertex(vert_mid)%boundary_vertex => new_bvert
            this_bvert => grid%vertex(vert_mid)%boundary_vertex
         endif
         nullify(new_bvert)
         exit
      endif
      bvert2 => bvert2%next
   end do
   bvert1 => bvert1%next
end do

end subroutine define_boundary_vertex

!          -------------------
subroutine new_boundary_vertex(grid,bvert1,bvert2,this_bvert,vert1,vert2, &
                               vert_mid)
!          -------------------

!----------------------------------------------------
! This routine defines a new boundary vertex in this_bvert as the point
! equidistant from bvert1 and bvert2 on their common surface piece.
! Also, the coordinate of vertex vert_mid is set.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout), target :: grid
type(bvert_type), intent(in), target :: bvert1, bvert2
type(bvert_type), pointer :: this_bvert
integer, intent(in) :: vert1, vert2, vert_mid
!----------------------------------------------------
! Local variables:

real(my_real) :: t, a1, a2
integer :: line
!----------------------------------------------------
! Begin executable code

! convenience variables

coord1 = grid%vertex(vert1)%coord
coord2 = grid%vertex(vert2)%coord

select case (grid%surface(bvert1%surface)%type)

case (SURF_PLANAR)

! planar surface

   line = common_line(bvert1,bvert2)

   if (line /= 0) then

      if (grid%line(line)%type == LINE_CIRCLE .or. &
          grid%line(line)%type == LINE_ELLIPSE) then

! bvert1 and bvert2 are on the same curved edge, find the "midpoint" of the edge

         bv1 => bvert1
         bv2 => bvert2
         hold_grid => grid
         hold_line = line
         if (bvert1%lines(1) == line) then
            a1 = bvert1%p1
         else
            a1 = bvert1%p2
         endif
         if (bvert2%lines(1) == line) then
            a2 = bvert2%p1
         else
            a2 = bvert2%p2
         endif
         if (a2 < a1) then
            t = a1
            a1 = a2
            a2 = t
         endif
         t = zeroin(a1,a2,diff_dist,100*epsilon(0.0_my_real))
         allocate(this_bvert)
         this_bvert%surface = bvert1%surface
         nullify(this_bvert%next)
         this_bvert%lines = (/line,0/)
         this_bvert%p1 = t
         this_bvert%p2 = 0
         grid%vertex(vert_mid)%coord = line_point(t,grid,line)

      else

! bvert1 and bvert2 are on the same straight edge, average to get the midpoint

         allocate(this_bvert)
         this_bvert%surface = bvert1%surface
         nullify(this_bvert%next)
         this_bvert%lines = (/line,0/)
         if (bvert1%lines(1) == line) then
            a1 = bvert1%p1
         else
            a1 = bvert1%p2
         endif
         if (bvert2%lines(1) == line) then
            a2 = bvert2%p1
         else
            a2 = bvert2%p2
         endif
         this_bvert%p1 = (a1+a2)/2
         this_bvert%p2 = 0
         grid%vertex(vert_mid)%coord = average_coords(coord1,coord2)

      endif

! bvert1 and bvert2 are not on a common line, average to get midpoint

   else

      allocate(this_bvert)
      this_bvert%surface = bvert1%surface
      nullify(this_bvert%next)
      this_bvert%lines = (/0,0/)
      this_bvert%p1 = 0
      this_bvert%p2 = 0
      grid%vertex(vert_mid)%coord = average_coords(coord1,coord2)

   endif

case (SURF_RULED_3, SURF_RULED_4)

! if it is a ruled surface, use a root finder on the transfinite interpolation
! that defines the surface piece to find the midpoint

   line = common_line(bvert1,bvert2)
   bv1 => bvert1
   bv2 => bvert2
   hold_grid => grid
   t = zeroin(0.0_my_real,1.0_my_real,diff_dist,100*epsilon(0.0_my_real))
   allocate(this_bvert)
   this_bvert%surface = bvert1%surface
   nullify(this_bvert%next)
   this_bvert%p1 = (1-t)*bvert1%p1 + t*bvert2%p1
   this_bvert%p2 = (1-t)*bvert1%p2 + t*bvert2%p2
   grid%vertex(vert_mid)%coord = surface_point(this_bvert%p1,this_bvert%p2, &
                                               bvert1%surface,grid)
   this_bvert%lines = (/line,0/)

case default

   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized surface type in new_boundary_vertex")
   stop

end select

contains

function average_coords(c1,c2)
type(point), intent(in) :: c1, c2
type(point) :: average_coords
average_coords = point((c1%x+c2%x)/2,(c1%y+c2%y)/2,(c1%z+c2%z)/2)
end function average_coords

function common_line(b1,b2)
type(bvert_type), intent(in) :: b1,b2
integer :: common_line
if (b1%lines(1) == b2%lines(1)) then
   common_line = b1%lines(1)
elseif (b1%lines(1) == b2%lines(2)) then
   common_line = b1%lines(1)
elseif (b1%lines(2) == b2%lines(1)) then
   common_line = b1%lines(2)
elseif (b1%lines(2) == b2%lines(2)) then
   common_line = b1%lines(2)
else
   common_line = 0
endif
end function common_line

end subroutine new_boundary_vertex

!        ------
function zeroin(ax,bx,f,tol)
!        ------

!----------------------------------------------------
! This is a root finding routine from netlib.
! My modifications to free source form and other modernizations

!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).

!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real) :: ax,bx,tol
real(my_real) :: zeroin
interface
   function f(x)
   use global
   real(my_real) :: x
   real(my_real) :: f
   end function f
end interface

!----------------------------------------------------
! Local variables:

real(my_real) ::  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s

!----------------------------------------------------
! Begin executable code

!   10 eps = d1mach(4)
   10 eps = epsilon(0.0_my_real)
      tol1 = eps+1.0_my_real

      a=ax
      b=bx
      fa=f(a)
      fb=f(b)
!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0_my_real .or. fb .eq. 0.0_my_real) go to 20
      if (fa * (fb/abs(fb)) .le. 0.0d0) go to 20
      ierr = PHAML_INTERNAL_ERROR
      call fatal("zeroin: f(ax) and f(bx) do not have different signs")
      stop
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (abs(fc).ge.abs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*abs(b)+0.5_my_real*tol
      xm = 0.5_my_real*(c-b)
      if ((abs(xm).le.tol1).or.(fb.eq.0.0_my_real)) go to 150

! see if a bisection is forced

      if ((abs(e).ge.tol1).and.(abs(fa).gt.abs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60

! linear interpolation

      p=2.0_my_real*xm*s
      q=1.0_my_real-s
      go to 70

! inverse quadratic interpolation

   60 q=fa/fc
      r=fb/fc
      p=s*(2.0_my_real*xm*q*(q-r)-(b-a)*(r-1.0_my_real))
      q=(q-1.0_my_real)*(r-1.0_my_real)*(s-1.0_my_real)
   70 if (p.le.0.0_my_real) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0_my_real*p).ge.(3.0_my_real*xm*q-abs(tol1*q))).or. &
          (p.ge.abs(0.5_my_real*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (abs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0_my_real) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
  140 fb=f(b)
      if ((fb*(fc/abs(fc))).gt.0.0_my_real) go to 20
      go to 30
  150 zeroin=b

end function zeroin

!        ------------
function diff_dist(t)
!        ------------

!----------------------------------------------------
! This routine gives the difference between the distances of the point
! t'th of the way from boundary vertex bv1 to bv2
! hold_grid, bv1, bv2, coord1 and coord2 need to be module variables
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real) :: t
real(my_real) :: diff_dist
!----------------------------------------------------
! Local variables:

real(my_real) :: p1, p2
type(point) :: coord_candidate
!----------------------------------------------------
! Begin executable code

if (hold_grid%surface(bv1%surface)%type == SURF_PLANAR) then
   coord_candidate = line_point(t,hold_grid,hold_line)
else ! SURF_RULED_3 or SURF_RULED_4
   p1 = (1-t)*bv1%p1 + t*bv2%p1
   p2 = (1-t)*bv1%p2 + t*bv2%p2
   coord_candidate = surface_point(p1,p2,bv1%surface,hold_grid)
endif
diff_dist = dist(coord_candidate,coord1) - dist(coord_candidate,coord2)

end function diff_dist

!        ----
function dist(p1,p2)
!        ----

!----------------------------------------------------
! This routine computes the distance between two points
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: p1, p2
real(my_real) :: dist
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

dist = sqrt((p1%x-p2%x)**2 + (p1%y-p2%y)**2 + (p1%z-p2%z)**2)

end function dist

!        -------------
function surface_point(p1,p2,surface,grid)
!        -------------

!----------------------------------------------------
! This routine computes the point whose transfinite interpolation with
! parameters p1, p2 on the given surface piece
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: p1, p2
integer, intent(in) :: surface
type(grid_type), intent(in), target :: grid
type(point) :: surface_point
!----------------------------------------------------
! Local variables:

type(lineloop_type), pointer :: lineloop
type(point) :: v_1_0_0, v_0_1_0, v_0_0_1, v_1mb2_b2_0, v_1mb3_0_b3, &
               v_0_1mb3_b3, v_b1_1mb1_0, v_b1_0_1mb1, v_0_b2_1mb2, &
               v_0_eta, v_1_eta, v_xi_0, v_xi_1, v_0_0, v_1_0, v_1_1, v_0_1, &
               center, p
real(my_real) :: b1, b2, b3, xi, eta
logical :: same_center
!----------------------------------------------------
! Begin executable code

! the line loop for this surface; convenience variable

lineloop => grid%lineloop(grid%surface(surface)%lineloop(1))

select case (grid%surface(surface)%type)

case (SURF_RULED_3)

! Ruled surface with 3 lines forms a triangular patch.
! In barycentric coordinates, the first line goes from v(1,0,0) to v(0,1,0).
! A point on the line, v(A,1-A,0), is given by
! line_point(1-A,grid,lineloop%line(1)), because, for example, if A=1 we want
! to go 0th of the way from the first endpoint to the second end point on line 1
! The second line goes from v(0,1,0) to v(0,0,1).
! A point on the line, v(0,A,1-A), is given by
! line_point(1-A,grid,lineloop%line(2)).
! The third line goes from v(0,0,1) to v(1,0,0).
! A point on the line, v(1-A,0,A), is given by
! line_point(1-A,grid,lineloop%line(3)).

! set the barycentric coordinates in the surface patch

   b1 = p1
   b2 = p2
   b3 = 1-p1-p2

! evaluate the line coordinates at the corners and edges

   v_1_0_0 = line_point(0.0_my_real,grid,lineloop%line(1))
   v_0_1_0 = line_point(0.0_my_real,grid,lineloop%line(2))
   v_0_0_1 = line_point(0.0_my_real,grid,lineloop%line(3))
   v_1mb2_b2_0 = line_point(b2,grid,lineloop%line(1))
   v_1mb3_0_b3 = line_point(1-b3,grid,lineloop%line(3))
   v_0_1mb3_b3 = line_point(b3,grid,lineloop%line(2))
   v_b1_1mb1_0 = line_point(1-b1,grid,lineloop%line(1))
   v_b1_0_1mb1 = line_point(b1,grid,lineloop%line(3))
   v_0_b2_1mb2 = line_point(1-b2,grid,lineloop%line(2))

! interpolate

   surface_point = b1*(v_1mb2_b2_0 + v_1mb3_0_b3 - v_1_0_0) + &
                   b2*(v_0_1mb3_b3 + v_b1_1mb1_0 - v_0_1_0) + &
                   b3*(v_b1_0_1mb1 + v_0_b2_1mb2 - v_0_0_1)

case (SURF_RULED_4)

! Ruled surface with 4 lines forms a quadrilateral patch.
! In (xi,eta) coordinates, the first line goes from v(0,0) to v(1,0).
! A point on the line, v(A,0) is given by
! line_point(A,grid,lineloop%line(1)).
! Similarly, the second, third and fourth lines go from v(1,0) to v(1,1) to
! v(0,1) to v(0,0) with points on the lines given by
! v(1,A); line_point(A,grid,lineloop%line(2))
! v(A,1); line_point(1-A,grid,lineloop%line(3))
! v(0,A); line_point(1-A,grid,lineloop%line(4))

! set (xi,eta), for clarity

   xi = p1
   eta = p2

! evaluate the line coordinates at the corners and edges

   v_0_0 = line_point(0.0_my_real,grid,lineloop%line(1))
   v_1_0 = line_point(0.0_my_real,grid,lineloop%line(2))
   v_1_1 = line_point(0.0_my_real,grid,lineloop%line(3))
   v_0_1 = line_point(0.0_my_real,grid,lineloop%line(4))
   v_0_eta = line_point(1-eta,grid,lineloop%line(4))
   v_1_eta = line_point(eta,grid,lineloop%line(2))
   v_xi_0 = line_point(xi,grid,lineloop%line(1))
   v_xi_1 = line_point(1-xi,grid,lineloop%line(3))

! interpolate

   surface_point = (1-xi)*v_0_eta + xi*v_1_eta + (1-eta)*v_xi_0 + eta*v_xi_1 &
                  -((1-xi)*(1-eta)*v_0_0 + xi*(1-eta)*v_1_0 + xi*eta*v_1_1 + &
                    (1-xi)*eta*v_0_1)

case (SURF_PLANAR)

   ierr = PHAML_INTERNAL_ERROR
   call fatal("cannot call surface_point with a planar surface piece")
   stop

case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized surface type in surface_point")
   stop

end select

! if we are applying the Gmsh sphere hack, then if all edges of this surface
! are circle arcs with the same center, move the point to the sphere with
! that center and radius of the circle arcs

if (gmsh_sphere_hack) then
   if (all(grid%line(abs(lineloop%line(:)))%type == LINE_CIRCLE)) then
      same_center = .true.
      if (grid%line(abs(lineloop%line(2)))%point(2) /= &
          grid%line(abs(lineloop%line(1)))%point(2)) same_center = .false.
      if (grid%line(abs(lineloop%line(3)))%point(2) /= &
          grid%line(abs(lineloop%line(1)))%point(2)) same_center = .false.
      if (grid%surface(surface)%type == SURF_RULED_4) then
         if (grid%line(abs(lineloop%line(4)))%point(2) /= &
             grid%line(abs(lineloop%line(1)))%point(2)) same_center = .false.
      endif
      if (same_center) then
         center = grid%point(grid%line(abs(lineloop%line(1)))%point(2))%coord
         p      = grid%point(grid%line(abs(lineloop%line(1)))%point(1))%coord
         surface_point = (dist(p,center)/dist(surface_point,center)) * &
                         (surface_point - center) + center
      endif
   endif
endif

end function surface_point

!        ----------
function line_point(A,grid,lineID)
!        ----------

!----------------------------------------------------
! This routine returns the point A'th of the way along the given line.
! If lineID is negative, it starts at the end of the line and goes backwards.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: A
type(grid_type), intent(in), target :: grid
integer, intent(in) :: lineID
type(point) :: line_point
!----------------------------------------------------
! Local variables:

real(my_real), parameter :: roundoff = 100*epsilon(0.0_my_real)
type(line_type), pointer :: line
type(point) :: P1, P2, C, M, u, v, tempp
real(my_real) :: r, theta, alpha, beta, temp, f1, f2, f3, f4, cost1, cost2, &
                 t1, t2, t
!----------------------------------------------------
! Begin executable code

! line, convenience variable

if (lineID < 0) then
   line => grid%line(-lineID)
else
   line => grid%line(lineID)
endif

select case(line%type)

case(LINE_LINE)

! A straight line.  Simple linear combination of the two points.

   if (lineID < 0) then
      line_point =    A *grid%point(line%point(1))%coord + &
                   (1-A)*grid%point(line%point(2))%coord
   else
      line_point =    A *grid%point(line%point(2))%coord + &
                   (1-A)*grid%point(line%point(1))%coord
   endif

case(LINE_CIRCLE)

! The point P needs to lie in the plane (P1,P2,C) where P1 is the starting point
! of the circular arc, P2 is the ending point and C is the center, and needs
! to have dist(P,C) = r := dist(P1,C) = dist(P2,C).  Being in that plane
! means P-C is in the span of {P1-C, P2-C}, i.e.
! P = C + alpha*(P1-C) + beta*(P2-C)
! The point P is found by starting at C, moving in the direction of P1-C
! a distance of alpha, then moving in the direction of P2-C a distance of beta.
! The triangle from C to alpha(P1-C) to P back to C has sides of length alpha*r,
! beta*r and r.  Let theta be the angle between lines P1-C and P2-C, i.e.
! theta = acos(<P1-C,P2-C>/r^2), 0 < theta < pi.
! We want to go A'th of the way from P1 to P2, so the angle between the line
! from C to P1 and the line from C to P (the first angle of our triangle) is
! A*theta.  The line from C+alpha*(P1-C) to P is parallel to the line from
! C to P2, so the second angle is pi-theta.  This makes the third angle
! (1-A)*theta.  Using the law of sines and the identity for sin(x-y) we find
! alpha = sin((1-A)*theta)/sin(theta) and beta = sin(A*theta)/sin(theta).

   if (lineID < 0) then
      P1 = grid%point(line%point(3))%coord
      P2 = grid%point(line%point(1))%coord
   else
      P1 = grid%point(line%point(1))%coord
      P2 = grid%point(line%point(3))%coord
   endif
   C = grid%point(line%point(2))%coord

   r = sqrt(dot_point(P1-C,P1-C))
   if (r == 0.0_my_real) then
      ierr = USER_INPUT_ERROR
      call fatal("line_point: looks like a point on a circle is the center of the circle")
      stop
   endif
   theta = acos(dot_point(P1-C,P2-C)/r**2)
   if (theta == 0.0_my_real) then
      ierr = USER_INPUT_ERROR
      call fatal("line_point: looks like the two endpoints of a circular arc are the same point")
      stop
   endif

   alpha = sin((1-A)*theta)/sin(theta)
   beta  = sin(A*theta)/sin(theta)

   line_point = C + alpha*(P1-C) + beta*(P2-C)

case(LINE_ELLIPSE)

! Given the starting point of the elliptical arc, P1, the end, P2, the center
! of the ellipse, C, and another point on the major axis, M.
! Let u be the unit vector in the direction of the major axis, i.e.
! u = (M-C)/||M-C||
! Let v be the unit vector in the direction of the minor axis on the same
! side as P1 (if P1 lies on the major axis, use P2).  v is in span{u,P1-C}
! and <v,u> = 0 (v is orthogonal to u), so
! v = (P1-C - <P1-C,u> u)/||its norm||
! An ellipse in this coordinate system in parametric form is given by
! P = C + a*cos(t)*u + b*sin(t)*v
! where a and b are the lengths of the major and minor axes.  We'll use
! alpha and beta for a and b.
! Taking inner products of P1-C and P2-C with u and v, we find
! a*cos(t1) = <P1-C,u>
! b*sin(t1) = <P1-C,v>
! a*cos(t2) = <P2-C,u>
! b*sin(t2) = <P2-C,v>
! where t1 and t2 are the parameters for P1 and P2.
! This system of 4 nonlinear equations in 4 unknowns a, b, cos(t1) and cos(t2)
! (replace sin(t1) with sqrt(1-cos^2(t1)) can be easily solved analytically.
! To get t1 (same for t2) from its cosine, use acos and test the result in
! the formula to see if it gives P1 back.  If not, t1 must be the other choice
! of arccosine, 2*pi-acos.

   if (lineID < 0) then
      P1 = grid%point(line%point(4))%coord
      P2 = grid%point(line%point(1))%coord
   else
      P1 = grid%point(line%point(1))%coord
      P2 = grid%point(line%point(4))%coord
   endif
   C = grid%point(line%point(2))%coord
   M = grid%point(line%point(3))%coord

   u = M-C
   temp = sqrt(dot_point(u,u))
   if (temp == 0.0_my_real) then
      ierr = USER_INPUT_ERROR
      call fatal("line_point: point on major axis of ellipse is the same as the center")
      stop
   endif
   u = u/temp

   f1 = dot_point(P1-C,u)
   f3 = dot_point(P2-C,u)
   v = P1-C - f1*u
   temp = sqrt(dot_point(v,v))
   if (temp == 0.0_my_real) then
      v = P2-C - f3*u
      temp = sqrt(dot_point(v,v))
      if (temp == 0) then
         ierr = USER_INPUT_ERROR
         call fatal("line_point: starting and ending points of elliptical arc are both on major axis")
         stop
      endif
   endif
   v = v/temp

   f2 = dot_point(P1-C,v)
   f4 = dot_point(P2-C,v)

   f1 = f1*f1
   f2 = f2*f2
   f3 = f3*f3
   f4 = f4*f4
   temp = f1*f4 - f2*f3

   alpha = sqrt(temp/(f4-f2))
   beta  = sqrt(temp/(f1-f3))
   cost1 = sqrt(f1)/a
   cost2 = sqrt(f3)/a

   t1 = acos(cost1)
   tempp = C + alpha*cos(t1)*u + beta*sin(t1)*v - P1
   temp = dot_point(tempp,tempp)
   if (temp > roundoff) then
      t1 = 8*atan(1.0_my_real) - t1
      tempp = C + alpha*cos(t1)*u + beta*sin(t1)*v - P1
      temp = dot_point(tempp,tempp)
      if (temp > roundoff) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("line_point: neither choice of acos worked")
      endif
   endif

   t2 = acos(cost2)
   tempp = C + alpha*cos(t2)*u + beta*sin(t2)*v - P2
   temp = dot_point(tempp,tempp)
   if (temp > roundoff) then
      t2 = 8*atan(1.0_my_real) - t2
      tempp = C + alpha*cos(t2)*u + beta*sin(t2)*v - P2
      temp = dot_point(tempp,tempp)
      if (temp > roundoff) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("line_point: neither choice of acos worked")
      endif
   endif

! to go A'th of the way from P1 to P2, use a root finder to find the zero of
! ||ellipse(t)-P1|| - A*||P2-P1||.  Use module variables to get the ellipse
! etc to the function to be zeroed.

   e_C = C
   e_alpha = alpha
   e_beta = beta
   e_u = u
   e_v = v
   e_p1 = P1
   e_p2 = P2
   e_A = A

   t = zeroin(t1,t2,ellipse_frac,100*epsilon(0.0_my_real))

   line_point = C + alpha*cos(t)*u + beta*sin(t)*v

case default

   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized line type in line_point")
   stop

end select

end function line_point

!        ------------
function ellipse_frac(t)
!        ------------

!----------------------------------------------------
! This routine is used when computing a point e_A'th of the way from
! e_P1 to e_P2 on an ellipse.  It returns
! ||ellipse(t) - P1|| - e_A*||P2-P1||
! Lots of stuff to define this is in module variables.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real) :: t
real(my_real) :: ellipse_frac
!----------------------------------------------------
! Local variables:

type(point) :: ellipse_point
!----------------------------------------------------
! Begin executable code

ellipse_point = e_C + e_alpha*cos(t)*e_u + e_beta*sin(t)*e_v
ellipse_frac = sqrt(dot_point(ellipse_point-e_P1,ellipse_point-e_P1)) - &
               e_A*sqrt(dot_point(e_P2-e_P1,e_P2-e_P1))

end function ellipse_frac

!          -------------------
subroutine deallocate_boundary(grid)
!          -------------------

!----------------------------------------------------
! This routine deallocates memory associated with a general boundary
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: i, lev, vert, elem, ivert
logical :: visited_vert(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

if (associated(grid%point)) then
   deallocate(grid%point)
endif

if (associated(grid%line)) then
   do i=1,grid%nline
      deallocate(grid%line(i)%point)
   end do
   deallocate(grid%line)
endif

if (associated(grid%lineloop)) then
   do i=1,grid%nlineloop
      deallocate(grid%lineloop(i)%line)
   end do
   deallocate(grid%lineloop)
endif

if (associated(grid%surface)) then
   do i=1,grid%nsurface
      deallocate(grid%surface(i)%lineloop)
   end do
   deallocate(grid%surface)
endif

visited_vert = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do ivert=1,VERTICES_PER_ELEMENT
         vert = grid%element(elem)%vertex(ivert)
         if (visited_vert(vert)) cycle
         visited_vert(vert) = .true.
         if (associated(grid%vertex(vert)%boundary_vertex)) then
            call deallocate_boundary_vertex_next(grid%vertex(vert)%boundary_vertex)
            deallocate(grid%vertex(vert)%boundary_vertex)
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

end subroutine deallocate_boundary

!                    -------------------------------
recursive subroutine deallocate_boundary_vertex_next(bvert)
!                    -------------------------------

!----------------------------------------------------
! This routine deallocates linked list of boundary vertices under bvert%next
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(bvert_type) :: bvert
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (associated(bvert%next)) then
   call deallocate_boundary_vertex_next(bvert%next)
   deallocate(bvert%next)
endif

end subroutine deallocate_boundary_vertex_next

!          ------------------
subroutine point_surface_dist(n,x,nf,fx,uiparm,urparm,ufparm)
!          ------------------

!----------------------------------------------------
! This routine returns the distance between coord1 and the point on surface
! uiparm(1).  The parameter list is dictated by subroutine dmnfb.
! module parameters coord1 and hold_grid must be set.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: n                   ! size of x
double precision :: x(*)       ! point at which to evaluate
integer :: nf                  ! invocation count
double precision :: fx         ! result of function evaluation
integer :: uiparm(*)           ! vector of integer parameters
double precision :: urparm(*)  ! vector of real parameters
external ufparm                ! function parameter
!----------------------------------------------------
! Local variables:

type(point) :: surf_point
real(my_real) :: p1, p2
!----------------------------------------------------
! Begin executable code

p1 = x(1)
p2 = x(2)
surf_point = surface_point(p1,p2,uiparm(1),hold_grid)
fx = dist(surf_point,coord1)

end subroutine point_surface_dist

!        -------------------
function coord_is_on_surface(coord,surf,grid,p1,p2,resid)
!        -------------------

!----------------------------------------------------
! This routine determines whether or not the point with the given coordinates
! is on surface surf.  If it is, (p1,p2) are the parameters that give the
! point, i.e. the first two barycentric coordinates for RULED_3 or (xi,eta)
! for RULED_4.  resid is the distance between the given point and the closest
! point on the surface.  Since the determination is made by an optimization
! routine that can get trapped in local minima, it is possible for this
! routine to fail if the surface patch is not convex.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: coord
integer, intent(in) :: surf
type(grid_type), intent(in), target :: grid
real(my_real), intent(out) :: p1, p2, resid
logical :: coord_is_on_surface
!----------------------------------------------------
! Local variables:

integer, parameter :: liv = 61, lv = 106
integer :: iv(liv), uiparm(1)
double precision :: x(2), v(lv), b(2,2), urparm(1), norm_point, norm_coord, &
                    bigger_norm
type(point) :: surf_point
external dummyf
!----------------------------------------------------
! Begin executable code

uiparm(1) = surf
urparm(1) = 0.0d0
hold_grid => grid
coord1 = coord

! initital guess

select case (grid%surface(surf)%type)
case (SURF_RULED_3)
   x = (/0.3333333d0,0.3333333d0/)
case (SURF_RULED_4)
   x = (/0.5d0,0.5d0/)
case (SURF_PLANAR)
   ierr = PHAML_INTERNAL_ERROR
   call fatal("cannot call coord_is_on_surface with a planar surface piece")
   stop
case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unrecognized surface type in coord_is_on_surface")
   stop
end select

b(1,1) = 0.0d0
b(2,1) = 1.0d0
b(1,2) = 0.0d0
b(2,2) = 1.0d0

call divset(2,iv,liv,lv,v)   ! set default values
iv(21) = 0 ! suppress all printing

call dmnfb(2, &              ! number of components of the solution
           (/1.0d0,1.0d0/), &! scale vector
           x, &              ! solution
           b, &              ! constraints
           point_surface_dist, & ! function to be minimized
           iv, &             ! integer control vector
           liv, &            ! size of iv
           lv, &             ! size of v
           v, &              ! real control vector
           uiparm, &         ! integers passed to the function to be minimized
           urparm, &         ! reals passed to the function to be minimized
           dummyf)           ! a function passed to the function to be minimized


p1 = x(1)
p2 = x(2)
surf_point = surface_point(p1,p2,surf,grid)
norm_point = dot_point(surf_point,surf_point)
norm_coord = dot_point(coord,coord)
bigger_norm = sqrt(max(norm_point,norm_coord))
resid = dist(coord,surf_point)
if (bigger_norm /= 0.0_my_real) resid = resid/bigger_norm
coord_is_on_surface = (resid < 1.0d-6)

end function coord_is_on_surface

end module boundary_util

!          ------
subroutine dummyf()
!          ------

!----------------------------------------------------
! This is a dummy routine to pass through dmnfb
!----------------------------------------------------

end subroutine dummyf
