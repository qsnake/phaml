/*
This file defines a cylinder with center up the z axis, bottom in the
x-y plane, and height and radius given in the next two lines.
The only boundary marker set is the bottom plane with bmark=1
*/

radius = 1.0;
height = 1.0;
size   = 1.0;  /* Min(radius,height), but PHAML doesn't parse functions in
                  expressions, and Gmsh doesn't have Min in its functions
                  anyway */

/* line up the center */

Point(1) = {0,0,0,size};
Point(2) = {0,0,height,size};
Line(1) = {1,2};

/* vertical lines on the +x, +y, -x and -y axes */

Point(3) = {radius,0,0,size};
Point(4) = {radius,0,height,size};
Line(2) = {3,4};

Point(5) = {0,radius,0,size};
Point(6) = {0,radius,height,size};
Line(3) = {5,6};

Point(7) = {-radius,0,0,size};
Point(8) = {-radius,0,height,size};
Line(4) = {7,8};

Point( 9) = {0,-radius,0,size};
Point(10) = {0,-radius,height,size};
Line(5) = {9,10};

/* lines along the axes on the bottom surface */

Line(6) = {1,3};
Line(7) = {1,5};
Line(8) = {1,7};
Line(9) = {1,9};

/* lines along the axes on the top surface */

Line(10) = {2,4};
Line(11) = {2,6};
Line(12) = {2,8};
Line(13) = {2,10};

/* circular lines on the bottom */

Circle(14) = {3,1,5};
Circle(15) = {5,1,7};
Circle(16) = {7,1,9};
Circle(17) = {9,1,3};

/* circular lines on the top */

Circle(18) = {4,2,6};
Circle(19) = {6,2,8};
Circle(20) = {8,2,10};
Circle(21) = {10,2,4};

/* faces on the bottom */

Line Loop(1) = {6,14,-7};
Plane Surface(1) = {1};

Line Loop(2) = {7,15,-8};
Plane Surface(2) = {2};

Line Loop(3) = {8,16,-9};
Plane Surface(3) = {3};

Line Loop(4) = {9,17,-6};
Plane Surface(4) = {4};

/* faces on the top */

Line Loop(5) = {10,18,-11};
Plane Surface(5) = {5};

Line Loop(6) = {11,19,-12};
Plane Surface(6) = {6};

Line Loop(7) = {12,20,-13};
Plane Surface(7) = {7};

Line Loop(8) = {13,21,-10};
Plane Surface(8) = {8};

/* circular faces that form the cylinder side */

Line Loop(9) = {14,3,-18,-2};
Ruled Surface(9) = {9};

Line Loop(10) = {15,4,-19,-3};
Ruled Surface(10) = {10};

Line Loop(11) = {16,5,-20,-4};
Ruled Surface(11) = {11};

Line Loop(12) = {17,2,-21,-5};
Ruled Surface(12) = {12};

/* final volume */

Surface Loop(1) = {1,2,3,4,12,11,10,6,7,8,5,9};
Volume(1) = {1};

/* boundary marker for the bottom surface */

/* this is causing the .msh file to only contain the bottom surface.
   need to track down that problem
Physical Surface(100) = {1,2,3,4};
*/
