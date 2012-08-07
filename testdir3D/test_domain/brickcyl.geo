/*
This file defines a brick with a cylinder on top
*/

/* extents of the brick */

xmin = -2;
xmax =  2;
ymin = -1;
ymax =  1;
zmin =  0;
zmax =  1.5;

/* radius and height of the cylinder */

radius = .5;
height = .5;

/* size of elements near cylinders and elsewhere */

size1 = .5;
size2 = 2;

/* the cylinder */

/* the points at the center of the circles at the bottom and top */

Point(1) = {0,0,zmax,size1};
Point(2) = {0,0,zmax+height,size1};

/* vertical lines on the +x, +y, -x and -y axes */

Point(3) = {radius,0,zmax,size1};
Point(4) = {radius,0,zmax+height,size1};
Line(2) = {3,4};

Point(5) = {0,radius,zmax,size1};
Point(6) = {0,radius,zmax+height,size1};
Line(3) = {5,6};

Point(7) = {-radius,0,zmax,size1};
Point(8) = {-radius,0,zmax+height,size1};
Line(4) = {7,8};

Point( 9) = {0,-radius,zmax,size1};
Point(10) = {0,-radius,zmax+height,size1};
Line(5) = {9,10};

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

/* line loop around the bottom */

Line Loop(1) = {14,15,16,17};

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

/* Vertices of the parallelepiped */

Point(11) = {xmin,ymin,zmin,size2};
Point(12) = {xmax,ymin,zmin,size2};
Point(13) = {xmax,ymax,zmin,size2};
Point(14) = {xmin,ymax,zmin,size2};
Point(15) = {xmin,ymin,zmax,size2};
Point(16) = {xmax,ymin,zmax,size2};
Point(17) = {xmax,ymax,zmax,size2};
Point(18) = {xmin,ymax,zmax,size2};

/* Lines around the bottom */

Line(31) = {11,12};
Line(32) = {12,13};
Line(33) = {13,14};
Line(34) = {14,11};

/* Vertical lines */

Line(35) = {11,15};
Line(36) = {12,16};
Line(37) = {13,17};
Line(38) = {14,18};

/* Lines around the top */

Line(39) = {15,16};
Line(40) = {16,17};
Line(41) = {17,18};
Line(42) = {18,15};

/* Bottom */

Line Loop(21) = {31,32,33,34};
Plane Surface(21) = {21};

/* Front, i.e. x-z plane with y=ymin */

Line Loop(22) = {31,36,-39,-35};
Plane Surface(22) = {22};

/* Right */

Line Loop(23) = {32,37,-40,-36};
Plane Surface(23) = {23};

/* Back */

Line Loop(24) = {33,38,-41,-37};
Plane Surface(24) = {24};

/* Left */

Line Loop(25) = {34,35,-42,-38};
Plane Surface(25) = {25};

/* Top, including a hole where the cylinder sits */

Line Loop(26) = {39,40,41,42};
Plane Surface(26) = {26,1};

/* Volume */

Surface Loop(1) = {21,22,23,24,25,26,12,11,10,6,7,8,5,9};
Volume(1) = {1};
