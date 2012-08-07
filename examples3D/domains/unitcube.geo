/*
This file defines a cube with all boundary entities given unique bmarks.
The bmarks start at the bottom and origin and work their way around
counterclockwise and upwards.  The values of bmark are given by the
Physical entities.
See the figure in the section on 3D problems in the PHAML user's guide.
*/

/* Vertices of the unit cube */

Point(1) = {0,0,0,1};
Physical Point(2) = {1};
Point(2) = {1,0,0,1};
Physical Point(4) = {2};
Point(3) = {1,1,0,1};
Physical Point(6) = {3};
Point(4) = {0,1,0,1};
Physical Point(8) = {4};
Point(5) = {0,0,1,1};
Physical Point(18) = {5};
Point(6) = {1,0,1,1};
Physical Point(20) = {6};
Point(7) = {1,1,1,1};
Physical Point(22) = {7};
Point(8) = {0,1,1,1};
Physical Point(24) = {8};

/* Lines around the bottom */

Line(1) = {1,2};
Physical Line(3) = {1};
Line(2) = {2,3};
Physical Line(5) = {2};
Line(3) = {3,4};
Physical Line(7) = {3};
Line(4) = {4,1};
Physical Line(9) = {4};

/* Vertical lines */

Line(5) = {1,5};
Physical Line(10) = {5};
Line(6) = {2,6};
Physical Line(12) = {6};
Line(7) = {3,7};
Physical Line(14) = {7};
Line(8) = {4,8};
Physical Line(16) = {8};

/* Lines around the top */

Line( 9) = {5,6};
Physical Line(19) = {9};
Line(10) = {6,7};
Physical Line(21) = {10};
Line(11) = {7,8};
Physical Line(23) = {11};
Line(12) = {8,5};
Physical Line(25) = {12};

/* Bottom */

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

/* Front, i.e. x-z plane with y=0 */

Line Loop(2) = {1,6,-9,-5};
Plane Surface(2) = {2};
Physical Surface(11) = {2};

/* Right */

Line Loop(3) = {2,7,-10,-6};
Plane Surface(3) = {3};
Physical Surface(13) = {3};

/* Back */

Line Loop(4) = {3,8,-11,-7};
Plane Surface(4) = {4};
Physical Surface(15) = {4};

/* Left */

Line Loop(5) = {4,5,-12,-8};
Plane Surface(5) = {5};
Physical Surface(17) = {5};

/* Top */

Line Loop(6) = {9,10,11,12};
Plane Surface(6) = {6};
Physical Surface(26) = {6};

/* Volume */

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};
Physical Volume(0) = {1};
