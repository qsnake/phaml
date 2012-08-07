/*
This file defines a unit tetrahedron.
*/

/* Vertices of the tetrahedron */

Point(1) = {0,0,0,1};
Point(2) = {1,0,0,1};
Point(3) = {0,1,0,1};
Point(4) = {0,0,1,1};

/* Edges of the tetrahedron */

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};
Line(4) = {1,4};
Line(5) = {2,4};
Line(6) = {3,4};

/* Faces of the tetrahedron */

Line Loop(1) = {1,2,3};
Plane Surface(1) = {1};
Line Loop(2) = {1,5,-4};
Plane Surface(2) = {2};
Line Loop(3) = {3,4,-6};
Plane Surface(3) = {3};
Line Loop(4) = {2,6,-5};
Plane Surface(4) = {4};

/* Volume */

Surface Loop(1) = {1,2,3,4};
Volume(1) = {1};
