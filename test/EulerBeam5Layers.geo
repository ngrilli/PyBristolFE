// Nicolò Grilli
// University of Bristol
// 15 Giugno 2023

// Euler beam with five layers of elements to avoid volumetric locking

lc = 0.2;

Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0, 0, lc};
Point(3) = {10, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Surface("left") = {4};
Physical Surface("right") = {2};
