// Parameters
R = 1;

// Define center and 4 points on the circle
Point(1) = {0, 0, 0};       // center
Point(2) = { R, 0, 0 };
Point(3) = { 0, R, 0 };
Point(4) = { -R, 0, 0 };
Point(5) = { 0, -R, 0 };

// Circle arcs
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Line loop and surface
Line Loop(10) = {1,2,3,4};
Plane Surface(20) = {10};

// Optional: physical groups
Physical Surface("Disk") = {20};
Physical Curve("Boundary") = {1,2,3,4};

