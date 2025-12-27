//---------------------------------------------------------------
// Example: unit disk (circle) in 2D for Gmsh
// - Domain: circle of radius R centered at the origin
//---------------------------------------------------------------
// Parameters
R  = 1.0;      // radius of the disk
ms = R / 8.0;  // characteristic mesh size 

//---------------------------------------------------------------
// Points
// Syntax: Point(ID) = {x, y, z, meshSize};
//   - meshSize optional; if ignored, we use 
//       Mesh.CharacteristicLength to control the mesh size
//---------------------------------------------------------------
Point(1) = {0, 0, 0};
Point(2) = { R, 0, 0 };
Point(3) = { 0, R, 0 };
Point(4) = { -R, 0, 0 };
Point(5) = { 0, -R, 0 };

//---------------------------------------------------------------
// Circle arcs
// Syntax: Circle(ID) = {startPoint, centerPoint, endPoint};
//---------------------------------------------------------------
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

//---------------------------------------------------------------
// Curve loop and surface
// Syntax:
//   Curve Loop(ID)    = {curveIDs in order};
//   Plane Surface(ID) = {loopID};
//---------------------------------------------------------------
Line Loop(10) = {1,2,3,4};
Plane Surface(20) = {10};

// Optional: physical groups
Physical Surface("Domain") = {20};
Physical Curve("Dirichlet") = {1,2};
Physical Curve("Neumann") = {3,4};

//---------------------------------------------------------------
// Mesh options
//---------------------------------------------------------------
Mesh.ElementOrder            = 1;    // 1 = linear elements; 2 = quadratic
Mesh.Algorithm               = 8;    // 8 = Delaunay
Mesh.CharacteristicLengthMax = ms;   // use ms as global size hint
Mesh.CharacteristicLengthMin = ms;
Mesh.Optimize                = 1;    // improve element quality

// EOF
