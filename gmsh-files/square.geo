//---------------------------------------------------------------
// Example: simple square geometry for Gmsh
// Demonstrates the grammar and structure of a .geo script
//---------------------------------------------------------------
// Define a parameter (variable) for the square's side length.
// Variables in Gmsh are declared with a name and value using '='.
// You can reuse them later by their name.
L  = 1.0; // domain length
ms = 1.0/100.0; // characteristic mesh size

//---------------------------------------------------------------
// Define corner points of the square.
// The general syntax is: Point(ID) = {x, y, z, meshSize};
// - ID is the point label (an integer, unique in the file)
// - x, y, z are coordinates in 3D
// - meshSize (optional) sets the target element size near this point
//---------------------------------------------------------------
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, L, 0};
Point(4) = {0, L, 0};

//---------------------------------------------------------------
// Define lines connecting the points.
// Syntax: Line(ID) = {startPoint, endPoint};
// The IDs must match previously defined points.
//---------------------------------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//---------------------------------------------------------------
// Define a closed loop of curves (lines) and then a surface.
// Syntax:
//   Curve Loop(ID) = {lineIDs in order};
//   Plane Surface(ID) = {loopID};
// The order of line IDs must follow the geometric boundary.
//---------------------------------------------------------------
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//---------------------------------------------------------------
// Define Physical Groups to assign names to boundaries or domains.
// These names are used by FEM solvers to identify regions or BCs.
// Syntax:
//   Physical Surface("name") = {surfaceIDs};
//   Physical Curve("name")   = {lineIDs};
//---------------------------------------------------------------
Physical Surface("Domain") = {1};
Physical Curve("Dirichlet")= {1,2,3};
Physical Curve("Neuman")   = {4};

//---------------------------------------------------------------
// Mesh options
//---------------------------------------------------------------
Mesh.ElementOrder = 1; // 1 = linear elements; 2 = quadratic
Mesh.Algorithm = 8;     // mesh generation algorithm (8 = Delaunay)
Mesh.CharacteristicLengthMax = ms; // uniform mesh sizing
Mesh.CharacteristicLengthMin = ms;
Mesh.Optimize = 1; // run optimization to improve element quality

// EOF
