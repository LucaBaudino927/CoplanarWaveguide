/*
// Coplanar Waveguide with Lumped Ports - Gmsh .geo version
SetFactory("OpenCASCADE");

// Parameters
trace_width = 90;
gap_width = 100;
boundary_distance = 765;
substrate_height = 25;
metal_height = 20;
length = 8000;

total_width = 2000;

// -------------------------
// Geometry
// -------------------------

// CPW cross-section (from bottom to top)
dy = 0;
Rectangle(1) = {0, dy, 0, length, boundary_distance};
dy += boundary_distance;
Rectangle(2) = {0, dy, 0, length, trace_width}; //gnd1
dy += trace_width;
Rectangle(3) = {0, dy, 0, length, gap_width};
dy += gap_width;
Rectangle(4) = {0, dy, 0, length, trace_width}; //signal
dy += trace_width;
Rectangle(5) = {0, dy, 0, length, gap_width};
dy += gap_width;
Rectangle(6) = {0, dy, 0, length, trace_width}; //gnd2
dy += trace_width;
Rectangle(7) = {0, dy, 0, length, boundary_distance};

// -------------------------
// Ports as rectangles
// -------------------------

//------------------------------------------------
//
//
//------------------------------------------------
//			gnd2
//------------------------------------------------
//| |p1b                                    p2b| |
//------------------------------------------------
//			signal
//------------------------------------------------
//| |p1a                                    p2a| |
//------------------------------------------------
//			gnd1
//------------------------------------------------
//
//
//------------------------------------------------
py = 0 + boundary_distance + trace_width;
Rectangle(8) = {0, py, 0, gap_width, gap_width};                       // p1a
Rectangle(9) = {length - gap_width, py, 0, gap_width, gap_width};      // p2a
py += gap_width + trace_width;
Rectangle(10) = {0, py, 0, gap_width, gap_width};                       // p1b
Rectangle(11) = {length - gap_width, py, 0, gap_width, gap_width};      // p2b


// kapton
Box(1) = {0, 0, -substrate_height, length, total_width, substrate_height};

// air
Box(2) = {-200, -200, -substrate_height-200, length + 2*200, total_width + 2*200, 2*(substrate_height + 200)};

// signal
Box(3) = {0, boundary_distance + trace_width + gap_width, 0, length, trace_width, metal_height};

// gnd1
Box(4) = {0, boundary_distance, 0, length, trace_width, metal_height};

// gnd2
Box(5) = {0, boundary_distance + trace_width + gap_width + trace_width + gap_width, 0, length, trace_width, metal_height};


// -------------------------
// Physical groups
// -------------------------

Physical Volume("kapton", 3001) = {1};
Physical Volume("air",    3002) = {2};
Physical Volume("signal", 3003) = {3};
Physical Volume("gnd1",   3004) = {4};
Physical Volume("gnd2",   3005) = {5};

Physical Surface("boundary_surface_signal", 2001) = {24, 25, 26, 27, 28, 29}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd1", 2002)   = {30, 31, 32, 33, 34, 35}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd2", 2003)   = {36, 37, 38, 39, 40, 41}; //NB: trovare un modo per automatizzare questa operazione

Physical Surface("p1a", 2004) = {8};
Physical Surface("p2a", 2005) = {9};
Physical Surface("p1b", 2006) = {10};
Physical Surface("p2b", 2007) = {11};

Physical Surface("farfield", 2008) = {18, 19, 20, 21, 22, 23}; //NB: trovare un modo per automatizzare questa operazione

// -------------------------
// Mesh controls
// -------------------------
Mesh.CharacteristicLengthMin = 30;
Mesh.CharacteristicLengthMax = 500;

// Refinement near CPW
Field[1] = Distance;
Field[1].SurfacesList = {24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41}; // gap and trace
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 30;
Field[2].SizeMax = 500;
Field[2].DistMin = trace_width;
Field[2].DistMax = 2 * gap_width;

Field[3] = Min;
Field[3].FieldsList = {2};

Background Field = 3;

Mesh.Algorithm3D = 10; // Delaunay

// -------------------------
// Mesh generation
// -------------------------
//Mesh 3;

*/

// Coplanar Waveguide with Lumped Ports - Gmsh .geo version
SetFactory("OpenCASCADE");

// Parameters
trace_width = 90;
gap_width = 100;
boundary_distance = 765;
substrate_height = 25;
metal_height = 20;
length = 8000;
total_width = 2000;
air_gap = 200;

// -------------------------
// Geometry
// -------------------------

// CPW cross-section (from bottom to top)
//kapton points
Point(1) = {0, 			0, 			0};
Point(2) = {length, 		0, 			0};
Point(3) = {length, 		total_width, 		0};
Point(4) = {0, 			total_width, 		0};
//air points
Point(5) = {-air_gap, 		-air_gap, 		-air_gap};
Point(6) = {length+air_gap,	-air_gap, 		-air_gap};
Point(7) = {length+air_gap, 	total_width+air_gap,	-air_gap};
Point(8) = {-air_gap, 		total_width+air_gap, 	-air_gap};
//gnd1
Point(9)  = {0, 		boundary_distance, 		0};
Point(10) = {length,		boundary_distance, 		0};
Point(11) = {length, 		boundary_distance+trace_width,	0};
Point(12) = {0, 		boundary_distance+trace_width, 	0};
//signal
Point(13) = {0, 		boundary_distance+trace_width+gap_width,		0};
Point(14) = {length,		boundary_distance+trace_width+gap_width, 		0};
Point(15) = {length, 		boundary_distance+trace_width+gap_width+trace_width,	0};
Point(16) = {0, 		boundary_distance+trace_width+gap_width+trace_width, 	0};
//gnd2
Point(17) = {0, 		boundary_distance+trace_width+gap_width+trace_width+gap_width,			0};
Point(18) = {length,		boundary_distance+trace_width+gap_width+trace_width+gap_width, 			0};
Point(19) = {length, 		boundary_distance+trace_width+gap_width+trace_width+gap_width+trace_width,	0};
Point(20) = {0, 		boundary_distance+trace_width+gap_width+trace_width+gap_width+trace_width, 	0};

// -------------------------
// Ports as rectangles
// -------------------------

//------------------------------------------------
//
//
//------------------------------------------------
//			gnd2
//------------------------------------------------
//| |p1b                                    p2b| |
//------------------------------------------------
//			signal
//------------------------------------------------
//| |p1a                                    p2a| |
//------------------------------------------------
//			gnd1
//------------------------------------------------
//
//
//------------------------------------------------

/*
//auxilliary points for p1a
Point(21) = {gap_width, 	boundary_distance+trace_width, 			0};
Point(22) = {gap_width,		boundary_distance+trace_width+gap_width, 	0};
//auxilliary points for p2a
Point(23) = {length-gap_width, 	boundary_distance+trace_width, 			0};
Point(24) = {length-gap_width,	boundary_distance+trace_width+gap_width, 	0};
//auxilliary points for p1b
Point(25) = {gap_width, 	boundary_distance+trace_width+gap_width+trace_width, 		0};
Point(26) = {gap_width,		boundary_distance+trace_width+gap_width+trace_width+gap_width, 	0};
//auxilliary points for p2b
Point(27) = {length-gap_width, 	boundary_distance+trace_width+gap_width+trace_width, 		0};
Point(28) = {length-gap_width,	boundary_distance+trace_width+gap_width+trace_width+gap_width, 	0};
*/



//kapton
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//air
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

//signal
Line(9)  = {13, 14};
Line(10) = {14, 15};
Line(11) = {15, 16};
Line(12) = {16, 13};

//gnd1
Line(13) = {9,  10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};

//gnd2
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 17};
/*
//p1a
Line(21) = {12, 21};
Line(22) = {21, 22};
Line(23) = {22, 13};
Line(24) = {13, 12};

//p2a
Line(25) = {23, 11};
Line(26) = {11, 14};
Line(27) = {14, 24};
Line(28) = {24, 23};

//p1b
Line(29) = {16, 25};
Line(30) = {25, 26};
Line(31) = {26, 17};
Line(32) = {17, 16};

//p2b
Line(33) = {27, 15};
Line(34) = {15, 18};
Line(35) = {18, 28};
Line(36) = {28, 27};
*/


//Kapton
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//air
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

//signal
Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3};

//gnd1
Curve Loop(4) = {13, 14, 15, 16};
Plane Surface(4) = {4};

//gnd2
Curve Loop(5) = {17, 18, 19, 20};
Plane Surface(5) = {5};
/*
//p1a
Curve Loop(6) = {21, 22, 23, 24};
Plane Surface(6) = {6};

//p2a
Curve Loop(7) = {25, 26, 27, 28};
Plane Surface(7) = {7};

//p1b
Curve Loop(8) = {29, 30, 31, 32};
Plane Surface(8) = {8};

//p2b
Curve Loop(9) = {33, 34, 35, 36};
Plane Surface(9) = {9};
*/


Rectangle(6) = {0, boundary_distance+trace_width, 0, gap_width,	gap_width, 0};



//Creo volume 1 (kapton)
Extrude {0, 0, -substrate_height} { Surface{1}; }

//Creo volume 2 (air)
Extrude {0, 0, 2*substrate_height + 2*air_gap} { Surface{2}; }

//Creo volume 3 (signal)
Extrude {0, 0, metal_height} { Surface{3}; }

//Creo volume 4 (gnd1)
Extrude {0, 0, metal_height} { Surface{4}; }

//Creo volume 5 (gnd2)
Extrude {0, 0, metal_height} { Surface{5}; }



// -------------------------
// Physical groups
// -------------------------

Physical Volume("kapton", 3001) = {1};
Physical Volume("air",    3002) = {2};
Physical Volume("signal", 3003) = {3};
Physical Volume("gnd1",   3004) = {4};
Physical Volume("gnd2",   3005) = {5};

Physical Surface("boundary_surface_signal", 2001) = {3, 20, 21, 22, 23, 24}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd1", 2002)   = {4, 25, 26, 27, 28, 29}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd2", 2003)   = {5, 30, 31, 32, 33, 34}; //NB: trovare un modo per automatizzare questa operazione

Physical Surface("p1a", 2004) = {6, 32, 33, 34, 35, 36};
//Physical Surface("p2a", 2005) = {7};
//Physical Surface("p1b", 2006) = {8};
//Physical Surface("p2b", 2007) = {9};

Physical Surface("farfield", 2008) = {2, 15, 16, 17, 18, 19}; //NB: trovare un modo per automatizzare questa operazione

// -------------------------
// Mesh controls
// -------------------------
Mesh.CharacteristicLengthMin = 30;
Mesh.CharacteristicLengthMax = 500;

// Refinement near CPW
Field[1] = Distance;
Field[1].SurfacesList = {3, 20, 21, 22, 23, 24, 4, 25, 26, 27, 28, 29, 5, 30, 31, 32, 33, 34}; // gap and trace
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 30;
Field[2].SizeMax = 500;
Field[2].DistMin = trace_width;
Field[2].DistMax = 2 * gap_width;

Field[3] = Min;
Field[3].FieldsList = {2};

Background Field = 3;

Mesh.Algorithm3D = 10; // Delaunay

// -------------------------
// Mesh generation
// -------------------------
//Mesh 3;
