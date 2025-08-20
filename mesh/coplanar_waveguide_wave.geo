SetFactory("OpenCASCADE");

//mesh size
large = 200; //200 um
small = 20;  //20 um
lc = 100; //100 um

//--------------------------Points--------------------------


//Kapton points (da estrudere verso il basso)
Point(1) = {0,		-1000,		0,		lc};
Point(2) = {8000,	-1000,		0,		lc};
Point(3) = {8000, 	1000,		0,		lc};
Point(4) = {0, 		1000,		0,		lc};

//Air (da estrudere verso l'alto)
Point(5) = {-200,	-1200,		-225,		lc};
Point(6) = {8200,	-1200,		-225,		lc};
Point(7) = {8200, 	1200,		-225,		lc};
Point(8) = {-200, 	1200,		-225,		lc};

//----------------------------------------------------------

//--------------------------Lines--------------------------

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


//----------------------------------------------------------

//--------------------Curves and Surfaces--------------------

//Kapton
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//air
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

//----------------------------------------------------------

//-------------------------Volumes-------------------------


//Creo volume 1 (kapton)
Extrude {0, 0, -25} { Surface{1}; }

//Creo volume 2 (air)
Extrude {0, 0, 550} { Surface{2}; }


//----------------------------------------------------------


//--------------------------Points--------------------------


p = newp;

//Signal points
Point(p+1) = {0,	0,		0,		lc};
Point(p+2) = {8000,	0,		0,		lc};
Point(p+3) = {8000, 	90,		0,		lc};
Point(p+4) = {0, 	90,		0,		lc};
//Point{p+1}  In Surface{1};
//Point{p+2} In Surface{1};
//Point{p+3} In Surface{1};
//Point{p+4} In Surface{1};

//Gnd1 points
Point(p+5) = {0,	-190,		0,		lc};
Point(p+6) = {8000,	-190,		0,		lc};
Point(p+7) = {8000, 	-100,		0,		lc};
Point(p+8) = {0, 	-100,		0,		lc};
//Point{p+5} In Surface{1};
//Point{p+6} In Surface{1};
//Point{p+7} In Surface{1};
//Point{p+8} In Surface{1};

//Gnd2 points
Point(p+9)  = {0,	190,		0,		lc};
Point(p+10) = {8000,	190,		0,		lc};
Point(p+11) = {8000, 	280,		0,		lc};
Point(p+12) = {0, 	280,		0,		lc};
//Point{p+9}  In Surface{1};
//Point{p+10} In Surface{1};
//Point{p+11} In Surface{1};
//Point{p+12} In Surface{1};


//----------------------------------------------------------

//--------------------------Lines--------------------------

l = newc;

//Signal lines
Line(l+1)  = {p+1,  p+2};
Line(l+2)  = {p+2,  p+3};
Line(l+3)  = {p+3,  p+4};
Line(l+4)  = {p+4,  p+1};
//Curve{l+1} In Surface{1};
//Curve{l+2} In Surface{1};
//Curve{l+3} In Surface{1};
//Curve{l+4} In Surface{1};

//Gnd1 lines
Line(l+5)  = {p+5,  p+6};
Line(l+6)  = {p+6,  p+7};
Line(l+7)  = {p+7,  p+8};
Line(l+8)  = {p+8,  p+5};
//Curve{l+5} In Surface{1};
//Curve{l+6} In Surface{1};
//Curve{l+7} In Surface{1};
//Curve{l+8} In Surface{1};

//Gnd2 lines
Line(l+9)   = {p+9,   p+10};
Line(l+10)  = {p+10,  p+11};
Line(l+11)  = {p+11,  p+12};
Line(l+12)  = {p+12,  p+9};
//Curve{l+9}  In Surface{1};
//Curve{l+10} In Surface{1};
//Curve{l+11} In Surface{1};
//Curve{l+12} In Surface{1};


//----------------------------------------------------------

//--------------------Curves and Surfaces-------------------

cl = newcl;
s = news;

Curve Loop(cl+1) = {l+1, l+2, l+3, l+4};
Plane Surface(s+1) = {cl+1};
//Surface{s+1} In Volume {1};

Curve Loop(cl+2) = {l+5, l+6, l+7, l+8};
Plane Surface(s+2) = {cl+2};
//Surface{s+2} In Volume {1};

Curve Loop(cl+3) = {l+9, l+10, l+11, l+12};
Plane Surface(s+3) = {cl+3};
//Surface{s+3} In Volume {1};


//----------------------------------------------------------

//-------------------------Volumes-------------------------

//Creo volume 3 (signal)
Extrude {0, 0, 20} { Surface{s+1}; }

//Creo volume 4 (gnd1)
Extrude {0, 0, 20} { Surface{s+2}; }

//Creo volume 5 (gnd2)
Extrude {0, 0, 20} { Surface{s+3}; }


//----------------------------------------------------------












///////////////PHYSICAL SURFACES AND VOLUMES/////////////////////

//--------------------------Volumes---------------------------

Physical Volume("kapton", 3001) = {1};
Physical Volume("air", 3002) = {2};
Physical Volume("signal", 3003) = {3};
Physical Volume("gnd1", 3004) = {4};
Physical Volume("gnd2", 3005) = {5};


//-------------------Signal, gdn buondaries------------------

//Definisco le superfici di contorno del signal e dei gnd escludendo le superfici assegnate alle porte

Physical Surface("boundary_surface_signal", 2001) = {39, 42, 44, 46}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd1", 2002) = {40, 47, 49, 51}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd2", 2003) = {41, 52, 54, 56}; //NB: trovare un modo per automatizzare questa operazione


//---------------------------Ports--------------------------

//Definisco le superfici delle porte

Physical Surface("signal_port_in", 2004)  = {45}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("signal_port_out", 2005) = {43}; //NB: trovare un modo per automatizzare questa operazione

Physical Surface("gnd1_port_in", 2006)  = {50}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("gnd1_port_out", 2007) = {48}; //NB: trovare un modo per automatizzare questa operazione

Physical Surface("gnd2_port_in", 2008)  = {55}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("gnd2_port_out", 2009) = {53}; //NB: trovare un modo per automatizzare questa operazione

//--------------------------Farfield-------------------------

Physical Surface("farfield", 2010) = {2, 8, 9, 10, 11, 12}; //NB: trovare un modo per automatizzare questa operazione




//----------------------------------------------------------







/*
//////////////////////////FIELD/////////////////////////////////////

Field[1] = Distance;
Field[1].CurvesList = {l+1, l+2, l+3, l+4, l+5, l+6, l+7, l+8, l+9, l+10, l+11, l+12};
Field[1].Sampling = 100;

// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = small;
Field[2].SizeMax = lc;
Field[2].DistMin = 100;
Field[2].DistMax = 120;

Field[3] = Min;
Field[3].FieldsList = {1, 2};
Background Field = 3;


Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// This will prevent over-refinement due to small mesh sizes on the boundary.

// Finally, while the default "Frontal-Delaunay" 2D meshing algorithm
// (Mesh.Algorithm = 6) usually leads to the highest quality meshes, the
// "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size
// fields better - in particular size fields with large element size gradients:

Mesh.Algorithm = 5;
*/
