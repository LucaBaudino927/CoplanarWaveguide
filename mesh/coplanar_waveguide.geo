SetFactory("OpenCASCADE");

lc = 1e+2; //um

//--------------------------Points--------------------------


//Kapton points (da estrudere verso il basso)
//Air (da estrudere verso l'alto)
Point(1) = {0,		-1000,		0,		lc};
Point(2) = {20000,	-1000,		0,		lc};
Point(3) = {20000, 	1000,		0,		lc};
Point(4) = {0, 		1000,		0,		lc};

//----------------------------------------------------------

//--------------------------Lines--------------------------


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


//----------------------------------------------------------

//--------------------Curves and Surfaces--------------------


Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};


//----------------------------------------------------------

//-------------------------Volumes-------------------------


//Creo volume 1
Extrude {0, 0, -100} { Surface{1}; }

//Creo volume 2
Extrude {0, 0, 100} { Surface{1}; }


//----------------------------------------------------------


//--------------------------Points--------------------------


p = newp;

//Signal points
Point(p+1) = {0,	0,		0,		lc};
Point(p+2) = {20000,	0,		0,		lc};
Point(p+3) = {20000, 	100,		0,		lc};
Point(p+4) = {0, 	100,		0,		lc};
//Point{p+1}  In Surface{1};
//Point{p+2} In Surface{1};
//Point{p+3} In Surface{1};
//Point{p+4} In Surface{1};

//Gnd1 points
Point(p+5) = {0,	-200,		0,		lc};
Point(p+6) = {20000,	-200,		0,		lc};
Point(p+7) = {20000, 	-100,		0,		lc};
Point(p+8) = {0, 	-100,		0,		lc};
//Point{p+5} In Surface{1};
//Point{p+6} In Surface{1};
//Point{p+7} In Surface{1};
//Point{p+8} In Surface{1};

//Gnd2 points
Point(p+9)  = {0,	200,		0,		lc};
Point(p+10) = {20000,	200,		0,		lc};
Point(p+11) = {20000, 	300,		0,		lc};
Point(p+12) = {0, 	300,		0,		lc};
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

//Creo volume 3
Extrude {0, 0, 20} { Surface{s+1}; }

//Creo volume 4
Extrude {0, 0, 20} { Surface{s+2}; }

//Creo volume 5
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

Physical Surface("boundary_surface_signal", 2001) = {35, 38, 40, 42}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd1", 2002) = {43, 45, 47, 36}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("boundary_surface_gnd2", 2003) = {48, 50, 52, 37}; //NB: trovare un modo per automatizzare questa operazione


//---------------------------Ports--------------------------

//Definisco le superfici delle porte

Physical Surface("signal_port_in", 2004)  = {41}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("signal_port_out", 2005) = {39}; //NB: trovare un modo per automatizzare questa operazione

Physical Surface("gnd1_port_in", 2006)  = {46}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("gnd1_port_out", 2007) = {44}; //NB: trovare un modo per automatizzare questa operazione

Physical Surface("gnd2_port_in", 2008)  = {51}; //NB: trovare un modo per automatizzare questa operazione
Physical Surface("gnd2_port_out", 2009) = {49}; //NB: trovare un modo per automatizzare questa operazione


//----------------------------------------------------------

