//gmsh geometry for branching problem
l  = 100.;
w  = 40.;
a  = 50.;  //notch length
b  = 0.5; //notch width 

//refinement
c1  = 2.0; //length of refinement
c2  = 15.0; //length of refinement
d  = 5.0; //height of refinement

// mesh size
l0  = 0.5;
h = 1.0;
hh = 0.1;

//points
Point(1)  = {     0,          0,  0, h};
Point(2)  = { 0.5*l - c1,   0,  0, hh};
Point(3)  = {     l,          0,  0, hh};
Point(4)  = {     l,          w,  0, hh};
Point(5)  = { 0.5*l - c1,   w,  0, hh};

Point(6)  = {     0,          w,  0, h};
Point(7)  = { 0,   0.5*w + b/2.,  0, h};
Point(8)  = { 0.5*l - c1,   0.5*w + b/2.,  0, hh};
Point(9)  = { 0.5*l,   0.5*w + b/2.,  0, hh};
Point(10)  = { 0.5*l,   0.5*w - b/2.,  0, hh};
Point(11)  = { 0.5*l - c1,   0.5*w - b/2.,  0, hh};
Point(12)  = { 0,   0.5*w - b/2.,  0, h};

Point(13)  = { 0.5*l - c1,   0.5*w + d/2.,  0, hh};
Point(14)  = { 0.5*l + c2,   0.5*w + d/2.,  0, hh};
Point(15)  = { 0.5*l + c2,   0.5*w - d/2.,  0, hh};
Point(16)  = { 0.5*l - c1,   0.5*w - d/2.,  0, hh};



//lines
Line(1)  = { 1,  2};
Line(2)  = { 2,  3};
Line(3)  = { 3,  4};
Line(4)  = { 4,  5};
Line(5)  = { 5,  6};
Line(6)  = { 6,  7};
Line(7)  = { 7,  8};
Line(8)  = { 8,  9};
Line(9) = {9, 10};
Line(10) = {10,  11};
Line(11) = {11,  12};
Line(12) = {12,  1};
Line(13) = {8,  13};
Line(14) = {13,  14};
Line(15) = {14,  15};
Line(16) = {15,  16};
Line(17) = {16,  11};

Line(18) = {16,  2};
Line(19) = {13,  5};


// -------------------
//  Surfaces
// -------------------

Line Loop(101) = {1,-18,17,11,12}; Plane Surface(101) = {101};
Line Loop(102) = {5,6,7,13,19}; Plane Surface(102) = {102};
Line Loop(103) = {2,3,4,-19,14,15,16,18}; Plane Surface(103) = {103};
Line Loop(104) = {13,14,15,16,17,-10,-9,-8}; Plane Surface(104) = {104};

//Point{121} In Surface{201}; 

//for Q4/Q8 elements
//Mesh.RandomFactor = 1e-8;//default = 1e-9
//Recombine Surface{101,102,103,104};//T3->T4
//Mesh.ElementOrder = 2;//T4->Q9 //Mesh.SecondOrderLinear = 0;
//Mesh.SecondOrderIncomplete=1;//Q9->Q8

// ----------------------
// Physical quantities
// ----------------------
Physical Surface("bulk") = {101, 102, 103, 104}; //bulk

Physical Line("bottom")    = {1,2}; // bottom
Physical Line("top")    = {4,5}; // top
