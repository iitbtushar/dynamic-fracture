//gmsh geometry for Zhou problem
l  = 32.;
w  = 16.;
a  = 1.0;  //notch length
b  = 0.1; //notch width 

//refinement
c1  = 0.125; //length of refinement
c2  = 10.0; //length of refinement
d  = 10.0; //height of refinement

// mesh size
l0  = 0.1;
h = 5*l0;
hh = l0/4;

//points
Point(1)  = {     0,          0,  0, h};
Point(2)  = { a - c1,   0,  0, h};
Point(3)  = {     l,          0,  0, h};
Point(4)  = {     l,          w,  0, h};
Point(5)  = { a - c1,   w,  0, h};

Point(6)  = {     0,          w,  0, h};
Point(7)  = { 0,   0.5*w + b/2.,  0, hh};
Point(8)  = { a - c1,   0.5*w + b/2.,  0, hh};
Point(9)  = { a,   0.5*w + b/2.,  0, hh};
Point(10)  = { a,   0.5*w - b/2.,  0, hh};
Point(101)  = { a,   0.5*w,  0, hh};
Point(102)  = { a+0.5*b,   0.5*w,  0, hh};

Point(11)  = { a - c1,   0.5*w - b/2.,  0, hh};
Point(12)  = { 0,   0.5*w - b/2.,  0, hh};

Point(13)  = { a - c1,   0.5*w + d/2.,  0, hh};
Point(14)  = { a + c2,   0.5*w + d/2.,  0, hh};
Point(15)  = { a + c2,   0.5*w - d/2.,  0, hh};
Point(16)  = { a - c1,   0.5*w - d/2.,  0, hh};

Point(21)  = {     l,          0.5*w + d/2.,  0, hh};
Point(22)  = {     l,          0.5*w - d/2.,  0, hh};


//lines
Line(1)  = { 1,  2};
Line(2)  = { 2,  3};
Line(3)  = { 3, 22};
Line(31)  = { 22, 21};
Line(32)  = { 21, 4};
Line(4)  = { 4,  5};
Line(5)  = { 5,  6};
Line(6)  = { 6,  7};
Line(7)  = { 7,  8};
Line(8)  = { 8,  9};
Line(91) = {9, 102};
Line(92) = {102, 10};
Line(93) = {7,12};
//Circle(9) = {10, 101, 9};
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

Line(21) = {14,  21};
Line(22) = {15,  22};




// -------------------
//  Surfaces
// -------------------

Line Loop(101) = {1,-18,17,11,12}; Plane Surface(101) = {101};
Line Loop(102) = {5,6,7,13,19}; Plane Surface(102) = {102};
Line Loop(103) = {22,31,-21,15}; Plane Surface(103) = {103};
Line Loop(104) = {13,14,15,16,17,-10,-91,-92,-8}; Plane Surface(104) = {104};
//Line Loop(104) = {13,14,15,16,17,-10,9,-8}; Plane Surface(104) = {104};
Line Loop(105) = {7,8,91,92,10,11,-93}; Plane Surface(105) = {105};
Line Loop(106) = {2,3,-22,16,18}; Plane Surface(106) = {106};
Line Loop(107) = {14,21,32,4,-19}; Plane Surface(107) = {107};

//Point{121} In Surface{201}; 

//for Q4/Q8 elements
//Mesh.RandomFactor = 1e-8;//default = 1e-9
//Recombine Surface{101,102,103,104};//T3->T4
//Mesh.ElementOrder = 2;//T4->Q9 //Mesh.SecondOrderLinear = 0;
//Mesh.SecondOrderIncomplete=1;//Q9->Q8

// ----------------------
// Physical quantities
// ----------------------
Physical Surface("bulk1") = {101, 102, 103, 104, 106, 107}; //bulk1
Physical Surface("notch") = {105}; //notch
/*
<ElementGroup name="bulk">
bulk1 | notch
</ElementGroup>
*/

Physical Line("bottom")    = {1,2}; // bottom
Physical Line("top")    = {4,5}; // top
