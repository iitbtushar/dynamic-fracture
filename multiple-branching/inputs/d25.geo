//Arakawa and Takashi(1991)
l  = 150.;
w  = 120.;
a  = 8.0; //notch length
b  = 1.0; //notch width 

//refinement
//length of refinement
c0  = 0.5; c1  = 15.0; c2 = 65;
//height of refinement
d  = 5.0; d1 = 70.0; d2 = 70.0;

bb = 60;//vetical distance betn pin-loading
aa = 25;//lever-arm
dd = 20.0;
v = 1.0;//v-notch

// mesh size
l0  = 0.1;
h = 30*l0; h1 = 10*l0;
hh = l0/4;

//points
Point(1)  = {     0,          0,  0, h};
Point(2)  = { a - c0,   0,  0, h};
Point(3)  = {     l,          0,  0, h};
Point(4)  = {     l,          w,  0, h};
Point(5)  = { a - c0,   w,  0, h};
Point(6)  = {     0,          w,  0, h};

Point(7)  = { 0,   0.5*w + b/2.+ v,  0, h1};
Point(8)  = { a - c0,   0.5*w + b/2.,  0, hh};
Point(9)  = { a,   0.5*w + b/2.,  0, hh};
Point(10)  = { a,   0.5*w - b/2.,  0, hh};
Point(101)  = { a,   0.5*w,  0, hh};
Point(102)  = { a+0.5*b,   0.5*w,  0, hh};

Point(11)  = { a - c0,   0.5*w - b/2.,  0, hh};
Point(12)  = { 0,   0.5*w - b/2.- v,  0, h1};

Point(13)  = { a - c0,   0.5*w + d/2.,  0, hh};
Point(14)  = { a + c1,   0.5*w + d/2.,  0, hh};
Point(15)  = { a + c1,   0.5*w - d/2.,  0, hh};
Point(16)  = { a - c0,   0.5*w - d/2.,  0, hh};

Point(21)  = {     l,          0.5*w + d1/2.,  0, h};
Point(22)  = {     l,          0.5*w - d1/2.,  0, h};

Point(40)  = {     a + aa,  w/2. - bb/2.,  0, h1};
Point(41)  = {     a + aa + dd/2.,  w/2. - bb/2.,  0, h1};
Point(42)  = {     a + aa - dd/2.,  w/2. - bb/2.,  0, h1};

Point(50)  = {     a + aa,  w/2. + bb/2.,  0, h1};
Point(51)  = {     a + aa + dd/2.,  w/2. + bb/2.,  0, h1};
Point(52)  = {     a + aa - dd/2.,  w/2. + bb/2.,  0, h1};

Point(61)  = { a + c2,   0.5*w + d1/2.,  0, hh};
Point(62)  = { a + c2,   0.5*w - d1/2.,  0, hh};


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
//Line(93) = {7,12};
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

Line(21) = {14,  61}; Line(210) = {61,  21};
Line(22) = {15,  62}; Line(220) = {62,  22};
Line(61) = {61,  62};

Circle(41) = {41, 40, 42}; Circle(42) = {42, 40, 41};
Circle(51) = {51, 50, 52}; Circle(52) = {52, 50, 51};

Line(43) = {42,  40}; Line(44) = {40,  41};
Line(53) = {52,  50}; Line(54) = {50,  51};


// -------------------
//  Surfaces
// -------------------

Line Loop(101) = {1,-18,17,11,12}; Plane Surface(101) = {101};
Line Loop(102) = {5,6,7,13,19}; Plane Surface(102) = {102};
Line Loop(103) = {220,31,-210,61}; Plane Surface(103) = {103};
Line Loop(1030) = {22,-61,-21,15}; Plane Surface(1030) = {1030};
Line Loop(104) = {13,14,15,16,17,-10,-91,-92,-8}; Plane Surface(104) = {104};
//Line Loop(104) = {13,14,15,16,17,-10,9,-8}; Plane Surface(104) = {104};
//Line Loop(105) = {7,8,91,92,10,11,-93}; Plane Surface(105) = {105};

Line Loop(106) = {2,3,-22,-220,16,18}; 
Line Loop(107) = {14,21,210,32,4,-19};
Line Loop(108) = {41,43,44}; Plane Surface(108) = {108};
Line Loop(109) = {-42,43,44}; Plane Surface(109) = {109};
Line Loop(110) = {51,53,54}; Plane Surface(110) = {110};
Line Loop(111) = {-52,53,54}; Plane Surface(111) = {111};

Plane Surface(106) = {106, 108, 109};
Plane Surface(107) = {107, 110, 111};

//for Q4/Q8 elements
Mesh.RandomFactor = 1e-8;//default = 1e-9
//Recombine Surface{101,102,103,104};//T3->T4
//Mesh.ElementOrder = 2;//T4->Q9 //Mesh.SecondOrderLinear = 0;
//Mesh.SecondOrderIncomplete=1;//Q9->Q8

// ----------------------
// Physical quantities
// ----------------------
Physical Surface("bulk1") = {101, 102, 103, 1030, 104, 106, 107}; //bulk1
Physical Surface("hinge") = {108, 109, 110, 111}; //hinge
/*
<ElementGroup name="bulk">
bulk1 | hinge
</ElementGroup>
*/

Physical Point("bottom")    = {40}; // bottom
Physical Point("top")    = {50}; // top
